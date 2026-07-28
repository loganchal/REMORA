/**
 * \file REMORA_NCTideData.cpp
 *
 * Reader for ROMS tidal forcing data (TIDENAME file) and for the ROMS grid
 * rotation angle, which the tidal current ellipses need.
 *
 * ROMS references:
 *   ROMS/Modules/mod_tides.F      TIDES(ng) storage
 *   ROMS/Nonlinear/get_idata.F    reads Tperiod / SSH_T* / UV_T* once at startup
 *   ROMS/External/varinfo.dat     unit scale factors (Fscale) applied on read
 *   ROMS/Utility/get_grid.F       reads "angle" into GRID(ng)%angler
 */

#include "REMORA_NCTideData.H"
#include "REMORA_NCFile.H"
#include "REMORA_Constants.H"
#include "REMORA.H"

#include "AMReX_ParallelDescriptor.H"

#ifdef REMORA_USE_NETCDF

using namespace amrex;

/**
 * @param[in] a_file_name  name of the ROMS tide file
 * @param[in] a_domain     level-0 problem domain
 */
NCTideData::NCTideData (const std::string& a_file_name,
                        const Box& a_domain)
{
    file_name = a_file_name;
    domain    = a_domain;
}

void
NCTideData::Initialize (const BoxArray& ba2d,
                        const DistributionMapping& dm,
                        const IntVect& ngrow,
                        const Geometry& geom)
{
    amrex::Print() << "Loading tidal constituents from NetCDF file " << file_name << std::endl;

    //
    // Constituent periods. varinfo.dat gives 'tide_period' in hours with
    // Fscale = 3600 (60*60), so the stored value is in seconds.
    //
    {
        using RARRAY = NDArray<Real>;
        Vector<RARRAY> array_per(1);
        ReadNetCDFFile(file_name, {"tide_period"}, array_per); // filled only on proc 0
        if (ParallelDescriptor::IOProcessor()) {
            ntide = static_cast<int>(array_per[0].get_vshape()[0]);
            for (int n(0); n < ntide; n++) {
                period_h.push_back( (*(array_per[0].get_data() + n)) * Real(3600.0) );
            }
        }
        int ioproc = ParallelDescriptor::IOProcessorNumber();
        ParallelDescriptor::Bcast(&ntide, 1, ioproc);
        if (!(ParallelDescriptor::IOProcessor())) {
            period_h.resize(ntide);
        }
        ParallelDescriptor::Bcast(period_h.data(), period_h.size(), ioproc);
    }

    if (ntide < 1) {
        Abort("No tidal constituents found in " + file_name);
    }

    period_d.resize(ntide);
    Gpu::copyAsync(Gpu::hostToDevice, period_h.begin(), period_h.end(), period_d.begin());
    Gpu::streamSynchronize();

    amrex::Print() << "Found " << ntide << " tidal constituent(s); periods [h]:";
    for (int n = 0; n < ntide; n++) {
        amrex::Print() << " " << period_h[n] / Real(3600.0);
    }
    amrex::Print() << std::endl;

    mf_Eamp   = std::make_unique<MultiFab>(ba2d, dm, ntide, ngrow);
    mf_Ephase = std::make_unique<MultiFab>(ba2d, dm, ntide, ngrow);
    mf_Cangle = std::make_unique<MultiFab>(ba2d, dm, ntide, ngrow);
    mf_Cphase = std::make_unique<MultiFab>(ba2d, dm, ntide, ngrow);
    mf_Cmax   = std::make_unique<MultiFab>(ba2d, dm, ntide, ngrow);
    mf_Cmin   = std::make_unique<MultiFab>(ba2d, dm, ntide, ngrow);

    // Scale factors are exactly the Fscale column of ROMS/External/varinfo.dat.
    // The degrees-to-radians factor is the literal that ROMS reads from that file
    // (0.017453292519943295), not a recomputed pi/180, so that the scaling is
    // bit-identical to the ROMS build.
    const Real deg2rad = Real(0.017453292519943295);

    read_field("tide_Eamp"  , Real(1.0), *mf_Eamp  , geom);
    read_field("tide_Ephase", deg2rad  , *mf_Ephase, geom);
    read_field("tide_Cangle", deg2rad  , *mf_Cangle, geom);
    read_field("tide_Cphase", deg2rad  , *mf_Cphase, geom);
    read_field("tide_Cmax"  , Real(1.0), *mf_Cmax  , geom);
    read_field("tide_Cmin"  , Real(1.0), *mf_Cmin  , geom);
}

/**
 * @param[in   ] var_name  name of the variable in the tide file
 * @param[in   ] scale     varinfo.dat Fscale for this variable
 * @param[inout] mf        MultiFab (ntide components) to fill
 * @param[in   ] geom      geometry, for FillBoundary
 */
void
NCTideData::read_field (const std::string& var_name, Real scale,
                        MultiFab& mf, const Geometry& geom)
{
    mf.setVal(Real(0.0));

    for (int itide = 0; itide < ntide; itide++)
    {
        FArrayBox NC_fab;
        Vector<FArrayBox*> NC_fabs;
        Vector<std::string> NC_names;
        Vector<enum NC_Data_Dims_Type> NC_dim_types;

        NC_fabs.push_back(&NC_fab); NC_names.push_back(var_name);
        // The leading dimension of the tide fields is the constituent index, which
        // has the same layout as a time dimension, so we read one "record" at a time.
        NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);

        BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, file_name, NC_names, NC_dim_types,
                                                NC_fabs, true, itide);

        // Don't tile: we are operating on full FABs here
        for (MFIter mfi(mf, false); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = mf[mfi];
            // FArrayBox to FArrayBox copy does "copy on intersection". This works
            // because the netCDF FAB has been broadcast to all ranks.
            fab.template copy<RunOn::Device>(NC_fab, 0, itide, 1);
        }
    }

    // Apply the varinfo.dat scale factor, including in the ghost cells that were
    // filled from the file's outer rho rows
    if (scale != Real(1.0)) {
        const int nt = ntide;
        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box gbx = mfi.growntilebox(mf.nGrowVect());
            const Array4<Real>& arr = mf.array(mfi);
            ParallelFor(gbx, nt, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                arr(i,j,k,n) *= scale;
            });
        }
    }

    mf.FillBoundary(geom.periodicity());
}

/**
 * \brief Read the ROMS grid rotation angle ("angle", radians) from the grid file.
 *
 * ROMS: GRID(ng)%angler, read in get_grid.F. The angle between the grid XI-axis
 * and east is needed to rotate tidal current ellipses (which are given in
 * geographic coordinates) into grid coordinates in set_tides.
 */
void
read_angle_from_netcdf (int /*lev*/,
                        const Box& domain,
                        const std::string& fname,
                        FArrayBox& NC_angle_fab)
{
    amrex::Print() << "Loading grid angle from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_angle_fab); NC_names.push_back("angle"); NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE);

    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param[in] lev level to fill the grid angle on
 */
void
REMORA::init_angler_from_netcdf (int lev)
{
    Vector<FArrayBox> NC_angle_fab; NC_angle_fab.resize(num_boxes_at_level[lev]);

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_angle_from_netcdf(lev, boxes_at_level[lev][idx], nc_grid_file[lev][idx],
                               NC_angle_fab[idx]);

        // Don't tile since we are operating on full FABs in this routine
        for (MFIter mfi(*vec_angler[lev], false); mfi.isValid(); ++mfi)
        {
            FArrayBox& angler_fab = (*vec_angler[lev])[mfi];
            angler_fab.template copy<RunOn::Device>(NC_angle_fab[idx]);
        }
    }

    vec_angler[lev]->FillBoundary(geom[lev].periodicity());
}

#endif // REMORA_USE_NETCDF
