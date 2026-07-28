#ifndef _REMORA_NCTIDEDATA_H_
#define _REMORA_NCTIDEDATA_H_

#ifdef REMORA_USE_NETCDF

#include <string>
#include <memory>

#include <AMReX_AmrCore.H>

/** \brief Static (time-invariant) tidal constituent data read from a ROMS tide file
 *
 * This is the REMORA analogue of ROMS's TIDES(ng) structure (ROMS/Modules/mod_tides.F),
 * populated by get_idata.F.  Unlike NCTimeSeries, tidal harmonic constants do not vary
 * in time: the file is read exactly once at initialization and the constituent fields
 * are held for the whole run.  The time dependence enters only through the harmonic
 * evaluation in REMORA::set_tides (ROMS/Nonlinear/set_tides.F).
 *
 * Variables and unit conventions follow ROMS/External/varinfo.dat exactly:
 *
 *   tide_period   hours        -> seconds  (Fscale = 3600)
 *   tide_Eamp     meter        -> meter    (Fscale = 1)
 *   tide_Ephase   degrees      -> radians  (Fscale = pi/180)
 *   tide_Cangle   degrees      -> radians  (Fscale = pi/180)
 *   tide_Cphase   degrees      -> radians  (Fscale = pi/180)
 *   tide_Cmax     meter/second -> m/s      (Fscale = 1)
 *   tide_Cmin     meter/second -> m/s      (Fscale = 1)
 *
 * The 2D fields are stored as MultiFabs with one component per constituent, on the
 * rho (cell-centered) grid, with ghost cells filled from the file's outer rho rows
 * (ROMS's IstrR..IendR, JstrR..JendR halo) in the same way NCTimeSeries fills its
 * grown tile boxes.
 */
class NCTideData
{
    public:
        /** \brief Constructor
         *
         * @param[in] a_file_name  name of the ROMS tide (TIDENAME) file
         * @param[in] a_domain     level-0 problem domain (cell-centered)
         */
        NCTideData (const std::string& a_file_name,
                    const amrex::Box& a_domain);

        /** \brief Read constituent periods and 2D harmonic fields from the file
         *
         * @param[in] ba2d   2D (k-slab) cell-centered BoxArray to store the data on
         * @param[in] dm     DistributionMapping matching ba2d
         * @param[in] ngrow  number of grow cells for the stored MultiFabs
         * @param[in] geom   geometry of the level (used for FillBoundary)
         */
        void Initialize (const amrex::BoxArray& ba2d,
                         const amrex::DistributionMapping& dm,
                         const amrex::IntVect& ngrow,
                         const amrex::Geometry& geom);

        /** \brief Number of tidal constituents found in the file (ROMS NTC) */
        int n_constituents () const { return ntide; }

        /** \brief Device pointer to the constituent periods [s] (ROMS Tperiod) */
        const amrex::Real* period_ptr () const { return period_d.dataPtr(); }

        /** \brief Constituent periods [s] on the host */
        const amrex::Vector<amrex::Real>& periods () const { return period_h; }

        const amrex::MultiFab& Eamp   () const { return *mf_Eamp; }    //!< SSH_Tamp    [m]
        const amrex::MultiFab& Ephase () const { return *mf_Ephase; }  //!< SSH_Tphase  [rad]
        const amrex::MultiFab& Cangle () const { return *mf_Cangle; }  //!< UV_Tangle   [rad]
        const amrex::MultiFab& Cphase () const { return *mf_Cphase; }  //!< UV_Tphase   [rad]
        const amrex::MultiFab& Cmax   () const { return *mf_Cmax; }    //!< UV_Tmajor   [m/s]
        const amrex::MultiFab& Cmin   () const { return *mf_Cmin; }    //!< UV_Tminor   [m/s]

    private:
        /** \brief Read one 2D-per-constituent field, apply the varinfo scale factor,
         *         and store it (one component per constituent) into mf */
        void read_field (const std::string& var_name, amrex::Real scale,
                         amrex::MultiFab& mf, const amrex::Geometry& geom);

        /// File to read from
        std::string file_name;
        /// Level-0 problem domain
        amrex::Box domain;
        /// Number of tidal constituents
        int ntide = 0;
        /// Constituent periods [s], host and device copies
        amrex::Vector<amrex::Real> period_h;
        amrex::Gpu::DeviceVector<amrex::Real> period_d;

        std::unique_ptr<amrex::MultiFab> mf_Eamp;
        std::unique_ptr<amrex::MultiFab> mf_Ephase;
        std::unique_ptr<amrex::MultiFab> mf_Cangle;
        std::unique_ptr<amrex::MultiFab> mf_Cphase;
        std::unique_ptr<amrex::MultiFab> mf_Cmax;
        std::unique_ptr<amrex::MultiFab> mf_Cmin;
};

#endif // REMORA_USE_NETCDF
#endif //_REMORA_NCTIDEDATA_H_
