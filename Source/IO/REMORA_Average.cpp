/**
 * \file REMORA_Average.cpp
 *
 * ROMS-equivalent time-averaged ("AVERAGES") output.
 *
 * This is a port of the plain running-sum averaging performed by ROMS
 * ROMS/Nonlinear/set_avg.F. See Tests/MoanaParity/avg_port_notes.md for the
 * mapping between the ROMS source line ranges and the code below.
 *
 * Semantics reproduced here:
 *   - one instantaneous sample is added per baroclinic time step
 *     (set_avg.F accumulate block, from line 1201 of the preprocessed Moana
 *      source: `avgzeta = avgzeta + zeta(Kout)`),
 *   - the window holds exactly nAVG samples (set_avg.F:674-677 reset logic),
 *   - at the end of the window the sums are multiplied by fac = 1/nAVG and the
 *     record is stamped with the model time at the end of the window
 *     (set_avg.F:1782-1800, `AVGtime = AVGtime + nAVG*dt`),
 *   - the accumulators are then zeroed for the next window.
 *
 * NOTE (CUDA): the two functions here that launch device lambdas
 * (REMORA::AverageNormalize and REMORA::mask_avg_arrays_for_write) are both
 * *public* members of REMORA, because nvcc rejects extended __device__ lambdas
 * inside private/protected member functions. Accumulation itself uses
 * amrex::MultiFab::Saxpy, which runs on the device without introducing any
 * lambda of our own.
 */

#include <REMORA.H>
#include "REMORA_IndexDefines.H"

using namespace amrex;

/**
 * Allocate the time-average accumulators (once) and zero them.
 *
 * Called lazily from AverageAccumulate so that a run with averaging turned off
 * (remora.avg_int <= 0, the default) allocates nothing and behaves exactly as
 * before.
 */
void
REMORA::AverageInit ()
{
    if (avg_initialized) return;

    // Which cell-centered scalars go into the avg file. We mirror the history
    // file: only scalars the user asked to plot, and only those for which the
    // NetCDF writer has a ROMS-named variable.
    avg_var_names_3d.clear();
    avg_cons_comp.clear();
    for (const auto& name : plot_var_names_3d) {
        if (name == "temp") {
            avg_var_names_3d.push_back(name);
            avg_cons_comp.push_back(Temp_comp);
        } else if (name == "salt") {
            avg_var_names_3d.push_back(name);
            avg_cons_comp.push_back(Salt_comp);
        } else if (name == "tracer") {
            avg_var_names_3d.push_back(name);
            avg_cons_comp.push_back(Tracer_comp);
        }
    }
    const int n_avg_cons = static_cast<int>(avg_var_names_3d.size());

    vec_avg_zeta.resize(finest_level+1);
    vec_avg_ubar.resize(finest_level+1);
    vec_avg_vbar.resize(finest_level+1);
    vec_avg_u.resize(finest_level+1);
    vec_avg_v.resize(finest_level+1);
    vec_avg_cons.resize(finest_level+1);
    vec_avg_sustr.resize(finest_level+1);
    vec_avg_svstr.resize(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        // Each accumulator mirrors the box array, distribution map and ghost
        // region of the field it accumulates, so that the ghost cells the
        // NetCDF writer reads at the domain boundary are averaged too (this is
        // the analogue of the ROMS IstrR:IendR / JstrR:JendR loop bounds).
        vec_avg_zeta[lev] = std::make_unique<MultiFab>(vec_Zt_avg1[lev]->boxArray(), dmap[lev],
                                                       1, vec_Zt_avg1[lev]->nGrowVect());
        vec_avg_ubar[lev] = std::make_unique<MultiFab>(vec_ubar[lev]->boxArray(), dmap[lev],
                                                       1, vec_ubar[lev]->nGrowVect());
        vec_avg_vbar[lev] = std::make_unique<MultiFab>(vec_vbar[lev]->boxArray(), dmap[lev],
                                                       1, vec_vbar[lev]->nGrowVect());
        vec_avg_u[lev]    = std::make_unique<MultiFab>(xvel_new[lev]->boxArray(), dmap[lev],
                                                       1, xvel_new[lev]->nGrowVect());
        vec_avg_v[lev]    = std::make_unique<MultiFab>(yvel_new[lev]->boxArray(), dmap[lev],
                                                       1, yvel_new[lev]->nGrowVect());
        vec_avg_cons[lev] = std::make_unique<MultiFab>(cons_new[lev]->boxArray(), dmap[lev],
                                                       std::max(n_avg_cons,1), cons_new[lev]->nGrowVect());
        vec_avg_sustr[lev] = std::make_unique<MultiFab>(vec_sustr[lev]->boxArray(), dmap[lev],
                                                        1, vec_sustr[lev]->nGrowVect());
        vec_avg_svstr[lev] = std::make_unique<MultiFab>(vec_svstr[lev]->boxArray(), dmap[lev],
                                                        1, vec_svstr[lev]->nGrowVect());
    }

    avg_initialized = true;
    AverageReset();

    Print() << "Time-averaged (avg) output enabled: nAVG = " << avg_int
            << " baroclinic steps per window, file prefix '" << avg_file_name << "'" << std::endl;
}

/**
 * Zero all accumulators and the per-window sample counter.
 *
 * ROMS instead *assigns* the first sample of a window (set_avg.F:674-690) and
 * adds the remaining nAVG-1; zeroing then adding all nAVG samples is
 * arithmetically identical and avoids a separate first-sample branch.
 */
void
REMORA::AverageReset ()
{
    if (!avg_initialized) return;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vec_avg_zeta[lev]->setVal(0.0_rt);
        vec_avg_ubar[lev]->setVal(0.0_rt);
        vec_avg_vbar[lev]->setVal(0.0_rt);
        vec_avg_u[lev]->setVal(0.0_rt);
        vec_avg_v[lev]->setVal(0.0_rt);
        vec_avg_cons[lev]->setVal(0.0_rt);
        vec_avg_sustr[lev]->setVal(0.0_rt);
        vec_avg_svstr[lev]->setVal(0.0_rt);
    }
    avg_nsamples = 0;
}

/**
 * Scale every accumulator, over its valid *and* ghost cells, by fac.
 *
 * This is ROMS set_avg.F:1782-1800 (`fac = 1/REAL(nAVG)` applied to each C-grid
 * variable type). We do it with an explicit ParallelFor over growntilebox()
 * rather than MultiFab::mult, because MultiFab::mult's int-ghost argument would
 * grow the 2D (slab) accumulators in the vertical, past their allocation.
 *
 * Public member: launches device lambdas, which nvcc forbids in
 * private/protected member functions.
 *
 * @param[in] fac  factor to multiply every accumulator by
 */
void
REMORA::AverageNormalize (Real fac)
{
    if (!avg_initialized) return;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        MultiFab* mfs[] = { vec_avg_zeta[lev].get(), vec_avg_ubar[lev].get(), vec_avg_vbar[lev].get(),
                            vec_avg_u[lev].get(),    vec_avg_v[lev].get(),    vec_avg_cons[lev].get(),
                            vec_avg_sustr[lev].get(), vec_avg_svstr[lev].get() };
        for (MultiFab* mf : mfs) {
            const int nc = mf->nComp();
            for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.growntilebox();
                Array4<Real> const& arr = mf->array(mfi);
                ParallelFor(bx, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                {
                    arr(i,j,k,n) *= fac;
                });
            }
        }
    }
    Gpu::streamSynchronize();
}

/**
 * Add one instantaneous sample of every averaged field to its running sum.
 *
 * This is the ROMS set_avg.F accumulation block (from line 1201 of the
 * preprocessed Moana source). It is called once per baroclinic step, at the
 * same point in the step cycle at which REMORA writes a history record, so
 * that the sampled state is identical to what a history file written every
 * step would contain.
 *
 * All work is a device Saxpy over valid + ghost cells; nothing is copied to the
 * host and no communication is required, so this is correct for any MPI
 * decomposition and on GPU.
 */
void
REMORA::AverageAccumulate ()
{
    if (avg_int <= 0) return;

    AverageInit();

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // zeta: REMORA writes vec_Zt_avg1 as "zeta" in the history file, so that
        // is what we accumulate (ROMS uses zeta(:,:,Kout) = zeta(:,:,kstp)).
        MultiFab::Saxpy(*vec_avg_zeta[lev], 1.0_rt, *vec_Zt_avg1[lev],
                        0, 0, 1, vec_avg_zeta[lev]->nGrowVect());

        // ubar / vbar: component 0 is the time level REMORA writes as ubar/vbar
        // (ROMS: ubar(:,:,Kout), vbar(:,:,Kout) with Kout = kstp).
        MultiFab::Saxpy(*vec_avg_ubar[lev], 1.0_rt, *vec_ubar[lev],
                        0, 0, 1, vec_avg_ubar[lev]->nGrowVect());
        MultiFab::Saxpy(*vec_avg_vbar[lev], 1.0_rt, *vec_vbar[lev],
                        0, 0, 1, vec_avg_vbar[lev]->nGrowVect());

        // u / v (ROMS: u(:,:,:,Nout), v(:,:,:,Nout) with Nout = nrhs)
        MultiFab::Saxpy(*vec_avg_u[lev], 1.0_rt, *xvel_new[lev],
                        0, 0, 1, vec_avg_u[lev]->nGrowVect());
        MultiFab::Saxpy(*vec_avg_v[lev], 1.0_rt, *yvel_new[lev],
                        0, 0, 1, vec_avg_v[lev]->nGrowVect());

        // surface momentum stress (ROMS: idUsms/idVsms averages)
        MultiFab::Saxpy(*vec_avg_sustr[lev], 1.0_rt, *vec_sustr[lev],
                        0, 0, 1, vec_avg_sustr[lev]->nGrowVect());
        MultiFab::Saxpy(*vec_avg_svstr[lev], 1.0_rt, *vec_svstr[lev],
                        0, 0, 1, vec_avg_svstr[lev]->nGrowVect());

        // temp / salt / tracer (ROMS: t(:,:,:,Nout,itrc))
        for (int n = 0; n < static_cast<int>(avg_cons_comp.size()); ++n) {
            MultiFab::Saxpy(*vec_avg_cons[lev], 1.0_rt, *cons_new[lev],
                            avg_cons_comp[n], n, 1, vec_avg_cons[lev]->nGrowVect());
        }
    }

    avg_nsamples++;
}

/**
 * Convert the accumulated sums into time averages, write the record, and reset.
 *
 * ROMS set_avg.F:1782-1800: fac = 1/REAL(nAVG) applied to every accumulator,
 * and AVGtime = AVGtime + nAVG*dt (or AVGtime = time when nAVG == 1). Because
 * we sample once per completed step, the model time when the window closes is
 * already t_start_of_window_run + nAVG*dt, so t_new[0] is that timestamp and no
 * separate AVGtime counter is needed.
 *
 * @param[in] which_step  step index at which the window closed (used for file naming)
 */
void
REMORA::AverageWriteAndReset (int which_step)
{
    if (avg_int <= 0 || !avg_initialized) return;

    const Real fac = 1.0_rt / static_cast<Real>(avg_nsamples);

    AverageNormalize(fac);

#ifdef REMORA_USE_NETCDF
    AMREX_ALWAYS_ASSERT(finest_level == 0);
    const int lev = 0;
    vec_avg_cons[lev]->FillBoundary(geom[lev].periodicity());
    WriteNCPlotFile(which_step, vec_avg_cons[lev].get(), true);
    avg_count++;
#else
    amrex::ignore_unused(which_step);
    amrex::Abort("remora.avg_int > 0 requires a NetCDF-enabled build (REMORA_ENABLE_PNETCDF=ON)");
#endif

    AverageReset();
}

/**
 * Accumulate one sample and, if the averaging window has just closed, write the
 * averaged record and start a new window.
 *
 * @param[in] step      number of completed coarse steps
 * @param[in] cur_time  model time after that step
 */
void
REMORA::AverageAtIntermediateTime (int step, Real cur_time)
{
    amrex::ignore_unused(cur_time);
    if (avg_int <= 0) return;

    AverageAccumulate();

    if (avg_nsamples == avg_int) {
        AverageWriteAndReset(step);
    }
}

/**
 * Apply the land/sea mask to the time-average accumulators before/after output,
 * mirroring REMORA::mask_arrays_for_write for the instantaneous state.
 *
 * Land cells accumulate exactly zero every step, so the round trip
 * (0 -> fill_value -> 0) is exact.
 *
 * Public member: launches device lambdas, which nvcc forbids in
 * private/protected member functions.
 *
 * @param[in] lev          level to mask
 * @param[in] fill_value   value written into masked cells
 * @param[in] fill_where   value a face-based cell must currently hold to be masked
 */
void
REMORA::mask_avg_arrays_for_write (int lev, Real fill_value, Real fill_where)
{
    if (!avg_initialized) return;

    const int n_avg_cons = vec_avg_cons[lev]->nComp();

    for (MFIter mfi(*vec_avg_cons[lev],false); mfi.isValid(); ++mfi) {
        Box gbx1 = mfi.growntilebox(IntVect(NGROW+1,NGROW+1,0));
        Box ubx  = mfi.grownnodaltilebox(0,IntVect(NGROW,NGROW,0));
        Box vbx  = mfi.grownnodaltilebox(1,IntVect(NGROW,NGROW,0));

        Array4<Real> const& zeta_a = vec_avg_zeta[lev]->array(mfi);
        Array4<Real> const& ubar_a = vec_avg_ubar[lev]->array(mfi);
        Array4<Real> const& vbar_a = vec_avg_vbar[lev]->array(mfi);
        Array4<Real> const& xvel_a = vec_avg_u[lev]->array(mfi);
        Array4<Real> const& yvel_a = vec_avg_v[lev]->array(mfi);
        Array4<Real> const& cons_a = vec_avg_cons[lev]->array(mfi);

        Array4<Real const> const& mskr = vec_mskr[lev]->array(mfi);
        Array4<Real const> const& msku = vec_msku[lev]->array(mfi);
        Array4<Real const> const& mskv = vec_mskv[lev]->array(mfi);

        ParallelFor(makeSlab(gbx1,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            if (mskr(i,j,0) == 0.0_rt) {
                zeta_a(i,j,0) = fill_value;
            }
        });
        ParallelFor(gbx1, n_avg_cons, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            if (mskr(i,j,0) == 0.0_rt) {
                cons_a(i,j,k,n) = fill_value;
            }
        });
        ParallelFor(makeSlab(ubx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            if (msku(i,j,0) == 0.0_rt && ubar_a(i,j,0) == fill_where) {
                ubar_a(i,j,0) = fill_value;
            }
        });
        ParallelFor(makeSlab(vbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            if (mskv(i,j,0) == 0.0_rt && vbar_a(i,j,0) == fill_where) {
                vbar_a(i,j,0) = fill_value;
            }
        });
        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (msku(i,j,0) == 0.0_rt && xvel_a(i,j,k) == fill_where) {
                xvel_a(i,j,k) = fill_value;
            }
        });
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (mskv(i,j,0) == 0.0_rt && yvel_a(i,j,k) == fill_where) {
                yvel_a(i,j,k) = fill_value;
            }
        });
    } // mfi
    Gpu::streamSynchronize();
}
