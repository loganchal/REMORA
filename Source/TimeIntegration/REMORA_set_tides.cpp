/**
 * \file REMORA_set_tides.cpp
 *
 * Harmonic evaluation of barotropic tidal forcing, ported from
 * ROMS/Nonlinear/set_tides.F (SSH_TIDES + UV_TIDES, as compiled for the
 * Moana hindcast: RAMP_TIDES undefined, ADD_FSOBC and ADD_M2OBC defined).
 *
 * The oracle used for this port is the preprocessed Moana build of that file
 * (external/moana_spec/set_tides.f90). Line references below are to that file.
 *
 * ROMS index convention vs REMORA: ROMS rho index i (0..Lm+1) corresponds to
 * REMORA cell index i-1 (-1..Nx), and ROMS u index i (1..Lm+1) corresponds to
 * REMORA x-face index i-1 (0..Nx). Because both shift by the same amount, every
 * stencil below keeps the algebraic form it has in ROMS.
 *
 * NOTE on the CLIMA blocks in the oracle (lines 679-691, 824-844): those add the
 * tide to CLIMA(ng)%ssh / ubarclm / vbarclm for 2D climatology nudging. The Moana
 * configuration has LsshCLM/Lm2CLM false (no M2/ssh climatology nudging), so those
 * blocks are intentionally not ported.
 */

#include <REMORA.H>

using namespace amrex;

namespace {
/**
 * Tidal current ellipse at a rho point for one constituent
 * (ROMS set_tides.f90:787-796).
 *
 *   angle = UV_Tangle - angler ; phase = omega - UV_Tphase
 *   Uwrk  = Tmajor*cos(angle)*cos(phase) - Tminor*sin(angle)*sin(phase)
 *   Vwrk  = Tmajor*sin(angle)*cos(phase) + Tminor*cos(angle)*sin(phase)
 *
 * A free function rather than a device lambda so that it is safe to call from
 * inside a GPU kernel under any compiler (nvcc restricts nested extended lambdas).
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void tide_ellipse (int i, int j, int itide, amrex::Real omega,
                   const amrex::Array4<const amrex::Real>& Cmax,
                   const amrex::Array4<const amrex::Real>& Cmin,
                   const amrex::Array4<const amrex::Real>& Cangle,
                   const amrex::Array4<const amrex::Real>& Cphase,
                   const amrex::Array4<const amrex::Real>& angler,
                   amrex::Real& Uwrk, amrex::Real& Vwrk)
{
    amrex::Real angle  = Cangle(i,j,0,itide) - angler(i,j,0);
    amrex::Real Cang   = std::cos(angle);
    amrex::Real Sang   = std::sin(angle);
    amrex::Real phase  = omega - Cphase(i,j,0,itide);
    amrex::Real Cpha   = std::cos(phase);
    amrex::Real Spha   = std::sin(phase);
    Uwrk = Cmax(i,j,0,itide) * Cang * Cpha - Cmin(i,j,0,itide) * Sang * Spha;
    Vwrk = Cmax(i,j,0,itide) * Sang * Cpha + Cmin(i,j,0,itide) * Cang * Spha;
}
} // anonymous namespace

/**
 * Evaluate the tidal elevation and tidal currents for the current baroclinic step.
 *
 * ROMS calls set_tides once per baroclinic step (main3d.F line 624), after
 * set_data has interpolated the open-boundary data to time(ng) and before the
 * barotropic (step2d) loop. REMORA's analogue is the top of REMORA::Advance, with
 * the same time argument (t_old[lev]) that the boundary data is interpolated to.
 *
 * @param[in] lev   level (tides are only supported at level 0)
 * @param[in] time  model time [s] on the ocean_time base; ROMS time(ng)
 */
void
REMORA::set_tides (int lev, Real time)
{
#ifdef REMORA_USE_NETCDF
    BL_PROFILE("REMORA::set_tides()");

    if (!solverChoice.use_tides) return;
    if (lev != 0) return;

    AMREX_ALWAYS_ASSERT(tide_data_from_file != nullptr);

    const int ntide = tide_data_from_file->n_constituents();
    const Real* Tperiod = tide_data_from_file->period_ptr();

    // ROMS set_tides.f90:650 -- RAMP_TIDES is undefined in the Moana build, so ramp = 1
    const Real ramp = one;

    // ROMS set_tides.f90:658, 781:  cff = 2*pi*(time(ng) - tide_start*day2sec)
    const Real cff_time = two * PI * (time - solverChoice.tide_start * Real(86400.0));

    const MultiFab& mf_Eamp   = tide_data_from_file->Eamp();
    const MultiFab& mf_Ephase = tide_data_from_file->Ephase();
    const MultiFab& mf_Cangle = tide_data_from_file->Cangle();
    const MultiFab& mf_Cphase = tide_data_from_file->Cphase();
    const MultiFab& mf_Cmax   = tide_data_from_file->Cmax();
    const MultiFab& mf_Cmin   = tide_data_from_file->Cmin();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_Etide[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // rho points, one ghost ring: ROMS IstrR-1..IendR, JstrR-1..JendR coverage
        Box bx_r = mfi.growntilebox(IntVect(1,1,0));
        // u faces of this grid, one ghost row in y: ROMS Istr..IendR, JstrR..JendR
        Box bx_u = mfi.grownnodaltilebox(0, IntVect(0,1,0));
        // v faces of this grid, one ghost column in x: ROMS IstrR..IendR, Jstr..JendR
        Box bx_v = mfi.grownnodaltilebox(1, IntVect(1,0,0));

        const Array4<Real      >& Etide  = vec_Etide[lev]->array(mfi);
        const Array4<Real      >& Utide  = vec_Utide[lev]->array(mfi);
        const Array4<Real      >& Vtide  = vec_Vtide[lev]->array(mfi);

        const Array4<const Real>& Eamp   = mf_Eamp.const_array(mfi);
        const Array4<const Real>& Ephase = mf_Ephase.const_array(mfi);
        const Array4<const Real>& Cangle = mf_Cangle.const_array(mfi);
        const Array4<const Real>& Cphase = mf_Cphase.const_array(mfi);
        const Array4<const Real>& Cmax   = mf_Cmax.const_array(mfi);
        const Array4<const Real>& Cmin   = mf_Cmin.const_array(mfi);

        const Array4<const Real>& angler = vec_angler[lev]->const_array(mfi);
        const Array4<const Real>& mskr   = vec_mskr[lev]->const_array(mfi);
        const Array4<const Real>& msku   = vec_msku[lev]->const_array(mfi);
        const Array4<const Real>& mskv   = vec_mskv[lev]->const_array(mfi);

        //
        // Tidal elevation (ROMS set_tides.f90:657-673)
        //
        //   Etide(i,j) = sum_c ramp*SSH_Tamp(i,j,c)*COS(omega_c - SSH_Tphase(i,j,c))
        //   masked by rmask inside the constituent loop
        //
        ParallelFor(makeSlab(bx_r,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real etide = zero;
            for (int itide = 0; itide < ntide; itide++) {
                if (Tperiod[itide] > zero) {
                    Real omega = cff_time / Tperiod[itide];
                    etide = etide + ramp * Eamp(i,j,0,itide) *
                            std::cos(omega - Ephase(i,j,0,itide));
                    etide = etide * mskr(i,j,0);
                }
            }
            Etide(i,j,0) = etide;
        });

        //
        // Tidal currents (ROMS set_tides.f90:779-818)
        //
        // Uwrk/Vwrk are the ellipse velocities at rho points:
        //   angle  = UV_Tangle(i,j,c) - angler(i,j)
        //   phase  = omega_c - UV_Tphase(i,j,c)
        //   Uwrk   = Tmajor*cos(angle)*cos(phase) - Tminor*sin(angle)*sin(phase)
        //   Vwrk   = Tmajor*sin(angle)*cos(phase) + Tminor*cos(angle)*sin(phase)
        // They are then averaged to u/v points and masked. Rather than storing the
        // rho-point work arrays, each u (v) point recomputes the two rho-point values
        // it needs; this is algebraically identical and avoids a scratch MultiFab and
        // any inter-thread dependence on the GPU.
        //
        ParallelFor(makeSlab(bx_u,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real utide = zero;
            for (int itide = 0; itide < ntide; itide++) {
                if (Tperiod[itide] > zero) {
                    Real omega = cff_time / Tperiod[itide];
                    Real Uwrk_im1, Vwrk_dummy, Uwrk_i;
                    tide_ellipse(i-1,j,itide,omega,Cmax,Cmin,Cangle,Cphase,angler,Uwrk_im1,Vwrk_dummy);
                    tide_ellipse(i  ,j,itide,omega,Cmax,Cmin,Cangle,Cphase,angler,Uwrk_i  ,Vwrk_dummy);
                    utide = utide + ramp * Real(0.5) * (Uwrk_im1 + Uwrk_i);
                    utide = utide * msku(i,j,0);
                }
            }
            Utide(i,j,0) = utide;
        });

        ParallelFor(makeSlab(bx_v,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real vtide = zero;
            for (int itide = 0; itide < ntide; itide++) {
                if (Tperiod[itide] > zero) {
                    Real omega = cff_time / Tperiod[itide];
                    Real Vwrk_jm1, Uwrk_dummy, Vwrk_j;
                    tide_ellipse(i,j-1,itide,omega,Cmax,Cmin,Cangle,Cphase,angler,Uwrk_dummy,Vwrk_jm1);
                    tide_ellipse(i,j  ,itide,omega,Cmax,Cmin,Cangle,Cphase,angler,Uwrk_dummy,Vwrk_j  );
                    vtide = (vtide + ramp * Real(0.5) * (Vwrk_jm1 + Vwrk_j));
                    vtide = vtide * mskv(i,j,0);
                }
            }
            Vtide(i,j,0) = vtide;
        });
    }
#else
    amrex::ignore_unused(lev, time);
#endif
}
