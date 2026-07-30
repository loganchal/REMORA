#include <REMORA.H>

using namespace amrex;

#ifdef REMORA_USE_FUNWAVE_FORT
#include <REMORA_funwave_Fortran_Interface.H>
#endif

/**
 * @param[in] lev            level of refinement
 * @param[in] time           simulation time at start of step
 * @param[in] dt_lev         baroclinic time step at level
 * @param[in] iteration      iteration in subcycling, if using
 * @param[in] ncycle         total number of subcycles, if using
 */
 void
REMORA::Advance (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("REMORA::Advance()");

    // ROMS main3d.F calls set_tides once per baroclinic step, on the same
    // time(ng) that set_data used, and holds the result fixed across all
    // NDTFAST barotropic sub-steps.
    //
    // Which absolute time that is, is NOT settled by reading the increment at
    // main3d.f90:522 alone. The order in the step body is
    //   iic++ ; time(ng)+=dt ; set_data ; set_tides ; output ; step2d-loop
    // and `output` writes the state BEFORE it is stepped. Since history record 0
    // carries ocean_time = T0, the loop must be entered with time(ng) = T0-dt,
    // so the step that advances T0 -> T0+dt runs set_tides at T0 = t_old, not
    // at t_new. That is the same conclusion the measured bdy_time_shift scan
    // reached for the open-boundary data (sharp minimum at exactly -dt), and
    // ROMS uses ONE time(ng) for both, so the tide and the boundary data cannot
    // be evaluated a step apart.
    //
    // MEASURED, and it is t_old. A four-point scan of remora.tide_time_shift
    // against the ROMS twin over 40 steps: the step-40 error is exactly LINEAR in
    // the offset and vanishes at -dt. Median |port-ROMS| at 80+ cells from land,
    // minus the non-tidal floor of 1.07e-04:
    //   shift  -50 -> 1.466e-03   (predicted offset 50 s  -> ratio 1)
    //   shift    0 -> 2.934e-03   (predicted 100 s -> 2)  measured 2.001
    //   shift +100 -> 5.910e-03   (predicted 200 s -> 4)  measured 4.031
    // and at -dt the error drops to 1.0698e-04 against 1.0701e-04 for the same
    // run with remora.tides=false -- i.e. with the tide evaluated at t_old the
    // tidal code contributes NOTHING to the divergence, to four figures.
    // Improvement over t_new at 80+ cells: zeta 28.4x, ubar 16.9x, vbar 41.0x,
    // which are precisely the factors the tide ablation had attributed to "the
    // tides". Step 1 is unchanged at every shift, because the tidal boundary
    // error has not yet propagated past the 8-cell edge cut.
    //
    // So the tide is evaluated at the time the step STARTS, exactly like the
    // open-boundary data, and for the same reason: ROMS uses one time(ng) for
    // both, and `output` writing the pre-step state forces that value to be T0.
    // The knob stays as a diagnostic and now defaults to no additional shift.
    set_tides(lev, time + solverChoice.tide_time_shift);

    setup_step(lev, time, dt_lev);

    if (solverChoice.use_barotropic)
    {
        int nfast_counter=nfast + 1;

        //***************************************************
        //Compute fast timestep from dt_lev and ratio
        //***************************************************
        Real dtfast_lev=dt_lev/Real(fixed_ndtfast_ratio);

        //***************************************************
        //Advance nfast_counter steps of the 2d integrator
        //***************************************************
        for (int my_iif = 0; my_iif < nfast_counter; my_iif++) {
            advance_2d_onestep(lev, dt_lev, dtfast_lev, my_iif, nfast_counter);
        }
    }

#ifdef REMORA_USE_FUNWAVE_FORT
    MultiFab* mf_rhoS = vec_rhoS[lev].get();
    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box bx = mfi.validbox();
        int ims = bx.smallEnd(0);
        int jms = bx.smallEnd(1);
        int kms = bx.smallEnd(2);
        int ime = bx.bigEnd(0);
        int jme = bx.bigEnd(1);
        int kme = bx.bigEnd(2);

        Array4<Real> const& rho_salt = mf_rhoS->array(mfi);

        funwave_advance_c(rho_salt.dataPtr(), ims, ime, jms, jme, kms, kme);
    }
#endif

    //***************************************************
    //Advance one step of the 3d integrator
    //***************************************************
    advance_3d_ml(lev, dt_lev);

}
