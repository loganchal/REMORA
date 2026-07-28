#include <REMORA.H>

using namespace amrex;

//
// Mellor-Yamada level 2.5 vertical turbulence closure.
//
// Port of ROMS 3.9 ROMS/Nonlinear/my25_prestep.F and my25_corstep.F as
// compiled for the Moana hindcast, i.e. with
//
//     MY25_MIXING, KANTHA_CLAYSON, N2S2_HORAVG, RI_SPLINES, MASKING
//
// active and
//
//     K_C2ADVECTION, K_C4ADVECTION, LIMIT_VDIFF, LIMIT_VVISC,
//     PERFECT_RESTART
//
// inactive.  Line references in the comments below are to the annotated
// originals in external/roms-3.9/ROMS/Nonlinear/; the preprocessed Moana
// build (all #ifdefs resolved) lives in external/moana_spec/.
//
// Index mapping between ROMS and REMORA (see Tests/MoanaParity/my25_port_notes.md):
//
//   N_ROMS   = N + 1                    (N is REMORA's top rho index)
//   w-points : k_REMORA = k_ROMS        (0 = bottom, N+1 = surface)
//   rho cells: k_REMORA = k_ROMS - 1    (so Hz_ROMS(k) == Hz(k-1))
//
// MY2.5 reuses the GLS storage slots: "tke" holds q^2 (twice the turbulent
// kinetic energy) and "gls" holds q^2*l.  Component 2 of each MultiFab is
// the ROMS time level "3" scratch slot written by the prestep.
//

/**
 * Predictor step for the turbulent kinetic energy variables (ROMS my25_prestep.F).
 *
 * A non-conservative but constancy-preserving auxiliary advective substep for
 * the tke and gls equations.  The result is stored in component 2 and is used
 * to build the advective terms in the corrector.  No dissipation here.
 *
 * NOTE: the algebra of ROMS my25_prestep.F is byte-identical to gls_prestep.F
 * (verified by diff of lines 180-424 of both files); MY2.5 differs from GLS
 * only in the corrector.
 *
 * @param[in   ] lev            level to operate on
 * @param[inout] mf_gls         q^2*l
 * @param[inout] mf_tke         q^2
 * @param[in   ] mf_W           vertical velocity
 * @param[in   ] mf_msku        land-sea mask on u points
 * @param[in   ] mf_mskv        land-sea mask on v points
 * @param[in   ] nstp           index of last time step in gls and tke MultiFabs
 * @param[in   ] nnew           index of time step to update in gls and tke MultiFabs
 * @param[in   ] iic            which time step we're on
 * @param[in   ] ntfirst        what is the first time step?
 * @param[in   ] N              top rho index
 * @param[in   ] dt_lev         time step at this level
 */
void
REMORA::my25_prestep (int lev, MultiFab* mf_gls, MultiFab* mf_tke,
                      MultiFab& mf_W, MultiFab* mf_msku, MultiFab* mf_mskv,
                      const int nstp, const int nnew,
                      const int iic, const int ntfirst, const int N, const Real dt_lev)
{
    BL_PROFILE("REMORA::my25_prestep()");

    for ( MFIter mfi(*mf_gls, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& gls = mf_gls->array(mfi);
        Array4<Real> const& tke = mf_tke->array(mfi);
        Array4<Real const> const& W = mf_W.const_array(mfi);

        Array4<Real const> const& Huon = vec_Huon[lev]->const_array(mfi);
        Array4<Real const> const& Hvom = vec_Hvom[lev]->const_array(mfi);
        Array4<Real const> const& Hz = vec_Hz[lev]->const_array(mfi);
        Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);
        Array4<Real const> const& msku = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv = mf_mskv->const_array(mfi);

        Box bx = mfi.tilebox();          // z-node box, k = 0..N+1
        Box xbx = surroundingNodes(bx,0);
        Box ybx = surroundingNodes(bx,1);

        Box xbx_hi = growHi(xbx,0,1);
        Box ybx_hi = growHi(ybx,0,1);

        const Box& domain = geom[lev].Domain();
        const auto dlo = amrex::lbound(domain);
        const auto dhi = amrex::ubound(domain);

        GeometryData const& geomdata = geom[0].data();
        bool is_periodic_in_x = geomdata.isPeriodic(0);
        bool is_periodic_in_y = geomdata.isPeriodic(1);

        FArrayBox fab_XF(xbx_hi, 1, amrex::The_Async_Arena()); fab_XF.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FX(xbx_hi, 1, amrex::The_Async_Arena()); fab_FX.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FXL(xbx_hi, 1, amrex::The_Async_Arena()); fab_FXL.template setVal<RunOn::Device>(zero);
        FArrayBox fab_EF(ybx_hi, 1, amrex::The_Async_Arena()); fab_EF.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FE(ybx_hi, 1, amrex::The_Async_Arena()); fab_FE.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FEL(ybx_hi, 1, amrex::The_Async_Arena()); fab_FEL.template setVal<RunOn::Device>(zero);
        FArrayBox fab_Hz_half(bx, 1, amrex::The_Async_Arena()); fab_Hz_half.template setVal<RunOn::Device>(zero);
        FArrayBox fab_CF(convert(bx,IntVect(0,0,0)), 1, amrex::The_Async_Arena()); fab_CF.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FC(convert(bx,IntVect(0,0,0)), 1, amrex::The_Async_Arena()); fab_FC.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FCL(convert(bx,IntVect(0,0,0)), 1, amrex::The_Async_Arena()); fab_FCL.template setVal<RunOn::Device>(zero);

        auto XF  = fab_XF.array();
        auto FX  = fab_FX.array();
        auto FXL = fab_FXL.array();
        auto EF  = fab_EF.array();
        auto FE  = fab_FE.array();
        auto FEL = fab_FEL.array();
        auto Hz_half = fab_Hz_half.array();
        auto CF  = fab_CF.array();
        auto FC  = fab_FC.array();
        auto FCL = fab_FCL.array();

        // ------------------------------------------------------------------
        // my25_prestep.F:208-252 -- fourth-order centered differences, xi flux.
        // grad(i) := tke(i)-tke(i-1) on u points, masked (:213-219); the
        // gradient is copied inward one point at a non-periodic west/east
        // boundary (:222-237), then FX/FXL are formed on u points (:239-252).
        // ------------------------------------------------------------------
        ParallelFor(grow(xbx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real grad_im1 = (tke(i-1,j,k,nstp) - tke(i-2,j,k,nstp)) * msku(i-1,j,0);
            Real grad_ip1 = (tke(i+1,j,k,nstp) - tke(i  ,j,k,nstp)) * msku(i+1,j,0);

            Real gradL_im1 = (gls(i-1,j,k,nstp) - gls(i-2,j,k,nstp)) * msku(i-1,j,0);
            Real gradL_ip1 = (gls(i+1,j,k,nstp) - gls(i  ,j,k,nstp)) * msku(i+1,j,0);

            // my25_prestep.F:222-229 grad(Istr-1,j)=grad(Istr,j).  The u face
            // ROMS calls Istr is REMORA's dlo.x, so the one-sided copy is felt
            // by the flux at i == dlo.x through its grad(i-1) stencil arm.
            if (i == dlo.x && !is_periodic_in_x) {
                grad_im1  = (tke(i,j,k,nstp) - tke(i-1,j,k,nstp)) * msku(i,j,0);
                gradL_im1 = (gls(i,j,k,nstp) - gls(i-1,j,k,nstp)) * msku(i,j,0);
            }
            // my25_prestep.F:230-237 grad(Iend+2,j)=grad(Iend+1,j); ROMS Iend+1
            // is REMORA's dhi.x+1.
            else if (i == dhi.x+1 && !is_periodic_in_x) {
                grad_ip1  = (tke(i,j,k,nstp) - tke(i-1,j,k,nstp)) * msku(i,j,0);
                gradL_ip1 = (gls(i,j,k,nstp) - gls(i-1,j,k,nstp)) * msku(i,j,0);
            }
            Real cff = one/Real(6.0);
            XF(i,j,k) = Real(0.5) * (Huon(i,j,k) + Huon(i,j,k-1));
            FX(i,j,k) = XF(i,j,k) * Real(0.5) * (tke(i-1,j,k,nstp) + tke(i,j,k,nstp) -
                cff * (grad_ip1 - grad_im1));
            FXL(i,j,k) = XF(i,j,k) * Real(0.5) * (gls(i-1,j,k,nstp) + gls(i,j,k,nstp) -
                cff * (gradL_ip1 - gradL_im1));
        });

        // ------------------------------------------------------------------
        // my25_prestep.F:253-291 -- same construction in the eta direction.
        // ------------------------------------------------------------------
        ParallelFor(grow(ybx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real grad_jm1 = (tke(i,j-1,k,nstp) - tke(i,j-2,k,nstp)) * mskv(i,j-1,0);
            Real grad_jp1 = (tke(i,j+1,k,nstp) - tke(i,j  ,k,nstp)) * mskv(i,j+1,0);

            Real gradL_jm1 = (gls(i,j-1,k,nstp) - gls(i,j-2,k,nstp)) * mskv(i,j-1,0);
            Real gradL_jp1 = (gls(i,j+1,k,nstp) - gls(i,j  ,k,nstp)) * mskv(i,j+1,0);

            // my25_prestep.F:263-270 grad(i,Jstr-1)=grad(i,Jstr)
            if (j == dlo.y && !is_periodic_in_y) {
                grad_jm1  = (tke(i,j,k,nstp) - tke(i,j-1,k,nstp)) * mskv(i,j,0);
                gradL_jm1 = (gls(i,j,k,nstp) - gls(i,j-1,k,nstp)) * mskv(i,j,0);
            }
            // my25_prestep.F:271-278 grad(i,Jend+2)=grad(i,Jend+1)
            else if (j == dhi.y+1 && !is_periodic_in_y) {
                grad_jp1  = (tke(i,j,k,nstp) - tke(i,j-1,k,nstp)) * mskv(i,j,0);
                gradL_jp1 = (gls(i,j,k,nstp) - gls(i,j-1,k,nstp)) * mskv(i,j,0);
            }
            Real cff = one/Real(6.0);
            EF(i,j,k) = Real(0.5) * (Hvom(i,j,k) + Hvom(i,j,k-1));
            FE(i,j,k) = EF(i,j,k) * Real(0.5) * (tke(i,j-1,k,nstp) + tke(i,j,k,nstp) -
                cff * (grad_jp1 - grad_jm1));
            FEL(i,j,k) = EF(i,j,k) * Real(0.5) * (gls(i,j-1,k,nstp) + gls(i,j,k,nstp) -
                cff * (gradL_jp1 - gradL_jm1));
        });

        // ------------------------------------------------------------------
        // my25_prestep.F:293-325 -- LF step plus AM3 interpolation back a half
        // step, fused into the weights cff1/cff2/cff3.  Also stashes
        // tke,gls(:,:,:,nnew) = Hz * tke,gls(:,:,:,nstp) while the old-time Hz
        // is still available.
        // ------------------------------------------------------------------
        Real gamma = one / Real(6.0);             // my25_prestep.F:139
        Real cff1, cff2, cff3;
        int indx;
        if (iic == ntfirst) {                     // my25_prestep.F:295-299
            cff1 = one;
            cff2 = zero;
            cff3 = Real(0.5) * dt_lev;
            indx = nstp;
        } else {                                  // my25_prestep.F:300-305
            cff1 = Real(0.5) + gamma;
            cff2 = Real(0.5) - gamma;
            cff3 = (one - gamma) * dt_lev;
            indx = 1 - nstp;
        }

        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = Real(0.5) * (Hz(i,j,k) + Hz(i,j,k-1));
            Real cff4 = cff3 * pm(i,j,0) * pn(i,j,0);
            Hz_half(i,j,k) = cff - cff4 * (XF(i+1,j,k)-XF(i,j,k)+EF(i,j+1,k)-EF(i,j,k));
            tke(i,j,k,2) = cff * (cff1*tke(i,j,k,nstp) + cff2*tke(i,j,k,indx)) -
                           cff4 * (FX(i+1,j,k)-FX(i,j,k)+FE(i,j+1,k)-FE(i,j,k));
            gls(i,j,k,2) = cff * (cff1 * gls(i,j,k,nstp) + cff2 * gls(i,j,k,indx)) -
                           cff4 * (FXL(i+1,j,k)-FXL(i,j,k)+FEL(i,j+1,k)-FEL(i,j,k));
            tke(i,j,k,nnew) = cff * tke(i,j,k,nstp);
            gls(i,j,k,nnew) = cff * gls(i,j,k,nstp);
        });

        // ------------------------------------------------------------------
        // my25_prestep.F:327-374 -- vertical advective flux at rho points.
        // Fourth-order interior (:339-354), one-sided at the two end cells
        // (:355-373).  Column kernel: no running k dependence, but keeping the
        // whole column in one thread avoids a separate boundary launch.
        // ------------------------------------------------------------------
        Box bxD = bx; bxD.makeSlab(2,0);
        ParallelFor(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            const Real c7o12 = Real(7.0)/Real(12.0);
            const Real c1o12 = one/Real(12.0);
            const Real c1o3  = one/Real(3.0);
            const Real c5o6  = Real(5.0)/Real(6.0);
            const Real c1o6  = one/Real(6.0);

            for (int k=1; k<=N-1; k++) {
                CF(i,j,k) = Real(0.5) * (W(i,j,k+1) + W(i,j,k));
                FC(i,j,k)  = CF(i,j,k) * (c7o12 * (tke(i,j,k  ,nstp) + tke(i,j,k+1,nstp)) -
                                          c1o12 * (tke(i,j,k-1,nstp) + tke(i,j,k+2,nstp)));
                FCL(i,j,k) = CF(i,j,k) * (c7o12 * (gls(i,j,k  ,nstp) + gls(i,j,k+1,nstp)) -
                                          c1o12 * (gls(i,j,k-1,nstp) + gls(i,j,k+2,nstp)));
            }
            CF(i,j,0) = Real(0.5) * (W(i,j,1) + W(i,j,0));
            FC(i,j,0)  = CF(i,j,0) * (c1o3*tke(i,j,0,nstp) + c5o6*tke(i,j,1,nstp) - c1o6*tke(i,j,2,nstp));
            FCL(i,j,0) = CF(i,j,0) * (c1o3*gls(i,j,0,nstp) + c5o6*gls(i,j,1,nstp) - c1o6*gls(i,j,2,nstp));

            CF(i,j,N) = Real(0.5) * (W(i,j,N+1) + W(i,j,N));
            FC(i,j,N)  = CF(i,j,N) * (c1o3*tke(i,j,N+1,nstp) + c5o6*tke(i,j,N,nstp) - c1o6*tke(i,j,N-1,nstp));
            FCL(i,j,N) = CF(i,j,N) * (c1o3*gls(i,j,N+1,nstp) + c5o6*gls(i,j,N,nstp) - c1o6*gls(i,j,N-1,nstp));
        });

        // ------------------------------------------------------------------
        // my25_prestep.F:376-393 -- time-step the vertical advective term and
        // divide through by the auxiliary layer thickness Hz_half.
        // ------------------------------------------------------------------
        Real cff3_v = (iic == ntfirst) ? Real(0.5) * dt_lev : (one - gamma) * dt_lev;
        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff4 = cff3_v * pm(i,j,0) * pn(i,j,0);
            Hz_half(i,j,k) = Hz_half(i,j,k) - cff4 * (CF(i,j,k)-CF(i,j,k-1));
            Real cff1_loc = one / Hz_half(i,j,k);
            tke(i,j,k,2) = cff1_loc * (tke(i,j,k,2) - cff4 * (FC (i,j,k) - FC (i,j,k-1)));
            gls(i,j,k,2) = cff1_loc * (gls(i,j,k,2) - cff4 * (FCL(i,j,k) - FCL(i,j,k-1)));
        });
    }

    // my25_prestep.F:396-411 -- tkebc_tile + exchange_w3d_tile.
    for (int icomp=0; icomp<3; icomp++) {
        FillPatch(lev, t_old[lev], *vec_tke[lev], GetVecOfPtrs(vec_tke), tke_bc(), BdyVars::null, icomp, false, false);
        FillPatch(lev, t_old[lev], *vec_gls[lev], GetVecOfPtrs(vec_gls), tke_bc(), BdyVars::null, icomp, false, false);
    }
}

/**
 * Corrector step for the Mellor-Yamada 2.5 closure (ROMS my25_corstep.F).
 *
 * @param[in   ] lev            level to operate on
 * @param[inout] mf_gls         q^2*l
 * @param[inout] mf_tke         q^2
 * @param[in   ] mf_W           vertical velocity
 * @param[inout] mf_Akv         vertical viscosity coefficient
 * @param[inout] mf_Akt         vertical diffusivity coefficients
 * @param[inout] mf_Akk         turbulent kinetic energy vertical diffusion coefficient
 * @param[in   ] mf_mskr        land-sea mask on rho points
 * @param[in   ] mf_msku        land-sea mask on u points
 * @param[in   ] mf_mskv        land-sea mask on v points
 * @param[in   ] nstp           index of last time step in gls and tke MultiFabs
 * @param[in   ] nnew           index of time step to update in gls and tke MultiFabs
 * @param[in   ] N              top rho index
 * @param[in   ] dt_lev         time step at this level
 */
void
REMORA::my25_corrector (int lev, MultiFab* mf_gls, MultiFab* mf_tke,
                        MultiFab& mf_W, MultiFab* mf_Akv, MultiFab* mf_Akt,
                        MultiFab* mf_Akk,
                        MultiFab* mf_mskr,
                        MultiFab* mf_msku, MultiFab* mf_mskv,
                        const int nstp, const int nnew,
                        const int N, const Real dt_lev)
{
    BL_PROFILE("REMORA::my25_corrector()");

    // ------------------------------------------------------------------
    // Closure constants.  Base values are mod_scalars.F:1754-1767; the
    // derived combinations are mod_scalars.F:3761-3775 with KANTHA_CLAYSON
    // defined (so my_Sm3 is *not* used and my_B1pm1o3 takes its place in Sm).
    // ------------------------------------------------------------------
    const Real my_A1 = solverChoice.my_A1;
    const Real my_A2 = solverChoice.my_A2;
    const Real my_B1 = solverChoice.my_B1;
    const Real my_B2 = solverChoice.my_B2;
    const Real my_C1 = solverChoice.my_C1;
    const Real my_C2 = solverChoice.my_C2;
    const Real my_C3 = solverChoice.my_C3;
    const Real my_E1 = solverChoice.my_E1;
    const Real my_E2 = solverChoice.my_E2;
    const Real my_Gh0  = solverChoice.my_Gh0;
    const Real my_Sq   = solverChoice.my_Sq;
    const Real my_lmax = solverChoice.my_lmax;
    const Real my_qmin = solverChoice.my_qmin;

    const Real my_B1p2o3  = std::pow(my_B1, Real(2.0)/Real(3.0));                 // mod_scalars.F:3761
    const Real my_B1pm1o3 = one / std::pow(my_B1, one/Real(3.0));                 // mod_scalars.F:3762
    const Real my_Sm2 = Real(9.0)*my_A1*my_A2;                                    // mod_scalars.F:3767
    const Real my_Sh1 = my_A2*(one - Real(6.0)*my_A1/my_B1);                      // mod_scalars.F:3768
    const Real my_Sh2 = Real(3.0)*my_A2*(Real(6.0)*my_A1 + my_B2*(one - my_C3));  // mod_scalars.F:3770 (KANTHA_CLAYSON)
    const Real my_Sm4 = Real(18.0)*my_A1*my_A1
                      + Real(9.0)*my_A1*my_A2*(one - my_C2);                      // mod_scalars.F:3771 (KANTHA_CLAYSON)
    // Wall-proximity prefactor, my25_corstep.F:598
    const Real wall_cff = my_E2 / (vonKar*vonKar);

    const Real Akv_bak = solverChoice.Akv_bak;
    const Real Akt_bak = solverChoice.Akt_bak;
    const Real Akk_bak = solverChoice.Akk_bak;

    const Real Gadv = one/Real(3.0);        // my25_corstep.F:189
    const Real eps  = Real(1.0e-10);        // my25_corstep.F:190

    const BoxArray&            ba = cons_old[lev]->boxArray();
    const DistributionMapping& dm = cons_old[lev]->DistributionMap();

    int ncomp_w = 0;
    int dU_comp = ncomp_w++;
    int dV_comp = ncomp_w++;
    int CF_comp = ncomp_w++;

    int ncomp = 0;
    int shear2_comp = ncomp++;
    int shear2_cache_comp = ncomp++;
    int buoy2_comp = ncomp++;

    MultiFab mf_w(convert(ba, IntVect(0,0,1)),dm,ncomp_w,IntVect(NGROW,NGROW,0));
    MultiFab mf(ba,dm,ncomp,IntVect(NGROW,NGROW,0));

    const Box& domain = geom[0].Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    GeometryData const& geomdata = geom[0].data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);

    // ==================================================================
    // my25_corstep.F:216-255 (RI_SPLINES branch) -- vertical velocity shear
    // at W points from a parabolic spline reconstruction of du/dz, dv/dz.
    // The forward elimination and back substitution are serial in k, so this
    // is a strict column-per-thread kernel.
    // ==================================================================
    for ( MFIter mfi(*mf_gls, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box   bx = mfi.tilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));

        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);

        Array4<Real> const& Hz = vec_Hz[lev]->array(mfi);
        Array4<Real> const& u = xvel_old[lev]->array(mfi);
        Array4<Real> const& v = yvel_old[lev]->array(mfi);

        auto dU = mf_w.array(mfi,dU_comp);
        auto dV = mf_w.array(mfi,dV_comp);
        auto CF = mf_w.array(mfi,CF_comp);
        auto shear2_cached = mf.array(mfi,shear2_cache_comp);

        ParallelFor(gbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            CF(i,j,0) = zero;                                   // my25_corstep.F:222-224
            dU(i,j,0) = zero;
            dV(i,j,0) = zero;
            for (int k=1; k<=N; k++) {                          // my25_corstep.F:226-238
                Real cff = one / (two * Hz(i,j,k) + Hz(i,j,k-1)*(two - CF(i,j,k-1)));
                CF(i,j,k) = cff * Hz(i,j,k);
                dU(i,j,k)=cff*(Real(3.0)*(u(i  ,j,k)-u(i,  j,k-1)+
                                          u(i+1,j,k)-u(i+1,j,k-1))-Hz(i,j,k-1)*dU(i,j,k-1));
                dV(i,j,k)=cff*(Real(3.0)*(v(i,j  ,k)-v(i,j  ,k-1)+
                                          v(i,j+1,k)-v(i,j+1,k-1))-Hz(i,j,k-1)*dV(i,j,k-1));
            }
            dU(i,j,N+1) = zero;                                 // my25_corstep.F:239-242
            dV(i,j,N+1) = zero;
            for (int k=N; k>=1; k--) {                          // my25_corstep.F:243-248
                dU(i,j,k) = dU(i,j,k) - CF(i,j,k) * dU(i,j,k+1);
                dV(i,j,k) = dV(i,j,k) - CF(i,j,k) * dV(i,j,k+1);
            }
            shear2_cached(i,j,0) = zero;
            for (int k=1; k<=N; k++) {                          // my25_corstep.F:249-253
                shear2_cached(i,j,k) = dU(i,j,k) * dU(i,j,k) + dV(i,j,k) * dV(i,j,k);
            }
        });
    }

    // my25_corstep.F:285-318 -- one-sided copy of shear2 onto the edge and
    // corner rows before the N2S2_HORAVG smoother.  ROMS does this at every
    // boundary, periodic included, hence foextrap here as in the GLS path.
    (*physbcs[lev])(mf,*mf_mskr,shear2_cache_comp,1,mf.nGrowVect(),t_new[lev],foextrap_bc());
    mf.setVal(zero,CF_comp,1);

    int ncomp_fab = 0;
    int tmp_buoy_comp  = ncomp_fab++;
    int tmp_shear_comp = ncomp_fab++;
    int curvK_comp = ncomp_fab++;
    int curvP_comp = ncomp_fab++;
    int FXK_comp = ncomp_fab++;
    int FXP_comp = ncomp_fab++;
    int FEK_comp = ncomp_fab++;
    int FEP_comp = ncomp_fab++;
    int FCK_comp = ncomp_fab++;
    int FCP_comp = ncomp_fab++;
    int BCK_comp = ncomp_fab++;
    int BCP_comp = ncomp_fab++;

    const int ncons_local = ncons;

    for ( MFIter mfi(*mf_gls, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box  bx = mfi.tilebox();                 // z-node box, k = 0..N+1
        Box xbx = surroundingNodes(bx,0);
        Box ybx = surroundingNodes(bx,1);
        Box gbx1 = grow(bx,IntVect(NGROW-1,NGROW-1,0));

        Box bx_growloxy = growLo(growLo(grow(bx,IntVect(0,0,-1)),0,1),1,1);

        Box bxD = bx;
        bxD.makeSlab(2,0);

        Array4<Real const> const& W = mf_W.const_array(mfi);
        Array4<Real> const& Hz = vec_Hz[lev]->array(mfi);
        Array4<Real> const& pm = vec_pm[lev]->array(mfi);
        Array4<Real> const& pn = vec_pn[lev]->array(mfi);
        Array4<Real> const& Lscale = vec_Lscale[lev]->array(mfi);

        Array4<Real> const& Huon = vec_Huon[lev]->array(mfi);
        Array4<Real> const& Hvom = vec_Hvom[lev]->array(mfi);
        Array4<Real> const& z_w = vec_z_w[lev]->array(mfi);

        Array4<Real> const& tke = mf_tke->array(mfi);
        Array4<Real> const& gls = mf_gls->array(mfi);

        Array4<Real const> const& sustr = vec_sustr[lev]->const_array(mfi);
        Array4<Real const> const& svstr = vec_svstr[lev]->const_array(mfi);
        Array4<Real const> const& bustr = vec_bustr[lev]->const_array(mfi);
        Array4<Real const> const& bvstr = vec_bvstr[lev]->const_array(mfi);
        Array4<Real const> const& msku = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv = mf_mskv->const_array(mfi);

        FArrayBox fab(gbx1,ncomp_fab, amrex::The_Async_Arena()); fab.template setVal<RunOn::Device>(zero);

        auto CF = mf_w.array(mfi,CF_comp);
        auto shear2 = mf.array(mfi,shear2_comp);
        auto shear2_cached = mf.array(mfi,shear2_cache_comp);
        auto buoy2 = mf.array(mfi,buoy2_comp);
        Array4<Real> const& bvf = vec_bvf[lev]->array(mfi);

        auto tmp_buoy = fab.array(tmp_buoy_comp);
        auto tmp_shear = fab.array(tmp_shear_comp);
        auto curvK = fab.array(curvK_comp);
        auto curvP = fab.array(curvP_comp);
        auto FXK = fab.array(FXK_comp);
        auto FXP = fab.array(FXP_comp);
        auto FEK = fab.array(FEK_comp);
        auto FEP = fab.array(FEP_comp);
        auto FCK = fab.array(FCK_comp);
        auto FCP = fab.array(FCP_comp);
        auto BCK = fab.array(BCK_comp);
        auto BCP = fab.array(BCP_comp);

        auto Akt = mf_Akt->array(mfi);
        auto Akv = mf_Akv->array(mfi);
        auto Akk = mf_Akk->array(mfi);

        // ==============================================================
        // my25_corstep.F:279-338 (N2S2_HORAVG) -- two passes of a 2x2 box
        // average of N^2 (from bvf, loaded at :270-278) and S^2.  The first
        // pass writes the upper-right corner average into scratch (:319-328),
        // the second averages the four surrounding corners back to rho
        // (:329-337).
        // ==============================================================
        ParallelFor(bx_growloxy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            tmp_buoy(i,j,k)=Real(0.25) * (bvf(i,j,k) + bvf(i+1,j,k) + bvf(i,j+1,k)+bvf(i+1,j+1,k));
            tmp_shear(i,j,k)=Real(0.25) * (shear2_cached(i,j,k) + shear2_cached(i+1,j,k)
                                         + shear2_cached(i,j+1,k) + shear2_cached(i+1,j+1,k));
        });

        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            buoy2(i,j,k)=Real(0.25) * (tmp_buoy(i,j,k) + tmp_buoy(i-1,j,k) + tmp_buoy(i,j-1,k)+tmp_buoy(i-1,j-1,k));
            shear2(i,j,k)=Real(0.25) * (tmp_shear(i,j,k) + tmp_shear(i-1,j,k) + tmp_shear(i,j-1,k)+tmp_shear(i-1,j-1,k));
        });

        // ==============================================================
        // my25_corstep.F:369-397 + 412-438 -- third-order upstream-biased
        // horizontal advection of tke,gls at time level "3" (component 2).
        // gradK/gradP are the masked one-sided differences (:373-380), copied
        // inward at non-periodic boundaries (:381-397), and curvK/curvP is
        // their second difference (:407-411 xi / :456-461 eta).
        // ==============================================================
        ParallelFor(growLo(grow(xbx,IntVect(0,0,-1)),0,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real gradK, gradK_ip1, gradP, gradP_ip1;

            if (i == dlo.x-1 && !is_periodic_in_x) {          // my25_corstep.F:381-388
                gradK_ip1 = (tke(i+1,j,k,2)-tke(i  ,j,k,2)) * msku(i+1,j,0);
                gradK = gradK_ip1;
                gradP_ip1 = (gls(i+1,j,k,2)-gls(i  ,j,k,2)) * msku(i+1,j,0);
                gradP = gradP_ip1;
            } else if (i == dhi.x+1 && !is_periodic_in_x) {   // my25_corstep.F:389-396
                gradK = (tke(i  ,j,k,2)-tke(i-1,j,k,2)) * msku(i,j,0);
                gradK_ip1 = gradK;
                gradP = (gls(i  ,j,k,2)-gls(i-1,j,k,2)) * msku(i,j,0);
                gradP_ip1 = gradP;
            } else {                                          // my25_corstep.F:373-380
                gradK     = (tke(i  ,j,k,2)-tke(i-1,j,k,2)) * msku(i  ,j,0);
                gradK_ip1 = (tke(i+1,j,k,2)-tke(i  ,j,k,2)) * msku(i+1,j,0);
                gradP     = (gls(i  ,j,k,2)-gls(i-1,j,k,2)) * msku(i  ,j,0);
                gradP_ip1 = (gls(i+1,j,k,2)-gls(i  ,j,k,2)) * msku(i+1,j,0);
            }

            curvK(i,j,k) = gradK_ip1 - gradK;
            curvP(i,j,k) = gradP_ip1 - gradP;
        });
        // my25_corstep.F:412-438 -- upstream pick of the curvature and the flux
        ParallelFor(grow(xbx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = Real(0.5) * (Huon(i,j,k) + Huon(i,j,k-1));
            Real cff1 = (cff > zero) ? curvK(i-1,j,k) : curvK(i,j,k);
            Real cff2 = (cff > zero) ? curvP(i-1,j,k) : curvP(i,j,k);

            FXK(i,j,k) = cff * Real(0.5) * (tke(i-1,j,k,2)+tke(i,j,k,2)-Gadv*cff1);
            FXP(i,j,k) = cff * Real(0.5) * (gls(i-1,j,k,2)+gls(i,j,k,2)-Gadv*cff2);
        });

        // my25_corstep.F:441-461 -- eta-direction gradients and curvature
        ParallelFor(growLo(grow(ybx,IntVect(0,0,-1)),1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real gradK     = (tke(i,j  ,k,2)-tke(i,j-1,k,2)) * mskv(i,j  ,0);
            Real gradK_jp1 = (tke(i,j+1,k,2)-tke(i,j  ,k,2)) * mskv(i,j+1,0);
            Real gradP     = (gls(i,j  ,k,2)-gls(i,j-1,k,2)) * mskv(i,j  ,0);
            Real gradP_jp1 = (gls(i,j+1,k,2)-gls(i,j  ,k,2)) * mskv(i,j+1,0);

            if (j == dlo.y-1 && !is_periodic_in_y) {          // my25_corstep.F:450-455 (S)
                gradK = gradK_jp1;
                gradP = gradP_jp1;
            }
            else if (j == dhi.y+1 && !is_periodic_in_y) {     // my25_corstep.F:456-461 (N)
                gradK_jp1 = gradK;
                gradP_jp1 = gradP;
            }

            curvK(i,j,k) = gradK_jp1 - gradK;
            curvP(i,j,k) = gradP_jp1 - gradP;
        });
        // my25_corstep.F:481-502
        ParallelFor(grow(ybx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = Real(0.5) * (Hvom(i,j,k) + Hvom(i,j,k-1));
            Real cff1 = (cff > zero) ? curvK(i,j-1,k) : curvK(i,j,k);
            Real cff2 = (cff > zero) ? curvP(i,j-1,k) : curvP(i,j,k);

            FEK(i,j,k) = cff * Real(0.5) * (tke(i,j-1,k,2)+tke(i,j,k,2)-Gadv*cff1);
            FEP(i,j,k) = cff * Real(0.5) * (gls(i,j-1,k,2)+gls(i,j,k,2)-Gadv*cff2);
        });

        // my25_corstep.F:505-521 -- time-step horizontal advection.  Unlike
        // gls_corstep.F, MY2.5 applies NO floor here; the only clamp in the
        // whole corrector is my_qmin at :706-707.
        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = dt_lev * pm(i,j,0) * pn(i,j,0);
            tke(i,j,k,nnew) = tke(i,j,k,nnew) - cff * (FXK(i+1,j  ,k)-FXK(i,j,k)+
                                                       FEK(i  ,j+1,k)-FEK(i,j,k));
            gls(i,j,k,nnew) = gls(i,j,k,nnew) - cff * (FXP(i+1,j  ,k)-FXP(i,j,k)+
                                                       FEP(i  ,j+1,k)-FEP(i,j,k));
        });

        // ==============================================================
        // my25_corstep.F:523-566 -- vertical advective flux of tke,gls at
        // rho points (fourth-order interior, one-sided end cells).
        // ==============================================================
        ParallelFor(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            const Real c7o12 = Real(7.0)/Real(12.0);
            const Real c1o12 = one/Real(12.0);
            const Real c1o3  = one/Real(3.0);
            const Real c5o6  = Real(5.0)/Real(6.0);
            const Real c1o6  = one/Real(6.0);

            for (int k=1; k<=N-1; k++) {                       // my25_corstep.F:534-545
                Real cff = Real(0.5) * (W(i,j,k+1)+W(i,j,k));
                FCK(i,j,k) = cff * (c7o12 * (tke(i,j,k  ,2)+tke(i,j,k+1,2))-
                                    c1o12 * (tke(i,j,k-1,2)+tke(i,j,k+2,2)));
                FCP(i,j,k) = cff * (c7o12 * (gls(i,j,k  ,2)+gls(i,j,k+1,2))-
                                    c1o12 * (gls(i,j,k-1,2)+gls(i,j,k+2,2)));
            }
            Real cff = Real(0.5) * (W(i,j,0)+W(i,j,1));        // my25_corstep.F:549-556
            FCK(i,j,0) = cff * (c1o3*tke(i,j,0,2)+c5o6*tke(i,j,1,2)-c1o6*tke(i,j,2,2));
            FCP(i,j,0) = cff * (c1o3*gls(i,j,0,2)+c5o6*gls(i,j,1,2)-c1o6*gls(i,j,2,2));
            cff = Real(0.5) * (W(i,j,N+1)+W(i,j,N));           // my25_corstep.F:557-564
            FCK(i,j,N) = cff * (c1o3*tke(i,j,N+1,2)+c5o6*tke(i,j,N,2)-c1o6*tke(i,j,N-1,2));
            FCP(i,j,N) = cff * (c1o3*gls(i,j,N+1,2)+c5o6*gls(i,j,N,2)-c1o6*gls(i,j,N-1,2));
        });

        // my25_corstep.F:568-578 -- time-step vertical advection (no floor).
        ParallelFor(grow(bx,2,-1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = dt_lev * pm(i,j,0) * pn(i,j,0);
            tke(i,j,k,nnew) = tke(i,j,k,nnew) - cff*(FCK(i,j,k  )-FCK(i,j,k-1));
            gls(i,j,k,nnew) = gls(i,j,k,nnew) - cff*(FCP(i,j,k  )-FCP(i,j,k-1));
        });

        // ==============================================================
        // my25_corstep.F:585-593 -- off-diagonal (vertical mixing) term of the
        // tridiagonal system.  MY2.5 uses Akk for BOTH the tke and the q^2*l
        // equation (there is no Akp), and unlike gls_corstep.F it does NOT
        // zero FCK at k=1 and k=N: the surface/bottom conditions are Dirichlet,
        // so those coefficients carry the boundary values into the RHS.
        // FCK here is on rho points: FCK_REMORA(k) == FCK_ROMS(k+1).
        // ==============================================================
        Real cff_diff = -Real(0.5) * dt_lev;
        ParallelFor(convert(bx,IntVect(0,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FCK(i,j,k) = cff_diff * (Akk(i,j,k) + Akk(i,j,k+1)) / Hz(i,j,k);
            CF(i,j,k) = zero;
        });

        // ==============================================================
        // my25_corstep.F:595-636 -- shear and buoyancy production, the
        // Mellor-Yamada dissipation, and the wall-proximity function.
        // ==============================================================
        ParallelFor(grow(bx,2,-1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // my25_corstep.F:601-612: ignore small negative buoyancy
            Real strat2 = ((buoy2(i,j,k) > Real(-5.0e-5)) && (buoy2(i,j,k) < zero))
                          ? zero : buoy2(i,j,k);
            Real Qprod = shear2(i,j,k) * (Akv(i,j,k)-Akv_bak) -
                         strat2 * (Akt(i,j,k,Temp_comp)-Akt_bak);

            // my25_corstep.F:613-617: old-time-step unlimited length scale
            Real Ls_unlmt = std::max(eps, gls(i,j,k,nstp)/std::max(tke(i,j,k,nstp),eps));

            // my25_corstep.F:618-625: time-step the production terms
            Real cff1 = Real(0.5) * (Hz(i,j,k-1) + Hz(i,j,k));
            tke(i,j,k,nnew) = tke(i,j,k,nnew) + dt_lev*cff1*Qprod*two;
            gls(i,j,k,nnew) = gls(i,j,k,nnew) + dt_lev*cff1*Qprod*my_E1*Ls_unlmt;

            // my25_corstep.F:626-636: dissipation + wall proximity, folded into
            // the tridiagonal diagonal.  ROMS z_w(N(ng)) is REMORA z_w(N+1).
            Real Qdiss = dt_lev*std::sqrt(tke(i,j,k,nstp))/(my_B1*Ls_unlmt);
            Real cffw = Ls_unlmt*(one/(z_w(i,j,N+1)-z_w(i,j,k)) +
                                  one/(z_w(i,j,k)  -z_w(i,j,0)));
            Real Wscale = one + wall_cff*cffw*cffw;
            // ROMS -FCK(i,k)-FCK(i,k+1) == -FCK(k-1)-FCK(k) in REMORA indices
            BCK(i,j,k) = cff1*(one + two*Qdiss)    - FCK(i,j,k-1) - FCK(i,j,k);
            BCP(i,j,k) = cff1*(one + Wscale*Qdiss) - FCK(i,j,k-1) - FCK(i,j,k);
        });

        // ==============================================================
        // my25_corstep.F:638-697 -- Dirichlet surface/bottom values followed
        // by the two tridiagonal solves.  Serial in k both ways, so this is a
        // single column-per-thread kernel; CF is per-column scratch and is
        // reused by the gls sweep after the tke sweep has finished.
        // ==============================================================
        ParallelFor(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            // my25_corstep.F:642-653.  tke = q^2 = B1^(2/3) * u*^2, gls = 0.
            tke(i,j,N+1,nnew) = my_B1p2o3*Real(0.5)*
                std::sqrt((sustr(i,j,0)+sustr(i+1,j,0))*(sustr(i,j,0)+sustr(i+1,j,0))+
                          (svstr(i,j,0)+svstr(i,j+1,0))*(svstr(i,j,0)+svstr(i,j+1,0)));
            gls(i,j,N+1,nnew) = zero;
            tke(i,j,0,nnew) = my_B1p2o3*Real(0.5)*
                std::sqrt((bustr(i,j,0)+bustr(i+1,j,0))*(bustr(i,j,0)+bustr(i+1,j,0))+
                          (bvstr(i,j,0)+bvstr(i,j+1,0))*(bvstr(i,j,0)+bvstr(i,j+1,0)));
            gls(i,j,0,nnew) = zero;

            // my25_corstep.F:655-675 -- "tke" tridiagonal.  ROMS index N(ng)-1
            // is REMORA N; ROMS FCK(k) is REMORA FCK(k-1).
            Real cff = one/BCK(i,j,N);
            CF(i,j,N) = cff*FCK(i,j,N-1);
            tke(i,j,N,nnew) = cff*(tke(i,j,N,nnew) - FCK(i,j,N)*tke(i,j,N+1,nnew));
            for (int k=N-1; k>=1; k--) {
                cff = one/(BCK(i,j,k) - CF(i,j,k+1)*FCK(i,j,k));
                CF(i,j,k) = cff*FCK(i,j,k-1);
                tke(i,j,k,nnew) = cff*(tke(i,j,k,nnew) - FCK(i,j,k)*tke(i,j,k+1,nnew));
            }
            for (int k=1; k<=N; k++) {
                tke(i,j,k,nnew) = tke(i,j,k,nnew) - CF(i,j,k)*tke(i,j,k-1,nnew);
            }

            // my25_corstep.F:677-697 -- "gls" tridiagonal.  Note it uses BCP
            // for the diagonal but the SAME FCK off-diagonals as tke.
            cff = one/BCP(i,j,N);
            CF(i,j,N) = cff*FCK(i,j,N-1);
            gls(i,j,N,nnew) = cff*(gls(i,j,N,nnew) - FCK(i,j,N)*gls(i,j,N+1,nnew));
            for (int k=N-1; k>=1; k--) {
                cff = one/(BCP(i,j,k) - CF(i,j,k+1)*FCK(i,j,k));
                CF(i,j,k) = cff*FCK(i,j,k-1);
                gls(i,j,k,nnew) = cff*(gls(i,j,k,nnew) - FCK(i,j,k)*gls(i,j,k+1,nnew));
            }
            for (int k=1; k<=N; k++) {
                gls(i,j,k,nnew) = gls(i,j,k,nnew) - CF(i,j,k)*gls(i,j,k-1,nnew);
            }
        });

        // ==============================================================
        // my25_corstep.F:699-753 -- length scale, Galperin limit, Kantha and
        // Clayson stability functions, and the mixing coefficients.
        // ==============================================================
        ParallelFor(grow(bx,2,-1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // my25_corstep.F:706-716
            tke(i,j,k,nnew) = std::max(tke(i,j,k,nnew), my_qmin);
            gls(i,j,k,nnew) = std::max(gls(i,j,k,nnew), my_qmin);
            Real Ls_unlmt = gls(i,j,k,nnew)/tke(i,j,k,nnew);
            // NOTE: Real(0.0) rather than the namespace-scope `zero`.
            // std::max takes const references, which ODR-uses its argument;
            // a host constexpr has no device address, so nvcc rejects it with
            // "identifier \"zero\" is undefined in device code". Passing a
            // temporary avoids the ODR-use. (Assigning `= zero` is fine --
            // that only reads the value.)
            Real Ls_lmt = std::min(Ls_unlmt,
                                   my_lmax*std::sqrt(tke(i,j,k,nnew)/
                                       (std::max(Real(0.0),buoy2(i,j,k))+eps)));

            // my25_corstep.F:717-731 (KANTHA_CLAYSON branch at :727-728)
            Real Gh = std::min(my_Gh0, -buoy2(i,j,k)*Ls_lmt*Ls_lmt/tke(i,j,k,nnew));
            Real cff = one - my_Sh2*Gh;
            Real Sh = my_Sh1/cff;
            Real Sm = (my_B1pm1o3 + Sh*Gh*my_Sm4)/(one - my_Sm2*Gh);

            // my25_corstep.F:733-742: average q*l over the two time steps
            Real ql = Real(0.5)*(Ls_lmt*std::sqrt(tke(i,j,k,nnew)) +
                                 Lscale(i,j,k)*std::sqrt(tke(i,j,k,nstp)));
            Akv(i,j,k) = Akv_bak + ql*Sm;
            for (int n=0; n<ncons_local; n++) {
                Akt(i,j,k,n) = Akt_bak + ql*Sh;
            }

            // my25_corstep.F:744-748
            Akk(i,j,k) = Akk_bak + ql*my_Sq;

            // my25_corstep.F:750-752
            Lscale(i,j,k) = Ls_lmt;
        });
        // NOTE: my25_corstep.F never touches Akv/Akt/Akk at k=0 or k=N(ng);
        // they keep the values set in mod_mixing.F:986-1008 (zero at the two
        // end faces).  Deliberately not writing them here.
    }

    // my25_corstep.F:755-772 -- tkebc_tile + exchange_w3d_tile
    for (int icomp=0; icomp<3; icomp++) {
        FillPatch(lev, t_old[lev], *mf_tke, GetVecOfPtrs(vec_tke), tke_bc(), BdyVars::null, icomp, false, false);
        FillPatch(lev, t_old[lev], *mf_gls, GetVecOfPtrs(vec_gls), tke_bc(), BdyVars::null, icomp, false, false);
    }
    for (int icomp=0; icomp<ncons; icomp++) {
        FillPatch(lev, t_old[lev], *mf_Akt, GetVecOfPtrs(vec_Akt), zvel_bc(), BdyVars::null, icomp, false, false);
    }
    FillPatchNoBC(lev, t_old[lev], *mf_Akv, GetVecOfPtrs(vec_Akv), BdyVars::null);
    FillPatchNoBC(lev, t_old[lev], *mf_Akk, GetVecOfPtrs(vec_Akk), BdyVars::null);
}
