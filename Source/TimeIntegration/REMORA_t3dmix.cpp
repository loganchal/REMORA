#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] bx         box to operate on
 * @param[inout] state      scalar data
 * @param[  out] state_rhs  scalar data RHS
 * @param[in   ] diff2      diffusivity
 * @param[in   ] Hz         vertical cell height
 * @param[in   ] z_r        rho-point z coordinate
 * @param[in   ] pden       potential density anomaly (only needed for isopycnal mixing)
 * @param[in   ] pm         1/dx
 * @param[in   ] pn         1/dy
 * @param[in   ] msku       land-sea mask on u-points
 * @param[in   ] mskv       land-sea mask on v-points
 * @param[in   ] dt_lev     time step at level
 * @param[in   ] ncomp      number of components to do this calculation on
 * @param[in   ] N          number of vertical levels
 */
void
REMORA::t3dmix2   (const Box& bx,
                   const Array4<Real      >& state,
                   const Array4<Real      >& state_rhs,
                   const Array4<Real const>& diff2,
                   const Array4<Real const>& Hz,
                   const Array4<Real const>& z_r,
                   const Array4<Real const>& pden,
                   const Array4<Real const>& pm,
                   const Array4<Real const>& pn,
                   const Array4<Real const>& msku,
                   const Array4<Real const>& mskv,
                   const Real dt_lev, const int ncomp,
                   const int N) {
    if (solverChoice.harmonic_mixing_type == HarmonicMixingType::s) {
        t3dmix2_s(bx, state, state_rhs, diff2, Hz, pm, pn, msku, mskv, dt_lev, ncomp);
    } else if (solverChoice.harmonic_mixing_type == HarmonicMixingType::geopotential) {
        t3dmix2_geo(bx, state, state_rhs, diff2, Hz, z_r, pm, pn, msku, mskv, dt_lev, ncomp, N);
    } else if (solverChoice.harmonic_mixing_type == HarmonicMixingType::isopycnal) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pden.dataPtr() != nullptr,
            "t3dmix2_iso requires potential density (pden); it was not allocated");
        t3dmix2_iso(bx, state, state_rhs, diff2, Hz, z_r, pden, pm, pn, msku, mskv, dt_lev, ncomp, N);
    }
}


/**
 * @param[in   ] bx         box to operate on
 * @param[inout] state      scalar data
 * @param[  out] state_rhs  scalar data RHS
 * @param[in   ] diff2      diffusivity
 * @param[in   ] Hz         vertical cell height
 * @param[in   ] pm         1/dx
 * @param[in   ] pn         1/dy
 * @param[in   ] msku       land-sea mask on u-points
 * @param[in   ] mskv       land-sea mask on v-points
 * @param[in   ] dt_lev     time step at level
 * @param[in   ] ncomp      number of components to do this calculation on
 */
void
REMORA::t3dmix2_s (const Box& bx,
                   const Array4<Real      >& state,
                   const Array4<Real      >& state_rhs,
                   const Array4<Real const>& diff2,
                   const Array4<Real const>& Hz,
                   const Array4<Real const>& pm,
                   const Array4<Real const>& pn,
                   const Array4<Real const>& msku,
                   const Array4<Real const>& mskv,
                   const Real dt_lev, const int ncomp)
{
    BL_PROFILE("REMORA::t3dmix2_s()");
    //-----------------------------------------------------------------------
    //  Add in harmonic diffusivity s terms.
    //-----------------------------------------------------------------------

    Box xbx(bx); xbx.surroundingNodes(0);
    Box ybx(bx); ybx.surroundingNodes(1);

    FArrayBox fab_FX(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_FE(ybx,ncomp,The_Async_Arena());

    auto FX=fab_FX.array();
    auto FE=fab_FE.array();

    ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real pmon_u = (pm(i-1,j,0)+pm(i,j,0))/(pn(i-1,j,0)+pn(i,j,0));

        const Real cff = Real(0.25) * (diff2(i,j,0,n) + diff2(i-1,j,0,n)) * pmon_u;
        FX(i,j,k,n) = cff * (Hz(i,j,k) + Hz(i-1,j,k)) * (state_rhs(i,j,k,n)-state_rhs(i-1,j,k,n));
        FX(i,j,k,n) *= msku(i,j,0);
    });

    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real pnom_v = (pn(i,j-1,0)+pn(i,j,0))/(pm(i,j-1,0)+pm(i,j,0));

        const Real cff = Real(0.25)*(diff2(i,j,0,n)+diff2(i,j-1,0,n)) * pnom_v;
        FE(i,j,k,n) = cff * (Hz(i,j,k) + Hz(i,j-1,k)) * (state_rhs(i,j,k,n) - state_rhs(i,j-1,k,n));
        FE(i,j,k,n) *= mskv(i,j,0);
    });

    /*
     Time-step harmonic, S-surfaces diffusion term.
    */
    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real cff = dt_lev*pm(i,j,0)*pn(i,j,0);

        state(i,j,k,n) += cff * ( (FX(i+1,j  ,k,n)-FX(i,j,k,n))
                                 +(FE(i  ,j+1,k,n)-FE(i,j,k,n)) );
    });
}

/**
 * @param[in   ] bx         box to operate on
 * @param[inout] state      scalar data
 * @param[  out] state_rhs  scalar data RHS
 * @param[in   ] diff2      diffusivity
 * @param[in   ] Hz         vertical cell height
 * @param[in   ] z_r        rho-point z coordinates
 * @param[in   ] pm         1/dx
 * @param[in   ] pn         1/dy
 * @param[in   ] msku       land-sea mask on u-points
 * @param[in   ] mskv       land-sea mask on v-points
 * @param[in   ] dt_lev     time step at level
 * @param[in   ] ncomp      number of components to do this calculation on
 * @param[in   ] N          number of vertical levels
 */
void
REMORA::t3dmix2_geo(const Box& bx,
                   const Array4<Real      >& state,
                   const Array4<Real      >& state_rhs,
                   const Array4<Real const>& diff2,
                   const Array4<Real const>& Hz,
                   const Array4<Real const>& z_r,
                   const Array4<Real const>& pm,
                   const Array4<Real const>& pn,
                   const Array4<Real const>& msku,
                   const Array4<Real const>& mskv,
                   const Real dt_lev, const int ncomp, const int N)
{
    BL_PROFILE("REMORA::t3dmix2_geo()");
    //-----------------------------------------------------------------------
    //  Add in harmonic diffusivity s terms.
    //-----------------------------------------------------------------------

    Box xbx(bx); xbx.surroundingNodes(0);
    Box ybx(bx); ybx.surroundingNodes(1);
    Box zbx(bx); zbx.surroundingNodes(2);

    FArrayBox fab_dZdx(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_dTdx(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_dZde(ybx,ncomp,The_Async_Arena());
    FArrayBox fab_dTde(ybx,ncomp,The_Async_Arena());
    FArrayBox fab_dTdz(grow(zbx,IntVect(1,1,0)),ncomp,The_Async_Arena());
    FArrayBox fab_FS(grow(zbx,IntVect(1,1,0)),ncomp,The_Async_Arena());
    FArrayBox fab_FX(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_FE(ybx,ncomp,The_Async_Arena());

    auto dZdx = fab_dZdx.array();
    auto dTdx = fab_dTdx.array();
    auto dZde = fab_dZde.array();
    auto dTde = fab_dTde.array();
    auto dTdz = fab_dTdz.array();
    auto FS   = fab_FS.array();
    auto FX   = fab_FX.array();
    auto FE   = fab_FE.array();

    ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = Real(0.5) * (pm(i,j,0) + pm(i-1,j,0)) * msku(i,j,0);
        dZdx(i,j,k,n) = cff * (z_r(i,j,k) - z_r(i-1,j,k));
        dTdx(i,j,k,n)=cff*(state_rhs(i  ,j,k,n)-state_rhs(i-1,j,k,n));
    });
    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = Real(0.5) * (pn(i,j,0) + pn(i,j-1,0)) * mskv(i,j,0);
        dZde(i,j,k,n) = cff * (z_r(i,j,k) - z_r(i,j-1,k));
        dTde(i,j,k,n)=cff*(state_rhs(i,j,k,n)-state_rhs(i,j-1,k,n));
    });
    ParallelFor(makeSlab(grow(bx,IntVect(1,1,0)),2,0), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
    {
        dTdz(i,j,0,n) = Real(0.0);
        dTdz(i,j,N+1,n) = Real(0.0);
        FS(i,j,0,n) = Real(0.0);
        FS(i,j,N+1,n) = Real(0.0);
    });
    ParallelFor(grow(zbx,IntVect(1,1,-1)), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = one / (z_r(i,j,k)-z_r(i,j,k-1));
        dTdz(i,j,k,n) = cff * (state_rhs(i,j,k,n) - state_rhs(i,j,k-1,n));
    });

    //  Compute components of the rotated tracer flux (T m3/s) along
    //  geopotential surfaces.
    ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real on_u = two / (pn(i,j,0) + pn(i-1,j,0));
        Real cff = Real(0.25) * (diff2(i,j,0,n)+diff2(i-1,j,0,n)) * on_u;
        FX(i,j,k,n) = cff *
                       (Hz(i,j,k)+Hz(i-1,j,k))*
                       (dTdx(i,j,k,n)-
                        Real(0.5)*(std::min(dZdx(i,j,k,n),Real(0.0))*
                                   (dTdz(i-1,j,k  ,n)+
                                    dTdz(i  ,j,k+1,n))+
                                std::max(dZdx(i,j,k,n),Real(0.0))*
                                   (dTdz(i-1,j,k+1,n)+
                                    dTdz(i  ,j,k  ,n))));
    });
    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real om_v = two / (pm(i,j,0) + pm(i,j-1,0));
        Real cff = Real(0.25) * (diff2(i,j,0,n)+diff2(i,j-1,0,n)) * om_v;
        FE(i,j,k,n) = cff *
                       (Hz(i,j,k)+Hz(i,j-1,k))*
                       (dTde(i,j,k,n)-
                        Real(0.5)*(std::min(dZde(i,j,k,n),Real(0.0))*
                                   (dTdz(i,j-1,k  ,n)+
                                    dTdz(i,j  ,k+1,n))+
                                std::max(dZde(i,j,k,n),Real(0.0))*
                                   (dTdz(i,j-1,k+1,n)+
                                    dTdz(i,j  ,k  ,n))));
    });
    ParallelFor(grow(zbx,IntVect(0,0,-1)), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = Real(0.5) * diff2(i,j,0,n);
        Real cff1=std::min(dZdx(i  ,j,k-1,n),Real(0.0));
        Real cff2=std::min(dZdx(i+1,j,k  ,n),Real(0.0));
        Real cff3=std::max(dZdx(i  ,j,k  ,n),Real(0.0));
        Real cff4=std::max(dZdx(i+1,j,k-1,n),Real(0.0));
        FS(i,j,k,n) = cff *
                        (cff1*(cff1*dTdz(i,j,k,n)-dTdx(i  ,j,k-1,n))+
                         cff2*(cff2*dTdz(i,j,k,n)-dTdx(i+1,j,k  ,n))+
                         cff3*(cff3*dTdz(i,j,k,n)-dTdx(i  ,j,k  ,n))+
                         cff4*(cff4*dTdz(i,j,k,n)-dTdx(i+1,j,k-1,n)));
        cff1=std::min(dZde(i,j  ,k-1,n),Real(0.0));
        cff2=std::min(dZde(i,j+1,k  ,n),Real(0.0));
        cff3=std::max(dZde(i,j  ,k  ,n),Real(0.0));
        cff4=std::max(dZde(i,j+1,k-1,n),Real(0.0));
        FS(i,j,k,n)=FS(i,j,k,n)+
                            cff*
                            (cff1*(cff1*dTdz(i,j,k,n)-dTde(i,j  ,k-1,n))+
                             cff2*(cff2*dTdz(i,j,k,n)-dTde(i,j+1,k  ,n))+
                             cff3*(cff3*dTdz(i,j,k,n)-dTde(i,j  ,k  ,n))+
                             cff4*(cff4*dTdz(i,j,k,n)-dTde(i,j+1,k-1,n)));
    });

    // Time-step harmonic, geopotential diffusion term (m Tunits).

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff=dt_lev*pm(i,j,0)*pn(i,j,0);
        Real cff1=cff*(FX(i+1,j  ,k,n)-FX(i,j,k,n));
        Real cff2=cff*(FE(i  ,j+1,k,n)-FE(i,j,k,n));
        Real cff3=dt_lev*(FS(i,j,k+1,n)-FS(i,j,k,n));
        Real cff4=cff1+cff2+cff3;
        state(i,j,k,n)=state(i,j,k,n)+cff4;
    });
}

/**
 * Harmonic diffusion of tracers along isopycnic (constant potential density)
 * surfaces. This is a term-for-term port of ROMS 3.9
 * ROMS/Nonlinear/t3dmix2_iso.h (cppdefs TS_DIF2 + MIX_ISO_TS), as compiled for
 * the Moana hindcast. The preprocessed Moana source used as the specification is
 * moana_spec/t3dmix.f90 lines 517-800 (MASKING on; WET_DRY, DIFF_3DCOEF,
 * TS_MIX_CLIMA, TS_MIX_STABILITY, TS_MIX_MAX_SLOPE, TS_MIX_MIN_STRAT and
 * DIAGNOSTICS_TS all off). The ROMS parameters `small`, `slope_max` and
 * `strat_min` (t3dmix2_iso.h:182-184) are declared but unreachable with those
 * cppdefs, so only `eps` is carried over here.
 *
 * Index conventions. ROMS carries the vertical dependence through the recursive
 * two-slot (k1,k2) blocking of its K_LOOP; here the same quantities are stored
 * as full 3D arrays, exactly as REMORA::t3dmix2_geo does:
 *   - dRdx, dTdx, dRde, dTde live at rho levels k (ROMS slot k1 at iteration k).
 *   - dTdr, FS live at w faces, with face k the interface between cells k-1 and
 *     k, so face k is ROMS slot k1 and face k+1 is ROMS slot k2 at level k.
 *     Faces 0 and N+1 are the bottom and top and are set to zero
 *     (t3dmix.f90:695-701).
 * There is no vertical running dependence in the kernel, so every loop below is
 * a fully fused (i,j,k,n) ParallelFor.
 *
 * @param[in   ] bx         box to operate on
 * @param[inout] state      scalar data
 * @param[  out] state_rhs  scalar data RHS
 * @param[in   ] diff2      diffusivity
 * @param[in   ] Hz         vertical cell height
 * @param[in   ] z_r        rho-point z coordinates
 * @param[in   ] pden       potential density anomaly referenced to the surface
 * @param[in   ] pm         1/dx
 * @param[in   ] pn         1/dy
 * @param[in   ] msku       land-sea mask on u-points
 * @param[in   ] mskv       land-sea mask on v-points
 * @param[in   ] dt_lev     time step at level
 * @param[in   ] ncomp      number of components to do this calculation on
 * @param[in   ] N          number of vertical levels
 */
void
REMORA::t3dmix2_iso(const Box& bx,
                    const Array4<Real      >& state,
                    const Array4<Real      >& state_rhs,
                    const Array4<Real const>& diff2,
                    const Array4<Real const>& Hz,
                    const Array4<Real const>& z_r,
                    const Array4<Real const>& pden,
                    const Array4<Real const>& pm,
                    const Array4<Real const>& pn,
                    const Array4<Real const>& msku,
                    const Array4<Real const>& mskv,
                    const Real dt_lev, const int ncomp, const int N)
{
    BL_PROFILE("REMORA::t3dmix2_iso()");
    //-----------------------------------------------------------------------
    //  Compute horizontal harmonic diffusion along isopycnic surfaces.
    //  ROMS t3dmix2_iso.h:200-202 / t3dmix.f90:648-650
    //-----------------------------------------------------------------------

    // ROMS t3dmix2_iso.h:181 / t3dmix.f90:565
    const Real eps = Real(0.5);

    Box xbx(bx); xbx.surroundingNodes(0);
    Box ybx(bx); ybx.surroundingNodes(1);
    Box zbx(bx); zbx.surroundingNodes(2);

    FArrayBox fab_dRdx(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_dTdx(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_dRde(ybx,ncomp,The_Async_Arena());
    FArrayBox fab_dTde(ybx,ncomp,The_Async_Arena());
    FArrayBox fab_dTdr(grow(zbx,IntVect(1,1,0)),ncomp,The_Async_Arena());
    FArrayBox fab_FS(grow(zbx,IntVect(1,1,0)),ncomp,The_Async_Arena());
    FArrayBox fab_FX(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_FE(ybx,ncomp,The_Async_Arena());

    auto dRdx = fab_dRdx.array();
    auto dTdx = fab_dTdx.array();
    auto dRde = fab_dRde.array();
    auto dTde = fab_dTde.array();
    auto dTdr = fab_dTdr.array();
    auto FS   = fab_FS.array();
    auto FX   = fab_FX.array();
    auto FE   = fab_FE.array();

    // Horizontal density and tracer gradients at rho levels, u points.
    // ROMS t3dmix2_iso.h:218-249 / t3dmix.f90:666-679
    ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = Real(0.5) * (pm(i,j,0) + pm(i-1,j,0));
        cff = cff * msku(i,j,0);
        dRdx(i,j,k,n) = cff * (pden(i  ,j,k) -
                               pden(i-1,j,k));
        dTdx(i,j,k,n) = cff * (state_rhs(i  ,j,k,n) -
                               state_rhs(i-1,j,k,n));
    });

    // Horizontal density and tracer gradients at rho levels, v points.
    // ROMS t3dmix2_iso.h:250-281 / t3dmix.f90:680-693
    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = Real(0.5) * (pn(i,j,0) + pn(i,j-1,0));
        cff = cff * mskv(i,j,0);
        dRde(i,j,k,n) = cff * (pden(i,j  ,k) -
                               pden(i,j-1,k));
        dTde(i,j,k,n) = cff * (state_rhs(i,j  ,k,n) -
                               state_rhs(i,j-1,k,n));
    });

    // Bottom and top boundary faces: dTdr = FS = 0.
    // ROMS t3dmix2_iso.h:283-289 / t3dmix.f90:695-701
    ParallelFor(makeSlab(grow(bx,IntVect(1,1,0)),2,0), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
    {
        dTdr(i,j,0,n) = Real(0.0);
        dTdr(i,j,N+1,n) = Real(0.0);
        FS(i,j,0,n) = Real(0.0);
        FS(i,j,N+1,n) = Real(0.0);
    });

    // Interior w faces: tracer gradient with respect to potential density, and
    // the (clipped) inverse vertical density gradient times the layer thickness.
    // The stratification is clipped with eps = 0.5 kg/m3 (the TS_MIX_MAX_SLOPE /
    // TS_MIX_MIN_STRAT variants are not compiled in Moana).
    // ROMS t3dmix2_iso.h:290-333 / t3dmix.f90:702-715
    ParallelFor(grow(zbx,IntVect(1,1,-1)), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff1 = std::max(pden(i,j,k-1)-pden(i,j,k), eps);
        Real cff = Real(-1.0)/cff1;

        dTdr(i,j,k,n) = cff * (state_rhs(i,j,k  ,n) -
                               state_rhs(i,j,k-1,n));

        FS(i,j,k,n) = cff * (z_r(i,j,k) - z_r(i,j,k-1));
    });

    //  Compute components of the rotated tracer flux (T m4/s) along
    //  isopycnic surfaces.
    //  ROMS t3dmix2_iso.h:335-336 / t3dmix.f90:717-718

    // ROMS t3dmix2_iso.h:339-358 / t3dmix.f90:721-737.
    // on_u(i,j) = 2/(pn(i-1,j)+pn(i,j)) (ROMS metrics.F:411).
    ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real on_u = two / (pn(i-1,j,0) + pn(i,j,0));
        Real cff = Real(0.25) * (diff2(i,j,0,n)+diff2(i-1,j,0,n)) * on_u;
        FX(i,j,k,n) = cff *
                       (Hz(i,j,k)+Hz(i-1,j,k))*
                       (dTdx(i,j,k,n)-
                        Real(0.5)*(std::max(dRdx(i,j,k,n),Real(0.0))*
                                      (dTdr(i-1,j,k  ,n)+
                                       dTdr(i  ,j,k+1,n))+
                                   std::min(dRdx(i,j,k,n),Real(0.0))*
                                      (dTdr(i-1,j,k+1,n)+
                                       dTdr(i  ,j,k  ,n))));
    });

    // ROMS t3dmix2_iso.h:359-378 / t3dmix.f90:738-754.
    // om_v(i,j) = 2/(pm(i,j-1)+pm(i,j)) (ROMS metrics.F:447).
    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real om_v = two / (pm(i,j-1,0) + pm(i,j,0));
        Real cff = Real(0.25) * (diff2(i,j,0,n)+diff2(i,j-1,0,n)) * om_v;
        FE(i,j,k,n) = cff *
                       (Hz(i,j,k)+Hz(i,j-1,k))*
                       (dTde(i,j,k,n)-
                        Real(0.5)*(std::max(dRde(i,j,k,n),Real(0.0))*
                                      (dTdr(i,j-1,k  ,n)+
                                       dTdr(i,j  ,k+1,n))+
                                   std::min(dRde(i,j,k,n),Real(0.0))*
                                      (dTdr(i,j-1,k+1,n)+
                                       dTdr(i,j  ,k  ,n))));
    });

    // Vertical cross-term. ROMS applies this for 1 <= k <= N-1 in its own
    // numbering, i.e. every interior w face; the bottom (0) and top (N+1) faces
    // keep the zero set above.
    // ROMS t3dmix2_iso.h:379-406 / t3dmix.f90:755-780
    ParallelFor(grow(zbx,IntVect(0,0,-1)), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff1=std::max(dRdx(i  ,j,k-1,n),Real(0.0));
        Real cff2=std::max(dRdx(i+1,j,k  ,n),Real(0.0));
        Real cff3=std::min(dRdx(i  ,j,k  ,n),Real(0.0));
        Real cff4=std::min(dRdx(i+1,j,k-1,n),Real(0.0));
        Real cff=cff1*(cff1*dTdr(i,j,k,n)-dTdx(i  ,j,k-1,n))+
                 cff2*(cff2*dTdr(i,j,k,n)-dTdx(i+1,j,k  ,n))+
                 cff3*(cff3*dTdr(i,j,k,n)-dTdx(i  ,j,k  ,n))+
                 cff4*(cff4*dTdr(i,j,k,n)-dTdx(i+1,j,k-1,n));
        cff1=std::max(dRde(i,j  ,k-1,n),Real(0.0));
        cff2=std::max(dRde(i,j+1,k  ,n),Real(0.0));
        cff3=std::min(dRde(i,j  ,k  ,n),Real(0.0));
        cff4=std::min(dRde(i,j+1,k-1,n),Real(0.0));
        cff=cff+
            cff1*(cff1*dTdr(i,j,k,n)-dTde(i,j  ,k-1,n))+
            cff2*(cff2*dTdr(i,j,k,n)-dTde(i,j+1,k  ,n))+
            cff3*(cff3*dTdr(i,j,k,n)-dTde(i,j  ,k  ,n))+
            cff4*(cff4*dTdr(i,j,k,n)-dTde(i,j+1,k-1,n));

        FS(i,j,k,n)=Real(0.5)*cff*diff2(i,j,0,n)*FS(i,j,k,n);
    });

    // Time-step harmonic, isopycnic diffusion term (m Tunits).
    // ROMS t3dmix2_iso.h:408-425 / t3dmix.f90:782-794
    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff=dt_lev*pm(i,j,0)*pn(i,j,0);
        Real cff1=cff*(FX(i+1,j  ,k,n)-FX(i,j,k,n));
        Real cff2=cff*(FE(i  ,j+1,k,n)-FE(i,j,k,n));
        Real cff3=dt_lev*(FS(i,j,k+1,n)-FS(i,j,k,n));
        Real cff4=cff1+cff2+cff3;
        state(i,j,k,n)=state(i,j,k,n)+cff4;
    });
}
