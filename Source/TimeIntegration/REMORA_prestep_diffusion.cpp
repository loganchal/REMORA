#include <REMORA.H>

using namespace amrex;

/**
 * Called from prestep. The tracer update is a bit different from the u,v updates so we test
 * for it, but checking if ioff=0 and joff=0. In some cases, though, we
 * can recover the tracer update from the generic one by setting those indices.
 * Setting icc and ntfirst identically for the tracers should be equivalent
 * to setting ioff=0 and joff=0
 *
 * @param[in   ] vel_bx   tile box
 * @param[in   ] gbx      grown tile box
 * @param[in   ] ioff     offset in x direction
 * @param[in   ] joff     offset in y direction
 * @param[  out] vel      velocity or scalar to update
 * @param[in   ] vel_old  velocity or scalar at last time
 * @param[inout] rvel     velocity or scalar RHS
 * @param[in   ] Hz       vertical cell height
 * @param[in   ] Akv      vertical viscosity coefficient
 * @param        FC       temporary
 * @param[in   ] sstr     surface flux
 * @param[in   ] bstr     bottom flux
 * @param[in   ] z_r      z coordinates at rho points
 * @param[in   ] pm       1/dx
 * @param[in   ] pn       1/dy
 * @param[in   ] iic      which time step we're on
 * @param[in   ] ntfirst  what is the first time step?
 * @param[in   ] nnew     index of time step to update
 * @param[in   ] nstp     index of last time step
 * @param[in   ] nrhs     index of RHS component
 * @param[in   ] N        number of vertical levels
 * @param[in   ] lambda   weighting coefficient for the newest (implicit) time step derivatives
 * @param[in   ] dt_lev   time step at this level
 */

void
REMORA::prestep_diffusion (const Box& vel_bx, const Box& gbx,
                          const int ioff, const int joff,
                          const Array4<Real      >& vel,
                          const Array4<Real const>& vel_old,
                          const Array4<Real      >& rvel,
                          const Array4<Real const>& Hz,
                          const Array4<Real const>& Akv,
                          const Array4<Real      >& FC,
                          const Array4<Real const>& sstr,
                          const Array4<Real const>& bstr,
                          const Array4<Real const>& z_r,
                          const Array4<Real const>& pm,
                          const Array4<Real const>& pn,
                          const int iic, const int ntfirst, const int nnew, int nstp, int nrhs, int N,
                          const Real lambda, const Real dt_lev,
                          const Array4<Real const>& srflx,
                          const Array4<Real const>& z_w,
                          const bool solar_source)
{
    BL_PROFILE("REMORA::prestep_diffusion()");

    BoxArray ba_gbxvel = intersect(BoxArray(vel_bx), gbx);
    AMREX_ASSERT((ba_gbxvel.size() == 1));
    Box gbxvel = ba_gbxvel[0];

    //
    //  Weighting coefficient for the newest (implicit) time step derivatives
    //  using either a Crank-Nicolson implicit scheme (lambda=Real(0.5)) or a
    //  backward implicit scheme (lambda=one).
    //

    Real oml_dt = dt_lev*(one-lambda);
    //N is one less than ROMS

    //  Except the commented out part means lambda is always one
    if (verbose > 1) {
        amrex::Print() << "in update_vel_3d with box " << vel_bx << std::endl;
        Print() << "vel old " << Box(vel_old) << std::endl;
        Print() << "Akv " << Box(Akv) << std::endl;
    }
    ParallelFor(grow(surroundingNodes(vel_bx,2),IntVect(0,0,-1)),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
            Real cff = one / ( z_r(i,j,k)+z_r(i-ioff,j-joff,k)
                              -z_r(i,j,k-1)-z_r(i-ioff,j-joff,k-1));
            FC(i,j,k) = oml_dt * cff * (vel_old(i,j,k,nstp)-vel_old(i,j,k-1,nstp)) *
                                           (Akv(i,j,k)     +Akv(i-ioff,j-joff,k));
    });
    // Shortwave penetration (ROMS SOLAR_SOURCE, pre_step3d.F + lmd_swfrac.F):
    // add dt*srflx*swdk at every interior W-face, where swdk is the
    // Paulson & Simpson (1977) two-band decay for the Jerlov water type.
    // ROMS's srflx is kinematic (degC m/s); REMORA's srflx MultiFab is W/m2,
    // so the 1/(rho0*Cp) scale is applied here. The surface face (below)
    // carries the full stflx, so the top cell absorbs srflx*(1-swdk(N)).
    if (solar_source) {
        static constexpr Real lmd_mu1[9] = {0.35, 0.6, 1.0, 1.5, 1.4,
                                            0.42, 0.37, 0.33, 0.00468592};
        static constexpr Real lmd_mu2[9] = {23.0, 20.0, 17.0, 14.0, 7.9,
                                            5.13, 3.54, 2.34, 1.51};
        static constexpr Real lmd_r1[9]  = {0.58, 0.62, 0.67, 0.77, 0.78,
                                            0.57, 0.57, 0.57, 0.55};
        const int jwt = solverChoice.jerlov_water_type - 1;
        // ROMS lmd_swfrac_tile with Zscale=-1 and Z = z_w(N)-z_w(k) >= 0
        const Real fac1 = -one/lmd_mu1[jwt];
        const Real fac2 = -one/lmd_mu2[jwt];
        const Real fac3 = lmd_r1[jwt];
        const Real Hscale2 = one/(solverChoice.rho0*Cp);
        ParallelFor(grow(surroundingNodes(vel_bx,2),IntVect(0,0,-1)),
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real zdep = z_w(i,j,N+1) - z_w(i,j,k);
            Real swdk = std::exp(zdep*fac1)*fac3 +
                        std::exp(zdep*fac2)*(one-fac3);
            FC(i,j,k) += dt_lev*(srflx(i,j,0)*Hscale2)*swdk;
        });
    }

    ParallelFor(makeSlab(vel_bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int ) {
        FC(i,j,0) = dt_lev*bstr(i,j,0);
        FC(i,j, N+1) = dt_lev*sstr(i,j,0);
    });

    if (iic==ntfirst || (ioff==0 and joff==0)) {
        ParallelFor(gbxvel,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff1=vel_old(i,j,k,nstp)*Real(0.5)*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
            Real cff2=FC(i,j,k+1)-FC(i,j,k);
            vel(i,j,k,nnew)=cff1+cff2;
        });
    } else if (iic==ntfirst+1) {
        ParallelFor(gbxvel,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff1=vel_old(i,j,k,nstp)*Real(0.5)*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
            Real cff2=FC(i,j,k+1)-FC(i,j,k);
            Real DC = Real(0.25) * dt_lev * (pm(i,j,0)+pm(i-ioff,j-joff,0))
                                  * (pn(i,j,0)+pn(i-ioff,j-joff,0));
            Real cff3 = DC * Real(0.5);
            int indx=nrhs ? 0 : 1;
            Real r_swap= rvel(i,j,k,indx);
            rvel(i,j,k,indx) = rvel(i,j,k,nrhs);
            rvel(i,j,k,nrhs) = r_swap;
            vel(i,j,k,nnew)=cff1- cff3*rvel(i,j,k,indx)+ cff2;
        });
    } else {
        ParallelFor(gbxvel,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff1 =  Real(5.0)/Real(12.0);
            Real cff2 = Real(16.0)/Real(12.0);
            Real cff3=vel_old(i,j,k,nstp)*Real(0.5)*(Hz(i,j,k)+Hz(i-ioff,j-joff,k));
            Real cff4=FC(i,j,k+1)-FC(i,j,k);
            Real DC = Real(0.25) * dt_lev * (pm(i,j,0)+pm(i-ioff,j-joff,0))
                                  * (pn(i,j,0)+pn(i-ioff,j-joff,0));

            int indx=nrhs ? 0 : 1;
            Real r_swap= rvel(i,j,k,indx);
            rvel(i,j,k,indx) = rvel(i,j,k,nrhs);
            rvel(i,j,k,nrhs) = r_swap;

            vel(i,j,k,nnew) = cff3 + DC*(cff1*rvel(i,j,k,nrhs)-
                                                cff2*rvel(i,j,k,indx))+cff4;
            rvel(i,j,k,nrhs) = zero;
        });
    }
}
