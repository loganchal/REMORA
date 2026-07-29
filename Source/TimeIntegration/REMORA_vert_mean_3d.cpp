#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] phi_bx     box to update on
 * @param[in   ] ioff       x-direction offset
 * @param[in   ] joff       y-direction offset
 * @param[inout] phi        velocity (u or v)
 * @param[in   ] Dphi_avg1  time average of barotropic velocity
 * @param        DC         temporary
 * @param        CF         temporary
 * @param[in   ] pm_or_pn   1/dx or 1/dy
 * @param[in   ] msk        land-sea mask
 * @param[in   ] nnew       index of time step to update
 * @param[in   ] N          number of vertical levels
 */
void
REMORA::vert_mean_3d (const Box& phi_bx, const int ioff, const int joff,
                     const Array4<Real      >& phi,
                     const Array4<Real const>& Hz,
                     const Array4<Real const>& Dphi_avg1,
                     const Array4<Real      >& DC,
                     const Array4<Real      >& CF,
                     const Array4<Real const>& pm_or_pn,
                     const Array4<Real const>& msk,
                     const int nnew, const int N)
{
    BL_PROFILE("REMORA::vert_mean_3d()");

    // Operation order follows ROMS step3d_uv.f90:985-998 exactly, because at
    // this stage the port's only remaining difference from ROMS is last-bit
    // rounding and Sum(c*a_k) != c*Sum(a_k) in floating point.
    //
    // ROMS folds the metric into EVERY LEVEL before summing:
    //     cff     = 0.5*on_u(i,j)
    //     DC(i,k) = cff*(Hz(i,j,k)+Hz(i-1,j,k))
    //     DC(i,0) = DC(i,0)+DC(i,k)
    //     CF(i,0) = CF(i,0)+DC(i,k)*u(i,j,k,nnew)
    // then takes the reciprocal of the already-scaled total:
    //     DC(i,0) = 1/DC(i,0)
    //     CF(i,0) = DC(i,0)*(CF(i,0)-DU_avg1(i,j))
    //
    // REMORA previously summed bare thicknesses and applied on_u to the total,
    // which is algebraically the same and numerically is not.
    ParallelFor(makeSlab(phi_bx,2,0),
    [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        const Real on_u_or_om_v = two / (pm_or_pn(i-ioff,j-joff,0) + pm_or_pn(i,j,0));
        const Real cff = Real(0.5) * on_u_or_om_v;

        Real DCk = cff*(Hz(i-ioff,j-joff,0)+Hz(i,j,0));
        CF(i,j,-1) = DCk;
        DC(i,j,-1) = DCk*phi(i,j,0,nnew);

        for (int k=1; k<=N; k++) {
            DCk = cff*(Hz(i-ioff,j-joff,k)+Hz(i,j,k));
            CF(i,j,-1) += DCk;
            DC(i,j,-1) += DCk*phi(i,j,k,nnew);
        }
    });

    ParallelFor(makeSlab(phi_bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        const Real recip = one/CF(i,j,-1);
        DC(i,j,-1) = recip*(DC(i,j,-1) - Dphi_avg1(i,j,0)); // recursive
    });

    ParallelFor(phi_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        phi(i,j,k) -= DC(i,j,-1);
        phi(i,j,k) *= msk(i,j,0);
    });
}
