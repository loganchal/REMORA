#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] bx             box to apply climatology on
 * @param[in   ] ioff           offset in x-direction
 * @param[in   ] joff           offset in y-direction
 * @param[inout] var            variable to update
 * @param[in   ] var_old        variable to compare against for nudging
 * @param[in   ] var_clim       climatology value to nudge towards
 * @param[in   ] clim_coeff     nudging time scale (1/s)
 * @param[in   ] Hz             vertical cell height
 * @param[in   ] pm             1/dx
 * @param[in   ] pn             1/dy
 * @param[in   ] dt_lev         time step
 */
void
REMORA::apply_clim_nudg (const Box& bx,
                         int ioff, int joff,
                         const Array4<Real      >& var,
                         const Array4<Real const>& var_old,
                         const Array4<Real const>& var_clim,
                         const Array4<Real const>& clim_coeff,
                         const Array4<Real const>& Hz,
                         const Array4<Real const>& pm,
                         const Array4<Real const>& pn,
                         const Real dt_lev)
{
    BL_PROFILE("REMORA::apply_clim_nudg()");
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (ioff==1 || joff==1) {
            // ROMS rhs3d.f90:334-341 (u) and :346-353 (v):
            //     cff = 0.25*(M(i-1)+M(i)) * om_u * on_u
            //     ru  = ru + cff*(Hz(i-1)+Hz(i)) * (uclm - u)
            // i.e. ((((0.25*csum)*om)*on) * Hzsum) * delta. Splitting the 0.25
            // into two halves and grouping the metrics with the thickness gives
            // a different product tree, and this term feeds ru -> rufrc ->
            // rhs_ubar on every step (LnudgeM3CLM is T in the Moana deck).
            const Real om = two / (pm(i-ioff,j-joff,0)+pm(i,j,0));
            const Real on = two / (pn(i-ioff,j-joff,0)+pn(i,j,0));
            const Real cff = Real(0.25) * (clim_coeff(i-ioff,j-joff,k) + clim_coeff(i,j,k)) * om * on;
            var(i,j,k) += cff * (Hz(i-ioff,j-joff,k) + Hz(i,j,k)) * (var_clim(i,j,k) - var_old(i,j,k));
        } else {
            Real cff = Real(0.5) * (clim_coeff(i-ioff,j-joff,k) + clim_coeff(i,j,k));
            cff *= dt_lev;
            var(i,j,k) += cff * (var_clim(i,j,k) - var_old(i,j,k));
        }
    });
}
