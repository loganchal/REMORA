#ifndef REMORA_TIMEINTERP_H_
#define REMORA_TIMEINTERP_H_

#include <AMReX_REAL.H>
#include <cmath>

/**
 * \file REMORA_TimeInterp.H
 *
 * ROMS's linear time-interpolation weights, reproduced bit-for-bit.
 *
 * Why this exists
 * ---------------
 * Every time-varying input ROMS reads -- open-boundary data (set_ngfld.F),
 * surface forcing (set_2dfld.F) and river sources (set_ngfld.F again, via
 * set_data.F:183) -- is interpolated with the SAME arithmetic:
 *
 *     SecScale = 1000                       ! seconds -> milliseconds
 *     fac1 = ANINT((t_after  - t) * SecScale)
 *     fac2 = ANINT((t - t_before) * SecScale)
 *     fac  = 1 / (fac1 + fac2)
 *     fac1 = fac * fac1
 *     fac2 = fac * fac2
 *     Fout = fac1*F_before + fac2*F_after
 *
 * REMORA originally used the algebraically equivalent slope form
 *
 *     Fout = F_before + (t - t_before) * (F_after - F_before) / dt
 *
 * which is a different floating-point expression tree and differs from ROMS in
 * the last bit or two. Measured on the Moana July 2009 boundary file at the
 * step-1 model time, the two forms disagree on ~27% of the western boundary
 * cells, by up to 5.55e-17 in zeta (~4.6e-16 relative, one to two ulp) and
 * 2.78e-17 in ubar and vbar.
 *
 * That is small, and it is also NOT harmless: these values are re-imposed at
 * the open boundary on every one of the 44 barotropic sub-steps of every
 * baroclinic step, and the Flather condition turns a boundary free-surface
 * error directly into a face velocity, hence a volume flux. The parity campaign
 * measured the port exporting volume through the open boundary monotonically
 * from step 1, with a per-side transport bias of ~5 ulp and a growth of about
 * 1.9x per step -- and it was independent of the compiler's arithmetic mode
 * (fast-math on and off agreed to six significant figures) and of the boundary
 * DATA (a full regeneration of the GLORYS files moved the step-1 volume error
 * by less than one part in a thousand). A deterministic wrong number in the
 * boundary values, reused every sub-step, is exactly what that fingerprint
 * describes, and a differing expression tree is exactly such a number.
 *
 * Three details of the ROMS form are load-bearing and must not be "simplified":
 *
 *   1. ANINT rounds HALF AWAY FROM ZERO, which is std::round, not std::rint
 *      (round-half-to-even under the default mode) and not std::nearbyint.
 *   2. The normalisation divides by (fac1 + fac2), i.e. by the ROUNDED span,
 *      not by (t_after - t_before). When the times are not whole milliseconds
 *      these differ.
 *   3. fac1 and fac2 are each rounded separately after scaling by fac, so they
 *      do NOT in general sum to exactly 1. ROMS's interpolant is therefore not
 *      exactly a convex combination, and reproducing it means reproducing that.
 *
 * ROMS does this arithmetic in `real(dp)` (kind 8) regardless of the model's
 * working precision. amrex::Real is double in every build of this port, so the
 * plain type is used here rather than forcing a separate double.
 */
namespace REMORATimeInterp {

/** ROMS Fortran ANINT: nearest whole number, halves away from zero. */
AMREX_FORCE_INLINE
amrex::Real anint (amrex::Real x) { return std::round(x); }

/**
 * ROMS set_ngfld.F:83-100 / set_2dfld.F:94-121.
 *
 * @param[in ] time         time to interpolate to
 * @param[in ] time_before  time of the earlier snapshot
 * @param[in ] time_after   time of the later snapshot
 * @param[out] fac1         weight on the EARLIER snapshot
 * @param[out] fac2         weight on the LATER snapshot
 */
AMREX_FORCE_INLINE
void weights (amrex::Real time, amrex::Real time_before, amrex::Real time_after,
              amrex::Real& fac1, amrex::Real& fac2)
{
    constexpr amrex::Real SecScale = amrex::Real(1000.0); // seconds to milliseconds
    fac1 = anint((time_after  - time) * SecScale);
    fac2 = anint((time - time_before) * SecScale);
    const amrex::Real fac = amrex::Real(1.0) / (fac1 + fac2);
    fac1 = fac * fac1;
    fac2 = fac * fac2;
}

} // namespace REMORATimeInterp

#endif
