#ifndef REMORA_LoopBounds_H
#define REMORA_LoopBounds_H

#include <AMReX_Box.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>

/**
 * \brief Shrink a face-centred tile box to ROMS's momentum-equation loop range.
 *
 * ROMS never solves the momentum equation on the domain-edge normal face. Every
 * loop in ROMS 3.9 that WRITES a u-momentum quantity -- `rhs_ubar`, `rufrc`,
 * `ru`, `ubar(knew)`, `rubar(krhs)`, `u(nnew)` -- is bounded
 * `DO j=Jstr,Jend ; DO i=IstrU,Iend`, and every v-momentum one is
 * `DO j=JstrV,Jend ; DO i=Istr,Iend`. Checked term by term against
 * `step2d_LF_AM3.h` (919, 1219, 1269, 1320, 1407, 1847, 1888, 1940, 2000,
 * 2074, 2129, 2174, 2382 and their v twins), `rhs3d.F` (491, 542, 947, 1536),
 * `step3d_uv.F` (431, 470, 557) and `prsgrd32.h` (303, 367).
 *
 * `get_bounds.F:1373` sets `IstrU = Istr+1` on the WESTERN-edge tile only
 * (non-periodic), and `Iend` on the eastern-edge tile is `Lm`, so the u faces
 * at `Istr` and `Iend+1` are left to `u2dbc`/`u3dbc` instead. Interior tiles
 * keep `IstrU = Istr`. This function reproduces exactly that: shrink by one at
 * the DOMAIN edge, leave internal box edges alone.
 *
 * \param[in] nodal_bx  face-centred box, e.g. `mfi.nodaltilebox(dir)`
 * \param[in] domain    the CELL-CENTRED problem domain, `Geom(lev).Domain()`
 * \param[in] dir       the direction in which `nodal_bx` is nodal (0 for u, 1 for v)
 */
AMREX_FORCE_INLINE
amrex::Box
roms_mom_box (const amrex::Box& nodal_bx, const amrex::Box& domain, int dir)
{
    AMREX_ASSERT(nodal_bx.type()[dir] == amrex::IndexType::NODE);

    amrex::Box bx(nodal_bx);

    // The low test is against the domain's low CELL index because a box nodal
    // in `dir` shares its low index with the cell. The high test is against
    // `bigEnd+1` because the same box's high index is one PAST the last cell --
    // an earlier version of this trim compared against `bigEnd` directly, which
    // is unreachable for a box touching the high domain edge, so the eastern and
    // northern trims silently never fired.
    if (bx.smallEnd(dir) == domain.smallEnd(dir)) {
        bx.growLo(dir,-1);
    }
    if (bx.bigEnd(dir) == domain.bigEnd(dir) + 1) {
        bx.growHi(dir,-1);
    }
    return bx;
}

#endif
