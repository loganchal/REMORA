# MY2.5 (Mellor–Yamada level 2.5) port notes

Port of ROMS 3.9's `MY25_MIXING` closure into REMORA as a third
`VertMixingType`, selected with `remora.vertical_mixing_type = MY25`.

## Sources

| Role | Path |
| --- | --- |
| Spec (preprocessed Moana build, all `#ifdef`s resolved) | `external/moana_spec/my25_prestep.f90`, `my25_corstep.f90`, `tkebc_im.f90` |
| Annotated originals | `external/roms-3.9/ROMS/Nonlinear/my25_prestep.F`, `my25_corstep.F` |
| Closure constants | `external/roms-3.9/ROMS/Modules/mod_scalars.F` |
| Initial values | `external/roms-3.9/ROMS/Modules/mod_mixing.F` (`initialize_mixing`) |
| Structural template | `Source/TimeIntegration/REMORA_gls.cpp` |
| Port | `Source/TimeIntegration/REMORA_my25.cpp` |

## CPP guards that survived preprocessing

Checked by diffing each `#ifdef` region of the annotated `.F` against the
preprocessed `.f90`:

| Guard | State in the Moana build | Consequence |
| --- | --- | --- |
| `MY25_MIXING` | **on** | the whole module compiles |
| `KANTHA_CLAYSON` | **on** | `my_Sh2`/`my_Sm4` take their KC forms; `Sm` numerator is `my_B1pm1o3`, *not* `my_Sm3` (`my25_corstep.F:727-728`) |
| `N2S2_HORAVG` | **on** | two-stage 2×2 horizontal smoothing of N² and S² (`my25_corstep.F:279-338`) |
| `RI_SPLINES` | **on** | shear at W-points from a parabolic spline solve (`my25_corstep.F:220-255`), not the simple finite difference in the `#else` at `:256-267` |
| `MASKING` | **on** | `umask`/`vmask` multiply every advective gradient |
| `K_C2ADVECTION` | off | — |
| `K_C4ADVECTION` | off | corrector uses **third-order upstream** horizontal advection (`my25_corstep.F:412-438`); the *prestep* separately uses fourth-order centred (`my25_prestep.F:208-291`), which is its `#else` branch and unrelated to `K_C4ADVECTION` |
| `LIMIT_VDIFF` / `LIMIT_VVISC` | off | no upper clamp on Akv/Akt (`my25_corstep.F:754-769` not compiled) |
| `PERFECT_RESTART` | off | tke/gls/Lscale are always re-seeded internally, never read from restart |

## Constants table

Base constants, `mod_scalars.F:1754-1767`. All nine `my_A*`/`my_B*`/`my_C*`/`my_E*`
already existed in `SolverChoice` with matching values, so no new inputs were added.
They are ROMS compile-time parameters, not `roms.in` entries, and are deliberately
**not** exposed as runtime options.

| Symbol | Value | `mod_scalars.F` |
| --- | --- | --- |
| `my_A1` | 0.92 | :1754 |
| `my_A2` | 0.74 | :1755 |
| `my_B1` | 16.6 | :1756 |
| `my_B2` | 10.1 | :1757 |
| `my_C1` | 0.08 | :1758 |
| `my_C2` | 0.7 | :1759 |
| `my_C3` | 0.2 | :1760 |
| `my_E1` | 1.8 | :1761 |
| `my_E2` | 1.33 | :1762 |
| `my_Gh0` | 0.0233 | :1763 |
| `my_Sq` | 0.2 | :1764 |
| `my_dtfac` | 0.05 | :1765 (Asselin filter — unused by prestep/corstep) |
| `my_lmax` | 0.53 | :1766 |
| `my_qmin` | 1.0e-8 | :1767 |

Derived combinations, `mod_scalars.F:3761-3775`:

| Symbol | Expression | Value | Line | Note |
| --- | --- | --- | --- | --- |
| `my_B1p2o3` | `B1^(2/3)` | 6.4938… | :3761 | surface/bottom `q²` |
| `my_B1pm1o3` | `1/B1^(1/3)` | 0.39238… | :3762 | KC `Sm` numerator |
| `my_E1o2` | `0.5*E1` | 0.9 | :3763 | **unused** in these two routines |
| `my_Sm1` | `A1*A2*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1))` | — | :3764-3766 | **unused** (Canuto path only) |
| `my_Sm2` | `9*A1*A2` | 6.1272 | :3767 | |
| `my_Sh1` | `A2*(1-6*A1/B1)` | 0.49347… | :3768 | |
| `my_Sh2` | `3*A2*(6*A1 + B2*(1-C3))` | 30.190… | :3770 | **KANTHA_CLAYSON** form |
| `my_Sm3` | `A1*(1-3*C1-6*A1/B1)` | — | :3774 | **not compiled** (KC branch) |
| `my_Sm4` | `18*A1² + 9*A1*A2*(1-C2)` | 17.078… | :3771 | **KANTHA_CLAYSON** form |

Local literals: `Gamma = 1/6` (`my25_prestep.F:139`), `Gadv = 1/3`
(`my25_corstep.F:189`), `eps = 1e-10` (`my25_corstep.F:190`),
wall prefactor `my_E2/vonKar²` (`my25_corstep.F:598`), buoyancy dead-band
`-5.0e-5` (`my25_corstep.F:604-608`).

The two KC entries are the crux of the port. REMORA's `GLS_StabilityType`
offers only `Canuto_A`, `Canuto_B`, `Galperin`; `Galperin` is the *other*
MY2.5 branch (`KANTHA_CLAYSON` undefined) and is **not** a drop-in
substitute — it uses `my_Sm3` and the non-KC `my_Sh2`. The KC forms are
therefore written out explicitly in `REMORA_my25.cpp` rather than routed
through the GLS stability machinery.

## Index mapping

`N` in REMORA is the top rho index, so there are `N+1` rho levels and `N+2`
z-nodes.

```
N_ROMS            = N + 1
w-point  index    k_REMORA = k_ROMS            (0 = bottom, N+1 = surface)
rho-cell index    k_REMORA = k_ROMS - 1        =>  Hz_ROMS(k) == Hz(k-1)
```

Consequences used throughout:

| ROMS | REMORA |
| --- | --- |
| `DO k=1,N(ng)-1` (interior W) | `grow(bx,2,-1)`, i.e. `k = 1..N` |
| `tke(i,j,N(ng),nnew)` (surface) | `tke(i,j,N+1,nnew)` |
| `tke(i,j,0,nnew)` (bottom) | `tke(i,j,0,nnew)` |
| `0.5*(Hz(i,j,k)+Hz(i,j,k+1))` at W-point k | `0.5*(Hz(i,j,k-1)+Hz(i,j,k))` |
| `0.5*(Huon(i,j,k)+Huon(i,j,k+1))` | `0.5*(Huon(i,j,k)+Huon(i,j,k-1))` |
| `z_w(i,j,N(ng))` | `z_w(i,j,N+1)` |
| `FCK(i,k)`, `k=1..N(ng)` (rho) | `FCK(i,j,k-1)`, `k_REMORA = 0..N` |
| `BCK(i,k) - FCK(i,k) - FCK(i,k+1)` | `BCK(i,j,k) - FCK(i,j,k-1) - FCK(i,j,k)` |
| `BCK(i,N(ng)-1)` (top of the solve) | `BCK(i,j,N)` |

Horizontal: ROMS rho `i=1..Lm` ↔ REMORA `i=0..nx-1`; ROMS u-face `i` ↔
REMORA u-face `i-1`. So ROMS `Istr` is REMORA `dlo.x`, ROMS `Iend+1` is
REMORA `dhi.x+1`.

## Structural line-by-line check

### `my25_prestep.F` → `REMORA_my25.cpp::my25_prestep`

`my25_prestep.F` is **byte-identical** to `gls_prestep.F` over the whole
computational body (verified: `diff <(sed -n '180,424p' gls_prestep.F)
<(sed -n '180,424p' my25_prestep.F)` differs only in a comment character and
the subroutine name). MY2.5 diverges from GLS only in the corrector.

| ROMS `my25_prestep.F` | `REMORA_my25.cpp` | Notes |
| --- | --- | --- |
| :208-219 | :121-127 | masked `grad`/`gradL` on u-points |
| :222-237 | :129-141 | west/east one-sided gradient copy — **see deviation 1** |
| :239-252 | :142-148 | `XF`, `FX`, `FXL` |
| :253-291 | :151-178 | same construction in η |
| :293-305 | :185-197 | AM3 weights `cff1/cff2/cff3`, `indx` |
| :306-325 | :200-212 | `Hz_half`, level-"3" update, `nnew = Hz*nstp` stash |
| :327-374 | :214-243 | vertical advective flux (4th-order interior, one-sided ends) |
| :376-393 | :245-256 | time-step vertical advection, divide by `Hz_half` |
| :396-411 | :259-263 | `tkebc_tile` + `exchange_w3d_tile` → `FillPatch` |

### `my25_corstep.F` → `REMORA_my25.cpp::my25_corrector`

| ROMS `my25_corstep.F` | `REMORA_my25.cpp` | Notes |
| --- | --- | --- |
| :220-255 (`RI_SPLINES`) | :356-401 | spline shear at W-points; strictly column-per-thread |
| :270-278 | folded into :490-500 | `buoy2 = bvf` load |
| :285-318 (`N2S2_HORAVG`) | :404-407 | edge/corner one-sided copy of `shear2` → `physbcs` foextrap |
| :319-337 | :484-501 | two-stage 2×2 average of `buoy2`, `shear2` |
| :373-397 | :510-532 | masked `gradK`/`gradP`, boundary copies |
| :407-438 | :528-543 | `curvK`/`curvP`, upstream pick, `FXK`/`FXP` |
| :441-461 | :545-563 | η gradients, curvature |
| :481-502 | :565-574 | `FEK`/`FEP` |
| :505-521 | :576-587 | time-step horizontal advection — **no floor** (deviation 2) |
| :523-566 | :589-613 | vertical advective flux |
| :568-578 | :615-621 | time-step vertical advection — **no floor** |
| :585-593 | :624-637 | `FCK = -0.5*dt*(Akk(k)+Akk(k-1))/Hz(k)`, `CF = 0` — **not** zeroed at the ends (deviation 3) |
| :601-612 | :644-649 | buoyancy dead-band, `Qprod` |
| :613-617 | :650-651 | `Ls_unlmt` from the old time step |
| :618-625 | :653-656 | production into `tke` (×2) and `gls` (×`E1*Ls_unlmt`) |
| :626-636 | :658-667 | `Qdiss`, wall function `Wscale`, `BCK`/`BCP` |
| :642-653 | :677-685 | Dirichlet surface/bottom: `q² = B1^(2/3)·½|τ|`, `q²l = 0` |
| :655-675 | :687-699 | `tke` tridiagonal |
| :677-697 | :701-713 | `gls` tridiagonal — diagonal `BCP`, off-diagonal **`FCK`** (deviation 4) |
| :706-716 | :722-728 | `my_qmin` floor, `Ls_unlmt`, Galperin length limit |
| :717-728 | :730-734 | `Gh`, `Sh`, KC `Sm` |
| :733-742 | :736-743 | two-time-level `ql` average, `Akv`, `Akt` |
| :744-748 | :744-745 | `Akk` |
| :750-752 | :747-748 | `Lscale` |
| :755-772 | :755-763 | `tkebc_tile` + exchanges |

## Deviations between spec and expectation

1. **Prestep boundary gradient (fixed relative to the GLS template).**
   `REMORA_gls.cpp:95` / `:122` test `i == dlo.x-1` / `j == dlo.y-1` inside a
   `ParallelFor` over `grow(xbx,…)` / `grow(ybx,…)`, whose lower face index is
   `dlo.x` / `dlo.y`. The branch is therefore **unreachable** and the west/south
   one-sided copy of `my25_prestep.F:222-229` / `:263-270` never happens — the
   flux at the first face falls through to a ghost-cell difference instead.
   (The high side, `i == dhi.x+1`, is correct because that face *is* in the box.)
   The MY2.5 port uses `i == dlo.x` / `j == dlo.y`, which reproduces ROMS.
   The masks `umask(Istr,j)` / `vmask(i,Jstr)` are applied to the copied
   gradient as ROMS does; the GLS template drops them on the boundary rows.
   `REMORA_gls.cpp` is left untouched (GLS stays bit-identical), so GLS and
   MY2.5 now differ here even though ROMS's two prestep routines are identical.
   Only reachable with a non-periodic boundary; harmless for the periodic-x
   Upwelling case, but it matters for Moana. **Worth reporting upstream.**

2. **No `gls_Kmin`/`gls_Pmin` floors in the corrector.** `gls_corstep.F`
   clamps `tke`/`gls` after each advective update; `my25_corstep.F` does not.
   The only clamp in MY2.5 is `my_qmin = 1e-8` at `:706-707`, applied once,
   after the tridiagonal solves. The port omits the extra clamps that
   `REMORA_gls.cpp:617,621,652,654,802,803` carry.

3. **`FCK` is not zeroed at the end faces.** `gls_corstep.F:700-703` sets
   `FCK(1) = FCK(N) = FCP(1) = FCP(N) = 0` because GLS uses flux (Neumann)
   surface/bottom conditions. MY2.5 uses **Dirichlet** conditions, so those
   coefficients are live and carry the boundary values into the RHS of the
   first and last interior rows. Copying the GLS zeroing would silently
   decouple the boundary.

4. **One diffusion coefficient, not two.** MY2.5 has no `Akp`: the `q²l`
   equation is diffused with the same `Akk` as `q²` (`my25_corstep.F:696`
   defines only `FCK`, and the `gls` solve at `:677-697` uses `FCK`
   throughout, with only the diagonal `BCP` differing from `BCK`).
   `my25_corrector` therefore takes no `mf_Akp` argument and never writes
   `vec_Akp`. `vec_Akp` is still allocated and checkpointed unconditionally
   by existing code, so nothing else changes.

5. **`Akv`/`Akt`/`Akk` are never written at `k=0` and `k=N+1`.**
   `my25_corstep.F` loops `k=1,N(ng)-1` only, unlike `gls_corstep.F:806-828`
   which sets surface/bottom values from the roughness lengths. In ROMS those
   two faces keep the `IniVal = 0` from `mod_mixing.F:986-1008` for the whole
   run. `init_my25_vmix` reproduces that (zero at both end faces) and the
   corrector deliberately leaves them alone. Confirmed harmless: ROMS's
   momentum/tracer vertical solves only read `k=1..N(ng)-1`.

6. **No `Zob`/`ZoBot` dependence.** MY2.5's boundary conditions come straight
   from the stress magnitude, so `vec_ZoBot` is not needed. Its allocation
   and checkpoint gating (`REMORA_make_new_level.cpp`, `REMORA_Checkpoint.cpp`)
   were left GLS-only on purpose.

7. **`my_dtfac`, `my_E1o2`, `my_Sm1`, `my_Sm3` are unused.** `my_dtfac` is the
   Asselin filter coefficient used elsewhere in ROMS; `my_E1o2` and `my_Sm1`
   are not referenced by either routine; `my_Sm3` belongs to the non-KC branch.
   Listed above for completeness and deliberately not wired in.

## GPU safety

Every running k-dependence is confined to one thread per column:

* spline shear (`REMORA_my25.cpp:378-401`) — forward elimination + back
  substitution, one `ParallelFor` over a `makeSlab(bx,2,0)` box with the `k`
  loops inside;
* both tridiagonal solves (`:675-714`) — Dirichlet BC assignment, `tke` sweep
  and `gls` sweep are all in a **single** column kernel, so the shared `CF`
  scratch is reused safely (per-column storage, sequential within the thread);
* vertical advective flux (`:592-613`, and `:220-243` in the prestep) — column
  kernels, purely for launch economy; no serial dependence.

All other kernels are pointwise in `k` over fused `(i,j,k)` `ParallelFor`s.
No fused index carries a serial dependence. `FCK`/`FCP` are reused between the
advection and diffusion stages, but only across separate `ParallelFor` launches
in the same stream, which are ordered.

## Verification

| Check | Result |
| --- | --- |
| Clean CPU rebuild (`cmake` + `make -j8`, AppleClang, Release) | pass, no warnings from the new file |
| `vertical_mixing_type = MY25`, Upwelling, 100 steps | runs, exit 0 |
| …500 steps | runs, exit 0 |
| `amrex_fnan` on the 100- and 500-step plotfiles | clean (temp, salt, u, v, w) |
| Akv/Akt/Akk positive and bounded | Akv ∈ [0, 2.0e-3] at step 120, [0, 2.7e-2] at step 500; Akt ∈ [0, 3.4e-2]; Akk ∈ [0, 1.2e-2]. Minimum 0 is the deliberately-unwritten `k=0`/`k=N+1` faces (deviation 5); interior values are ≥ background |
| `Lscale` bounded | 0 → 2.23 m over 500 steps, monotone, no runaway |
| `tke`, `gls` bounded | `tke` ≤ 1.3e-3 m²/s², floored at `my_qmin = 1e-8`; `gls` ≤ 2.3e-3 |
| Physical plausibility | MY2.5 tracks the GLS solution closely (u_min −0.109 vs −0.115; temp_max 21.728 vs 21.736 at step 100) and both differ markedly from constant-Akv analytic mixing (−0.033, 21.931) — the expected signature of a real closure |
| Default (analytic) run bit-identical vs `moana/dev` build | `PLOTFILE AGREE`, absolute and relative error exactly 0 on all five variables |
| GLS run bit-identical vs `moana/dev` build | `PLOTFILE AGREE`, absolute and relative error exactly 0 on all five variables |

### Pre-existing test-suite state (not caused by this port)

`ctest` reports **15 of 16 regression tests failing on bare `moana/dev`
(8f95ac3)**, before any MY2.5 change. Confirmed by building `moana/dev` in the
same worktree and re-running: identical 15/16 failure set. The cause is
`amrex.fpe_trap_invalid = 1` in the test `.i` files tripping a floating-point
invalid inside `Problem::init_analytic_bathymetry` → `REMORA::set_bathymetry`
→ `REMORA::init_only`, i.e. during initialization, before any time step and in
a call path this port does not touch.

Running the same `Upwelling.i` with `amrex.fpe_trap_invalid=0` completes and
lands close to the stored gold file (temp 3.7e-8 abs / 1.7e-9 rel, salt 7.8e-14,
u 2.1e-9, v 1.8e-8), the residual being `moana/dev`'s intentional parity
divergence from upstream. So the suite is informative again once the trap is
off, and this port leaves it exactly where it was.

### Not verified

* **No GPU build or run.** Machine is CPU-only (AppleClang, no CUDA/HIP).
  GPU safety is by construction and inspection (section above), not measured.
  In particular the column kernels have not been checked for register pressure
  or the CPU/GPU reproducibility caveat flagged at `REMORA_gls.cpp:813-814`.
* **No numerical comparison against ROMS output.** No MY2.5 ROMS reference run
  was available, so the port is verified structurally (line-by-line against the
  preprocessed spec) and for physical plausibility, not term-by-term numerically.
  A single-column or short Moana-grid A/B against ROMS is the obvious next step.
* **MPI decomposition.** Runs were single-rank. The tile-boundary behaviour of
  the `physbcs` foextrap on `shear2` and of the `dlo.x`/`dhi.x+1` gradient
  tests has not been exercised across ranks.
* **Restart/checkpoint round-trip** with `MY25` selected was not exercised.
  `tke`/`gls`/`Lscale`/`Akk`/`Akv`/`Akt` are checkpointed unconditionally by
  existing code, so this should work, but it is untested.
* **`Akp`** is left at its initial value under MY2.5 and written to checkpoints
  as such. Harmless, but a restart that later switches to GLS would inherit it.

## Open questions

1. **Lateral BC type for `tke`/`gls`.** Moana's `roms.in` sets
   `LBC(isMtke) = Gra` (zero gradient) on all four edges, i.e. the
   `tkebc_im.F` gradient branch. The port calls `FillPatch(..., zvel_bc(), ...)`
   exactly as the GLS path does, which resolves to whatever the run configures
   for z-velocity rather than an unconditional zero gradient. This is inherited
   from the GLS plumbing and is fine for the Upwelling test, but the Moana
   configuration needs a `foextrap`-style BC on `tke`/`gls` specifically.
   Flagging rather than changing, since it would alter the shared BC path.
2. **Initial `tke`/`gls` seed.** ROMS `mod_mixing.F:1041-1046` seeds MY2.5 from
   `gls_Kmin`/`gls_Pmin` even though MY2.5's own floor is `my_qmin = 1e-8`.
   Reproduced as-is (`init_my25_vmix`), and `remora.gls_Kmin`/`gls_Pmin` remain
   readable under `MY25` for that reason. Worth confirming against Moana's
   `roms.in` values.
3. **Deviation 1 upstream.** Should the equivalent off-by-one be fixed in
   `REMORA_gls.cpp`? Doing so would break GLS bit-identity, so it was left
   alone here, but it looks like a genuine bug for non-periodic domains.
4. **`Akt` per-tracer.** ROMS loops `itrc=1,NAT` (temperature and salinity);
   the port loops over all `ncons`. Identical for a two-tracer setup, and
   Moana sets `Akt(temp) == Akt(salt)` anyway, but it would differ if passive
   tracers were added.
