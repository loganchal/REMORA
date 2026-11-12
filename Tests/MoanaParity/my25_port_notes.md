# Mellor–Yamada level 2.5 vertical mixing

Select `remora.vertical_mixing_type = MY25` to use the ROMS 3.9 closure with Kantha–Clayson stability functions. The [implementation](../../Source/TimeIntegration/REMORA_my25.cpp) advances `tke` (q²) and `gls` (q²ℓ) using the Moana configuration of `my25_prestep.F` and `my25_corstep.F`. ROMS defines the closure constants in `mod_scalars.F`, initial values in `mod_mixing.F` and lateral boundary conditions in `tkebc_im.F`.

## ROMS options

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
| `PERFECT_RESTART` | off | ROMS history/restart initialisation reseeds the closure variables. Native REMORA checkpoints store the evolving turbulence fields. |

## Closure constants

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
| `my_B1p2o3` | `B1^(2/3)` | ≈6.50737 | :3761 | surface/bottom `q²` |
| `my_B1pm1o3` | `1/B1^(1/3)` | ≈0.39201 | :3762 | KC `Sm` numerator |
| `my_E1o2` | `0.5*E1` | 0.9 | :3763 | **unused** in these two routines |
| `my_Sm1` | `A1*A2*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1))` | — | :3764-3766 | **unused** (Canuto path only) |
| `my_Sm2` | `9*A1*A2` | 6.1272 | :3767 | |
| `my_Sh1` | `A2*(1-6*A1/B1)` | ≈0.49393 | :3768 | |
| `my_Sh2` | `3*A2*(6*A1 + B2*(1-C3))` | 30.192 | :3770 | **KANTHA_CLAYSON** form |
| `my_Sm3` | `A1*(1-3*C1-6*A1/B1)` | — | :3774 | **not compiled** (KC branch) |
| `my_Sm4` | `18*A1² + 9*A1*A2*(1-C2)` | 17.07336 | :3771 | **KANTHA_CLAYSON** form |

Local literals: `Gamma = 1/6` (`my25_prestep.F:139`), `Gadv = 1/3`
(`my25_corstep.F:189`), `eps = 1e-10` (`my25_corstep.F:190`),
wall prefactor `my_E2/vonKar²` (`my25_corstep.F:598`), buoyancy dead-band
`-5.0e-5` (`my25_corstep.F:604-608`).

The Kantha–Clayson stability functions use the expressions for `my_Sh2` and `my_Sm4` above and `my_B1pm1o3` in the numerator of `Sm`.

## Index mapping

`N` in REMORA is the top rho index, so there are `N+1` rho levels and `N+2`
z-nodes.

```
N_ROMS            = N + 1
w-point  index    k_REMORA = k_ROMS            (0 = bottom, N+1 = surface)
rho-cell index    k_REMORA = k_ROMS - 1        =>  Hz_ROMS(k) == Hz(k-1)
```

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

## ROMS routines and REMORA kernels

The following line numbers refer to the original [MY2.5 implementation](https://github.com/loganchal/REMORA/blob/916d7e9/Source/TimeIntegration/REMORA_my25.cpp). The ROMS prestep uses the same advection equations as GLS. The corrector implements the MY2.5 production, dissipation and stability functions.

| ROMS `my25_prestep.F` | `REMORA_my25.cpp` | Notes |
| --- | --- | --- |
| :208-219 | :121-127 | masked `grad`/`gradL` on u-points |
| :222-237 | :129-141 | west/east one-sided gradient copy |
| :239-252 | :142-148 | `XF`, `FX`, `FXL` |
| :253-291 | :151-178 | same construction in η |
| :293-305 | :185-197 | AM3 weights `cff1/cff2/cff3`, `indx` |
| :306-325 | :200-212 | `Hz_half`, intermediate time-level update, store `nnew = Hz*nstp` |
| :327-374 | :214-243 | vertical advective flux (4th-order interior, one-sided ends) |
| :376-393 | :245-256 | time-step vertical advection, divide by `Hz_half` |
| :396-411 | :259-263 | `tkebc_tile` + `exchange_w3d_tile` → `FillPatch` |

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
| :505-521 | :576-587 | time-step horizontal advection — **no floor** |
| :523-566 | :589-613 | vertical advective flux |
| :568-578 | :615-621 | time-step vertical advection — **no floor** |
| :585-593 | :624-637 | `FCK = -0.5*dt*(Akk(k)+Akk(k-1))/Hz(k)`, `CF = 0` — **not** zeroed at the ends |
| :601-612 | :644-649 | buoyancy dead-band, `Qprod` |
| :613-617 | :650-651 | `Ls_unlmt` from the old time step |
| :618-625 | :653-656 | production into `tke` (×2) and `gls` (×`E1*Ls_unlmt`) |
| :626-636 | :658-667 | `Qdiss`, wall function `Wscale`, `BCK`/`BCP` |
| :642-653 | :677-685 | Dirichlet surface/bottom: `q² = B1^(2/3)·½|τ|`, `q²l = 0` |
| :655-675 | :687-699 | `tke` tridiagonal |
| :677-697 | :701-713 | `gls` tridiagonal — diagonal `BCP`, off-diagonal **`FCK`** |
| :706-716 | :722-728 | `my_qmin` floor, `Ls_unlmt`, Galperin length limit |
| :717-728 | :730-734 | `Gh`, `Sh`, KC `Sm` |
| :733-742 | :736-743 | two-time-level `ql` average, `Akv`, `Akt` |
| :744-748 | :744-745 | `Akk` |
| :750-752 | :747-748 | `Lscale` |
| :755-772 | :755-763 | `tkebc_tile` + exchanges |

## Boundary conditions and vertical solves

* At non-periodic west and south edges the prestep takes one-sided gradients at `dlo.x` and `dlo.y`. The copied gradients retain the u- or v-face mask. The high-side x face is `dhi.x+1`.
* The corrector applies the `my_qmin` floor after the tridiagonal solves. Both equations use `Akk` and the same off-diagonal coefficients (`FCK`), with separate diagonals (`BCK`, `BCP`). The end coefficients carry the prescribed surface and bottom values into the first and last interior equations. Boundary q² depends on stress magnitude, q²ℓ is zero, and the closure does not require bottom roughness.
* Each GPU thread solves one vertical column. The spline shear calculation and both tridiagonal sweeps preserve the vertical order. The q² and q²ℓ sweeps share column-local scratch storage. Separate kernels handle pointwise updates.
* [Initialisation](../../Source/Initialization/REMORA_init.cpp) sets `tke` and `gls` to `gls_Kmin` and `gls_Pmin`, respectively. It sets `Akv`, `Akt` for every tracer, and `Akk` to zero at the surface and bottom. The corrector updates interior faces. `Akp` is unused by MY2.5.
* The lateral fills for both turbulence fields use `tke_bc()`. The Moana case specifies `remora.bc.tke.type = outflow` on each side for its zero-gradient turbulence boundary. The tracer diffusivity update applies to all `ncons` components; the ROMS Moana case uses temperature and salinity.

## Component tests

The original tests used an AppleClang Release build on one CPU rank.

| Test | Result |
|---|---|
| Upwelling, 100 and 500 steps | Both runs completed. `amrex_fnan` found no NaNs in temperature, salinity or the three velocity components in either plotfile. |
| Mixing coefficients | Maximum `Akv` was 2.0 × 10⁻³ at step 120 and 2.7 × 10⁻² at step 500. Maximum `Akt` and `Akk` were 3.4 × 10⁻² and 1.2 × 10⁻² m²/s. All three were zero on the end faces and at least their background values in the interior. |
| Turbulence variables | `Lscale` ranged from 0 to 2.23 m over 500 steps. Maximum `tke` was 1.3 × 10⁻³ m²/s² and maximum `gls` was 2.3 × 10⁻³. |
| Step 100, minimum u and maximum temperature | MY2.5: −0.109 and 21.728. GLS: −0.115 and 21.736. Constant viscosity: −0.033 and 21.931. |
| Existing analytic and GLS cases | All five plotfile variables matched the baseline build bit for bit. |

With invalid floating-point exceptions trapped, 15 of 16 CTest cases failed during analytic bathymetry initialisation on both the baseline (`8f95ac3`) and the MY2.5 build. With trapping disabled, Upwelling completed. Absolute differences from the stored reference were 3.7 × 10⁻⁸ for temperature (1.7 × 10⁻⁹ relative), 7.8 × 10⁻¹⁴ for salinity, 2.1 × 10⁻⁹ for u and 1.8 × 10⁻⁸ for v.

These component tests did not compare the individual closure terms with ROMS, exercise a checkpoint restart, or test multiple ranks or GPUs. The [5 km model comparisons](../../Exec/Moana/model-notes.rst) cover the integrated GPU configuration.
