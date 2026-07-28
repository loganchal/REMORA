# `t3dmix2_iso` port: structural self-check

Port of ROMS 3.9 rotated isoneutral harmonic tracer diffusion (`TS_DIF2` + `MIX_ISO_TS`)
into REMORA as `REMORA::t3dmix2_iso`, selected with `remora.harmonic_mixing_type = isopycnal`
(synonyms `isopycnic`, `iso`). Default (`s`) is unchanged.

## Sources

| Role | File |
| --- | --- |
| Specification (as compiled for Moana; all `#ifdef`s resolved) | `external/moana_spec/t3dmix.f90`, `t3dmix2_tile` at lines 517-800 |
| Annotated original | `external/roms-3.9/ROMS/Nonlinear/t3dmix2_iso.h` |
| Port | `Source/TimeIntegration/REMORA_t3dmix.cpp`, `REMORA::t3dmix2_iso` at lines 263-471 |

Guards that survived preprocessing in the Moana build: `MASKING` only. Off: `WET_DRY`,
`DIFF_3DCOEF`, `TS_MIX_CLIMA`, `TS_MIX_STABILITY`, `TS_MIX_MAX_SLOPE`, `TS_MIX_MIN_STRAT`,
`DIAGNOSTICS_TS`, `PROFILE`. The `.f90` and the `.h` agree term for term once those are
resolved; no discrepancy was found between spec and original.

Unused-in-Moana ROMS declarations deliberately not carried over (they are dead code with
these cppdefs): `small = 1.0e-14` (`t3dmix2_iso.h:182`, used only under `TS_MIX_MAX_SLOPE`),
`slope_max = 0.0001` (`:183`, ditto), `strat_min = 0.1` (`:184`, used only under
`TS_MIX_MIN_STRAT`). Only `eps = 0.5` (`:181` / `t3dmix.f90:565`) is live.

## Index-convention mapping

ROMS carries the vertical dependence through the two-slot `(k1,k2)` recursive blocking of
its `K_LOOP` (`t3dmix.f90:660-664`); this is a scratch-memory device, not a serial
dependence. The port stores the same quantities as full 3D arrays, exactly as the existing
`t3dmix2_geo` does. Because there is no vertical running dependence, every loop is a fully
fused `(i,j,k,n)` `ParallelFor` — GPU-safe with no k-loop-inside-column construction needed.

REMORA cells are `k = 0..N` (`AMREX_ASSERT(bx.smallEnd(2)==0 && bx.bigEnd(2)==N)` in
`lin_eos`), bottom-up; ROMS rho levels are `1..N(ng)`. So REMORA cell `c` = ROMS rho level
`c+1`, and `N(ng) = N+1`.

| ROMS | Meaning | Port |
| --- | --- | --- |
| `dRdx/dTdx(:,:,k1)` at iteration `k` | rho level `k` | `dRdx/dTdx(i,j,c)` with `c = k-1` |
| `dRdx/dTdx(:,:,k2)` at iteration `k` | rho level `k+1` | `dRdx/dTdx(i,j,c+1)` |
| `dTdr/FS(:,:,k1)` at iteration `k` | w face `k-1/2` (below rho `k`) | `dTdr/FS(i,j,c)` (face index = upper cell) |
| `dTdr/FS(:,:,k2)` at iteration `k` | w face `k+1/2` (above rho `k`) | `dTdr/FS(i,j,c+1)` |
| ROMS faces `0` and `N(ng)` (zeroed) | bottom, top | faces `0` and `N+1` |
| `IF (k.lt.N(ng))` around FS-final | interior faces only | `grow(zbx,IntVect(0,0,-1))` = `k = 1..N` |

This is the same face convention `t3dmix2_geo` uses for `dTdz` (face `k` = lower face of
cell `k`; `dTdz(i,j,0)` and `dTdz(i,j,N+1)` set to zero).

## Density field: `pden`, not `rho`

Traced in the spec: `t3dmix.f90:506` passes `OCEAN(ng) % pden` (not `OCEAN(ng) % rho`) —
i.e. **surface-referenced potential density anomaly**, not in-situ density. Set in
`ROMS/Utility/../Nonlinear/rho_eos.F`:

* nonlinear EOS — `rho_eos.F:479-482`: `rho = den` (in-situ), `pden = (den1 - 1000)*rmask`
  (`den1` = zero-pressure density, i.e. before the secant bulk modulus compression).
* linear EOS — `rho_eos.F:722-743`: `pden = rho`, identical.

REMORA computed only the in-situ `rho`, so `pden` was plumbed through the same path as
`rho`:

* `REMORA::rho_eos` / `lin_eos` / `nonlin_eos` gain a `pden` output argument
  (`Source/REMORA.H`, `Source/TimeIntegration/REMORA_rho_eos.cpp`).
  `REMORA_rho_eos.cpp:107-109` mirrors `rho_eos.F:743`; `REMORA_rho_eos.cpp:259-262` mirrors
  `rho_eos.F:479-482`.
* `Source/TimeIntegration/REMORA_setup_step.cpp:54-57` allocates a local `mf_pden` with the
  same `BoxArray`/`DistributionMapping`/ghost width as the existing local `mf_rho`, **only
  when isopycnal mixing is selected**; otherwise an empty `Array4` is passed and the EOS
  kernels skip the store (`calc_pden`). No extra memory or work in the default path.
* `rho_eos` is already called on `gbx2` (grown by `NGROW`), so `pden` has the `i±1, j±1`
  halo the kernel needs.
* `pden` is derived from `state_old` (`nrhs`), matching ROMS, where `pden` comes from
  `rho_eos` on `t(:,:,:,nrhs,:)` while `t3dmix` updates `nnew`.

## Line-by-line flux comparison

`SPEC` = `moana_spec/t3dmix.f90`; `H` = `roms-3.9/.../t3dmix2_iso.h`;
`PORT` = `Source/TimeIntegration/REMORA_t3dmix.cpp`.

### 1. Horizontal gradients at u points

* SPEC:668-676 / H:220-246
  ```
  cff=0.5*(pm(i,j)+pm(i-1,j)); cff=cff*umask(i,j)
  dRdx(i,j,k2)=cff*(pden(i,j,k+1)-pden(i-1,j,k+1))
  dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-t(i-1,j,k+1,nrhs,itrc))
  ```
* PORT:349-354 — same, with `k+1 -> k` under the cell-index shift, `umask -> msku(i,j,0)`,
  `t(...,nrhs,...) -> state_rhs`. Mask is applied to `cff` (so to both `dRdx` and `dTdx`),
  as a separate statement, exactly as in ROMS — not folded into the `0.5*(pm+pm)` product,
  which would change rounding. Loop box `xbx` = `bx` node-centred in x, matching ROMS
  `j=Jstr..Jend, i=Istr..Iend+1`.

### 2. Horizontal gradients at v points

* SPEC:682-690 / H:252-278; PORT:361-366. Same substitutions, `vmask -> mskv`, loop box
  `ybx` matching ROMS `j=Jstr..Jend+1, i=Istr..Iend`.

### 3. Boundary faces

* SPEC:695-701 / H:283-289: `dTdr(i,j,k2)=0`, `FS(i,j,k2)=0` for `k=0` and `k=N(ng)`, over
  `i=Istr-1..Iend+1, j=Jstr-1..Jend+1`.
* PORT:371-377: same values at faces `0` and `N+1` over `grow(bx,(1,1,0))` slabbed in z.

### 4. Interior faces — the stratification clip

* SPEC:706-712 / H:308-309, 327-330
  ```
  cff1=MAX(pden(i,j,k)-pden(i,j,k+1),eps)
  cff=-1.0/cff1
  dTdr(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc))
  FS(i,j,k2)=cff*(z_r(i,j,k+1)-z_r(i,j,k))
  ```
* PORT:386-392
  ```
  Real cff1 = std::max(pden(i,j,k-1)-pden(i,j,k), eps);
  Real cff = Real(-1.0)/cff1;
  dTdr(i,j,k,n) = cff * (state_rhs(i,j,k,n) - state_rhs(i,j,k-1,n));
  FS(i,j,k,n)   = cff * (z_r(i,j,k) - z_r(i,j,k-1));
  ```
  Face `k` has lower cell `k-1` and upper cell `k`, so ROMS `(k, k+1)` becomes `(k-1, k)`.
  The clip is the plain `eps = 0.5` kg/m³ branch (H:307-310), i.e. the `#else` of
  `TS_MIX_MAX_SLOPE` / `TS_MIX_MIN_STRAT`, which is what the Moana build compiled
  (SPEC:706-707 has no alternative branch). `-1.0/cff1` is kept as a division, not
  rewritten as a reciprocal multiply.

### 5. Rotated x flux `FX`

* SPEC:724-735 / H:342-356
  ```
  cff=0.25*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*on_u(i,j)
  FX(i,j)=cff*(Hz(i,j,k)+Hz(i-1,j,k))*
          (dTdx(i,j,k1)-
           0.5*(MAX(dRdx(i,j,k1),0)*(dTdr(i-1,j,k1)+dTdr(i,j,k2))+
                MIN(dRdx(i,j,k1),0)*(dTdr(i-1,j,k2)+dTdr(i,j,k1))))
  ```
* PORT:403-413 — identical, with `k1 -> k`, `k2 -> k+1` for `dTdr`. The `MAX` term is kept
  first in the sum (ROMS order; the geopotential variant has `MIN` first because its slope
  has the opposite sign convention). `on_u` is expanded to `two/(pn(i-1,j,0)+pn(i,j,0))`,
  which is verbatim ROMS `metrics.F:411`, including operand order. Curvilinear metrics
  therefore enter exactly as ROMS precomputes them.

### 6. Rotated y flux `FE`

* SPEC:741-752 / H:362-376; PORT:420-430. Same, with `om_v = two/(pm(i,j-1,0)+pm(i,j,0))`
  from `metrics.F:447`.

### 7. Vertical cross-term `FS`

* SPEC:758-776 / H:382-402
  ```
  cff1=MAX(dRdx(i  ,j,k1),0); cff2=MAX(dRdx(i+1,j,k2),0)
  cff3=MIN(dRdx(i  ,j,k2),0); cff4=MIN(dRdx(i+1,j,k1),0)
  cff = cff1*(cff1*dTdr(i,j,k2)-dTdx(i  ,j,k1)) + cff2*(...) + cff3*(...) + cff4*(...)
  ... same four terms for dRde/dTde, accumulated as cff=cff+...
  FS(i,j,k2)=0.5*cff*diff2(i,j,itrc)*FS(i,j,k2)
  ```
* PORT:439-457 — identical with `k1 -> k-1`, `k2 -> k` for the cell-level fields and
  `dTdr(:,:,k2) -> dTdr(i,j,k)`. The x-part is accumulated into a local `cff`, then the
  y-part is added with `cff = cff + ...`, then `FS` is scaled — i.e. ROMS's association is
  preserved (the geopotential port instead writes `FS` twice, which associates differently;
  that difference is intentional and follows each ROMS variant). `FS` is read-modify-written
  in a separate kernel from the one that filled it, so the kernel boundary provides the
  ordering; each thread touches only its own `(i,j,k,n)`.
  ROMS restricts this to `1 <= k <= N(ng)-1`, i.e. interior faces; the port's loop box
  `grow(zbx,IntVect(0,0,-1))` is `k = 1..N`, the same set. Faces `0`/`N+1` keep the zero
  from step 3. Note the halo values of `FS` (from step 4, `i,j` in `Istr-1..Iend+1`) are
  left un-scaled in both ROMS and the port, and are never read.

### 8. Time step

* SPEC:786-791 / H:412-417
  ```
  cff=dt(ng)*pm(i,j)*pn(i,j)
  cff1=cff*(FX(i+1,j)-FX(i,j)); cff2=cff*(FE(i,j+1)-FE(i,j))
  cff3=dt(ng)*(FS(i,j,k2)-FS(i,j,k1)); cff4=cff1+cff2+cff3
  t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff4
  ```
* PORT:464-469 — identical, `FS(k2)-FS(k1) -> FS(i,j,k+1)-FS(i,j,k)`, `t(...,nnew,...) ->
  state`. The four `cff*` temporaries are retained so the sum associates as
  `(cff1+cff2)+cff3`, as in ROMS.

## Masking placement

ROMS applies `umask`/`vmask` only to the `cff` metric factor of the horizontal gradients
(H:221-223, 253-255) — nowhere else. The port does the same. In particular the mask is
**not** applied to `FX`/`FE` after the fact (which is what `t3dmix2_s` does, following its
own ROMS variant), and `dTdr`/`FS` carry no mask. `WET_DRY` masks (H:224-226, 256-258) are
absent from the preprocessed spec and are not ported.

## Boundary/edge handling

Loop extents follow the ROMS tile bounds one-for-one: `xbx`/`ybx` are the x/y node-centred
versions of `bx` (= `Istr..Iend+1` / `Jstr..Jend+1`); `dTdr`/`FS` are computed over `bx`
grown by one cell in x and y (= `Istr-1..Iend+1`, `Jstr-1..Jend+1`) because `FX`/`FE` read
`dTdr(i-1,·)` / `dTdr(·,j-1)`; the final update is over `bx`. Vertical boundaries are the
zeroed faces `0` and `N+1`. No physical BC is applied inside the kernel, matching ROMS
(the caller's fill-patching supplies the halo).

## Verification

1. **Build** — full CPU rebuild of `build/` (`make -j8`) is clean; no warnings or errors
   from the touched translation units.
2. **Smoke, new scheme** —
   `remora_exec Exec/Upwelling/inputs remora.max_step=20 remora.harmonic_mixing_type=isopycnal`
   runs 20 steps and writes `plt00020`.
3. **Smoke, default** — `remora_exec Exec/Upwelling/inputs remora.max_step=20` runs 20 steps
   and writes `plt00020`. Its plotfile is **bit-identical** (`diff -r`, excluding `job_info`)
   to the same run from a `moana/dev` build, confirming the `pden` plumbing does not perturb
   the default (`s`) path.
4. **The kernel actually does something** — the stock Upwelling inputs leave `tnu2 = 0`, so
   all three schemes trivially agree there. Re-running with
   `remora.horizontal_mixing_type=constant remora.tnu2_temp=25 remora.tnu2_salt=25` (the
   Moana value) gives three mutually distinct plotfiles for `s`, `geopotential` and `iso`,
   i.e. the isopycnal branch is reached and produces its own tendency.
5. **Stability / no NaN** — 100 steps with `iso` and `tnu2 = 25` complete; a raw scan of the
   `plt00100` cell data (262,400 doubles) finds no NaN or Inf, with values in a physical
   range.

Notes on coverage:

* The Upwelling case uses the **linear** EOS, so these runs exercise the kernel with
  `pden == rho`. The nonlinear-EOS `pden` path (`den1 - 1000`) is compiled and plumbed but
  is not exercised here.
* AMReX FPE trapping (`amrex.fpe_trap_invalid=1`) cannot be used as an extra check on this
  case: it fires in the pre-existing `Problem::init_analytic_bathymetry` during setup,
  before any of this code runs.

## Open questions

* No numerical parity test against ROMS output yet — this is a structural port only. A
  Moana-parity gate (single-column or small-box, comparing against a ROMS `t3dmix2_iso`
  reference) would be the natural next step.
* `diff2` is 2D-per-tracer in REMORA (matching ROMS `MIXING(ng)%diff2` with `DIFF_3DCOEF`
  off), so the `DIFF_3DCOEF` branch (`diff3d_r(i,j,k)`) is not ported. If a 3-D diffusivity
  is ever wanted, it is H:341-347, 361-367, 399-403.
