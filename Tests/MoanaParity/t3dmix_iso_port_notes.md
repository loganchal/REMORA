# Isoneutral tracer diffusion

`remora.harmonic_mixing_type = isopycnal` selects ROMS harmonic tracer diffusion along surfaces of constant potential density (`TS_DIF2` with `MIX_ISO_TS`). `isopycnic` and `iso` are aliases. The default, `s`, diffuses tracers along terrain-following surfaces.

## Density and grid indexing

The [diffusion kernel](../../Source/TimeIntegration/REMORA_t3dmix.cpp) uses the surface-referenced potential density anomaly from the [equation of state](../../Source/TimeIntegration/REMORA_rho_eos.cpp). For the nonlinear equation of state, `pden = (den1 − 1000) × rmask`, where `den1` is the density at zero pressure. For the linear equation of state, `pden = rho`. The field is computed from `state_old` at the ROMS `nrhs` time level. Isoneutral mixing allocates this field with the horizontal ghost cells needed for the adjacent-cell gradients.

Moana uses horizontal masks and a two-dimensional diffusivity for each tracer. The optional ROMS wetting/drying, three-dimensional diffusivity, climatology, stability and alternative slope-limiter branches are disabled. The active denominator limit is `eps = 0.5` kg/m³.

ROMS reuses two work-array slots, `k1` and `k2`, in a vertical loop. REMORA stores full three-dimensional work arrays so each kernel can operate on all cells and tracer components in parallel. REMORA rho cells are numbered 0 through N from the bottom, corresponding to ROMS levels 1 through N+1.

| ROMS | Meaning | Port |
| --- | --- | --- |
| `dRdx/dTdx(:,:,k1)` at iteration `k` | rho level `k` | `dRdx/dTdx(i,j,c)` with `c = k-1` |
| `dRdx/dTdx(:,:,k2)` at iteration `k` | rho level `k+1` | `dRdx/dTdx(i,j,c+1)` |
| `dTdr/FS(:,:,k1)` at iteration `k` | w face `k-1/2` (below rho `k`) | `dTdr/FS(i,j,c)` (face index = upper cell) |
| `dTdr/FS(:,:,k2)` at iteration `k` | w face `k+1/2` (above rho `k`) | `dTdr/FS(i,j,c+1)` |
| ROMS faces `0` and `N(ng)` (zeroed) | bottom, top | faces `0` and `N+1` |
| `IF (k.lt.N(ng))` around FS-final | interior faces only | `grow(zbx,IntVect(0,0,-1))` = `k = 1..N` |

## Flux calculation

The references below identify the preprocessed Moana ROMS routine (`SPEC`, `t3dmix.f90`), the original ROMS `t3dmix2_iso.h` (`H`), and the [REMORA implementation](../../Source/TimeIntegration/REMORA_t3dmix.cpp) (`PORT`).

### 1. Horizontal gradients at u points

* SPEC:668-676 / H:220-246
  ```
  cff=0.5*(pm(i,j)+pm(i-1,j)); cff=cff*umask(i,j)
  dRdx(i,j,k2)=cff*(pden(i,j,k+1)-pden(i-1,j,k+1))
  dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)-t(i-1,j,k+1,nrhs,itrc))
  ```
* PORT:349-354 — same, with `k+1 -> k` under the cell-index shift, `umask -> msku(i,j,0)`,
  `t(...,nrhs,...) -> state_rhs`. Mask is applied to `cff` (so to both `dRdx` and `dTdx`),
  as a separate statement, in ROMS order. Folding it into the `0.5*(pm+pm)` product
  would change rounding. Loop box `xbx` = `bx` node-centred in x, matching ROMS
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
* PORT:403-413 — with `k1 -> k`, `k2 -> k+1` for `dTdr`. The `MAX` term is kept
  first in the sum (ROMS order; the geopotential variant has `MIN` first because its slope
  has the opposite sign convention). `on_u` is expanded to `two/(pn(i-1,j,0)+pn(i,j,0))`,
  from ROMS `metrics.F:411`, including operand order. Curvilinear metrics
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
* PORT:439-457 — with `k1 -> k-1`, `k2 -> k` for the cell-level fields and
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
* PORT:464-469 — `FS(k2)-FS(k1) -> FS(i,j,k+1)-FS(i,j,k)`, `t(...,nnew,...) ->
  state`. The four `cff*` temporaries are retained so the sum associates as
  `(cff1+cff2)+cff3`, as in ROMS.

The u- and v-face masks enter the horizontal gradient metric. The fluxes use those masked gradients without a second mask. The face derivative and vertical flux extend one cell beyond the horizontal tile so adjacent faces can read them. The caller supplies the halo values. The divergence updates only the original cell box.

## Component tests

The CPU Upwelling case uses a linear equation of state.

| Test | Result |
|---|---|
| 20 steps with `remora.harmonic_mixing_type=isopycnal` | Completed and wrote `plt00020`. |
| 20 steps with the default scheme | The plotfile data matched the baseline build bit for bit. Run metadata was excluded. |
| Constant horizontal mixing with `remora.tnu2_temp=25` and `remora.tnu2_salt=25` m²/s | The terrain-following, geopotential and isoneutral runs produced distinct fields. The default Upwelling inputs have zero diffusivity. |
| 100 steps with the isoneutral scheme at 25 m²/s | All 262,400 double-precision values examined in the output were finite. |

The component test exercised `pden = rho`. It did not test the nonlinear equation of state or compare individual fluxes numerically with a ROMS calculation. Enabling `amrex.fpe_trap_invalid` in this case caused an exception during analytic bathymetry initialisation, before the diffusion kernel ran.
