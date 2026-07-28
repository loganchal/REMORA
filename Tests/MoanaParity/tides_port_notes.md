# Barotropic tidal boundary forcing: ROMS → REMORA port notes

Port of the Moana hindcast's tidal forcing: ROMS 3.9 `SSH_TIDES` + `UV_TIDES` with
`ADD_FSOBC` + `ADD_M2OBC` (and `RAMP_TIDES` **undefined**), i.e. exactly what
`external/moana_hindcast/configuration/roms3d.h` lines 78–82 compile.

Oracle: `external/moana_spec/set_tides.f90` — the preprocessed Moana build of
`ROMS/Nonlinear/set_tides.F`. Line numbers below refer to that preprocessed file.
Supporting ROMS sources: `ROMS/Modules/mod_tides.F` (storage),
`ROMS/Nonlinear/get_idata.F` (one-time read), `ROMS/External/varinfo.dat` (units),
`ROMS/Utility/get_grid.F` (grid angle), `ROMS/Nonlinear/main3d.F` (call site).

Everything is gated on `remora.tides` (default `false`). With the gate off the code
path is bit-identical to `moana/dev` (see *Verification* below).

---

## 1. New inputs

| input | default | meaning | ROMS equivalent |
|---|---|---|---|
| `remora.tides` | `false` | enable barotropic tidal boundary forcing | `SSH_TIDES`+`UV_TIDES`+`ADD_FSOBC`+`ADD_M2OBC` |
| `remora.tide_start` | `0.0` | reference time origin for tidal phase, **days** | `TIDE_START` (Moana: `0.0d0`) |
| `remora.nc_tide_file` | `""` | tidal harmonic constants file | `TIDENAME` |

Requirements enforced at initialization (`Source/REMORA.cpp`, in `init_only`):
`remora.tides=true` needs a NetCDF build, `remora.ic_bc_type=netcdf`,
NetCDF open-boundary data (`remora.nc_bdry_file`), a tide file, and a single level.

## 2. Files

| file | contents |
|---|---|
| `Source/IO/REMORA_NCTideData.{H,cpp}` | `NCTideData`: static (read-once) tide-file reader; also `REMORA::init_angler_from_netcdf` (grid `angle`) |
| `Source/TimeIntegration/REMORA_set_tides.cpp` | `REMORA::set_tides(lev,time)` — the harmonic evaluation |
| `Source/BoundaryConditions/REMORA_BoundaryConditions_netcdf.cpp` | the ADD_FSOBC / ADD_M2OBC additions inside `fill_from_bdyfiles` |
| `Source/TimeIntegration/REMORA_Advance.cpp` | call site (once per baroclinic step) |
| `Source/Initialization/REMORA_make_new_level.cpp` | allocation of `vec_angler`, `vec_Etide`, `vec_Utide`, `vec_Vtide` |
| `Source/REMORA_DataStruct.H` | `use_tides`, `tide_start` inputs |
| `Tests/MoanaParity/make_tide_test_files.py` | synthetic ROMS-format grid/init/bry/tide files for the functional test |
| `Tests/MoanaParity/inputs.tidebox` | the functional test case |
| `Tests/MoanaParity/check_tide_boundaries.py` | independent numpy re-implementation of `set_tides.F`, compared against model output |

## 3. ROMS line ranges ↔ ported code

| ROMS `set_tides.f90` | what | REMORA |
|---|---|---|
| 484–492 | `TIDES(ng)%SSH_Tamp/SSH_Tphase/UV_T*/Tperiod` | `NCTideData` members (`Eamp/Ephase/Cangle/Cphase/Cmax/Cmin`, `period_h/period_d`) |
| 478 (`GRID(ng)%angler`) | grid rotation angle | `vec_angler`, read by `init_angler_from_netcdf` |
| 645 `LprocessTides` | master switch | `solverChoice.use_tides` |
| 650 | `ramp = 1` (RAMP_TIDES undefined) | `const Real ramp = one;` in `set_tides` |
| 657–673 | `Etide` accumulation + `rmask` | first `ParallelFor` in `set_tides` |
| 679–691 | `CLIMA%ssh += Etide` (LsshCLM) | **not ported** — see §7 |
| 704–770 | ADD_FSOBC: `zeta_west/east/south/north += 0.5*(Etide,Etide)` | `tide_zeta_val` in each of the four edge kernels of `fill_from_bdyfiles` |
| 779–818 | `Uwrk/Vwrk` ellipse, averaging to u/v points, `umask/vmask` | `tide_ellipse()` + second/third `ParallelFor` in `set_tides` |
| 824–844 | `CLIMA%ubarclm/vbarclm += Utide/Vtide` (Lm2CLM) | **not ported** — see §7 |
| 849–931 | ADD_M2OBC: `ubar_*/vbar_* += Utide/Vtide` | `tide_val` in each of the four edge kernels of `fill_from_bdyfiles` |

The mathematics, verbatim:

```
omega_c = 2*pi*(time - tide_start*86400) / Tperiod(c)

Etide(i,j) = sum_c ramp * SSH_Tamp(i,j,c) * cos(omega_c - SSH_Tphase(i,j,c))     [* rmask]

angle   = UV_Tangle(i,j,c) - angler(i,j)
phase   = omega_c - UV_Tphase(i,j,c)
Uwrk    = UV_Tmajor*cos(angle)*cos(phase) - UV_Tminor*sin(angle)*sin(phase)
Vwrk    = UV_Tmajor*sin(angle)*cos(phase) + UV_Tminor*cos(angle)*sin(phase)
Utide(i,j) = sum_c ramp*0.5*(Uwrk(i-1,j) + Uwrk(i,j))                            [* umask]
Vtide(i,j) = sum_c ramp*0.5*(Vwrk(i,j-1) + Vwrk(i,j))                            [* vmask]
```

Two deliberate implementation choices, both algebraically identical to the Fortran:

* `Uwrk`/`Vwrk` are **recomputed** at the two rho points a u (v) point needs instead of
  being stored in a tile-sized scratch array. This removes a scratch MultiFab and any
  inter-thread dependence; the cost is one extra pair of sin/cos per point per
  constituent on a 2D field evaluated once per baroclinic step.
* The masks are applied inside the constituent loop exactly as ROMS does (`mask` is
  0 or 1, so this is also bit-identical to masking once at the end).

## 4. Unit conventions (from `ROMS/External/varinfo.dat`)

ROMS applies the `Fscale` column when it reads each variable (`get_idata.F` →
`get_2dfld`/`get_ngfld` → `nf_fread*`). The port applies exactly the same factors in
`NCTideData::Initialize`.

| variable | file units | Fscale | stored units |
|---|---|---|---|
| `tide_period` | hours | `3600.0` | seconds |
| `tide_Eamp` | meter | `1.0` | m |
| `tide_Ephase` | degrees | `0.017453292519943295` | radians |
| `tide_Cangle` | degrees | `0.017453292519943295` | radians |
| `tide_Cphase` | degrees | `0.017453292519943295` | radians |
| `tide_Cmax` | m s⁻¹ | `1.0` | m/s |
| `tide_Cmin` | m s⁻¹ | `1.0` | m/s |
| `angle` (grid file) | radians | — (no scale; `get_grid.F`) | radians |

The degrees→radians factor is coded as the literal `0.017453292519943295` that ROMS
reads from `varinfo.dat`, not as a recomputed `pi/180`.

Constituents with `Tperiod <= 0` are skipped, as in ROMS (`set_tides.f90:660, 783`).

## 5. Index map and the boundary additions

ROMS rho index `i` runs `0..Lm+1` (0 and `Lm+1` are the ghost ring); ROMS u index `i`
runs `1..Lm+1`. REMORA's cell index is `ROMS_rho_index - 1` (`-1..Nx`) and REMORA's
x-face index is `ROMS_u_index - 1` (`0..Nx`). **Both shift by the same amount**, so
every ROMS stencil keeps its algebraic form; only the loop bounds change.

`dlo`/`dhi` below are the bounds of the **cell-centered** domain (`0..Nx-1`, `0..Ny-1`).

| ROMS boundary array | ROMS index (set_tides.f90) | REMORA index | REMORA array holding it |
|---|---|---|---|
| `zeta_west(j)`  | `0.5*(Etide(Istr-1,j)+Etide(Istr,j))`   (712–713) | `0.5*(E(dlo.x-1,j)+E(dlo.x,j))`   | zeta bdry column `i=dlo.x-1` |
| `zeta_east(j)`  | `0.5*(Etide(Iend,j)+Etide(Iend+1,j))`   (729–730) | `0.5*(E(dhi.x,j)+E(dhi.x+1,j))`   | zeta bdry column `i=dhi.x+1` |
| `zeta_south(i)` | `0.5*(Etide(i,Jstr-1)+Etide(i,Jstr))`   (746–747) | `0.5*(E(i,dlo.y-1)+E(i,dlo.y))`   | zeta bdry row `j=dlo.y-1` |
| `zeta_north(i)` | `0.5*(Etide(i,Jend)+Etide(i,Jend+1))`   (763–764) | `0.5*(E(i,dhi.y)+E(i,dhi.y+1))`   | zeta bdry row `j=dhi.y+1` |
| `ubar_west(j)`  | `Utide(Istr,j)`    (856) | `Utide(dlo.x,j)`   | u face at the west boundary |
| `ubar_east(j)`  | `Utide(Iend+1,j)`  (877) | `Utide(dhi.x+1,j)` | u face at the east boundary |
| `ubar_south(i)` | `Utide(i,Jstr-1)`  (898) | `Utide(i,dlo.y-1)` | u row one below the domain |
| `ubar_north(i)` | `Utide(i,Jend+1)`  (919) | `Utide(i,dhi.y+1)` | u row one above the domain |
| `vbar_west(j)`  | `Vtide(Istr-1,j)`  (862) | `Vtide(dlo.x-1,j)` | v column one west of the domain |
| `vbar_east(j)`  | `Vtide(Iend+1,j)`  (883) | `Vtide(dhi.x+1,j)` | v column one east of the domain |
| `vbar_south(i)` | `Vtide(i,Jstr)`    (904) | `Vtide(i,dlo.y)`   | v face at the south boundary |
| `vbar_north(i)` | `Vtide(i,Jend+1)`  (925) | `Vtide(i,dhi.y+1)` | v face at the north boundary |

These indices coincide exactly with the locations at which REMORA's
`NCTimeSeriesBoundary` stores its interpolated boundary FABs (see the box definitions
in `REMORA_NCTimeSeriesBoundary::Initialize`), which is why the addition can be done
where the boundary value is consumed.

**Where the addition happens.** Rather than mutating the stored interpolated boundary
data, the tide is added to the local `bry_val` inside `fill_from_bdyfiles`. This is
equivalent to the ROMS "add to `BOUNDARY(ng)%*`" step and is idempotent: REMORA calls
`fill_from_bdyfiles` several times per step (and per fast sub-step), so accumulating
into the stored arrays would add the tide many times.

**Flather.** ROMS's Flather condition for `ubar`/`vbar` uses `BOUNDARY%zeta_*`, which
carries the tide under `ADD_FSOBC`; so the same zeta tide is added to `bry_val_zeta` in
the Flather branch. With Moana's actual configuration (`LBC(isFsur)=Cha`,
`LBC(isUbar)=LBC(isVbar)=Fla`) this is the *only* route by which the tidal elevation
enters the solution, since the Chapman condition ignores the boundary value.

## 6. Time base

ROMS `main3d.F`: `iic` and `time(ng)` are advanced at the top of the step
(lines 518–522), then `get_data`/`set_data` interpolate the open-boundary data to
`time(ng)`, then `set_tides` is called (line 624), then `step2d`. Because `initial.F`
sets `time(ng) = time(ng) - dt(ng)` before the loop starts (line 741), `time(ng)`
during a step is the time of the *current* state, i.e. the start of the step being
computed. The tidal signal is therefore frozen across the barotropic sub-steps.

REMORA: `set_tides(lev, time)` is called at the top of `REMORA::Advance`, whose `time`
argument is `t_old[lev]` — precisely the time to which `advance_2d` interpolates the
boundary data (`FillPatch(lev, t_old[lev], ...)` for `ubar`, `vbar`, `zeta`). So the
same "evaluate once per baroclinic step, at the step's start time, hold across the
fast loop" semantics hold.

The absolute origin: REMORA converts the boundary file's `ocean_time` (days) to
seconds and asserts the model time lies inside it, so `t_old` is on the same base as
the forcing files' `ocean_time`, which is ROMS's `time(ng)` base. `tide_start` (days)
is subtracted from it exactly as ROMS does. For Moana, `TIDE_START = 0.0`, so the
tidal phases are referenced to the `ocean_time` epoch of the run.

At initialization, `Etide/Utide/Vtide` are zero, so the boundary fill done by
`init_zeta_from_netcdf` carries no tide — matching ROMS, where `ini_zeta` runs before
the first `set_tides` call.

## 7. Deliberately not ported

* **`CLIMA` blocks** (`set_tides.f90:679–691` and `824–844`): these add `Etide`/
  `Utide`/`Vtide` to `CLIMA(ng)%ssh`, `%ubarclm`, `%vbarclm` for 2D climatology
  nudging, and are guarded by `LsshCLM(ng)` / `Lm2CLM(ng)`. The Moana hindcast has
  no free-surface or 2D-momentum climatology nudging (`Lm2CLM == F`), so those blocks
  are inert in the oracle build and are not ported. If 2D climatology nudging is ever
  switched on together with tides, they must be added (`REMORA_apply_clim_nudg.cpp`
  is where they would go).
* **`RAMP_TIDES`**: undefined in `roms3d.h`, so `ramp = 1` is hard-coded (with a
  comment). If a ramp is ever wanted, it is ROMS's
  `ramp = TANH((tdays(ng)-dstart)/1.0_r8)`.
* **Detiding / harmonic accumulation** (`AVERAGES_DETIDE`, the `CosW/SinW/...`
  machinery in `mod_tides.F`): not compiled in the Moana build.

## 8. Verification

All runs used `build_tides/` (CPU, Release, MPI + PnetCDF) in the `moana/tides`
worktree; `remora_exec_base` is the same configuration built from `moana/dev`.

**(a) Clean build.** `cmake` + `make -j8` with `REMORA_ENABLE_PNETCDF=ON`,
`REMORA_ENABLE_MPI=ON`: no warnings/errors from the new files.

**(b) Default-off bit-identity** (`remora.tides` absent/false):

| case | result |
|---|---|
| `Exec/Upwelling/inputs`, 10 steps, plotfiles at 5 and 10 | plotfiles byte-identical to the `moana/dev` build (only `job_info` differs: timings, path, git hash, and the two new default-valued inputs) |
| `Tests/MoanaParity/inputs.tidebox` with `remora.tides=false` (NetCDF grid + boundary path, clamped BCs), 40 steps | plotfiles byte-identical to the `moana/dev` build |

**(c) Functional test.** `make_tide_test_files.py` writes a 10×16×4, 100 m deep box
with **identically zero** boundary data, so everything seen at the open boundaries is
the tidal addition. `check_tide_boundaries.py` re-implements `set_tides.F` in numpy
(independently, from the Fortran) and compares against the model's own output.

| configuration | result |
|---|---|
| 1 constituent, T=12 h, `Eamp`=1 m, `Ephase`=90° (so `zeta_bry = sin(2πt/T)`), clamped BCs, 288 steps = 1 period | boundary `zeta` on all four edges matches the analytic tide to **4.4e-16** |
| 1 constituent, `Cmax`=0.5 m/s, `Cangle`=0, grid `angle`=0 | `ubar`/`vbar` on all 8 boundary arrays match to **≤1.7e-16** |
| same but grid `angle`=30° (ellipse rotated by `-angler`: `u` amp 0.43301, `v` amp −0.25) | all 8 match to **≤2.3e-16** |
| 2 constituents (12 h + 6 h), `tide_start`=0.25 d, 4 MPI ranks, `max_grid_size=8` | boundary `zeta` matches the analytic superposition to **4.4e-16** |
| 2 constituents (12 h + 6.2 h), **spatially varying** `Eamp/Ephase/Cmax/Cmin/Cangle/Cphase`, grid `angle`=25°, `tide_start`=0.3 d | all 12 boundary arrays match the numpy ROMS oracle to **≤6.7e-16**, serial and with 4 MPI ranks / `max_grid_size=6` |

The spatially varying case is the one that actually exercises the staggering: with
uniform fields the `0.5*(E(i-1)+E(i))` and `0.5*(Uwrk(i-1)+Uwrk(i))` averages are
indistinguishable from point values.

**(d) Moana boundary configuration.** The same box run with `zeta`=Chapman and
`ubar`/`vbar`=Flather (Moana's `roms.in` setting) is stable and fills/empties the basin
at the forcing period with interior amplitude 1.02 m for a 1 m boundary tide; with
`remora.tides=false` the same run stays exactly at rest (`max|zeta| = 0`). Since
Chapman ignores the boundary elevation, this specifically exercises the
`bry_val_zeta` addition in the Flather branch.

**What is not yet verified**

* No comparison against an actual ROMS run (no Moana tide file or ROMS executable
  available locally); the oracle here is the Fortran source re-implemented in numpy.
* Land masking of `Etide/Utide/Vtide` is ported verbatim but untested — the synthetic
  grid is all water. (Note that in ROMS a boundary point next to land gets *half* the
  tidal elevation, because `Etide` is masked before the two-point average.)
* GPU: the code is written for `ParallelFor` with no serial dependence and no nested
  extended lambdas, but only a CPU build was compiled here (no CUDA toolchain on this
  machine).
* AMR: tides are restricted to single-level runs (`max_level > 0` aborts). ROMS itself
  only processes tides on the coarsest grid in refinement applications.
* Restart: the tide file and grid angle are re-read on restart (the checkpoint path
  goes through `MakeNewLevelFromScratch` → `init_only`), but a restart run was not
  exercised.

## 9. Open questions

1. **Boundary-file dependency.** The addition is implemented inside
   `fill_from_bdyfiles`, so tides currently require `remora.nc_bdry_file`. ROMS would
   happily run with tides and otherwise-zero boundary arrays. For Moana this is moot
   (the hindcast has boundary files), but if a tide-only configuration is ever wanted,
   the additions would have to be hoisted into a path that runs without boundary data.
2. **`Etide` on the whole grid.** ROMS computes `Etide/Utide/Vtide` over the full tile
   and only uses the boundary rows. This port does the same (cheap, and it keeps the
   door open for `LsshCLM`-style interior use), but if profiling ever shows it matters,
   the kernels could be restricted to boxes touching the physical boundary.
3. **Reference to `angler`.** `vec_angler` is now read from the grid file only when
   `remora.tides=true`, because existing REMORA test grid files are not required to
   contain an `angle` variable. If other ROMS features that need `angler` are ported
   (e.g. rotating vector forcing onto the grid), that gate should be widened.
