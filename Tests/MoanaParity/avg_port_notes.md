# ROMS AVERAGES (avg file) port notes

Port of ROMS time-averaged output into REMORA, for the Moana-ROMS GPU port.
Moana runs `AVERAGES` with `NAVG = 864` (daily at `DT = 100 s`) and the downstream
atlas analysis pipeline consumes ROMS `avg` files, so this is v1-scope.

Spec: `external/moana_spec/set_avg.f90` — the preprocessed `ROMS/Nonlinear/set_avg.F`
from the Moana build (line numbers below refer to that preprocessed file).

## Inputs

| input | meaning | default |
|---|---|---|
| `remora.avg_int` | baroclinic steps per averaging window; ROMS `NAVG`. `<= 0` disables | `-1` (OFF) |
| `remora.avg_file` | prefix of the avg NetCDF file; output is `<prefix>_d01.nc` | `avg` |

Default is OFF, so an existing input deck is unaffected: no accumulators are
allocated, no extra kernels run, and no file is written.
`avg_int > 0` requires a PnetCDF-enabled build (checked in `ReadParameters`, aborts otherwise).

## Files touched

| file | change |
|---|---|
| `Source/IO/REMORA_Average.cpp` | **new** — all accumulate / normalize / reset / mask logic |
| `Source/REMORA.H` | new public `Average*` / `mask_avg_arrays_for_write` methods; accumulator MultiFabs and counters; `is_avg` parameter on the two NetCDF writer declarations |
| `Source/REMORA.cpp` | `avg_int` / `avg_file` parsing; `AverageAtIntermediateTime()` call in `Evolve()`; non-NetCDF-build abort |
| `Source/IO/REMORA_NCPlotFile.cpp` | `is_avg` mode: file naming, record index, `nt` sizing, source-MultiFab selection, masking |
| `CMake/BuildREMORAExe.cmake`, `Source/IO/Make.package` | register the new source file |
| `Docs/sphinx_doc/Plotfiles.rst` | document the two new inputs |

## Spec line ranges to code

| ROMS `set_avg.f90` | semantics | REMORA code |
|---|---|---|
| 476–482, 500–502, 535–537 | `Kout = kstp(ng)` (2D time index), `Nout = nrhs(ng)` (3D time index) — the time levels sampled | `REMORA_Average.cpp: AverageAccumulate()` — see "time index" below |
| 629 | `IF (nAVG(ng).eq.0) RETURN` — averaging disabled | `AverageAccumulate()`: `if (avg_int <= 0) return;` |
| 674–677 | window start / reset test: `MOD(iic-1,nAVG)==1` (or `nAVG==1`, or restart) | `AverageAtIntermediateTime()`: window closes when `avg_nsamples == avg_int`, then `AverageReset()` starts the next window |
| 686, 695, 703 | *initialize* the accumulator with the first sample of the window (`avgzeta = zeta(Kout)`) | folded into "zero, then add all `nAVG` samples" — see "reset" below |
| 728, 738, 797 | initialize 3D accumulators (`avgu3d`, `avgv3d`, `avgt`) | same |
| **1201–1360 onward** | **accumulate block**: `avgzeta = avgzeta + zeta(Kout)`, `avgu2d += ubar(Kout)`, `avgv2d += vbar(Kout)`, `avgu3d += u(Nout)`, `avgv3d += v(Nout)`, `avgt += t(Nout,itrc)` — a plain running **sum**, one term per baroclinic step | `AverageAccumulate()`, one `amrex::MultiFab::Saxpy(dst, 1.0, src, ...)` per field |
| 1782–1790 | window-end test and timestamp: `AVGtime = AVGtime + REAL(nAVG)*dt` (or `AVGtime = time` when `nAVG == 1`) | `AverageWriteAndReset()`; timestamp written is `t_new[0]` in `WriteNCPlotFile_which` — see "timestamp" below |
| 1793–1800 | `fac = 1.0_r8/REAL(nAVG,r8)`, applied via `pfac/rfac/ufac/vfac` to every C-grid variable type | `AverageNormalize(fac)` with `fac = 1/avg_nsamples` |

## Semantics as implemented

**Accumulation.** Plain running sum, one sample per *baroclinic* step, added on
device via `MultiFab::Saxpy` over valid **and** ghost cells. Nothing is copied to
the host, and no MPI communication is needed — every accumulator shares the
`BoxArray`, `DistributionMapping` and ghost region of the field it accumulates, so
each rank sums exactly its own data. The 4-rank MPI avg file is bit-identical to
the serial one (verified, see below).

Accumulating the ghost layer is the analogue of the ROMS `IstrR:IendR` /
`JstrR:JendR` loop bounds: the NetCDF writer emits one ghost row at the domain
boundary, and that row must be averaged too.

**Which time index each variable is accumulated at.** ROMS samples
`zeta/ubar/vbar` at `Kout = kstp(ng)` and `u/v/t` at `Nout = nrhs(ng)`
(`set_avg.f90:476-482` passes `nrhs(ng)` as `Nout` and `kstp(ng)` as `Kout`), and
`set_avg` is called once per `iic` on the state left by the previous step. REMORA
has no multi-time-level state arrays with the same indexing; the equivalent is
"the state at the end of a completed baroclinic step", so `AverageAccumulate()` is
called from `Evolve()` immediately after `timeStep()` returns and `cur_time` has
been advanced. Concretely it reads:

| avg variable | REMORA source | ROMS equivalent |
|---|---|---|
| `zeta` | `vec_Zt_avg1[lev]`, comp 0 | `zeta(:,:,Kout)` — see deviation D1 |
| `ubar` | `vec_ubar[lev]`, comp 0 | `ubar(:,:,Kout)` |
| `vbar` | `vec_vbar[lev]`, comp 0 | `vbar(:,:,Kout)` |
| `u` | `xvel_new[lev]` | `u(:,:,:,Nout)` |
| `v` | `yvel_new[lev]` | `v(:,:,:,Nout)` |
| `temp` | `cons_new[lev]`, `Temp_comp` | `t(:,:,:,Nout,itemp)` |
| `salt` | `cons_new[lev]`, `Salt_comp` | `t(:,:,:,Nout,isalt)` |
| `tracer` | `cons_new[lev]`, `Tracer_comp` (only if `tracer` is in `remora.plot_vars_3d`) | `t(:,:,:,Nout,itrc)` |
| `sustr` | `vec_sustr[lev]` | `avgsus` (ROMS `idUsms`) |
| `svstr` | `vec_svstr[lev]` | `avgsvs` (ROMS `idVsms`) |

These are *exactly* the MultiFabs the REMORA history writer reads, at the same
point in the step cycle, which is what makes the "avg == mean of per-step history"
test below a true test of the averaging.

**Window bookkeeping.** ROMS uses modular arithmetic on `iic`: with
`ntsAVG = 1`, `nAVG = N` and `ntstart = 1`, the reset test `MOD(iic-1,N)==1` fires
at `iic = 2, N+2, 2N+2, ...` and the normalization test `MOD(iic-1,N)==0` at
`iic = N+1, 2N+1, ...`, so window *k* holds exactly `N` samples — the states after
steps `(k-1)N+1 ... kN`. The initial condition (state after 0 steps) is **not** in
any window. This port reproduces that directly with a counter: sample after every
step, close the window when `avg_nsamples == avg_int`. `ntsAVG` is effectively 1
(averaging starts at the first step of the run); there is no separate input for it.

**Reset.** ROMS *assigns* the first sample of a window (`avgzeta = zeta(Kout)`,
line 686) and adds the remaining `N-1`. This port zeroes the accumulators and adds
all `N` samples. Arithmetically identical (`0 + x1 + x2 + ...` in the same order),
and it removes the first-sample branch. Reset happens after every write, so an avg
window never straddles two output records.

**Normalization.** `AverageNormalize(1/avg_nsamples)` multiplies every accumulator
(valid + ghost) by `fac`, i.e. ROMS `fac = 1/REAL(nAVG)`. `avg_nsamples` always
equals `avg_int` at that point, so this is `1/nAVG`; using the counter rather than
`avg_int` makes a short final window (if one is ever written) self-consistent
rather than silently biased.

**Timestamp.** ROMS sets `AVGtime = AVGtime + nAVG*dt` at window end, starting
from the run's initial time, i.e. the model time of the *last* sample in the
window (not the window centre). Because this port samples after each completed
step, `t_new[0]` at window close is already that value, so the writer's existing
`ocean_time` put is correct with no extra counter. Verified: with `dt = 300 s` and
`avg_int = 4`, history `ocean_time = 0,300,...,3600` and avg `ocean_time =
1200, 2400, 3600`. The `nAVG == 1` special case (`AVGtime = time`) reduces to the
same expression.

**Writer.** No new writer. `WriteNCPlotFile` / `WriteNCPlotFile_which` gained a
`bool is_avg` (default `false`); in avg mode they select the accumulators instead
of the instantaneous MultiFabs, use `avg_file_name` and `avg_count`, always behave
as a multi-record ("history-style") file, and size the `ocean_time` dimension as
`max_step / avg_int`. Everything else — dimension names (`xi_rho`, `eta_u`,
`s_rho`, ...), variable names, units, `_FillValue`, `grid`/`location`/`coordinates`
attributes, the static grid fields, the per-box offsets — is the shared code path.
Verified: the avg file's schema (names, dtypes, dimensions, every attribute) is
identical to the history file's.

## Variables written to the avg file

Time-averaged: `zeta`, `ubar`, `vbar`, `u`, `v`, `temp`, `salt`, `tracer` (if
plotted), `sustr`, `svstr`.
Static (written once, from the header path): `h`, `pm`, `pn`, `f`, `x_*`, `y_*`,
`s_rho`, `s_w`, `Cs_r`, `Cs_w`, `hc`, `theta_s`, `theta_b`, grid topology.

## Deviations from the literal spec

- **D1 — `zeta` source.** ROMS accumulates `zeta(:,:,kstp)`, the instantaneous free
  surface. REMORA's history writer emits `vec_Zt_avg1` (ROMS `Zt_avg1`, the
  fast-time-average of the free surface over the barotropic substeps of the step)
  under the name `zeta`, so that is what is accumulated here. Keeping avg
  consistent with REMORA's own `zeta` was preferred over matching the ROMS
  variable literally; for `NAVG = 864` the difference is a sub-baroclinic-step
  detail. Change one line in `AverageAccumulate()` if literal parity is wanted.
- **D2 — higher-order and derived ROMS averages not ported.** `set_avg.f90` also
  accumulates quadratic terms (`avgZZ`, `avgU2`, `avgUV`, `avgHuonT`, ...),
  `omega`/`w`, density, vorticity and the rotated-to-east/north velocities
  (`avgu3dE`/`avgv3dN`). None are ported; the atlas pipeline does not use them.
- **D3 — surface-flux group.** When `remora.output_forcing` is on, the history file
  carries `Tair`, `Pair`, `qnet`, `ssflux`, `latent`, `sensible`, `lwrad`, `swrad`,
  `evaporation`, `rain`. These are **not** averaged, and are deliberately neither
  defined nor written in the avg file, so nobody can mistake an end-of-window
  snapshot for a mean.
- **D4 — chunking.** `remora.chunk_history_file` is ignored for the avg file; all
  avg records go into one file.
- **D5 — restart.** The accumulators are not checkpointed. ROMS re-initializes the
  averages on restart when `nrrec > 0` (the third clause of `set_avg.f90:674-677`);
  here a restart simply starts a fresh window at the first step after the restart.
  If the restart is not window-aligned, avg record boundaries shift. Unverified
  (see below).
- **D6 — land masking of `temp`/`salt`.** In the avg file, land points of every
  written field including `temp`/`salt` carry `netcdf_fill_value` (ROMS-like). In
  the *history* file `temp`/`salt` are **not** masked, because `WritePlotFile`
  copies `cons_new` into `plotMF` before `mask_arrays_for_write` runs — a
  pre-existing REMORA quirk, not something this port changed. `zeta`, `u`, `v`,
  `ubar`, `vbar` are masked identically in both files (verified: mask arrays match
  exactly).

## CUDA compatibility (inspection only — no GPU available here)

Every function added or modified by this port that contains an
`amrex::ParallelFor` / extended `__device__` lambda:

| function | file | access | notes |
|---|---|---|---|
| `REMORA::AverageNormalize(amrex::Real)` | `Source/IO/REMORA_Average.cpp` | **public** | one `ParallelFor(bx, nc, ...)`; body is `arr(i,j,k,n) *= fac` |
| `REMORA::mask_avg_arrays_for_write(int, Real, Real)` | `Source/IO/REMORA_Average.cpp` | **public** | six `ParallelFor`s; bodies are `if (msk == 0.0_rt) x = fill_value` |

Functions added by this port that contain **no** device lambda:
`AverageInit`, `AverageReset`, `AverageAccumulate` (uses `MultiFab::Saxpy`),
`AverageWriteAndReset`, `AverageAtIntermediateTime`.

Rules applied:

1. **No device lambda in a private/protected member.** Both lambda-carrying
   functions above are declared `public` in `REMORA.H`, next to `Evolve()`.
2. **No namespace-scope `constexpr` passed by const reference inside device code.**
   Neither kernel calls `std::max`/`std::min`/`std::swap` at all. The only
   constants read inside the lambdas are the prvalue literal `0.0_rt` and the
   by-value captured `fac` / `fill_value` / `fill_where` function arguments. The
   `REMORA_Constants.H` symbols `zero` / `one` appear only in host code in this
   port (`WriteNCPlotFile_which`'s masking calls), never inside a lambda.
3. **No serial dependence across a fused `ParallelFor` index.** Both kernels are
   pointwise; the running sum over the averaging window is across *calls*, in the
   host-side `MultiFab::Saxpy`, not across a kernel index.
4. `mask_avg_arrays_for_write` grows only in `x`/`y` (`IntVect(NGROW+1,NGROW+1,0)`
   / `IntVect(NGROW,NGROW,0)`, matching the existing `mask_arrays_for_write`), and
   `AverageNormalize` uses `mfi.growntilebox()` with no argument so the box is
   clamped to each MultiFab's own allocated ghost region. `MultiFab::mult` was
   deliberately **not** used for the normalization: its `int nghost` argument
   becomes `IntVect(nghost)`, which would grow the 2D slab accumulators
   (`zeta`, `ubar`, `vbar`, `sustr`, `svstr`, allocated with `(n,n,0)` ghosts) in
   the vertical, past their allocation.

## Verification performed (AppleClang, PnetCDF 1.15.0, `-DREMORA_ENABLE_MPI=ON -DREMORA_ENABLE_PNETCDF=ON`)

1. **Clean build.** No warnings or errors introduced.
2. **Default-off bit identity.** Upwelling, 20 steps, AMReX plotfile: `plt00000`
   and `plt00020` are byte-identical to a build of `moana/dev` (only `job_info`
   differs, in timestamp / build dir / git hash).
   NetCDF history, 12 steps, `plot_int = 2`, avg off: every variable's data is
   identical to the `moana/dev` build. See "Known pre-existing issue" below for
   the one differing header byte.
3. **Averaging accuracy — the key test.** Upwelling with
   `plotfile_type = netcdf`, `plot_int = 1` (a history record *every* step),
   `avg_int = 5`, `max_step = 15`. For each of `zeta, ubar, vbar, u, v, temp,
   salt, sustr, svstr` and each of the 3 avg records, the avg field was compared
   against the arithmetic mean of the 5 per-step history fields in that window
   (read with `netCDF4` in double precision):

   > **max relative error = 2.03e-16** (max absolute error 7.1e-15 on `salt`,
   > field scale 35.0), against `float64` eps = 2.22e-16 — i.e. 1 ulp.

   With `avg_int = 4` (a power of two, so `1/N` is exact in binary) over 12 steps,
   the error is **exactly 0.0** for every field and window.
   Within-window spread was non-trivial for the fields that matter (`u`
   ~3e-4 m/s, `temp` ~2e-4 degC, `zeta` ~1e-5 m over a window), so the agreement
   is not the trivial "nothing changed" result.
4. **MPI correctness.** The same `avg_int = 5` case under `mpirun -n 4` gives an
   avg file bit-identical to the serial run (max |serial − 4-rank| = 0 for every
   averaged field), and reproduces the same 2.03e-16 agreement with its own
   per-step history.
5. **Land masking.** IdealMiniGrid with `idmini_grd_masked_v1_classic64.nc`,
   `avg_int = 5`: the avg file's `_FillValue` masks for `zeta`, `ubar`, `vbar`,
   `u`, `v` match the history file's exactly; averaging over unmasked points
   agrees to 2.09e-16. (`temp`/`salt` mask difference is deviation D6 above.)
6. **Schema parity.** The avg file and the history file have identical variable
   sets, dtypes, dimension lists and attributes; avg dimensions are the ROMS names
   (`xi_rho`, `eta_rho`, `s_rho`, `xi_u`, `eta_v`, `s_w`, `ocean_time`, ...).

## Not verified

- **CUDA / GPU.** No GPU on this machine. CUDA compliance is by inspection only,
  per the table above.
- **Restart behaviour (D5).** Not exercised. The accumulators are not
  checkpointed, so a restart begins a fresh window.
- **Multi-level (AMR).** The accumulate/normalize/reset loops are written over
  `0..finest_level`, but the writer asserts `finest_level == 0`, exactly as the
  existing NetCDF history writer does. Only single-level runs were tested.
- **Very long runs / file size.** No 2 GB-limit handling or chunking for the avg
  file. At Moana resolution with `NAVG = 864`, one avg record is the same size as
  one history record, so a multi-year run will need either chunking support or
  per-period runs. Worth revisiting before a production hindcast.
- **`plot_int_time`-style time-based averaging windows.** Only step-count windows
  (`avg_int`) are supported, matching ROMS `NAVG`.

## Known pre-existing issue found while verifying

`REMORA::start_bdy_time` (`Source/REMORA.H`) has **no initializer** and is only
assigned when initializing from a NetCDF boundary file. It is nevertheless written
into every NetCDF output as the global attribute `start_time`
(`REMORA_NCPlotFile.cpp`, `ncf.put_attr("start_time", ...)`). For a run that does
not read NetCDF boundary data (e.g. analytic Upwelling) this writes uninitialized
memory; the value is stable for a given binary but changes when the class layout
changes, which is why the avg-off NetCDF header differs by exactly those 4 bytes
between the `moana/dev` build and this one. All *variable data* is identical. This
is not a regression from this port; a one-word fix (`amrex::Real start_bdy_time =
amrex::Real(0.0);`) would remove it.
