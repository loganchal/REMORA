# Barotropic tidal forcing at open boundaries

The Moana ROMS configuration adds harmonic tidal elevations and depth-averaged currents to open-boundary data (`SSH_TIDES`, `UV_TIDES`, `ADD_FSOBC`, `ADD_M2OBC`). The calculation follows `ROMS/Nonlinear/set_tides.F`. Line numbers below refer to the Moana build of that routine after preprocessing.

## Inputs and units

| input | default | meaning | ROMS equivalent |
|---|---|---|---|
| `remora.tides` | `false` | enable barotropic tidal boundary forcing | `SSH_TIDES`+`UV_TIDES`+`ADD_FSOBC`+`ADD_M2OBC` |
| `remora.tide_start` | `0.0` | reference time origin for tidal phase, **days** | `TIDE_START` (Moana: `0.0d0`) |
| `remora.nc_tide_file` | `""` | tidal harmonic constants file | `TIDENAME` |
| `remora.time_ref` | required with tides | model reference epoch for input times | `TIME_REF` |

Tides require one grid level, a Parallel-NetCDF build, `remora.ic_bc_type=netcdf`, a boundary file, a grid file with `angle`, and a tidal constituent file. The [NetCDF reader](../../Source/IO/REMORA_NCTideData.cpp) loads the grid angle and harmonic constants at initialisation. Conversions follow the ROMS `varinfo.dat` scale factors.

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

## Harmonic calculation

[`set_tides`](../../Source/TimeIntegration/REMORA_set_tides.cpp) evaluates the constituents at the start of each baroclinic step by default and holds the values across its barotropic substeps. The diagnostic `remora.tide_time_shift` offsets the evaluation time in seconds. The phase origin is `tide_start` in days relative to the model epoch. Constituents with non-positive periods are skipped. Moana uses the full tidal amplitude from the first step (`ramp = 1`).

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

The elevation sum follows ROMS lines 657–673. For velocity (lines 779–818), each u or v face evaluates the ellipse at its two adjacent rho points, then averages and masks the result. The land mask is applied within the constituent loop. The fields cover the grid; the boundary condition uses the edge values.

## Boundary values

ROMS includes the outer rho-point ring in its grid indices. Subtracting one from the ROMS rho or face index gives the corresponding REMORA index. Here `dlo` and `dhi` are the bounds of the cell-centred domain.

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

[`fill_from_bdyfiles`](../../Source/BoundaryConditions/REMORA_BoundaryConditions_netcdf.cpp) adds the tide to the local interpolated boundary value. Repeated boundary fills in a step therefore leave the stored input intact. Under Moana’s Chapman free-surface and Flather velocity conditions, the Flather calculation uses the tidal elevation through `bry_val_zeta`. Before the first tidal calculation, the three tidal fields are zero.

Moana omits the ROMS tidal ramp, detided averages and tidal additions to sea-level or barotropic velocity climatology. The implementation requires open-boundary data even for a synthetic case with zero background boundary values.

## Component tests

The tests used a CPU Release build with MPI and Parallel-NetCDF. With tidal forcing disabled, both cases reproduced the original model output.

| case | result |
|---|---|
| `Exec/Upwelling/inputs`, 10 steps, plotfiles at 5 and 10 | Plotfile data matched the baseline bit for bit, excluding run metadata. |
| `Tests/MoanaParity/inputs.tidebox` with `remora.tides=false` (NetCDF grid + boundary path, clamped BCs), 40 steps | Plotfiles matched the baseline bit for bit. |

With tides enabled, the [generator](make_tide_test_files.py) produced a 10 × 16 × 4 grid, 100 m deep, with zero background boundary values. The [boundary comparison](check_tide_boundaries.py) evaluates the ROMS tidal equations independently in NumPy. The [test case](inputs.tidebox) predates the model epoch requirement; supply `remora.time_ref` when running it with the release.

| configuration | result |
|---|---|
| 1 constituent, T=12 h, `Eamp`=1 m, `Ephase`=90° (so `zeta_bry = sin(2πt/T)`), clamped BCs, 288 steps = 1 period | boundary `zeta` on all four edges matches the analytic tide to **4.4e-16** |
| 1 constituent, `Cmax`=0.5 m/s, `Cangle`=0, grid `angle`=0 | `ubar`/`vbar` on all 8 boundary arrays match to **≤1.7e-16** |
| Single current constituent, grid `angle`=30° (ellipse rotated by `-angler`: `u` amp 0.43301, `v` amp −0.25) | all 8 match to **≤2.3e-16** |
| 2 constituents (12 h + 6 h), `tide_start`=0.25 d, 4 MPI ranks, `max_grid_size=8` | boundary `zeta` matches the analytic superposition to **4.4e-16** |
| 2 constituents (12 h + 6.2 h), **spatially varying** `Eamp/Ephase/Cmax/Cmin/Cangle/Cphase`, grid `angle`=25°, `tide_start`=0.3 d | all 12 boundary arrays agree with the independent NumPy calculation to **≤6.7e-16**, serial and with 4 MPI ranks / `max_grid_size=6` |

The spatially varying case tests the interpolation from rho points to velocity faces. With Chapman and Flather boundaries, a 1 m boundary tide produced a 1.02 m interior amplitude at the forcing period. With tides disabled, the basin remained at rest (maximum absolute elevation zero).

These tests used an all-water grid on CPU and compared against the harmonic calculation. They did not exercise shoreline masking or a checkpoint restart. For the GPU comparisons with ROMS on the Moana grid, see the [model notes](../../Exec/Moana/model-notes.rst).
