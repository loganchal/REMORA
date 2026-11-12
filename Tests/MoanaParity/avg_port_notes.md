# Time-averaged NetCDF output

Moana uses ROMS `AVERAGES` with `NAVG = 864` and a 100 s baroclinic step, giving one mean per day. REMORA accumulates one sample after each completed baroclinic step and writes a record for each complete window.

| Input | Default | Meaning |
|---|---:|---|
| `remora.avg_int` | −1 | Baroclinic steps per window. A value ≤ 0 disables averaging and accumulator allocation. |
| `remora.avg_file` | `avg` | Output filename prefix. |
| `remora.chunk_avg_file` | true | Split average records across numbered files. |
| `remora.steps_per_avg_file` | −1 | Average **records** per file. With splitting enabled, the default estimates the count from the grid size. |

Averaging requires a Parallel-NetCDF build. The NetCDF writer supports a single grid level.

## Accumulation and output

| ROMS operation in the Moana build of `set_avg.F` | REMORA implementation |
|---|---|
| Select `kstp` for 2-D fields and `nrhs` for 3-D fields (lines 476–502) | Sample the arrays used by the history writer after `timeStep()`. The initial condition is outside the first window. |
| Initialise and sum the fields (lines 674–797 and 1201 onward) | Zero the accumulators and add one sample per step with `MultiFab::Saxpy`. Each accumulator has the layout and ghost cells of its source field, including the boundary row written to NetCDF. No inter-rank sum is needed. |
| Close the window and scale by `1/NAVG` (lines 1782–1800) | When the count reaches `avg_int`, multiply each field by the reciprocal count, write and reset. Normalisation uses the allocated ghost widths, including zero vertical ghost cells on 2-D arrays. |

[Accumulation and masking](../../Source/IO/REMORA_Average.cpp) run on the distributed arrays. The [NetCDF writer](../../Source/IO/REMORA_NCPlotFile.cpp) stamps each record at the **centre of the window**, `t_end − 0.5 × avg_int × dt`. In the July 2009 ROMS average file, the month starts at −4416 h and the first daily record is at −4404 h. Restarting begins a new window because checkpoints omit partial sums. A restart partway through a window therefore changes the subsequent averaging intervals.

| Output | Sample |
|---|---|
| `zeta` | `vec_Zt_avg1`, the free surface averaged over the barotropic substeps and used in REMORA history output. ROMS samples instantaneous `zeta(:,:,kstp)`. |
| `ubar`, `vbar` | Component 0 of the depth-averaged velocity fields. |
| `u`, `v` | The updated face velocities. |
| `temp`, `salt`, `tracer` | The corresponding components of `cons_new` selected for output. |
| `sustr`, `svstr` | Surface momentum stresses, written in N/m². |

Static grid coordinates, metrics and vertical stretching parameters use the history writer's schema. The mean fields exclude the atmospheric forcing group, quadratic products and derived fields such as vorticity. The window length is a step count.

## Component tests

The tests used CPU Upwelling and IdealMiniGrid cases built with AppleClang, MPI and Parallel-NetCDF 1.15.0.

| Test | Result |
|---|---|
| Upwelling, 15 steps, history every step, five steps per mean | Across three windows and nine fields (the fields above except the optional tracer), maximum relative difference from the arithmetic mean of history was 2.03 × 10⁻¹⁶. Maximum absolute difference was 7.1 × 10⁻¹⁵ for salinity near 35. |
| Four steps per mean, 12 steps | Exact agreement for every tested field and window. |
| Five steps per mean, four MPI ranks | All averaged fields were bit-identical to the serial result. |
| Masked IdealMiniGrid | The fill masks for sea level and velocities matched history. Means at unmasked points agreed to 2.09 × 10⁻¹⁶ relative error. In this test the average file masked land temperatures and salinities while history retained their values. |
| NetCDF schema | Variable names, data types, dimension lists and attributes matched history for the tested field set. |
| Averaging disabled | Upwelling plotfile data matched at steps 0 and 20. Variable data also matched in the 12-step NetCDF history test. A difference in the original test's global `start_time` attribute was traced to an uninitialised boundary start time, now initialised to zero. |

The component tests covered single-level CPU runs. They did not exercise partial-window restarts or GPU averaging.
