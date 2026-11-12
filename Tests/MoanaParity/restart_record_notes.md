# Initialisation from ROMS history and restart files

The Moana ROMS configuration uses `PERFECT_RESTART` disabled. It reads sea level, horizontal velocities, temperature and salinity from one time record in `ININAME`. For monthly runs, `ININAME = nz5km_his_YYYYMM.nc` and `NRREC = -1` select the last record in the history file.

## Record and time selection

| REMORA input | Meaning |
|---|---|
| `remora.nc_init_record` | Zero-based record index. The default is 0. Negative indices count from the end of the file, so −1 selects the last record. |
| `remora.init_ubar_from_file` | Read `ubar` and `vbar` from the selected record and populate the barotropic velocity time slots. Defaults to false. |
| `remora.start_time` | Model start time in seconds. Set it to the time of the selected record. Selecting a record does not set the model clock from its `ocean_time`. |

The [NetCDF initial-state reader](../../Source/IO/REMORA_ReadFromInitNetcdf.cpp) obtains the number of records, resolves the index and rejects indices outside the file. The same record selection applies to `zeta`, `temp`, `salt`, `u`, `v` and, when requested, `ubar` and `vbar`. The reader requires Parallel-NetCDF and MPI. CMake can locate Parallel-NetCDF under `PNETCDF_DIR` when pkg-config is unavailable.

## Record selection test

The IdealMiniGrid test used a three-record file derived from `idmini_ini_v1_classic64.nc`. Each successive record increased temperature by 10 K and `ubar` by 0.5 m/s. The initial temperature in the output identified the record actually read.

| Requested index | Record read (zero-based) | Output temperature (°C) |
|---:|---:|---:|
| 0 | 0 | 20.0 |
| 1 | 1 | 30.0 |
| −1 | 2 | 40.0 |

This test covered initial field selection. It did not measure continuity of a subsequent integration across a restart.
