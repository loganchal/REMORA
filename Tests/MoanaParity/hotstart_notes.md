# ROMS his/rst hot start (G9) — port notes

Moana runs `# undef PERFECT_RESTART`, so ROMS initializes from
zeta/ubar/vbar/u/v/temp/salt of one record of an ININAME file, selected by
`NRREC` (with `NRREC = -1` meaning "last record"). The extension-run protocol
sets `ININAME = nz5km_his_YYYYMM.nc`, `NRREC = -1`. Upstream REMORA always
read record 0 and never read ubar/vbar (it started them at zero), so it could
not reproduce a Moana month-chain start.

## Changes

- `remora.nc_init_record` (default 0). Negative counts back from the end, so
  `-1` is ROMS `NRREC = -1`. Resolved against the file's record count by
  `get_num_time_records()` (new helper in `Source/IO/REMORA_NCFile.H`), then
  passed to the existing `one_time`/`fill_time` slab-read path that
  `NCTimeSeries` already used — no new I/O machinery.
- `remora.init_ubar_from_file` (default false). When true, reads `ubar`/`vbar`
  from the same record into every time slot of `vec_ubar`/`vec_vbar`
  (`read_ubar_from_netcdf`), matching the ROMS non-PERFECT_RESTART set.
- `CMake/FindPNetCDF.cmake`: fall back to `find_path`/`find_library` (honoring
  `PNETCDF_DIR`) when pkg-config is unavailable. Portability only; upstream-worthy.

Defaults reproduce upstream behavior exactly (record 0, ubar/vbar zeroed).

## Verification

Built with `-DREMORA_ENABLE_PNETCDF=ON -DREMORA_ENABLE_MPI=ON` (the NetCDF path
is `#ifdef REMORA_USE_NETCDF`, so a non-NetCDF build does not exercise it).

Test: `IdealMiniGrid` with a synthetic 3-record init file derived from
`idmini_ini_v1_classic64.nc` (remora-data), offsetting temp by +10 K and ubar by
+0.5 m/s per record so the selected record is identifiable.

| `nc_init_record` | record used (log) | initial temp in plotfile |
|---|---|---|
| 0 | 0 of 3 | 20.0 |
| 1 | 1 of 3 | 30.0 |
| -1 | 2 of 3 | 40.0 |

Confirms both the record arithmetic and that the data actually loaded is the
selected record.

## Open

- Model start time is NOT taken from the record's `ocean_time`; the run still
  starts at `t=0` unless set by inputs. ROMS sets `time` from the ini record.
  Needed for the month-chain protocol — follow-up.
- `ubar`/`vbar` are written into every time slot; ROMS's `ini_fields` sets the
  equivalent slots. Fine for a cold-ish start; revisit if exact restart matters.
