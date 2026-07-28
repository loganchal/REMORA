#!/usr/bin/env python3
"""Generate a minimal synthetic ROMS-format test case for the REMORA tidal port.

Writes four classic-format NetCDF files into the output directory:

  tidebox_grd.nc   grid (h, pm/pn, coordinates, masks, Coriolis, angle)
  tidebox_ini.nc   initial state (constant T/S, zero velocity, flat free surface)
  tidebox_bry.nc   open-boundary data, identically zero at two times
  tidebox_tide.nc  one tidal constituent (ROMS TIDENAME file)

Because the boundary file is all zeros, everything the model sees at the open
boundaries is the tidal signal added by the ADD_FSOBC/ADD_M2OBC port, which makes
the response easy to check: with tide_Eamp = A, tide_Ephase = 90 deg and period T,
the boundary free surface is A*sin(2*pi*t/T).

Usage:
    python3 make_tide_test_files.py [outdir] [--eamp A] [--ephase P]
                                    [--cmax U] [--cangle ANG] [--cphase PH]
                                    [--period HOURS] [--angle DEG]
"""

import argparse
import os

import numpy as np
from netCDF4 import Dataset

FMT = "NETCDF3_64BIT_OFFSET"  # CDF-2: readable by PnetCDF, which REMORA uses

# Domain: 10 x 16 cells of 10 km, 100 m deep
NX, NY, NZ = 10, 16, 4
DX, DY = 10000.0, 10000.0
DEPTH = 100.0


def _def(nc, name, dims, data, **attrs):
    var = nc.createVariable(name, "f8", dims)
    var[:] = data
    for k, v in attrs.items():
        setattr(var, k, v)
    return var


def write_grid(path, angle_deg):
    nc = Dataset(path, "w", format=FMT)
    nc.createDimension("xi_rho", NX + 2)
    nc.createDimension("eta_rho", NY + 2)
    nc.createDimension("xi_u", NX + 1)
    nc.createDimension("eta_u", NY + 2)
    nc.createDimension("xi_v", NX + 2)
    nc.createDimension("eta_v", NY + 1)
    nc.createDimension("xi_psi", NX + 1)
    nc.createDimension("eta_psi", NY + 1)

    # rho points sit at cell centers; index 0 is the ghost row/column, as in ROMS
    xr = (np.arange(NX + 2) - 0.5) * DX
    yr = (np.arange(NY + 2) - 0.5) * DY
    xu = np.arange(NX + 1) * DX
    yv = np.arange(NY + 1) * DY

    Xr, Yr = np.meshgrid(xr, yr)
    Xu, Yu = np.meshgrid(xu, yr)
    Xv, Yv = np.meshgrid(xr, yv)
    Xp, Yp = np.meshgrid(xu, yv)

    _def(nc, "h", ("eta_rho", "xi_rho"), np.full_like(Xr, DEPTH), units="meter")
    _def(nc, "pm", ("eta_rho", "xi_rho"), np.full_like(Xr, 1.0 / DX), units="meter-1")
    _def(nc, "pn", ("eta_rho", "xi_rho"), np.full_like(Xr, 1.0 / DY), units="meter-1")
    _def(nc, "x_rho", ("eta_rho", "xi_rho"), Xr, units="meter")
    _def(nc, "y_rho", ("eta_rho", "xi_rho"), Yr, units="meter")
    _def(nc, "x_u", ("eta_u", "xi_u"), Xu, units="meter")
    _def(nc, "y_u", ("eta_u", "xi_u"), Yu, units="meter")
    _def(nc, "x_v", ("eta_v", "xi_v"), Xv, units="meter")
    _def(nc, "y_v", ("eta_v", "xi_v"), Yv, units="meter")
    _def(nc, "x_psi", ("eta_psi", "xi_psi"), Xp, units="meter")
    _def(nc, "y_psi", ("eta_psi", "xi_psi"), Yp, units="meter")
    _def(nc, "f", ("eta_rho", "xi_rho"), np.zeros_like(Xr), units="second-1")
    _def(nc, "mask_rho", ("eta_rho", "xi_rho"), np.ones_like(Xr))
    _def(nc, "mask_u", ("eta_u", "xi_u"), np.ones_like(Xu))
    _def(nc, "mask_v", ("eta_v", "xi_v"), np.ones_like(Xv))
    # ROMS grid angle is stored in RADIANS (get_grid.F applies no scale factor)
    _def(nc, "angle", ("eta_rho", "xi_rho"), np.full_like(Xr, np.deg2rad(angle_deg)),
         units="radians", long_name="angle between XI-axis and EAST")
    nc.close()


def write_init(path):
    nc = Dataset(path, "w", format=FMT)
    nc.createDimension("xi_rho", NX + 2)
    nc.createDimension("eta_rho", NY + 2)
    nc.createDimension("xi_u", NX + 1)
    nc.createDimension("eta_u", NY + 2)
    nc.createDimension("xi_v", NX + 2)
    nc.createDimension("eta_v", NY + 1)
    nc.createDimension("s_rho", NZ)
    nc.createDimension("ocean_time", 1)

    _def(nc, "ocean_time", ("ocean_time",), np.zeros(1), units="days")
    _def(nc, "temp", ("ocean_time", "s_rho", "eta_rho", "xi_rho"),
         np.full((1, NZ, NY + 2, NX + 2), 10.0), units="Celsius")
    _def(nc, "salt", ("ocean_time", "s_rho", "eta_rho", "xi_rho"),
         np.full((1, NZ, NY + 2, NX + 2), 32.0))
    _def(nc, "u", ("ocean_time", "s_rho", "eta_u", "xi_u"),
         np.zeros((1, NZ, NY + 2, NX + 1)), units="meter second-1")
    _def(nc, "v", ("ocean_time", "s_rho", "eta_v", "xi_v"),
         np.zeros((1, NZ, NY + 1, NX + 2)), units="meter second-1")
    _def(nc, "zeta", ("ocean_time", "eta_rho", "xi_rho"),
         np.zeros((1, NY + 2, NX + 2)), units="meter")
    nc.close()


def write_bry(path, ndays=10.0):
    """All-zero boundary data at two times, so the whole boundary signal is tidal."""
    nc = Dataset(path, "w", format=FMT)
    nc.createDimension("xi_rho", NX + 2)
    nc.createDimension("eta_rho", NY + 2)
    nc.createDimension("xi_u", NX + 1)
    nc.createDimension("eta_v", NY + 1)
    nc.createDimension("s_rho", NZ)
    nc.createDimension("ocean_time", 2)
    nt = 2

    _def(nc, "ocean_time", ("ocean_time",), np.array([0.0, ndays]), units="days")

    for side, dim2d, n2d in (("west", "eta_rho", NY + 2), ("east", "eta_rho", NY + 2),
                             ("south", "xi_rho", NX + 2), ("north", "xi_rho", NX + 2)):
        _def(nc, "zeta_" + side, ("ocean_time", dim2d), np.zeros((nt, n2d)), units="meter")
        _def(nc, "temp_" + side, ("ocean_time", "s_rho", dim2d),
             np.full((nt, NZ, n2d), 10.0), units="Celsius")
        _def(nc, "salt_" + side, ("ocean_time", "s_rho", dim2d), np.full((nt, NZ, n2d), 32.0))

    # ubar/u live on eta_rho for west/east and on xi_u for south/north
    for side, dim2d, n2d in (("west", "eta_rho", NY + 2), ("east", "eta_rho", NY + 2),
                             ("south", "xi_u", NX + 1), ("north", "xi_u", NX + 1)):
        _def(nc, "ubar_" + side, ("ocean_time", dim2d), np.zeros((nt, n2d)), units="meter second-1")
        _def(nc, "u_" + side, ("ocean_time", "s_rho", dim2d), np.zeros((nt, NZ, n2d)),
             units="meter second-1")

    # vbar/v live on eta_v for west/east and on xi_rho for south/north
    for side, dim2d, n2d in (("west", "eta_v", NY + 1), ("east", "eta_v", NY + 1),
                             ("south", "xi_rho", NX + 2), ("north", "xi_rho", NX + 2)):
        _def(nc, "vbar_" + side, ("ocean_time", dim2d), np.zeros((nt, n2d)), units="meter second-1")
        _def(nc, "v_" + side, ("ocean_time", "s_rho", dim2d), np.zeros((nt, NZ, n2d)),
             units="meter second-1")
    nc.close()


def varying_fields(nconst):
    """Spatially varying harmonic constants.

    Constant fields cannot distinguish the ROMS staggering (the 0.5 averages of
    Etide across an edge and of Uwrk/Vwrk onto u/v points), so this profile makes
    every field vary with both i and j.
    """
    jj, ii = np.meshgrid(np.arange(NY + 2), np.arange(NX + 2), indexing="ij")
    out = {}
    for name, base in (("Eamp", 1.0), ("Ephase", 90.0), ("Cmax", 0.5),
                       ("Cmin", 0.2), ("Cangle", 10.0), ("Cphase", 90.0)):
        flds = []
        for c in range(nconst):
            if name == "Eamp":
                f = 1.0 + 0.10 * ii + 0.05 * jj + 0.5 * c
            elif name == "Ephase":
                f = 90.0 + 3.0 * ii + 2.0 * jj + 17.0 * c
            elif name == "Cmax":
                f = 0.5 + 0.02 * ii + 0.01 * jj + 0.1 * c
            elif name == "Cmin":
                f = 0.2 - 0.01 * ii + 0.005 * jj
            elif name == "Cangle":
                f = 10.0 + 2.0 * ii + 1.0 * jj
            else:
                f = 90.0 + 3.0 * ii - 2.0 * jj + 11.0 * c
            flds.append(f)
        out[name] = np.stack(flds)
    return out


def write_tide(path, period_hours, eamp, ephase, cmax, cmin, cangle, cphase, vary=False):
    """ROMS TIDENAME file. Units follow varinfo.dat: period in HOURS, phases and
    angles in DEGREES, amplitudes in m and m/s."""
    nc = Dataset(path, "w", format=FMT)
    nc.createDimension("xi_rho", NX + 2)
    nc.createDimension("eta_rho", NY + 2)
    nc.createDimension("tide_period", len(period_hours))
    shape = (len(period_hours), NY + 2, NX + 2)

    def const(v):
        return np.full(shape, v, dtype="f8") if np.isscalar(v) else \
            np.stack([np.full(shape[1:], x) for x in v])

    if vary:
        fields = varying_fields(len(period_hours))
    else:
        fields = {"Eamp": const(eamp), "Ephase": const(ephase), "Cmax": const(cmax),
                  "Cmin": const(cmin), "Cangle": const(cangle), "Cphase": const(cphase)}

    _def(nc, "tide_period", ("tide_period",), np.asarray(period_hours), units="hours")
    _def(nc, "tide_Eamp", ("tide_period", "eta_rho", "xi_rho"), fields["Eamp"], units="meter")
    _def(nc, "tide_Ephase", ("tide_period", "eta_rho", "xi_rho"), fields["Ephase"], units="degrees")
    _def(nc, "tide_Cmax", ("tide_period", "eta_rho", "xi_rho"), fields["Cmax"], units="meter second-1")
    _def(nc, "tide_Cmin", ("tide_period", "eta_rho", "xi_rho"), fields["Cmin"], units="meter second-1")
    _def(nc, "tide_Cangle", ("tide_period", "eta_rho", "xi_rho"), fields["Cangle"], units="degrees")
    _def(nc, "tide_Cphase", ("tide_period", "eta_rho", "xi_rho"), fields["Cphase"], units="degrees")
    nc.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir", nargs="?", default=".")
    ap.add_argument("--period", type=float, nargs="+", default=[12.0], help="hours")
    ap.add_argument("--eamp", type=float, nargs="+", default=[1.0], help="meter")
    ap.add_argument("--ephase", type=float, nargs="+", default=[90.0], help="degrees")
    ap.add_argument("--cmax", type=float, nargs="+", default=[0.0], help="m/s")
    ap.add_argument("--cmin", type=float, nargs="+", default=[0.0], help="m/s")
    ap.add_argument("--cangle", type=float, nargs="+", default=[0.0], help="degrees")
    ap.add_argument("--cphase", type=float, nargs="+", default=[90.0], help="degrees")
    ap.add_argument("--angle", type=float, default=0.0, help="grid angle, degrees")
    ap.add_argument("--prefix", default="tidebox")
    ap.add_argument("--vary", action="store_true",
                    help="spatially varying harmonic constants (tests the ROMS staggering)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    p = os.path.join(args.outdir, args.prefix)
    write_grid(p + "_grd.nc", args.angle)
    write_init(p + "_ini.nc")
    write_bry(p + "_bry.nc")
    write_tide(p + "_tide.nc", args.period, args.eamp, args.ephase,
               args.cmax, args.cmin, args.cangle, args.cphase, args.vary)
    print("wrote", p + "_{grd,ini,bry,tide}.nc")


if __name__ == "__main__":
    main()
