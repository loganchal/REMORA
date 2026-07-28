#!/usr/bin/env python3
"""Oracle check for the REMORA barotropic tide port.

Independently re-implements ROMS set_tides.F (SSH_TIDES + UV_TIDES + ADD_FSOBC +
ADD_M2OBC) in numpy, straight from the Fortran, and compares it against the
boundary values that a REMORA run actually produced.

The comparison relies on the open-boundary data file being identically zero and
the boundary conditions being clamped, so that the ghost/edge value written by
REMORA IS the tidal contribution:

    zeta ghost row/column   == 0.5*(Etide straddling the edge)
    ubar/vbar boundary face == Utide/Vtide at the ROMS boundary index

ROMS -> file index map used below (file arrays carry the ROMS ghost ring, so file
index = ROMS index; REMORA index = ROMS index - 1):

    zeta_west (j) = 0.5*(E[j, 0]    + E[j, 1])       ROMS Etide(Istr-1), Etide(Istr)
    zeta_east (j) = 0.5*(E[j, Nx]   + E[j, Nx+1])    ROMS Etide(Iend),   Etide(Iend+1)
    zeta_south(i) = 0.5*(E[0, i]    + E[1, i])       ROMS Etide(Jstr-1), Etide(Jstr)
    zeta_north(i) = 0.5*(E[Ny, i]   + E[Ny+1, i])    ROMS Etide(Jend),   Etide(Jend+1)
    ubar_west (j) = Utide[j, 0]  (u face 0)          ROMS Utide(Istr)
    ubar_east (j) = Utide[j, Nx] (u face Nx)         ROMS Utide(Iend+1)
    ubar_south(i) = Utide[0, i]  (u row j=-1)        ROMS Utide(Jstr-1)
    ubar_north(i) = Utide[Ny+1, i]                   ROMS Utide(Jend+1)
    vbar_west (j) = Vtide[j, 0]  (v col i=-1)        ROMS Vtide(Istr-1)
    vbar_east (j) = Vtide[j, Nx+1]                   ROMS Vtide(Iend+1)
    vbar_south(i) = Vtide[0, i]  (v face 0)          ROMS Vtide(Jstr)
    vbar_north(i) = Vtide[Ny, i] (v face Ny)         ROMS Vtide(Jend+1)

Usage:
    python3 check_tide_boundaries.py <rundir> [--dt 150] [--tide-start 0.0]
"""

import argparse
import sys

import numpy as np
from netCDF4 import Dataset

DEG2RAD = 0.017453292519943295  # ROMS varinfo.dat Fscale for degrees


def roms_set_tides(tide, angler, time, tide_start):
    """Return (Etide, Utide, Vtide) on the ROMS rho/u/v grids for one time.

    Etide has shape (eta_rho, xi_rho); Utide (eta_rho, xi_u); Vtide (eta_v, xi_rho).
    """
    period = tide["tide_period"][:] * 3600.0            # hours -> s
    Eamp = tide["tide_Eamp"][:]                          # m
    Epha = tide["tide_Ephase"][:] * DEG2RAD              # deg -> rad
    Cmax = tide["tide_Cmax"][:]                          # m/s
    Cmin = tide["tide_Cmin"][:]                          # m/s
    Cang = tide["tide_Cangle"][:] * DEG2RAD              # deg -> rad
    Cpha = tide["tide_Cphase"][:] * DEG2RAD              # deg -> rad

    cff = 2.0 * np.pi * (time - tide_start * 86400.0)

    E = np.zeros(Eamp.shape[1:])
    Uwrk = np.zeros(Eamp.shape[1:])
    Vwrk = np.zeros(Eamp.shape[1:])
    for c in range(len(period)):
        if period[c] <= 0.0:
            continue
        omega = cff / period[c]
        E += Eamp[c] * np.cos(omega - Epha[c])
        angle = Cang[c] - angler
        phase = omega - Cpha[c]
        Uwrk += Cmax[c] * np.cos(angle) * np.cos(phase) - Cmin[c] * np.sin(angle) * np.sin(phase)
        Vwrk += Cmax[c] * np.sin(angle) * np.cos(phase) + Cmin[c] * np.cos(angle) * np.sin(phase)

    # Average the rho-point work arrays onto u and v points (ROMS 0.5 stencils)
    Utide = 0.5 * (Uwrk[:, :-1] + Uwrk[:, 1:])   # (eta_rho, xi_u)
    Vtide = 0.5 * (Vwrk[:-1, :] + Vwrk[1:, :])   # (eta_v,   xi_rho)
    return E, Utide, Vtide


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--history", default="plt_his_d01.nc")
    ap.add_argument("--tide", default="tidebox_tide.nc")
    ap.add_argument("--grid", default="tidebox_grd.nc")
    ap.add_argument("--dt", type=float, default=150.0)
    ap.add_argument("--tide-start", type=float, default=0.0)
    ap.add_argument("--tol", type=float, default=1.0e-12)
    args = ap.parse_args()

    his = Dataset(f"{args.rundir}/{args.history}")
    tide = Dataset(f"{args.rundir}/{args.tide}").variables
    angler = Dataset(f"{args.rundir}/{args.grid}")["angle"][:]

    t = his["ocean_time"][:]
    zeta = his["zeta"][:]
    ubar = his["ubar"][:]
    vbar = his["vbar"][:]
    ny_r, nx_r = zeta.shape[1], zeta.shape[2]
    Nx, Ny = nx_r - 2, ny_r - 2

    worst = {}
    # skip record 0: that is the initial condition, written before any set_tides call
    for it in range(1, len(t)):
        # REMORA evaluates the tide at the start of the baroclinic step (ROMS time(ng)),
        # and history is written at the end of the step
        E, U, V = roms_set_tides(tide, angler, t[it] - args.dt, args.tide_start)

        checks = {
            "zeta_west":  (zeta[it, 1:-1, 0],    0.5 * (E[1:-1, 0] + E[1:-1, 1])),
            "zeta_east":  (zeta[it, 1:-1, -1],   0.5 * (E[1:-1, Nx] + E[1:-1, Nx + 1])),
            "zeta_south": (zeta[it, 0, 1:-1],    0.5 * (E[0, 1:-1] + E[1, 1:-1])),
            "zeta_north": (zeta[it, -1, 1:-1],   0.5 * (E[Ny, 1:-1] + E[Ny + 1, 1:-1])),
            "ubar_west":  (ubar[it, 1:-1, 0],    U[1:-1, 0]),
            "ubar_east":  (ubar[it, 1:-1, -1],   U[1:-1, Nx]),
            "ubar_south": (ubar[it, 0, 1:-1],    U[0, 1:-1]),
            "ubar_north": (ubar[it, -1, 1:-1],   U[Ny + 1, 1:-1]),
            "vbar_west":  (vbar[it, 1:-1, 0],    V[1:-1, 0]),
            "vbar_east":  (vbar[it, 1:-1, -1],   V[1:-1, Nx + 1]),
            "vbar_south": (vbar[it, 0, 1:-1],    V[0, 1:-1]),
            "vbar_north": (vbar[it, -1, 1:-1],   V[Ny, 1:-1]),
        }
        for name, (got, want) in checks.items():
            err = float(np.abs(np.asarray(got) - np.asarray(want)).max())
            worst[name] = max(worst.get(name, 0.0), err)

    ok = True
    print(f"{len(t)-1} output times compared against the numpy ROMS oracle")
    for name, err in worst.items():
        flag = "OK  " if err <= args.tol else "FAIL"
        if err > args.tol:
            ok = False
        print(f"   {flag} {name:11s} max|REMORA - ROMS oracle| = {err:.3e}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
