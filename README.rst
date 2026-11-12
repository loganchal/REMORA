Moana REMORA
============

This fork of REMORA runs the 5 km New Zealand Moana Ocean Hindcast on NVIDIA
GPUs. We ported the existing ROMS configuration so we could use the same
physics on finer grids.

What changed
------------

The fork adds these Moana options.

* Moana’s vertical grid stretching (Vstretching=5).
* Mellor–Yamada 2.5 vertical mixing.
* Isoneutral tracer diffusion.
* Tidal elevation forcing at open boundaries.
* Tidal current forcing at open boundaries.
* Inverse-barometer terms in the pressure gradient and at open boundaries.
* Shortwave radiation penetration through the water column.
* Hot starts from a selected ROMS restart record.
* Time-averaged NetCDF output.

Three existing routines now follow ROMS timing and operation order.

* River transport is interpolated in time before it is distributed over depth.
* Surface forcing is sampled at the start of the time step.
* Open-boundary data are sampled at the start of the time step.

How we tested it
----------------

For July 2009, we compared the GPU run and eight ROMS runs against the same
unperturbed ROMS 3.9 control. Each of the eight runs started with round-off-scale
changes to the initial fields (about one ULP). At month end, 17 of 18 interior diagnostics showed a GPU–control difference
no larger than the biggest ROMS–control difference.
Mixed-layer depth was the exception.

In a two-month drift test on an earlier build, the surface temperature RMS difference grew
from 0.016 to 0.030 °C while the domain-mean bias fell from +0.00025 to
+0.00008 °C. On these tests we consider the 5 km port functionally equivalent
to ROMS and numerically stable over the periods tested.

`Build instructions and test results <Exec/Moana/model-notes.rst>`_ ·
`5 km case file <Exec/Moana/inputs>`_
