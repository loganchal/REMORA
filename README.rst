Moana REMORA
============

Moana REMORA runs the 5 km New Zealand Moana Ocean Hindcast configuration on
NVIDIA GPUs. We ported it so we could build higher-resolution models of New
Zealand waters using the same physics. The existing 5 km hindcast gave us a
direct comparison with ROMS before changing the grid.

What changed
------------

REMORA already had the GPU ocean solver. To run Moana, we added:

* Moana’s vertical grid (Vstretching=5), Mellor–Yamada 2.5 vertical mixing
  and isoneutral tracer diffusion.
* Tidal elevations and currents at open boundaries, atmospheric pressure
  forcing and penetration of sunlight below the surface.
* Support for starting from a selected ROMS restart record and writing
  time-averaged NetCDF output.

We also corrected river transport weighting and the timing of forcing and
boundary updates to match ROMS.

How we tested it
----------------

We compared the GPU model with ROMS 3.9 on the Moana grid. For July 2009,
we also ran ROMS eight times with round-off-scale perturbations to the initial
conditions (about one ULP). At month end, 17 of 18 interior diagnostics fell
at or below the largest difference between the ROMS runs. Mixed-layer depth
was the exception. In a separate two-month drift test on an earlier build,
the domain-mean surface temperature bias stayed below 0.0003 °C at both
month ends.

On these tests, we consider the 5 km port functionally equivalent to ROMS
and numerically stable over the periods tested.

`Build instructions and test results <Exec/Moana/model-notes.rst>`_ ·
`Moana configuration <Exec/Moana/inputs>`_
