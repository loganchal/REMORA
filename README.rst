Moana REMORA
============

Moana REMORA runs the 5 km New Zealand Moana Ocean Hindcast configuration on
NVIDIA GPUs. We ported it so we could build higher-resolution models of New
Zealand waters using the same physics. The existing 5 km hindcast gave us a
direct comparison with ROMS before changing the grid.

What changed
------------

* Added Moana’s vertical grid stretching (Vstretching=5).
* Ported Mellor–Yamada 2.5 vertical mixing.
* Ported isoneutral tracer diffusion.
* Added tidal elevation forcing at open boundaries.
* Added tidal current forcing at open boundaries.
* Added atmospheric pressure forcing.
* Added shortwave radiation penetration below the surface.
* Enabled hot starts from a selected ROMS restart record.
* Added time-averaged NetCDF output.
* Corrected river transport weighting.
* Matched surface forcing times to ROMS.
* Matched open-boundary update times to ROMS.

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
