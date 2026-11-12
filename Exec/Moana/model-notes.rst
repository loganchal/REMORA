Moana implementation and tests
=============================

What Moana needed
-----------------

Upstream REMORA already supplied the GPU ocean solver, the curvilinear grid,
split-explicit time stepping and the Chapman and Flather open boundaries.
At the fork point (``e8d2b6d``), it could read a ROMS grid but could not
put Moana's 50 layers at the right depths, apply its tidal constituents at the
boundaries, or use its vertical and lateral mixing schemes. The reference for
the port was the Fortran selected by the Moana build options and the settings
in its run file. `Original Moana hindcast configuration <https://github.com/joaometocean/moana_hindcast>`_.

* **Vertical grid.** Moana uses ROMS Vtransform=2 with Vstretching=5.
  Base REMORA implemented stretching option 4. The new option reproduced
  the reference stretching curves bit for bit. Reconstructed depths agreed
  to within 8 × 10⁻¹⁶ in relative terms.
* **Tides.** The fork adds harmonic sea-level and depth-averaged current
  forcing to the existing open boundaries. We tested a rotated grid with
  two spatially varying constituents. All twelve boundary arrays agreed
  with an independent calculation to within 6.7 × 10⁻¹⁶ in their respective
  units, in serial and on four MPI ranks. `Tidal boundary tests <../../Tests/MoanaParity/tides_port_notes.md>`_.
* **Vertical mixing.** We ported the ROMS Mellor–Yamada 2.5 closure with
  Kantha–Clayson stability functions and horizontal averaging of buoyancy
  and shear terms. The upstream GLS option used different equations.
  A 500-step upwelling test checked bounded mixing coefficients and turbulent
  kinetic energy. We also compared the mixing fields with ROMS near the coast,
  where the differences were largest. `MY2.5 implementation and tests <../../Tests/MoanaParity/my25_port_notes.md>`_.
* **Lateral tracer mixing.** Moana diffuses temperature and salinity along
  density surfaces. Base REMORA offered mixing along terrain-following
  surfaces or at constant depth. The new rotated isoneutral operator gave a
  distinct tracer tendency at the Moana diffusivity of 25 m²/s and ran for
  100 steps with finite output. This component test used a linear equation
  of state. `Diffusion implementation and tests <../../Tests/MoanaParity/t3dmix_iso_port_notes.md>`_.
* **Atmospheric forcing.** We added penetrating shortwave radiation and
  atmospheric pressure in both the interior pressure gradient and the
  Flather boundary condition. We corrected relative-humidity scaling and
  the time at which forcing is sampled. REMORA records the forcing applied
  during a step, while ROMS records the forcing for the next step. After
  aligning the records by one 100 s step, shortwave and rain agreed exactly
  in the forcing comparison.
* **ROMS inputs and output.** Base REMORA read the first initial record
  and started depth-averaged velocities at zero. The fork can select a ROMS
  history or restart record and read its barotropic velocity. A three-record
  file with distinct temperatures verified first, middle and last record
  selection. We added per-variable forcing time axes, CF epoch conversion
  and ROMS-style time averages alongside instantaneous history. Across nine
  fields the five-step averages agreed with the means calculated from history
  to 2.03 × 10⁻¹⁶ relative error. `Hot-start test <../../Tests/MoanaParity/hotstart_notes.md>`_ ·
  `Averaging tests <../../Tests/MoanaParity/avg_port_notes.md>`_.

Rivers, nudging and radiation boundaries were already present in REMORA.
Matching the hindcast also required fixing river transport weighting and the
substep at which river momentum is applied, aligning boundary and climatology
interpolation with ROMS, and correcting halo fills and GPU column dependencies.

Comparison with ROMS
--------------------

We consider the fork functionally equivalent to ROMS and numerically stable
for the 5 km Moana configuration over the periods tested. The comparisons used the Moana solver at revision
`bd85f0d <https://github.com/loganchal/REMORA/commit/bd85f0d959da244bb07335417a3dca21533df733>`_
except where an earlier build is identified.

**Short-step differences.** We perturbed the ROMS initial fields at wet points
by roughly one unit in the last place (ULP) and ran ROMS again. At step 40,
the GPU–ROMS differences in sea level and eastward depth-averaged velocity
were below the ROMS perturbation differences in all six distance-from-land
bands. In the band one to two cells from land, the sea-level differences were
2.22 × 10⁻¹² m for the port and 1.42 × 10⁻¹¹ m for the perturbed ROMS run.
Mixing was more sensitive at the coast: relative differences in viscosity
and turbulent kinetic energy were of order 10⁻¹¹ there.

**July 2009, 26,784 steps.** We compared the GPU run with a ROMS control
and eight ROMS runs with independent round-off-scale perturbations to the
same initial state. We scored 18 month-end diagnostics, including mixed-layer
depth and stratification, over the interior outside a 40-cell rim. Only
mixed-layer depth exceeded the largest corresponding ROMS ensemble difference.
The prespecified failure threshold was five exceedances. Two GPU runs with
identical inputs also reproduced sea level and both depth-averaged velocities
bit for bit at the first and last records checked. July trajectories had largely
decorrelated by month end, which limits how small a defect this test can detect.

**Two-month drift.** In an earlier build, REMORA carried its own July
checkpoint into August for a total of 53,568 steps. The surface temperature
RMS difference from ROMS grew from 0.016 to 0.030 °C. Over the same interval
the domain-mean bias fell from +0.00025 to +0.00008 °C. The August difference
includes a restart difference. ROMS reinitialised its turbulence prognostics
while REMORA carried them across the month boundary.

**Matched restarts.** Separate August and September 2009 runs began from
ROMS restarts with matching forcing. Month-end interior temperature RMS
differences at model level 20 were 0.0040 and 0.0028 °C, respectively.
We excluded October 2009–June 2010 from this comparison after finding that
the reference and GPU runs used different generations of the boundary and
climatology files.

Running the 5 km case
---------------------

The `Moana case file <inputs>`_ documents the ROMS options and
names the eight external NetCDF inputs: grid, initial state, open boundary
fields, meteorology, tides, climatology, nudging coefficients and rivers.
The example uses January 2021. Supply files readable by Parallel-NetCDF
(for example CDF-5), with times covering the run. The river file needs the
monthly climatology repeated over the simulation dates. For a ROMS comparison,
match the river tracer switches and the forcing files to the reference run.

This double-precision build targets an NVIDIA H100. It requires CMake, CUDA,
MPI and Parallel-NetCDF. Set ``PNETCDF_DIR`` to the Parallel-NetCDF installation
prefix and change ``CMAKE_CUDA_ARCHITECTURES`` for a different GPU.

.. code-block:: sh

   git clone --branch moana/dev --recurse-submodules https://github.com/loganchal/REMORA.git
   cd REMORA
   cmake -S . -B build-moana -DCMAKE_BUILD_TYPE=Release \
     -DREMORA_DIM=3 -DREMORA_PRECISION=DOUBLE -DREMORA_ENABLE_CUDA=ON \
     -DCMAKE_CUDA_ARCHITECTURES=90 -DREMORA_ENABLE_MPI=ON \
     -DREMORA_ENABLE_PNETCDF=ON -DREMORA_ENABLE_TESTS=OFF \
     -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx
   cmake --build build-moana --parallel

From the directory containing the input files, run a ten-step check:

.. code-block:: sh

   /path/to/REMORA/build-moana/Exec/remora_exec /path/to/REMORA/Exec/Moana/inputs remora.max_step=10

For a finer grid, prepare the bathymetry, initial and boundary fields,
nudging, rivers and tides on the new grid. Revisit the time step, barotropic
subcycling, viscosity and diffusivity for the new spacing. The comparison
above provides a baseline at 5 km. A 2 km hindcast is a separate experiment.

Upstream REMORA
---------------

REMORA is a GPU ocean model built on AMReX by the SEAHORCE SciDAC team.
See the `upstream repository <https://github.com/seahorce-scidac/REMORA>`_
and `REMORA documentation <https://remora-ocean.readthedocs.io/en/latest/>`_
for the solver, build options and other example configurations.

License
~~~~~~~~~

REMORA Copyright (c) 2023, The Regents of the University of California,
through Lawrence Berkeley National Laboratory, Triad National Security, LLC
through Los Alamos National Laboratory, UChicago Argonne LLC through
Argonne National Laboratory, Battelle Memorial Institute through Pacific
Northwest National Laboratory (subject to receipt of any required approvals
from the U.S. Dept. of Energy) and TechTrans International, Inc. All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative
works, and perform publicly and display publicly, and to permit others to do so.

The license for REMORA can be found in the LICENSE.md file.

Citation
~~~~~~~~~

To cite REMORA, please use the following reference:

.. image:: https://joss.theoj.org/papers/10.21105/joss.07958/status.svg
   :target: https://doi.org/10.21105/joss.07958

::

  @article{REMORA_JOSS,
      doi = {10.21105/joss.07958},
      url = {https://doi.org/10.21105/joss.07958},
      year = {2025},
      publisher = {The Open Journal},
      volume = {10},
      number = {110},
      pages = {7958},
      author = {Hannah Klion and Robert Hetland and Jean Sexton and Ann Almgren and Iulian Grindeanu and Kyle Hinson and Vijay Mahadevan},
      title = {REMORA: Regional Modeling of Oceans Refined Adaptively (built on AMReX)},
      journal = {Journal of Open Source Software}
  }
