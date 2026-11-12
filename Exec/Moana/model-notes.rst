Moana build and test notes
==========================

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

   git clone https://github.com/loganchal/REMORA.git
   cd REMORA
   git checkout --detach bd85f0d959da244bb07335417a3dca21533df733
   git submodule update --init --recursive
   cmake -S . -B build-moana -DCMAKE_BUILD_TYPE=Release \
     -DREMORA_DIM=3 -DREMORA_PRECISION=DOUBLE -DREMORA_ENABLE_CUDA=ON \
     -DCMAKE_CUDA_ARCHITECTURES=90 -DREMORA_ENABLE_MPI=ON \
     -DREMORA_ENABLE_PNETCDF=ON -DREMORA_ENABLE_TESTS=OFF \
     -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx
   cmake --build build-moana --parallel

From the directory containing the input files, run a ten-step check:

.. code-block:: sh

   /path/to/REMORA/build-moana/Exec/remora_exec /path/to/REMORA/Exec/Moana/inputs remora.max_step=10

Comparison with ROMS 3.9
------------------------

The tests used solver revision ``bd85f0d`` unless an earlier build is specified.

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
the domain-mean bias fell from +0.00025 to +0.00008 °C. At the July restart, ROMS reinitialised its turbulence prognostics
while REMORA carried them into August.

**Matched restarts.** Separate August and September 2009 runs began from
ROMS restarts with matching forcing. Month-end interior temperature RMS
differences at model level 20 were 0.0040 and 0.0028 °C, respectively.
We excluded October 2009–June 2010 from this comparison after finding that
the reference and GPU runs used different generations of the boundary and
climatology files.

Component tests
---------------

* **Vertical coordinates.** With Vstretching=5 the reference stretching curves
  agreed bit for bit. Reconstructed depths agreed to within 8 × 10⁻¹⁶ relative
  error. `Vertical grid comparison <../../Tests/MoanaParity/v5_compare.py>`_.
* **Tides.** For two spatially varying constituents on a rotated grid, all
  twelve boundary arrays agreed with the independent tidal calculation to
  within 6.7 × 10⁻¹⁶ in their respective units, in serial and on four MPI
  ranks. `Tidal boundary tests <../../Tests/MoanaParity/tides_port_notes.md>`_.
* **Vertical mixing.** Mixing coefficients and turbulent kinetic energy
  remained bounded in a 500-step upwelling test.
  `MY2.5 tests <../../Tests/MoanaParity/my25_port_notes.md>`_.
* **Isoneutral diffusion.** A 100-step test at 25 m²/s with a linear equation
  of state produced finite tracer fields. The tracer tendency differed from
  both terrain-following and constant-depth diffusion.
  `Diffusion tests <../../Tests/MoanaParity/t3dmix_iso_port_notes.md>`_.
* **Atmospheric forcing.** ROMS writes the forcing for the next step. REMORA
  writes the forcing it just applied. After shifting the records by one 100 s
  step, shortwave radiation and rain agreed exactly.
* **Hot starts.** A file with three distinct temperature records verified
  selection of the first, middle and last records.
  `Hot-start test <../../Tests/MoanaParity/hotstart_notes.md>`_.
* **Time averages.** For nine fields, five-step averages agreed with means
  calculated from the history files to within 2.03 × 10⁻¹⁶ relative error.
  `Averaging tests <../../Tests/MoanaParity/avg_port_notes.md>`_.

Model sources
-------------

`Moana hindcast configuration <https://github.com/joaometocean/moana_hindcast>`_ ·
`Upstream REMORA <https://github.com/seahorce-scidac/REMORA>`_ ·
`REMORA paper <https://doi.org/10.21105/joss.07958>`_ ·
`License <../../LICENSE.md>`_

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative
works, and perform publicly and display publicly, and to permit others to do so.
