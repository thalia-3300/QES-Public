# Physical 3D VTK input: measured results

Measured on 2026-10-07 using the Release build, GCC 13, CPU serial solver and
CUDA global-memory solver. Input: [`winds_WindNinja_GhostTown_3D.xml`](../../Capstone/QES_RunScripts/GhostTown/winds_WindNinja_GhostTown_3D.xml).
See [input documentation](WIND_VTK_INPUT.md) for constraints and reproduction commands.

## Existing paths inspected

- `WRFInput` reads the specialised WRF/SFIRE arrays; `WindProfilerWRF` maps
  horizontal components through its existing profile path. That path does not
  provide a general legacy-VTK importer preserving physical 3D w.
- `WindShapefileInput`, `WindShapefileGrid` and `WindProfilerShapefile` import
  horizontal point winds at one reference height and generate a profile.
- `WINDSGeneralData::applyWindProfile` creates u0/v0/w0;
  `applyParametrizations` changes those arrays before the solver consumes them.
- `Solver::buildRHS` and `recoverVelocity` already handle nonzero w. Their
  continuity equation and physical-face correction are retained. New hooks
  refresh unused padding and report all-air/boundary divergence for VTK mode.

## Manufactured interpolation and buildings

A 4 × 4 × 4 binary VTK fixture has bilinear sloping terrain and nonuniform
terrain-following vertical layers. Its three velocities are known affine
functions of physical x/y/z, with nonzero vertical velocity. The QES grid has
12 × 10 × 18 cells of 2 × 2 × 5 m.

| Check | Measured result |
| --- | ---: |
| Maximum double-precision interpolation error, 200 interior points × 3 components | 1.21529 × 10⁻¹¹ m/s |
| Maximum error at physical staggered faces after float storage | 2.28882 × 10⁻⁷ m/s |
| Nonzero imported w faces, no building / one building | 2,034 / 2,022 |
| No-building parameterisation changes | Exactly zero |
| No-building convergence | 96 iterations, tolerance 10⁻⁹ |
| One-building convergence | 97 iterations, tolerance 10⁻⁹ |
| No-building solved-cell maximum divergence | 1.37091 × 10⁻⁷ s⁻¹ |
| One-building solved-cell maximum divergence | 1.54972 × 10⁻⁷ s⁻¹ |

For the 4 × 4 × 15 m test building, existing upwind, wake, canyon, sidewall and
rooftop code executes. Changes from import to solver input, excluding padding:

| Component | Changed physical faces | Maximum change | RMS change |
| --- | ---: | ---: | ---: |
| u | 33 / 2,340 | 4.98570752 m/s | 0.312049740 m/s |
| v | 29 / 2,376 | 1.00788608 m/s | 0.066693042 m/s |
| w | 27 / 2,280 | 0.528999984 m/s | 0.053353517 m/s |

The test checks immutable imported snapshots, exact NetCDF round trips,
physical top w, corrected ghost copies, final solid masks, and rejected malformed
inputs. These are numerical implementation tests; they do not validate the
physical accuracy of the building parameterisations.

## Ghost Town import

Source: binary legacy `STRUCTURED_GRID`, 68 × 57 × 20, point-data
`wind_vectors`, three big-endian doubles per vector. The declared projected CRS
is EPSG:32612; the VTK itself does not encode a CRS or velocity units.

Reference terrain agreement: **5.52973 × 10⁻¹⁰ m maximum error** at 3,752 VTK
bottom-layer nodes with TIFF pixel-centre interpolation support. The other 124
outer nodes have no such support; the example QES footprint is cropped inside
both sources.

| Component | Sampled air faces | Solid/terrain faces | Imported range including zero solids | Maximum float rounding |
| --- | ---: | ---: | --- | ---: |
| u | 93,418 | 59,582 | −1.55053475 … 0.233978413 m/s | 5.01736 × 10⁻⁸ m/s |
| v | 93,895 | 59,855 | −0.294480834 … 2.30106067 m/s | 6.43023 × 10⁻⁸ m/s |
| w | 92,567 | 59,433 | −0.588798505 … 1.04269030 m/s | 5.84278 × 10⁻⁸ m/s |

All 92,567 sampled physical w faces remain nonzero in the imported field.
Ghost Town's example supplies no buildings/canopies, so its parameterisation
max/RMS deltas are **exactly zero** for all three components. The separate
manufactured building case above exercises nonzero parameterisation changes.

The terrain-registration and float-rounding numbers are not estimates of the
unknown interpolation error relative to the true Ghost Town wind. Analytic
interpolation error is measured with the manufactured field above.

## Ghost Town mass-conservation solve

50 × 40 × 75 QES cells, spacing 200 × 200 × 20 m, tolerance 10⁻⁸, cap 30,000
iterations. There are 82,654 solved air cells and 92,567 total air cells.

| Diagnostic | CPU serial | GPU global-memory |
| --- | ---: | ---: |
| Status | CONVERGED | CONVERGED |
| Iterations | 1,858 | 1,850 |
| Maximum potential update | 9.92566 × 10⁻⁹ | 9.99991 × 10⁻⁹ |
| Maximum equation residual | 1.05351 × 10⁻¹¹ | 2.29297 × 10⁻¹¹ |
| Solved-cell maximum divergence | 3.24100 × 10⁻⁹ s⁻¹ | 3.24100 × 10⁻⁹ s⁻¹ |
| Solved-cell RMS divergence | 3.91281 × 10⁻¹⁰ s⁻¹ | 3.91323 × 10⁻¹⁰ s⁻¹ |
| Net flux over solved cells | 0.0259535 m³/s | 0.0308011 m³/s |
| All-air maximum divergence | 0.0440675 s⁻¹ | 0.0440675 s⁻¹ |
| All-air RMS divergence | 0.000549101 s⁻¹ | 0.000549101 s⁻¹ |

**The large remaining divergence is in QES's prescribed outer ring/top layer,
which is excluded from the Poisson solve.** The field is mass-converged in the
solved region; the entire saved volume is not divergence-free. This retained
solver boundary behaviour is reported explicitly. Changing it would require a
separate boundary-condition design and solver validation.

## Interpreting vertical MATLAB slices

The 40-cell y domain has a prescribed outer ring at MATLAB `j=40`. At that
slice, all 2,775 air cell centres have an out-of-plane `|v|` greater than the
plotted x–z magnitude `hypot(u,w)`. Thus x–z `streamslice` lines are projections,
not physical 3D paths. Near the terrain on this slice, median `|v|` is about
1.17 m/s while median `|u|` is about 0.017 m/s. Cell-centre `w` can be nonzero
just above the zero-flux stair-step surface; local upward flow does not by
itself indicate flow through a solid face.

Directly checking the staggered workspace arrays with their saved cell flags
found zero velocity at all 2,000 vertical air/terrain interfaces, including the
50 interfaces in the displayed `j=40` slice. Every internal u/v/w face next to
a flagged solid cell is also zero. The MATLAB viewer under `Capstone/QES Matlab
Scripts/Slice_optimized.m` now labels prescribed boundary slices and suppresses
projected streamlines where the omitted component dominates. The `x` key
switches to a y–z slice, which is more representative of this flow.

CPU/GPU `u0/v0/w0` match exactly. Maximum final CPU/GPU differences are
7.45058 × 10⁻⁹, 1.19209 × 10⁻⁷ and 1.49012 × 10⁻⁸ m/s for u/v/w.
The independent workspace mass checker also confirms zero velocities on solid
faces and the solved-cell divergence above.

## Regression status and saved artifacts

All four CTest tests pass:

- `wind_shapefile_regression`
- `wind_solver_regression`
- `wind_vtk_manufactured`
- `wind_vtk_ghosttown`

A separate legacy single-sensor 20³-cell CLI smoke test converges in one
iteration with zero divergence and finite nonzero wind output. WRF integration
was inspected but not run with a WRF dataset in this change.

Saved run artifacts on this machine:

- `/home/thalia/Capstone/QES_Data/GhostTownVTK/ghosttown_cpu_windsWk.nc`
- `/home/thalia/Capstone/QES_Data/GhostTownVTK/ghosttown_gpu_windsWk.nc`
- `/home/thalia/Capstone/QES_Data/GhostTownVTK/ghosttown_gpu_windsOut.nc`
- Corresponding `ghosttown_cpu.log` and `ghosttown_gpu.log` in the same directory.

Both workspace files contain the imported snapshots, parameterised solver input,
final velocities, import masks and physical-origin metadata.

Source SHA-256 hashes:

```text
Ghosttown.vtk
74c8e6f2bf608be1607a2b8e68f9bdab193f94f8da108b76fae178c109832993
Ghosttown.tiff.tif
f0fb766a27c0bbd0a2b6abf68787e1f390c02f21af32611e6fba3690a2b273d0
```
