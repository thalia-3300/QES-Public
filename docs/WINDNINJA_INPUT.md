# Direct WindNinja point-shapefile input (version 1)

QES can initialize its wind field directly from one regular WindNinja point
shapefile. The importer uses the existing GDAL dependency, not a grid of XML
sensors. Rebuild QES once after this source change; thereafter select a different
shapefile in XML without rebuilding. Existing sensor/WRF/HRRR configurations
continue to use their original paths.

## Build and run with your GPU

From the QES-Public root, using your existing CUDA-enabled build directory:

```bash
cmake -S . -B build -DENABLE_WIND_SHAPEFILE_TESTS=ON
cmake --build build --target qesWinds windShapefileTests -j4
ctest --test-dir build -R 'wind_(shapefile|solver)_regression' --output-on-failure
./build/tests/windShapefile/windShapefileTests --solver-tests-gpu
mkdir -p output
./build/qesWinds/qesWinds -q docs/examples/winds_WindNinja_Testing.xml \
  -s 3 -t -w -o output/windninja_testing
```

`-s 3` selects the existing global-memory GPU solver. `-s 4` selects the existing
shared-memory GPU solver; `-s 1` selects CPU. Reading and interpolating the
shapefile is a one-time host operation for either solver: it fills the same
`u0/v0/w0` arrays. The mass-conservation solver includes convergence fixes shared
by all input types: double-precision correction values, corrected GPU indexing
and reduction, and consistent final face-velocity updates. SOR iterations still
run on the GPU; RHS construction, final velocity recovery and diagnostics run
on the host. Outputs remain single precision. Internal doubles increase working
memory and may reduce GPU arithmetic throughput; performance has not been
benchmarked on large production domains. No new dependency is required.
`-t` produces turbulence output; `-w` also saves workspace fields including the
initial wind.

## XML configuration

Use this **instead of** `sensor` or `sensorName` inside `metParams`:

```xml
<metParams>
  <z0_domain_flag>0</z0_domain_flag>
  <windShapefile>
    <file>wind-output.shp</file>
    <speedField>speed</speedField>
    <directionField>dir</directionField>
    <speedUnits>m/s</speedUnits>
    <directionUnits>degrees</directionUnits>
    <originX>374610</originX>
    <originY>5835000</originY>
    <height>10</height>
    <siteZ0>0.1</siteZ0>
    <reciprocal>0</reciprocal>
    <boundaryLayerFlag>1</boundaryLayerFlag>
    <timeStamp>2026-09-10T21:09:00</timeStamp>
  </windShapefile>
</metParams>
```

The numerical settings above are **test assumptions**, not recovered simulation
metadata. Set the actual WindNinja output height before production use. The `7m`
and `20m` filename suffixes refer to output resolution, not wind measurement height.

| Option | Meaning |
| --- | --- |
| `file` | Required `.shp` path; keep `.shx`, `.dbf` and `.prj` companions. Relative paths are resolved from the XML directory. |
| `speedField`, `directionField` | Numeric attributes; defaults `speed`, `dir`. Do not use WindNinja's display-angle fields `AM_dir` or `QGIS_dir`. |
| `speedUnits` | Required: `m/s`, `mps`, `mph`, `km/h`, `kph`, `knots`, or `kts`. Converted to m/s. |
| `directionUnits` | `degrees` (default, 0–360) or `radians` (0–2π). Meteorological **FROM**, clockwise from projected grid north. |
| `originX`, `originY` | Required double-precision projected coordinates of QES local **(0,0)**, in the shapefile's CRS. Not the first QES cell centre. |
| `height` | Required positive reference/output wind height in **metres above ground**. No automatic vegetation/displacement-height correction. |
| `siteZ0` | Required uniform aerodynamic roughness length, metres, positive. Also used by QES turbulence/canopy consumers. |
| `reciprocal` | Required 1/L, m⁻¹: zero neutral, positive stable, negative unstable. Not inferred from WindNinja's stability filename. |
| `boundaryLayerFlag` | `1` (default): QES-style stability-corrected logarithmic profile; `0`: vertically uniform above ground. Other profiles are not supported in v1. |
| `timeStamp` | Required single timestamp, `YYYY-MM-DDThh:mm:ss`; set its intended time basis explicitly. |

Unknown or repeated options are rejected. `simulationParameters/totalTimeIncrements`
must be omitted, 0, or 1. There is only one wind record, regardless of timestamp.

## Coordinates, coverage and missing cells

The source must form a complete axis-aligned rectangular lattice with constant
spacing on each axis (x/y spacings may differ), at least 2×2 points, and exactly
one point at each location. Feature order does not matter. Coordinate tolerance
is 0.00001 m. Holes, duplicates, irregular spacing, non-point geometry, missing
values, negative speeds and invalid directions fail with an error.

The `.prj` must describe a projected, metre-based CRS. QES x/y axes must align with
its easting/northing axes: `domainRotation=0`, `halo_x=halo_y=0`. There is no automatic
reprojection or true-north/grid-north angle conversion. Reprojecting point
coordinates alone may destroy regular spacing and does not rotate wind vectors.

For a QES domain `Nx Ny Nz` with cell sizes `dx dy dz`, the **whole** rectangle
`[originX, originX+Nx*dx] × [originY, originY+Ny*dy]` must lie inside the source
point bounds. Raster outer edges often extend half a cell beyond point centres:
they do **not** count as wind coverage. Crop the QES domain or generate a larger
WindNinja grid; this version never extrapolates or silently falls back to Barnes.

With a DEM, its CRS must match the wind grid, and the wind origin must exactly
match the DEM crop selected by QES (`originFlag=0` with `DEMDistancex/y`, or
`originFlag=1` with `UTMx/y`). The crop must fall on DEM pixel boundaries. Use the
actual elevation DEM, not a WindNinja speed/direction TIFF. With building
shapefiles, their CRS must also match, and simulation `UTMx/UTMy` must match wind
`originX/originY` so QES's existing building importer uses the same origin.
XML rectangular buildings remain in QES local coordinates.

Source angles/speeds are converted to `u=-speed*sin(direction)` and
`v=-speed*cos(direction)` **before** bilinear interpolation. This avoids the
359°/1° wrap-around error and correctly permits cancellation of opposing winds.
Interpolation is evaluated at QES's staggered u/v face positions. Only internal
padding rows/columns copy their adjacent physical values; no physical boundary
is filled from outside coverage. Initial `w=0`.

Each horizontal face receives a vertical profile normalized at `height`.
Neutral scaling is `log(1+z/siteZ0) / log(1+height/siteZ0)`; stable/unstable
corrections and the free-convection limit follow QES's existing log-profile
convention. Height above ground is measured from adjacent-cell-averaged terrain
at each face. Below-ground and bottom-ghost velocities are zero. The QES grid
must currently have uniform vertical spacing and at least two cells per axis.
Buildings, wakes, terrain masking and mass conservation are applied afterward.

This is a diagnostic initialization, not import of WindNinja's full 3-D flow,
turbulence or stratification. It does not guarantee closer agreement with
observations. Existing QES parameterizations and its solve can change the input
wind substantially, and applying terrain influences again needs scientific review.

## Testing and Clare examples

- [Testing](examples/winds_WindNinja_Testing.xml): 142×142 source points, spacing
  about 6.56 m, WGS84 / UTM 55S (EPSG:32755); flat 880×880 m QES domain inside coverage.
- [Clare](examples/winds_WindNinja_Clare.xml): 51×51 source points, spacing 20 m,
  GDA2020 / MGA 54 (EPSG:7854); 960×960 m QES domain with the supplied Clare DEM.

These example paths assume `Capstone` and `QES-Public` are sibling directories.
They use m/s based on the supplied settings, but **10 m AGL, z0=0.1 m and neutral
QES stability are explicitly unverified test assumptions**, including for the
Clare file named `non_neutral_stability`. If WindNinja output height was specified
above vegetation, determine the effective above-ground reference height first.

The standalone CTest regression generates synthetic shapefiles and tests affine
u/v interpolation, wrap-around, cardinal directions, calm wind, units, vertical
profiles, staggered positions, terrain-relative heights, padding, source defects,
CRS and coverage errors. It needs no Catch2 download.

Run the private-data smoke suite using the actual paths to your datasets:

```bash
python3 tests/windShapefile/run_real_cases.py --build build --solver 3 --require-convergence \
  --testing "/path/to/Testing/wind.shp" \
  --clare "/path/to/Clare/wind.shp" \
  --clare-dem "/path/to/Clare/elevation.tif"
```

The suite samples both source grids, runs Testing flat/building and Clare
flat/terrain cases plus the original sensor path, checks finite wind/turbulence
outputs, and exercises rejected XML configurations. With `--require-convergence`,
it additionally requires a converged solver log and independently checks saved
face velocities for mass conservation and zero flow through solid faces.
Each run keeps XML, logs
and NetCDF results in a fresh `/tmp/qes-windninja-s...` directory. Repeat with
`--solver 1` for CPU or `--solver 4` for shared-memory GPU. Compare CPU/GPU wind
outputs using `windShapefileTests --compare first_windsOut.nc second_windsOut.nc`.
Compare complete initial fields using `--compare-initial first_windsWk.nc
second_windsWk.nc`. See [the implementation test record](WINDNINJA_TEST_RESULTS.md)
for measured convergence and CPU/GPU comparisons. The manufactured-solution
regression tests a known flow with unequal x/y spacing and nonzero initial w;
the GPU variant exercises both GPU modes as well as CPU solvers.

Cases default to a 1,000-iteration cap (override with `--iterations`). All four
current WindNinja regression cases converge below that cap at tolerance 1e-7.
The solver now prints `CONVERGED` or `NOT CONVERGED`, the maximum potential update,
an independent equation residual, and direct final-velocity divergence. A finite
nonconverged run still writes output, but prints a warning: do not use that output
as a converged result. The convergence flag requires both algebraic measures to
meet the configured tolerance; the separately reported velocity divergence has
units s^-1 and is limited by single-precision output. That diagnostic covers
solved regular/stair-step air cells, excluding prescribed outer-ring/top cells
and cut cells (which require area/volume-weighted checks).
Numerical convergence is **not** physical validation. Confirm source metadata
and compare against measurements or a suitable reference before scientific use.

Not implemented: irregular/rotated source grids, multiple files/times/heights,
automatic unit/height inference, spatial roughness/stability, automatic CRS
conversion, or assimilation/blending with sensors/WRF/HRRR.

For WindNinja binary legacy VTK volumes with u/v/w at multiple physical heights, see [full 3D VTK input](WIND_VTK_INPUT.md). That input preserves vertical velocity and does not use this shapefile mode's single-height vertical profile.
