# Full 3D VTK wind input

QES-Winds accepts a single binary legacy VTK `STRUCTURED_GRID` through
`metParams/windVTK`. All three point-data `wind_vectors` components are imported
at their physical 3D locations. The existing building/canopy parameterisations
then operate on that field, followed by the existing mass-conservation solve.

## Ghost Town example

Keep the runnable XML and launcher together in the labelled
[`GhostTown` run folder](../../Capstone/QES_RunScripts/GhostTown). This follows
the `QES_RunScripts/ClaireValley` layout. Run the GPU case with:

```sh
/home/thalia/Capstone/QES_RunScripts/GhostTown/RunCodes_GhostTown_3D_Winds.sh
```

Pass `1` for the serial CPU solver; the default is CUDA solver `3`. Outputs and
the run log go to `Capstone/QES_Data/GhostTownVTK` with prefix
`ghosttown_full3d_50m`. The launcher includes `-w` to save the staggered workspace
file with the original imported field and masks. The visualisation file alone
contains cell-centred corrected velocities.

The example uses 200 × 160 × 125 cells of 50 × 50 × 12 m, with physical origin
(636700, 4781600, 1770) m. Its top is 3270 m, below the source top of about
3345.528 m. The horizontal footprint is deliberately cropped inside both source
VTK coverage and the TIFF's pixel-centre interpolation support. This coarse
example demonstrates import and mass conservation; it does not resolve individual
Ghost Town buildings. No buildings are supplied in this example.

The file paths in the example are machine-specific. The timestamp is an example
label, **not a timestamp recovered from the VTK**. The `m/s` velocity units and
EPSG:32612 are explicit assertions in the XML. Confirm those assertions when
using a different source file.

## XML options

Place exactly one `windVTK` under `metParams`:

```xml
<metParams>
  <windVTK>
    <file>wind.vtk</file>
    <terrainFile>terrain.tif</terrainFile>
    <epsg>32612</epsg>
    <velocityUnits>m/s</velocityUnits>
    <originX>636700</originX>
    <originY>4781600</originY>
    <originZ>1770</originZ>
    <timeStamp>2020-01-01T00:00:00</timeStamp>
    <siteZ0>0.1</siteZ0>
    <reciprocal>0</reciprocal>
    <terrainTolerance>0.1</terrainTolerance>
  </windVTK>
</metParams>
```

| Option | Meaning |
| --- | --- |
| `file` | Required VTK path; relative paths resolve from the XML directory. |
| `terrainFile` | Required single-band north-up reference elevation TIFF, in metres and the same vertical datum as VTK z. It must have a projected metre-based CRS and finite, non-nodata elevations. |
| `epsg` | Required explicit declaration of the VTK's projected CRS; must match the TIFF and any building shapefile. No automatic reprojection is performed. |
| `velocityUnits` | Required `m/s`. Components are Cartesian velocities along projected x, projected y and physical upward z. |
| `originX/Y/Z` | Required finite physical coordinates of local QES (0,0,0). Absolute height is `originZ + local z`; it is not height above local terrain. |
| `timeStamp` | Required one-timestep label in QES timestamp format. |
| `siteZ0` | Optional uniform roughness, default 0.1 m, positive. Used by existing parameterisations/turbulence; it does not generate an imported wind profile. |
| `reciprocal` | Optional inverse Monin–Obukhov length, default 0 m⁻¹. Used by existing canopy parameterisations. |
| `terrainTolerance` | Optional maximum VTK/TIFF terrain registration error, default 0.1 m. Do not increase it to conceal a datum or alignment error. |

Unknown/repeated options and mixed wind-input modes are errors. Only one time,
unrotated uniform QES grids, zero halos, uniform roughness, and stair-step
geometry (`meshTypeFlag=0`) are supported. Require `logLawFlag=0` and
`readCoefficientsFlag=0`, zero `DEMDistancex/y`, and no WRF coupling. Cut cells, vertical stretching, sensors, wind shapefiles,
WRF and HRRR cannot be combined with this mode. Use `terrainFile`, not a separate
`simulationParameters/DEM`, because the VTK bottom layer supplies the terrain.
Existing modes retain their original parsing and profile paths.

Optional simulation `UTMx/UTMy/UTMZone` must agree with the declared origin/CRS.
For building shapefiles, the legacy building reader requires nonzero originX/Y
that can be represented exactly as floats; otherwise this mode rejects the case
rather than introducing a hidden offset. Inline XML buildings use local x/y.
Buildings outside wind coverage also fail when their terrain heights are sampled.

## Physical interpolation

The reader supports big-endian `float` or `double` points/vectors, with x varying
fastest, then y, then z. Horizontal x/y coordinates must be separable and strictly
increasing; spacing may be nonuniform. Vertical columns must have strictly
increasing z. ASCII, arbitrary horizontally curvilinear grids, cell-data winds,
other vector names, extra data arrays and folded layers are rejected.

For a query `(x,y,z)`:

1. Locate the containing horizontal rectangle and its four bilinear weights.
2. Evaluate each source layer's **physical height** using those weights.
3. Locate the two surfaces that enclose z. Compute the vertical fraction from
   their physical heights at `(x,y)`.
4. Apply the eight resulting hexahedral weights to each Cartesian velocity
   component. This is trilinear isoparametric interpolation on the source mesh.

It reproduces affine vector fields in physical XYZ to floating-point precision,
including on sloping, nonuniform terrain-following layers. It preserves nonzero
vertical velocity. There is no height collapse, logarithmic profile generation,
nearest-neighbour fallback, or extrapolation. The remap itself is not a
flux-conservative remap; QES's subsequent solve enforces its discrete continuity
constraint.

For `Nx,Ny,Nz` configured cells, local physical face locations are:

| Component | Location | Physical indices |
| --- | --- | --- |
| u | `(i dx, (j+½)dy, z[k])` | i=0…Nx, j=0…Ny−1, k=1…Nz |
| v | `((i+½)dx, j dy, z[k])` | i=0…Nx−1, j=0…Ny, k=1…Nz |
| w | `((i+½)dx, (j+½)dy, z_face[k])` | i=0…Nx−1, j=0…Ny−1, k=1…Nz+1 |

Every physical air face must lie within the source. The complete horizontal
footprint and top surface are checked before terrain construction. Missing
coverage, non-finite source data, inconsistent array sizes, truncated data,
invalid reference terrain and declared CRS mismatches produce errors.

## Terrain, solids and boundaries

The VTK bottom layer is the continuous terrain. Its heights are checked against
bilinear interpolation of the reference TIFF **at pixel centres**, at all source
nodes with TIFF interpolation support. At least four shared nodes must exist.
The entire QES footprint must be inside that support. Outer half-pixel regions
are never extrapolated. Ghost Town has 3,752 checked nodes and 124 outer nodes
without TIFF interpolation support.

QES terrain heights and building base-height queries use the same continuous
VTK bottom surface, translated by originZ. There is no automatic minimum-height
subtraction. A triangle mesh is also constructed for existing terrain-mesh
consumers; that mesh approximates the surface at QES horizontal grid vertices.

For this mode, stair-step solids conservatively cover every cell intersecting
the continuous terrain: the maximum terrain over each cell footprint is used,
including source-grid break lines, and the solid stack is rounded up to a QES
vertical face. Faces adjoining terrain/building solids have zero velocity. This
can remove near-ground detail, especially with coarse horizontal/vertical cells,
but prevents the importer or solver from flowing through terrain. Refinement
reduces this stair-step geometry approximation.

QES allocates all components to one padded array. The k=0 below-ground ghost
layer is zero. Unused tangential last rows/columns copy the last physical row or
column; upper u/v padding copies the last physical u/v layer. **The upper w face
is physical and is sampled at its actual height**, not copied. Padding is updated
after parameterisations and after the solve. It is not extrapolated source data.

The existing solver's prescribed outer horizontal cell ring and top cell layer
are retained. Their divergence is reported separately and is not guaranteed to
vanish. The solved interior has a separate maximum/RMS divergence diagnostic.
A converged interior solve must not be described as a divergence-free result over
the entire saved array. The imported field is an initial condition, not a field
that QES promises to leave unchanged at all boundaries or near buildings.

## Auditable output and diagnostics

`<basename>_windsWk.nc` includes:

- `u0_imported/v0_imported/w0_imported`: immutable remapped initial field, after
  solid masking and padding, **before** parameterisations.
- `u0/v0/w0`: the usual QES initial arrays, after parameterisations and wall
  constraints, used as input to the solver.
- `u/v/w`: final mass-corrected velocities.
- `vtk_u_mask/vtk_v_mask/vtk_w_mask`: 0=physical solid/ground face, 1=interpolated
  physical air face, 2=unused padding/below-ground ghost. Exclude mask=2 from
  physical statistics. No raw source samples inside solids are claimed.
- Global attributes recording source paths, EPSG, double-precision physical
  origins, the CRS/units assertion and staggered-coordinate conventions.

The common array dimensions are `[t,z_face,y_face,x_face]` for storage only.
Use the component's staggered coordinates from the table above. In this new
mode `x_face=i dx` and `y_face=j dy`; `x/y/z` contain cell-centre coordinates.
Add the saved origins to obtain projected x/y and absolute z. The legacy modes'
workspace coordinate output is unchanged.

## MATLAB inspection

`Capstone/QES Matlab Scripts/Slice_optimized.m` can show the final corrected
field, or the imported field before parameterisations and the mass-conservation
solve. Set its editable `dataFolder` and `fileName` variables. The script
selects imported mode when the filename contains `windsWk` and final mode when
it contains `windsOut`. It then reads
`u0_imported/v0_imported/w0` from the workspace NetCDF file, then averages the
staggered faces to cell centres for the slice plot. A filename containing
`windsOut` reads the visualisation file's corrected cell-centre field. The
imported view is
an audit view; its plotted lines are still 2D projections and may omit lines
when through-plane flow dominates.

For the original WindNinja source without QES interpolation, run
`Capstone/QES Matlab Scripts/Inspect_GhostTown_VTK.m`. It reads the big-endian
legacy VTK points and vectors and plots the selected terrain-following y layer
with physical z coordinates and quiver vectors. It does not re-grid the source
field, so it does not introduce a second interpolation step for visualization.

The run log reports import counts/ranges, maximum double-to-float rounding
error, per-component parameterisation max/RMS changes, solver iteration count,
maximum update, equation residual, solved-cell divergence and all-air/boundary
divergence. Rounding error is not an estimate of the unknown physical error in
WindNinja's field. The manufactured test supplies a known interpolation truth.

## Tests and measured results

```sh
cmake -S . -B build -DENABLE_WIND_VTK_TESTS=ON \
  -DENABLE_WIND_SHAPEFILE_TESTS=ON \
  -DWIND_VTK_GHOSTTOWN_XML="$PWD/../Capstone/QES_RunScripts/GhostTown/winds_WindNinja_GhostTown_3D.xml"
cmake --build build --target qesWinds windVTKTests windShapefileTests -j 6
ctest --test-dir build --output-on-failure
```

Omit `WIND_VTK_GHOSTTOWN_XML` (or set it to an empty string) on machines without
these private GIS files. The manufactured test creates its own VTK/TIFF/XML
fixtures, checks physical affine interpolation and all face locations, nonzero w,
ghosts, solid masking, a building parameterisation case, saved snapshots, solver
convergence, and rejection of invalid coordinates/data/coverage/CRS/options.
The Ghost Town integration test executes the complete import/parameterisation/
CPU-solve/workspace-output pipeline and checks the saved arrays and solid faces.

For a saved integration run, optionally choosing GPU solver 3:

```sh
./build/tests/windVTK/windVTKTests --integration \
  /home/thalia/Capstone/QES_RunScripts/GhostTown/winds_WindNinja_GhostTown_3D.xml \
  /tmp/ghosttown_cpu_windsWk.nc 1
```

See [measured test results](WIND_VTK_TEST_RESULTS.md) for errors, parameterisation
changes, CPU/GPU convergence and the boundary-divergence limitation.
