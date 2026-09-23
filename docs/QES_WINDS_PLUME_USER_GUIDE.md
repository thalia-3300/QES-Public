# QES-Winds and QES-Plume user guide

This guide explains how to build and run QES-Winds, QES-Turb, and QES-Plume, how their XML files are structured, and what files they produce. It is based on the current `main` branch at QES v2.3.1 (commit `dced7ec571487be7538eca17bfc25bb3262fb6d0`, 2025-12-11), not on older example syntax that remains in parts of the repository.

The two ready-to-run inputs used throughout this guide are:

- [`examples/winds_minimal.xml`](examples/winds_minimal.xml)
- [`examples/plume_minimal.xml`](examples/plume_minimal.xml)

Both examples were exercised end to end. The primary commands below assume an accessible NVIDIA GPU and a CUDA-enabled build. Together they produce a wind visualisation file, the wind and turbulence files required by standalone Plume, and a concentration file.

## 1. What each executable does

QES contains three useful entry points for a Winds/Plume workflow:

| Executable | Purpose | Typical use |
|---|---|---|
| `qesWinds` | Builds an initial wind field, applies terrain/building/vegetation parameterisations, solves the mass-conservation problem, and optionally computes QES-Turb fields. | Generate wind and turbulence files, inspect Winds results, or prepare inputs for standalone Plume. |
| `qesPlume` | Reads precomputed QES-Winds and QES-Turb NetCDF files and advects Lagrangian particles. | Re-run dispersion with different sources or collection grids without recomputing Winds. |
| `qes` | Runs Winds, Turb, and optionally Plume in one process. | A coupled one-command run where intermediate NetCDF input files are not needed by Plume. |

The normal data flow is:

```text
Winds XML
   |
   v
qesWinds --turbcomp --wkout
   |------------------------------|
   v                              v
<base>_windsWk.nc          <base>_turbOut.nc
   |                              |
   |------------ qesPlume --------|
                    |
                    v
          <base>_<tag>_plumeOut.nc
```

`qesWinds` always creates `<base>_windsOut.nc` when `-o` is supplied. The extra `-w` and `-t` switches create the two files that standalone Plume requires.

## 2. Requirements and building

QES requires C++17, CMake, Boost Program Options and Date Time, NetCDF C/C++, and GDAL. CUDA is optional at compile time. A CUDA-enabled binary checks for a usable NVIDIA device at program startup, even before displaying `--help`.

On Debian/Ubuntu, the project README lists these packages:

```bash
sudo apt install \
  libgdal-dev \
  libnetcdf-c++4-dev \
  libnetcdf-cxx-legacy-dev \
  libnetcdf-dev \
  netcdf-bin \
  libboost-all-dev \
  cmake \
  cmake-curses-gui
```

### CUDA build (recommended)

For a CUDA build:

```bash
cd /path/to/QES-Public
cmake -S . -B build-gpu \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_ARCHITECTURES=native
cmake --build build-gpu --parallel
```

The repository README states that GPU acceleration requires NVIDIA compute capability 7.0 or newer and has been tested with CUDA 11.8. The current build files also contain work for CUDA 12.8. If OptiX is installed outside a standard location, add `-DOptiX_INSTALL_DIR=/path/to/optix`.

The main binaries will be:

```text
build-gpu/qesWinds/qesWinds
build-gpu/qesPlume/qesPlume
build-gpu/qes/qes
```

### CPU-only fallback

For a machine without an accessible NVIDIA GPU:

```bash
cd /path/to/QES-Public
cmake -S . -B build-cpu \
  -DCMAKE_BUILD_TYPE=Release \
  -DHAS_CUDA_SUPPORT=OFF \
  -DHAS_OPTIX_SUPPORT=OFF
cmake --build build-cpu --parallel
```

In a CPU-only build, any `--solvetype` value is forced to solver type `1`.

### Verify the build

```bash
build-gpu/qesWinds/qesWinds --help
build-gpu/qesPlume/qesPlume --help
build-gpu/qes/qes --help
```

If a CUDA-enabled build reports `cudaGetDeviceCount` or `no CUDA-capable device is detected`, use a GPU node or make a CPU-only build. Selecting `-s 1` does not bypass the startup GPU check in a CUDA-enabled binary.

## 3. Quick start: recommended two-stage workflow

Run all commands from the repository root. Output directories must already exist.

### Step 1: compute Winds and Turb

```bash
mkdir -p results

build-gpu/qesWinds/qesWinds \
  -q docs/examples/winds_minimal.xml \
  -o results/minimal \
  -w \
  -t \
  -s 3
```

This creates:

```text
results/minimal_windsOut.nc   cell-centred wind fields for analysis/visualisation
results/minimal_windsWk.nc    full Winds workspace consumed by qesPlume
results/minimal_turbOut.nc    turbulence fields consumed by qesPlume
```

Important: `-w` and `-t` are Boolean switches. Write `-w -t`, not `-w 1 -t 1`.

### Step 2: run standalone Plume

The project-basename form is the least error-prone:

```bash
build-gpu/qesPlume/qesPlume \
  -p docs/examples/plume_minimal.xml \
  -f results/minimal
```

`-f results/minimal` expands to:

```text
Winds input:       results/minimal_windsWk.nc
Turbulence input:  results/minimal_turbOut.nc
Output basename:   results/minimal
```

The example has a particle tag of `tracer`, so its concentration output is:

```text
results/minimal_tracer_plumeOut.nc
```

You can instead name each input explicitly:

```bash
build-gpu/qesPlume/qesPlume \
  -p docs/examples/plume_minimal.xml \
  -w results/minimal_windsWk.nc \
  -t results/minimal_turbOut.nc \
  -o results/plume_run_02
```

## 4. One-command integrated workflow

The combined executable keeps Winds and Turb in memory and passes them directly to Plume:

```bash
build-gpu/qes/qes \
  -q docs/examples/winds_minimal.xml \
  -p docs/examples/plume_minimal.xml \
  -o results/integrated \
  -s 3
```

Supplying `-p` automatically enables the turbulence computation required by Plume. Add `-t` if you also want `<base>_turbOut.nc`, and add `-w` if you also want `<base>_windsWk.nc` for later standalone runs.

The two-stage workflow is generally easier to debug and is more efficient when testing several Plume configurations against one wind field.

## 5. Command-line reference

Both short and long forms are accepted, for example `-q file.xml` and `--qesWindsParamFile file.xml`.

### `qesWinds`

| Option | Meaning |
|---|---|
| `-h`, `--help` | Print help. |
| `-v`, `--verbose` | Enable verbose QES logging. |
| `-q`, `--qesWindsParamFile FILE` | Winds XML input. Required. |
| `-o`, `--outbasename BASE` | Output path and basename. If omitted, file output is disabled. |
| `-s`, `--solvetype N` | Solver: `1` CPU, `2` dynamic-parallelism request (currently redirected to global-memory GPU), `3` global-memory GPU, `4` shared-memory GPU. Default `1`. |
| `-w`, `--wkout` | Write `<base>_windsWk.nc`, required for standalone Plume. |
| `-t`, `--turbcomp` | Compute turbulence and write `<base>_turbOut.nc`. Requires `<turbParams>` in the Winds XML. |
| `-m`, `--terrainout` | Write `<base>_terrainOut.obj`. Requires a DEM-backed terrain. |
| `-x`, `--windsolveroff` | Registered as “wind solver off”, but the current main loop still constructs and runs a solver. Do not rely on this switch in v2.3.1. |
| `-f`, `--firemode` | Registered and reported, but the normal WRF write-back path is controlled by `<WRFCoupling>` in the XML. |

### `qesPlume`

| Option | Meaning |
|---|---|
| `-h`, `--help` | Print help. |
| `-v`, `--verbose` | Enable verbose QES logging. |
| `-p`, `--qesPlumeParamFile FILE` | Plume XML input. Required. |
| `-f`, `--projectQESFiles BASE` | Use `<base>_windsWk.nc` and `<base>_turbOut.nc`, and write outputs using the same base. |
| `-w`, `--inputWINDSFile FILE` | Explicit Winds workspace input. Required when `-f` is not used. |
| `-t`, `--inputTURBFile FILE` | Explicit turbulence input. Required when `-f` is not used. |
| `-o`, `--outputbasename BASE` | Explicit output basename. Required when `-f` is not used. |
| `-l`, `--particleOutput` | Write periodic per-particle CSV snapshots. This requires `<particleOutputParameters>` in the Plume XML. |

### `qes`

The integrated executable accepts `-q`, optional `-p`, `-o`, `-s`, `-t`, `-w`, `-m`, `-l`, `-v`, and `-h`. Its `-x` switch has the same current limitation as `qesWinds`. The registered `-z/--visuout` switch is not consulted; a Winds visualisation file is already enabled whenever an output basename is supplied.

## 6. Paths, coordinates, and time

### Input paths

Relative data paths inside a Winds XML, such as `<DEM>../GISFiles/terrain.tif</DEM>`, are resolved relative to the XML file’s directory. The XML file itself may be passed as an absolute or working-directory-relative path.

Output basenames are interpreted relative to the process working directory. QES creates files but does not create missing parent directories.

### Wind direction

The input uses meteorological wind direction in degrees. The existing examples use values such as `270` for a wind arriving from the west and directed toward the east.

### Domain dimensions

`<domain>nx ny nz</domain>` describes base cell counts. QES internally adds staggered-grid faces and vertical ghost cells. `<cellSize>dx dy dz</cellSize>` is in metres.

For v2.3.1, use equal horizontal spacing (`dx == dy`). The current `Domain::getDomainSize()` implementation returns `dx` for both horizontal dimensions, which can make downstream calculations inconsistent for a rectangular horizontal cell.

### Time stamps

Use ISO-like timestamps in the form:

```xml
<timeStamp>2026-01-01T00:00:00</timeStamp>
```

`<timeEpoch>` is also accepted. If both are supplied, they must represent the same time. If neither is supplied, QES uses the current UTC time, which makes runs less reproducible.

## 7. QES-Winds XML

The root element must be exactly:

```xml
<QESWindsParameters>
  ...
</QESWindsParameters>
```

XML tags are case-sensitive. For example, the parser reads `UTMx` and `UTMy`; older files containing `UTMX` and `UTMY` do not set those values.

### 7.1 `simulationParameters`

This section is required by the root parser. In practice, `domain`, `cellSize`, and `totalTimeIncrements` should always be specified for a non-WRF case even though the low-level parser marks them optional.

| Tag | Required in practice | Meaning/default |
|---|---:|---|
| `domain` | Yes for non-WRF | Three integer base cell counts: `nx ny nz`. |
| `cellSize` | Yes for non-WRF | Three cell sizes in metres: `dx dy dz`. |
| `verticalStretching` | No | `0` uniform, `1` custom; default `0`. See current limitation below. |
| `dz_value` | No | Repeat this element for custom vertical cell sizes. |
| `totalTimeIncrements` | Yes | Number of Winds time records to process. |
| `logLawFlag` | No | `0` off, `1` apply the wall log-law iteration; default `0`. |
| `maxIterations` | No | Solver iteration limit; default `500`. |
| `tolerance` | No | Solver convergence tolerance; default `1e-9`. |
| `meshTypeFlag` | No | `0` stair-step, `1` cut-cell; default `0`. |
| `domainRotation` | No | Domain rotation relative to true north, in degrees; default `0`. |
| `originFlag` | No | `0` DEM/local origin, `1` UTM origin; default `0`. |
| `UTMx`, `UTMy` | Conditional | UTM origin coordinates in metres. |
| `UTMZone`, `UTMZoneLetter` | No | UTM zone metadata. |
| `DEMDistancex`, `DEMDistancey` | No | Offset into a DEM when `originFlag=0`; defaults `0`. |
| `halo_x`, `halo_y` | No | Halo distance added on each horizontal side, in metres; defaults `0`. |
| `heightFactor` | No | Scale applied to imported building heights; default `1`. |
| `readCoefficientsFlag` | No | `0` compute solver coefficients, `1` read them; default `0`. |
| `COEFF` | Conditional | Coefficient NetCDF path. |
| `DEM` | No | Terrain raster path. |
| `WRF` | No | WRF NetCDF path. If used without a DEM, WRF also supplies terrain/domain data. |
| `WRFSensorSample` | No | WRF sensor sampling interval; default `1`. |
| `WRFCoupling` | No | Enable iterative WRF write-back/re-read; default `false`. |

Current limitation: `verticalStretching` and `dz_value` are parsed, but the current `qes::Domain` construction path always calls the uniform-grid implementation. Treat custom vertical stretching as inactive in v2.3.1 unless that code path is updated and validated.

### 7.2 `metParams`, `sensor`, and `timeSeries`

Wind profiles can come from inline sensors, external sensor XML files, WRF, or HRRR.

```xml
<metParams>
  <z0_domain_flag>0</z0_domain_flag>
  <sensor>
    <site_coord_flag>1</site_coord_flag>
    <site_xcoord>10.0</site_xcoord>
    <site_ycoord>10.0</site_ycoord>
    <timeSeries>
      <timeStamp>2026-01-01T00:00:00</timeStamp>
      <boundaryLayerFlag>1</boundaryLayerFlag>
      <siteZ0>0.1</siteZ0>
      <reciprocal>0.0</reciprocal>
      <height>10.0</height>
      <speed>5.0</speed>
      <direction>270.0</direction>
    </timeSeries>
  </sensor>
</metParams>
```

`site_coord_flag` is required for every sensor:

| Value | Coordinates to provide |
|---:|---|
| `1` | `site_xcoord`, `site_ycoord` in local QES metres. |
| `2` | `site_UTM_x`, `site_UTM_y`, and `site_UTM_zone`. |
| `3` | `site_lon` and `site_lat`. |

Each `timeSeries` requires `siteZ0`, `reciprocal`, and one or more corresponding `height`, `speed`, and `direction` values. `boundaryLayerFlag` defaults to `1`. Optional `canopyHeight` and `attenuationCoefficient` values support canopy profiles.

Repeat `<timeSeries>` for multiple times and repeat `<sensor>` for multiple sites. The number and ordering of time records must be consistent with `totalTimeIncrements` and across sensors.

External sensor files can be listed with repeated `<sensorName>path.xml</sensorName>` entries under `metParams`.

### 7.3 `buildingsParams`

If this section is present, the following parameterisation controls are required by the parser:

| Tag | Meaning |
|---|---|
| `wallRoughness` | Building-wall roughness length. |
| `rooftopFlag` | `0` none, `1` log-profile, `2` vortex. |
| `upwindCavityFlag` | `0` none, `1` Röckle, `2` MVP, `3` HMVP. |
| `streetCanyonFlag` | `0` none, `1` canyon parameterisation. |
| `streetIntersectionFlag` | `0` off, `1` on. |
| `wakeFlag` | `0` none, `1` Röckle, `2` modified Röckle. Some old comments also mention `3`; validate before production use. |
| `highRiseFlag` | `0` off, `1` on. |
| `sidewallFlag` | `0` off, `1` on. |

Buildings can be inline rectangles:

```xml
<rectangularBuilding>
  <height>10.0</height>
  <baseHeight>0.0</baseHeight>
  <xStart>20.0</xStart>
  <yStart>20.0</yStart>
  <length>15.0</length>
  <width>10.0</width>
  <buildingRotation>0.0</buildingRotation>
</rectangularBuilding>
```

or inline polygons:

```xml
<QUICBuilding>
  <height>10.0</height>
  <baseHeight>0.0</baseHeight>
  <xVertex>20.0</xVertex>
  <yVertex>20.0</yVertex>
  <xVertex>35.0</xVertex>
  <yVertex>20.0</yVertex>
  <xVertex>35.0</xVertex>
  <yVertex>30.0</yVertex>
</QUICBuilding>
```

For shapefile input, use `SHPFile`, with optional `SHPBuildingLayer` (default `buildings`), `SHPHeightField` (default `H`), and `heightFactor` (default `1`). `numBuildings` and `numPolygonNodes` are parsed for compatibility but the vectors and shapefile features are the operative inputs.

### 7.4 `vegetationParams`

`wakeFlag` is optional. Repeat any of these canopy elements:

- `Homogeneous`: requires `attenuationCoefficient`, `height`, and `baseHeight`; define either rectangle fields (`xStart`, `yStart`, `length`, `width`, `canopyRotation`) or repeated polygon `xVertex`/`yVertex` values.
- `IsolatedTree`: requires `attenuationCoefficient`, `height`, `zMaxLAI`, `baseHeight`, `xCenter`, `yCenter`, and `width`.
- `Windbreak`: requires `height`, `understroyHeight` (the parser uses this misspelling), `xStart`, `yStart`, `length`, `width`, `canopyRotation`, and `opticalPorosity`; `baseHeight`, `fenceModel`, and `fetch` are optional.
- `ROC`: a specialised row-of-canopy model with additional required row geometry, leaf-area, porosity, and TKE values.

Shapefile vegetation input uses `SHPFile` and optional `SHPTreeLayer`.

### 7.5 `turbParams`

This section is required whenever `-t` is used or Plume is enabled in the integrated executable.

| Tag | Required | Meaning/default |
|---|---:|---|
| `method` | Yes | Mixing length: `0` height above terrain, `1` serial ray tracing, `2` ray tracing, `3` OptiX, `4` read from file. |
| `samples` | No | Samples per air cell for ray methods. |
| `save` | No | Save mixing length to `LMfile`; default `false`. |
| `LMfile` | Conditional | Mixing-length NetCDF file. Required for method `4` or saving. |
| `varname` | No | NetCDF variable name; default `mixlength`. |
| `nonLocalMixing` | No | Enable non-local mixing; default `false`. |
| `terrainWallFlag` | No | Terrain wall treatment; default `2`. |
| `buildingWallFlag` | No | Building wall treatment; default `2`. |
| `turbUpperBound` | No | Upper bound for turbulence fields; default `20`. |
| `backgroundMixing` | No | Background mixing value; default `0`. |
| `sigmaConst` | No | Three constant sigma values in one element. |

If method `3` is requested in a build without OptiX, QES warns and falls back to method `1`. Method `0` is the portable starting point.

### 7.6 `fileOptions`

When an output basename is supplied, keep this section in the Winds XML. The current visualisation writer dereferences `fileOptions` directly.

```xml
<fileOptions>
  <outputFlag>1</outputFlag>
  <outputFields>all</outputFields>
</fileOptions>
```

`outputFlag` is required when the section exists. Valid visualisation fields defined by the current writer are `u`, `v`, `w`, `mag`, `icell`, `icellInitial`, and `terrain`; `all` selects all of them. `massConservedFlag`, `sensorVelocityFlag`, and the intentionally misspelled `staggerdVelocityFlag` are also parsed but are not the main output-selection controls in the current writer.

## 8. QES-Plume XML

The root must be exactly:

```xml
<QESPlumeParameters>
  ...
</QESPlumeParameters>
```

The current v2.3.1 root parser requires `plumeParameters`, `collectionParameters`, and `boundaryConditions`. `particleParameters` is technically optional, but a run without it has no particle model and produces no concentration file.

### 8.1 `plumeParameters`

| Tag | Required | Meaning/constraint |
|---|---:|---|
| `simDur` | Yes | Duration to advance when the Winds input has one time record, seconds; must be greater than `0`. |
| `timeStep` | Yes | Overall Plume integration step, seconds; `0 < timeStep <= simDur`. |
| `CourantNumber` | Yes | Particle sub-stepping control in `[0, 1]`. `0` disables Courant subdivision. |
| `invarianceTol` | Yes | Stress-tensor realisability tolerance; must be greater than `0`. |
| `updateFrequency_particleLoop` | Yes | Particle-loop reporting interval; integer at least `1`. |
| `updateFrequency_timeLoop` | Yes | Simulation-time reporting interval; must be non-negative. |
| `interpolationMethod` | No | `triLinear` (default), `nearestCell`, or `analyticalPowerLaw`. Case-sensitive. |

`C_0` appears in older XML files but is not read by the current parser.

For multiple Winds time records, Plume loads each Winds/Turb record in turn and advances to the next Winds timestamp. For the last record, it reuses the interval between the final two Winds timestamps. In that multi-record path, `simDur` is still required by the parser but is not used to choose the final end time.

### 8.2 `collectionParameters`

This defines the Eulerian concentration sampling grid, which may be coarser than the Winds grid.

| Tag | Meaning/constraint |
|---|---|
| `timeAvgStart` | Seconds after simulation start before collection begins; at least `0`. |
| `timeAvgFreq` | Averaging and output period in seconds; greater than `0`. |
| `boxBoundsX1`, `boxBoundsX2` | Lower and upper x bounds in metres. |
| `boxBoundsY1`, `boxBoundsY2` | Lower and upper y bounds in metres. |
| `boxBoundsZ1`, `boxBoundsZ2` | Lower and upper z bounds in metres. |
| `nBoxesX`, `nBoxesY`, `nBoxesZ` | Number of collection boxes; each at least `1`. |

Lower bounds must not exceed upper bounds. Keep the collection grid within the usable Winds domain. Approximate memory grows with `nBoxesX * nBoxesY * nBoxesZ`, so start coarse and refine only where required.

### 8.3 `particleParameters`

The current schema groups sources by particle tag:

```xml
<particleParameters>
  <baseParticle>
    <tag>tracer</tag>
    <source>
      ...exactly one release type...
      ...exactly one source geometry...
    </source>
  </baseParticle>
</particleParameters>
```

Repeat `source` to combine sources into one output species/tag. Repeat `baseParticle` to generate separate concentration files. Tags should be unique and filename-safe because each output is named `<base>_<tag>_plumeOut.nc`.

The parser also accepts `depositionFlag`, `c1`, `c2`, and `decayConst` under `baseParticle`, but deposition and decay are not connected to the current `ParticleModel` advection path: deposition construction and execution are commented out. Do not interpret those settings as active physics in v2.3.1.

### 8.4 Release types

Each source must contain exactly one of the following case-sensitive elements.

Continuous release:

```xml
<releaseType_continuous>
  <particlePerTimestep>100</particlePerTimestep>
  <releaseStartTime>0.0</releaseStartTime>
  <massPerSec>1.0</massPerSec>
</releaseType_continuous>
```

- `particlePerTimestep` is required.
- `releaseStartTime` defaults to `0`.
- Set at most one of `massPerSec` and `massPerTimestep`.
- If neither mass is set, particles carry no configured mass and concentration remains zero even though `p_count` is populated.

Instantaneous release:

```xml
<releaseType_instantaneous>
  <releaseTime>0.0</releaseTime>
  <numPar>1000</numPar>
  <totalMass>10.0</totalMass>
</releaseType_instantaneous>
```

`numPar` is required and must be positive. `releaseTime` and `totalMass` default to `0`.

Finite-duration release:

```xml
<releaseType_duration>
  <releaseStartTime>0.0</releaseStartTime>
  <releaseDuration>60.0</releaseDuration>
  <particlePerTimestep>100</particlePerTimestep>
  <massPerSec>1.0</massPerSec>
</releaseType_duration>
```

- `releaseStartTime` and `particlePerTimestep` are required.
- Set exactly one of `releaseEndTime` and `releaseDuration`.
- Set at most one of `totalMass` and `massPerSec`.

Use a zero release start in v2.3.1 unless you have validated timing carefully. The controller construction adds the stored end value to the already shifted start time, which can extend non-zero-start instantaneous and finite-duration releases beyond the intuitive end time.

### 8.5 Source geometries

Each source must contain exactly one geometry:

| Element | Required children |
|---|---|
| `sourceGeometry_Point` | `posX`, `posY`, `posZ` |
| `sourceGeometry_Line` | `posX_0`, `posY_0`, `posZ_0`, `posX_1`, `posY_1`, `posZ_1` |
| `sourceGeometry_Cube` | `minX`, `minY`, `minZ`, `maxX`, `maxY`, `maxZ` |
| `sourceGeometry_SphereShell` | `posX`, `posY`, `posZ`, `radius` |
| `sourceGeometry_FullDomain` | No children; uses interpolation-domain bounds. |

Keep point, line, cube, and sphere coordinates inside the interpolation domain and out of solid terrain/building cells. Much of the older explicit geometry validation is commented out, so a malformed geometry may fail later or simply create inactive particles.

`particleDensity` and `particleDiameter` children belong directly under `source`. Although the XML parser marks them optional, set both explicitly: the current source builder always copies both local values into generated particles. They do not currently enable settling/deposition, because that physics is inactive in the present advection path. Use zero values for an ideal tracer if those properties are not otherwise relevant.

### 8.6 Boundary conditions

All three domain boundary types are required and case-sensitive:

```xml
<boundaryConditions>
  <xBCtype>exiting</xBCtype>
  <yBCtype>periodic</yBCtype>
  <zBCtype>reflection</zBCtype>
  <wallReflection>stairstepReflection</wallReflection>
</boundaryConditions>
```

For `xBCtype`, `yBCtype`, and `zBCtype`, choose:

- `exiting`: deactivate a particle that leaves the domain.
- `periodic`: wrap it to the opposite side.
- `reflection`: reflect position and velocity fluctuation at the domain boundary.

`wallReflection` is optional and defaults to `doNothing`. Supported values are:

- `doNothing`
- `setInactive`
- `stairstepReflection`
- `meshReflection` (falls back to stair-step if no triangle mesh is available)

Domain boundary conditions and wall reflection are different: the first handles outer domain edges; the second handles terrain/building solids.

### 8.7 Particle CSV output

To use `qesPlume -l` or `qes -l`, include:

```xml
<particleOutputParameters>
  <outputStartTime>0.0</outputStartTime>
  <outputFrequency>10.0</outputFrequency>
</particleOutputParameters>
```

`outputFrequency` is required. `outputStartTime` defaults to the first simulation time plus one frequency interval. `outputEndTime` and repeated `outputFields` are parsed but not used by the current CSV writer. Each snapshot contains all active particles and fixed columns, with a filename like:

```text
<base>_<tag>_particleData_2026-01-01T00:00:10.csv
```

Particle snapshots can become very large; enable them only for diagnostics.

## 9. Output files and variables

### Winds visualisation: `<base>_windsOut.nc`

The cell-centred output defines dimensions `t`, `z`, `y`, and `x`. Available fields are:

| Variable | Meaning | Units |
|---|---|---|
| `u`, `v`, `w` | Cell-centred velocity components | `m s-1` |
| `mag` | Velocity magnitude | `m s-1` |
| `terrain` | Terrain height | `m` |
| `icell` | Final cell flag | dimensionless |
| `icellInitial` | Initial cell flag | dimensionless |

### Winds workspace: `<base>_windsWk.nc`

This is a model interchange/restart-style file, not just a plotting product. It contains staggered `u`, `v`, `w`, initial velocities `u0`, `v0`, `w0`, cell flags, terrain and roughness, mixing length, vertical spacing, cut-cell solver coefficients, and building/terrain volume fractions.

### Turbulence: `<base>_turbOut.nc`

The turbulence file includes:

- velocity-gradient tensor components `Gxx` through `Gzz`;
- mixing length `L`;
- stress components `txx`, `tyy`, `tzz`, `txy`, `txz`, `tyz`;
- turbulent kinetic energy `tke`;
- dissipation-related field `CoEps`;
- stress-divergence components `div_tau_x`, `div_tau_y`, `div_tau_z`;
- `iturbflag`.

### Plume concentration: `<base>_<tag>_plumeOut.nc`

Each particle tag receives its own NetCDF file. The main variables are:

| Variable | Meaning | Units |
|---|---|---|
| `x_c`, `y_c`, `z_c` | Collection-box centres | `m` |
| `t_avg` | Actual averaging duration | `s` |
| `p_count` | Particle count accumulated per box | particle count |
| `c` | Time-averaged mass concentration | `g m-3` |

Inspect any NetCDF header with:

```bash
ncdump -h results/minimal_windsOut.nc
ncdump -h results/minimal_turbOut.nc
ncdump -h results/minimal_tracer_plumeOut.nc
```

## 10. Troubleshooting

### `cudaGetDeviceCount` fails before `--help`

The binary was compiled with CUDA and performs a startup GPU check. Run on a GPU host or build with `-DHAS_CUDA_SUPPORT=OFF -DHAS_OPTIX_SUPPORT=OFF`.

### `qesPlume` cannot open a Winds or Turb file

Confirm that the Winds stage used both `-w` and `-t`, and that the basenames match exactly:

```bash
ls -l results/minimal_windsWk.nc results/minimal_turbOut.nc
```

With `qesPlume -f results/minimal`, do not include `_windsWk.nc` in the argument; `-f` expects the shared base only.

### A command fails after `-w 1`, `-t 1`, or `-l 1`

These are switches and accept no value. Use `-w`, `-t`, or `-l` alone.

### `Turbulence model is turned on without turbParams`

Add `<turbParams><method>0</method></turbParams>` to the Winds XML, or remove the operation that requires Turb. Standalone Plume always requires a turbulence input file.

### Plume produces particle counts but zero concentration

Set a positive mass using one supported mass field in the release type. Particle count and concentration are separate: `p_count` can be non-zero while `c` is zero when no mass was configured.

### Plume runs but writes no concentration file

Use the current `<particleParameters><baseParticle>...` schema with a unique `<tag>` and at least one valid source. Older top-level `<sourceParameters>` blocks are not read by the v2.3.1 `PlumeInputData` root parser.

### `-l` crashes or fails during initialisation

Add `<particleOutputParameters>` with `outputFrequency` to the Plume XML. The command-line flag assumes that section exists.

### A relative DEM, WRF, shapefile, coefficient, sensor, or mixing-length path fails

Resolve it from the directory containing the Winds XML, not necessarily from the shell’s current directory. The path’s parent directory must already exist so QES can canonicalise it.

### An older XML example is ignored or behaves strangely

Check these schema changes and exact spellings:

| Legacy/stale form | Current v2.3.1 form |
|---|---|
| Plume `<simulationParameters>` | `<plumeParameters>` |
| Top-level Plume `<sourceParameters>` or `<sources>` | `<particleParameters><baseParticle><source>...` |
| `parPerTimestep` | `particlePerTimestep` |
| `ReleaseType_continuous` | `releaseType_continuous` |
| `SourcePoint` | `sourceGeometry_Point` |
| `UTMX`, `UTMY` | `UTMx`, `UTMy` |
| `<doDepositions>` | Parsed nowhere in the current root schema; deposition is currently inactive. |
| `C_0` | Not parsed by the current Plume parameter class. |

## 11. Pre-run checklist

Before a production run, verify:

1. The executable reports the expected QES version and CPU/CUDA/OptiX capabilities.
2. Winds and Plume XML roots and tag capitalisation match this guide.
3. `dx == dy` for the current domain implementation.
4. Sensor/source/collection coordinates lie inside the domain.
5. Every Winds time step has consistent sensor data.
6. `-t` has a matching `turbParams` section.
7. Standalone Plume has both `_windsWk.nc` and `_turbOut.nc` from the same Winds run.
8. Each Plume particle tag is unique and filename-safe.
9. Every Plume source has exactly one release type and one geometry.
10. A positive release mass is configured if physical concentration is required.
11. Collection-grid size is affordable in memory and storage.
12. A small case completes and its NetCDF headers, timestamps, dimensions, and ranges look plausible before scaling up.

## 12. Code locations used for this guide

The primary user-facing behaviour comes from:

- `qesWinds/handleWINDSArgs.cpp` and `qesWinds/qesWindsMain.cpp`
- `qesPlume/handlePlumeArgs.cpp` and `qesPlume/qesPlumeMain.cpp`
- `qes/handleQESArgs.cpp` and `qes/qesMain.cpp`
- `src/winds/WINDSInputData.h`, `SimulationParameters.h`, `Sensor.h`, `TimeSeries.h`, `BuildingsParams.h`, `VegetationParams.h`, `TURBParams.h`, and the Winds/Turb output classes
- `src/plume/PLUMEInputData.h`, all `PI_*` input classes, `PLUMEGeneralData.cpp`, `ParticleModel.cpp`, `Concentration.cpp`, and `ParticleOutput.cpp`
- `src/util/ParseInterface.h`, `QESFileSystemHandler.cpp`, and the NetCDF output classes

The Doxygen target (`windsdoc`) generates an API reference when Doxygen is installed. It complements this operational guide but does not replace it.
