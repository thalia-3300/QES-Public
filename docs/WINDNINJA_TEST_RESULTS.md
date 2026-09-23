# WindNinja input v1 — convergence regression record

Tested locally on 2026-09-23. These are numerical/software regression tests,
not physical validation of WindNinja-to-QES downscaling.

## Outcome

The previously observed convergence failure and large CPU/GPU disagreement are
resolved in the tested configurations, without relaxing the 1e-7 tolerance.
All four WindNinja cases converge below the examples' 1,000-iteration cap.
The GPU suites used a 5,000-iteration safety cap, but stopped at the counts below;
the separate OpenMP CPU suite passed with a 1,000-iteration cap.

| Case | GPU iterations | Equation residual | Maximum final divergence (s^-1) | Maximum CPU/GPU velocity difference (m/s) |
| --- | ---: | ---: | ---: | ---: |
| Testing flat | 536 | 1.59e-9 | 4.84e-8 | 9.54e-7 |
| Testing building | 532 | 1.57e-9 | 4.95e-8 | 9.54e-7 |
| Clare flat | 544 | 1.61e-9 | 2.42e-7 | 3.81e-6 |
| Clare terrain | 464 | 1.29e-9 | 2.38e-7 | 3.81e-6 |

Global-memory (`-s 3`) and shared-memory (`-s 4`) GPU final fields match exactly
in all four cases. The legacy uniform sensor case converges in one iteration
with zero divergence. Complete initial CPU/GPU `u0/v0/w0` arrays match exactly.

Divergence is independently recomputed from saved staggered face velocities,
not inferred from the iteration stopping metric. It covers solved regular and
stair-step air cells, excluding solids, prescribed outer-ring/top cells, and
cut cells. Float32 final velocities impose a rounding floor, distinct from
double-precision potential convergence. Solid faces are separately checked for
zero velocity. This is not a full-domain or cut-cell mass-balance validation.

For comparison, the previous GPU outputs had maximum divergence of 0.515 s^-1
(Testing flat) and 1.325 s^-1 (Clare terrain) on the same solved-cell subset.
The previous maximum CPU/GPU u difference reached 0.747 m/s in Clare terrain.
Those old outputs should not be used as converged results: rerun with the rebuilt
binary to produce corrected wind and downstream turbulence/plume inputs.

## Defects corrected

- GPU RHS construction read beyond the vertical-spacing array at the top.
  NVIDIA memcheck reproduced the invalid read before the fix.
- The GPU ground mirror used a face-sized plane instead of a cell-sized plane,
  overwriting part of the first physical potential layer.
- The shared-memory error reduction relied on unsafe cross-block synchronization.
- Single-precision potentials and coefficient-sum accumulation could stall
  strict convergence. Potential, RHS, updates and reductions now use doubles;
  coefficient and output arrays remain floats.
- Final recovery skipped a top correction face. Its coefficient must come from
  the scaled row below, not the unscaled prescribed top row.
- The initial vertical divergence contained an extra dz factor. A manufactured
  case with nonzero initial w now guards this path.
- OpenMP reset its maximum error from multiple threads, and solid-cell-owned
  velocity masking could overlap neighboring writes. These now have single-writer
  operations.
- `Domain::getDomainSize()` returned dx in place of dy; unequal-spacing tests
  now exercise the corrected accessor.

All supported solvers use common RHS construction, velocity recovery, solid-face
masking and convergence diagnostics. SOR remains GPU-accelerated in modes 3/4.
Finite nonconvergence is explicitly logged, not a process-exit failure by itself;
the real-case runner's `--require-convergence` makes it a test failure.

## Verification

- CUDA 12.4 / GCC 13: rebuilt `qesWinds`, combined `qes`, `qesPlume`, and tests.
- Separate CPU-only GCC 15 / OpenMP build: both CTest regressions and real-case
  suite pass with four threads.
- Importer regression: 416 checks for interpolation, angles/units, profiles,
  staggered positions, terrain-relative heights, padding, CRS/DEM alignment,
  coverage and invalid source data.
- Manufactured-solution regression: 8,461 checks across serial CPU, red/black
  CPU and both GPU solvers; known corrected velocity, unequal x/y spacing,
  nonzero initial w, top-face correction, and deliberate nonconvergence at a
  one-iteration cap.
- NVIDIA memcheck and synccheck: zero errors for the manufactured GPU suite.
- Clare terrain example run directly with a 1,000-iteration cap under NVIDIA
  memcheck: converged at iteration 464, zero memory errors.
- Both real source grids accepted: Testing 142×142, Clare 51×51, including
  10,201 sampled interpolation positions per grid.
- CPU and both GPU modes: Testing flat/building, Clare flat/DEM terrain and
  legacy sensors; convergence, mass conservation, one output time and finite
  wind/turbulence outputs checked on an NVIDIA GeForce RTX 3060 Laptop GPU.
- Twelve input configurations rejected: incomplete coverage, multiple times,
  unknown speed units, missing attribute, nonzero halo, rotated domain, missing
  height, unknown option, repeated height/source, mixed sensors/source and
  DEM-origin mismatch.

## Remaining scientific limitations

The private shapefiles do not encode the required output height, roughness or
QES Monin–Obukhov length. Tests use m/s per the supplied clarification and
**unverified assumptions of 10 m AGL, z0=0.1 m, reciprocal=0**. This includes the
Clare filename containing `non_neutral_stability`. Confirm actual WindNinja
settings, including any above-vegetation height convention, before production.

Convergence and agreement between implementations do not establish agreement
with real wind. Terrain influence applied in both models, profile assumptions,
domain size/resolution and building parameterizations still need scientific
assessment. These tests do not validate irregular/multi-time input, cut-cell
mass balance, stretched vertical grids, or large-domain performance. Internal
double precision increases solver working memory and may affect GPU throughput.

## Reproduction and local artifacts

Follow [the input guide](WINDNINJA_INPUT.md). Run `run_real_cases.py` with
`--require-convergence`; default iteration cap is 1,000. Use `--solver 1`,
`--solver 3` and `--solver 4` to repeat each backend. Use the test executable's
`--compare`, `--compare-initial`, `--check-mass` and `--solver-tests-gpu` modes
for additional comparisons. NVIDIA checks can be repeated with:

```bash
compute-sanitizer --tool memcheck --error-exitcode 99 \
  ./build/tests/windShapefile/windShapefileTests --solver-tests-gpu
compute-sanitizer --tool synccheck --error-exitcode 99 \
  ./build/tests/windShapefile/windShapefileTests --solver-tests-gpu
```

Temporary outputs (not committed):

- Serial CPU: `/tmp/qes-windninja-s1-3g0ul9og`
- Global GPU: `/tmp/qes-windninja-s3-tnxjx_wo`
- Shared GPU: `/tmp/qes-windninja-s4-ylm3bnu7`
- OpenMP CPU-only: `/tmp/qes-windninja-s1-wrudmah8`
- CPU-only build: `/tmp/qes-windninja-cpu-build-OtYgsy`
- Sanitizer logs: `/tmp/qes-after-memcheck.log`, `/tmp/qes-after-synccheck.log`
  and `/tmp/qes-clare-after-memcheck.log`

Run directories contain XML, logs, workspace/visualization NetCDF and turbulence
outputs. They may be cleared by the system; regenerate them as needed.
