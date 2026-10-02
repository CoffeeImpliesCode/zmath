# zmath

`zmath` is a Zig 0.16.0 mathematics module for `f16`, `f32`, and `f64`.

## Implemented

All nine module paths named by the [API contract](API.md) exist and are exported from the root.

| Module | Source | Surface |
| --- | --- | --- |
| `zmath.ieee` | `src/ieee.zig` | The C23 IEEE primitives |
| `zmath.scalar` | `src/scalar.zig` | The accepted scalar kernels |
| `zmath.vector` | `src/vector.zig` plus six siblings and `src/vector_value.zig` | 417 `@Vector` entry points and the six `Vec*` value types |
| `zmath.complex` | `src/complex.zig` | `Complexf`, `Complexd`, and 54 suffixed functions |
| `zmath.quaternion` | `src/quaternion.zig` | `Quatf`, `Quatd`, `UnitQuatf`, `UnitQuatd` |
| `zmath.matrix` | `src/matrix.zig` | The `Mat` generic, eight size aliases, 19 module functions |
| `zmath.linalg` | `src/linalg.zig` | Six LU and six Cholesky factorization types, 10 module functions |
| `zmath.geometry` | `src/geometry.zig` | 2D/3D/4D point, vector, line, ray, segment, flat, ball, AABB families for `f` and `d`, 59 types |
| `zmath.algebra` | `src/algebra.zig` | `Cl3`, `Cl4`, `Pga3`, `Pga4`, `Cga3`, `Cga4` in `f` and `d`, with a published normalization table |

- Implementation and tests are flat under `src/`: `src/root.zig`, `src/ieee.zig`, `src/scalar*.zig`, `src/vector*.zig`, the six non-scalar modules, `src/*_test.zig`, and C kernels in `src/coremath_*.[ch]`.
- The root exports 97 top-level `pub fn` operations and 25 public types and module aliases.
- Each `@Vector` lane applies the identical accepted scalar kernel, so a vector result is bit-for-bit equal to the accepted scalar result for that lane and the `max_ulp = 0` contract is preserved. Zig forbids a struct as a `@Vector` element type, so the thirteen result-struct operations and the `pown`, `compoundn`, and `rootn` power family return a lane-ordered `[N]Result` array rather than a `@Vector`; that is a language limitation, and the array shape buys no throughput over a scalar loop.
- `zmath.vector` carries two unrelated layers that share the module path: the portable `@Vector` elementwise surface, and the named-component `Vec2f` through `Vec4d` value types from `src/vector_value.zig`. The elementwise operations use `v` or `vLanes` suffixes; no `Vec*` type name does.
- `Mat` is the one public structural generic. Every other public type in the library is concretely named, with private comptime generics behind it so the `f` and `d` forms cannot drift.
- `zmath.geometry` takes an explicit `PredicatePolicy`. The `adaptive_exact` path evaluates the predicate in exact integer arithmetic on a fixed `[132]u64` limb array, with no epsilon and no fallback to the `fast` answer.
- The runtime has no external oracle or math-library dependency. The vector layer uses portable Zig `@Vector` and is lowered by the compiler per target; it has no hand-written per-ISA kernels, no runtime dispatch, and no SSE2 fallback.
- `build.zig` makes the `zmath` module available to consumers and installs the `zig-out/bin/zmath` demonstration executable, not a compiled library artifact.
- C kernels require IEEE binary32/binary64 storage, a known little- or big-endian byte order, and `_BitInt(128)` or native `__int128`. The binary16 C ABI is bit-encoded `uint16_t`, with optional native conversion and direct portable RNE packing. C and Zig share a compile-time native-FMA policy and use bounded integer binary32/binary64 fusion where hardware fusion is absent. The x86/x86_64 `+fma` floor remains.

## Not covered by evidence

The non-scalar modules have no general accuracy proof against an independent oracle.
The corpus has ten `definitional` rows for four Cl(3) involutions and quaternion multiplication in `f32` and `f64`.
These rows check 176 stored convention cases and do not establish either complete module's numerical accuracy.
Ten further numerical rows are five selected operations, each with one `f32` and one `f64` row: `complex` addition and multiplication, `matrix.determinantMat2`, `linalg.lu` with `solveLu` on a 2x2 system, and the 2D orientation predicate under both policies.
Those ten rows check 324 stored bounded cases against an independent exact-rational reference.
Each row follows the rounding of the graph it checks, so the rule is per primitive step for `complexMul`, `matrixDet2`, and the `fast` predicate, and one round per component for the exact rational `luSolve2` solution.
This is the selected sample the user accepted: the other operations of complex, matrix, linalg, and geometry stay outside the manifest as a scope limit.
Their unit tests do not substitute for oracle evidence.

## Numerical evidence

The [scalar contract](spec/core-math-2026-08-22.md) defines the corpus schema and target profile. The [checked-in manifest](../zmath-testing/test/corpus/manifest.toml) in sibling `zmath-testing` is the acceptance authority. Its oracle generation, corpus replay, and benchmarks are separate from this runtime. A source export or selected corpus does not prove whole-domain accuracy.

## Commands

From this repository, use `zig build` to install the demonstration and `zig build run` to run it. Use `zig build test` for local tests. This build graph has no `bench` step. Run qualified oracle and benchmark steps from `zmath-testing` as documented in [API.md: Tests and benchmarks](API.md#tests-and-benchmarks).

Use `cross-test` to attempt the complete 53-platform Zig 0.16.0 Tier 1/2
matrix inherited by `cross-test.profile.sh`; append `-- -Doptimize=ReleaseSafe`
for optimized correctness checks. `cross-test --list` shows detected runners.
Missing foreign-OS execution environments remain build-only, and
`--require-runners` makes incomplete execution coverage fail.
`zig build test-build` installs the same test executable without executing it.
These checks run no corpus or profiling workload and do not replace the gate.
