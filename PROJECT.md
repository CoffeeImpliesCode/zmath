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

- Implementation and tests are flat under `src/`: `src/root.zig`, `src/ieee.zig`, `src/scalar*.zig`, `src/vector*.zig`, the six non-scalar modules, `src/*_test.zig`, and C kernels in `src/coremath_*.[ch]`. 208 `.zig` files.
- The root exports 97 top-level `pub fn` operations and 25 public types and module aliases.
- Each `@Vector` lane applies the identical accepted scalar kernel, so a vector result is bit-for-bit equal to the accepted scalar result for that lane and the `max_ulp = 0` contract is preserved. Zig forbids a struct as a `@Vector` element type, so the thirteen result-struct operations and the `pown`, `compoundn`, and `rootn` power family return a lane-ordered `[N]Result` array rather than a `@Vector`; that is a language limitation, and the array shape buys no throughput over a scalar loop.
- `zmath.vector` carries two unrelated layers that share the module path: the portable `@Vector` elementwise surface, and the named-component `Vec2f` through `Vec4d` value types from `src/vector_value.zig`. Every `@Vector` entry point ends in `v` and no `Vec*` name does, so the call sites are unambiguous.
- `Mat` is the one public structural generic. Every other public type in the library is concretely named, with private comptime generics behind it so the `f` and `d` forms cannot drift.
- `zmath.geometry` takes an explicit `PredicatePolicy`. The `adaptive_exact` path evaluates the predicate in exact integer arithmetic on a fixed `[132]u64` limb array, with no epsilon and no fallback to the `fast` answer.
- The runtime has no external oracle or math-library dependency. The vector layer uses portable Zig `@Vector` and is lowered by the compiler per target; it has no hand-written per-ISA kernels, no runtime dispatch, and no SSE2 fallback.
- `build.zig` makes the `zmath` module available to consumers and installs the `zig-out/bin/zmath` demonstration executable, not a compiled library artifact.
- C kernels require `+fma` on x86/x86_64 or a little-endian aarch64 target. A big-endian target is a hard error, both at configure time and through a `_Static_assert` in the generated ABI guard.

## Not covered by evidence

The non-scalar modules are numerically unproven against any independent oracle. `zmath.matrix`, `zmath.linalg`, `zmath.geometry`, and `zmath.algebra` are exact-arithmetic-shaped operations whose accuracy rests on their own unit tests, not on the corpus. Only the IEEE and scalar surfaces inherit the `max_ulp = 0` contract from the corpus.

## Numerical evidence

The [scalar contract](spec/core-math-2026-08-22.md) defines the corpus schema and target profile. The [checked-in manifest](../zmath-testing/test/corpus/manifest.toml) in sibling `zmath-testing` is the acceptance authority. Its oracle generation, corpus replay, and benchmarks are separate from this runtime. A source export or selected corpus does not prove whole-domain accuracy.

## Commands

From this repository, use `zig build` to install the demonstration and `zig build run` to run it. Use `zig build test` for local tests. This build graph has no `bench` step. Run qualified oracle and benchmark steps from `zmath-testing` as documented in [API.md: Tests and benchmarks](API.md#tests-and-benchmarks).
