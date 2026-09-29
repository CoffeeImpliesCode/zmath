# zmath

`zmath` is a Zig 0.16.0 mathematics module for `f16`, `f32`, and `f64`.

## Implemented

- The public root exposes IEEE-754 and scalar math through `zmath.ieee`, `zmath.scalar`, and scalar root operations, plus a portable `@Vector` layer through `zmath.vector`: 417 public entry points over the elementwise scalar operations, covering the float-to-float, result-struct, integer-argument, integer-returning, boolean-returning, and narrowing forms.
- Implementation and tests are flat under `src/`: `src/root.zig`, `src/ieee.zig`, `src/scalar*.zig`, `src/vector*.zig`, `src/*_test.zig`, and C kernels in `src/coremath_*.[ch]`.
- Each vector lane applies the identical accepted scalar kernel, so a vector result is bit-for-bit equal to the accepted scalar result for that lane and the `max_ulp = 0` contract is preserved. Zig forbids a struct as a `@Vector` element type, so the result-struct operations return a lane-ordered `[N]Result` array rather than a `@Vector`; that is a language limitation, and the array shape buys no throughput over a scalar loop.
- The runtime has no external oracle or math-library dependency. The vector layer uses portable Zig `@Vector` and is lowered by the compiler per target; it has no hand-written per-ISA kernels, no runtime dispatch, and no SSE2 fallback.
- `build.zig` makes the `zmath` module available to consumers and installs the `zig-out/bin/zmath` demonstration executable, not a compiled library artifact.
- C kernels require `+fma` on x86/x86_64 or a little-endian aarch64 target.

## Planned

The [API contract](API.md) specifies complex numbers, quaternions, matrices, linear algebra, geometry, and geometric algebra. None of those modules is implemented. The contract also places the linear-algebra `Vec2f` through `Vec4d` value types under the same `zmath.vector` module path; those are target names and are not implemented either. The implemented `zmath.vector` module is the portable `@Vector` layer over the elementwise scalar operations, not the linear-algebra `Vec*` types. The API distinguishes target names from current source exports.

## Numerical evidence

The [scalar contract](spec/core-math-2026-08-22.md) defines the corpus schema and target profile. The [checked-in manifest](../zmath-testing/test/corpus/manifest.toml) in sibling `zmath-testing` is the acceptance authority. Its oracle generation, corpus replay, and benchmarks are separate from this runtime. A source export or selected corpus does not prove whole-domain accuracy.

## Commands

From this repository, use `zig build` to install the demonstration and `zig build run` to run it. Use `zig build test` for local tests. This build graph has no `bench` step. Run qualified oracle and benchmark steps from `zmath-testing` as documented in [API.md](API.md#tests-and-benchmarks).
