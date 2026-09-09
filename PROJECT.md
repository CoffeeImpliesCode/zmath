# zmath

`zmath` is a Zig 0.16.0 mathematics module for `f16`, `f32`, and `f64`.

## Implemented

- The public root exposes IEEE-754 and scalar math through `zmath.ieee`, `zmath.scalar`, and scalar root operations.
- Implementation and tests are flat under `src/`: `src/root.zig`, `src/ieee.zig`, `src/scalar*.zig`, `src/*_test.zig`, and C kernels in `src/coremath_*.[ch]`.
- The runtime has no external oracle or math-library dependency. There is no SIMD implementation or SSE2 fallback.
- `build.zig` makes the `zmath` module available to consumers and installs the `zig-out/bin/zmath` demonstration executable, not a compiled library artifact.
- C kernels require `+fma` on x86/x86_64 or a little-endian aarch64 target.

## Planned

The [API contract](API.md) specifies complex numbers, vectors, quaternions, matrices, linear algebra, geometry, and geometric algebra. These modules and SIMD kernels are not implemented. The API distinguishes target names from current source exports.

## Numerical evidence

The [scalar contract](spec/core-math-2026-08-22.md) defines the corpus schema and target profile. The [checked-in manifest](../zmath-testing/test/corpus/manifest.toml) in sibling `zmath-testing` is the acceptance authority. Its oracle generation, corpus replay, and benchmarks are separate from this runtime. A source export or selected corpus does not prove whole-domain accuracy.

## Commands

From this repository, use `zig build` to install the demonstration and `zig build run` to run it. Use `zig build test` for local tests. This build graph has no `bench` step. Run qualified oracle and benchmark steps from `zmath-testing` as documented in [API.md](API.md#tests-and-benchmarks).
