# zmath

Zig 0.16.0 IEEE and scalar math for `f16`, `f32`, and `f64`, plus portable `@Vector` overloads for the elementwise operations. The repository uses `jj` with colocated Git.

## Layout

- `src/root.zig` is the root module. `src/ieee.zig` and `src/scalar.zig` implement the current IEEE and scalar APIs. `src/vector.zig` implements the portable `@Vector` layer for the elementwise operations.
- Source and tests are flat under `src/`: scalar kernels use `src/scalar*.zig`, the vector layer uses `src/vector*.zig`, C kernels use `src/coremath_*.[ch]`, and tests use `src/*_test.zig`.
- The vector layer uses Zig's portable `@Vector` and is lowered by the compiler per target; it has no hand-written per-ISA kernels and no runtime dispatch.
- There is no `src/scalar/`, `src/simd/`, or `src/tests/` directory.
- `build.zig` exposes the `zmath` module, installs a demonstration executable, and defines `run`, `test`, and `docs` steps.
- The sibling [`zmath-testing`](../zmath-testing/) project owns oracle generation, corpus replay, and benchmarks.

## Numerical contracts

- Handle IEEE-754 NaNs, signed zero, infinities, and subnormals according to [API.md](API.md) and the [scalar contract](spec/core-math-2026-08-22.md).
- Use the [schema-3 corpus manifest](../zmath-testing/test/corpus/manifest.toml) for acceptance status. An export or a sampled corpus does not prove whole-domain accuracy.
- `build.zig` requires `+fma` for x86/x86_64 C kernels and supports little-endian aarch64. The vector layer adds no dispatch; the compiler lowers `@Vector` for the selected target. There is no baseline SSE2 fallback and no hand-written SIMD kernel selection.
- Do not allocate dynamically in evaluation loops. Keep scratch storage statically bounded or caller-supplied.

## Commands

Run at the repository root:

- `zig build` installs the demonstration executable at `zig-out/bin/zmath` and makes the `zmath` module available to consumers. It does not install a compiled library artifact.
- `zig build run` runs the demonstration. `zig build test` runs the local unit tests. `zig build docs` generates documentation from the demonstration executable.
- `zig fmt --check .` checks formatting. `gate` runs the house gate.
- `jj st` shows working-copy status. Use `jj`, not raw Git commands.

Run oracle and benchmark commands in sibling `zmath-testing`, not here. Its accepted corpus profile is `x86_64-linux-gnu.2.38`, CPU `haswell`, and `ReleaseSafe`. See [API.md: Tests and benchmarks](API.md#tests-and-benchmarks) for its current build steps.

This NixOS project declares development dependencies and tools in `devenv.nix`. Do not install system or development libraries imperatively.
Use direnv or `devenv shell`; run `direnv allow` after changing the environment.

