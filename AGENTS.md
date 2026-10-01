# zmath

Zig 0.16.0 IEEE and scalar math for `f16`, `f32`, and `f64`, plus a portable `@Vector` layer over the elementwise scalar operations, including the result-struct, integer-argument, integer-returning, boolean-returning, and narrowing forms. The repository uses `jj` with colocated Git.

## Layout

- `src/root.zig` is the root module and exports `zmath.ieee`, `zmath.scalar`, `zmath.vector`, `zmath.complex`, `zmath.quaternion`, `zmath.matrix`, `zmath.linalg`, `zmath.geometry`, and `zmath.algebra`. `src/ieee.zig` and `src/scalar.zig` implement the current IEEE and scalar APIs. `src/vector.zig` holds the float-to-float `@Vector` surface and re-exports `src/vector_core.zig` (shared per-lane helpers), `src/vector_struct.zig`, `src/vector_intarg.zig`, `src/vector_narrow.zig`, and `src/vector_ieee.zig` for the other lane-to-non-lane forms. It also re-exports the six `Vec*` value types from `src/vector_value.zig`, which is a different layer that shares the module path.
- Source and tests are flat under `src/`: scalar kernels use `src/scalar*.zig`, the vector layer uses `src/vector*.zig`, C kernels use `src/coremath_*.[ch]`, and tests use `src/*_test.zig`.
- The non-scalar modules are `src/complex.zig`, `src/quaternion.zig`, `src/matrix.zig`, `src/linalg.zig`, `src/geometry.zig`, and `src/algebra.zig`. `src/matrix.zig` holds the one sanctioned public structural generic, `Mat`; every other public type is concretely named. `src/geometry.zig` carries the `PredicatePolicy` enum and its `adaptive_exact` path evaluates predicates in exact integer arithmetic on a fixed `[132]u64` limb array.
- The vector layer uses Zig's portable `@Vector` and is lowered by the compiler per target; it has no hand-written per-ISA kernels and no runtime dispatch.
- There is no `src/scalar/`, `src/simd/`, or `src/tests/` directory.
- `build.zig` exposes the `zmath` module, installs a demonstration executable, and defines `run`, `test`, and `docs` steps.
- The sibling [`zmath-testing`](../zmath-testing/) project owns oracle generation, corpus replay, and benchmarks.

## Numerical contracts

- Handle IEEE-754 NaNs, signed zero, infinities, and subnormals according to [API.md](API.md) and the [scalar contract](spec/core-math-2026-08-22.md).
- Use the [schema-5 corpus manifest](../zmath-testing/test/corpus/manifest.toml) for acceptance status. An export or a sampled corpus does not prove whole-domain accuracy.
- Every `@Vector` lane runs the identical accepted scalar kernel, so a vector result is bit-for-bit the scalar result for that lane; never re-derive a polynomial for the vector path. Zig forbids a struct as a `@Vector` element type, so the thirteen result-struct operations and `pown`, `compoundn`, and `rootn` return a lane-ordered `[N]Result` array instead of a `@Vector`. That is a language limitation, not a placeholder, and the array shape buys no throughput over a scalar loop.
- `build.zig` requires `+fma` for x86/x86_64 C kernels and supports little-endian aarch64. The vector layer adds no dispatch; the compiler lowers `@Vector` for the selected target. There is no baseline SSE2 fallback and no hand-written SIMD kernel selection.
- Do not allocate dynamically in evaluation loops. Keep scratch storage statically bounded or caller-supplied.

## Commands

Run at the repository root:

- `zig build` installs the demonstration executable at `zig-out/bin/zmath` and makes the `zmath` module available to consumers. It does not install a compiled library artifact.
- `zig build run` runs the demonstration. `zig build test` runs the local unit tests. `zig build docs` generates documentation from the demonstration executable.
- `zig fmt --check .` checks formatting. `gate` runs the house gate: format, build, test, AST check, and three linters. `gate` is a devenv script; `.envrc` synthesizes a `gate` shim into the gitignored `.devenv/host-bin` because a devenv script is a store derivation and this host has no `/nix/store`. Run it as `direnv exec . gate`, or `direnv exec . zig-gate` if the shim is missing.
- `jj st` shows working-copy status. Use `jj`, not raw Git commands.

Run oracle and benchmark commands in sibling `zmath-testing`, not here. Its accepted corpus profile is `x86_64-linux-gnu.2.38`, CPU `haswell`, and `ReleaseSafe`. See [API.md: Tests and benchmarks](API.md#tests-and-benchmarks) for its current build steps.

This NixOS project declares development dependencies and tools in `devenv.nix`. Do not install system or development libraries imperatively.
Use direnv or `devenv shell`; run `direnv allow` after changing the environment. Nix is not
installed natively on this host: `nix` and `devenv` run inside an apptainer sandbox, so
`/nix/store` does not exist and every store path a devshell exports is unreachable. That is why
`.envrc` puts `~/.local/bin` on `PATH` before `use devenv`, and why a `devenv direnvrc` run prints
one harmless `readlink: Datei oder Verzeichnis nicht gefunden` from the generated shell hook.

