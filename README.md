# zmath

IEC 60559 (`IEEE 754`) elementary functions for Zig, across `f16`, `f32`, `f64`
and `@Vector`s of those.

`zmath` implements the standard's recommended operations rather than whatever
the platform's libm happens to provide, so results are the same on every target
Zig supports. Correctness is pinned against the standard's own specification
(`spec/core-math-2026-08-22.md`) and checked against a differential oracle
(see `zmath-testing`).

## Modules

Nine namespace modules hold the full surface, and each carries its own contract:

| module | contents |
|---|---|
| `ieee` | encoding, classification, payload operations |
| `scalar` | the per-format elementary functions |
| `vector` | the portable `@Vector` layer and the `Vec*` value types |
| `complex` | complex arithmetic |
| `quaternion` | quaternion arithmetic |
| `matrix` | matrix types and operations |
| `geometry` | geometric primitives |
| `linalg` | decompositions and solvers |
| `algebra` | polynomial and related operations |

`root.zig` adds generic convenience functions over `scalar` and `ieee` and
re-exports their shared result types.

## Usage

The floating-point argument type selects the format:

```zig
const std = @import("std");
const zmath = @import("zmath");

pub fn main() void {
    // Both calls use binary32: the scalar and the four-lane vector.
    const scalar_result = zmath.sin(@as(f32, 0.5));
    const lane_result = zmath.sin(@as(@Vector(4, f32), @splat(0.5)));
    std.debug.print("{d} {d}\n", .{ scalar_result, lane_result[0] });
}
```

Generic operations accept supported scalar floats or `@Vector`s of those, except
`nan`, which takes an explicit scalar type and a payload string. A vector
argument keeps its lane count.

Operations not named in `root.zig` are reached through their namespace:

```zig
const zmath = @import("zmath");

const bits = zmath.ieee.toBits(@as(f32, 1.0));   // encoding
const c = zmath.complex.exp(zmath.complex.Complex.f32.init(0, 0));
```

## Building

```sh
zig build            # compile the library
zig build test       # run the unit suite
zig build test-build # build the test binaries without running them
zig build run        # run the demonstration
zig build docs       # generate the autodoc site into zig-out/docs
```

Requires Zig 0.16.0. The repository has no third-party dependencies; the
`coremath_*` translation units are vendored reference sources from the
CORE-MATH project, not a dependency.

## Correctness

The suite is organised around what can actually go wrong in a floating-point
library, and the files name the concern: `ieee_edge_test.zig`,
`ieee_exhaustive_h_test.zig` (every binary16 value), `ieee_scale_widths_test.zig`
(the widening and narrowing boundaries), `root_pow_test.zig`,
`root_sinpi_test.zig`, `root_tgamma_contract_test.zig`,
`c23_public_surface_test.zig`, and `api_suffix_test.zig` for the naming
convention that selects a format by suffix.

Differential testing against an independent oracle, corpus replay, fuzzing and
benchmark harnesses live in `zmath-testing`.

## License

See [LICENSE](LICENSE).