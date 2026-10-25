# zmath

Experimental IEC 60559 (IEEE 754) elementary functions for Zig `f16`, `f32`,
`f64`, and vectors of those formats. The library also provides complex numbers,
quaternions, matrices, linear algebra, geometry, and geometric algebras.

## Requirements

- Zig 0.17.0.
- FMA support on x86 and x86-64 (`-Dcpu=baseline+fma` when cross-compiling).
- A C target with IEEE binary32/binary64, known byte order, and `_BitInt(128)`
  or `__int128` support.

No external runtime libraries are required. Private C kernels derived from
[CORE-MATH](https://core-math.gitlabpages.inria.fr/) are included in `src/`.
Public fixed-size operations do not allocate.

## Usage

Fetch a revision into your application's dependencies:

```sh
zig fetch --save=zmath git+https://github.com/CoffeeImpliesCode/zmath
```

Import its `zmath` module from your `build.zig`:

```zig
const dependency = b.dependency("zmath", .{
    .target = target,
    .optimize = optimize,
});
exe.root_module.addImport("zmath", dependency.module("zmath"));
```

The module includes its private C sources and link settings. In application code:

```zig
const std = @import("std");
const zmath = @import("zmath");

pub fn main() void {
    const scalar_result = zmath.sin(@as(f32, 0.5));
    const lane_result = zmath.sin(@as(@Vector(4, f32), @splat(0.5)));
    const bits = zmath.ieee.toBitsh(@as(f16, 1.0));
    const c = zmath.complex.expf(.{ .re = 0, .im = 0 });
    std.debug.print("{d} {d} {x} {d}\n", .{
        scalar_result, lane_result[0], bits, c.re,
    });
}
```

The argument type selects the format for generic scalar and vector operations.
Suffixed names use `h`, `f`, and `d` for `f16`, `f32`, and `f64`.
See [API.md](API.md) for result types, rounding rules, and error contracts.

## Build and test

```sh
zig build            # install the demonstration executable
zig build test       # run the numerical unit suite
zig build test-build # install the test executable without running it
zig build run        # run the demonstration
zig build docs       # install the autodoc site in zig-out/docs
```

Serve `zig-out/docs` over HTTP to view the generated documentation.

## Layout

| Module | Contents |
|---|---|
| `ieee` | Encoding, classification, payloads, and rounding |
| `scalar` | Elementary functions for each format |
| `vector` | Elementwise `@Vector` operations and `Vec*` value types |
| `complex` | Complex arithmetic |
| `quaternion` | Raw and unit quaternions |
| `matrix` | Fixed-size matrices and operations |
| `linalg` | Decompositions and solvers |
| `geometry` | Geometric primitives and predicates |
| `algebra` | Clifford, projective, and conformal geometric algebras |

Implementation and unit tests are under `src/`. The root module adds generic
convenience functions over `scalar` and `ieee`.

Independent oracle generation, corpus replay, fuzzing, and benchmarks live in
[zmath-testing](https://github.com/CoffeeImpliesCode/zmath-testing).
The [scalar corpus contract](spec/core-math-2026-08-22.md) defines that harness's
schema. Unit tests and selected corpus cases do not prove whole-domain accuracy;
the non-scalar modules have no general independent accuracy proof.

Big-endian vector code retains compiler-specific workarounds. Tests exclude
some boolean mask widths, and returned-vector register ordering has known
compiler limits on ARM big-endian targets. Cross-target support requires tests
on the selected compiler and target.

## License

First-party code is licensed under MIT; see [LICENSE](LICENSE).
Derived CORE-MATH files retain their upstream copyright and MIT license notices.
