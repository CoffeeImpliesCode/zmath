# zmath Public API

## Purpose

`zmath` provides dependency-free mathematical primitives for Zig 0.16.

The API favors compile-time evaluation, fixed storage, and named semantic types. It uses one structural generic, `Mat`, because matrix dimensions and scalar type define a matrix. It does not expose a generic algebra generator.

No public operation requires an allocator.

## Public root

```zig
const zmath = @import("zmath");

pub const ieee = @import("ieee.zig");
pub const scalar = @import("scalar.zig");
pub const complex = @import("complex.zig");
pub const vector = @import("vector.zig");
pub const quaternion = @import("quaternion.zig");
pub const matrix = @import("matrix.zig");
pub const linalg = @import("linalg.zig");
pub const geometry = @import("geometry.zig");
pub const algebra = @import("algebra.zig");
```

The root module exports stable public modules and selected named types. It does not re-export implementation helpers.

## Module tree

```text
src/
  root.zig

  ieee.zig
  ieee/
    binary16.zig
    bfloat16.zig
    binary32.zig
    binary64.zig
    binary80.zig
    binary128.zig
    classify.zig
    order.zig
    rounding.zig
    scale.zig

  scalar.zig
  scalar/
    f16.zig
    bf16.zig
    f32.zig
    f64.zig
    f80.zig
    f128.zig
    constants.zig
    trig.zig
    exp_log.zig
    power.zig
    hyperbolic.zig
    special.zig

  complex.zig
  complex/
    complex32.zig
    complex64.zig
    complex80.zig
    complex128.zig
    arithmetic.zig
    elementary.zig
    branch.zig

  vector.zig
  vector/
    vec2.zig
    vec3.zig
    vec4.zig
    simd.zig

  quaternion.zig
  quaternion/
    quatf.zig
    quatd.zig

  matrix.zig
  matrix/
    mat.zig
    affine3.zig
    affine4.zig
    simd.zig

  linalg.zig
  linalg/
    solve.zig
    lu.zig
    cholesky.zig
    qr.zig
    eigen.zig

  geometry.zig
  geometry/
    primitives2.zig
    primitives3.zig
    primitives4.zig
    intersections2.zig
    intersections3.zig
    intersections4.zig
    predicates.zig

  algebra.zig
  algebra/
    cl3.zig
    cl4.zig
    pga3.zig
    pga4.zig
    cga3.zig
    cga4.zig
```

The tree describes the target public shape. A file is added only when its concrete API is implemented.

## Dependency order

```text
ieee -> scalar -> complex
             +-> vector -> quaternion
complex + vector + quaternion -> matrix -> linalg
vector + quaternion + matrix -> geometry
complex + vector + quaternion + matrix -> algebra
```

`geometry` can use `scalar`, `vector`, `quaternion`, and `matrix`.

`algebra` can use `scalar`, `complex`, `vector`, `quaternion`, and `matrix`.

`geometry` and `algebra` must not depend on each other. A later adapter module can combine their concrete types when a real use case requires it.

SIMD code is an implementation detail of scalar, vector, quaternion, and matrix operations. It must not define a second public mathematical model.

## IEEE module

`zmath.ieee` owns exact IEEE-754 bit semantics.

It supports binary16, bfloat16, binary32, binary64, binary80, and binary128 formats. `BFloat16` is an explicit sixteen-bit value type. It is not treated as an alias for a Zig language float type.

Public functions use an explicit format suffix:

```zig
pub fn classifyF32(x: f32) Class;
pub fn classifyF64(x: f64) Class;
pub fn totalOrderF32(a: f32, b: f32) std.math.Order;
pub fn totalOrderF64(a: f64, b: f64) std.math.Order;
pub fn nextAfterF32(x: f32, y: f32) f32;
pub fn nextAfterF64(x: f64, y: f64) f64;
pub fn copySignF32(magnitude: f32, sign: f32) f32;
pub fn copySignF64(magnitude: f64, sign: f64) f64;
```

The module defines `Class` for zero, subnormal, normal, infinity, quiet NaN, and signaling NaN.

This module owns signed zero, NaN payload, subnormal, overflow, underflow, and rounding edge behavior. Higher modules must use it instead of duplicating bit tests.

## Scalar module

`zmath.scalar` provides real elementary functions with separate entry points for every supported format.

```zig
pub fn sinF32(x: f32) f32;
pub fn sinF64(x: f64) f64;
pub fn expF32(x: f32) f32;
pub fn expF64(x: f64) f64;
pub fn logF32(x: f32) f32;
pub fn logF64(x: f64) f64;
```

The suffix rule also applies to binary16, bfloat16, binary80, and binary128 entry points.

### CORE-MATH coverage

The scalar target is the complete published CORE-MATH function table for every format that CORE-MATH supports:

- binary16
- bfloat16
- binary32
- binary64
- binary80
- binary128

The implementation includes all published CORE-MATH real function families. This includes trigonometric, inverse trigonometric, hyperbolic, inverse hyperbolic, exponential, logarithmic, power, root, error, gamma, and pi-scaled functions.

Each exported function documents one accuracy contract:

- correctly rounded for every supported rounding mode, or
- a maximum error in ULP for round-to-nearest-even.

A function must not promise both contracts unless its implementation proves both.

`scalar.constants` exports named constants for every supported format. It does not expose a run-time precision configuration.

## Complex module

`zmath.complex` provides the complete C99 and C23 complex function family.

```zig
pub const Complex32 = complex.Complex32;
pub const Complex64 = complex.Complex64;
pub const Complex80 = complex.Complex80;
pub const Complex128 = complex.Complex128;
```

Each complex type uses explicit `re` and `im` fields. It has no allocator and no run-time precision field.

```zig
pub const Complex64 = struct {
    re: f64,
    im: f64,

    pub fn add(a: Complex64, b: Complex64) Complex64;
    pub fn mul(a: Complex64, b: Complex64) Complex64;
    pub fn div(a: Complex64, b: Complex64) Complex64;
    pub fn conjugate(value: Complex64) Complex64;
    pub fn abs(value: Complex64) f64;
    pub fn arg(value: Complex64) f64;
    pub fn exp(value: Complex64) Complex64;
    pub fn log(value: Complex64) Complex64;
    pub fn sqrt(value: Complex64) Complex64;
    pub fn pow(a: Complex64, b: Complex64) Complex64;
};
```

The full family includes arithmetic, conjugate, projection, absolute value, argument, exponential, logarithmic, power, root, trigonometric, inverse trigonometric, hyperbolic, and inverse hyperbolic functions.

Each branch-sensitive function documents its branch cut, signed-zero behavior, NaN behavior, and infinity behavior.

`Complex32` and `Complex64` are first. `Complex80` and `Complex128` follow the same API when their underlying real scalar contract exists.

## Vector module

`zmath.vector` exports only fixed, named vector types.

```zig
pub const Vec2f = vector.Vec2f;
pub const Vec3f = vector.Vec3f;
pub const Vec4f = vector.Vec4f;
pub const Vec2d = vector.Vec2d;
pub const Vec3d = vector.Vec3d;
pub const Vec4d = vector.Vec4d;
```

Each type is a value type with named components and no allocation.

```zig
pub const Vec4f = struct {
    x: f32,
    y: f32,
    z: f32,
    w: f32,

    pub fn add(a: Vec4f, b: Vec4f) Vec4f;
    pub fn scale(a: Vec4f, factor: f32) Vec4f;
    pub fn dot(a: Vec4f, b: Vec4f) f32;
    pub fn lengthSquared(a: Vec4f) f32;
    pub fn normalized(a: Vec4f) ?Vec4f;
};
```

The API does not export `Vec(comptime dimensions, comptime Scalar)`.

Private helpers can use `@Vector` or narrow comptime helpers. These helpers must reduce implementation duplication. Public callers must not need them to name a type or call an operation.

## Quaternion module

`zmath.quaternion` provides explicit Hamilton quaternion types.

```zig
pub const Quatf = quaternion.Quatf;
pub const Quatd = quaternion.Quatd;
```

```zig
pub const Quatf = struct {
    x: f32,
    y: f32,
    z: f32,
    w: f32,

    pub fn identity() Quatf;
    pub fn add(a: Quatf, b: Quatf) Quatf;
    pub fn mul(a: Quatf, b: Quatf) Quatf;
    pub fn conjugate(value: Quatf) Quatf;
    pub fn inverse(value: Quatf) ?Quatf;
    pub fn normalized(value: Quatf) ?Quatf;
    pub fn rotate(value: Quatf, vector: Vec3f) Vec3f;
};
```

`Quatf` and `Quatd` are semantic rotation types. They are not aliases for `Cl3` rotors.

`Cl3` provides explicit conversion functions such as `rotorFromQuatf` and `quatFromRotorF`.

## Matrix module

`zmath.matrix` exports one structural generic with static storage.

```zig
pub fn Mat(
    comptime rows: usize,
    comptime columns: usize,
    comptime Scalar: type,
) type;

pub const Mat2f = Mat(2, 2, f32);
pub const Mat3f = Mat(3, 3, f32);
pub const Mat4f = Mat(4, 4, f32);
pub const Mat2d = Mat(2, 2, f64);
pub const Mat3d = Mat(3, 3, f64);
pub const Mat4d = Mat(4, 4, f64);
pub const Mat4c = Mat(4, 4, Complex64);
```

`Mat` stores `[rows][columns]Scalar` in one documented storage order. Its dimensions and scalar type are compile-time values. It does not allocate.
`Scalar` must be a zmath real or complex scalar type with the operations required by the selected matrix operation. `Mat` does not accept arbitrary user types.


`Mat` supports element access, construction, transpose, addition, scalar multiplication, and multiplication where dimensions and scalar type match.

Concrete aliases own semantic operations. For example, `Mat3f` owns `determinant`, `inverse`, and `mulVec3`. `Mat4f` owns homogeneous 3D and 4D transforms.

The API does not include a dynamic matrix type, expression templates, or a generic tensor API.

## Linear algebra module

`zmath.linalg` contains factorization and solve algorithms over `Mat` specializations.

```zig
pub fn solveMat2f(a: Mat2f, b: Vec2f) ?Vec2f;
pub fn solveMat3f(a: Mat3f, b: Vec3f) ?Vec3f;
pub fn solveMat4f(a: Mat4f, b: Vec4f) ?Vec4f;
pub fn factorLuMat4f(a: Mat4f) ?LuMat4f;
pub fn factorCholeskyMat4f(a: Mat4f) ?CholeskyMat4f;
```

A function returns `null` for a normal mathematically undefined result. Examples include a singular solve and normalization of a zero vector.

A function returns an error only for a recoverable exceptional condition. Fixed-size operations do not allocate, so they normally do not return allocation errors.

## Geometry module

`zmath.geometry` provides named Euclidean primitives, intersections, and predicates in two, three, and four dimensions.

```zig
pub const Point2f = geometry.Point2f;
pub const Point3f = geometry.Point3f;
pub const Point4f = geometry.Point4f;
pub const Ray3f = geometry.Ray3f;
pub const Ray4f = geometry.Ray4f;
pub const Plane3f = geometry.Plane3f;
pub const Hyperplane4f = geometry.Hyperplane4f;
pub const Sphere3f = geometry.Sphere3f;
pub const Hypersphere4f = geometry.Hypersphere4f;
pub const Aabb3f = geometry.Aabb3f;
pub const Aabb4f = geometry.Aabb4f;
```

A point and a vector use separate types. This prevents accidental point-plus-point operations.

`predicates` owns orientation and containment predicates. It exposes an explicit accuracy contract. It does not silently treat near-zero values as zero.

## Algebra module

`zmath.algebra` provides named foundational algebras. Each algebra has a fixed dimension, metric, basis order, storage layout, and public types.

```zig
pub const Cl3f = algebra.Cl3f;
pub const Cl3d = algebra.Cl3d;
pub const Cl4f = algebra.Cl4f;
pub const Cl4d = algebra.Cl4d;
pub const Pga3f = algebra.Pga3f;
pub const Pga3d = algebra.Pga3d;
pub const Pga4f = algebra.Pga4f;
pub const Pga4d = algebra.Pga4d;
pub const Cga3f = algebra.Cga3f;
pub const Cga3d = algebra.Cga3d;
pub const Cga4f = algebra.Cga4f;
pub const Cga4d = algebra.Cga4d;
```

### Cl3 and Cl4

`Cl3` is `Cl(3, 0, 0)` with eight coefficients.

`Cl4` is `Cl(4, 0, 0)` with sixteen coefficients.

Each module provides fixed-dimension vectors, bivectors, rotors, pseudoscalars, and full multivectors.

### Pga3 and Pga4

`Pga3` is `P(3, 0, 1)` with sixteen coefficients.

`Pga4` is `P(4, 0, 1)` with thirty-two coefficients.

Each module provides fixed-dimension projective points, hyperplanes, lines, flats, motors, and full multivectors.

PGA operations use explicit names for concrete operand pairs. Examples include `joinPoints`, `meetHyperplanes`, and `applyMotorToPoint`.

### Cga3 and Cga4

`Cga3` is `Cl(4, 1, 0)` with thirty-two coefficients.

`Cga4` is `Cl(5, 1, 0)` with sixty-four coefficients.

Each module provides conformal points, spheres, hyperspheres, flats, rounds, versors, and full multivectors.

`Cga3` and `Cga4` are added after the corresponding Cl and PGA modules define their basis, dual, reverse, normalization, and degeneracy conventions.

## Generic policy

`Mat` is the only public structural generic in the first API version.

The library can add a generic function only when it models one stable concept and does not create a type family that callers must inspect.

The API does not expose these forms:

```zig
pub fn Algebra(comptime p: usize, comptime q: usize, comptime r: usize, comptime T: type) type;
pub fn Vec(comptime dimensions: usize, comptime T: type) type;
pub fn Quat(comptime T: type) type;
pub fn product(a: anytype, b: anytype) @TypeOf(a);
```

Private modules can use comptime tables, `@Vector`, or narrow helpers. These tools must not make the public API open-ended.

## Error and allocation policy

Public fixed-size operations do not allocate.

A public operation returns:

- a value for a total operation,
- `?T` for a normal mathematically undefined result,
- a named error set only for a recoverable exceptional condition.

The library does not use `anyerror` in public signatures.

## Oracle and benchmark boundary

`zmath` has no MPFR, GMP, FLINT, Arb, BLIS, Eigen, or C++ dependency.

`~/projects/zmath-testing` owns integration tests, oracle tests, and benchmarks.

- MPFR checks scalar correctly rounded results.
- Arb and FLINT check interval and special-function results.
- External linear algebra libraries provide differential-test and benchmark baselines.
- Each test and benchmark names the exact public function, scalar type, matrix specialization, and concrete algebra.

## Addition order

1. IEEE format primitives and conversions
2. CORE-MATH real scalar coverage
3. `Complex32` and `Complex64` C complex coverage
4. `Vec2f`, `Vec3f`, `Vec4f`, then `f64` variants
5. `Quatf` and `Quatd`
6. `Mat` and named matrix aliases
7. concrete linear-algebra functions
8. Euclidean geometry through four dimensions
9. `Cl3` and `Cl4`
10. `Pga3` and `Pga4`
11. `Cga3` and `Cga4`

A later module must not force an open-ended generic abstraction into an earlier module.
