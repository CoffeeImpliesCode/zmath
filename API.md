# zmath Public API

## Purpose

`zmath` provides fast, dependency-free mathematical primitives for Zig 0.16.

The public API uses explicit scalar types and concrete algebra types. It does not use a public algebra generator.

The library must support compile-time evaluation when a function uses only compile-time inputs. The API must not require an allocator.

## Public root

```zig
const zmath = @import("zmath");

pub const ieee = @import("ieee.zig");
pub const scalar = @import("scalar.zig");
pub const vector = @import("vector.zig");
pub const matrix = @import("matrix.zig");
pub const linalg = @import("linalg.zig");
pub const geometry = @import("geometry.zig");
pub const algebra = @import("algebra.zig");
```

The root module exports only stable public modules and selected concrete types. It does not re-export implementation helpers.

## Module tree

```text
src/
  root.zig

  ieee.zig
  ieee/
    binary32.zig
    binary64.zig
    classify.zig
    order.zig
    rounding.zig
    scale.zig

  scalar.zig
  scalar/
    f32.zig
    f64.zig
    constants.zig
    trig.zig
    exp_log.zig
    power.zig
    special.zig

  vector.zig
  vector/
    vec2.zig
    vec3.zig
    vec4.zig
    simd.zig

  matrix.zig
  matrix/
    mat2.zig
    mat3.zig
    mat4.zig
    affine3.zig

  linalg.zig
  linalg/
    solve2.zig
    solve3.zig
    factor3.zig
    transform.zig

  geometry.zig
  geometry/
    primitives2.zig
    primitives3.zig
    intersections2.zig
    intersections3.zig
    predicates.zig

  algebra.zig
  algebra/
    cl3.zig
    pga3.zig
    cga3.zig
```

The file tree describes the target public shape. A file is added only when its concrete API is implemented.

## Dependency order

```text
ieee -> scalar -> vector -> matrix -> linalg
                   |         |
                   +---------+-> geometry
                   |
                   +------------> algebra
```

`geometry` can use `vector`, `matrix`, and `scalar`.

`algebra` can use `vector`, `matrix`, and `scalar`.

`geometry` and `algebra` must not depend on each other. A later adapter module can combine their concrete types when a real use case requires it.

`simd` is an implementation detail of `scalar`, `vector`, and `matrix`. It must not define a second public mathematical model.

## IEEE module

`zmath.ieee` owns IEEE-754 operations that need exact bit semantics.

Public functions include these families:

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

The module defines `Class` for zero, subnormal, normal, infinity, quiet NaN, and signaling NaN classification.

This module owns signed zero, NaN payload, subnormal, overflow, and rounding edge behavior. Higher modules must use it instead of duplicating bit tests.

## Scalar module

`zmath.scalar` provides scalar mathematical functions with separate `f32` and `f64` entry points.

```zig
pub fn sinF32(x: f32) f32;
pub fn sinF64(x: f64) f64;
pub fn cosF32(x: f32) f32;
pub fn cosF64(x: f64) f64;
pub fn expF32(x: f32) f32;
pub fn expF64(x: f64) f64;
pub fn logF32(x: f32) f32;
pub fn logF64(x: f64) f64;
```

Other families use the same suffix rule:

- `asin`, `acos`, `atan`, `atan2`, `sincos`
- `sinh`, `cosh`, `tanh`
- `exp2`, `expm1`, `log2`, `log10`, `log1p`
- `pow`, `sqrt`, `cbrt`, `hypot`
- `erf`, `erfc`, `gamma`, and other special functions

Each function documents one accuracy contract:

- correctly rounded for every supported rounding mode, or
- a maximum error in ULP for round-to-nearest-even.

A function must not promise both contracts unless its implementation proves both.

`scalar.constants` exports named `f32` and `f64` constants. It does not expose a run-time precision configuration.

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
pub const Vec3f = struct {
    x: f32,
    y: f32,
    z: f32,

    pub fn add(a: Vec3f, b: Vec3f) Vec3f;
    pub fn sub(a: Vec3f, b: Vec3f) Vec3f;
    pub fn scale(a: Vec3f, factor: f32) Vec3f;
    pub fn dot(a: Vec3f, b: Vec3f) f32;
    pub fn cross(a: Vec3f, b: Vec3f) Vec3f;
    pub fn lengthSquared(a: Vec3f) f32;
    pub fn normalized(a: Vec3f) ?Vec3f;
};
```

The API does not export `Vec(comptime dimensions, comptime Scalar)`.

Private helpers can use `@Vector` or small comptime helpers when they remove repeated implementation code. Public callers must not need those helpers to name a type or call an operation.

## Matrix module

`zmath.matrix` exports fixed-size matrices for the dimensions that have stable operations.

```zig
pub const Mat2f = matrix.Mat2f;
pub const Mat3f = matrix.Mat3f;
pub const Mat4f = matrix.Mat4f;
pub const Mat2d = matrix.Mat2d;
pub const Mat3d = matrix.Mat3d;
pub const Mat4d = matrix.Mat4d;
pub const Affine3f = matrix.Affine3f;
pub const Affine3d = matrix.Affine3d;
```

A matrix type specifies its storage order in its documentation. All operations use that order.

Each type owns only operations that match its dimension. For example, `Mat3f` owns `determinant`, `transpose`, `inverse`, and `mulVec3`. `Mat4f` owns homogeneous transform operations.

The API does not export `Mat(comptime rows, comptime columns, comptime Scalar)` in the first public version. Dynamic dense matrices are a separate future package, not a fallback representation for fixed matrices.

## Linear algebra module

`zmath.linalg` contains concrete algorithms that operate on the fixed matrix types.

```zig
pub fn solveMat2f(a: Mat2f, b: Vec2f) ?Vec2f;
pub fn solveMat3f(a: Mat3f, b: Vec3f) ?Vec3f;
pub fn factorLuMat3f(a: Mat3f) ?LuMat3f;
pub fn factorCholeskyMat3f(a: Mat3f) ?CholeskyMat3f;
```

A function returns `null` for a mathematically invalid input when the invalid state is part of normal use. Examples include a singular solve or normalization of a zero vector.

A function returns an error only when the caller can recover from an external condition. Fixed-size operations do not allocate, so they normally do not return allocation errors.

## Geometry module

`zmath.geometry` provides named Euclidean primitives and predicates.

```zig
pub const Point2f = geometry.Point2f;
pub const Point3f = geometry.Point3f;
pub const Ray3f = geometry.Ray3f;
pub const Plane3f = geometry.Plane3f;
pub const Sphere3f = geometry.Sphere3f;
pub const Aabb3f = geometry.Aabb3f;
```

A point and a vector use separate types. This prevents accidental point-plus-point operations.

`predicates` owns orientation and containment predicates. It exposes an explicit accuracy contract. It does not silently treat near-zero values as zero.

## Algebra module

`zmath.algebra` provides named foundational algebras. Each algebra has fixed dimension, metric, basis order, storage, and public types.

```zig
pub const Cl3f = algebra.Cl3f;
pub const Cl3d = algebra.Cl3d;
pub const Pga3f = algebra.Pga3f;
pub const Pga3d = algebra.Pga3d;
pub const Cga3f = algebra.Cga3f;
pub const Cga3d = algebra.Cga3d;
```

### Cl3

`Cl3` is `Cl(3, 0, 0)` with eight coefficients.

It provides Euclidean vectors, bivectors, rotors, pseudoscalars, and full multivectors.

```zig
pub const Cl3f = struct {
    pub const Vector = Vec3f;
    pub const Bivector = Bivector3f;
    pub const Rotor = Rotor3f;
    pub const Multivector = Cl3Multivectorf;

    pub fn geometricProduct(a: Multivector, b: Multivector) Multivector;
    pub fn outerProduct(a: Multivector, b: Multivector) Multivector;
    pub fn reverse(value: Multivector) Multivector;
    pub fn dual(value: Multivector) Multivector;
};
```

### Pga3

`Pga3` is `P(3, 0, 1)` with sixteen coefficients.

It provides projective points, planes, lines, motors, and full multivectors.

```zig
pub const Pga3f = struct {
    pub const Point = PgaPoint3f;
    pub const Plane = PgaPlane3f;
    pub const Line = PgaLine3f;
    pub const Motor = PgaMotor3f;
    pub const Multivector = Pga3Multivectorf;

    pub fn joinPoints(a: Point, b: Point) Line;
    pub fn meetPlanes(a: Plane, b: Plane) Line;
    pub fn applyMotorToPoint(motor: Motor, point: Point) Point;
};
```

Zig has no overloads. The public API uses explicit named functions for each concrete operand pair.

### Cga3

`Cga3` is `Cl(4, 1, 0)` with thirty-two coefficients.

It provides conformal points, spheres, circles, lines, planes, versors, and full multivectors.

`Cga3` is added after `Cl3` and `Pga3` define their basis, dual, reverse, normalization, and degeneracy conventions.

## Generic policy

The public API uses generics only where a single type parameter expresses one stable concept without creating a type family that callers must inspect.

Allowed public generic examples:

```zig
pub fn lerp(comptime T: type, a: T, b: T, t: T) T;
pub fn clamp(comptime T: type, value: T, min: T, max: T) T;
```

These functions are added only after their valid type set is explicit and tested.

Disallowed public generic examples:

```zig
pub fn Algebra(comptime p: usize, comptime q: usize, comptime r: usize, comptime T: type) type;
pub fn Vec(comptime dimensions: usize, comptime T: type) type;
pub fn Mat(comptime rows: usize, comptime columns: usize, comptime T: type) type;
pub fn product(a: anytype, b: anytype) @TypeOf(a);
```

Private modules can use comptime tables, `@Vector`, or narrowly scoped helpers. These tools must reduce implementation duplication. They must not make the public API open-ended.

## Error and allocation policy

Public fixed-size operations do not allocate.

A public operation returns:

- a value for total operations,
- `?T` for a normal mathematically undefined result,
- a named error set only for a recoverable exceptional condition.

The library does not use `anyerror` in public signatures.

## Oracle and benchmark boundary

`zmath` has no MPFR, GMP, FLINT, Arb, BLIS, Eigen, or C++ dependency.

`~/projects/zmath-testing` owns integration tests, oracle tests, and benchmarks.

- MPFR checks scalar correctly-rounded results.
- Arb and FLINT check interval and special-function results.
- External linear algebra libraries provide differential-test and benchmark baselines.
- Each test and benchmark names the exact public function, scalar type, and concrete algebra.

## Addition order

1. `ieee`
2. `scalar` basic functions
3. `Vec2f`, `Vec3f`, `Vec4f`, then `f64` variants
4. `Mat2f`, `Mat3f`, `Mat4f`, then `f64` variants
5. concrete linear-algebra functions
6. Euclidean geometry and predicates
7. `Cl3`
8. `Pga3`
9. `Cga3`

A later module must not force a generic abstraction into an earlier module.
