# zmath API Contract

## Status

This document defines the target API and identifies the currently implemented public subset.

The current public root exports `add`, `printAnotherMessage`, `ieee`, and `scalar`. The implemented `ieee` and `scalar` symbols are limited to the subsets listed below. Other IEEE and scalar functions remain planned.

## Scope

`zmath` is a dependency-free Zig 0.16 mathematics library. Public fixed-size operations do not allocate.

Public real formats are `f16`, `f32`, and `f64` for the current IEEE and scalar primitives. Public type and function names use lower-case `h`, `f`, and `d` suffixes. The `q` suffix remains reserved for a future `f128` API. No bfloat16 API exists.

The API uses one structural generic, `Mat`. Vectors, quaternions, geometry, factorizations, and geometric algebras use named concrete types.

## Planned module paths

The module path is canonical. Root re-exports are optional convenience names after an implementation becomes public.

```text
zmath.ieee
zmath.scalar
zmath.complex
zmath.vector
zmath.quaternion
zmath.matrix
zmath.linalg
zmath.geometry
zmath.algebra
```

```text
src/
  ieee.zig
  scalar.zig
  complex.zig
  vector.zig
  quaternion.zig
  matrix.zig
  linalg.zig
  geometry.zig
  algebra.zig
```

The dependency direction is:

```text
ieee -> scalar -> complex
             +-> vector -> quaternion
complex + vector + quaternion -> matrix -> linalg
vector + quaternion + matrix -> geometry
complex + vector + quaternion + matrix -> algebra
```

Geometry and algebra do not depend on each other.

## IEEE and scalar API

The complete planned `zmath.ieee` API owns classification, total order, bit conversion, successor, predecessor, scaling, and rounding for the public real formats.

`h` is implemented for the current IEEE and scalar primitives. These primitives implement the `h`, `f`, and `d` suffix families:

```zig
pub const RoundMode = enum {
    nearest_even,
    toward_zero,
    toward_positive,
    toward_negative,
};

pub const Class = enum {
    zero,
    subnormal,
    normal,
    infinity,
    quiet_nan,
    signaling_nan,
};

pub fn classifyh(value: f16) Class;
pub fn classifyf(value: f32) Class;
pub fn classifyd(value: f64) Class;
pub fn copySignh(value: f16, sign: f16) f16;
pub fn copySignf(value: f32, sign: f32) f32;
pub fn copySignd(value: f64, sign: f64) f64;
pub fn nextAfterh(from: f16, to: f16) f16;
pub fn nextAfterf(from: f32, to: f32) f32;
pub fn nextAfterd(from: f64, to: f64) f64;
pub fn totalOrderh(left: f16, right: f16) bool;
pub fn totalOrderf(left: f32, right: f32) bool;
pub fn totalOrderd(left: f64, right: f64) bool;
```

`nextAfterh`, `nextAfterf`, and `nextAfterd` choose the first NaN operand. If it is a signaling NaN, they set its quiet bit and preserve its sign and payload. They return a quiet NaN unchanged.

`totalOrderh`, `totalOrderf`, and `totalOrderd` are Boolean predicates. They return `true` when `left` precedes or equals `right` in the IEEE total order.

The implemented public scalar subset is:

```zig
pub fn absh(value: f16) f16;
pub fn absf(value: f32) f32;
pub fn absd(value: f64) f64;
pub fn ulpDistanceh(a: f16, b: f16) u16;
pub fn ulpDistancef(a: f32, b: f32) u32;
pub fn ulpDistanced(a: f64, b: f64) u64;
```

A planned scalar operation takes an explicit `RoundMode` when its result depends on rounding. The API does not read or change the process floating-point environment.

The CORE-MATH manifest-gated transcendental functions below remain planned. `zmath` does not currently export them.

```zig
pub fn sinh(mode: RoundMode, x: f16) f16;
pub fn sinf(mode: RoundMode, x: f32) f32;
pub fn sind(mode: RoundMode, x: f64) f64;
pub fn exph(mode: RoundMode, x: f16) f16;
pub fn expf(mode: RoundMode, x: f32) f32;
pub fn expd(mode: RoundMode, x: f64) f64;
pub fn logh(mode: RoundMode, x: f16) f16;
pub fn logf(mode: RoundMode, x: f32) f32;
pub fn logd(mode: RoundMode, x: f64) f64;
```

### CORE-MATH manifest

Before any manifest-gated transcendental function becomes public, the repository adds a versioned manifest at `spec/core-math-<commit>.md`.

The manifest pins one CORE-MATH commit. It lists every required f32 and f64 function-format pair, its zmath symbol, domain, special values, supported rounding modes, and accuracy class.

A manifest entry is one of:

- `correct`: the result is correctly rounded in the requested `RoundMode`.
- `bounded`: the result has a stated maximum ULP error in `nearest_even` mode only.
- `unsupported`: zmath does not export the pair.

A ULP comparison applies only to finite results of the same sign. The manifest defines the exact expected bit pattern for zero, subnormal, infinity, and NaN results. A NaN comparison checks class, sign when required, and payload policy when required.

## Complex API

`zmath.complex` owns the canonical `Complexf` and `Complexd` types.

```zig
pub const Complexf = struct { re: f32, im: f32 };
pub const Complexd = struct { re: f64, im: f64 };
```

The complete C99/C23 complex surface is required for both types. The planned module functions use an explicit type suffix:

```zig
pub fn realf(value: Complexf) f32;
pub fn reald(value: Complexd) f64;
pub fn imagf(value: Complexf) f32;
pub fn imagd(value: Complexd) f64;
pub fn conjf(value: Complexf) Complexf;
pub fn conjd(value: Complexd) Complexd;
pub fn projf(value: Complexf) Complexf;
pub fn projd(value: Complexd) Complexd;
pub fn absf(value: Complexf) f32;
pub fn absd(value: Complexd) f64;
pub fn argf(value: Complexf) f32;
pub fn argd(value: Complexd) f64;
```

For both suffixes, the manifest enumerates `add`, `sub`, `mul`, `div`, `exp`, `log`, `log10`, `pow`, `sqrt`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, and `atanh`.

Each entry defines its Zig name, inputs, result, branch cut, principal range, signed-zero result, NaN result, infinity result, pole result, overflow result, and underflow result.

Finite complex results use the scalar manifest accuracy contract for each component after range reduction. A branch-cut or special-value case uses an exact expected component-bit pattern. Complex functions do not return `null` for C mathematical poles.

## Vectors and quaternions

`zmath.vector` owns `Vec2f`, `Vec3f`, `Vec4f`, `Vec2d`, `Vec3d`, and `Vec4d`. Each is a value type with named components.

`zmath.quaternion` owns raw Hamilton types and validated rotation types:

```zig
pub const Quatf = struct { x: f32, y: f32, z: f32, w: f32 };
pub const Quatd = struct { x: f64, y: f64, z: f64, w: f64 };
pub const UnitQuatf = struct { raw: Quatf };
pub const UnitQuatd = struct { raw: Quatd };
```

`w` is the scalar component. `Quatf.mul(a, b)` is Hamilton multiplication. A `UnitQuatf` represents an active, right-handed rotation of a column vector.

`UnitQuatf.rotate(v)` applies `q * (v, 0) * conjugate(q)`. `a.then(b)` means apply `a` first, then `b`, and equals `b * a`.

`raw` is private to `zmath.quaternion`. `UnitQuatf.init(raw)` returns `?UnitQuatf`. It rejects zero, non-finite, and non-unit values after normalization. Raw `Quatf` supports algebraic operations. It does not expose a total rotation operation.

## Matrices

`zmath.matrix` owns the sole public structural generic.

```zig
pub fn Mat(
    comptime rows: usize,
    comptime columns: usize,
    comptime Scalar: type,
) type;

pub const Mat2f = Mat(2, 2, f32);
pub const Mat3f = Mat(3, 3, f32);
pub const Mat4f = Mat(4, 4, f32);
pub const Mat5f = Mat(5, 5, f32);
pub const Mat2d = Mat(2, 2, f64);
pub const Mat3d = Mat(3, 3, f64);
pub const Mat4d = Mat(4, 4, f64);
pub const Mat5d = Mat(5, 5, f64);
```

`Scalar` is one of `f32`, `f64`, `Complexf`, or `Complexd`. Other user-defined scalar types are unsupported.

`Mat(rows, columns, Scalar)` stores `data: [columns][rows]Scalar`. The first index is the column. The second index is the row. `at(row, column)` returns `data[column][row]`.

Vectors are columns. `M * v` applies matrix `M` to vector `v`. Composition is right-to-left. `(A * B) * v` applies `B` first and `A` second.

```zig
pub fn mul(
    comptime output_columns: usize,
    left: Mat(rows, columns, Scalar),
    right: Mat(columns, output_columns, Scalar),
) Mat(rows, output_columns, Scalar);
```

Specialized operations are module functions, not alias methods.

```zig
pub fn determinantMat3f(value: Mat3f) f32;
pub fn inverseMat4d(value: Mat4d) ?Mat4d;
pub fn mulMat4fVec4f(matrix: Mat4f, vector: Vec4f) Vec4f;
```

`Mat4f` and `Mat4d` represent homogeneous 3D transforms and linear 4D transforms. `Mat5f` and `Mat5d` represent homogeneous 4D affine or projective transforms.

The API has no dynamic matrix, expression template, or tensor type.

## Linear algebra

`zmath.linalg` owns named factorization types for named matrix sizes.

```zig
pub const LuMat2f = struct { lu: Mat2f, permutation: [2]u8, parity: i1 };
pub const LuMat3f = struct { lu: Mat3f, permutation: [3]u8, parity: i1 };
pub const LuMat4f = struct { lu: Mat4f, permutation: [4]u8, parity: i1 };
pub const CholeskyMat4f = struct { lower: Mat4f };
```

Each factorization type has a same-size `solve`, `reconstruct`, and `determinant` operation. A factorization returns `null` for a singular, non-positive-definite, or non-finite input, as applicable.

## Geometry

`zmath.geometry` owns concrete f/d families in each dimension.

- 2D: point, vector, line, ray, segment, circle, and AABB.
- 3D: point, vector, line, ray, segment, plane, sphere, and AABB.
- 4D: point, vector, line, ray, segment, hyperplane, hypersphere, and AABB.

A point and a vector are distinct types. A line, ray, plane, or hyperplane constructor normalizes its direction or normal. A circle, sphere, or hypersphere requires a finite nonnegative radius. An AABB requires `minimum[i] <= maximum[i]`.

Each dimension and precision has a concrete intersection union. It distinguishes `none`, `point`, `two_points`, `segment`, `line`, `plane`, `hyperplane`, `coincident`, and `contained` when applicable. It does not use `null` to combine these cases.

Predicates use one explicit policy:

```zig
pub const PredicatePolicy = enum { fast, adaptive_exact };
```

`fast` reports the sign of the scalar computation. `adaptive_exact` reports the exact sign for finite f32 or f64 inputs. Predicates do not use an implicit epsilon. Every geometry operation rejects non-finite input with `error.NonFiniteInput`.

## Concrete geometric algebras

`zmath.algebra` owns `Cl3f/d`, `Cl4f/d`, `Pga3f/d`, `Pga4f/d`, `Cga3f/d`, and `Cga4f/d`.

| Module | Ordered basis | Metric | Coefficients |
| --- | --- | --- | --- |
| Cl3 | e1, e2, e3 | +1, +1, +1 | 8 |
| Cl4 | e1, e2, e3, e4 | +1, +1, +1, +1 | 16 |
| Pga3 | e0, e1, e2, e3 | 0, +1, +1, +1 | 16 |
| Pga4 | e0, e1, e2, e3, e4 | 0, +1, +1, +1, +1 | 32 |
| Cga3 | e1, e2, e3, e_plus, e_minus | +1, +1, +1, +1, -1 | 32 |
| Cga4 | e1, e2, e3, e4, e_plus, e_minus | +1, +1, +1, +1, +1, -1 | 64 |

The basis order assigns a bit to each basis vector. A full multivector stores coefficients in ascending unsigned blade-mask order. The coefficient at index `mask` multiplies the ascending canonical blade represented by `mask`.

Reverse maps grade `k` to `(-1)^(k * (k - 1) / 2)`. Cl and CGA use the right dual `A * inverse(I)`. PGA uses the signed complementary Poincare dual.

For CGA, `n_o = (e_minus - e_plus) / 2` and `n_inf = e_minus + e_plus`. Their inner product is `-1`. A Euclidean vector `x` embeds as `n_o + x + (x dot x) * n_inf / 2`.

Each algebra module publishes a normalization table before it becomes public. The table defines scale, sign, preconditions, and failure behavior for every named rotor, motor, versor, point, flat, round, and full multivector operation.

No public `Algebra(p, q, r, T)`, generic vector, generic quaternion, or `anytype` product exists.

## Tests and benchmarks

`zmath` has no MPFR, GMP, FLINT, Arb, BLIS, Eigen, or C++ runtime dependency.

`~/projects/zmath-testing` owns integration tests, oracle tests, and benchmarks. The current oracle confirms its zmath import and external-library environment only. It is not a numerical differential test.

A numerical oracle test is added with each public scalar, complex, matrix, geometry, or algebra contract. It names the exact zmath symbol, type, manifest entry, input corpus, expected result, and acceptance rule.

An external differential test or benchmark names and explicitly configures its baseline library, version, build options, scalar type, dimensions, and operation. The current scaffold has no configured external baseline.
