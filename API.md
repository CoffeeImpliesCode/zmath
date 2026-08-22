# zmath API Contract

## Status

This document defines the target API and identifies the currently implemented public subset.

The current public root exports `add`, `printAnotherMessage`, `ieee`, and `scalar`. The implemented `ieee` and `scalar` symbols are limited to the subsets listed below. Other IEEE and scalar functions remain planned.

## Scope

`zmath` is a dependency-free Zig 0.16 mathematics library. Public fixed-size operations do not allocate.

Public real formats are `f16`, `f32`, and `f64` for the current IEEE and scalar primitives. Public type and function names use lower-case `h`, `f`, and `d` suffixes. `q`, `f128`, and bfloat16 are outside this contract and have no public API.

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

### Status terms

This contract uses three independent status terms:

- **Target** names a planned public symbol.
- **Implemented** means that the current `zmath` root or module exports the symbol.
- **Manifest-accepted** means that a `correct` or `bounded` manifest row accepts one operation-format pair.

A target symbol is not necessarily implemented. An implemented symbol is not necessarily manifest-accepted.

The current root exports the `zmath.ieee` and `zmath.scalar` modules. A target `zmath.<op>` generic is a separate root export and is implemented only when it appears in an implementation-source snapshot.

The implemented `zmath.ieee` exports are:

```zig
pub const RoundMode = enum { nearest_even, toward_zero, toward_positive, toward_negative };
pub const Class = enum { zero, subnormal, normal, infinity, quiet_nan, signaling_nan };
pub fn toBitsh(value: f16) u16;
pub fn fromBitsh(bits: u16) f16;
pub fn classifyh(value: f16) Class;
pub fn classifyf(value: f32) Class;
pub fn classifyd(value: f64) Class;
pub fn copySignh(value: f16, sign: f16) f16;
pub fn copySignf(value: f32, sign: f32) f32;
pub fn copySignd(value: f64, sign: f64) f64;
pub fn nextUph(value: f16) f16;
pub fn nextDownh(value: f16) f16;
pub fn nextAfterh(from: f16, to: f16) f16;
pub fn nextAfterf(from: f32, to: f32) f32;
pub fn nextAfterd(from: f64, to: f64) f64;
pub fn scaleBh(mode: RoundMode, value: f16, exponent: i32) f16;
pub fn roundIntegralh(mode: RoundMode, value: f16) f16;
pub fn totalOrderh(left: f16, right: f16) bool;
pub fn totalOrderf(left: f32, right: f32) bool;
pub fn totalOrderd(left: f64, right: f64) bool;
```

The implemented `zmath.scalar` exports are:

```zig
pub const logd = logd_impl.logd;
pub const sind = @import("scalar_sind.zig").sind;
pub const cosd = @import("scalar_sind.zig").cosd;
pub fn exph(value: f16) f16;
pub fn logh(value: f16) f16;
pub fn sinh(value: f16) f16;
pub fn cosh(value: f16) f16;
pub fn logf(value: f32) f32;
pub fn absh(value: f16) f16;
pub fn absf(value: f32) f32;
pub fn absd(value: f64) f64;
pub fn ulpDistanceh(a: f16, b: f16) u16;
pub fn ulpDistancef(a: f32, b: f32) u32;
pub fn ulpDistanced(a: f64, b: f64) u64;
pub fn expf(value: f32) f32;
pub fn expd(value: f64) f64;
pub fn sinf(value: f32) f32;
pub fn cosf(value: f32) f32;
```

These are source exports. They do not imply a target generic root export or manifest acceptance.

The implemented target scalar families are:

```text
abs{h,f,d}
ulpDistance{h,f,d}
exp{h,f,d}
log{h,f,d}
sin{h,f,d}
cos{h,f,d}
```

`ulpDistance` is a zmath extension. It is not a C23 scalar operation.

This implemented-symbol list does not grant manifest acceptance. The versioned corpus manifest is the only acceptance authority.


### Generic dispatch

The target root API uses one generic operation name for `f16`, `f32`, and `f64`.

```zig
zmath.abs(value)
zmath.sin(value)
zmath.sinh(value)
zmath.atan2(y, x)
zmath.compoundn(value, exponent)
```

The floating input type selects one exact scalar implementation:

| Input type | Scalar suffix |
| --- | --- |
| `f16` | `h` |
| `f32` | `f` |
| `f64` | `d` |

For example, `zmath.sin(x)` selects `zmath.scalar.sinh`, `sinf`, or `sind`.
`zmath.sinh(x)` selects `zmath.scalar.sinhh`, `sinhf`, or `sinhd`.

The scalar name is the operation name followed by the type suffix. The suffix rule produces `sinhh`, `coshh`, `tanhh`, `asinhh`, `acoshh`, and `atanhh` for binary16 hyperbolic operations.

All floating arguments to a normal binary or ternary generic must have the same type. Dispatch does not widen, narrow, or coerce a floating argument.

The first floating argument selects the family when an operation also takes an integer or an output pointer. Integer-result operations keep the suffix of their floating input.

`zmath.abs` is the C23 `fabs` equivalent. Its scalar symbols remain `absh`, `absf`, and `absd`.
The API does not add `fabs`, `fabsh`, `fabsf`, or `fabsd` aliases.

The target has no `q`, `f128`, bfloat16, or bfloat suffix family. These types do not participate in generic dispatch.

### Rounding environment

Callers must keep IEEE round-to-nearest, ties-to-even active for each scalar call. Scalar functions do not take a general `RoundMode` argument.

An operation-defined direction parameter, such as the `fromfp` direction, is part of that operation. It does not change the process environment.

`zmath` does not read, validate, set, or restore the process floating-point environment (`fenv`). Behavior with another active rounding mode is unsupported.

This caller precondition also applies to narrowing operations. They compute the exact mathematical operation and round once to the declared destination type.

### Target C23 scalar families

In the tables below, `op{h,f,d}` names three exact `zmath.scalar` symbols. For example, `exp{h,f,d}` means `exph`, `expf`, and `expd`.

Each operation also has one target root generic named `zmath.<op>`, unless a table states an exception. These tables define target names only.

#### Basic and rounding

| Root generic | Scalar implementation family |
| --- | --- |
| `abs` | `abs{h,f,d}` |
| `ceil` | `ceil{h,f,d}` |
| `floor` | `floor{h,f,d}` |
| `nearbyint` | `nearbyint{h,f,d}` |
| `rint` | `rint{h,f,d}` |
| `round` | `round{h,f,d}` |
| `roundeven` | `roundeven{h,f,d}` |
| `trunc` | `trunc{h,f,d}` |
| `lrint` | `lrint{h,f,d}` |
| `llrint` | `llrint{h,f,d}` |
| `lround` | `lround{h,f,d}` |
| `llround` | `llround{h,f,d}` |
| `fromfp` | `fromfp{h,f,d}` |
| `ufromfp` | `ufromfp{h,f,d}` |
| `fromfpx` | `fromfpx{h,f,d}` |
| `ufromfpx` | `ufromfpx{h,f,d}` |

#### Exponential and logarithmic

| Root generic | Scalar implementation family |
| --- | --- |
| `exp` | `exp{h,f,d}` |
| `exp2` | `exp2{h,f,d}` |
| `exp10` | `exp10{h,f,d}` |
| `expm1` | `expm1{h,f,d}` |
| `exp2m1` | `exp2m1{h,f,d}` |
| `exp10m1` | `exp10m1{h,f,d}` |
| `log` | `log{h,f,d}` |
| `log2` | `log2{h,f,d}` |
| `log10` | `log10{h,f,d}` |
| `log1p` | `log1p{h,f,d}` |
| `logp1` | `logp1{h,f,d}` |
| `log2p1` | `log2p1{h,f,d}` |
| `log10p1` | `log10p1{h,f,d}` |

`log1p` and `logp1` are distinct C23 public spellings with the same mathematical operation.

#### Trigonometric, inverse, and pi-scaled

| Root generic | Scalar implementation family |
| --- | --- |
| `sin` | `sin{h,f,d}` |
| `cos` | `cos{h,f,d}` |
| `tan` | `tan{h,f,d}` |
| `asin` | `asin{h,f,d}` |
| `acos` | `acos{h,f,d}` |
| `atan` | `atan{h,f,d}` |
| `atan2` | `atan2{h,f,d}` |
| `sinpi` | `sinpi{h,f,d}` |
| `cospi` | `cospi{h,f,d}` |
| `tanpi` | `tanpi{h,f,d}` |
| `asinpi` | `asinpi{h,f,d}` |
| `acospi` | `acospi{h,f,d}` |
| `atanpi` | `atanpi{h,f,d}` |
| `atan2pi` | `atan2pi{h,f,d}` |

#### Hyperbolic and inverse hyperbolic

| Root generic | Scalar implementation family |
| --- | --- |
| `sinh` | `sinh{h,f,d}` |
| `cosh` | `cosh{h,f,d}` |
| `tanh` | `tanh{h,f,d}` |
| `asinh` | `asinh{h,f,d}` |
| `acosh` | `acosh{h,f,d}` |
| `atanh` | `atanh{h,f,d}` |

The family notation expands mechanically. Thus, `sinh{h,f,d}` means `sinhh`, `sinhf`, and `sinhd`.

#### Power, root, and norm

| Root generic | Scalar implementation family |
| --- | --- |
| `pow` | `pow{h,f,d}` |
| `pown` | `pown{h,f,d}` |
| `powr` | `powr{h,f,d}` |
| `compoundn` | `compoundn{h,f,d}` |
| `sqrt` | `sqrt{h,f,d}` |
| `cbrt` | `cbrt{h,f,d}` |
| `rsqrt` | `rsqrt{h,f,d}` |
| `rootn` | `rootn{h,f,d}` |
| `hypot` | `hypot{h,f,d}` |

`compoundn`, `pown`, and `rootn` take a signed integer exponent or degree. The target fixed-width Zig integer is `i64`.

The C23 `compoundn(x, n)` operation computes `(1 + x)^n`. A CORE-MATH `compound` source with a floating exponent is not a direct implementation of this operation.

#### Error and gamma

| Root generic | Scalar implementation family |
| --- | --- |
| `erf` | `erf{h,f,d}` |
| `erfc` | `erfc{h,f,d}` |
| `lgamma` | `lgamma{h,f,d}` |
| `tgamma` | `tgamma{h,f,d}` |

`lgamma` returns `LgammaResult(T) { value: T, sign: i2 }`. The sign field is `-1`, `0`, or `1`.
The root generic returns the same result shape for its input type.

The API has no global `signgam` value and no hidden mutable gamma-sign state.

#### Remainder, decomposition, and scaling

| Root generic | Scalar implementation family |
| --- | --- |
| `fmod` | `fmod{h,f,d}` |
| `remainder` | `remainder{h,f,d}` |
| `remquo` | `remquo{h,f,d}` |
| `frexp` | `frexp{h,f,d}` |
| `modf` | `modf{h,f,d}` |
| `ilogb` | `ilogb{h,f,d}` |
| `llogb` | `llogb{h,f,d}` |
| `logb` | `logb{h,f,d}` |
| `ldexp` | `ldexp{h,f,d}` |
| `scalbn` | `scalbn{h,f,d}` |
| `scalbln` | `scalbln{h,f,d}` |

`frexp`, `modf`, and `remquo` keep their C23 decomposition semantics. Their target Zig result containers replace C output-pointer syntax.

#### Fused operation, difference, maximum, and minimum

| Root generic | Scalar implementation family |
| --- | --- |
| `fma` | `fma{h,f,d}` |
| `fdim` | `fdim{h,f,d}` |
| `fmax` | `fmax{h,f,d}` |
| `fmin` | `fmin{h,f,d}` |
| `fmaximum` | `fmaximum{h,f,d}` |
| `fminimum` | `fminimum{h,f,d}` |
| `fmaximum_mag` | `fmaximum_mag{h,f,d}` |
| `fminimum_mag` | `fminimum_mag{h,f,d}` |
| `fmaximum_num` | `fmaximum_num{h,f,d}` |
| `fminimum_num` | `fminimum_num{h,f,d}` |
| `fmaximum_mag_num` | `fmaximum_mag_num{h,f,d}` |
| `fminimum_mag_num` | `fminimum_mag_num{h,f,d}` |

#### Narrowing operations

Narrowing operations are the exception to one-operation-name dispatch. C23 encodes the source and destination format in each IEC 60559 symbol.

C23 7.12.14 defines these public binary-format signatures:

```zig
pub fn fadd(x: f64, y: f64) f32;
pub fn fsub(x: f64, y: f64) f32;
pub fn fmul(x: f64, y: f64) f32;
pub fn fdiv(x: f64, y: f64) f32;
pub fn ffma(x: f64, y: f64, z: f64) f32;
pub fn fsqrt(x: f64) f32;
```

The target also includes the C23 Annex F `_FloatM` forms for the supported public formats:

| Destination | Source | Add | Subtract | Multiply | Divide | Fused multiply-add | Square root |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `f16` | `f32` | `f16addf32(x: f32, y: f32) f16` | `f16subf32(x: f32, y: f32) f16` | `f16mulf32(x: f32, y: f32) f16` | `f16divf32(x: f32, y: f32) f16` | `f16fmaf32(x: f32, y: f32, z: f32) f16` | `f16sqrtf32(x: f32) f16` |
| `f16` | `f64` | `f16addf64(x: f64, y: f64) f16` | `f16subf64(x: f64, y: f64) f16` | `f16mulf64(x: f64, y: f64) f16` | `f16divf64(x: f64, y: f64) f16` | `f16fmaf64(x: f64, y: f64, z: f64) f16` | `f16sqrtf64(x: f64) f16` |
| `f32` | `f64` | `f32addf64(x: f64, y: f64) f32` | `f32subf64(x: f64, y: f64) f32` | `f32mulf64(x: f64, y: f64) f32` | `f32divf64(x: f64, y: f64) f32` | `f32fmaf64(x: f64, y: f64, z: f64) f32` | `f32sqrtf64(x: f64) f32` |

These names are target scalar implementation names. They are not aliases for a same-type operation.
The target has no `f128`, `q`, bfloat16, decimal, or long-double source or destination form.

The legacy `fadd`, `fsub`, `fmul`, `fdiv`, `ffma`, and `fsqrt` forms are root-level C23 names with the signatures above.
The format-pair forms are available from `zmath.scalar` and use the exact names in the table.

Each narrowing operation computes as if to infinite precision and rounds once to the declared destination type.

### C23 operations owned by `zmath.ieee`

The scalar module does not duplicate representation and successor operations. The target `zmath.ieee` module owns these C23 equivalents:

```text
canonicalize{h,f,d}
copySign{h,f,d}
nan{h,f,d}
nextAfter{h,f,d}
nextToward{h,f,d}
nextUp{h,f,d}
nextDown{h,f,d}
totalOrder{h,f,d}
```

Classification predicates, signaling tests, ordered comparisons, bit conversion, and total order also remain in `zmath.ieee`.

### Optional IEC 60559 extension family

C23 Annex F defines these optional IEC 60559 operations. They are target names, but they are not baseline C23 7.12 scalar names:

| Root generic | `zmath.ieee` implementation family |
| --- | --- |
| `totalOrderMag` | `totalOrderMag{h,f,d}` |
| `getPayload` | `getPayload{h,f,d}` |
| `setPayload` | `setPayload{h,f,d}` |
| `setPayloadSig` | `setPayloadSig{h,f,d}` |

No symbol in this optional family is currently implemented or manifest-accepted.

### Exact special-value policy

Each manifest-accepted operation-format pair must define exact output bits for zeros, subnormals, infinities, and NaNs.

The format encodings are:

| Format | `+0` | `-0` | `+1` | `+infinity` | `-infinity` | Canonical quiet NaN |
| --- | --- | --- | --- | --- | --- | --- |
| `h` | `0000` | `8000` | `3c00` | `7c00` | `fc00` | `7e00` |
| `f` | `00000000` | `80000000` | `3f800000` | `7f800000` | `ff800000` | `7fc00000` |
| `d` | `0000000000000000` | `8000000000000000` | `3ff0000000000000` | `7ff0000000000000` | `fff0000000000000` | `7ff8000000000000` |

For a NaN input, preserve the sign and payload and set the format quiet bit.
The quiet-bit masks are `0200`, `00400000`, and `0008000000000000` for `h`, `f`, and `d`.
A quiet NaN remains unchanged. A signaling NaN becomes a quiet NaN.

Unless an operation has a stricter C23 rule, the first NaN argument in parameter order selects the result.
An invalid non-NaN input returns the positive canonical quiet NaN for its format.
An operation-specific C23 rule takes precedence. This includes the number-preferring maximum and minimum families.

Signed-zero and infinity results follow the operation-specific C23 and IEC 60559 rules.
Implementations must not replace those results with a host-library choice.

The manifest must record every operation-specific exception before acceptance.
Target status or source availability cannot substitute for this evidence.

### CORE-MATH manifest

Before a manifest-gated scalar pair becomes accepted, the repository adds or updates a versioned contract at `spec/core-math-<date>.md`.

The machine-readable manifest pins one CORE-MATH commit. It lists each accepted or unsupported function-format pair separately for `f16`, `f32`, and `f64`.
At this contract revision, the schema version 1 manifest has nine `exp`, `log`, and `sin` rows.
All nine rows have status `unsupported`. Other target pairs have no manifest acceptance in schema version 1.

Each entry defines the zmath symbol, domain, special values, supported rounding policy, and accuracy class.

A manifest entry is one of:

- `correct`: the result is correctly rounded in `nearest_even`.
- `bounded`: the result has a stated maximum ULP error in `nearest_even`.
- `unsupported`: the manifest has no accepted accuracy claim for the pair. It does not state that an implementation is absent.

A ULP comparison applies only to finite results of the same sign. Special-value rows require exact expected bit patterns.

CORE-MATH source availability does not prove API implementation or manifest acceptance. Each format requires independent implementation and oracle evidence.

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
