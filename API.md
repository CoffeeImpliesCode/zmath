# zmath API Contract

## Status

This document defines the target API and identifies the currently implemented public subset.

The current public root exports the `zmath.ieee`, `zmath.scalar`, and `zmath.vector` modules, their shared result types, and the C23 source surface defined in `spec/core-math-2026-08-22.md`.
It also exports `add` and `printAnotherMessage`.
The `zmath.complex`, `zmath.quaternion`, `zmath.matrix`, `zmath.linalg`, `zmath.geometry`, and `zmath.algebra` module paths below remain planned.
`zmath.vector` is implemented as the SIMD `@Vector` layer described in [Vector API](#vector-api); it is distinct from the planned linear-algebra `Vec*` value types in [Vectors and quaternions](#vectors-and-quaternions).

## Scope

`zmath` is a dependency-free Zig 0.16 mathematics library. Public fixed-size operations do not allocate.

Public real formats are `f16`, `f32`, and `f64` for the current IEEE and scalar primitives. Public type and function names use lower-case `h`, `f`, and `d` suffixes. `q`, `f128`, and bfloat16 are outside this contract and have no public API.

The API uses one structural generic, `Mat`. Vectors, quaternions, geometry, factorizations, and geometric algebras use named concrete types.

## Planned module paths

The module path is canonical. Root re-exports are optional convenience names after an implementation becomes public.

```text
zmath.ieee
zmath.scalar
zmath.vector
zmath.complex
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
  vector.zig
  complex.zig
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

## Vector API

`zmath.vector` provides SIMD-style `@Vector` overloads for the elementwise
float-to-float operations. It is a distinct layer from the planned
linear-algebra `Vec*` value types described in [Vectors and
quaternions](#vectors-and-quaternions); those remain unimplemented.

Every vector operation applies the identical accepted scalar kernel to each
lane. Lane `i` of a vector result is therefore bit-for-bit equal to the scalar
result for lane `i` of the input. This preserves the `max_ulp = 0` contract and
the immutable corpus: the vector path is a per-lane application of the same
scalar code, never a re-derived polynomial, and it requires no new acceptance
evidence.

Each elementwise operation is available in three forms:

- As a top-level generic overload: `zmath.sin(x)` where `x: @Vector(N, T)` with
  `T` one of `f16`, `f32`, or `f64`. The element type selects the scalar kernel,
  and the vector width is inferred from the argument.
- As a generic vector entry point: `zmath.vector.sinv(x)`, which selects the
  scalar kernel from the lane type.
- As specialized entry points that pin the lane type: `zmath.vector.sinfv` and
  `zmath.vector.sindv` for `f32` and `f64`, plus `zmath.vector.sinhv` for
  `f16` where that name is not already used by another generic entry point.

The specialized name is the operation name, the format letter, and `v`
(for example `sinfv` for the sine of a `@Vector(N, f32)`).
Where a specialized `f16` name would collide with another operation's generic
entry point, the specialized name is omitted and the generic `opv` entry point
serves binary16 as well. This affects the bare trigonometry operations
`sin`, `cos`, `tan`, `asin`, `acos`, and `atan`, because their binary16
specialization would otherwise collide with the generic entry point of the
corresponding hyperbolic operation (`sinhv`, `coshv`, `tanhv`, `asinhhv`,
`acoshhv`, `atanhhv`). No operation loses coverage; binary16 callers use the
generic `opv` form.

```text
zmath.sin(@Vector(4, f32))          // top-level generic overload
zmath.vector.sinv(@Vector(4, f32))   // generic vector entry point
zmath.vector.sinfv(@Vector(4, f32))  // binary32 specialization
zmath.vector.sindv(@Vector(2, f64))  // binary64 specialization
zmath.vector.powv(a, b)              // elementwise binary
```

Vector operations cover the elementwise float-to-float operations only.
Operations that return a struct (`atan2`, `atan2pi`, `lgamma`, `frexp`,
`modf`, `remquo`, `canonicalize`, and the integer-power family), that return an
integer (`ilogb`, `llogb`, `lrint`, `llrint`), that take an integer argument
(`ldexp`, `scalbn`, `scalbln`, `pown`, `rootn`, `compoundn`), or that are
ternary or narrowing (`fma`, `fadd`, `fsub`, `fmul`, `fdiv`, `ffma`, `fsqrt`,
`fromfp`) are not vectorized.

## IEEE and scalar API

### Status terms

These are independent status terms:

- **Target** names a planned public symbol.
- **Implemented** means that the current `zmath` root or module exports the symbol.
- **Manifest-accepted** means that a `correct` or `bounded` manifest row accepts one operation-format pair.
- **Corpus-evidenced** means that immutable rows for one pair passed the stated oracle checks.
- **Whole-domain-proved** means that exhaustive evidence or a separate proof covers every input in the operation-format domain.

Implementation, manifest acceptance, corpus coverage, and proof scope are separate facts.
Export presence does not change manifest status.
Manifest status does not widen a corpus beyond its stored inputs.
Sampled corpus evidence does not establish whole-domain proof.

The implemented C23 inventory is the `Implemented-source snapshot` in `spec/core-math-2026-08-22.md`.
This document does not maintain a second exhaustive implemented-export list.

The current source exports the `zmath.ieee` and `zmath.scalar` families described in that contract, their shared result types, and the corresponding generic root operations.
`src/c23_public_surface_test.zig` validates those module and root signatures for `f16`, `f32`, and `f64`.

The source also exports the six C23 root narrowing operations (`fadd`, `fsub`, `fmul`, `fdiv`, `ffma`, and `fsqrt`) and the Annex H format-pair forms from `zmath.scalar`.
These source exports do not imply manifest acceptance.
The versioned corpus manifest is the only acceptance authority.

`ulpDistance` is a zmath extension. It is not a C23 scalar operation.

`nextAfterh`, `nextAfterf`, and `nextAfterd` choose the first NaN operand.
If it is signaling, they set its quiet bit and preserve its sign and payload.
They return a quiet NaN unchanged.

`totalOrderh`, `totalOrderf`, and `totalOrderd` return `true` when the left operand precedes or equals the right operand in IEEE total order.


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
`zmath.cos(x)` selects `zmath.scalar.cosh`, `cosf`, or `cosd`.
`zmath.cosh(x)` selects `zmath.scalar.coshh`, `coshf`, or `coshd`.
`zmath.tan(x)` selects `tanh`, `tanf`, or `tand`.
`zmath.tanh(x)` selects `tanhh`, `tanhf`, or `tanhd`.

The scalar name is the operation name followed by the type suffix. The suffix rule produces `sinhh`, `coshh`, `tanhh`, `asinhh`, `acoshh`, and `atanhh` for binary16 hyperbolic operations.

All floating arguments to a normal binary or ternary generic must have the same type. Dispatch does not widen, narrow, or coerce a floating argument.

The first floating argument selects the family when an operation also takes an integer or an output pointer. Integer-result operations keep the suffix of their floating input.

`zmath.abs` is the C23 `fabs` equivalent. Its scalar symbols remain `absh`, `absf`, and `absd`.
The API does not add `fabs`, `fabsh`, `fabsf`, or `fabsd` aliases.

The target has no `q`, `f128`, bfloat16, or bfloat suffix family. These types do not participate in generic dispatch.

### Pure-value and rounding profile

The scalar and IEEE APIs use a pure-value profile.
A call communicates only through its return value, including an explicit result struct when required.
No result contract uses `errno`, floating-point exception flags, or hidden mutable status.
`zmath` does not read or write `errno`.
It does not inspect, clear, report, or restore floating-point exception flags.
Incidental processor flags are outside this contract.

Callers must keep IEEE round-to-nearest, ties-to-even active for each scalar call.
Scalar functions do not take a general `RoundMode` argument.
An operation-defined direction parameter, such as the `fromfp` direction, does not change the process environment.
`zmath` does not read, validate, set, or restore the process rounding mode.
Behavior with another active rounding mode is unsupported.

This caller precondition also applies to narrowing operations.
They compute the exact mathematical operation and round once to the declared destination type.

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

The `lrint`, `llrint`, `lround`, and `llround` families return this shared result type:

```zig
pub const IntegerResult = struct { value: i64, ok: bool };
```

A valid conversion returns `.{ .value = result, .ok = true }`.
A NaN, either infinity, or a rounded value outside the `i64` range returns `.{ .value = 0, .ok = false }`.
Invalid input does not trap and does not use a numeric sentinel.

The `fromfp`, `ufromfp`, `fromfpx`, and `ufromfpx` families return an integral value in the input floating type `T`.
Each family has the target signature `(value: T, round: FromFpRound, width: u32) T`.
The `u32` width maps the C `unsigned int` argument.

```zig
pub const FromFpRound = enum {
    nearest_even,
    nearest_from_zero,
    toward_zero,
    toward_positive,
    toward_negative,
};
```

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
`atan2` and `atan2pi` return `Atan2Result(T) { value: T, resolved: bool }`.
The root generics and all suffix variants use this result shape.
The public symbols are wrappers, and direct CORE-MATH kernel symbols remain private.
The wrappers own all NaN, axis, zero, and infinity results.
The `value` field is an operation result only when `resolved=true`.
Every ordinary result and every documented special result sets `resolved=true`.
This includes a quiet NaN selected from a NaN input.
If a finite kernel cannot prove one nearest-even result, the call sets `resolved=false`.
The `value` field is then the format's positive canonical quiet NaN.
This false status reports an unresolved computation, not a mathematical NaN or a domain error.
No valid input may trap or terminate the process.
The wrappers do not guess a result bit when a kernel cannot decide.
They do not call host libm or use `fenv` or `errno`.


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

`compoundn`, `pown`, and `rootn` take an `i64` exponent or degree.
They return `PowerResult(T) { value: T, ok: bool }`.
The root generics and all suffix variants use this result shape.
The accepted integer range is `-64` through `+64`, inclusive.
When `n` is below `-64` or above `+64`, the range check precedes all input-value rules.
The `value` field is then the format's positive canonical quiet NaN, and `ok=false`.
Every call with `-64 <= n` and `n <= 64` sets `ok=true`.
The `value` field then holds the specified result, including an in-range domain or NaN result.
Thus, `ok` reports only whether `n` is in the accepted range.

The C23 `compoundn(x, n)` operation computes `(1 + x)^n`.
A CORE-MATH `compound` source with a floating exponent is not a direct implementation of this operation.
The power families use the ordered IEC 60559 and glibc value precedence below.
An earlier rule wins.
For the three `PowerResult` families, each result below describes `value` when `ok=true`.
This project adopts the return values, but it does not adopt `fenv` or `errno` effects.

For `pow`, signaling-NaN quieting precedes an identity that would otherwise hide that signaling NaN.
Quiet NaNs do not override `pow(x, +/-0)=+1` or `pow(+1, y)=+1`.
`pow(-1, +/-infinity)` is `+1`.
After these identities, the first remaining NaN in `(x, y)` parameter order wins.
For an infinite exponent, `abs(x)` relative to one selects `+0` or `+infinity`.
For a zero or infinite base and a finite nonzero exponent, the exponent sign selects zero or infinity.
Such a result is negative only when the base is negative and the exponent is an odd integer.
A finite negative base with a finite noninteger exponent returns the positive canonical quiet NaN.

For `pown`, a signaling NaN is quieted before the `n=0` identity.
Then `n=0` returns `+1`, including for a quiet NaN base.
A remaining NaN is quieted.
Zero and infinity results use the exponent sign, and a negative sign occurs only for an odd `n`.

For `powr`, four domain cases precede ordinary NaN selection.
They are a negative non-NaN base, `0^0`, `1^infinity`, and `infinity^0`.
Each domain case returns the positive canonical quiet NaN.
After these cases, the first NaN in `(x, y)` parameter order wins.

For `compoundn`, a signaling NaN is quieted first.
A non-NaN `x < -1` then returns the positive canonical quiet NaN.
Next, `n=0` returns `+1`, including for a quiet NaN.
A remaining NaN is quieted.
For `x=-1`, a positive `n` returns `+0`, and a negative `n` returns `+infinity`.
For `x=+infinity`, a positive `n` returns `+infinity`, and a negative `n` returns `+0`.

For `rootn`, `n=0` returns the positive canonical quiet NaN before NaN selection.
A remaining NaN is quieted.
A negative nonzero value with an even `n` returns the positive canonical quiet NaN.
Zero and infinity results preserve a negative sign only for odd `n`.
The sign of `n` selects the zero or infinity magnitude.

Each finite direct result from a power, root, or pi-scaled family must be correctly rounded in `nearest_even`.
For a result-container function, this rule applies to a finite `float` only when its status field is true.
For `pown`, `compoundn`, and `rootn`, the rule applies only when `-64 <= n` and `n <= 64`.
A false status supplies no mathematical result and makes no accuracy claim for its sentinel `float`.
This requirement does not claim that a current implementation satisfies it.
A composition of separately rounded operations is not final evidence.
This includes post-scaling by an approximation of `1/pi`, reciprocal-exponent power, and target-format power-by-squaring.

The implemented bounded power kernels use fixed, format-specific integer storage.
They use 32 `u64` limbs for `f16`, 256 limbs for `f32`, and 2,048 limbs for `f64`.
`pown` and `compoundn` keep finite dyadic arithmetic exact and round once.
`rootn` uses exact candidate and midpoint comparisons instead of a rounded reciprocal-power composition.
These storage bounds prevent allocation. They do not specify a universal execution-time bound.

#### Error and gamma

| Root generic | Scalar implementation family |
| --- | --- |
| `erf` | `erf{h,f,d}` |
| `erfc` | `erfc{h,f,d}` |
| `lgamma` | `lgamma{h,f,d}` |
| `tgamma` | `tgamma{h,f,d}` |

`lgamma` returns `LgammaResult(T) { value: T, sign: i2 }`.
The root generic and all suffix variants use this result shape.
For a finite non-pole input, `value` is the rounded value of `log(abs(Gamma(x)))`.
The `sign` field is the sign of `Gamma(x)` and is either `-1` or `+1`.
No public result uses `sign=0`.

At each finite negative-integer pole and at either zero, `value` is `+infinity`.
The `sign` field is `-1` at each finite negative-integer pole and at `-0`.
It is `+1` at `+0`, either infinity, and every NaN.
For either infinity, `value` is `+infinity`.
For a NaN, `value` follows the common same-width NaN policy.

`lgamma` follows the pure-value and caller-RNE profile.
It returns only `value` and `sign` and exposes no `errno` or floating-point exception status.
It does not inspect, clear, report, set, or restore floating-point exception flags.
Incidental processor flags remain outside this contract.
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

`frexp` returns `FrexpResult(T) { fraction: T, exponent: i32 }`.
`modf` returns `ModfResult(T) { fractional: T, integral: T }`.
`remquo` returns `RemquoResult(T) { value: T, quotient: i32 }`.
These structs replace C output pointers.

`remquo` uses the nearest integer quotient, with ties to even.
Let `n` be that quotient.
Its quotient field is `sign(x/y) * (abs(n) mod 8)`.
Thus, the contract exposes exactly three low quotient bits.
A selected NaN, an infinite dividend, or a zero divisor sets `quotient=0`.
The NaN result follows the common NaN policy.
An infinite dividend or zero divisor without a NaN returns the positive canonical quiet NaN.
A zero dividend with a nonzero divisor, or a finite dividend with an infinite divisor, returns the dividend and `quotient=0`.

`frexp` returns `exponent=0` for zero, infinity, and NaN.
Its fraction preserves a zero or infinity bit pattern.
For a NaN, its fraction preserves the sign and payload and sets the quiet bit.

`modf` returns the input zero in both fields.
For a finite value with magnitude below one, the integral field is a zero with the input sign.
For an integral finite value, the fractional field is a zero with the input sign.
For infinity, the fractional field is a zero with the input sign, and the integral field is the input infinity.
For a NaN, both fields contain the same sign-preserving, payload-preserving quiet NaN.

`ilogb` returns `minInt(i32)` for either zero and `maxInt(i32)` for infinity or NaN.
`llogb` uses the corresponding `i64` limits.
`logb` returns `-infinity` for either zero, `+infinity` for either infinity, and a quieted input NaN for NaN.

The `ldexp`, `scalbn`, and `scalbln` families preserve a NaN sign and payload.
A quiet NaN is bit-identical.
A signaling NaN sets only the format quiet bit.
These rules freeze the current deterministic remquo and decomposition outputs.
They are not delegated to host library behavior.

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

The binary16 `fmah` kernel multiplies and adds exact integer dyadics.
It performs one final round-to-nearest, ties-to-even packing step.
It does not use a widened host FMA followed by a second target-format rounding.

The legacy `fmax` and `fmin` families use a deterministic local policy.
This policy is a zmath contract, not a portable C guarantee.
With exactly one NaN operand, they return the numeric operand.
With two NaN operands, they return the first operand after preserving its sign and payload and setting its quiet bit.
For an opposite signed-zero tie, `fmax` returns `+0` and `fmin` returns `-0`.

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

The target also includes the C23 Annex H `_FloatM` forms for the supported public formats:

| Destination | Source | Add | Subtract | Multiply | Divide | Fused multiply-add | Square root |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `f16` | `f32` | `f16addf32(x: f32, y: f32) f16` | `f16subf32(x: f32, y: f32) f16` | `f16mulf32(x: f32, y: f32) f16` | `f16divf32(x: f32, y: f32) f16` | `f16fmaf32(x: f32, y: f32, z: f32) f16` | `f16sqrtf32(x: f32) f16` |
| `f16` | `f64` | `f16addf64(x: f64, y: f64) f16` | `f16subf64(x: f64, y: f64) f16` | `f16mulf64(x: f64, y: f64) f16` | `f16divf64(x: f64, y: f64) f16` | `f16fmaf64(x: f64, y: f64, z: f64) f16` | `f16sqrtf64(x: f64) f16` |
| `f32` | `f64` | `f32addf64(x: f64, y: f64) f32` | `f32subf64(x: f64, y: f64) f32` | `f32mulf64(x: f64, y: f64) f32` | `f32divf64(x: f64, y: f64) f32` | `f32fmaf64(x: f64, y: f64, z: f64) f32` | `f32sqrtf64(x: f64) f32` |

These format-pair declarations are Annex H forms. They are not Annex F declarations.

These names are target scalar implementation names. They are not aliases for a same-type operation.
The target has no `f128`, `q`, bfloat16, decimal, or long-double source or destination form.

The legacy `fadd`, `fsub`, `fmul`, `fdiv`, `ffma`, and `fsqrt` forms are root-level C23 names with the signatures above.
The format-pair forms are available from `zmath.scalar` and use the exact names in the table.

Each narrowing operation computes as if to infinite precision and rounds once to the declared destination type.
The implementation must form the exact intermediate in portable software and round only at the destination.
Target-dependent extended precision, contraction, or a direct narrowing shortcut is not final correctness evidence.

For every cross-width narrowing operation, the first NaN in parameter order selects the result.
The result preserves its sign and the high-order diagnostic payload bits that fit the destination.
It discards lower payload bits and then sets the destination quiet bit.


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
```

For input type `T`, `canonicalize{h,f,d}(value: T)` returns `CanonicalizeResult(T) { value: T, ok: bool }`.
This value-returning local adaptation replaces the C output pointer and integer success result.
Each supported `h`, `f`, or `d` input returns `ok=true`.
A non-NaN `value` is bit-identical to the input, including zero and infinity.
A quiet NaN is bit-identical to the input.
A signaling NaN preserves its sign and diagnostic payload and sets only the format quiet bit.
The supported fixed formats never return `ok=false`.
No supported fixed-format call can produce a failure result, so no fixed-format failure value exists.
That status is reserved for a future representation which the API cannot canonicalize.

The `nan` family has the target signature `nan{h,f,d}(tag: []const u8) T`.
The tag grammar is one or more ASCII hexadecimal digits: `[0-9A-Fa-f]+`.
The parser accepts no sign, `0x` prefix, separator, whitespace, empty tag, or non-ASCII byte.
It interprets the tag as an unsigned diagnostic-payload integer.
The value must fit 9 bits for `h`, 22 bits for `f`, or 51 bits for `d`.
A valid tag produces a positive quiet NaN with that diagnostic payload.
An invalid, empty, or out-of-range tag produces the positive canonical quiet NaN.

The `nextToward` family has the target signature `nextToward{h,f,d}(value: T, toward: f64) T`.
Its comparator is `f64` because `f64` is the widest public format.

C leaves this cross-format NaN mapping implementation-defined.
When `toward` is an `f64` NaN, `nextTowardh` and `nextTowardf` preserve its sign.
They retain the high-order diagnostic payload bits that fit the target, discard lower bits, and set the target quiet bit.
For `f64` to `f16`, source payload bits 50 through 42 become destination payload bits 8 through 0.
For `f64` to `f32`, source payload bits 50 through 29 become destination payload bits 21 through 0.
This mapping is a zmath contract, not a portable C guarantee.

Classification predicates, signaling tests, ordered comparisons, and bit conversion also remain in `zmath.ieee`.

### Optional IEC 60559 extension family

C23 Annex F defines these optional IEC 60559 operations. They are target names, but they are not baseline C23 7.12 scalar names:

| Root generic | `zmath.ieee` implementation family |
| --- | --- |
| `totalOrder` | `totalOrder{h,f,d}` |
| `totalOrderMag` | `totalOrderMag{h,f,d}` |
| `getPayload` | `getPayload{h,f,d}` |
| `setPayload` | `setPayload{h,f,d}` |
| `setPayloadSig` | `setPayloadSig{h,f,d}` |

`setPayload` and `setPayloadSig` return `PayloadResult(T) { value: T, ok: bool }`.
This struct replaces the C result pointer with a value and success flag.
These optional operations remain a by-value API fixed to the `h`, `f`, and `d` formats.
They do not expose C output pointers or a dynamic format parameter.

C leaves the accepted NaN payload range implementation-defined.
This project accepts a fixed-format floating argument only when it represents an allowed numeric integer payload.
For `setPayload`, both signed zeros represent `q=0`.
The exact success ranges and result bits are:

| Format | Accepted `q` | `setPayload` result bits | Accepted `s` | `setPayloadSig` result bits |
| --- | --- | --- | --- | --- |
| `h` | `0` through `511` | `0x7e00 OR q` | `1` through `511` | `0x7c00 OR s` |
| `f` | `0` through `4,194,303` | `0x7fc00000 OR q` | `1` through `4,194,303` | `0x7f800000 OR s` |
| `d` | `0` through `2,251,799,813,685,247` | `0x7ff8000000000000 OR q` | `1` through `2,251,799,813,685,247` | `0x7ff0000000000000 OR s` |

Each successful `setPayload` call returns the exact positive quiet NaN in the table with `ok=true`.
Each successful `setPayloadSig` call returns the exact positive signaling NaN in the table with `ok=true`.
Zero, including `-0`, is invalid for `setPayloadSig` because it would encode infinity.
A negative nonzero, fractional, subnormal, infinite, NaN, or out-of-range argument is invalid for either setter.
Each failure returns the exact `PayloadResult(T) { value=+0, ok=false }`.
This policy is a zmath contract, not a portable C guarantee.

### Exact special-value policy

Each manifest-accepted operation-format pair must define exact output bits for zeros, subnormals, infinities, and NaNs.

The format encodings are:

| Format | `+0` | `-0` | `+1` | `+infinity` | `-infinity` | Canonical quiet NaN |
| --- | --- | --- | --- | --- | --- | --- |
| `h` | `0000` | `8000` | `3c00` | `7c00` | `fc00` | `7e00` |
| `f` | `00000000` | `80000000` | `3f800000` | `7f800000` | `ff800000` | `7fc00000` |
| `d` | `0000000000000000` | `8000000000000000` | `3ff0000000000000` | `7ff0000000000000` | `fff0000000000000` | `7ff8000000000000` |

When an operation returns a same-width input NaN, it preserves the sign and diagnostic payload.
It then sets the format quiet bit.
The quiet-bit masks are `0200`, `00400000`, and `0008000000000000` for `h`, `f`, and `d`.
A quiet NaN remains unchanged.
A signaling NaN becomes a quiet NaN.

When an operation returns a narrower NaN, let `S` and `D` be the source and destination diagnostic-payload widths.
The destination payload is `source_payload >> (S - D)`.
The operation preserves the sign and then sets the destination quiet bit.

Unless an operation has a stricter rule, the first NaN argument in parameter order selects the result.
An invalid non-NaN input returns the positive canonical quiet NaN for its format.
An operation-specific documented rule takes precedence.
This includes the local `fmax` and `fmin` policy and the C23 number-preferring variants.

Signed-zero and infinity results follow each operation's documented policy.
That policy can be local, C23, or IEC 60559.
Implementations must not replace those results with a host-library choice.

The manifest must record every operation-specific exception before acceptance.
Target status or source availability cannot substitute for this evidence.

### Typed scalar corpus manifest

Before a manifest-gated scalar pair becomes accepted, the repository updates a versioned contract at `spec/core-math-<date>.md`.

The schema-3 manifest is `zmath-testing/test/corpus/manifest.toml`.
It contains exactly 315 sorted typed function rows for `f16`, `f32`, and `f64`.
It pins the CORE-MATH revision, 512-bit MPFR/Arb settings for evaluator paths that use them, `nearest_even`, and the special-value policy.
The checked-in manifest is the only acceptance authority.
A successful candidate campaign does not change a row's manifest status.

Each row defines the operation, provenance, status, arity, typed inputs, typed result, oracle, special-value policy, and maximum ULP bound.
The root `operation_spelling` value is `"c23_hfd_lowercase_v1"`.
For a standard H/F/D family, the typed descriptors select the suffix.
Dedicated encoded names, including narrowing operations, remain unchanged.

The provenance vocabulary is exactly `["zmath", "core_math"]`.
Each row uses `["zmath"]`, `["core_math"]`, or `["zmath", "core_math"]`.
Scalar results use `float:<format>`, `integer:<type>`, or `bool`.
Structure results use `struct:<kind>:<format>`.
Status-container descriptors include:

- `struct:atan2:<format>` for `Atan2Result(T) { value: T, resolved: bool }`.
- `struct:power:<format>` for `PowerResult(T) { value: T, ok: bool }`.
- `struct:integer:<format>` for `IntegerResult { value: i64, ok: bool }`.

These descriptors record the public return shape.
They do not grant implementation, corpus evidence, or manifest acceptance.

The status values are `correct`, `bounded`, and `unsupported`.
An `unsupported` row uses `oracle = "none"` and `max_ulp = 0`.
An accepted row uses the evaluator class assigned by the authoritative C operation table.
A `correct` row uses `max_ulp = 0`.
A `bounded` row uses an approved non-negative `max_ulp` value.
Schema 3 accepts `none`, `exact_bits`, `algebraic`, and `mpfr_arb`.

Status and oracle snapshots are not part of numerical descriptor identity.
This permits review of candidate artifacts before a manifest promotion.
The loader still checks that each snapshot has a valid status-oracle relation.
It hashes the profile, numerical descriptors, special policy, ULP bound, and authoritative C signature.

The current checked-in manifest accepts all 315 descriptors as `correct` with `max_ulp = 0`, backed by full publication and replay of 63,160,084 records and benchmarks for every descriptor.
Do not infer acceptance from source exports, evaluator availability, or a filtered campaign.

#### Oracle evidence classes

The authoritative C operation table selects one evaluator for each of the 315 descriptors.
The table assigns 150 descriptors to `exact_bits`, 126 to `mpfr_arb`, and 39 to `algebraic`.
Operation spelling and manifest status do not select an evaluator.

- `exact_bits` uses direct IEEE bit and integer models. It preserves bit-sensitive NaN and status behavior without numerical approximation.
- `algebraic` uses exact GMP rational arithmetic for finite dyadic operations. Irrational algebraic roots use 512-bit MPFR candidates and 512-bit Arb enclosures.
- `mpfr_arb` uses a 512-bit MPFR candidate and a 512-bit Arb proof, with operation-specific certified reductions where required.

Exact-bit and algebraic boundary checks use exact integer or rational arithmetic.
The MPFR/Arb oracle uses fixed 1,024-bit MPFR boundary containers where their capacity is sufficient for exact binary values.
This boundary storage is not adaptive numerical evaluation.
The numerical MPFR and Arb precision remains fixed at 512 bits.

The generator fails closed when an evaluator reports unsupported, error, uncertainty, a descriptor mismatch, or an unresolved rounding cell.
It does not increase numerical precision or emit a best-effort record.

Fixed MPFR512 evaluation followed by target RNE can double-round
`atan2((2^52-1)*2^-1074, 2)`.
The generator may correct the initial candidate by one adjacent lower-magnitude value only under the following proof:

1. Exact dyadic arithmetic proves that `t = abs(y/x)` equals the midpoint above the corrected candidate.
2. The proof establishes `0 < t < 1`.
3. Arb512 proves `t^3` is less than the gap from `t` to the next lower rounding midpoint.
4. The analytic bound `t - t^3 < atan(t) < t` places the exact result inside the corrected candidate's rounding cell.

The correction preserves the result sign.
If any condition is unresolved, the generator rejects the record.

For exponentials and powers, Arb512 can prove positive overflow in the log domain.
This proof does not construct the enormous result.
For negative noninteger `lgamma` inputs, the Arb proof uses the reflection formula and evaluates `lgamma` at a positive argument.

#### Typed corpus publication

The generator uses oracle ABI v2.
Every request declares the operation bytes, arity, typed input fields, result layout, structure kind, and typed result fields.
The generator compares the zmath runtime result with independent evaluator evidence before it writes a record.

Published records use the length-framed binary `ZCDR` codec version 1.
A record stores operation, ordered provenance, typed arguments, result descriptor, typed result fields, and `max_ulp`.
Floating fields store the exact format, bits, and IEEE class.
Other fields store their integer width, Boolean value, tag bytes, or rounding enumeration.

A full publication covers all 315 descriptors.
A filtered publication covers the exact `--operation` and `--format` selection.
Metadata records the campaign class, seed, cardinality, selection rule, profile, and artifact digest.
The 61 exhaustive unary H campaigns each enumerate all 65,536 input bit patterns.
Only those campaigns claim exhaustive coverage of their one-argument binary16 domain.
Selected Cartesian, explicit special-tuple, selected unary, bounded-integer, and exact-status campaigns cover only their declared tuples.

The publication loader rejects `.incomplete`, requires `COMPLETE`, and verifies the metadata SHA-256.
It also verifies the full or filtered scope, target profile, descriptor digest, record count, record order, and artifact SHA-256.
Unavailable per-campaign source metadata uses the exact `source-metadata-unavailable` digest sentinel and `none` source fields.

Generated raw records can be multi-gigabyte artifacts.
Keep `generated/records/`, `generated/probes/`, `.incomplete`, and `*.tmp` local and ignored.
Commit publication metadata, `COMPLETE`, `replay.txt`, and `benchmark.txt`.
Regenerate the raw records with `generate-corpus` in a new checkout.

The accepted corpus profile is target `x86_64-linux-gnu.2.38`, CPU `haswell`, and mode `ReleaseSafe`.
No other target, CPU, or optimization mode is manifest-qualified.
Corpus evidence applies only to stored tuples unless an exhaustive campaign or separate proof covers the whole operation-format domain.
Benchmark output measures performance only and does not change manifest acceptance.
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
`~/projects/zmath-testing` owns integration tests, oracle generation, typed corpus replay, and benchmarks.

The scalar tools use only target `x86_64-linux-gnu.2.38`, CPU `haswell`, and mode `ReleaseSafe` for manifest-qualified evidence.
These build steps install the typed tools:

```text
zig build oracle-generator -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
zig build corpus-replay -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
zig build bench-scalar -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
zig build bench-helpers -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
```

The installed executables are:

```text
zig-out/bin/zmath-testing-corpus-producer
zig-out/bin/zmath-corpus-replay
zig-out/bin/zmath-scalar-matrix
zig-out/bin/zmath-public-helpers
```

`zmath-public-helpers` measures the public evaluators that no corpus descriptor
covers: `classify`, `ulpDistance`, all four `scaleB` rounding modes, all four
`roundIntegral` rounding modes, and integer `add`. It takes no corpus root and
no operation filter, and accepts only `--repetitions=1..100`. The typed
benchmark accepts `--repetitions=1..1024` and `--with-zig`, which times
same-input Zig intrinsic references beside zmath.

The corpus producer, corpus replay, and scalar-matrix executables accept an
optional first positional corpus root. The default root is `test/corpus`.
These three executables also accept `--operation=<manifest-op>` and
`--format=h|f|d`. The filters can be used independently or together.

The standard full-publication commands are:

```text
zig build generate-corpus -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
zig build replay-corpus -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
zig build run-bench-scalar -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe
```

The build steps already pass `test/corpus` as the positional root.
Append filters after `--`, for example:

```text
zig build generate-corpus -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe -- --operation=atan2 --format=d
zig build replay-corpus -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe -- --operation=atan2 --format=d
zig build run-bench-scalar -Dtarget=x86_64-linux-gnu.2.38 -Dcpu=haswell -Doptimize=ReleaseSafe -- --operation=atan2 --format=d
```

The generator owns independent oracle evaluation and proof.
Replay reads stored typed outputs and compares them with current zmath results.
The replay executable does not link MPFR, Arb, GMP, or FLINT.
The scalar benchmark validates publication metadata and artifact digests before measurement.
It reports measurements only and does not grant manifest acceptance or whole-domain proof.

A numerical oracle record names the exact zmath symbol, typed inputs, typed outputs, manifest descriptor, campaign, and acceptance rule.
An external differential test names its baseline library, version, build options, scalar type, dimensions, and operation.
