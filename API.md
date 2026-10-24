# zmath API Contract

The public root exports `ieee`, `scalar`, `vector`, `complex`, `quaternion`,
`matrix`, `linalg`, `geometry`, and `algebra`, their shared result types, and
the C23 source surface in [the scalar contract](spec/core-math-2026-08-22.md).

## Scope

`zmath` is a dependency-free Zig 0.17.0 mathematics library. Public fixed-size operations do not allocate.

Public real formats are `f16`, `f32`, and `f64` for the current IEEE and scalar primitives. Public type and function names use lower-case `h`, `f`, and `d` suffixes. `q`, `f128`, and bfloat16 are outside this contract and have no public API.

The API uses one structural generic, `Mat`. Complex values, vectors, quaternions, geometry, factorizations, and geometric algebras use named concrete types.

## Module paths

Root re-exports provide convenience names for scalar and IEEE operations.

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
  vector_value.zig
  complex.zig
  quaternion.zig
  matrix.zig
  linalg.zig
  geometry.zig
  algebra.zig
```

`zmath.vector` is implemented across `src/vector*.zig`: `src/vector.zig` plus the five modules it re-exports, as described in [Vector API](#vector-api).

The dependency direction is:

```text
ieee -> scalar -> complex
             +-> vector -> quaternion
complex + vector + quaternion -> matrix -> linalg
vector + quaternion + matrix -> geometry
complex + vector + quaternion + matrix -> algebra
```

Geometry and algebra do not depend on each other.

## Dependency and build setup

See [README.md](README.md#usage) for dependency setup and build commands.
The module carries its private C kernels and link settings.
The C headers are not a standalone C API.

Checked geometry returns `NonFiniteInput` for invalid supplied coordinates
and `NonFiniteResult` for arithmetic range failure. Fixed-size values need no cleanup.
Other numerical APIs use IEEE values, status fields, or `null` as documented below.


## Vector API

`zmath.vector` provides `@Vector` overloads for the elementwise scalar
operations: 417 public entry points across six modules. It is a distinct
layer from the linear-algebra `Vec*` value types described in
[Vectors and quaternions](#vectors-and-quaternions), which the same module
also owns and which are implemented in `src/vector_value.zig`. The two layers
share this module path and nothing else.

Every vector operation applies the scalar kernel to each lane. Lane `i` of
a vector result has the scalar result for lane `i` of the input.

`src/vector_surface_test.zig` guards this surface, and `src/vector_test.zig` carries the per-lane bit-identity tests, including an exhaustive pass over all 65,536 binary16 bit patterns.

The layer is portable. `@Vector` is lowered by the compiler for the selected
target. There is no runtime CPU detection, no hand-written per-ISA kernel, and
no per-format selection beyond the lane type itself. What the layer provides
is a uniform lane-parallel API, not a numerical speedup over a scalar loop.

### Modules

| Module | Source | Covers |
| --- | --- | --- |
| `zmath.vector` | `src/vector.zig` | The 214 float-to-float entry points, where a lane maps to a lane of the same float type |
| `zmath.vector.struct_family` | `src/vector_struct.zig` | The 13 operations whose scalar kernel returns a result struct |
| `zmath.vector.intarg` | `src/vector_intarg.zig` | The 11 operations that take a non-float argument |
| `zmath.vector.narrow` | `src/vector_narrow.zig` | The ternary `fma` and the 6 C23 narrowing operations |
| `zmath.vector.ieee_family` | `src/vector_ieee.zig` | The remaining IEEE and miscellaneous operations, including the integer-returning and boolean-returning ones |
| `zmath.vector.core` | `src/vector_core.zig` | The shared per-lane mapping core behind `intarg`, `narrow`, and `ieee_family` |

The four families and `core` are re-exported from `zmath.vector`, so the whole `@Vector`
surface is reachable as `zmath.vector.*`. `core` is public so that a caller
can write a further lane-parallel entry point from the same helpers, with the
same per-lane bit-identity. `struct_family` does not use it, because every
`core` helper returns a `@Vector`; it carries its own per-lane helpers instead.

The six `Vec*` value types are re-exported from `zmath.vector` as well, so the module path serves
both layers. The table above is the `@Vector` surface only; `src/vector_value.zig` is not a
lane-parallel module and its operations carry no per-lane bit-identity argument.

The two counts are disjoint and are taken the same way. 417 is the number of
`pub fn` declarations across the six `@Vector` modules in the table above. 158 is
the number of `pub fn` declarations the six `Vec*` structs publish in
`src/vector_value.zig`, counting their instance methods, their constructors, and
their static factories. The two numbers share a module path and nothing else, and
neither is a subset of the other.

The root of `zmath` publishes 103 runtime callables: 97 `pub fn` declarations and
6 `pub const` function aliases, which are `fmaximum_mag`, `fminimum_mag`,
`fmaximum_num`, `fminimum_num`, `fmaximum_mag_num`, and `fminimum_mag_num`.
96 of the 103 accept a `@Vector` and route to this layer, which is 90 of the 97
`pub fn` declarations plus all 6 aliases. The seven that do not are:

- `nan`, whose type argument is a `comptime` type and therefore has no lane to
  map over. The vector spelling of a NaN is `@splat` at the call site.
- `fadd`, `fsub`, `fmul`, `fdiv`, `ffma`, and `fsqrt`, whose root signatures take
  concrete `f64` arguments and return `f32`. A concrete `f64` parameter cannot
  accept a `@Vector`, so those six are scalar-only at the root. Their vector
  forms exist and belong to this layer, as `zmath.vector.narrow.faddv`, `fsubv`,
  `fmulv`, `fdivv`, `ffmav`, and `fsqrtv`, which take a `@Vector(N, f64)` and
  return `@Vector(N, f32)`.

A fixed `f64` parameter by itself does not make a root operation scalar-only.
`nextToward(value, toward: f64)` still takes a `@Vector` value and broadcasts its
`f64` comparator, and it is one of the 96.

### Naming

Each operation is available in these forms:

- As a top-level generic overload: `zmath.sin(x)` where `x: @Vector(N, T)` with
  `T` one of `f16`, `f32`, or `f64`. The element type selects the scalar kernel,
  and the vector width is inferred from the argument.
- As a generic vector entry point: `zmath.vector.sinv(x)`, which selects the
  scalar kernel from the lane type.
- As specialized entry points that pin the lane type: `zmath.vector.sinfv` and
  `zmath.vector.sindv` for `f32` and `f64`, plus `zmath.vector.exphv` for `f16`,
  which is the spelling wherever the `f16` name is not already used by another
  generic entry point.

The specialized name is the operation name, the format letter, and `v`
(for example `sinfv` for the sine of a `@Vector(N, f32)`).
Where a specialized `f16` name would collide with another operation's generic
entry point, the specialized name is omitted and the generic `opv` entry point
serves binary16 as well. This affects the bare trigonometry operations
`sin`, `cos`, `tan`, `asin`, `acos`, and `atan`, because their binary16
specialization would otherwise collide with the generic entry point of the
corresponding hyperbolic operation (`sinhv`, `coshv`, `tanhv`, `asinhv`,
`acoshv`, `atanhv`). No operation loses coverage; binary16 callers use the
generic `opv` form.

`zmath.vector.sinhv` is therefore the binary16 HYPERBOLIC sine applied to
`@Vector(N, f16)` lanes, not the binary16 SINE. The bare trigonometric
generics `sinv`, `cosv`, `tanv`, `asinv`, `acosv`, and `atanv` keep the f16
binding of the root dispatch, where `sin` maps to `zmath.scalar.sinh` (the
binary16 sine), `cos` to `zmath.scalar.cosh`, and `tan` to `zmath.scalar.tanh`,
while `sinh` maps to `zmath.scalar.sinhh`. See
[Generic dispatch](#generic-dispatch).

```text
zmath.sin(@Vector(4, f32))           // top-level generic overload
zmath.vector.sinv(@Vector(4, f32))    // generic vector entry point
zmath.vector.sinfv(@Vector(4, f32))   // binary32 specialization
zmath.vector.sindv(@Vector(2, f64))   // binary64 specialization
zmath.vector.powv(a, b)               // elementwise binary
zmath.vector.sinhv(@Vector(4, f16))   // binary16 HYPERBOLIC sine
```

### Result shapes

Most results are a real `@Vector(N, R)`:

| Operations | Result |
| --- | --- |
| Float-to-float, and the narrowing operations | `@Vector(N, T)`, and `@Vector(N, f32)` for `@Vector(N, f64)` narrowing input |
| `ilogb` and `llogb` | `@Vector(N, i32)` and `@Vector(N, i64)` |
| `totalOrder` and `totalOrderMag` | `@Vector(N, bool)` |

The 13 result-struct operations cannot be a `@Vector`. Zig forbids a
struct as a vector element type and reports "expected integer, float, bool, or
pointer for the vector element type", so `atan2`, `atan2pi`, `lgamma`,
`frexp`, `modf`, `remquo`, `canonicalize`, `setPayload`, `setPayloadSig`,
`lrint`, `llrint`, `lround`, and `llround` return a lane-ordered array
`[N]Result` instead, and so do `pown`, `compoundn`, and `rootn`, which also
return a result struct. Element `i` of that array is the value the scalar
kernel itself returned for lane `i`, moved whole and in lane order, so the
per-lane bit-identity guarantee is unchanged. This is a language limitation,
not an oversight.

The array shape buys no throughput. The kernel is applied once per lane,
exactly as a scalar loop would, so the array form is a uniform lane-parallel
API and not SIMD acceleration.

Result structs in `scalar_c23_types.zig` are bare `struct` values with auto
layout. They are never `extern` and never `packed`, so `@bitCast` and
`@offsetOf` do not apply to them. Every field is assigned explicitly and
nothing is reinterpreted.

### Integer-argument operations

`zmath.vector.intarg` holds the 11 operations that take a non-float argument.
Each extra argument has two forms:

- A broadcast runtime scalar applied to every lane, for example
  `ldexpv(x, exponent: i32)` and `scalbnv(x, exponent: i32)`,
  `scalblnv(x, exponent: i64)`, the power family, and
  `nextTowardv(x, toward: f64)`. The `fromfp` family takes two broadcast
  scalars, a `FromFpRound` direction and a `u32` width.
- A per-lane integer vector under a `Lanes` suffix, holding a `@Vector(N, i32)`
  or `@Vector(N, i64)` of exponents, for example `ldexpvLanes(x, e)` and
  `pownvLanes(x, e)`.

`nextToward` has no `Lanes` form, because its `toward` argument is always
`f64` and each per-format kernel applies its own internal tolerance.

`pownv`, `compoundnv`, and `rootnv` return the lane-ordered
`[N]PowerResult(T)` array described above, in both the broadcast and the
`Lanes` form.

```text
zmath.vector.ldexpv(x, 3)               // one broadcast exponent
zmath.vector.ldexpvLanes(x, e)           // one exponent per lane
zmath.vector.pownvLanes(x, e)            // [N]PowerResult(T), lane-ordered
zmath.vector.ilogbv(x)                   // @Vector(N, i32)
zmath.vector.totalOrderv(a, b)           // @Vector(N, bool)
zmath.vector.atan2v(y, x)                // [N]Atan2Result(T), lane-ordered
```

## IEEE and scalar API


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
The only public functions that do are the `zmath.ieee` `scaleB` and `roundIntegral` families, which name the rounding of their own result.
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

The general bounded power cases use fixed, format-specific integer storage.
They use 32 `u64` limbs for `f16`, 256 limbs for `f32`, and 2,048 limbs for `f64`.
The exact `pown` and `compoundn` paths keep finite dyadic arithmetic exact and round once.
The exact `rootn` path uses candidate and midpoint comparisons instead of a rounded reciprocal-power composition.
These storage bounds prevent allocation. They do not specify a universal execution-time bound.

For binary64 exact rationals with a unit denominator, rounding extracts the
significand bits and uses guard, sticky, and parity bits instead of long division.
Binary64 scaled-integer comparisons read virtually shifted limbs without
materializing shifted operands. Both paths retain the same exact rounding and
candidate decisions; the binary16 and binary32 integer paths are unchanged.

Four binary64 cases use a single strict IEEE operation instead of the integer path:
`pownd(x, -1)` and `rootnd(x, -1)` use `1 / x`,
`rootnd(x, 2)` uses `@sqrt(x)`, and `compoundnd(x, 1)` uses `1 + x`.
The exponent bounds and exceptional-value rules above precede these paths.
Each finite result rounds once in `nearest_even`.
The general reciprocal-square-root case, `rootnd(x, -2)`, retains exact comparisons;
it is not a separately rounded `1 / @sqrt(x)` composition.
The `f16` and `f32` paths are unchanged.

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
It does not inspect, clear, report, or restore floating-point exception flags.
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

### Mode-parameterised scaling and integer rounding

`zmath.ieee` owns the two families that name the rounding of their own result instead of reading the process rounding mode.
`scaleB` takes any `i32` exponent, and both families take the same four directions at every width.

| Family | Signatures |
| --- | --- |
| `scaleB` | `scaleBh(mode: RoundMode, value: f16, exponent: i32) f16`, `scaleBf(mode: RoundMode, value: f32, exponent: i32) f32`, `scaleBd(mode: RoundMode, value: f64, exponent: i32) f64` |
| `roundIntegral` | `roundIntegralh(mode: RoundMode, value: f16) f16`, `roundIntegralf(mode: RoundMode, value: f32) f32`, `roundIntegrald(mode: RoundMode, value: f64) f64` |

`RoundMode` declares exactly four values: `nearest_even`, `toward_zero`, `toward_positive`, and `toward_negative`.
Those four are the modes both families accept, with no default and no fifth direction.
`nearest_even` is ties to an even least-significant significand bit.

`scaleB` returns `value` scaled by `2^exponent`, rounded once in the named mode.
An exact scale never rounds: the significand is unchanged, so only a right shift into the subnormal range can discard anything.
An exponent beyond the representable range is not an error, because the whole `i32` range is a valid argument.
A result that leaves the format saturates at its boundary: infinity when the mode lets the magnitude grow, and the largest finite value when it does not.
A result that rounds to zero keeps the sign of the input, so the underflow of a negative value returns `-0`.

`roundIntegral` returns `value` rounded to an integer in the named mode, with the same tie, direction, and sign rules.
An integral value is returned bit-identical.
A finite fractional value is rounded in the named mode, including values with magnitude greater than one.
A result of zero retains the input sign.

Both families use gradual underflow.
A subnormal input is normalized before it is scaled, a subnormal result is never flushed to zero, and a zero and an infinity pass through at any exponent.
A NaN input preserves its sign and diagnostic payload and sets only the format quiet bit.

This contract previously named no signature for either family, and only the binary16 forms existed.
The `f` and `d` forms are now implemented with the signatures in the table, over the same four `RoundMode` values.
The three widths share one width-generic implementation, so a change to one width is a change to all three.

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

## Complex API

`zmath.complex` owns the canonical `Complexf` and `Complexd` types.

```zig
pub const Complexf = struct { re: f32, im: f32 };
pub const Complexd = struct { re: f64, im: f64 };
```

The complete C99/C23 complex surface is implemented for both types. The module functions use an explicit type suffix:

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

Both suffixes implement `add`, `sub`, `mul`, `div`, `exp`, `log`, `log10`, `pow`, `sqrt`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, and `atanh`.
The manifest names only `complexAdd` and `complexMul`, so the other nineteen pairs are implemented without oracle evidence.

Each entry defines its Zig name, inputs, result, branch cut, principal range, signed-zero result, NaN result, infinity result, pole result, overflow result, and underflow result.

Finite complex results use the scalar manifest accuracy contract for each component after range reduction. A branch-cut or special-value case uses an exact expected component-bit pattern. Complex functions do not return `null` for C mathematical poles.

## Vectors and quaternions

`zmath.vector` owns `Vec2f`, `Vec3f`, `Vec4f`, `Vec2d`, `Vec3d`, and `Vec4d`. Each is a value type with named components.

These are implemented in `src/vector_value.zig` and re-exported from `zmath.vector`. They are not
the `@Vector` layer of [Vector API](#vector-api), which reuses the same module path for a different
purpose.

`zmath.quaternion` owns raw Hamilton types and validated rotation types:

```zig
pub const Quatf = struct { x: f32, y: f32, z: f32, w: f32 };
pub const Quatd = struct { x: f64, y: f64, z: f64, w: f64 };
pub const UnitQuatf = struct { raw: Quatf };
pub const UnitQuatd = struct { raw: Quatd };
```

`w` is the scalar component. `Quatf.mul(a, b)` is Hamilton multiplication. A `UnitQuatf` represents an active, right-handed rotation of a column vector.

`UnitQuatf.rotate(v)` applies `q * (v, 0) * conjugate(q)`. `a.then(b)` means apply `a` first, then `b`, and equals `b * a`.

`UnitQuatf.init(raw)` returns `?UnitQuatf`. It rejects a zero or non-finite input, normalizes a
finite nonzero one, and validates that the normalized result lands in the unit band, so `init` is
the constructor that maintains the unit norm. Raw `Quatf` supports algebraic operations and does
not expose a total rotation operation.

`length` on `Vec2f` through `Vec4d` and on `Quatf` and `Quatd` is an ordered fold
through the accepted `hypot` kernel over the components rather than a square root
of a sum of squares, so it avoids the spurious overflow and the flush to zero that
squaring a component can cause. The first component enters the fold as the exact
value `hypot(+0, x)` returns for it: the magnitude of `x`, or, when `x` is a NaN,
that NaN quieted with its sign and payload preserved. Every later component goes
through the same kernel in the same order. Nesting the kernel is not a
single-rounded multicomponent norm. It is the accepted kernel applied again, never
a re-derived polynomial, and it carries no acceptance of its own here.

`fromAxisAngle` measures its axis once. It rejects a non-finite axis, rejects a
zero or non-finite length, and divides the three components by that length. It
does not normalize the quotient and measure it a second time: a largest component
divided by the accepted finite, nonzero length cannot round to zero, so the
quotient is not the zero vector and a second norm check decides nothing the length
has not already decided. The accepted and rejected calls are the same calls.

`raw` is a public field. Zig cannot make a struct field private, and
`UnitQuatf` is a plain struct with no opaque handle, so `UnitQuatf{ .raw = q }` compiles for any
consumer and can store any quaternion at all. The unit norm is therefore an invariant the
constructors maintain and the consumer preserves, not one the language or the type enforces:
`init` is the only validating constructor, `fromAxisAngle` is the second way in and routes through
`init`, and `raw` is the read path back to the algebra. A consumer must not assign arbitrary values
to `raw` and must not bypass validation to obtain a `UnitQuatf`. `rotate` and the composition
methods rely on that unit-norm precondition, so a forged non-unit value has no rotation guarantee:
the rotation it performs scales by the norm squared, and every `rotate` result from it is
meaningless rather than merely inaccurate. There is no `raw()` accessor method: Zig puts fields and
declarations in one namespace per struct, so a field named `raw` and a method named `raw` cannot
both exist.

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

`Mat.mulVec` applies the matrix to a column vector: it takes exactly `columns`
lanes and returns exactly `rows`. The input is a named vector when a named vector of
that lane count and entry type exists, so `Vec2f`, `Vec3f`, or `Vec4f` for a real
`f32` matrix and `Vec2d`, `Vec3d`, or `Vec4d` for a real `f64` matrix. Every other
real column count, and every `Complexf` or `Complexd` column count, takes the plain
array `[columns]Scalar`. The result is `Mat.vector_type`, which stays row-sized: a
three-row `f32` matrix returns `Vec3f`, a five-row one returns `[5]f32`, and a
complex matrix returns `[rows]Complexf` or `[rows]Complexd`. A square `Mat` is
unchanged by that split, so `Mat3f.mulVec` still takes and returns `Vec3f` and
`Mat4d.mulVec` still takes and returns `Vec4d`, and `vector_type` is still the
right-hand side type `solveLu` and `solveCholesky` accept.

```zig
const rect = Mat(2, 4, f32).init(.{
    .{ 1, 3, 5, 7 },
    .{ 2, 4, 6, 8 },
});
const applied: zmath.vector.Vec2f =
    rect.mulVec(zmath.vector.Vec4f.init(1, 1, 1, 1)); // x = 16, y = 20
```

A rectangular matrix no longer takes a row-sized argument. Handing `Mat(5, 4, f32)`
or `Mat(2, 3, f32)` the five-lane or two-lane vector the old signature accepted is
rejected by type checking at compile time rather than by a checked-build panic, and
there is no compatibility overload for that spelling. The arithmetic does not move
with the shape: each row still seeds from its column-zero product and adds the
remaining columns in increasing order, strictly, with no reassociation and no fused
multiply-add.

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

The implementation supplies the size-paired forms of the three named operations rather than only
`Mat3f`, `Mat4d`, and `Mat4f`: `determinantMat2f` through `determinantMat5d`, `inverseMat2f`
through `inverseMat4d`, and `mulMat4fVec4f` with `mulMat4dVec4d`.
The fixed relative pivot thresholds are `inverse_pivot_threshold_f = 0x1p-21`
and `inverse_pivot_threshold_d = 0x1p-50`.
They derive from unit roundoff and the reduction's element-growth bound.
Callers can inspect these constants. The API has no per-call threshold option.

`inverse` writes only the cells a later step or the result can still read. The
reduction holds `[ value | identity ]` in one row-major work array of a comptime
extent, so nothing is allocated. Step `k` reads column `k` of the left block for
the pivot search and for the multiplier, and writes every column from `k + 1` onward.
Column `j` of the left block is therefore read by the pivot search of step `j` and
by nothing else: after step `j` no cell of column `j` is read again whatever the
entries are, and the result loop reads only the right block. The normalization and
the elimination skip the settled left-block columns, and a row interchange moves
only the live columns, because the columns before `k` are settled in every row.
This is not a value-dependent claim. An update skipped for a settled column would
have written a cell that nothing reads, including a NaN that an overflow would
otherwise have left behind. Every live operation keeps its position in the
sequence, the reduction stays in strict float mode, and `null` still means exactly
the three things it means above: the pre-scan found a non-finite entry, the pivot
magnitude failed the threshold test against `inverse_pivot_threshold_f` or
`inverse_pivot_threshold_d` times the largest entry magnitude, or the reduction
left a non-finite entry in the right block. The threshold, the pivot order, and
the returned value are not changed by the skipped cells.

`zmath.matrix.mul` is declared with `anytype` parameters so that both operand shapes are checked
while the output width is inferred. Zig has no function overloading, so a separate literal
declaration per size pair is not expressible; the returned type is still exactly
`Mat(rows, output_columns, Scalar)`.

## Linear algebra

`zmath.linalg` owns named factorization types for named matrix sizes.

```zig
pub const LuMat2f = struct { lu: Mat2f, permutation: [2]u8, parity: i1 };
pub const LuMat3f = struct { lu: Mat3f, permutation: [3]u8, parity: i1 };
pub const LuMat4f = struct { lu: Mat4f, permutation: [4]u8, parity: i1 };
pub const CholeskyMat4f = struct { lower: Mat4f };
```

The factorization types are plain value structs with no methods. Their operations are module
functions: `lu` and `cholesky` produce a factorization, `solveLu` and `solveCholesky` return
`?Vec2f` through `?Vec4d`, `solveLuMultiple` and `solveCholeskyMultiple` carry a right-hand side
matrix, and `reconstructLu`, `determinantLu`, `reconstructCholesky`, and `determinantCholesky`
read a factorization back. `lu` and `cholesky` return `null` for a singular, non-positive-definite,
or non-finite input, as applicable. A solve returns `null` for a zero or non-finite
stored pivot. It does not reject non-finite computed entries: finite right-hand sides
can overflow to a successful result containing infinity or NaN.
Single right-hand sides must have the matrix's exact `vector_type`.
Multiple right-hand sides must have the factorization's row count and scalar type.
The result is a new value. A failure writes no partial result into caller storage.

The implementation extends this to every named size in both precisions: `LuMat2f` through
`LuMat4d` and `CholeskyMat2f` through `CholeskyMat4d`, twelve factorization types rather than the
four named above. `lu` and `cholesky` take `anytype` for the same reason `zmath.matrix.mul` does.
There is also a multiple-right-hand-side solve for each factorization family. `lu` is Doolittle
with partial pivoting, so the stored diagonal of U is the pivot sequence and the determinant is
`parity` times the product of that diagonal; `zmath.linalg` reuses the pivot thresholds
`zmath.matrix` publishes instead of declaring a second spelling of the same constant.

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

`fast` reports the sign of the rounded scalar polynomial: infinity retains its sign,
while NaN returns `error.NonFiniteResult`.
`adaptive_exact` reports the exact polynomial sign for finite `f32` or `f64` inputs.
It can resolve a sign when `fast` rounds to zero or fails.
Floating constructions can still fail under either policy.
Predicates do not use an implicit epsilon.
Checked geometry rejects non-finite arguments and receiver fields with `error.NonFiniteInput`.
It rejects non-finite computed results or required construction intermediates with `error.NonFiniteResult`.
These checks apply in every build mode and return no partial result.
Infallible queries and conversions (`isFinite`, `eq`, `toArray`) and `Bivector4.neg`
do not perform checked validation and can inspect or preserve non-finite fields.

Geometry fields are writable. Constructors establish finite coordinates, unit directions or normals,
nonnegative radii, and ordered bounds where applicable.
Callers must preserve the geometric invariants after direct field writes.
Checked calls validate finiteness. They do not re-run every constructor invariant check.

The implementation adds three tags beyond the nine listed above, because the nine do not name
every result the families can produce. `box` carries the common AABB of two AABBs.
`circle` and `hypersphere` carry center/radius summaries of a flat's ball cross-section.
Use the original flat to recover its placement and orientation; these payloads alone
do not represent an oriented disk or three-dimensional ball in ambient space.
Applicability follows the dimension. The 2D union has no flat tags.
The 3D union has no three-flat tag.
A coincident pair of two-sided infinite flats returns its span tag.
Identical rays, segments, balls, and boxes return `coincident`.
A ray's supporting line is not its whole common set.
`contained` identifies which entire operand is the common set.
A solid ball intersected by a ray returns its entire closed chord as `segment`,
including an inside ray's origin; a zero-length chord is `point`.
Proper ball/ball overlaps retain `two_points` boundary summaries. In 2D these are
the circumference crossings, not the common solid lens. In 3D and 4D they are
selected boundary points, not the complete boundary or solid intersection.

`adaptive_exact` is exact integer arithmetic, not an error-bounded approximation with a fast path.
Every finite `f64` is `m * 2^e` with `|m| < 2^53` and `-1074 <= e <= 971`, so scaling by `2^1074`
reaches an exact integer; the predicate is then evaluated with schoolbook multiplication in a
fixed `[132]u64` limb magnitude plus a separate sign, 8,448 bits on the stack, and the sign of
that integer is the sign of the real value. An `f32` widens to `f64` exactly, so one
representation covers both precisions. There is no epsilon, no value-dependent limb count, and
no fallback to the `fast` answer.
Coordinate subtraction happens after integer conversion, so rounding or overflow in the corresponding
floating difference cannot change the exact sign. Mixed-degree ball and ray comparisons align their
power-of-two scales before comparison. The fixed scratch bound covers these rescaled predicates.

One limit is stated rather than hidden. In 3D and 4D the affine intersections decide coplanarity
with the exact predicate on the ORIGINAL coordinates: in 3D that gate is the exact orientation
determinant, and in 4D it is an exact rank test on three difference vectors. They then project
into a plane and run the exact 2D predicate on the PROJECTED coordinates. The coplanarity gate,
which is the decision a plain evaluation gets wrong, is exact on the inputs as given.
Intersection PARAMETERS are `f32`/`f64` and their range tests are ordinary comparisons.

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

The signed complementary Poincare dual is implemented as: a grade `k` blade of mask `m` maps to
the blade of mask `(top & ~m)` with the sign `(-1)^(k * (n - 1 - k))`, where `n` is the basis
count. The double dual then carries the sign
`(-1)^(k(n-1-k) + (n-k)(k-1)) = (-1)^(2k(n-k) - n) = (-1)^n`, so `dual(dual(A))` is `A` for Pga3
and `-A` for Pga4. That involution is the test that pins the sign rule.

Blade product tables are generated at comptime from the metric vector rather than written out, so
the `Cga4` table is not a hand-maintained 64 by 64 constant. Each of the twelve types is a distinct
concrete type built from one private comptime generic, so no public signature carries `anytype` or
a type parameter. Storage is `[2^n]T` and is published per type as `storage_bytes`: 32 bytes for
`Cl3f` up to 512 for `Cga4d`. Those are stack value types, not allocations.

`dot` is `scalarProduct(a, b.reverse())`, and `scalarProduct` reads the diagonal
of the generated product table rather than all `coefficient_count` squared entries
of it. A blade pair lands on the scalar coefficient only when it multiplies with
itself, so every other pair lands above grade `0` and is discarded by the scalar
form. The surviving terms are the same terms, in the same ascending mask order,
with the same `+0` seed, and with the same `coeff * a * b` operand order that
`geometricProduct` uses when it accumulates into its scalar slot, so the value is
what `a.geometricProduct(b).c[0]` is. Both skips of the full product are kept. A
zero left coefficient contributes nothing, and a zero metric entry contributes
nothing, which is the degenerate-metric rule: the `Pga` basis vector `e0` has a
zero entry on the diagonal, so a `Pga` dot product still ignores it exactly as the
full product does. `dot` stays symmetric and `norm` stays `dot(self, self)`.

`inverseChecked` uses the even Clifford conjugate, forms its product with `self`,
and checks that product for a finite, nonzero scalar sandwich. It rejects non-finite
input and a non-scalar sandwich. This is not a general invertibility test: an
odd-only blade returns `null`. Final scaled coefficients are not checked and can
overflow. The unchecked `inverse` requires the same scalar-sandwich precondition.

`evenPart` retains all even grades, including grade 6 where supported.
`unitRotor`, `unitVersor`, and `unitFlat` normalize candidate values; they do not
prove rotor membership, versor membership, or blade decomposability. The caller
must establish those properties. A single grade or an even grade alone is not
such a proof.

## Tests and benchmarks

`zig build test` runs the numerical unit suite. `zig build test-build` installs
the same executable without running it. The optional `cross-test` host tool
uses [cross-test.profile.sh](cross-test.profile.sh) to build that executable
for its selected targets. x86 targets use `baseline+fma`.

Oracle generation, stored-corpus replay, fuzzing, and benchmarks are in
[zmath-testing](https://github.com/CoffeeImpliesCode/zmath-testing).
Its [scalar contract](spec/core-math-2026-08-22.md) defines the corpus schema
and target profile. Source exports and unit tests do not establish
whole-domain accuracy. See the harness documentation for its commands.
