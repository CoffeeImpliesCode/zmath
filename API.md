# zmath API Contract

## Status

This document defines the target API and identifies the currently implemented public subset.

The current public root exports the `zmath.ieee`, `zmath.scalar`, `zmath.vector`, `zmath.complex`,
`zmath.quaternion`, `zmath.matrix`, `zmath.linalg`, `zmath.geometry`, and `zmath.algebra` modules,
their shared result types, and the C23 source surface defined in
`spec/core-math-2026-08-22.md`. Every module path listed under
[Planned module paths](#planned-module-paths) is now implemented.
`zmath.vector` is implemented as the portable `@Vector` layer described in [Vector API](#vector-api):
417 public entry points covering the elementwise scalar operations, including the result-struct,
integer-argument, integer-returning, boolean-returning, and narrowing ones. The same module also
owns the linear-algebra `Vec*` value types from [Vectors and quaternions](#vectors-and-quaternions);
the two layers share the module path and nothing else, and every `@Vector` entry point ends in `v`
while no `Vec*` name does.

## Scope

`zmath` is a dependency-free Zig 0.16 mathematics library. Public fixed-size operations do not allocate.

Public real formats are `f16`, `f32`, and `f64` for the current IEEE and scalar primitives. Public type and function names use lower-case `h`, `f`, and `d` suffixes. `q`, `f128`, and bfloat16 are outside this contract and have no public API.

The API uses one structural generic, `Mat`. Complex values, vectors, quaternions, geometry, factorizations, and geometric algebras use named concrete types.

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

## Vector API

`zmath.vector` provides `@Vector` overloads for the elementwise scalar
operations: 417 public entry points across six modules. It is a distinct
layer from the linear-algebra `Vec*` value types described in
[Vectors and quaternions](#vectors-and-quaternions), which the same module
also owns and which are implemented in `src/vector_value.zig`. The two layers
share this module path and nothing else.

Every vector operation applies the identical accepted scalar kernel to each
lane. Lane `i` of a vector result is therefore bit-for-bit equal to the scalar
result for lane `i` of the input. That is what preserves the `max_ulp = 0`
contract and the immutable corpus, so the vector surface requires no new
acceptance evidence under the independent [status terms](#status-terms): the
vector path is a per-lane application of the same scalar code, never a
re-derived polynomial.

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

The 13 result-struct operations cannot be a `@Vector`. Zig 0.16 forbids a
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

### Typed scalar corpus manifest

Before a manifest-gated scalar pair becomes accepted, the repository updates a versioned contract at `spec/core-math-<date>.md`.

The schema-5 manifest is `zmath-testing/test/corpus/manifest.toml`.
It contains exactly 347 sorted typed function rows for `f16`, `f32`, and `f64`.
It pins the CORE-MATH revision, 512-bit MPFR/Arb settings for evaluator paths that use them, `nearest_even`, and the special-value policy.
The checked-in manifest is the only acceptance authority.
A successful candidate campaign does not change a row's manifest status.

Each row defines the operation, provenance, status, arity, typed inputs, typed result, oracle, special-value policy, and maximum ULP bound.
The root `operation_spelling` value is `"c23_hfd_lowercase_v1"`.
For a standard H/F/D family, the typed descriptors select the suffix.
Dedicated encoded names, including narrowing operations, remain unchanged.

The provenance vocabulary is exactly `["zmath", "core_math"]`.
Each row uses `["zmath"]`, `["core_math"]`, or `["zmath", "core_math"]`.
Scalar results use `float:<format>`, `integer:<type>`, `bool`, or `classification`.
Structure results use `struct:<kind>:<format>`.
Status-container descriptors include:

- `struct:atan2:<format>` for `Atan2Result(T) { value: T, resolved: bool }`.
- `struct:power:<format>` for `PowerResult(T) { value: T, ok: bool }`.
- `struct:integer:<format>` for `IntegerResult { value: i64, ok: bool }`.
- `struct:quaternion:<format>` for four quaternion components in `x`, `y`, `z`, `w` order.
- `struct:algebra:<format>` for eight Cl(3) coefficients in basis-mask order.
- `struct:complex:<format>` for the two components of `Complexf` or `Complexd`, `re` then `im`.
- `struct:luSolve:<format>` for three stored channels, `x0`, `x1`, and `ok`.

These descriptors name the layout a stored record uses.
Most of them also carry the public return shape.
Two of them name no public type of their own.
`zmath.complex` has no generic `Complex(T)`: its only public types are the concrete
`Complexf` and `Complexd`, each with `re` and `im`, so `struct:complex` mirrors
those two fields and promises no third type.
`zmath.linalg` publishes no `LuSolve2Result(T)`: `lu` returns `?LuMat2f` through
`?LuMat4d` and `solveLu` returns an optional vector, `?Vec2f` through `?Vec4d`.
`struct:luSolve` is the adapter that turns that optional vector into stored
fields. A refused factorization or solve publishes the canonical `+0`, `+0`, and
`false`, because a `null` has no record encoding, and an accepted solve whose
exact solution leaves the format may publish an infinite component with
`ok=true`.
No descriptor grants implementation, corpus evidence, or manifest acceptance.

Ten rows cover selected non-scalar operations: five operations, each with one `f`
and one `d` row. `complexAdd` and `complexMul` take four same-format floats `a.re`,
`a.im`, `b.re`, `b.im`; `matrixDet2` takes four row-major `a00`, `a01`, `a10`,
`a11`; `luSolve2` takes `a00`, `a01`, `a10`, `a11`, `b0`, `b1`; and
`geometryOrient2` takes six coordinates `ax`, `ay`, `bx`, `by`, `cx`, `cy` plus a
policy tag of exactly `adaptive_exact` or `fast`, and returns `integer:i32`
carrying `-1`, `0`, or `+1` rather than a sign-enum ordinal. All ten use
`oracle = "algebraic"`, `status = "correct"`, and `max_ulp = 0`.
The independent exact-rational reference follows the arithmetic each public graph
performs, so the ten rows do not share one rounding rule: `complexAdd` rounds each
exact component sum once, `complexMul` rounds each product and each component sum,
`matrixDet2` rounds two products and the final subtraction, `luSolve2` computes the
exact rational Cramer solution and rounds once per component, and `geometryOrient2`
uses the exact determinant of the stored coordinates under `adaptive_exact` and the
target-rounded graph under `fast`.
These rows hold 324 stored records over bounded named cases, and they are the
approved selected coverage rather than a claim about the remaining operations of
those modules.

The status values are `correct`, `bounded`, and `unsupported`.
An `unsupported` row uses `oracle = "none"` and `max_ulp = 0`.
An accepted row uses the evaluator class assigned by the authoritative C operation table.
A `correct` row uses `max_ulp = 0`.
A `bounded` row uses an approved non-negative `max_ulp` value.
Schema 5 accepts `none`, `exact_bits`, `algebraic`, `mpfr_arb`, and `definitional`.

Status and oracle snapshots are not part of numerical descriptor identity.
This permits review of candidate artifacts before a manifest promotion.
The loader still checks that each snapshot has a valid status-oracle relation.
It hashes the profile, numerical descriptors, special policy, ULP bound, and authoritative C signature.

The current manifest accepts all 347 descriptors as `correct` with `max_ulp = 0`.
Full publication and replay verified 72,035,642 records on the scalar lane and on the scalar-backed vector lane.
The scalar benchmark receipt covers all 347 descriptors of this publication.
It retains only repetition zero for each descriptor.
A benchmark result is a timing record, not acceptance or whole-domain proof.
Do not infer acceptance from source exports, evaluator availability, or a filtered campaign.

#### Oracle evidence classes

The authoritative C operation table selects one evaluator for each of the 347 descriptors.
It assigns 150 descriptors to `exact_bits`, 126 to `mpfr_arb`, 49 to `algebraic`, and 22 to `definitional`.
Operation spelling and manifest status do not select an evaluator.

- `exact_bits` uses direct IEEE bit and integer models. It preserves bit-sensitive NaN and status behavior without numerical approximation.
- `algebraic` uses exact GMP rational arithmetic for finite dyadic operations. Irrational algebraic roots use 512-bit MPFR candidates and 512-bit Arb enclosures.
- `mpfr_arb` uses a 512-bit MPFR candidate and a 512-bit Arb proof, with operation-specific certified reductions where required.
- `definitional` checks IEEE representations or published algebraic conventions without an MPFR/Arb enclosure.
  Its 22 rows cover the four public helpers in three formats and five algebra/quaternion operations in two formats.
  These rows do not establish general non-scalar numerical accuracy.

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
The generator compares the zmath result with the declared evaluator before it writes a record.

Published records use the length-framed binary `ZCDR` codec version 1.
A record stores operation, ordered provenance, typed arguments, result descriptor, typed result fields, and `max_ulp`.
Floating fields store the exact format, bits, and IEEE class.
Other fields store their integer width, Boolean value, tag bytes, or rounding enumeration.

A full publication covers all 347 descriptors.
A filtered publication covers the exact `--operation` and `--format` selection.
Metadata records the campaign class, seed, cardinality, selection rule, profile, and artifact digest.
The 62 exhaustive unary H campaigns each enumerate all 65,536 input bit patterns.
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

`raw` is a public field, not a private one. Zig 0.16 cannot make a struct field private, and
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
through `inverseMat4d`, and `mulMat4fVec4f` with `mulMat4dVec4d`. The singularity threshold for
`inverse` is published as `inverse_pivot_threshold_f` and `inverse_pivot_threshold_d`, derived as
the unit roundoff of the format times the element-growth bound of the reduction, so a caller can
see and change it rather than discover an implicit epsilon.

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
or non-finite input, as applicable, and a solve returns `null` for a zero or non-finite stored
pivot.

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

`fast` reports the sign of the scalar computation. `adaptive_exact` reports the exact sign for finite f32 or f64 inputs. Predicates do not use an implicit epsilon. Every geometry operation rejects non-finite input with `error.NonFiniteInput`.

The implementation adds three tags beyond the nine listed above, because the nine do not name
every result the families can produce. `box` carries an AABB against an AABB, which is an AABB
and not a point set. `circle` and `hypersphere` carry a flat cutting a ball, which is again a ball
rather than a finite point set. Applicability still follows the dimension: 2D has no flat, so its
union has neither `plane` nor `hyperplane` nor `circle`; 3D has no three-flat, so it has no
`hyperplane`. A coincident pair of INFINITE flats returns the tag naming its own span, because
"the intersection is this line" is the more useful answer than "coincident"; `coincident` is
reserved for identical segments, balls, and boxes.

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

`inverseChecked` reuses the product it has already formed. It computes the sandwich
conjugate, forms `conjugate_part.geometricProduct(self)` once, scans that product
for a nonzero coefficient above grade `0`, and then reads its grade-0 part as the
sandwich scalar instead of building a second product to read the same scalar from.
The residue scan has to read the whole product anyway, so the grade-0 part is
already in hand. `conjugate_part.scale(1 / sandwich)` is then the same expression
`inverse` builds from its own second product, and the refusals are unchanged: a
non-finite coefficient in `self`, a nonzero coefficient above grade `0` in the
product, and a zero or non-finite sandwich scalar each return `null`. The
precondition `inverse` documents is untouched, and `inverseChecked` remains the
test for it.

## Tests and benchmarks

`zmath` has no MPFR, GMP, FLINT, Arb, BLIS, Eigen, or C++ runtime dependency.
The sibling `zmath-testing` project owns integration tests, oracle generation, typed corpus replay, and benchmarks.
It is the acceptance authority and is colocated as `../zmath-testing`.
Run the following commands there, not in this repository.
The current publication has 347 descriptors and 72,035,642 records.
See [PLAN.md](PLAN.md) for measured evidence and remaining scope.

That authority covers the scalar manifest. Evidence gathered for the implemented
consumer surfaces is a separate thing and carries none of it. The optimization
round that reworked those surfaces is gated by paired reference comparison over
bounded named inputs rather than by a publication, and its gate is bounded
globally at 4096 evaluations per canonical f32 operation counted across every lane,
consumer, repetition, and reference together, not once per policy and not once per
fixture. Bounded consumer evidence is not a publication under
[Typed scalar corpus manifest](#typed-scalar-corpus-manifest): it changes no
descriptor row, status, oracle, ULP bound, or stored record, it is not a campaign
the corpus loader accepts, and it establishes neither corpus coverage nor
whole-domain proof, whatever such a comparison reports. At its strongest it says
that a reworked path reproduced the reference it was compared against over the
inputs it was compared on, and it says nothing about the rest of an
operation-format domain. The bound is what makes the round checkable without a
publication, and it relaxes no threshold stated elsewhere in this document.

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

`zmath-public-helpers` measures `classify`, `ulpDistance`, and the binary16
`scaleB` and `roundIntegral` entry points in all four rounding modes.
The manifest now covers these helpers under `definitional` rows.
The helper benchmark takes no corpus root and no
operation filter, and accepts only `--repetitions=1..100`. The typed benchmark
accepts `--repetitions=1..100` and `--with-zig`, which times same-input Zig
intrinsic references beside zmath. 1024 is the typed benchmark's internal
sample-loop cap, not a command-line range.

The corpus producer, corpus replay, and scalar-matrix executables accept an
optional first positional corpus root. The default root is `test/corpus`.
These three executables also accept `--operation=<manifest-op>` and
`--format=h|f|d`. The filters can be used independently or together.
Each filtered command uses its own `generated/probes/<op>-<fmt>` publication.
An omitted filter uses `all` in that path.
Generate that probe before filtered replay or measurement.

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
