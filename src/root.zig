//! By convention, root.zig is the root source file when making a package.
const std = @import("std");

pub const ieee = @import("ieee.zig");
pub const scalar = @import("scalar.zig");
pub const vector = @import("vector.zig");
pub const complex = @import("complex.zig");
pub const quaternion = @import("quaternion.zig");
pub const algebra = @import("algebra.zig");
pub const matrix = @import("matrix.zig");
pub const geometry = @import("geometry.zig");
pub const linalg = @import("linalg.zig");
const generic_scalar_test = @import("generic_scalar_test.zig");
const c23_public_surface_test = @import("c23_public_surface_test.zig");
const scalar_c23_types = @import("scalar_c23_types.zig");
const scalar_narrow = @import("scalar_narrow.zig");

pub const FromFpRound = scalar_c23_types.FromFpRound;
pub const IntegerResult = scalar_c23_types.IntegerResult;
pub const Atan2Result = scalar_c23_types.Atan2Result;
pub const PowerResult = scalar_c23_types.PowerResult;
pub const FrexpResult = scalar_c23_types.FrexpResult;
pub const ModfResult = scalar_c23_types.ModfResult;
pub const RemquoResult = scalar_c23_types.RemquoResult;
pub const LgammaResult = scalar_c23_types.LgammaResult;
pub const CanonicalizeResult = scalar_c23_types.CanonicalizeResult;
pub const PayloadResult = scalar_c23_types.PayloadResult;

// ---------------------------------------------------------------------------
// Lane-parallel result shapes
//
// The operations below route a `@Vector` argument to the lane-parallel entry
// point of their family in `vector.zig`, so each of them must spell BOTH the
// scalar result and the vector result. These three helpers are the whole
// adjustment: every one of them leaves the scalar result type exactly as it
// was and derives the vector result type from the argument alone.
// ---------------------------------------------------------------------------

/// The lane float type of `X`: the `@Vector` element type for a vector
/// argument, and `X` itself for a scalar one.
///
/// A result struct such as `FrexpResult(T)` is parameterized by the lane type,
/// never by the vector type, so a `@Vector(N, f32)` argument must bind
/// `FrexpResult(f32)`.
fn FloatLane(comptime X: anytype) type {
    return switch (@typeInfo(X)) {
        .vector => |v| v.child,
        else => X,
    };
}

/// The result type of an operation whose per-lane result is the struct `R`:
/// `R` itself for a scalar argument, and a lane-ordered `[N]R` for a vector
/// one.
///
/// Zig forbids a struct as a `@Vector` element type and reports "expected
/// integer, float, bool, or pointer for the vector element type", so the
/// vector result here is an array, never a `@Vector`. Element `i` is the value
/// the scalar kernel returned for lane `i`.
fn LaneResult(comptime R: type, comptime X: anytype) type {
    return switch (@typeInfo(X)) {
        .vector => |v| [v.len]R,
        else => R,
    };
}

/// The result type of an operation whose per-lane result is the non-struct
/// value `R`: `R` itself for a scalar argument, and a real `@Vector(N, R)` for
/// a vector one.
///
/// This is the shape of the integer-returning `ilogb` and `llogb` and of the
/// boolean-returning `totalOrder` and `totalOrderMag`.
fn LaneVector(comptime R: type, comptime X: anytype) type {
    return switch (@typeInfo(X)) {
        .vector => |v| @Vector(v.len, R),
        else => R,
    };
}

/// Computes the absolute value of `value` with the scalar implementation for
/// its type.
pub fn abs(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.absv(value);
    if (T == f16) return scalar.absh(value);
    if (T == f32) return scalar.absf(value);
    if (T == f64) return scalar.absd(value);
    @compileError("zmath.abs supports only f16, f32, and f64");
}

/// Computes e raised to `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn exp(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.expv(value);
    if (T == f16) return scalar.exph(value);
    if (T == f32) return scalar.expf(value);
    if (T == f64) return scalar.expd(value);
    @compileError("zmath.exp supports only f16, f32, and f64");
}

/// Computes `base` raised to `exponent` with the scalar implementation for
/// their common type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn pow(base: anytype, exponent: anytype) @TypeOf(base) {
    const T = @TypeOf(base);
    if (@typeInfo(T) == .vector) return vector.powv(base, exponent);
    if (@TypeOf(exponent) != T) {
        @compileError("zmath.pow requires base and exponent to have the same type");
    }
    if (T == f16) return scalar.powh(base, exponent);
    if (T == f32) return scalar.powf(base, exponent);
    if (T == f64) return scalar.powd(base, exponent);
    @compileError("zmath.pow supports only f16, f32, and f64");
}

/// Computes 10 raised to `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn exp10(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.exp10v(value);
    if (T == f16) return scalar.exp10h(value);
    if (T == f32) return scalar.exp10f(value);
    if (T == f64) return scalar.exp10d(value);
    @compileError("zmath.exp10 supports only f16, f32, and f64");
}

/// Computes the base-2 exponential biased by one with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn exp2m1(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.exp2m1v(value);
    if (T == f16) return scalar.exp2m1h(value);
    if (T == f32) return scalar.exp2m1f(value);
    if (T == f64) return scalar.exp2m1d(value);
    @compileError("zmath.exp2m1 supports only f16, f32, and f64");
}

/// Computes the base-10 exponential biased by one with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn exp10m1(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.exp10m1v(value);
    if (T == f16) return scalar.exp10m1h(value);
    if (T == f32) return scalar.exp10m1f(value);
    if (T == f64) return scalar.exp10m1d(value);
    @compileError("zmath.exp10m1 supports only f16, f32, and f64");
}

/// Computes `a` multiplied by `b` then added to `c` with the scalar implementation
/// for their common type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn fma(a: anytype, b: anytype, c: anytype) @TypeOf(a) {
    const T = @TypeOf(a);
    if (@typeInfo(T) == .vector) return vector.narrow.fmav(a, b, c);
    if (@TypeOf(b) != T or @TypeOf(c) != T) {
        @compileError("zmath.fma requires all arguments to have the same type");
    }
    if (T == f16) return scalar.fmah(a, b, c);
    if (T == f32) return scalar.fmaf(a, b, c);
    if (T == f64) return scalar.fmad(a, b, c);
    @compileError("zmath.fma supports only f16, f32, and f64");
}

/// Computes the complementary error function of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn erfc(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.erfcv(value);
    if (T == f16) return scalar.erfch(value);
    if (T == f32) return scalar.erfcf(value);
    if (T == f64) return scalar.erfcd(value);
    @compileError("zmath.erfc supports only f16, f32, and f64");
}

/// Computes the true gamma function of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn tgamma(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.tgammav(value);
    if (T == f16) return scalar.tgammah(value);
    if (T == f32) return scalar.tgammaf(value);
    if (T == f64) return scalar.tgammad(value);
    @compileError("zmath.tgamma supports only f16, f32, and f64");
}

/// Computes the natural logarithm of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn log(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.logv(value);
    if (T == f16) return scalar.logh(value);
    if (T == f32) return scalar.logf(value);
    if (T == f64) return scalar.logd(value);
    @compileError("zmath.log supports only f16, f32, and f64");
}
/// Computes the base-2 logarithm of `value` with the scalar implementation for
/// its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn log2(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.log2v(value);
    if (T == f16) return scalar.log2h(value);
    if (T == f32) return scalar.log2f(value);
    if (T == f64) return scalar.log2d(value);
    @compileError("zmath.log2 supports only f16, f32, and f64");
}
/// Computes the base-2 logarithm of one plus `value` with the scalar
/// implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn log2p1(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.log2p1v(value);
    if (T == f16) return scalar.log2p1h(value);
    if (T == f32) return scalar.log2p1f(value);
    if (T == f64) return scalar.log2p1d(value);
    @compileError("zmath.log2p1 supports only f16, f32, and f64");
}

/// Computes the base-10 logarithm of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn log10(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.log10v(value);
    if (T == f16) return scalar.log10h(value);
    if (T == f32) return scalar.log10f(value);
    if (T == f64) return scalar.log10d(value);
    @compileError("zmath.log10 supports only f16, f32, and f64");
}
/// Computes the base-10 logarithm of one plus `value` with the scalar
/// implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn log10p1(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.log10p1v(value);
    if (T == f16) return scalar.log10p1h(value);
    if (T == f32) return scalar.log10p1f(value);
    if (T == f64) return scalar.log10p1d(value);
    @compileError("zmath.log10p1 supports only f16, f32, and f64");
}

/// Computes the error function of `value` with the scalar implementation for
/// its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn erf(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.erfv(value);
    if (T == f16) return scalar.erfh(value);
    if (T == f32) return scalar.erff(value);
    if (T == f64) return scalar.erfd(value);
    @compileError("zmath.erf supports only f16, f32, and f64");
}

/// Computes the hyperbolic sine of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn sinh(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.sinhv(value);
    if (T == f16) return scalar.sinhh(value);
    if (T == f32) return scalar.sinhf(value);
    if (T == f64) return scalar.sinhd(value);
    @compileError("zmath.sinh supports only f16, f32, and f64");
}

/// Computes the sine of `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn sin(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.sinv(value);
    if (T == f16) return scalar.sinh(value);
    if (T == f32) return scalar.sinf(value);
    if (T == f64) return scalar.sind(value);
    @compileError("zmath.sin supports only f16, f32, and f64");
}
/// Computes the sine of `value` in half-revolutions with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn sinpi(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.sinpiv(value);
    if (T == f16) return scalar.sinpih(value);
    if (T == f32) return scalar.sinpif(value);
    if (T == f64) return scalar.sinpid(value);
    @compileError("zmath.sinpi supports only f16, f32, and f64");
}

/// Computes the cosine of `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn cos(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.cosv(value);
    if (T == f16) return scalar.cosh(value);
    if (T == f32) return scalar.cosf(value);
    if (T == f64) return scalar.cosd(value);
    @compileError("zmath.cos supports only f16, f32, and f64");
}

/// Computes cos(pi * `value`) with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn cospi(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.cospiv(value);
    if (T == f16) return scalar.cospih(value);
    if (T == f32) return scalar.cospif(value);
    if (T == f64) return scalar.cospid(value);
    @compileError("zmath.cospi supports only f16, f32, and f64");
}
/// Computes the hyperbolic cosine of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn cosh(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.coshv(value);
    if (T == f16) return scalar.coshh(value);
    if (T == f32) return scalar.coshf(value);
    if (T == f64) return scalar.coshd(value);
    @compileError("zmath.cosh supports only f16, f32, and f64");
}
/// Computes the inverse hyperbolic cosine of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn acosh(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.acoshv(value);
    if (T == f16) return scalar.acoshh(value);
    if (T == f32) return scalar.acoshf(value);
    if (T == f64) return scalar.acoshd(value);
    @compileError("zmath.acosh supports only f16, f32, and f64");
}

/// Computes the tangent of `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn tan(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.tanv(value);
    if (T == f16) return scalar.tanh(value);
    if (T == f32) return scalar.tanf(value);
    if (T == f64) return scalar.tand(value);
    @compileError("zmath.tan supports only f16, f32, and f64");
}

/// Computes the tangent of pi times `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn tanpi(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.tanpiv(value);
    if (T == f16) return scalar.tanpih(value);
    if (T == f32) return scalar.tanpif(value);
    if (T == f64) return scalar.tanpid(value);
    @compileError("zmath.tanpi supports only f16, f32, and f64");
}

/// Computes the hyperbolic tangent of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn tanh(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.tanhv(value);
    if (T == f16) return scalar.tanhh(value);
    if (T == f32) return scalar.tanhf(value);
    if (T == f64) return scalar.tanhd(value);
    @compileError("zmath.tanh supports only f16, f32, and f64");
}

/// Computes the reciprocal square root of `value` with the scalar
/// implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn rsqrt(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.rsqrtv(value);
    if (T == f16) return scalar.rsqrth(value);
    if (T == f32) return scalar.rsqrtf(value);
    if (T == f64) return scalar.rsqrtd(value);
    @compileError("zmath.rsqrt supports only f16, f32, and f64");
}

/// Computes the square root of `value` with the scalar implementation for its
/// type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn sqrt(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.sqrtv(value);
    if (T == f16) return scalar.sqrth(value);
    if (T == f32) return scalar.sqrtf(value);
    if (T == f64) return scalar.sqrtd(value);
    @compileError("zmath.sqrt supports only f16, f32, and f64");
}

/// Computes the hypotenuse of `x` and `y` with the scalar implementation for
/// their common type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn hypot(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.hypotv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.hypot requires both arguments to have the same type");
    }
    if (T == f16) return scalar.hypoth(x, y);
    if (T == f32) return scalar.hypotf(x, y);
    if (T == f64) return scalar.hypotd(x, y);
    @compileError("zmath.hypot supports only f16, f32, and f64");
}

/// Computes the inverse cosine of `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn acos(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.acosv(value);
    if (T == f16) return scalar.acosh(value);
    if (T == f32) return scalar.acosf(value);
    if (T == f64) return scalar.acosd(value);
    @compileError("zmath.acos supports only f16, f32, and f64");
}

/// Computes the inverse sine of `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn asin(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.asinv(value);
    if (T == f16) return scalar.asinh(value);
    if (T == f32) return scalar.asinf(value);
    if (T == f64) return scalar.asind(value);
    @compileError("zmath.asin supports only f16, f32, and f64");
}

/// Computes the inverse hyperbolic sine of `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn asinh(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.asinhv(value);
    if (T == f16) return scalar.asinhh(value);
    if (T == f32) return scalar.asinhf(value);
    if (T == f64) return scalar.asinhd(value);
    @compileError("zmath.asinh supports only f16, f32, and f64");
}

/// Computes the inverse hyperbolic tangent of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn atanh(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.atanhv(value);
    if (T == f16) return scalar.atanhh(value);
    if (T == f32) return scalar.atanhf(value);
    if (T == f64) return scalar.atanhd(value);
    @compileError("zmath.atanh supports only f16, f32, and f64");
}

pub fn atan(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.atanv(value);
    if (T == f16) return scalar.atanh(value);
    if (T == f32) return scalar.atanf(value);
    if (T == f64) return scalar.atand(value);
    @compileError("zmath.atan supports only f16, f32, and f64");
}

pub fn cbrt(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.cbrtv(value);
    if (T == f16) return scalar.cbrth(value);
    if (T == f32) return scalar.cbrtf(value);
    if (T == f64) return scalar.cbrtd(value);
    @compileError("zmath.cbrt supports only f16, f32, and f64");
}

pub fn exp2(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.exp2v(value);
    if (T == f16) return scalar.exp2h(value);
    if (T == f32) return scalar.exp2f(value);
    if (T == f64) return scalar.exp2d(value);
    @compileError("zmath.exp2 supports only f16, f32, and f64");
}

pub fn expm1(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.expm1v(value);
    if (T == f16) return scalar.expm1h(value);
    if (T == f32) return scalar.expm1f(value);
    if (T == f64) return scalar.expm1d(value);
    @compileError("zmath.expm1 supports only f16, f32, and f64");
}

pub fn log1p(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.log1pv(value);
    if (T == f16) return scalar.log1ph(value);
    if (T == f32) return scalar.log1pf(value);
    if (T == f64) return scalar.log1pd(value);
    @compileError("zmath.log1p supports only f16, f32, and f64");
}

pub fn logp1(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.logp1v(value);
    if (T == f16) return scalar.logp1h(value);
    if (T == f32) return scalar.logp1f(value);
    if (T == f64) return scalar.logp1d(value);
    @compileError("zmath.logp1 supports only f16, f32, and f64");
}

pub fn ceil(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.ceilv(value);
    if (T == f16) return scalar.ceilh(value);
    if (T == f32) return scalar.ceilf(value);
    if (T == f64) return scalar.ceild(value);
    @compileError("zmath.ceil supports only f16, f32, and f64");
}

pub fn floor(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.floorv(value);
    if (T == f16) return scalar.floorh(value);
    if (T == f32) return scalar.floorf(value);
    if (T == f64) return scalar.floord(value);
    @compileError("zmath.floor supports only f16, f32, and f64");
}

pub fn nearbyint(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.nearbyintv(value);
    if (T == f16) return scalar.nearbyinth(value);
    if (T == f32) return scalar.nearbyintf(value);
    if (T == f64) return scalar.nearbyintd(value);
    @compileError("zmath.nearbyint supports only f16, f32, and f64");
}

pub fn rint(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.rintv(value);
    if (T == f16) return scalar.rinth(value);
    if (T == f32) return scalar.rintf(value);
    if (T == f64) return scalar.rintd(value);
    @compileError("zmath.rint supports only f16, f32, and f64");
}

pub fn round(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.ieee_family.roundv(value);
    if (T == f16) return scalar.roundh(value);
    if (T == f32) return scalar.roundf(value);
    if (T == f64) return scalar.roundd(value);
    @compileError("zmath.round supports only f16, f32, and f64");
}

pub fn roundeven(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.roundevenv(value);
    if (T == f16) return scalar.roundevenh(value);
    if (T == f32) return scalar.roundevenf(value);
    if (T == f64) return scalar.roundevend(value);
    @compileError("zmath.roundeven supports only f16, f32, and f64");
}

pub fn trunc(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.truncv(value);
    if (T == f16) return scalar.trunch(value);
    if (T == f32) return scalar.truncf(value);
    if (T == f64) return scalar.truncd(value);
    @compileError("zmath.trunc supports only f16, f32, and f64");
}

pub fn lrint(value: anytype) LaneResult(IntegerResult, @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.lrintv(value);
    if (T == f16) return scalar.lrinth(value);
    if (T == f32) return scalar.lrintf(value);
    if (T == f64) return scalar.lrintd(value);
    @compileError("zmath.lrint supports only f16, f32, and f64");
}

pub fn llrint(value: anytype) LaneResult(IntegerResult, @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.llrintv(value);
    if (T == f16) return scalar.llrinth(value);
    if (T == f32) return scalar.llrintf(value);
    if (T == f64) return scalar.llrintd(value);
    @compileError("zmath.llrint supports only f16, f32, and f64");
}

pub fn lround(value: anytype) LaneResult(IntegerResult, @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.lroundv(value);
    if (T == f16) return scalar.lroundh(value);
    if (T == f32) return scalar.lroundf(value);
    if (T == f64) return scalar.lroundd(value);
    @compileError("zmath.lround supports only f16, f32, and f64");
}

pub fn llround(value: anytype) LaneResult(IntegerResult, @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.llroundv(value);
    if (T == f16) return scalar.llroundh(value);
    if (T == f32) return scalar.llroundf(value);
    if (T == f64) return scalar.llroundd(value);
    @compileError("zmath.llround supports only f16, f32, and f64");
}

pub fn fromfp(value: anytype, direction: FromFpRound, width: u32) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.fromfpv(value, direction, width);
    if (T == f16) return scalar.fromfph(value, direction, width);
    if (T == f32) return scalar.fromfpf(value, direction, width);
    if (T == f64) return scalar.fromfpd(value, direction, width);
    @compileError("zmath.fromfp supports only f16, f32, and f64");
}

pub fn ufromfp(value: anytype, direction: FromFpRound, width: u32) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.ufromfpv(value, direction, width);
    if (T == f16) return scalar.ufromfph(value, direction, width);
    if (T == f32) return scalar.ufromfpf(value, direction, width);
    if (T == f64) return scalar.ufromfpd(value, direction, width);
    @compileError("zmath.ufromfp supports only f16, f32, and f64");
}

pub fn fromfpx(value: anytype, direction: FromFpRound, width: u32) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.fromfpxv(value, direction, width);
    if (T == f16) return scalar.fromfpxh(value, direction, width);
    if (T == f32) return scalar.fromfpxf(value, direction, width);
    if (T == f64) return scalar.fromfpxd(value, direction, width);
    @compileError("zmath.fromfpx supports only f16, f32, and f64");
}

pub fn ufromfpx(value: anytype, direction: FromFpRound, width: u32) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.ufromfpxv(value, direction, width);
    if (T == f16) return scalar.ufromfpxh(value, direction, width);
    if (T == f32) return scalar.ufromfpxf(value, direction, width);
    if (T == f64) return scalar.ufromfpxd(value, direction, width);
    @compileError("zmath.ufromfpx supports only f16, f32, and f64");
}

pub fn atan2(y: anytype, x: anytype) LaneResult(Atan2Result(FloatLane(@TypeOf(y))), @TypeOf(y)) {
    const T = @TypeOf(y);
    if (@typeInfo(T) == .vector) return vector.struct_family.atan2v(y, x);
    if (@TypeOf(x) != T) {
        @compileError("zmath.atan2 requires both arguments to have the same type");
    }
    if (T == f16) return scalar.atan2h(y, x);
    if (T == f32) return scalar.atan2f(y, x);
    if (T == f64) return scalar.atan2d(y, x);
    @compileError("zmath.atan2 supports only f16, f32, and f64");
}

pub fn asinpi(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.asinpiv(value);
    if (T == f16) return scalar.asinpih(value);
    if (T == f32) return scalar.asinpif(value);
    if (T == f64) return scalar.asinpid(value);
    @compileError("zmath.asinpi supports only f16, f32, and f64");
}

pub fn acospi(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.acospiv(value);
    if (T == f16) return scalar.acospih(value);
    if (T == f32) return scalar.acospif(value);
    if (T == f64) return scalar.acospid(value);
    @compileError("zmath.acospi supports only f16, f32, and f64");
}

pub fn atanpi(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.atanpiv(value);
    if (T == f16) return scalar.atanpih(value);
    if (T == f32) return scalar.atanpif(value);
    if (T == f64) return scalar.atanpid(value);
    @compileError("zmath.atanpi supports only f16, f32, and f64");
}

pub fn atan2pi(y: anytype, x: anytype) LaneResult(Atan2Result(FloatLane(@TypeOf(y))), @TypeOf(y)) {
    const T = @TypeOf(y);
    if (@typeInfo(T) == .vector) return vector.struct_family.atan2piv(y, x);
    if (@TypeOf(x) != T) {
        @compileError("zmath.atan2pi requires both arguments to have the same type");
    }
    if (T == f16) return scalar.atan2pih(y, x);
    if (T == f32) return scalar.atan2pif(y, x);
    if (T == f64) return scalar.atan2pid(y, x);
    @compileError("zmath.atan2pi supports only f16, f32, and f64");
}

pub fn pown(value: anytype, exponent: i64) LaneResult(PowerResult(FloatLane(@TypeOf(value))), @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.pownv(value, exponent);
    if (T == f16) return scalar.pownh(value, exponent);
    if (T == f32) return scalar.pownf(value, exponent);
    if (T == f64) return scalar.pownd(value, exponent);
    @compileError("zmath.pown supports only f16, f32, and f64");
}

pub fn powr(base: anytype, exponent: anytype) @TypeOf(base) {
    const T = @TypeOf(base);
    if (@typeInfo(T) == .vector) return vector.powrv(base, exponent);
    if (@TypeOf(exponent) != T) {
        @compileError("zmath.powr requires base and exponent to have the same type");
    }
    if (T == f16) return scalar.powrh(base, exponent);
    if (T == f32) return scalar.powrf(base, exponent);
    if (T == f64) return scalar.powrd(base, exponent);
    @compileError("zmath.powr supports only f16, f32, and f64");
}

pub fn compoundn(value: anytype, exponent: i64) LaneResult(PowerResult(FloatLane(@TypeOf(value))), @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.compoundnv(value, exponent);
    if (T == f16) return scalar.compoundnh(value, exponent);
    if (T == f32) return scalar.compoundnf(value, exponent);
    if (T == f64) return scalar.compoundnd(value, exponent);
    @compileError("zmath.compoundn supports only f16, f32, and f64");
}

pub fn rootn(value: anytype, degree: i64) LaneResult(PowerResult(FloatLane(@TypeOf(value))), @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.rootnv(value, degree);
    if (T == f16) return scalar.rootnh(value, degree);
    if (T == f32) return scalar.rootnf(value, degree);
    if (T == f64) return scalar.rootnd(value, degree);
    @compileError("zmath.rootn supports only f16, f32, and f64");
}

pub fn lgamma(value: anytype) LaneResult(LgammaResult(FloatLane(@TypeOf(value))), @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.lgammav(value);
    if (T == f16) return scalar.lgammah(value);
    if (T == f32) return scalar.lgammaf(value);
    if (T == f64) return scalar.lgammad(value);
    @compileError("zmath.lgamma supports only f16, f32, and f64");
}

pub fn fmod(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.fmodv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fmod requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fmodh(x, y);
    if (T == f32) return scalar.fmodf(x, y);
    if (T == f64) return scalar.fmodd(x, y);
    @compileError("zmath.fmod supports only f16, f32, and f64");
}

pub fn remainder(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.remainderv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.remainder requires both arguments to have the same type");
    }
    if (T == f16) return scalar.remainderh(x, y);
    if (T == f32) return scalar.remainderf(x, y);
    if (T == f64) return scalar.remainderd(x, y);
    @compileError("zmath.remainder supports only f16, f32, and f64");
}

pub fn remquo(x: anytype, y: anytype) LaneResult(RemquoResult(FloatLane(@TypeOf(x))), @TypeOf(x)) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.struct_family.remquov(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.remquo requires both arguments to have the same type");
    }
    if (T == f16) return scalar.remquoh(x, y);
    if (T == f32) return scalar.remquof(x, y);
    if (T == f64) return scalar.remquod(x, y);
    @compileError("zmath.remquo supports only f16, f32, and f64");
}

pub fn frexp(value: anytype) LaneResult(FrexpResult(FloatLane(@TypeOf(value))), @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.frexpv(value);
    if (T == f16) return scalar.frexph(value);
    if (T == f32) return scalar.frexpf(value);
    if (T == f64) return scalar.frexpd(value);
    @compileError("zmath.frexp supports only f16, f32, and f64");
}

pub fn modf(value: anytype) LaneResult(ModfResult(FloatLane(@TypeOf(value))), @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.modfv(value);
    if (T == f16) return scalar.modfh(value);
    if (T == f32) return scalar.modff(value);
    if (T == f64) return scalar.modfd(value);
    @compileError("zmath.modf supports only f16, f32, and f64");
}

pub fn ilogb(value: anytype) LaneVector(i32, @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.ieee_family.ilogbv(value);
    if (T == f16) return scalar.ilogbh(value);
    if (T == f32) return scalar.ilogbf(value);
    if (T == f64) return scalar.ilogbd(value);
    @compileError("zmath.ilogb supports only f16, f32, and f64");
}

pub fn llogb(value: anytype) LaneVector(i64, @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.ieee_family.llogbv(value);
    if (T == f16) return scalar.llogbh(value);
    if (T == f32) return scalar.llogbf(value);
    if (T == f64) return scalar.llogbd(value);
    @compileError("zmath.llogb supports only f16, f32, and f64");
}

pub fn logb(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.logbv(value);
    if (T == f16) return scalar.logbh(value);
    if (T == f32) return scalar.logbf(value);
    if (T == f64) return scalar.logbd(value);
    @compileError("zmath.logb supports only f16, f32, and f64");
}

pub fn ldexp(value: anytype, exponent: i32) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.ldexpv(value, exponent);
    if (T == f16) return scalar.ldexph(value, exponent);
    if (T == f32) return scalar.ldexpf(value, exponent);
    if (T == f64) return scalar.ldexpd(value, exponent);
    @compileError("zmath.ldexp supports only f16, f32, and f64");
}

pub fn scalbn(value: anytype, exponent: i32) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.scalbnv(value, exponent);
    if (T == f16) return scalar.scalbnh(value, exponent);
    if (T == f32) return scalar.scalbnf(value, exponent);
    if (T == f64) return scalar.scalbnd(value, exponent);
    @compileError("zmath.scalbn supports only f16, f32, and f64");
}

pub fn scalbln(value: anytype, exponent: i64) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.scalblnv(value, exponent);
    if (T == f16) return scalar.scalblnh(value, exponent);
    if (T == f32) return scalar.scalblnf(value, exponent);
    if (T == f64) return scalar.scalblnd(value, exponent);
    @compileError("zmath.scalbln supports only f16, f32, and f64");
}

pub fn fdim(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.fdimv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fdim requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fdimh(x, y);
    if (T == f32) return scalar.fdimf(x, y);
    if (T == f64) return scalar.fdimd(x, y);
    @compileError("zmath.fdim supports only f16, f32, and f64");
}

pub fn fmax(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.fmaxv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fmax requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fmaxh(x, y);
    if (T == f32) return scalar.fmaxf(x, y);
    if (T == f64) return scalar.fmaxd(x, y);
    @compileError("zmath.fmax supports only f16, f32, and f64");
}

pub fn fmin(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.fminv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fmin requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fminh(x, y);
    if (T == f32) return scalar.fminf(x, y);
    if (T == f64) return scalar.fmind(x, y);
    @compileError("zmath.fmin supports only f16, f32, and f64");
}

pub fn fmaximum(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.fmaximumv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fmaximum requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fmaximumh(x, y);
    if (T == f32) return scalar.fmaximumf(x, y);
    if (T == f64) return scalar.fmaximumd(x, y);
    @compileError("zmath.fmaximum supports only f16, f32, and f64");
}

pub fn fminimum(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.fminimumv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fminimum requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fminimumh(x, y);
    if (T == f32) return scalar.fminimumf(x, y);
    if (T == f64) return scalar.fminimumd(x, y);
    @compileError("zmath.fminimum supports only f16, f32, and f64");
}

fn fmaximumMag(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.ieee_family.fmaximumMagv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fmaximum_mag requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fmaximum_magh(x, y);
    if (T == f32) return scalar.fmaximum_magf(x, y);
    if (T == f64) return scalar.fmaximum_magd(x, y);
    @compileError("zmath.fmaximum_mag supports only f16, f32, and f64");
}

pub const fmaximum_mag = fmaximumMag;

fn fminimumMag(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.ieee_family.fminimumMagv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fminimum_mag requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fminimum_magh(x, y);
    if (T == f32) return scalar.fminimum_magf(x, y);
    if (T == f64) return scalar.fminimum_magd(x, y);
    @compileError("zmath.fminimum_mag supports only f16, f32, and f64");
}

pub const fminimum_mag = fminimumMag;

fn fmaximumNum(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.ieee_family.fmaximumNumv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fmaximum_num requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fmaximum_numh(x, y);
    if (T == f32) return scalar.fmaximum_numf(x, y);
    if (T == f64) return scalar.fmaximum_numd(x, y);
    @compileError("zmath.fmaximum_num supports only f16, f32, and f64");
}

pub const fmaximum_num = fmaximumNum;

fn fminimumNum(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.ieee_family.fminimumNumv(x, y);
    if (@TypeOf(y) != T) {
        @compileError("zmath.fminimum_num requires both arguments to have the same type");
    }
    if (T == f16) return scalar.fminimum_numh(x, y);
    if (T == f32) return scalar.fminimum_numf(x, y);
    if (T == f64) return scalar.fminimum_numd(x, y);
    @compileError("zmath.fminimum_num supports only f16, f32, and f64");
}

pub const fminimum_num = fminimumNum;

fn fmaximumMagNum(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.ieee_family.fmaximumMagNumv(x, y);
    if (@TypeOf(y) != T) {
        @compileError(
            "zmath.fmaximum_mag_num requires both arguments to have the same type",
        );
    }
    if (T == f16) return scalar.fmaximum_mag_numh(x, y);
    if (T == f32) return scalar.fmaximum_mag_numf(x, y);
    if (T == f64) return scalar.fmaximum_mag_numd(x, y);
    @compileError("zmath.fmaximum_mag_num supports only f16, f32, and f64");
}

pub const fmaximum_mag_num = fmaximumMagNum;

fn fminimumMagNum(x: anytype, y: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.ieee_family.fminimumMagNumv(x, y);
    if (@TypeOf(y) != T) {
        @compileError(
            "zmath.fminimum_mag_num requires both arguments to have the same type",
        );
    }
    if (T == f16) return scalar.fminimum_mag_numh(x, y);
    if (T == f32) return scalar.fminimum_mag_numf(x, y);
    if (T == f64) return scalar.fminimum_mag_numd(x, y);
    @compileError("zmath.fminimum_mag_num supports only f16, f32, and f64");
}

pub const fminimum_mag_num = fminimumMagNum;

pub fn canonicalize(value: anytype) LaneResult(CanonicalizeResult(FloatLane(@TypeOf(value))), @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.canonicalizev(value);
    if (T == f16) return ieee.canonicalizeh(value);
    if (T == f32) return ieee.canonicalizef(value);
    if (T == f64) return ieee.canonicalized(value);
    @compileError("zmath.canonicalize supports only f16, f32, and f64");
}

pub fn copySign(value: anytype, sign: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.ieee_family.copySignv(value, sign);
    if (@TypeOf(sign) != T) {
        @compileError("zmath.copySign requires both arguments to have the same type");
    }
    if (T == f16) return ieee.copySignh(value, sign);
    if (T == f32) return ieee.copySignf(value, sign);
    if (T == f64) return ieee.copySignd(value, sign);
    @compileError("zmath.copySign supports only f16, f32, and f64");
}

pub fn nan(comptime T: type, tag: []const u8) T {
    if (T == f16) return ieee.nanh(tag);
    if (T == f32) return ieee.nanf(tag);
    if (T == f64) return ieee.nand(tag);
    @compileError("zmath.nan supports only f16, f32, and f64");
}

pub fn nextAfter(from: anytype, to: anytype) @TypeOf(from) {
    const T = @TypeOf(from);
    if (@typeInfo(T) == .vector) return vector.ieee_family.nextAfterv(from, to);
    if (@TypeOf(to) != T) {
        @compileError("zmath.nextAfter requires both arguments to have the same type");
    }
    if (T == f16) return ieee.nextAfterh(from, to);
    if (T == f32) return ieee.nextAfterf(from, to);
    if (T == f64) return ieee.nextAfterd(from, to);
    @compileError("zmath.nextAfter supports only f16, f32, and f64");
}

pub fn nextToward(value: anytype, toward: f64) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.intarg.nextTowardv(value, toward);
    if (T == f16) return ieee.nextTowardh(value, toward);
    if (T == f32) return ieee.nextTowardf(value, toward);
    if (T == f64) return ieee.nextTowardd(value, toward);
    @compileError("zmath.nextToward supports only f16, f32, and f64");
}

pub fn nextUp(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.ieee_family.nextUpv(value);
    if (T == f16) return ieee.nextUph(value);
    if (T == f32) return ieee.nextUpf(value);
    if (T == f64) return ieee.nextUpd(value);
    @compileError("zmath.nextUp supports only f16, f32, and f64");
}

pub fn nextDown(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.ieee_family.nextDownv(value);
    if (T == f16) return ieee.nextDownh(value);
    if (T == f32) return ieee.nextDownf(value);
    if (T == f64) return ieee.nextDownd(value);
    @compileError("zmath.nextDown supports only f16, f32, and f64");
}

pub fn totalOrder(left: anytype, right: anytype) LaneVector(bool, @TypeOf(left)) {
    const T = @TypeOf(left);
    if (@typeInfo(T) == .vector) return vector.ieee_family.totalOrderv(left, right);
    if (@TypeOf(right) != T) {
        @compileError("zmath.totalOrder requires both arguments to have the same type");
    }
    if (T == f16) return ieee.totalOrderh(left, right);
    if (T == f32) return ieee.totalOrderf(left, right);
    if (T == f64) return ieee.totalOrderd(left, right);
    @compileError("zmath.totalOrder supports only f16, f32, and f64");
}

pub fn totalOrderMag(left: anytype, right: anytype) LaneVector(bool, @TypeOf(left)) {
    const T = @TypeOf(left);
    if (@typeInfo(T) == .vector) return vector.ieee_family.totalOrderMagv(left, right);
    if (@TypeOf(right) != T) {
        @compileError("zmath.totalOrderMag requires both arguments to have the same type");
    }
    if (T == f16) return ieee.totalOrderMagh(left, right);
    if (T == f32) return ieee.totalOrderMagf(left, right);
    if (T == f64) return ieee.totalOrderMagd(left, right);
    @compileError("zmath.totalOrderMag supports only f16, f32, and f64");
}

pub fn getPayload(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.ieee_family.getPayloadv(value);
    if (T == f16) return ieee.getPayloadh(value);
    if (T == f32) return ieee.getPayloadf(value);
    if (T == f64) return ieee.getPayloadd(value);
    @compileError("zmath.getPayload supports only f16, f32, and f64");
}

pub fn setPayload(value: anytype) LaneResult(PayloadResult(FloatLane(@TypeOf(value))), @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.setPayloadv(value);
    if (T == f16) return ieee.setPayloadh(value);
    if (T == f32) return ieee.setPayloadf(value);
    if (T == f64) return ieee.setPayloadd(value);
    @compileError("zmath.setPayload supports only f16, f32, and f64");
}

pub fn setPayloadSig(value: anytype) LaneResult(PayloadResult(FloatLane(@TypeOf(value))), @TypeOf(value)) {
    const T = @TypeOf(value);
    if (@typeInfo(T) == .vector) return vector.struct_family.setPayloadSigv(value);
    if (T == f16) return ieee.setPayloadSigh(value);
    if (T == f32) return ieee.setPayloadSigf(value);
    if (T == f64) return ieee.setPayloadSigd(value);
    @compileError("zmath.setPayloadSig supports only f16, f32, and f64");
}

pub fn fadd(x: f64, y: f64) f32 {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.narrow.faddv(x, y);
    return scalar_narrow.fadd(x, y);
}

pub fn fsub(x: f64, y: f64) f32 {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.narrow.fsubv(x, y);
    return scalar_narrow.fsub(x, y);
}

pub fn fmul(x: f64, y: f64) f32 {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.narrow.fmulv(x, y);
    return scalar_narrow.fmul(x, y);
}

pub fn fdiv(x: f64, y: f64) f32 {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.narrow.fdivv(x, y);
    return scalar_narrow.fdiv(x, y);
}

pub fn ffma(x: f64, y: f64, z: f64) f32 {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.narrow.ffmav(x, y, z);
    return scalar_narrow.ffma(x, y, z);
}

pub fn fsqrt(x: f64) f32 {
    const T = @TypeOf(x);
    if (@typeInfo(T) == .vector) return vector.narrow.fsqrtv(x);
    return scalar_narrow.fsqrt(x);
}

test {
    _ = @import("ieee_test.zig");
    _ = @import("ieee_edge_test.zig");
    _ = @import("scalar_test.zig");
    _ = @import("scalar_fmah_test.zig");
    _ = @import("scalar_fmaf_test.zig");
    _ = @import("scalar_fmad_test.zig");
    _ = @import("scalar_erfc_test.zig");
    _ = @import("scalar_log10_test.zig");
    _ = @import("scalar_log10p1_test.zig");
    _ = @import("scalar_cosh_test.zig");
    _ = @import("scalar_boundary_test.zig");
    _ = @import("scalar_rsqrt_test.zig");
    _ = @import("scalar_sqrt_test.zig");
    _ = @import("api_suffix_test.zig");
    _ = @import("ieee_h_complete_test.zig");
    _ = @import("ieee_exhaustive_h_test.zig");
    _ = @import("ieee_scale_widths_test.zig");
    _ = @import("scalar_exp_test.zig");
    _ = @import("scalar_expm1_test.zig");
    _ = @import("scalar_exp2m1_test.zig");
    _ = @import("scalar_exp10m1_test.zig");
    _ = @import("scalar_exp10_test.zig");
    _ = @import("scalar_log_test.zig");
    _ = @import("scalar_log2_test.zig");
    _ = @import("scalar_log1p_test.zig");
    _ = @import("scalar_log2p1_test.zig");
    _ = @import("scalar_atanh_test.zig");
    _ = @import("scalar_sin_test.zig");
    _ = @import("scalar_acosh_test.zig");
    _ = @import("scalar_acos_test.zig");
    _ = @import("scalar_cos_test.zig");
    _ = @import("scalar_cospi_test.zig");
    _ = @import("scalar_sind_test.zig");
    _ = @import("scalar_pow_test.zig");
    _ = @import("root_pow_test.zig");
    _ = @import("scalar_tan_test.zig");
    _ = @import("scalar_tanpi_test.zig");
    _ = @import("scalar_stdlib_agreement_test.zig");
    _ = @import("scalar_exp2_test.zig");
    _ = @import("scalar_cbrt_test.zig");
    _ = @import("scalar_hypot_test.zig");
    _ = @import("scalar_asin_test.zig");
    _ = @import("scalar_asinh_test.zig");
    _ = @import("scalar_atan_test.zig");
    _ = @import("scalar_sinh_test.zig");
    _ = @import("scalar_erf_test.zig");
    _ = @import("scalar_tanh_test.zig");
    _ = @import("scalar_tgamma_test.zig");
    _ = @import("scalar_tgamma_contract_test.zig");
    _ = @import("scalar_sinpi_test.zig");
    _ = @import("root_sinpi_test.zig");
    _ = @import("scalar_c23_common_test.zig");
    _ = @import("scalar_c23_types_test.zig");
    _ = @import("scalar_rounding_test.zig");
    _ = @import("scalar_integer_rounding_test.zig");
    _ = @import("scalar_fromfp_test.zig");
    _ = @import("scalar_remainder_test.zig");
    _ = @import("scalar_decompose_test.zig");
    _ = @import("scalar_scaling_test.zig");
    _ = @import("scalar_extrema_test.zig");
    _ = @import("scalar_lgamma_test.zig");
    _ = @import("scalar_atan2_test.zig");
    _ = @import("scalar_pi_inverse_test.zig");
    _ = @import("scalar_power_c23_test.zig");
    _ = @import("scalar_logp1_alias_test.zig");
    _ = @import("scalar_narrow_test.zig");
    _ = @import("ieee_c23_core_test.zig");
    _ = @import("ieee_c23_optional_test.zig");
    _ = @import("c23_public_surface_test.zig");
    _ = @import("vector_test.zig");
    _ = @import("vector_surface_test.zig");
    _ = @import("vector_value_test.zig");
    _ = @import("complex_test.zig");
    _ = @import("quaternion_test.zig");
    _ = @import("algebra_test.zig");
    _ = @import("matrix_test.zig");
    _ = @import("geometry_test.zig");
    _ = @import("linalg_test.zig");
}

test "root exports nearest-even Atan H, F, and D" {
    try std.testing.expectEqual(
        @as(u16, 0x3a48),
        @as(u16, @bitCast(scalar.atanh(@as(f16, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f49_0fdb),
        @as(u32, @bitCast(scalar.atanf(@as(f32, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3fe9_21fb_5444_2d18),
        @as(u64, @bitCast(scalar.atand(@as(f64, 1.0)))),
    );
}

test "root exports nearest-even Atanh family" {
    try std.testing.expectEqual(
        @as(u16, 0x3865),
        @as(u16, @bitCast(atanh(@as(f16, 0.5)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f0c_9f54),
        @as(u32, @bitCast(atanh(@as(f32, 0.5)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3fe1_93ea_7aad_030b),
        @as(u64, @bitCast(atanh(@as(f64, 0.5)))),
    );
}

test "root exports nearest-even Asin family" {
    try std.testing.expectEqual(
        @as(u16, 0x3e48),
        @as(u16, @bitCast(asin(@as(f16, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3fc9_0fdb),
        @as(u32, @bitCast(asin(@as(f32, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff9_21fb_5444_2d18),
        @as(u64, @bitCast(asin(@as(f64, 1.0)))),
    );
}
test "root exports nearest-even Acos family" {
    try std.testing.expectEqual(
        @as(u16, 0x3c30),
        @as(u16, @bitCast(acos(@as(f16, 0.5)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f86_0a92),
        @as(u32, @bitCast(acos(@as(f32, 0.5)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff0_c152_382d_7366),
        @as(u64, @bitCast(acos(@as(f64, 0.5)))),
    );
}

test "root exports nearest-even Asinh family" {
    try std.testing.expectEqual(
        @as(u16, 0x3b0d),
        @as(u16, @bitCast(asinh(@as(f16, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f61_a1b3),
        @as(u32, @bitCast(asinh(@as(f32, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3fec_3436_6179_d427),
        @as(u64, @bitCast(asinh(@as(f64, 1.0)))),
    );
}

test "generic scalar dispatch covers all formats" {
    try generic_scalar_test.run(@This());
}

test "C23 public root surface" {
    try c23_public_surface_test.run(@This());
}

test "root exports nearest-even binary64 sine" {
    const input: f64 = 1.0;
    try std.testing.expectEqual(@as(u64, 0x3fea_ed54_8f09_0cee), @as(u64, @bitCast(scalar.sind(input))));
}
test "root exports nearest-even binary64 cosine" {
    const input: f64 = 1.0;
    try std.testing.expectEqual(@as(u64, 0x3fe1_4a28_0fb5_068c), @as(u64, @bitCast(scalar.cosd(input))));
}
test "root exports nearest-even binary16/32/64 tangent" {
    try std.testing.expectEqual(@as(u16, 0x3e3b), @as(u16, @bitCast(tan(@as(f16, 1.0)))));
    try std.testing.expectEqual(@as(u32, 0x3fc7_5923), @as(u32, @bitCast(tan(@as(f32, 1.0)))));
    try std.testing.expectEqual(@as(u64, 0x3ff8_eb24_5cbe_e3a6), @as(u64, @bitCast(tan(@as(f64, 1.0)))));
}

test "root exports nearest-even Tanpi family" {
    try std.testing.expectEqual(
        @as(u16, 0x3c00),
        @as(u16, @bitCast(tanpi(@as(f16, 0.25)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f80_0000),
        @as(u32, @bitCast(tanpi(@as(f32, 0.25)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff0_0000_0000_0000),
        @as(u64, @bitCast(tanpi(@as(f64, 0.25)))),
    );
}

test "root exports nearest-even hyperbolic Sinh" {
    try std.testing.expectEqual(
        @as(u16, 0x3cb3),
        @as(u16, @bitCast(sinh(@as(f16, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f96_6cfe),
        @as(u32, @bitCast(sinh(@as(f32, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff2_cd9f_c44e_b982),
        @as(u64, @bitCast(sinh(@as(f64, 1.0)))),
    );
}

test "root exports nearest-even hyperbolic Cosh family" {
    try std.testing.expectEqual(
        @as(u16, 0x3c83),
        @as(u16, @bitCast(cosh(@as(f16, 0.5)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f90_560c),
        @as(u32, @bitCast(cosh(@as(f32, 0.5)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff2_0ac1_862a_e8d0),
        @as(u64, @bitCast(cosh(@as(f64, 0.5)))),
    );
}
test "root exports nearest-even inverse hyperbolic Acosh family" {
    try std.testing.expectEqual(
        @as(u16, 0x3bb3),
        @as(u16, @bitCast(acosh(@as(f16, 1.5)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f76_6165),
        @as(u32, @bitCast(acosh(@as(f32, 1.5)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3fee_cc2c_aec5_160a),
        @as(u64, @bitCast(acosh(@as(f64, 1.5)))),
    );
}

test "root exports nearest-even binary64 hyperbolic tangent" {
    const input: f64 = 1.0;
    try std.testing.expectEqual(
        @as(u64, 0x3fe8_5efa_b514_f394),
        @as(u64, @bitCast(tanh(input))),
    );
}
test "root exports nearest-even Exp2 family" {
    try std.testing.expectEqual(
        @as(u16, 0x4000),
        @as(u16, @bitCast(scalar.exp2h(1.0))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x4000_0000),
        @as(u32, @bitCast(scalar.exp2f(1.0))),
    );

    try std.testing.expectEqual(
        @as(u64, 0x4000_0000_0000_0000),
        @as(u64, @bitCast(scalar.exp2d(1.0))),
    );
}

test "root exports nearest-even Log2 family" {
    try std.testing.expectEqual(@as(u16, 0x3c00), @as(u16, @bitCast(log2(@as(f16, 2.0)))));
    try std.testing.expectEqual(@as(u32, 0x3f80_0000), @as(u32, @bitCast(log2(@as(f32, 2.0)))));
    try std.testing.expectEqual(
        @as(u64, 0x3ff0_0000_0000_0000),
        @as(u64, @bitCast(log2(@as(f64, 2.0)))),
    );
}
test "root exports nearest-even Rsqrt family" {
    try std.testing.expectEqual(
        @as(u16, 0x39a8),
        @as(u16, @bitCast(rsqrt(@as(f16, 2.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f35_04f3),
        @as(u32, @bitCast(rsqrt(@as(f32, 2.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3fe6_a09e_667f_3bcd),
        @as(u64, @bitCast(rsqrt(@as(f64, 2.0)))),
    );
}

test "root exports nearest-even Sqrt family" {
    try std.testing.expectEqual(@as(u16, 0x3da8), @as(u16, @bitCast(sqrt(@as(f16, 2.0)))));
    try std.testing.expectEqual(@as(u32, 0x3fb5_04f3), @as(u32, @bitCast(sqrt(@as(f32, 2.0)))));
    try std.testing.expectEqual(@as(u64, 0x3ff6_a09e_667f_3bcd), @as(u64, @bitCast(sqrt(@as(f64, 2.0)))));
}
test "root exports nearest-even Cbrt H, F, and D" {
    try std.testing.expectEqual(
        @as(u16, 0x3d0a),
        @as(u16, @bitCast(scalar.cbrth(@as(f16, 2.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3fa1_4518),
        @as(u32, @bitCast(scalar.cbrtf(@as(f32, 2.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff4_28a2_f98d_728b),
        @as(u64, @bitCast(scalar.cbrtd(@as(f64, 2.0)))),
    );
}

test "root exports nearest-even Log1p H, F, and D" {
    try std.testing.expectEqual(@as(u16, 0x398c), @as(u16, @bitCast(scalar.log1ph(@as(f16, 1.0)))));
    try std.testing.expectEqual(@as(u32, 0x3f31_7218), @as(u32, @bitCast(scalar.log1pf(@as(f32, 1.0)))));
    try std.testing.expectEqual(
        @as(u64, 0x3fe6_2e42_fefa_39ef),
        @as(u64, @bitCast(scalar.log1pd(@as(f64, 1.0)))),
    );
}

test "root exports nearest-even Hypot H, F, and D" {
    try std.testing.expectEqual(
        @as(u16, 0x4500),
        @as(u16, @bitCast(scalar.hypoth(@as(f16, 3.0), @as(f16, 4.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x40a0_0000),
        @as(u32, @bitCast(scalar.hypotf(3.0, 4.0))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x4014_0000_0000_0000),
        @as(u64, @bitCast(scalar.hypotd(3.0, 4.0))),
    );
}

test "root exports nearest-even Exp10m1 family" {
    try std.testing.expectEqual(
        @as(u16, 0x4880),
        @as(u16, @bitCast(exp10m1(@as(f16, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x4110_0000),
        @as(u32, @bitCast(exp10m1(@as(f32, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x4022_0000_0000_0000),
        @as(u64, @bitCast(exp10m1(@as(f64, 1.0)))),
    );
}

test "root exports nearest-even Exp2m1 family" {
    try std.testing.expectEqual(
        @as(u16, 0x3c00),
        @as(u16, @bitCast(exp2m1(@as(f16, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f80_0000),
        @as(u32, @bitCast(exp2m1(@as(f32, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff0_0000_0000_0000),
        @as(u64, @bitCast(exp2m1(@as(f64, 1.0)))),
    );
}

test "root exports nearest-even Log10 family" {
    try std.testing.expectEqual(
        @as(u16, 0x3c00),
        @as(u16, @bitCast(log10(@as(f16, 10.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f80_0000),
        @as(u32, @bitCast(log10(@as(f32, 10.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff0_0000_0000_0000),
        @as(u64, @bitCast(log10(@as(f64, 10.0)))),
    );
}

test "root exports nearest-even Log10p1 family" {
    try std.testing.expectEqual(
        @as(u16, 0x3c00),
        @as(u16, @bitCast(log10p1(@as(f16, 9.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f80_0000),
        @as(u32, @bitCast(log10p1(@as(f32, 9.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ff0_0000_0000_0000),
        @as(u64, @bitCast(log10p1(@as(f64, 9.0)))),
    );
}

test "root exports nearest-even Exp10 family" {
    try std.testing.expectEqual(
        @as(u16, 0x4900),
        @as(u16, @bitCast(exp10(@as(f16, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x4120_0000),
        @as(u32, @bitCast(exp10(@as(f32, 1.0)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x4024_0000_0000_0000),
        @as(u64, @bitCast(exp10(@as(f64, 1.0)))),
    );
}
