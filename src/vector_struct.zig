//! Vector (`@Vector`) overloads for the elementwise operations whose scalar
//! kernel returns a result struct.
//!
//! Every operation here applies the *identical* accepted scalar kernel to each
//! lane, so lane `i` of the result is bit-for-bit equal to the scalar result for
//! lane `i` of the input. That preserves the schema-3 `max_ulp = 0` contract and
//! the immutable corpus without new acceptance evidence: the vector path is a
//! per-lane application of the same scalar code, never a re-derived polynomial.
//!
//! `vector.zig` covers the float-to-float operations, where a lane maps to a
//! lane. This module covers the operations where a lane maps to a RESULT STRUCT:
//! `IntegerResult`, `LgammaResult(T)`, `FrexpResult(T)`, `ModfResult(T)`,
//! `Atan2Result(T)`, `RemquoResult(T)`, `PayloadResult(T)`, and
//! `CanonicalizeResult(T)`.
//!
//! The result is a lane-ordered ARRAY `[N]Result`, NOT a `@Vector(N, Result)`.
//! Zig forbids a struct as a vector element type and reports "expected integer,
//! float, bool, or pointer for the vector element type". Every other vector
//! operation in this library returns a real `@Vector`; these thirteen cannot.
//! That is a language limitation, not an oversight, and it is why the two
//! mapping helpers below are private instead of the `@Vector`-returning helpers
//! in `vector_core.zig`.
//!
//! Element `i` of the array is the value the scalar kernel itself returned for
//! lane `i`, moved as a whole, in lane order. The per-lane bit-identity
//! guarantee is therefore unchanged, and the `max_ulp = 0` contract and the
//! immutable corpus still hold with no new acceptance evidence. Note what this
//! shape does NOT buy: the kernel runs once per lane, exactly as a scalar loop
//! would, so there is no throughput advantage over a scalar loop and no SIMD
//! acceleration. The value here is a uniform lane-parallel API, nothing more.
//!
//! Those structs are bare `struct` values with auto layout. They are never
//! `extern` and never `packed`, so `@bitCast` and `@offsetOf` do not apply to
//! them. Every field was assigned by the kernel, and nothing in this module
//! reinterprets a result or reads a field through an assumed offset.
//!
//! Each operation provides the generic entry point `OPv`, which selects the
//! kernel from the lane type, plus `OPhv`, `OPfv`, and `OPdv`, which pin the
//! lane type to f16, f32, and f64. The lane count is read from the argument in
//! every form, so one entry point serves every width. None of these names
//! collides with the generic hyperbolic entry points `sinhv`, `coshv`, `tanhv`,
//! `asinhv`, `acoshv`, and `atanhv`, so every operation here carries all four.
//!
//! The kernel for each format is the one bound by the top-level dispatch in
//! `root.zig`. `canonicalize`, `setPayload`, and `setPayloadSig` come from the
//! `ieee` module; the rest come from `scalar`.

const std = @import("std");
const ieee = @import("ieee.zig");
const scalar = @import("scalar.zig");
const types = @import("scalar_c23_types.zig");

// ---------------------------------------------------------------------------
// Per-lane application of a result struct
// ---------------------------------------------------------------------------

/// Returns the lane count of the vector type `X`.
fn widthOf(comptime X: type) usize {
    return switch (@typeInfo(X)) {
        .vector => |v| v.len,
        else => @compileError("zmath vector operations require a @Vector argument"),
    };
}

/// The `[N]R` lane-ordered result array for the lane count of the vector type `X`.
fn ResultArray(comptime R: type, comptime X: anytype) type {
    return [widthOf(X)]R;
}

/// Applies the result-returning kernel `f` to every lane of `x`.
///
/// Element `i` of the result is the `R` that `f` itself returned for lane `i`.
fn mapUnaryStruct(
    comptime T: type,
    comptime R: type,
    comptime f: fn (T) R,
    x: anytype,
) ResultArray(R, @TypeOf(x)) {
    const n = comptime widthOf(@TypeOf(x));
    if (comptime std.meta.Child(@TypeOf(x)) != T) {
        @compileError("zmath vector lane type does not match the selected scalar kernel");
    }
    var out: ResultArray(R, @TypeOf(x)) = undefined;
    inline for (0..n) |i| out[i] = f(x[i]);
    return out;
}

/// Applies the result-returning kernel `f` to corresponding lanes of `a` and `b`.
///
/// Element `i` of the result is the `R` that `f` itself returned for lanes `i`.
fn mapBinaryStruct(
    comptime T: type,
    comptime R: type,
    comptime f: fn (T, T) R,
    a: anytype,
    b: anytype,
) ResultArray(R, @TypeOf(a)) {
    const n = comptime widthOf(@TypeOf(a));
    if (comptime std.meta.Child(@TypeOf(a)) != T) {
        @compileError("zmath vector lane type does not match the selected scalar kernel");
    }
    if (@TypeOf(b) != @TypeOf(a)) {
        @compileError("zmath vector binary operations require both arguments to have the same type");
    }
    var out: ResultArray(R, @TypeOf(a)) = undefined;
    inline for (0..n) |i| out[i] = f(a[i], b[i]);
    return out;
}

// ---------------------------------------------------------------------------
// Integer rounding, one `IntegerResult` per lane
// ---------------------------------------------------------------------------

/// Computes the lrint of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds the nearest `i64` to lane `i`, ties to even, and whether
/// that value is representable. The result is a lane-ordered `[N]` array, not a
/// `@Vector`.
pub fn lrintv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.IntegerResult, scalar.lrinth, x);
    if (T == f32) return mapUnaryStruct(f32, types.IntegerResult, scalar.lrintf, x);
    if (T == f64) return mapUnaryStruct(f64, types.IntegerResult, scalar.lrintd, x);
    @compileError("zmath.lrintv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the lrint of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn lrinthv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f16, types.IntegerResult, scalar.lrinth, x);
}

/// Computes the lrint of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn lrintfv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f32, types.IntegerResult, scalar.lrintf, x);
}

/// Computes the lrint of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn lrintdv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f64, types.IntegerResult, scalar.lrintd, x);
}

/// Computes the llrint of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds the nearest `i64` to lane `i`, ties to even, and whether
/// that value is representable. The result is a lane-ordered `[N]` array, not a
/// `@Vector`.
pub fn llrintv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.IntegerResult, scalar.llrinth, x);
    if (T == f32) return mapUnaryStruct(f32, types.IntegerResult, scalar.llrintf, x);
    if (T == f64) return mapUnaryStruct(f64, types.IntegerResult, scalar.llrintd, x);
    @compileError("zmath.llrintv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the llrint of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn llrinthv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f16, types.IntegerResult, scalar.llrinth, x);
}

/// Computes the llrint of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn llrintfv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f32, types.IntegerResult, scalar.llrintf, x);
}

/// Computes the llrint of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn llrintdv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f64, types.IntegerResult, scalar.llrintd, x);
}

/// Computes the lround of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds the nearest `i64` to lane `i`, ties away from zero, and
/// whether that value is representable. The result is a lane-ordered `[N]`
/// array, not a `@Vector`.
pub fn lroundv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.IntegerResult, scalar.lroundh, x);
    if (T == f32) return mapUnaryStruct(f32, types.IntegerResult, scalar.lroundf, x);
    if (T == f64) return mapUnaryStruct(f64, types.IntegerResult, scalar.lroundd, x);
    @compileError("zmath.lroundv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the lround of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn lroundhv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f16, types.IntegerResult, scalar.lroundh, x);
}

/// Computes the lround of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn lroundfv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f32, types.IntegerResult, scalar.lroundf, x);
}

/// Computes the lround of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn lrounddv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f64, types.IntegerResult, scalar.lroundd, x);
}

/// Computes the llround of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds the nearest `i64` to lane `i`, ties away from zero, and
/// whether that value is representable. The result is a lane-ordered `[N]`
/// array, not a `@Vector`.
pub fn llroundv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.IntegerResult, scalar.llroundh, x);
    if (T == f32) return mapUnaryStruct(f32, types.IntegerResult, scalar.llroundf, x);
    if (T == f64) return mapUnaryStruct(f64, types.IntegerResult, scalar.llroundd, x);
    @compileError("zmath.llroundv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the llround of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn llroundhv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f16, types.IntegerResult, scalar.llroundh, x);
}

/// Computes the llround of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn llroundfv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f32, types.IntegerResult, scalar.llroundf, x);
}

/// Computes the llround of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn llrounddv(x: anytype) ResultArray(types.IntegerResult, @TypeOf(x)) {
    return mapUnaryStruct(f64, types.IntegerResult, scalar.llroundd, x);
}

// ---------------------------------------------------------------------------
// lgamma, one `LgammaResult(T)` per lane
// ---------------------------------------------------------------------------

/// Computes the lgamma of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds the rounded logarithm of the absolute value at lane `i` and
/// the sign of gamma there. The result is a lane-ordered `[N]` array, not a
/// `@Vector`.
pub fn lgammav(x: anytype) ResultArray(types.LgammaResult(std.meta.Child(@TypeOf(x))), @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.LgammaResult(f16), scalar.lgammah, x);
    if (T == f32) return mapUnaryStruct(f32, types.LgammaResult(f32), scalar.lgammaf, x);
    if (T == f64) return mapUnaryStruct(f64, types.LgammaResult(f64), scalar.lgammad, x);
    @compileError("zmath.lgammav supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the lgamma of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn lgammahv(x: anytype) ResultArray(types.LgammaResult(f16), @TypeOf(x)) {
    return mapUnaryStruct(f16, types.LgammaResult(f16), scalar.lgammah, x);
}

/// Computes the lgamma of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn lgammafv(x: anytype) ResultArray(types.LgammaResult(f32), @TypeOf(x)) {
    return mapUnaryStruct(f32, types.LgammaResult(f32), scalar.lgammaf, x);
}

/// Computes the lgamma of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn lgammadv(x: anytype) ResultArray(types.LgammaResult(f64), @TypeOf(x)) {
    return mapUnaryStruct(f64, types.LgammaResult(f64), scalar.lgammad, x);
}

// ---------------------------------------------------------------------------
// Decomposition, one `FrexpResult(T)` or `ModfResult(T)` per lane
// ---------------------------------------------------------------------------

/// Computes the frexp of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds a signed fraction and a binary exponent for lane `i`. The
/// result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn frexpv(x: anytype) ResultArray(types.FrexpResult(std.meta.Child(@TypeOf(x))), @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.FrexpResult(f16), scalar.frexph, x);
    if (T == f32) return mapUnaryStruct(f32, types.FrexpResult(f32), scalar.frexpf, x);
    if (T == f64) return mapUnaryStruct(f64, types.FrexpResult(f64), scalar.frexpd, x);
    @compileError("zmath.frexpv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the frexp of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn frexphv(x: anytype) ResultArray(types.FrexpResult(f16), @TypeOf(x)) {
    return mapUnaryStruct(f16, types.FrexpResult(f16), scalar.frexph, x);
}

/// Computes the frexp of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn frexpfv(x: anytype) ResultArray(types.FrexpResult(f32), @TypeOf(x)) {
    return mapUnaryStruct(f32, types.FrexpResult(f32), scalar.frexpf, x);
}

/// Computes the frexp of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn frexpdv(x: anytype) ResultArray(types.FrexpResult(f64), @TypeOf(x)) {
    return mapUnaryStruct(f64, types.FrexpResult(f64), scalar.frexpd, x);
}

/// Computes the modf of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds the fractional and integral parts of lane `i`. The result
/// is a lane-ordered `[N]` array, not a `@Vector`.
pub fn modfv(x: anytype) ResultArray(types.ModfResult(std.meta.Child(@TypeOf(x))), @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.ModfResult(f16), scalar.modfh, x);
    if (T == f32) return mapUnaryStruct(f32, types.ModfResult(f32), scalar.modff, x);
    if (T == f64) return mapUnaryStruct(f64, types.ModfResult(f64), scalar.modfd, x);
    @compileError("zmath.modfv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the modf of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn modfhv(x: anytype) ResultArray(types.ModfResult(f16), @TypeOf(x)) {
    return mapUnaryStruct(f16, types.ModfResult(f16), scalar.modfh, x);
}

/// Computes the modf of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn modffv(x: anytype) ResultArray(types.ModfResult(f32), @TypeOf(x)) {
    return mapUnaryStruct(f32, types.ModfResult(f32), scalar.modff, x);
}

/// Computes the modf of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn modfdv(x: anytype) ResultArray(types.ModfResult(f64), @TypeOf(x)) {
    return mapUnaryStruct(f64, types.ModfResult(f64), scalar.modfd, x);
}

// ---------------------------------------------------------------------------
// IEEE canonicalization and payloads, one result struct per lane
// ---------------------------------------------------------------------------

/// Computes the canonicalize of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds the canonicalized lane `i` and whether the canonical form
/// was produced. The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn canonicalizev(x: anytype) ResultArray(types.CanonicalizeResult(std.meta.Child(@TypeOf(x))), @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.CanonicalizeResult(f16), ieee.canonicalizeh, x);
    if (T == f32) return mapUnaryStruct(f32, types.CanonicalizeResult(f32), ieee.canonicalizef, x);
    if (T == f64) return mapUnaryStruct(f64, types.CanonicalizeResult(f64), ieee.canonicalized, x);
    @compileError("zmath.canonicalizev supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the canonicalize of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn canonicalizehv(x: anytype) ResultArray(types.CanonicalizeResult(f16), @TypeOf(x)) {
    return mapUnaryStruct(f16, types.CanonicalizeResult(f16), ieee.canonicalizeh, x);
}

/// Computes the canonicalize of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn canonicalizefv(x: anytype) ResultArray(types.CanonicalizeResult(f32), @TypeOf(x)) {
    return mapUnaryStruct(f32, types.CanonicalizeResult(f32), ieee.canonicalizef, x);
}

/// Computes the canonicalize of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn canonicalizedv(x: anytype) ResultArray(types.CanonicalizeResult(f64), @TypeOf(x)) {
    return mapUnaryStruct(f64, types.CanonicalizeResult(f64), ieee.canonicalized, x);
}

/// Computes the setPayload of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds the quiet NaN carrying lane `i` as its payload, and whether
/// that payload is an admissible integer. The result is a lane-ordered `[N]`
/// array, not a `@Vector`.
pub fn setPayloadv(x: anytype) ResultArray(types.PayloadResult(std.meta.Child(@TypeOf(x))), @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.PayloadResult(f16), ieee.setPayloadh, x);
    if (T == f32) return mapUnaryStruct(f32, types.PayloadResult(f32), ieee.setPayloadf, x);
    if (T == f64) return mapUnaryStruct(f64, types.PayloadResult(f64), ieee.setPayloadd, x);
    @compileError("zmath.setPayloadv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the setPayload of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn setPayloadhv(x: anytype) ResultArray(types.PayloadResult(f16), @TypeOf(x)) {
    return mapUnaryStruct(f16, types.PayloadResult(f16), ieee.setPayloadh, x);
}

/// Computes the setPayload of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn setPayloadfv(x: anytype) ResultArray(types.PayloadResult(f32), @TypeOf(x)) {
    return mapUnaryStruct(f32, types.PayloadResult(f32), ieee.setPayloadf, x);
}

/// Computes the setPayload of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn setPayloaddv(x: anytype) ResultArray(types.PayloadResult(f64), @TypeOf(x)) {
    return mapUnaryStruct(f64, types.PayloadResult(f64), ieee.setPayloadd, x);
}

/// Computes the setPayloadSig of every lane in `x` with the scalar kernel for its lane type.
///
/// Result `i` holds the signaling NaN carrying lane `i` as its payload, and
/// whether that payload is an admissible non-zero integer. The result is a
/// lane-ordered `[N]` array, not a `@Vector`.
pub fn setPayloadSigv(x: anytype) ResultArray(types.PayloadResult(std.meta.Child(@TypeOf(x))), @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapUnaryStruct(f16, types.PayloadResult(f16), ieee.setPayloadSigh, x);
    if (T == f32) return mapUnaryStruct(f32, types.PayloadResult(f32), ieee.setPayloadSigf, x);
    if (T == f64) return mapUnaryStruct(f64, types.PayloadResult(f64), ieee.setPayloadSigd, x);
    @compileError("zmath.setPayloadSigv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the setPayloadSig of every lane in `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn setPayloadSighv(x: anytype) ResultArray(types.PayloadResult(f16), @TypeOf(x)) {
    return mapUnaryStruct(f16, types.PayloadResult(f16), ieee.setPayloadSigh, x);
}

/// Computes the setPayloadSig of every lane in `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn setPayloadSigfv(x: anytype) ResultArray(types.PayloadResult(f32), @TypeOf(x)) {
    return mapUnaryStruct(f32, types.PayloadResult(f32), ieee.setPayloadSigf, x);
}

/// Computes the setPayloadSig of every lane in `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn setPayloadSigdv(x: anytype) ResultArray(types.PayloadResult(f64), @TypeOf(x)) {
    return mapUnaryStruct(f64, types.PayloadResult(f64), ieee.setPayloadSigd, x);
}

// ---------------------------------------------------------------------------
// Binary operations, one result struct per lane
// ---------------------------------------------------------------------------

/// Computes the atan2 of corresponding lanes of `y` and `x` with the scalar kernel for their lane type.
///
/// Result `i` holds the angle for lanes `i` and whether it is resolved. The
/// result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn atan2v(y: anytype, x: anytype) ResultArray(types.Atan2Result(std.meta.Child(@TypeOf(y))), @TypeOf(y)) {
    const T = std.meta.Child(@TypeOf(y));
    if (T == f16) return mapBinaryStruct(f16, types.Atan2Result(f16), scalar.atan2h, y, x);
    if (T == f32) return mapBinaryStruct(f32, types.Atan2Result(f32), scalar.atan2f, y, x);
    if (T == f64) return mapBinaryStruct(f64, types.Atan2Result(f64), scalar.atan2d, y, x);
    @compileError("zmath.atan2v supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the atan2 of corresponding lanes of `y` and `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn atan2hv(y: anytype, x: anytype) ResultArray(types.Atan2Result(f16), @TypeOf(y)) {
    return mapBinaryStruct(f16, types.Atan2Result(f16), scalar.atan2h, y, x);
}

/// Computes the atan2 of corresponding lanes of `y` and `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn atan2fv(y: anytype, x: anytype) ResultArray(types.Atan2Result(f32), @TypeOf(y)) {
    return mapBinaryStruct(f32, types.Atan2Result(f32), scalar.atan2f, y, x);
}

/// Computes the atan2 of corresponding lanes of `y` and `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn atan2dv(y: anytype, x: anytype) ResultArray(types.Atan2Result(f64), @TypeOf(y)) {
    return mapBinaryStruct(f64, types.Atan2Result(f64), scalar.atan2d, y, x);
}

/// Computes the atan2pi of corresponding lanes of `y` and `x` with the scalar kernel for their lane type.
///
/// Result `i` holds the angle in half-revolutions for lanes `i` and whether it
/// is resolved. The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn atan2piv(y: anytype, x: anytype) ResultArray(types.Atan2Result(std.meta.Child(@TypeOf(y))), @TypeOf(y)) {
    const T = std.meta.Child(@TypeOf(y));
    if (T == f16) return mapBinaryStruct(f16, types.Atan2Result(f16), scalar.atan2pih, y, x);
    if (T == f32) return mapBinaryStruct(f32, types.Atan2Result(f32), scalar.atan2pif, y, x);
    if (T == f64) return mapBinaryStruct(f64, types.Atan2Result(f64), scalar.atan2pid, y, x);
    @compileError("zmath.atan2piv supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the atan2pi of corresponding lanes of `y` and `x` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn atan2pihv(y: anytype, x: anytype) ResultArray(types.Atan2Result(f16), @TypeOf(y)) {
    return mapBinaryStruct(f16, types.Atan2Result(f16), scalar.atan2pih, y, x);
}

/// Computes the atan2pi of corresponding lanes of `y` and `x` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn atan2pifv(y: anytype, x: anytype) ResultArray(types.Atan2Result(f32), @TypeOf(y)) {
    return mapBinaryStruct(f32, types.Atan2Result(f32), scalar.atan2pif, y, x);
}

/// Computes the atan2pi of corresponding lanes of `y` and `x` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn atan2pidv(y: anytype, x: anytype) ResultArray(types.Atan2Result(f64), @TypeOf(y)) {
    return mapBinaryStruct(f64, types.Atan2Result(f64), scalar.atan2pid, y, x);
}

/// Computes the remquo of corresponding lanes of `x` and `y` with the scalar kernel for their lane type.
///
/// Result `i` holds the rounded remainder of lanes `i` and the low bits of the
/// quotient. The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn remquov(x: anytype, y: anytype) ResultArray(types.RemquoResult(std.meta.Child(@TypeOf(x))), @TypeOf(x)) {
    const T = std.meta.Child(@TypeOf(x));
    if (T == f16) return mapBinaryStruct(f16, types.RemquoResult(f16), scalar.remquoh, x, y);
    if (T == f32) return mapBinaryStruct(f32, types.RemquoResult(f32), scalar.remquof, x, y);
    if (T == f64) return mapBinaryStruct(f64, types.RemquoResult(f64), scalar.remquod, x, y);
    @compileError("zmath.remquov supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the remquo of corresponding lanes of `x` and `y` with the binary16 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn remquohv(x: anytype, y: anytype) ResultArray(types.RemquoResult(f16), @TypeOf(x)) {
    return mapBinaryStruct(f16, types.RemquoResult(f16), scalar.remquoh, x, y);
}

/// Computes the remquo of corresponding lanes of `x` and `y` with the binary32 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn remquofv(x: anytype, y: anytype) ResultArray(types.RemquoResult(f32), @TypeOf(x)) {
    return mapBinaryStruct(f32, types.RemquoResult(f32), scalar.remquof, x, y);
}

/// Computes the remquo of corresponding lanes of `x` and `y` with the binary64 scalar kernel.
///
/// The result is a lane-ordered `[N]` array, not a `@Vector`.
pub fn remquodv(x: anytype, y: anytype) ResultArray(types.RemquoResult(f64), @TypeOf(x)) {
    return mapBinaryStruct(f64, types.RemquoResult(f64), scalar.remquod, x, y);
}

comptime {
    _ = lrintv;
    _ = lrinthv;
    _ = lrintfv;
    _ = lrintdv;
    _ = llrintv;
    _ = llrinthv;
    _ = llrintfv;
    _ = llrintdv;
    _ = lroundv;
    _ = lroundhv;
    _ = lroundfv;
    _ = lrounddv;
    _ = llroundv;
    _ = llroundhv;
    _ = llroundfv;
    _ = llrounddv;
    _ = lgammav;
    _ = lgammahv;
    _ = lgammafv;
    _ = lgammadv;
    _ = frexpv;
    _ = frexphv;
    _ = frexpfv;
    _ = frexpdv;
    _ = modfv;
    _ = modfhv;
    _ = modffv;
    _ = modfdv;
    _ = canonicalizev;
    _ = canonicalizehv;
    _ = canonicalizefv;
    _ = canonicalizedv;
    _ = setPayloadv;
    _ = setPayloadhv;
    _ = setPayloadfv;
    _ = setPayloaddv;
    _ = setPayloadSigv;
    _ = setPayloadSighv;
    _ = setPayloadSigfv;
    _ = setPayloadSigdv;
    _ = atan2v;
    _ = atan2hv;
    _ = atan2fv;
    _ = atan2dv;
    _ = atan2piv;
    _ = atan2pihv;
    _ = atan2pifv;
    _ = atan2pidv;
    _ = remquov;
    _ = remquohv;
    _ = remquofv;
    _ = remquodv;
}
