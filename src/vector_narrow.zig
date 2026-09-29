//! Vector (`@Vector`) overloads for the narrowing and ternary operations.
//!
//! Every operation here applies the *identical* accepted scalar kernel to each
//! lane, so lane `i` of the result is bit-for-bit equal to the scalar result for
//! lane `i` of the input. This preserves the schema-3 `max_ulp = 0` contract and
//! the immutable corpus without new evidence: the vector path is a per-lane
//! application of the same scalar code, never a re-derived polynomial.
//!
//! Two surfaces live here.
//!
//! - The six C23 narrowing operations `fadd`, `fsub`, `fmul`, `fdiv`, `ffma` and
//!   `fsqrt` take `f64` lanes and produce `f32` lanes. The lane type and the
//!   result element type are therefore fixed, which is why they expose no
//!   per-format entry points: there is no `T` to select. The surface here
//!   mirrors exactly the root surface in `root.zig`, and the adjacent
//!   `f16addf32` to `f32sqrtf64` kernels are deliberately not promoted.
//! - The ternary float-to-float operation `fma` keeps its lane type, so it does
//!   get per-format entry points: `fmav` selects the kernel from the lane type,
//!   and `fmahv`, `fmafv` and `fmadv` pin binary16, binary32 and binary64.
//!
//! The narrowing operations are the first users of a result element type that
//! differs from the lane type, which is what the explicit `R` parameter of the
//! shared core exists for. Every lane value is assigned whole by that core, so
//! nothing here is reinterpreted: this family produces no result struct, and
//! where a family does, its auto-layout fields are assigned, never bit-cast.

const std = @import("std");
const core = @import("vector_core.zig");
const scalar = @import("scalar.zig");
const narrow = @import("scalar_narrow.zig");

// ---------------------------------------------------------------------------
// Shared width and result types
// ---------------------------------------------------------------------------

/// Returns the lane count of the vector type `X`.
fn widthOf(comptime X: type) usize {
    return switch (@typeInfo(X)) {
        .vector => |v| v.len,
        else => @compileError("zmath vector operations require a @Vector argument"),
    };
}

/// The `@Vector(N, f32)` result type for the width of the vector type `X`.
///
/// Every narrowing operation in this module shares this result element type,
/// so none of them returns its argument type.
fn Narrowed(comptime X: anytype) type {
    return @Vector(widthOf(X), f32);
}

// ---------------------------------------------------------------------------
// Narrowing operations: f64 lanes in, f32 lanes out
// ---------------------------------------------------------------------------

/// Computes the narrowing add of corresponding lanes of `a` and `b` with the
/// binary64 scalar kernel.
pub fn faddv(a: anytype, b: anytype) Narrowed(@TypeOf(a)) {
    return core.mapBinaryR(f64, f32, narrow.fadd, a, b);
}

/// Computes the narrowing subtract of corresponding lanes of `a` and `b` with
/// the binary64 scalar kernel.
pub fn fsubv(a: anytype, b: anytype) Narrowed(@TypeOf(a)) {
    return core.mapBinaryR(f64, f32, narrow.fsub, a, b);
}

/// Computes the narrowing multiply of corresponding lanes of `a` and `b` with
/// the binary64 scalar kernel.
pub fn fmulv(a: anytype, b: anytype) Narrowed(@TypeOf(a)) {
    return core.mapBinaryR(f64, f32, narrow.fmul, a, b);
}

/// Computes the narrowing divide of corresponding lanes of `a` and `b` with the
/// binary64 scalar kernel.
pub fn fdivv(a: anytype, b: anytype) Narrowed(@TypeOf(a)) {
    return core.mapBinaryR(f64, f32, narrow.fdiv, a, b);
}

/// Computes the narrowing fused multiply-add of corresponding lanes of `a`, `b`
/// and `c` with the binary64 scalar kernel.
pub fn ffmav(a: anytype, b: anytype, c: anytype) Narrowed(@TypeOf(a)) {
    return core.mapTernaryR(f64, f32, narrow.ffma, a, b, c);
}

/// Computes the narrowing square root of every lane in `x` with the binary64
/// scalar kernel.
pub fn fsqrtv(x: anytype) Narrowed(@TypeOf(x)) {
    return core.mapUnaryR(f64, f32, narrow.fsqrt, x);
}

// ---------------------------------------------------------------------------
// Ternary float-to-float operation
// ---------------------------------------------------------------------------

/// Computes the fused multiply-add of corresponding lanes of `a`, `b` and `c`
/// with the scalar kernel for their lane type.
pub fn fmav(a: anytype, b: anytype, c: anytype) @TypeOf(a) {
    const T = std.meta.Child(@TypeOf(a));
    if (T == f16) return core.mapTernaryR(f16, f16, scalar.fmah, a, b, c);
    if (T == f32) return core.mapTernaryR(f32, f32, scalar.fmaf, a, b, c);
    if (T == f64) return core.mapTernaryR(f64, f64, scalar.fmad, a, b, c);
    @compileError("zmath.fmav supports only @Vector(N, f16), @Vector(N, f32), and @Vector(N, f64)");
}

/// Computes the fused multiply-add of corresponding lanes of `a`, `b` and `c`
/// with the binary16 scalar kernel.
pub fn fmahv(a: anytype, b: anytype, c: anytype) @TypeOf(a) {
    return core.mapTernaryR(f16, f16, scalar.fmah, a, b, c);
}

/// Computes the fused multiply-add of corresponding lanes of `a`, `b` and `c`
/// with the binary32 scalar kernel.
pub fn fmafv(a: anytype, b: anytype, c: anytype) @TypeOf(a) {
    return core.mapTernaryR(f32, f32, scalar.fmaf, a, b, c);
}

/// Computes the fused multiply-add of corresponding lanes of `a`, `b` and `c`
/// with the binary64 scalar kernel.
pub fn fmadv(a: anytype, b: anytype, c: anytype) @TypeOf(a) {
    return core.mapTernaryR(f64, f64, scalar.fmad, a, b, c);
}

comptime {
    _ = faddv;
    _ = fsubv;
    _ = fmulv;
    _ = fdivv;
    _ = ffmav;
    _ = fsqrtv;
    _ = fmav;
    _ = fmahv;
    _ = fmafv;
    _ = fmadv;
}
