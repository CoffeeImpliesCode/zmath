//! By convention, root.zig is the root source file when making a package.
const std = @import("std");
const Io = std.Io;

pub const ieee = @import("ieee.zig");
pub const scalar = @import("scalar.zig");
const generic_scalar_test = @import("generic_scalar_test.zig");

/// Computes e raised to `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn exp(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (T == f16) return scalar.exph(value);
    if (T == f32) return scalar.expf(value);
    if (T == f64) return scalar.expd(value);
    @compileError("zmath.exp supports only f16, f32, and f64");
}

/// Computes the natural logarithm of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn log(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (T == f16) return scalar.logh(value);
    if (T == f32) return scalar.logf(value);
    if (T == f64) return scalar.logd(value);
    @compileError("zmath.log supports only f16, f32, and f64");
}

/// Computes the hyperbolic sine of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn sinh(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
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
    if (T == f16) return scalar.sinh(value);
    if (T == f32) return scalar.sinf(value);
    if (T == f64) return scalar.sind(value);
    @compileError("zmath.sin supports only f16, f32, and f64");
}

/// Computes the cosine of `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn cos(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (T == f16) return scalar.cosh(value);
    if (T == f32) return scalar.cosf(value);
    if (T == f64) return scalar.cosd(value);
    @compileError("zmath.cos supports only f16, f32, and f64");
}
/// Computes the hyperbolic cosine of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn cosh(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (T == f16) return scalar.coshh(value);
    if (T == f32) return scalar.coshf(value);
    if (T == f64) return scalar.coshd(value);
    @compileError("zmath.cosh supports only f16, f32, and f64");
}

/// Computes the tangent of `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn tan(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (T == f16) return scalar.tanh(value);
    if (T == f32) return scalar.tanf(value);
    if (T == f64) return scalar.tand(value);
    @compileError("zmath.tan supports only f16, f32, and f64");
}

/// Computes the hyperbolic tangent of `value` with the scalar implementation
/// for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn tanh(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
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
    if (T == f16) return scalar.rsqrth(value);
    if (T == f32) return scalar.rsqrtf(value);
    if (T == f64) return scalar.rsqrtd(value);
    @compileError("zmath.rsqrt supports only f16, f32, and f64");
}

/// Computes the inverse sine of `value` with the scalar implementation for its type.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn asin(value: anytype) @TypeOf(value) {
    const T = @TypeOf(value);
    if (T == f16) return scalar.asinh(value);
    if (T == f32) return scalar.asinf(value);
    if (T == f64) return scalar.asind(value);
    @compileError("zmath.asin supports only f16, f32, and f64");
}

/// This is a documentation comment to explain the `printAnotherMessage` function below.
///
/// Accepting an `Io.Writer` instance is a handy way to write reusable code.
pub fn printAnotherMessage(writer: *Io.Writer) Io.Writer.Error!void {
    try writer.print("Run `zig build test` to run the tests.\n", .{});
}

pub fn add(a: i32, b: i32) i32 {
    return a + b;
}

test {
    _ = @import("ieee_test.zig");
    _ = @import("ieee_edge_test.zig");
    _ = @import("scalar_test.zig");
    _ = @import("scalar_cosh_test.zig");
    _ = @import("scalar_boundary_test.zig");
    _ = @import("scalar_rsqrt_test.zig");
    _ = @import("api_suffix_test.zig");
    _ = @import("ieee_h_complete_test.zig");
    _ = @import("ieee_exhaustive_h_test.zig");
    _ = @import("scalar_exp_test.zig");
    _ = @import("scalar_expm1_test.zig");
    _ = @import("scalar_log_test.zig");
    _ = @import("scalar_log1p_test.zig");
    _ = @import("scalar_sin_test.zig");
    _ = @import("scalar_cos_test.zig");
    _ = @import("scalar_sind_test.zig");
    _ = @import("scalar_tan_test.zig");
    _ = @import("scalar_stdlib_agreement_test.zig");
    _ = @import("scalar_exp2_test.zig");
    _ = @import("scalar_cbrt_test.zig");
    _ = @import("scalar_hypot_test.zig");
    _ = @import("scalar_asin_test.zig");
    _ = @import("scalar_atan_test.zig");
    _ = @import("scalar_sinh_test.zig");
    _ = @import("scalar_tanh_test.zig");
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

test "basic add functionality" {
    try std.testing.expect(add(3, 7) == 10);
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

test "generic scalar dispatch covers all formats" {
    try generic_scalar_test.run(@This());
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
