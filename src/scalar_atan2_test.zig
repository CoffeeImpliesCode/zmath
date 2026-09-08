const std = @import("std");
const scalar = @import("scalar_atan2.zig");
const c23 = @import("scalar_c23_types.zig");

extern fn zmath_core_atan2d_test_subtract_two_limbs(
    h: *u64,
    m: *u64,
    l: *u64,
) callconv(.c) void;

fn expectAtan2dSubtraction(
    h_bits: u64,
    m_bits: u64,
    l_bits: u64,
    expected_h: u64,
    expected_m: u64,
    expected_l: u64,
) !void {
    var h = h_bits;
    var m = m_bits;
    var l = l_bits;
    zmath_core_atan2d_test_subtract_two_limbs(&h, &m, &l);
    try std.testing.expectEqual(expected_h, h);
    try std.testing.expectEqual(expected_m, m);
    try std.testing.expectEqual(expected_l, l);
}

fn expectAtan2hBits(y_bits: u16, x_bits: u16, expected_bits: u16) !void {
    const y: f16 = @bitCast(y_bits);
    const x: f16 = @bitCast(x_bits);
    const actual = scalar.atan2h(y, x);
    try std.testing.expect(actual.resolved);
    const actual_bits: u16 = @bitCast(actual.value);
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAtan2fBits(y_bits: u32, x_bits: u32, expected_bits: u32) !void {
    const y: f32 = @bitCast(y_bits);
    const x: f32 = @bitCast(x_bits);
    const actual = scalar.atan2f(y, x);
    try std.testing.expect(actual.resolved);
    const actual_bits: u32 = @bitCast(actual.value);
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAtan2dBits(y_bits: u64, x_bits: u64, expected_bits: u64) !void {
    const y: f64 = @bitCast(y_bits);
    const x: f64 = @bitCast(x_bits);
    const actual = scalar.atan2d(y, x);
    try std.testing.expect(actual.resolved);
    const actual_bits: u64 = @bitCast(actual.value);
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "Atan2 family exposes lower-case suffixes" {
    comptime {
        std.debug.assert(@TypeOf(scalar.atan2h(@as(f16, 0.0), @as(f16, 1.0))) == c23.Atan2Result(f16));
        std.debug.assert(@TypeOf(scalar.atan2f(@as(f32, 0.0), @as(f32, 1.0))) == c23.Atan2Result(f32));
        std.debug.assert(@TypeOf(scalar.atan2d(@as(f64, 0.0), @as(f64, 1.0))) == c23.Atan2Result(f64));
    }
}

test "Atan2 binary16 signed-zero quadrants" {
    try expectAtan2hBits(0x0000, 0x0000, 0x0000);
    try expectAtan2hBits(0x8000, 0x0000, 0x8000);
    try expectAtan2hBits(0x0000, 0x8000, 0x4248);
    try expectAtan2hBits(0x8000, 0x8000, 0xc248);
}

test "Atan2 binary32 signed-zero quadrants" {
    try expectAtan2fBits(0x0000_0000, 0x0000_0000, 0x0000_0000);
    try expectAtan2fBits(0x8000_0000, 0x0000_0000, 0x8000_0000);
    try expectAtan2fBits(0x0000_0000, 0x8000_0000, 0x4049_0fdb);
    try expectAtan2fBits(0x8000_0000, 0x8000_0000, 0xc049_0fdb);
}

test "Atan2 binary64 signed-zero quadrants" {
    try expectAtan2dBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectAtan2dBits(0x8000_0000_0000_0000, 0x0000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectAtan2dBits(0x0000_0000_0000_0000, 0x8000_0000_0000_0000, 0x4009_21fb_5444_2d18);
    try expectAtan2dBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0xc009_21fb_5444_2d18);
}

test "Atan2 binary16 finite axes" {
    try expectAtan2hBits(0x0000, 0x3c00, 0x0000);
    try expectAtan2hBits(0x8000, 0x3c00, 0x8000);
    try expectAtan2hBits(0x0000, 0xbc00, 0x4248);
    try expectAtan2hBits(0x8000, 0xbc00, 0xc248);
    try expectAtan2hBits(0x3c00, 0x0000, 0x3e48);
    try expectAtan2hBits(0xbc00, 0x0000, 0xbe48);
}

test "Atan2 binary32 finite axes" {
    try expectAtan2fBits(0x0000_0000, 0x3f80_0000, 0x0000_0000);
    try expectAtan2fBits(0x8000_0000, 0x3f80_0000, 0x8000_0000);
    try expectAtan2fBits(0x0000_0000, 0xbf80_0000, 0x4049_0fdb);
    try expectAtan2fBits(0x8000_0000, 0xbf80_0000, 0xc049_0fdb);
    try expectAtan2fBits(0x3f80_0000, 0x0000_0000, 0x3fc9_0fdb);
    try expectAtan2fBits(0xbf80_0000, 0x0000_0000, 0xbfc9_0fdb);
}

test "Atan2 binary64 finite axes" {
    try expectAtan2dBits(0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectAtan2dBits(0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectAtan2dBits(0x0000_0000_0000_0000, 0xbff0_0000_0000_0000, 0x4009_21fb_5444_2d18);
    try expectAtan2dBits(0x8000_0000_0000_0000, 0xbff0_0000_0000_0000, 0xc009_21fb_5444_2d18);
    try expectAtan2dBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000, 0x3ff9_21fb_5444_2d18);
    try expectAtan2dBits(0xbff0_0000_0000_0000, 0x0000_0000_0000_0000, 0xbff9_21fb_5444_2d18);
}

test "Atan2 binary16 double-infinity quadrants" {
    try expectAtan2hBits(0x7c00, 0x7c00, 0x3a48);
    try expectAtan2hBits(0x7c00, 0xfc00, 0x40b6);
    try expectAtan2hBits(0xfc00, 0x7c00, 0xba48);
    try expectAtan2hBits(0xfc00, 0xfc00, 0xc0b6);
}

test "Atan2 binary32 double-infinity quadrants" {
    try expectAtan2fBits(0x7f80_0000, 0x7f80_0000, 0x3f49_0fdb);
    try expectAtan2fBits(0x7f80_0000, 0xff80_0000, 0x4016_cbe4);
    try expectAtan2fBits(0xff80_0000, 0x7f80_0000, 0xbf49_0fdb);
    try expectAtan2fBits(0xff80_0000, 0xff80_0000, 0xc016_cbe4);
}

test "Atan2 binary64 double-infinity quadrants" {
    try expectAtan2dBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0x3fe9_21fb_5444_2d18);
    try expectAtan2dBits(0x7ff0_0000_0000_0000, 0xfff0_0000_0000_0000, 0x4002_d97c_7f33_21d2);
    try expectAtan2dBits(0xfff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0xbfe9_21fb_5444_2d18);
    try expectAtan2dBits(0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000, 0xc002_d97c_7f33_21d2);
}

test "Atan2 binary16 NaN precedence preserves payload and sign" {
    try expectAtan2hBits(0x7c01, 0x3c00, 0x7e01);
    try expectAtan2hBits(0x3c00, 0xfc01, 0xfe01);
    try expectAtan2hBits(0x7c2a, 0xfc55, 0x7e2a);
}

test "Atan2 binary32 NaN precedence preserves payload and sign" {
    try expectAtan2fBits(0x7f80_0001, 0x3f80_0000, 0x7fc0_0001);
    try expectAtan2fBits(0x3f80_0000, 0xff80_1234, 0xffc0_1234);
    try expectAtan2fBits(0x7f80_00a5, 0xff80_0055, 0x7fc0_00a5);
}

test "Atan2 binary64 NaN precedence preserves payload and sign" {
    try expectAtan2dBits(0x7ff0_0000_0000_0001, 0x3ff0_0000_0000_0000, 0x7ff8_0000_0000_0001);
    try expectAtan2dBits(0x3ff0_0000_0000_0000, 0xfff0_0000_0000_1234, 0xfff8_0000_0000_1234);
    try expectAtan2dBits(0x7ff0_0000_0000_00a5, 0xfff0_0000_0000_0055, 0x7ff8_0000_0000_00a5);
}

test "Atan2 binary16 single infinities and finite quadrants" {
    try expectAtan2hBits(0x7c00, 0x3c00, 0x3e48);
    try expectAtan2hBits(0xfc00, 0xbc00, 0xbe48);
    try expectAtan2hBits(0x3c00, 0x7c00, 0x0000);
    try expectAtan2hBits(0xbc00, 0x7c00, 0x8000);
    try expectAtan2hBits(0x3c00, 0xfc00, 0x4248);
    try expectAtan2hBits(0xbc00, 0xfc00, 0xc248);
    try expectAtan2hBits(0x3c00, 0x3c00, 0x3a48);
    try expectAtan2hBits(0x3c00, 0xbc00, 0x40b6);
    try expectAtan2hBits(0xbc00, 0xbc00, 0xc0b6);
    try expectAtan2hBits(0xbc00, 0x3c00, 0xba48);
}

test "Atan2 binary32 single infinities and finite quadrants" {
    try expectAtan2fBits(0x7f80_0000, 0x3f80_0000, 0x3fc9_0fdb);
    try expectAtan2fBits(0xff80_0000, 0xbf80_0000, 0xbfc9_0fdb);
    try expectAtan2fBits(0x3f80_0000, 0x7f80_0000, 0x0000_0000);
    try expectAtan2fBits(0xbf80_0000, 0x7f80_0000, 0x8000_0000);
    try expectAtan2fBits(0x3f80_0000, 0xff80_0000, 0x4049_0fdb);
    try expectAtan2fBits(0xbf80_0000, 0xff80_0000, 0xc049_0fdb);
    try expectAtan2fBits(0x3f80_0000, 0x3f80_0000, 0x3f49_0fdb);
    try expectAtan2fBits(0x3f80_0000, 0xbf80_0000, 0x4016_cbe4);
    try expectAtan2fBits(0xbf80_0000, 0xbf80_0000, 0xc016_cbe4);
    try expectAtan2fBits(0xbf80_0000, 0x3f80_0000, 0xbf49_0fdb);
}

test "Atan2 binary64 single infinities and finite quadrants" {
    try expectAtan2dBits(0x7ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x3ff9_21fb_5444_2d18);
    try expectAtan2dBits(0xfff0_0000_0000_0000, 0xbff0_0000_0000_0000, 0xbff9_21fb_5444_2d18);
    try expectAtan2dBits(0x3ff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectAtan2dBits(0xbff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectAtan2dBits(0x3ff0_0000_0000_0000, 0xfff0_0000_0000_0000, 0x4009_21fb_5444_2d18);
    try expectAtan2dBits(0xbff0_0000_0000_0000, 0xfff0_0000_0000_0000, 0xc009_21fb_5444_2d18);
    try expectAtan2dBits(0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x3fe9_21fb_5444_2d18);
    try expectAtan2dBits(0x3ff0_0000_0000_0000, 0xbff0_0000_0000_0000, 0x4002_d97c_7f33_21d2);
    try expectAtan2dBits(0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000, 0xc002_d97c_7f33_21d2);
    try expectAtan2dBits(0xbff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0xbfe9_21fb_5444_2d18);
}

test "Atan2 binary16 axis-adjacent and magnitude boundaries" {
    try expectAtan2hBits(0x0001, 0x3c00, 0x0001);
    try expectAtan2hBits(0x8001, 0x3c00, 0x8001);
    try expectAtan2hBits(0x3c00, 0x0001, 0x3e48);
    try expectAtan2hBits(0x0001, 0x7bff, 0x0000);
    try expectAtan2hBits(0x7bff, 0x7bff, 0x3a48);
}

test "Atan2 binary32 axis-adjacent and magnitude boundaries" {
    try expectAtan2fBits(0x0000_0001, 0x3f80_0000, 0x0000_0001);
    try expectAtan2fBits(0x8000_0001, 0x3f80_0000, 0x8000_0001);
    try expectAtan2fBits(0x3f80_0000, 0x0000_0001, 0x3fc9_0fdb);
    try expectAtan2fBits(0x0000_0001, 0x7f7f_ffff, 0x0000_0000);
    try expectAtan2fBits(0x7f7f_ffff, 0x7f7f_ffff, 0x3f49_0fdb);
}

test "Atan2 binary64 axis-adjacent and magnitude boundaries" {
    try expectAtan2dBits(0x0000_0000_0000_0001, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0001);
    try expectAtan2dBits(0x8000_0000_0000_0001, 0x3ff0_0000_0000_0000, 0x8000_0000_0000_0001);
    try expectAtan2dBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0001, 0x3ff9_21fb_5444_2d18);
    try expectAtan2dBits(0x0000_0000_0000_0001, 0x7fef_ffff_ffff_ffff, 0x0000_0000_0000_0000);
    try expectAtan2dBits(0x7fef_ffff_ffff_ffff, 0x7fef_ffff_ffff_ffff, 0x3fe9_21fb_5444_2d18);
}

test "Atan2 binary64 forced ambiguity reports unresolved" {
    const positive = scalar.testing.forcedAmbiguity(1.0, 1.0);
    const negative = scalar.testing.forcedAmbiguity(-1.0, -1.0);
    const positive_bits: u64 = @bitCast(positive.value);
    const negative_bits: u64 = @bitCast(negative.value);

    try std.testing.expect(!positive.resolved);
    try std.testing.expect(!negative.resolved);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), positive_bits);
    try std.testing.expectEqual(@as(u64, 0x7ff8_0000_0000_0000), negative_bits);
}

test "Atan2d subtract-two helper propagates low and middle borrows" {
    const max_limb: u64 = 0xffff_ffff_ffff_ffff;

    try expectAtan2dSubtraction(
        0x123,
        0,
        0,
        0x122,
        max_limb,
        0xffff_ffff_ffff_fffe,
    );
    try expectAtan2dSubtraction(
        0x123,
        0,
        1,
        0x122,
        max_limb,
        max_limb,
    );
    try expectAtan2dSubtraction(
        0x123,
        7,
        0,
        0x123,
        6,
        0xffff_ffff_ffff_fffe,
    );
    try expectAtan2dSubtraction(
        0x123,
        7,
        1,
        0x123,
        6,
        max_limb,
    );
    try expectAtan2dSubtraction(0x123, 0, 2, 0x123, 0, 0);
    try expectAtan2dSubtraction(0x123, 7, 2, 0x123, 7, 0);
}
