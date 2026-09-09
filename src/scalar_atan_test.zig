const std = @import("std");
const scalar = @import("scalar.zig");

fn expectAtanhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.atanh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAtanfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.atanf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAtandBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.atand(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "Atan family exposes lower-case suffixes" {
    comptime {
        std.debug.assert(@TypeOf(scalar.atanh(@as(f16, 0.0))) == f16);
        std.debug.assert(@TypeOf(scalar.atanf(@as(f32, 0.0))) == f32);
        std.debug.assert(@TypeOf(scalar.atand(@as(f64, 0.0))) == f64);
    }
}

test "Atan family preserves signed zero" {
    try expectAtanhBits(0x0000, 0x0000);
    try expectAtanhBits(0x8000, 0x8000);
    try expectAtanfBits(0x0000_0000, 0x0000_0000);
    try expectAtanfBits(0x8000_0000, 0x8000_0000);
    try expectAtandBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectAtandBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
}

test "Atan family rounds one and infinities" {
    try expectAtanhBits(0x3c00, 0x3a48);
    try expectAtanhBits(0xbc00, 0xba48);
    try expectAtanhBits(0x7c00, 0x3e48);
    try expectAtanhBits(0xfc00, 0xbe48);

    try expectAtanfBits(0x3f80_0000, 0x3f49_0fdb);
    try expectAtanfBits(0xbf80_0000, 0xbf49_0fdb);
    try expectAtanfBits(0x7f80_0000, 0x3fc9_0fdb);
    try expectAtanfBits(0xff80_0000, 0xbfc9_0fdb);

    try expectAtandBits(0x3ff0_0000_0000_0000, 0x3fe9_21fb_5444_2d18);
    try expectAtandBits(0xbff0_0000_0000_0000, 0xbfe9_21fb_5444_2d18);
    try expectAtandBits(0x7ff0_0000_0000_0000, 0x3ff9_21fb_5444_2d18);
    try expectAtandBits(0xfff0_0000_0000_0000, 0xbff9_21fb_5444_2d18);
}

test "binary16 Atan exceptional finite vectors" {
    try expectAtanhBits(0x32b5, 0x329d);
    try expectAtanhBits(0xb2b5, 0xb29d);
    try expectAtanhBits(0x48ae, 0x3ddb);
    try expectAtanhBits(0xc8ae, 0xbddb);
    try expectAtanhBits(0x5619, 0x3e3e);
    try expectAtanhBits(0xd619, 0xbe3e);
}

test "Atan family matches CORE-MATH rounding witnesses" {
    try expectAtanfBits(0x3880_0000, 0x3880_0000);
    try expectAtanfBits(0xb880_0000, 0xb880_0000);

    try expectAtandBits(0x434d_0296_7c31_cdb5, 0x3ff9_21fb_5444_2d18);
    try expectAtandBits(0xc34d_0296_7c31_cdb5, 0xbff9_21fb_5444_2d18);
    try expectAtandBits(0x3fe5_f457_e4f4_812e, 0x3fe3_3e04_6c84_3286);
    try expectAtandBits(0xbfe5_f457_e4f4_812e, 0xbfe3_3e04_6c84_3286);
}

test "binary64 atan refines small non-table arguments" {
    // MPFR-512 and CORE-MATH agree on this refinement-path result.
    try expectAtandBits(0x3f7d_7065_405b_4000, 0x3f7d_7044_0819_2854);
}

test "Atan family quiets NaNs and preserves payload and sign" {
    try expectAtanhBits(0x7c01, 0x7e01);
    try expectAtanhBits(0xfe55, 0xfe55);
    try expectAtanfBits(0x7f80_0001, 0x7fc0_0001);
    try expectAtanfBits(0xffc1_2345, 0xffc1_2345);
    try expectAtandBits(0x7ff0_0000_0000_0001, 0x7ff8_0000_0000_0001);
    try expectAtandBits(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}
