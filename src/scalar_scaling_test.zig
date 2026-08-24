// These tests compare raw IEEE encodings across the full H/F/D scaling surface.
// They cover identities, signed zeros, gradual underflow, ties-to-even, overflow,
// infinities, and both quiet and signaling NaN payload preservation.

const std = @import("std");
const scaling = @import("scalar_scaling.zig");

fn expectHBits(expected_bits: u16, actual: f16) !void {
    const actual_bits: u16 = @bitCast(actual);
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectFBits(expected_bits: u32, actual: f32) !void {
    const actual_bits: u32 = @bitCast(actual);
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectDBits(expected_bits: u64, actual: f64) !void {
    const actual_bits: u64 = @bitCast(actual);
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAllH(input_bits: u16, exponent: i32, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    try expectHBits(expected_bits, scaling.ldexph(input, exponent));
    try expectHBits(expected_bits, scaling.scalbnh(input, exponent));
    try expectHBits(expected_bits, scaling.scalblnh(input, @as(i64, exponent)));
}

fn expectAllF(input_bits: u32, exponent: i32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    try expectFBits(expected_bits, scaling.ldexpf(input, exponent));
    try expectFBits(expected_bits, scaling.scalbnf(input, exponent));
    try expectFBits(expected_bits, scaling.scalblnf(input, @as(i64, exponent)));
}

fn expectAllD(input_bits: u64, exponent: i32, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    try expectDBits(expected_bits, scaling.ldexpd(input, exponent));
    try expectDBits(expected_bits, scaling.scalbnd(input, exponent));
    try expectDBits(expected_bits, scaling.scalblnd(input, @as(i64, exponent)));
}

test "scaling family signatures use C23 exponent widths" {
    comptime {
        std.debug.assert(@TypeOf(scaling.ldexph(@as(f16, 1.0), @as(i32, 0))) == f16);
        std.debug.assert(@TypeOf(scaling.ldexpf(@as(f32, 1.0), @as(i32, 0))) == f32);
        std.debug.assert(@TypeOf(scaling.ldexpd(@as(f64, 1.0), @as(i32, 0))) == f64);
        std.debug.assert(@TypeOf(scaling.scalbnh(@as(f16, 1.0), @as(i32, 0))) == f16);
        std.debug.assert(@TypeOf(scaling.scalbnf(@as(f32, 1.0), @as(i32, 0))) == f32);
        std.debug.assert(@TypeOf(scaling.scalbnd(@as(f64, 1.0), @as(i32, 0))) == f64);
        std.debug.assert(@TypeOf(scaling.scalblnh(@as(f16, 1.0), @as(i64, 0))) == f16);
        std.debug.assert(@TypeOf(scaling.scalblnf(@as(f32, 1.0), @as(i64, 0))) == f32);
        std.debug.assert(@TypeOf(scaling.scalblnd(@as(f64, 1.0), @as(i64, 0))) == f64);
    }
}

test "binary16 scaling preserves identities and special encodings" {
    try expectAllH(0x0001, 0, 0x0001);
    try expectAllH(0x3c00, 0, 0x3c00);
    try expectAllH(0x7bff, 0, 0x7bff);
    try expectAllH(0x0000, 0, 0x0000);
    try expectAllH(0x8000, 0, 0x8000);
    try expectAllH(0x0000, std.math.maxInt(i32), 0x0000);
    try expectAllH(0x8000, std.math.minInt(i32), 0x8000);
    try expectAllH(0x7c00, std.math.minInt(i32), 0x7c00);
    try expectAllH(0xfc00, std.math.maxInt(i32), 0xfc00);
    try expectAllH(0x7e55, 37, 0x7e55);
    try expectAllH(0xfd55, -37, 0xff55);
}

test "binary16 scaling crosses range boundaries with ties to even" {
    try expectAllH(0x0400, -1, 0x0200);
    try expectAllH(0x0001, 10, 0x0400);
    try expectAllH(0x03ff, 1, 0x07fe);
    try expectAllH(0x07ff, -1, 0x0400);

    try expectAllH(0x0400, -11, 0x0000);
    try expectAllH(0x0401, -11, 0x0001);
    try expectAllH(0x0600, -10, 0x0002);
    try expectAllH(0x8001, -1, 0x8000);

    try expectAllH(0x7bff, 1, 0x7c00);
    try expectAllH(0xfbff, 1, 0xfc00);
    try expectHBits(0x7c00, scaling.scalblnh(@bitCast(@as(u16, 0x3c00)), std.math.maxInt(i64)));
    try expectHBits(0x8000, scaling.scalblnh(@bitCast(@as(u16, 0xbc00)), std.math.minInt(i64)));
}

test "binary32 scaling preserves identities and special encodings" {
    try expectAllF(0x0000_0001, 0, 0x0000_0001);
    try expectAllF(0x3f80_0000, 0, 0x3f80_0000);
    try expectAllF(0x7f7f_ffff, 0, 0x7f7f_ffff);
    try expectAllF(0x0000_0000, 0, 0x0000_0000);
    try expectAllF(0x8000_0000, 0, 0x8000_0000);
    try expectAllF(0x0000_0000, std.math.maxInt(i32), 0x0000_0000);
    try expectAllF(0x8000_0000, std.math.minInt(i32), 0x8000_0000);
    try expectAllF(0x7f80_0000, std.math.minInt(i32), 0x7f80_0000);
    try expectAllF(0xff80_0000, std.math.maxInt(i32), 0xff80_0000);
    try expectAllF(0x7fc0_1234, 137, 0x7fc0_1234);
    try expectAllF(0xff80_1234, -137, 0xffc0_1234);
}

test "binary32 scaling crosses range boundaries with ties to even" {
    try expectAllF(0x0080_0000, -1, 0x0040_0000);
    try expectAllF(0x0000_0001, 23, 0x0080_0000);
    try expectAllF(0x007f_ffff, 1, 0x00ff_fffe);
    try expectAllF(0x00ff_ffff, -1, 0x0080_0000);

    try expectAllF(0x0080_0000, -24, 0x0000_0000);
    try expectAllF(0x0080_0001, -24, 0x0000_0001);
    try expectAllF(0x00c0_0000, -23, 0x0000_0002);
    try expectAllF(0x8000_0001, -1, 0x8000_0000);

    try expectAllF(0x7f7f_ffff, 1, 0x7f80_0000);
    try expectAllF(0xff7f_ffff, 1, 0xff80_0000);
    try expectFBits(
        0x7f80_0000,
        scaling.scalblnf(@bitCast(@as(u32, 0x3f80_0000)), std.math.maxInt(i64)),
    );
    try expectFBits(
        0x8000_0000,
        scaling.scalblnf(@bitCast(@as(u32, 0xbf80_0000)), std.math.minInt(i64)),
    );
}

test "binary64 scaling preserves identities and special encodings" {
    try expectAllD(0x0000_0000_0000_0001, 0, 0x0000_0000_0000_0001);
    try expectAllD(0x3ff0_0000_0000_0000, 0, 0x3ff0_0000_0000_0000);
    try expectAllD(0x7fef_ffff_ffff_ffff, 0, 0x7fef_ffff_ffff_ffff);
    try expectAllD(0x0000_0000_0000_0000, 0, 0x0000_0000_0000_0000);
    try expectAllD(0x8000_0000_0000_0000, 0, 0x8000_0000_0000_0000);
    try expectAllD(0x0000_0000_0000_0000, std.math.maxInt(i32), 0x0000_0000_0000_0000);
    try expectAllD(0x8000_0000_0000_0000, std.math.minInt(i32), 0x8000_0000_0000_0000);
    try expectAllD(0x7ff0_0000_0000_0000, std.math.minInt(i32), 0x7ff0_0000_0000_0000);
    try expectAllD(0xfff0_0000_0000_0000, std.math.maxInt(i32), 0xfff0_0000_0000_0000);
    try expectAllD(0x7ff8_0000_0000_1234, 1037, 0x7ff8_0000_0000_1234);
    try expectAllD(0xfff0_0000_0000_1234, -1037, 0xfff8_0000_0000_1234);
}

test "binary64 scaling crosses range boundaries with ties to even" {
    try expectAllD(0x0010_0000_0000_0000, -1, 0x0008_0000_0000_0000);
    try expectAllD(0x0000_0000_0000_0001, 52, 0x0010_0000_0000_0000);
    try expectAllD(0x000f_ffff_ffff_ffff, 1, 0x001f_ffff_ffff_fffe);
    try expectAllD(0x001f_ffff_ffff_ffff, -1, 0x0010_0000_0000_0000);

    try expectAllD(0x0010_0000_0000_0000, -53, 0x0000_0000_0000_0000);
    try expectAllD(0x0010_0000_0000_0001, -53, 0x0000_0000_0000_0001);
    try expectAllD(0x0018_0000_0000_0000, -52, 0x0000_0000_0000_0002);
    try expectAllD(0x8000_0000_0000_0001, -1, 0x8000_0000_0000_0000);

    try expectAllD(0x7fef_ffff_ffff_ffff, 1, 0x7ff0_0000_0000_0000);
    try expectAllD(0xffef_ffff_ffff_ffff, 1, 0xfff0_0000_0000_0000);
    try expectDBits(
        0x7ff0_0000_0000_0000,
        scaling.scalblnd(
            @bitCast(@as(u64, 0x3ff0_0000_0000_0000)),
            std.math.maxInt(i64),
        ),
    );
    try expectDBits(
        0x8000_0000_0000_0000,
        scaling.scalblnd(
            @bitCast(@as(u64, 0xbff0_0000_0000_0000)),
            std.math.minInt(i64),
        ),
    );
}
