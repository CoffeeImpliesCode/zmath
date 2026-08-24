const std = @import("std");
const c23 = @import("scalar_c23_types.zig");
const decompose = @import("scalar_decompose.zig");

fn expectFrexpH(input_bits: u16, fraction_bits: u16, exponent: i32) !void {
    const result = decompose.frexph(@bitCast(input_bits));
    try std.testing.expectEqual(fraction_bits, @as(u16, @bitCast(result.fraction)));
    try std.testing.expectEqual(exponent, result.exponent);
}

fn expectFrexpF(input_bits: u32, fraction_bits: u32, exponent: i32) !void {
    const result = decompose.frexpf(@bitCast(input_bits));
    try std.testing.expectEqual(fraction_bits, @as(u32, @bitCast(result.fraction)));
    try std.testing.expectEqual(exponent, result.exponent);
}

fn expectFrexpD(input_bits: u64, fraction_bits: u64, exponent: i32) !void {
    const result = decompose.frexpd(@bitCast(input_bits));
    try std.testing.expectEqual(fraction_bits, @as(u64, @bitCast(result.fraction)));
    try std.testing.expectEqual(exponent, result.exponent);
}

fn expectModfH(input_bits: u16, fractional_bits: u16, integral_bits: u16) !void {
    const result = decompose.modfh(@bitCast(input_bits));
    try std.testing.expectEqual(fractional_bits, @as(u16, @bitCast(result.fractional)));
    try std.testing.expectEqual(integral_bits, @as(u16, @bitCast(result.integral)));
}

fn expectModfF(input_bits: u32, fractional_bits: u32, integral_bits: u32) !void {
    const result = decompose.modff(@bitCast(input_bits));
    try std.testing.expectEqual(fractional_bits, @as(u32, @bitCast(result.fractional)));
    try std.testing.expectEqual(integral_bits, @as(u32, @bitCast(result.integral)));
}

fn expectModfD(input_bits: u64, fractional_bits: u64, integral_bits: u64) !void {
    const result = decompose.modfd(@bitCast(input_bits));
    try std.testing.expectEqual(fractional_bits, @as(u64, @bitCast(result.fractional)));
    try std.testing.expectEqual(integral_bits, @as(u64, @bitCast(result.integral)));
}

fn expectLogbH(input_bits: u16, expected_bits: u16) !void {
    try std.testing.expectEqual(expected_bits, @as(u16, @bitCast(decompose.logbh(@bitCast(input_bits)))));
}

fn expectLogbF(input_bits: u32, expected_bits: u32) !void {
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(decompose.logbf(@bitCast(input_bits)))));
}

fn expectLogbD(input_bits: u64, expected_bits: u64) !void {
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(decompose.logbd(@bitCast(input_bits)))));
}

test "decomposition functions expose shared result containers" {
    comptime {
        std.debug.assert(@TypeOf(decompose.frexph(0)) == c23.FrexpResult(f16));
        std.debug.assert(@TypeOf(decompose.frexpf(0)) == c23.FrexpResult(f32));
        std.debug.assert(@TypeOf(decompose.frexpd(0)) == c23.FrexpResult(f64));
        std.debug.assert(@TypeOf(decompose.modfh(0)) == c23.ModfResult(f16));
        std.debug.assert(@TypeOf(decompose.modff(0)) == c23.ModfResult(f32));
        std.debug.assert(@TypeOf(decompose.modfd(0)) == c23.ModfResult(f64));
    }
}

test "frexp binary16 covers zeros, subnormal boundaries, normals, and infinities" {
    try expectFrexpH(0x0000, 0x0000, 0);
    try expectFrexpH(0x8000, 0x8000, 0);
    try expectFrexpH(0x0001, 0x3800, -23);
    try expectFrexpH(0x8001, 0xb800, -23);
    try expectFrexpH(0x03ff, 0x3bfe, -14);
    try expectFrexpH(0x83ff, 0xbbfe, -14);
    try expectFrexpH(0x0400, 0x3800, -13);
    try expectFrexpH(0x8400, 0xb800, -13);
    try expectFrexpH(0x3e00, 0x3a00, 1);
    try expectFrexpH(0xbe00, 0xba00, 1);
    try expectFrexpH(0x7c00, 0x7c00, 0);
    try expectFrexpH(0xfc00, 0xfc00, 0);
}

test "frexp binary32 covers zeros, subnormal boundaries, normals, and infinities" {
    try expectFrexpF(0x0000_0000, 0x0000_0000, 0);
    try expectFrexpF(0x8000_0000, 0x8000_0000, 0);
    try expectFrexpF(0x0000_0001, 0x3f00_0000, -148);
    try expectFrexpF(0x8000_0001, 0xbf00_0000, -148);
    try expectFrexpF(0x007f_ffff, 0x3f7f_fffe, -126);
    try expectFrexpF(0x807f_ffff, 0xbf7f_fffe, -126);
    try expectFrexpF(0x0080_0000, 0x3f00_0000, -125);
    try expectFrexpF(0x8080_0000, 0xbf00_0000, -125);
    try expectFrexpF(0x3fc0_0000, 0x3f40_0000, 1);
    try expectFrexpF(0xbfc0_0000, 0xbf40_0000, 1);
    try expectFrexpF(0x7f80_0000, 0x7f80_0000, 0);
    try expectFrexpF(0xff80_0000, 0xff80_0000, 0);
}

test "frexp binary64 covers zeros, subnormal boundaries, normals, and infinities" {
    try expectFrexpD(0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0);
    try expectFrexpD(0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0);
    try expectFrexpD(0x0000_0000_0000_0001, 0x3fe0_0000_0000_0000, -1073);
    try expectFrexpD(0x8000_0000_0000_0001, 0xbfe0_0000_0000_0000, -1073);
    try expectFrexpD(0x000f_ffff_ffff_ffff, 0x3fef_ffff_ffff_fffe, -1022);
    try expectFrexpD(0x800f_ffff_ffff_ffff, 0xbfef_ffff_ffff_fffe, -1022);
    try expectFrexpD(0x0010_0000_0000_0000, 0x3fe0_0000_0000_0000, -1021);
    try expectFrexpD(0x8010_0000_0000_0000, 0xbfe0_0000_0000_0000, -1021);
    try expectFrexpD(0x3ff8_0000_0000_0000, 0x3fe8_0000_0000_0000, 1);
    try expectFrexpD(0xbff8_0000_0000_0000, 0xbfe8_0000_0000_0000, 1);
    try expectFrexpD(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0);
    try expectFrexpD(0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000, 0);
}

test "frexp quiets signaling NaNs and preserves sign and payload" {
    try expectFrexpH(0x7e55, 0x7e55, 0);
    try expectFrexpH(0xfe55, 0xfe55, 0);
    try expectFrexpH(0x7d55, 0x7f55, 0);
    try expectFrexpH(0xfd55, 0xff55, 0);

    try expectFrexpF(0x7fc1_2345, 0x7fc1_2345, 0);
    try expectFrexpF(0xffc1_2345, 0xffc1_2345, 0);
    try expectFrexpF(0x7f81_2345, 0x7fc1_2345, 0);
    try expectFrexpF(0xff81_2345, 0xffc1_2345, 0);

    try expectFrexpD(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc, 0);
    try expectFrexpD(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc, 0);
    try expectFrexpD(0x7ff1_2345_6789_abcd, 0x7ff9_2345_6789_abcd, 0);
    try expectFrexpD(0xfff1_2345_6789_abcd, 0xfff9_2345_6789_abcd, 0);
}

test "modf binary16 preserves signs and splits fractional values" {
    try expectModfH(0x0000, 0x0000, 0x0000);
    try expectModfH(0x8000, 0x8000, 0x8000);
    try expectModfH(0x0001, 0x0001, 0x0000);
    try expectModfH(0x8001, 0x8001, 0x8000);
    try expectModfH(0x03ff, 0x03ff, 0x0000);
    try expectModfH(0x83ff, 0x83ff, 0x8000);
    try expectModfH(0x0400, 0x0400, 0x0000);
    try expectModfH(0x8400, 0x8400, 0x8000);
    try expectModfH(0x3800, 0x3800, 0x0000);
    try expectModfH(0xb800, 0xb800, 0x8000);
    try expectModfH(0x3e00, 0x3800, 0x3c00);
    try expectModfH(0xbe00, 0xb800, 0xbc00);
    try expectModfH(0x4300, 0x3800, 0x4200);
    try expectModfH(0xc300, 0xb800, 0xc200);
    try expectModfH(0x4200, 0x0000, 0x4200);
    try expectModfH(0xc200, 0x8000, 0xc200);
    try expectModfH(0x6400, 0x0000, 0x6400);
    try expectModfH(0x7c00, 0x0000, 0x7c00);
    try expectModfH(0xfc00, 0x8000, 0xfc00);
}

test "modf binary32 preserves signs and splits fractional values" {
    try expectModfF(0x0000_0000, 0x0000_0000, 0x0000_0000);
    try expectModfF(0x8000_0000, 0x8000_0000, 0x8000_0000);
    try expectModfF(0x0000_0001, 0x0000_0001, 0x0000_0000);
    try expectModfF(0x8000_0001, 0x8000_0001, 0x8000_0000);
    try expectModfF(0x007f_ffff, 0x007f_ffff, 0x0000_0000);
    try expectModfF(0x807f_ffff, 0x807f_ffff, 0x8000_0000);
    try expectModfF(0x0080_0000, 0x0080_0000, 0x0000_0000);
    try expectModfF(0x8080_0000, 0x8080_0000, 0x8000_0000);
    try expectModfF(0x3f00_0000, 0x3f00_0000, 0x0000_0000);
    try expectModfF(0xbf00_0000, 0xbf00_0000, 0x8000_0000);
    try expectModfF(0x3fc0_0000, 0x3f00_0000, 0x3f80_0000);
    try expectModfF(0xbfc0_0000, 0xbf00_0000, 0xbf80_0000);
    try expectModfF(0x4060_0000, 0x3f00_0000, 0x4040_0000);
    try expectModfF(0xc060_0000, 0xbf00_0000, 0xc040_0000);
    try expectModfF(0x4040_0000, 0x0000_0000, 0x4040_0000);
    try expectModfF(0xc040_0000, 0x8000_0000, 0xc040_0000);
    try expectModfF(0x4b00_0000, 0x0000_0000, 0x4b00_0000);
    try expectModfF(0x7f80_0000, 0x0000_0000, 0x7f80_0000);
    try expectModfF(0xff80_0000, 0x8000_0000, 0xff80_0000);
}

test "modf binary64 preserves signs and splits fractional values" {
    try expectModfD(0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectModfD(0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectModfD(0x0000_0000_0000_0001, 0x0000_0000_0000_0001, 0x0000_0000_0000_0000);
    try expectModfD(0x8000_0000_0000_0001, 0x8000_0000_0000_0001, 0x8000_0000_0000_0000);
    try expectModfD(0x3fe0_0000_0000_0000, 0x3fe0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectModfD(0xbfe0_0000_0000_0000, 0xbfe0_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectModfD(0x3ff8_0000_0000_0000, 0x3fe0_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectModfD(0xbff8_0000_0000_0000, 0xbfe0_0000_0000_0000, 0xbff0_0000_0000_0000);
    try expectModfD(0x000f_ffff_ffff_ffff, 0x000f_ffff_ffff_ffff, 0x0000_0000_0000_0000);
    try expectModfD(0x800f_ffff_ffff_ffff, 0x800f_ffff_ffff_ffff, 0x8000_0000_0000_0000);
    try expectModfD(0x0010_0000_0000_0000, 0x0010_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectModfD(0x8010_0000_0000_0000, 0x8010_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectModfD(0x400c_0000_0000_0000, 0x3fe0_0000_0000_0000, 0x4008_0000_0000_0000);
    try expectModfD(0xc00c_0000_0000_0000, 0xbfe0_0000_0000_0000, 0xc008_0000_0000_0000);
    try expectModfD(0x4008_0000_0000_0000, 0x0000_0000_0000_0000, 0x4008_0000_0000_0000);
    try expectModfD(0xc008_0000_0000_0000, 0x8000_0000_0000_0000, 0xc008_0000_0000_0000);
    try expectModfD(0x4330_0000_0000_0000, 0x0000_0000_0000_0000, 0x4330_0000_0000_0000);
    try expectModfD(0x7ff0_0000_0000_0000, 0x0000_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectModfD(0xfff0_0000_0000_0000, 0x8000_0000_0000_0000, 0xfff0_0000_0000_0000);
}

test "modf quiets signaling NaNs and preserves both result payloads" {
    try expectModfH(0x7e55, 0x7e55, 0x7e55);
    try expectModfH(0xfe55, 0xfe55, 0xfe55);
    try expectModfH(0x7d55, 0x7f55, 0x7f55);
    try expectModfH(0xfd55, 0xff55, 0xff55);

    try expectModfF(0x7fc1_2345, 0x7fc1_2345, 0x7fc1_2345);
    try expectModfF(0xffc1_2345, 0xffc1_2345, 0xffc1_2345);
    try expectModfF(0x7f81_2345, 0x7fc1_2345, 0x7fc1_2345);
    try expectModfF(0xff81_2345, 0xffc1_2345, 0xffc1_2345);

    try expectModfD(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectModfD(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectModfD(0x7ff1_2345_6789_abcd, 0x7ff9_2345_6789_abcd, 0x7ff9_2345_6789_abcd);
    try expectModfD(0xfff1_2345_6789_abcd, 0xfff9_2345_6789_abcd, 0xfff9_2345_6789_abcd);
}

test "ilogb and llogb return exact finite exponents and sentinels" {
    const quiet_nan_h: f16 = @bitCast(@as(u16, 0x7e55));
    const quiet_nan_f: f32 = @bitCast(@as(u32, 0x7fc1_2345));
    const quiet_nan_d: f64 = @bitCast(@as(u64, 0x7ff8_1234_5678_9abc));
    const positive_zero_h: f16 = @bitCast(@as(u16, 0x0000));
    const negative_zero_h: f16 = @bitCast(@as(u16, 0x8000));
    const positive_infinity_h: f16 = @bitCast(@as(u16, 0x7c00));
    const positive_zero_f: f32 = @bitCast(@as(u32, 0x0000_0000));
    const negative_zero_f: f32 = @bitCast(@as(u32, 0x8000_0000));
    const positive_infinity_f: f32 = @bitCast(@as(u32, 0x7f80_0000));
    const positive_zero_d: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0000));
    const negative_zero_d: f64 = @bitCast(@as(u64, 0x8000_0000_0000_0000));
    const positive_infinity_d: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0000));
    const min_subnormal_h: f16 = @bitCast(@as(u16, 0x0001));
    const max_subnormal_h: f16 = @bitCast(@as(u16, 0x03ff));
    const min_normal_h: f16 = @bitCast(@as(u16, 0x0400));
    const one_h: f16 = @bitCast(@as(u16, 0x3c00));
    const min_subnormal_f: f32 = @bitCast(@as(u32, 0x0000_0001));
    const max_subnormal_f: f32 = @bitCast(@as(u32, 0x007f_ffff));
    const min_normal_f: f32 = @bitCast(@as(u32, 0x0080_0000));
    const one_f: f32 = @bitCast(@as(u32, 0x3f80_0000));
    const min_subnormal_d: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0001));
    const max_subnormal_d: f64 = @bitCast(@as(u64, 0x000f_ffff_ffff_ffff));
    const min_normal_d: f64 = @bitCast(@as(u64, 0x0010_0000_0000_0000));
    const one_d: f64 = @bitCast(@as(u64, 0x3ff0_0000_0000_0000));

    try std.testing.expectEqual(std.math.minInt(i32), decompose.ilogbh(positive_zero_h));
    try std.testing.expectEqual(std.math.minInt(i32), decompose.ilogbh(negative_zero_h));
    try std.testing.expectEqual(-24, decompose.ilogbh(min_subnormal_h));
    try std.testing.expectEqual(-15, decompose.ilogbh(max_subnormal_h));
    try std.testing.expectEqual(-14, decompose.ilogbh(min_normal_h));
    try std.testing.expectEqual(0, decompose.ilogbh(one_h));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbh(positive_infinity_h));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbh(quiet_nan_h));

    try std.testing.expectEqual(std.math.minInt(i32), decompose.ilogbf(positive_zero_f));
    try std.testing.expectEqual(std.math.minInt(i32), decompose.ilogbf(negative_zero_f));
    try std.testing.expectEqual(-149, decompose.ilogbf(min_subnormal_f));
    try std.testing.expectEqual(-127, decompose.ilogbf(max_subnormal_f));
    try std.testing.expectEqual(-126, decompose.ilogbf(min_normal_f));
    try std.testing.expectEqual(0, decompose.ilogbf(one_f));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbf(positive_infinity_f));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbf(quiet_nan_f));

    try std.testing.expectEqual(std.math.minInt(i32), decompose.ilogbd(positive_zero_d));
    try std.testing.expectEqual(std.math.minInt(i32), decompose.ilogbd(negative_zero_d));
    try std.testing.expectEqual(-1074, decompose.ilogbd(min_subnormal_d));
    try std.testing.expectEqual(-1023, decompose.ilogbd(max_subnormal_d));
    try std.testing.expectEqual(-1022, decompose.ilogbd(min_normal_d));
    try std.testing.expectEqual(0, decompose.ilogbd(one_d));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbd(positive_infinity_d));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbd(quiet_nan_d));

    try std.testing.expectEqual(std.math.minInt(i64), decompose.llogbh(positive_zero_h));
    try std.testing.expectEqual(std.math.minInt(i64), decompose.llogbh(negative_zero_h));
    try std.testing.expectEqual(-24, decompose.llogbh(min_subnormal_h));
    try std.testing.expectEqual(-14, decompose.llogbh(min_normal_h));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbh(positive_infinity_h));
    try std.testing.expectEqual(std.math.minInt(i64), decompose.llogbf(positive_zero_f));
    try std.testing.expectEqual(std.math.minInt(i64), decompose.llogbf(negative_zero_f));
    try std.testing.expectEqual(-149, decompose.llogbf(min_subnormal_f));
    try std.testing.expectEqual(-126, decompose.llogbf(min_normal_f));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbf(positive_infinity_f));
    try std.testing.expectEqual(std.math.minInt(i64), decompose.llogbd(positive_zero_d));
    try std.testing.expectEqual(std.math.minInt(i64), decompose.llogbd(negative_zero_d));
    try std.testing.expectEqual(-1074, decompose.llogbd(min_subnormal_d));
    try std.testing.expectEqual(-1022, decompose.llogbd(min_normal_d));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbd(positive_infinity_d));
}

test "ilogb and llogb classify negative infinities and signaling NaNs" {
    const positive_signaling_nan_h: f16 = @bitCast(@as(u16, 0x7d55));
    const negative_quiet_nan_h: f16 = @bitCast(@as(u16, 0xfe55));
    const negative_signaling_nan_h: f16 = @bitCast(@as(u16, 0xfd55));
    const positive_signaling_nan_f: f32 = @bitCast(@as(u32, 0x7f81_2345));
    const negative_quiet_nan_f: f32 = @bitCast(@as(u32, 0xffc1_2345));
    const negative_signaling_nan_f: f32 = @bitCast(@as(u32, 0xff81_2345));
    const positive_signaling_nan_d: f64 = @bitCast(@as(u64, 0x7ff1_2345_6789_abcd));
    const negative_quiet_nan_d: f64 = @bitCast(@as(u64, 0xfff8_1234_5678_9abc));
    const negative_signaling_nan_d: f64 = @bitCast(@as(u64, 0xfff1_2345_6789_abcd));
    const negative_infinity_h: f16 = @bitCast(@as(u16, 0xfc00));
    const negative_infinity_f: f32 = @bitCast(@as(u32, 0xff80_0000));
    const negative_infinity_d: f64 = @bitCast(@as(u64, 0xfff0_0000_0000_0000));

    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbh(negative_infinity_h));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbh(positive_signaling_nan_h));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbh(negative_quiet_nan_h));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbh(negative_signaling_nan_h));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbf(negative_quiet_nan_f));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbf(negative_signaling_nan_f));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbd(negative_quiet_nan_d));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbd(negative_signaling_nan_d));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbh(negative_quiet_nan_h));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbh(negative_signaling_nan_h));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbf(negative_quiet_nan_f));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbf(negative_signaling_nan_f));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbd(negative_quiet_nan_d));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbd(negative_signaling_nan_d));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbf(negative_infinity_f));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbf(positive_signaling_nan_f));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbd(negative_infinity_d));
    try std.testing.expectEqual(std.math.maxInt(i32), decompose.ilogbd(positive_signaling_nan_d));

    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbh(negative_infinity_h));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbh(positive_signaling_nan_h));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbf(negative_infinity_f));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbf(positive_signaling_nan_f));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbd(negative_infinity_d));
    try std.testing.expectEqual(std.math.maxInt(i64), decompose.llogbd(positive_signaling_nan_d));
}

test "logb returns negative infinity for both signed zeros" {
    try expectLogbH(0x0000, 0xfc00);
    try expectLogbH(0x8000, 0xfc00);
    try expectLogbF(0x0000_0000, 0xff80_0000);
    try expectLogbF(0x8000_0000, 0xff80_0000);
    try expectLogbD(0x0000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLogbD(0x8000_0000_0000_0000, 0xfff0_0000_0000_0000);
}

test "logb returns exact finite exponents at all subnormal and normal boundaries" {
    try expectLogbH(0x0001, 0xce00);
    try expectLogbH(0x03ff, 0xcb80);
    try expectLogbH(0x0400, 0xcb00);
    try expectLogbH(0x3c00, 0x0000);
    try expectLogbH(0x7bff, 0x4b80);

    try expectLogbF(0x0000_0001, 0xc315_0000);
    try expectLogbF(0x007f_ffff, 0xc2fe_0000);
    try expectLogbF(0x0080_0000, 0xc2fc_0000);
    try expectLogbF(0x3f80_0000, 0x0000_0000);
    try expectLogbF(0x7f7f_ffff, 0x42fe_0000);

    try expectLogbD(0x0000_0000_0000_0001, 0xc090_c800_0000_0000);
    try expectLogbD(0x000f_ffff_ffff_ffff, 0xc08f_f800_0000_0000);
    try expectLogbD(0x0010_0000_0000_0000, 0xc08f_f000_0000_0000);
    try expectLogbD(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectLogbD(0x7fef_ffff_ffff_ffff, 0x408f_f800_0000_0000);
}

test "logb maps infinities to positive infinity and quiets NaNs" {
    try expectLogbH(0x7c00, 0x7c00);
    try expectLogbH(0xfc00, 0x7c00);
    try expectLogbH(0x7e55, 0x7e55);
    try expectLogbH(0xfe55, 0xfe55);
    try expectLogbH(0x7d55, 0x7f55);
    try expectLogbH(0xfd55, 0xff55);

    try expectLogbF(0x7f80_0000, 0x7f80_0000);
    try expectLogbF(0xff80_0000, 0x7f80_0000);
    try expectLogbF(0x7fc1_2345, 0x7fc1_2345);
    try expectLogbF(0xffc1_2345, 0xffc1_2345);
    try expectLogbF(0x7f81_2345, 0x7fc1_2345);
    try expectLogbF(0xff81_2345, 0xffc1_2345);

    try expectLogbD(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectLogbD(0xfff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectLogbD(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLogbD(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectLogbD(0x7ff1_2345_6789_abcd, 0x7ff9_2345_6789_abcd);
    try expectLogbD(0xfff1_2345_6789_abcd, 0xfff9_2345_6789_abcd);
}
