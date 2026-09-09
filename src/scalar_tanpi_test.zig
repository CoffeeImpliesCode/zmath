const std = @import("std");
const tanpi = @import("scalar_tanpi.zig");

fn expectBits(comptime T: type, expected: @Int(.unsigned, @bitSizeOf(T)), actual: T) !void {
    try std.testing.expectEqual(expected, @as(@Int(.unsigned, @bitSizeOf(T)), @bitCast(actual)));
}

test "Tanpi H preserves exact C23 special values" {
    try expectBits(f16, 0x0000, tanpi.tanpih(@bitCast(@as(u16, 0x0000))));
    try expectBits(f16, 0x8000, tanpi.tanpih(@bitCast(@as(u16, 0x8000))));
    try expectBits(f16, 0x8000, tanpi.tanpih(@as(f16, 1.0)));
    try expectBits(f16, 0x0000, tanpi.tanpih(@as(f16, -1.0)));
    try expectBits(f16, 0x7c00, tanpi.tanpih(@as(f16, 0.5)));
    try expectBits(f16, 0xfc00, tanpi.tanpih(@as(f16, -0.5)));
    try expectBits(f16, 0xfc00, tanpi.tanpih(@as(f16, 1.5)));
    try expectBits(f16, 0x7c00, tanpi.tanpih(@as(f16, -1.5)));
    try expectBits(f16, 0x7e00, tanpi.tanpih(@bitCast(@as(u16, 0x7c00))));
    try expectBits(f16, 0x7e00, tanpi.tanpih(@bitCast(@as(u16, 0xfc00))));
    try expectBits(f16, 0x7f55, tanpi.tanpih(@bitCast(@as(u16, 0x7d55))));
    try expectBits(f16, 0xff55, tanpi.tanpih(@bitCast(@as(u16, 0xff55))));
}

test "Tanpi F preserves exact C23 special values" {
    try expectBits(f32, 0x0000_0000, tanpi.tanpif(@bitCast(@as(u32, 0x0000_0000))));
    try expectBits(f32, 0x8000_0000, tanpi.tanpif(@bitCast(@as(u32, 0x8000_0000))));
    try expectBits(f32, 0x8000_0000, tanpi.tanpif(@as(f32, 1.0)));
    try expectBits(f32, 0x0000_0000, tanpi.tanpif(@as(f32, -1.0)));
    try expectBits(f32, 0x7f80_0000, tanpi.tanpif(@as(f32, 0.5)));
    try expectBits(f32, 0xff80_0000, tanpi.tanpif(@as(f32, -0.5)));
    try expectBits(f32, 0x7fc0_0000, tanpi.tanpif(@bitCast(@as(u32, 0x7f80_0000))));
    try expectBits(f32, 0x7fc0_0000, tanpi.tanpif(@bitCast(@as(u32, 0xff80_0000))));
    try expectBits(f32, 0x7fe1_2345, tanpi.tanpif(@bitCast(@as(u32, 0x7fa1_2345))));
    try expectBits(f32, 0xffe1_2345, tanpi.tanpif(@bitCast(@as(u32, 0xffe1_2345))));
}

test "Tanpi D preserves exact C23 special values" {
    try expectBits(f64, 0x0000_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0x0000_0000_0000_0000))));
    try expectBits(f64, 0x8000_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0x8000_0000_0000_0000))));
    try expectBits(f64, 0x8000_0000_0000_0000, tanpi.tanpid(@as(f64, 1.0)));
    try expectBits(f64, 0x0000_0000_0000_0000, tanpi.tanpid(@as(f64, -1.0)));
    try expectBits(f64, 0x7ff0_0000_0000_0000, tanpi.tanpid(@as(f64, 0.5)));
    try expectBits(f64, 0xfff0_0000_0000_0000, tanpi.tanpid(@as(f64, -0.5)));
    try expectBits(f64, 0x7ff8_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0x7ff0_0000_0000_0000))));
    try expectBits(f64, 0x7ff8_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0xfff0_0000_0000_0000))));
    try expectBits(f64, 0x7ff9_2345_6789_abcd, tanpi.tanpid(@bitCast(@as(u64, 0x7ff1_2345_6789_abcd))));
    try expectBits(f64, 0xfff9_2345_6789_abcd, tanpi.tanpid(@bitCast(@as(u64, 0xfff9_2345_6789_abcd))));
}

test "Tanpi H F and D round representative finite values" {
    try expectBits(f16, 0x3c00, tanpi.tanpih(@as(f16, 0.25)));
    try expectBits(f32, 0x3f80_0000, tanpi.tanpif(@as(f32, 0.25)));
    try expectBits(f64, 0x3ff0_0000_0000_0000, tanpi.tanpid(@as(f64, 0.25)));

    try expectBits(f16, 0x36a1, tanpi.tanpih(@as(f16, 0.125)));
    try expectBits(f32, 0x3ed4_13cd, tanpi.tanpif(@as(f32, 0.125)));
    try expectBits(f64, 0x3fda_8279_99fc_ef32, tanpi.tanpid(@as(f64, 0.125)));

    try expectBits(f16, 0xbc00, tanpi.tanpih(@as(f16, 0.75)));
    try expectBits(f32, 0xbf80_0000, tanpi.tanpif(@as(f32, 0.75)));
    try expectBits(f64, 0xbff0_0000_0000_0000, tanpi.tanpid(@as(f64, 0.75)));
}

test "Tanpi F uses the exponent field for the large integral cutoff" {
    try expectBits(f32, 0x8000_0000, tanpi.tanpif(@bitCast(@as(u32, 0x4b00_0001))));
    try expectBits(f32, 0x0000_0000, tanpi.tanpif(@bitCast(@as(u32, 0xcb00_0001))));
}

test "Tanpi D reflects the large-argument correction-table index" {
    try expectBits(
        f64,
        0xc034_5aff_ed20_1b55,
        tanpi.tanpid(@bitCast(@as(u64, 0x42d0_0000_0000_0021))),
    );
    try expectBits(
        f64,
        0x4034_5aff_ed20_1b55,
        tanpi.tanpid(@bitCast(@as(u64, 0xc2d0_0000_0000_0021))),
    );
}

test "Tanpi preserves exact quarter and half boundaries at wide exponents" {
    try expectBits(f32, 0xff80_0000, tanpi.tanpif(@bitCast(@as(u32, 0x4a00_0006))));
    try expectBits(f32, 0x7f80_0000, tanpi.tanpif(@bitCast(@as(u32, 0xca00_0006))));
    try expectBits(f32, 0x3f80_0000, tanpi.tanpif(@bitCast(@as(u32, 0x4a00_0001))));
    try expectBits(f32, 0xbf80_0000, tanpi.tanpif(@bitCast(@as(u32, 0xca00_0001))));

    try expectBits(f64, 0x0000_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0x42d0_0000_0000_0000))));
    try expectBits(f64, 0x3ff0_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0x42d0_0000_0000_0010))));
    try expectBits(f64, 0x7ff0_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0x42d0_0000_0000_0020))));
    try expectBits(f64, 0xbff0_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0x42d0_0000_0000_0030))));
    try expectBits(f64, 0x8000_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0x42d0_0000_0000_0040))));
    try expectBits(f64, 0xbff0_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0xc2d0_0000_0000_0010))));
    try expectBits(f64, 0xfff0_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0xc2d0_0000_0000_0020))));
    try expectBits(f64, 0x3ff0_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0xc2d0_0000_0000_0030))));
    try expectBits(f64, 0x0000_0000_0000_0000, tanpi.tanpid(@bitCast(@as(u64, 0xc2d0_0000_0000_0040))));
}

test "Tanpi preserves the tiny binary64 product boundaries" {
    try expectBits(
        f64,
        0x0000_0000_0000_0003,
        tanpi.tanpid(@bitCast(@as(u64, 0x0000_0000_0000_0001))),
    );
    try expectBits(
        f64,
        0x8000_0000_0000_0003,
        tanpi.tanpid(@bitCast(@as(u64, 0x8000_0000_0000_0001))),
    );
    try expectBits(
        f64,
        0x0029_21fb_5444_2d18,
        tanpi.tanpid(@bitCast(@as(u64, 0x0010_0000_0000_0000))),
    );
    try expectBits(
        f64,
        0x8029_21fb_5444_2d18,
        tanpi.tanpid(@bitCast(@as(u64, 0x8010_0000_0000_0000))),
    );
}
