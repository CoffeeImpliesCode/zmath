const std = @import("std");
const common = @import("scalar_c23_common.zig");
const rounding = @import("scalar_rounding.zig");

// Each case stores input, ceil, floor, trunc, ties-even, and ties-away bits.
fn expectH(entry: [6]u16) !void {
    const value = common.H.fromBits(entry[0]);
    try std.testing.expectEqual(entry[1], common.H.toBits(rounding.ceilh(value)));
    try std.testing.expectEqual(entry[2], common.H.toBits(rounding.floorh(value)));
    try std.testing.expectEqual(entry[3], common.H.toBits(rounding.trunch(value)));
    try std.testing.expectEqual(entry[4], common.H.toBits(rounding.nearbyinth(value)));
    try std.testing.expectEqual(entry[4], common.H.toBits(rounding.rinth(value)));
    try std.testing.expectEqual(entry[5], common.H.toBits(rounding.roundh(value)));
    try std.testing.expectEqual(entry[4], common.H.toBits(rounding.roundevenh(value)));
}

fn expectF(entry: [6]u32) !void {
    const value = common.F.fromBits(entry[0]);
    try std.testing.expectEqual(entry[1], common.F.toBits(rounding.ceilf(value)));
    try std.testing.expectEqual(entry[2], common.F.toBits(rounding.floorf(value)));
    try std.testing.expectEqual(entry[3], common.F.toBits(rounding.truncf(value)));
    try std.testing.expectEqual(entry[4], common.F.toBits(rounding.nearbyintf(value)));
    try std.testing.expectEqual(entry[4], common.F.toBits(rounding.rintf(value)));
    try std.testing.expectEqual(entry[5], common.F.toBits(rounding.roundf(value)));
    try std.testing.expectEqual(entry[4], common.F.toBits(rounding.roundevenf(value)));
}

fn expectD(entry: [6]u64) !void {
    const value = common.D.fromBits(entry[0]);
    try std.testing.expectEqual(entry[1], common.D.toBits(rounding.ceild(value)));
    try std.testing.expectEqual(entry[2], common.D.toBits(rounding.floord(value)));
    try std.testing.expectEqual(entry[3], common.D.toBits(rounding.truncd(value)));
    try std.testing.expectEqual(entry[4], common.D.toBits(rounding.nearbyintd(value)));
    try std.testing.expectEqual(entry[4], common.D.toBits(rounding.rintd(value)));
    try std.testing.expectEqual(entry[5], common.D.toBits(rounding.roundd(value)));
    try std.testing.expectEqual(entry[4], common.D.toBits(rounding.roundevend(value)));
}

test "f16 rounding uses exact target-format bits at all named boundaries" {
    const cases = [_][6]u16{
        .{ 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000 }, // +zero
        .{ 0x8000, 0x8000, 0x8000, 0x8000, 0x8000, 0x8000 }, // -zero
        .{ 0x0001, 0x3c00, 0x0000, 0x0000, 0x0000, 0x0000 }, // Minimum +subnormal.
        .{ 0x8001, 0x8000, 0xbc00, 0x8000, 0x8000, 0x8000 }, // Minimum -subnormal.
        .{ 0x03ff, 0x3c00, 0x0000, 0x0000, 0x0000, 0x0000 }, // Maximum +subnormal.
        .{ 0x83ff, 0x8000, 0xbc00, 0x8000, 0x8000, 0x8000 }, // Maximum -subnormal.
        .{ 0x0400, 0x3c00, 0x0000, 0x0000, 0x0000, 0x0000 }, // Minimum +normal.
        .{ 0x8400, 0x8000, 0xbc00, 0x8000, 0x8000, 0x8000 }, // Minimum -normal.
        .{ 0x37ff, 0x3c00, 0x0000, 0x0000, 0x0000, 0x0000 }, // Below +0.5.
        .{ 0x3800, 0x3c00, 0x0000, 0x0000, 0x0000, 0x3c00 }, // +0.5.
        .{ 0x3801, 0x3c00, 0x0000, 0x0000, 0x3c00, 0x3c00 }, // Above +0.5.
        .{ 0xb7ff, 0x8000, 0xbc00, 0x8000, 0x8000, 0x8000 }, // Above -0.5.
        .{ 0xb800, 0x8000, 0xbc00, 0x8000, 0x8000, 0xbc00 }, // -0.5.
        .{ 0xb801, 0x8000, 0xbc00, 0x8000, 0xbc00, 0xbc00 }, // Below -0.5.
        .{ 0x3bff, 0x3c00, 0x0000, 0x0000, 0x3c00, 0x3c00 }, // Below +1.
        .{ 0x3c00, 0x3c00, 0x3c00, 0x3c00, 0x3c00, 0x3c00 }, // +1.
        .{ 0x3c01, 0x4000, 0x3c00, 0x3c00, 0x3c00, 0x3c00 }, // Above +1.
        .{ 0xbbff, 0x8000, 0xbc00, 0x8000, 0xbc00, 0xbc00 }, // Above -1.
        .{ 0xbc00, 0xbc00, 0xbc00, 0xbc00, 0xbc00, 0xbc00 }, // -1.
        .{ 0xbc01, 0xbc00, 0xc000, 0xbc00, 0xbc00, 0xbc00 }, // Below -1.
        .{ 0x3e00, 0x4000, 0x3c00, 0x3c00, 0x4000, 0x4000 }, // +1.5.
        .{ 0xbe00, 0xbc00, 0xc000, 0xbc00, 0xc000, 0xc000 }, // -1.5.
        .{ 0x4100, 0x4200, 0x4000, 0x4000, 0x4000, 0x4200 }, // +2.5.
        .{ 0xc100, 0xc000, 0xc200, 0xc000, 0xc000, 0xc200 }, // -2.5.
        .{ 0x4300, 0x4400, 0x4200, 0x4200, 0x4400, 0x4400 }, // +3.5.
        .{ 0xc300, 0xc200, 0xc400, 0xc200, 0xc400, 0xc400 }, // -3.5.
        .{ 0x63ff, 0x6400, 0x63fe, 0x63fe, 0x6400, 0x6400 }, // Below 2^10.
        .{ 0xe3ff, 0xe3fe, 0xe400, 0xe3fe, 0xe400, 0xe400 }, // Above -2^10.
        .{ 0x6400, 0x6400, 0x6400, 0x6400, 0x6400, 0x6400 }, // First integral binade.
        .{ 0x6401, 0x6401, 0x6401, 0x6401, 0x6401, 0x6401 }, // Next integral value.
        .{ 0x7bff, 0x7bff, 0x7bff, 0x7bff, 0x7bff, 0x7bff }, // Maximum finite.
        .{ 0xfbff, 0xfbff, 0xfbff, 0xfbff, 0xfbff, 0xfbff }, // Minimum finite.
        .{ 0x7c00, 0x7c00, 0x7c00, 0x7c00, 0x7c00, 0x7c00 }, // +infinity.
        .{ 0xfc00, 0xfc00, 0xfc00, 0xfc00, 0xfc00, 0xfc00 }, // -infinity.
        .{ 0x7c01, 0x7e01, 0x7e01, 0x7e01, 0x7e01, 0x7e01 }, // Minimum +sNaN.
        .{ 0xfc01, 0xfe01, 0xfe01, 0xfe01, 0xfe01, 0xfe01 }, // Minimum -sNaN.
        .{ 0x7d55, 0x7f55, 0x7f55, 0x7f55, 0x7f55, 0x7f55 }, // +sNaN payload.
        .{ 0xfd55, 0xff55, 0xff55, 0xff55, 0xff55, 0xff55 }, // -sNaN payload.
        .{ 0x7e55, 0x7e55, 0x7e55, 0x7e55, 0x7e55, 0x7e55 }, // +qNaN payload.
        .{ 0xfe55, 0xfe55, 0xfe55, 0xfe55, 0xfe55, 0xfe55 }, // -qNaN payload.
        .{ 0x7fff, 0x7fff, 0x7fff, 0x7fff, 0x7fff, 0x7fff }, // Maximum +qNaN.
        .{ 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff }, // Maximum -qNaN.
    };

    for (cases) |entry| try expectH(entry);
}

test "f32 rounding covers ties, format boundaries, and special values" {
    const cases = [_][6]u32{
        .{ 0x0000_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000 },
        .{ 0x8000_0000, 0x8000_0000, 0x8000_0000, 0x8000_0000, 0x8000_0000, 0x8000_0000 },
        .{ 0x0000_0001, 0x3f80_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000 },
        .{ 0x8000_0001, 0x8000_0000, 0xbf80_0000, 0x8000_0000, 0x8000_0000, 0x8000_0000 },
        .{ 0x007f_ffff, 0x3f80_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000 },
        .{ 0x807f_ffff, 0x8000_0000, 0xbf80_0000, 0x8000_0000, 0x8000_0000, 0x8000_0000 },
        .{ 0x0080_0000, 0x3f80_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000 },
        .{ 0x8080_0000, 0x8000_0000, 0xbf80_0000, 0x8000_0000, 0x8000_0000, 0x8000_0000 },
        .{ 0x3eff_ffff, 0x3f80_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000 },
        .{ 0x3f00_0000, 0x3f80_0000, 0x0000_0000, 0x0000_0000, 0x0000_0000, 0x3f80_0000 },
        .{ 0x3f00_0001, 0x3f80_0000, 0x0000_0000, 0x0000_0000, 0x3f80_0000, 0x3f80_0000 },
        .{ 0xbeff_ffff, 0x8000_0000, 0xbf80_0000, 0x8000_0000, 0x8000_0000, 0x8000_0000 },
        .{ 0xbf00_0000, 0x8000_0000, 0xbf80_0000, 0x8000_0000, 0x8000_0000, 0xbf80_0000 },
        .{ 0xbf00_0001, 0x8000_0000, 0xbf80_0000, 0x8000_0000, 0xbf80_0000, 0xbf80_0000 },
        .{ 0x3f7f_ffff, 0x3f80_0000, 0x0000_0000, 0x0000_0000, 0x3f80_0000, 0x3f80_0000 },
        .{ 0x3f80_0000, 0x3f80_0000, 0x3f80_0000, 0x3f80_0000, 0x3f80_0000, 0x3f80_0000 },
        .{ 0x3f80_0001, 0x4000_0000, 0x3f80_0000, 0x3f80_0000, 0x3f80_0000, 0x3f80_0000 },
        .{ 0xbf7f_ffff, 0x8000_0000, 0xbf80_0000, 0x8000_0000, 0xbf80_0000, 0xbf80_0000 },
        .{ 0xbf80_0000, 0xbf80_0000, 0xbf80_0000, 0xbf80_0000, 0xbf80_0000, 0xbf80_0000 },
        .{ 0xbf80_0001, 0xbf80_0000, 0xc000_0000, 0xbf80_0000, 0xbf80_0000, 0xbf80_0000 },
        .{ 0x3fc0_0000, 0x4000_0000, 0x3f80_0000, 0x3f80_0000, 0x4000_0000, 0x4000_0000 },
        .{ 0xbfc0_0000, 0xbf80_0000, 0xc000_0000, 0xbf80_0000, 0xc000_0000, 0xc000_0000 },
        .{ 0x4020_0000, 0x4040_0000, 0x4000_0000, 0x4000_0000, 0x4000_0000, 0x4040_0000 },
        .{ 0xc020_0000, 0xc000_0000, 0xc040_0000, 0xc000_0000, 0xc000_0000, 0xc040_0000 },
        .{ 0x4aff_ffff, 0x4b00_0000, 0x4aff_fffe, 0x4aff_fffe, 0x4b00_0000, 0x4b00_0000 },
        .{ 0xcaff_ffff, 0xcaff_fffe, 0xcb00_0000, 0xcaff_fffe, 0xcb00_0000, 0xcb00_0000 },
        .{ 0x4b00_0000, 0x4b00_0000, 0x4b00_0000, 0x4b00_0000, 0x4b00_0000, 0x4b00_0000 },
        .{ 0x4b00_0001, 0x4b00_0001, 0x4b00_0001, 0x4b00_0001, 0x4b00_0001, 0x4b00_0001 },
        .{ 0x7f7f_ffff, 0x7f7f_ffff, 0x7f7f_ffff, 0x7f7f_ffff, 0x7f7f_ffff, 0x7f7f_ffff },
        .{ 0x7f80_0000, 0x7f80_0000, 0x7f80_0000, 0x7f80_0000, 0x7f80_0000, 0x7f80_0000 },
        .{ 0xff80_0000, 0xff80_0000, 0xff80_0000, 0xff80_0000, 0xff80_0000, 0xff80_0000 },
        .{ 0xff7f_ffff, 0xff7f_ffff, 0xff7f_ffff, 0xff7f_ffff, 0xff7f_ffff, 0xff7f_ffff },
        .{ 0x7f80_0001, 0x7fc0_0001, 0x7fc0_0001, 0x7fc0_0001, 0x7fc0_0001, 0x7fc0_0001 },
        .{ 0xff80_0001, 0xffc0_0001, 0xffc0_0001, 0xffc0_0001, 0xffc0_0001, 0xffc0_0001 },
        .{ 0x7f81_2345, 0x7fc1_2345, 0x7fc1_2345, 0x7fc1_2345, 0x7fc1_2345, 0x7fc1_2345 },
        .{ 0xff81_2345, 0xffc1_2345, 0xffc1_2345, 0xffc1_2345, 0xffc1_2345, 0xffc1_2345 },
        .{ 0x7fc1_2345, 0x7fc1_2345, 0x7fc1_2345, 0x7fc1_2345, 0x7fc1_2345, 0x7fc1_2345 },
        .{ 0xffc1_2345, 0xffc1_2345, 0xffc1_2345, 0xffc1_2345, 0xffc1_2345, 0xffc1_2345 },
    };

    for (cases) |entry| try expectF(entry);
}

test "f64 rounding covers ties, format boundaries, and special values" {
    const cases = [_][6]u64{
        .{
            0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000,
            0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000,
        },
        .{
            0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000,
            0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000,
        },
        .{
            0x0000_0000_0000_0001, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000,
            0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000,
        },
        .{
            0x8000_0000_0000_0001, 0x8000_0000_0000_0000, 0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000,
        },
        .{
            0x000f_ffff_ffff_ffff, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000,
            0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000,
        },
        .{
            0x800f_ffff_ffff_ffff, 0x8000_0000_0000_0000, 0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000,
        },
        .{
            0x0010_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000,
            0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000,
        },
        .{
            0x8010_0000_0000_0000, 0x8000_0000_0000_0000, 0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000,
        },
        .{
            0x3fdf_ffff_ffff_ffff, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000,
            0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x0000_0000_0000_0000,
        },
        .{
            0x3fe0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000,
            0x0000_0000_0000_0000, 0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000,
        },
        .{
            0x3fe0_0000_0000_0001, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000,
            0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000,
        },
        .{
            0xbfdf_ffff_ffff_ffff, 0x8000_0000_0000_0000, 0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0x8000_0000_0000_0000,
        },
        .{
            0xbfe0_0000_0000_0000, 0x8000_0000_0000_0000, 0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000, 0x8000_0000_0000_0000, 0xbff0_0000_0000_0000,
        },
        .{
            0xbfe0_0000_0000_0001, 0x8000_0000_0000_0000, 0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000, 0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000,
        },
        .{
            0x3fef_ffff_ffff_ffff, 0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000,
            0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000,
        },
        .{
            0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000,
            0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000,
        },
        .{
            0x3ff0_0000_0000_0001, 0x4000_0000_0000_0000, 0x3ff0_0000_0000_0000,
            0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000,
        },
        .{
            0xbfef_ffff_ffff_ffff, 0x8000_0000_0000_0000, 0xbff0_0000_0000_0000,
            0x8000_0000_0000_0000, 0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000,
        },
        .{
            0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000,
            0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000,
        },
        .{
            0xbff0_0000_0000_0001, 0xbff0_0000_0000_0000, 0xc000_0000_0000_0000,
            0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000, 0xbff0_0000_0000_0000,
        },
        .{
            0x3ff8_0000_0000_0000, 0x4000_0000_0000_0000, 0x3ff0_0000_0000_0000,
            0x3ff0_0000_0000_0000, 0x4000_0000_0000_0000, 0x4000_0000_0000_0000,
        },
        .{
            0xbff8_0000_0000_0000, 0xbff0_0000_0000_0000, 0xc000_0000_0000_0000,
            0xbff0_0000_0000_0000, 0xc000_0000_0000_0000, 0xc000_0000_0000_0000,
        },
        .{
            0x4004_0000_0000_0000, 0x4008_0000_0000_0000, 0x4000_0000_0000_0000,
            0x4000_0000_0000_0000, 0x4000_0000_0000_0000, 0x4008_0000_0000_0000,
        },
        .{
            0xc004_0000_0000_0000, 0xc000_0000_0000_0000, 0xc008_0000_0000_0000,
            0xc000_0000_0000_0000, 0xc000_0000_0000_0000, 0xc008_0000_0000_0000,
        },
        .{
            0x432f_ffff_ffff_ffff, 0x4330_0000_0000_0000, 0x432f_ffff_ffff_fffe,
            0x432f_ffff_ffff_fffe, 0x4330_0000_0000_0000, 0x4330_0000_0000_0000,
        },
        .{
            0xc32f_ffff_ffff_ffff, 0xc32f_ffff_ffff_fffe, 0xc330_0000_0000_0000,
            0xc32f_ffff_ffff_fffe, 0xc330_0000_0000_0000, 0xc330_0000_0000_0000,
        },
        .{
            0x4330_0000_0000_0000, 0x4330_0000_0000_0000, 0x4330_0000_0000_0000,
            0x4330_0000_0000_0000, 0x4330_0000_0000_0000, 0x4330_0000_0000_0000,
        },
        .{
            0x4330_0000_0000_0001, 0x4330_0000_0000_0001, 0x4330_0000_0000_0001,
            0x4330_0000_0000_0001, 0x4330_0000_0000_0001, 0x4330_0000_0000_0001,
        },
        .{
            0x7fef_ffff_ffff_ffff, 0x7fef_ffff_ffff_ffff, 0x7fef_ffff_ffff_ffff,
            0x7fef_ffff_ffff_ffff, 0x7fef_ffff_ffff_ffff, 0x7fef_ffff_ffff_ffff,
        },
        .{
            0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000,
            0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000,
        },
        .{
            0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000,
            0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000,
        },
        .{
            0xffef_ffff_ffff_ffff, 0xffef_ffff_ffff_ffff, 0xffef_ffff_ffff_ffff,
            0xffef_ffff_ffff_ffff, 0xffef_ffff_ffff_ffff, 0xffef_ffff_ffff_ffff,
        },
        .{
            0x7ff0_0000_0000_0001, 0x7ff8_0000_0000_0001, 0x7ff8_0000_0000_0001,
            0x7ff8_0000_0000_0001, 0x7ff8_0000_0000_0001, 0x7ff8_0000_0000_0001,
        },
        .{
            0xfff0_0000_0000_0001, 0xfff8_0000_0000_0001, 0xfff8_0000_0000_0001,
            0xfff8_0000_0000_0001, 0xfff8_0000_0000_0001, 0xfff8_0000_0000_0001,
        },
        .{
            0x7ff0_0000_0001_2345, 0x7ff8_0000_0001_2345, 0x7ff8_0000_0001_2345,
            0x7ff8_0000_0001_2345, 0x7ff8_0000_0001_2345, 0x7ff8_0000_0001_2345,
        },
        .{
            0xfff0_0000_0001_2345, 0xfff8_0000_0001_2345, 0xfff8_0000_0001_2345,
            0xfff8_0000_0001_2345, 0xfff8_0000_0001_2345, 0xfff8_0000_0001_2345,
        },
        .{
            0x7ff8_0000_0001_2345, 0x7ff8_0000_0001_2345, 0x7ff8_0000_0001_2345,
            0x7ff8_0000_0001_2345, 0x7ff8_0000_0001_2345, 0x7ff8_0000_0001_2345,
        },
        .{
            0xfff8_0000_0001_2345, 0xfff8_0000_0001_2345, 0xfff8_0000_0001_2345,
            0xfff8_0000_0001_2345, 0xfff8_0000_0001_2345, 0xfff8_0000_0001_2345,
        },
    };

    for (cases) |entry| try expectD(entry);
}

const OracleMode = enum {
    nearest_even,
    nearest_from_zero,
    toward_zero,
    toward_positive,
    toward_negative,
};

fn encodeHInteger(sign: u16, magnitude: u32) u16 {
    if (magnitude == 0) return sign;

    const exponent: u5 = @intCast(31 - @clz(magnitude));
    std.debug.assert(exponent <= 10);
    const leading = @as(u32, 1) << exponent;
    const fraction_shift: u5 = @intCast(10 - exponent);
    const fraction: u16 = @intCast((magnitude - leading) << fraction_shift);
    const exponent_field: u16 = @intCast(@as(u32, exponent) + 15);
    return sign | (exponent_field << 10) | fraction;
}

// The oracle decodes an exact dyadic rational, rounds its integer quotient,
// and encodes the integer again. It does not share the production mask-and-add path.
fn oracleH(bits: u16, mode: OracleMode) u16 {
    const sign = bits & 0x8000;
    const magnitude = bits & 0x7fff;
    const exponent_field = (magnitude >> 10) & 0x1f;
    const fraction = magnitude & 0x03ff;

    if (exponent_field == 0x1f) {
        if (fraction == 0) return bits;
        return bits | 0x0200;
    }
    if (magnitude == 0) return bits;

    const significand: u32 = if (exponent_field == 0)
        fraction
    else
        @as(u32, 0x0400 | fraction);
    const shift: i32 = if (exponent_field == 0)
        24
    else
        25 - @as(i32, exponent_field);
    if (shift <= 0) return bits;

    const shift_amount: u5 = @intCast(shift);
    const denominator = @as(u32, 1) << shift_amount;
    const quotient = significand >> shift_amount;
    const remainder = significand & (denominator - 1);
    const doubled_remainder = remainder * 2;
    const increment = switch (mode) {
        .nearest_even => doubled_remainder > denominator or
            (doubled_remainder == denominator and quotient & 1 != 0),
        .nearest_from_zero => doubled_remainder >= denominator,
        .toward_zero => false,
        .toward_positive => sign == 0 and remainder != 0,
        .toward_negative => sign != 0 and remainder != 0,
    };
    const rounded = if (increment) quotient + 1 else quotient;
    return encodeHInteger(sign, rounded);
}

test "all f16 encodings agree with an exact dyadic oracle" {
    // Exhaustion covers every zero, subnormal, normal-binade boundary,
    // infinity, signaling-NaN, quiet-NaN, payload, and sign encoding.
    var raw: u32 = 0;
    while (raw < 1 << 16) : (raw += 1) {
        const bits: u16 = @intCast(raw);
        const value = common.H.fromBits(bits);

        try std.testing.expectEqual(
            oracleH(bits, .toward_positive),
            common.H.toBits(rounding.ceilh(value)),
        );
        try std.testing.expectEqual(
            oracleH(bits, .toward_negative),
            common.H.toBits(rounding.floorh(value)),
        );
        try std.testing.expectEqual(
            oracleH(bits, .toward_zero),
            common.H.toBits(rounding.trunch(value)),
        );
        try std.testing.expectEqual(
            oracleH(bits, .nearest_even),
            common.H.toBits(rounding.nearbyinth(value)),
        );
        try std.testing.expectEqual(
            oracleH(bits, .nearest_even),
            common.H.toBits(rounding.rinth(value)),
        );
        try std.testing.expectEqual(
            oracleH(bits, .nearest_from_zero),
            common.H.toBits(rounding.roundh(value)),
        );
        try std.testing.expectEqual(
            oracleH(bits, .nearest_even),
            common.H.toBits(rounding.roundevenh(value)),
        );
    }
}
