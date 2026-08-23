const std = @import("std");
const scalar = @import("scalar.zig");

fn expectLog2hBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.log2h(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectLog2fBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.log2f(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectLog2dBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.log2d(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "log2 exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.log2h) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.log2f) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.log2d) == fn (f64) f64);
    }
}

test "log2 returns exact powers of two" {
    try expectLog2hBits(0x3800, 0xbc00);
    try expectLog2hBits(0x3c00, 0x0000);
    try expectLog2hBits(0x4000, 0x3c00);
    try expectLog2hBits(0x7000, 0x4a80);

    try expectLog2fBits(0x3f00_0000, 0xbf80_0000);
    try expectLog2fBits(0x3f80_0000, 0x0000_0000);
    try expectLog2fBits(0x4000_0000, 0x3f80_0000);
    try expectLog2fBits(0x7f00_0000, 0x42fe_0000);

    try expectLog2dBits(0x3fe0_0000_0000_0000, 0xbff0_0000_0000_0000);
    try expectLog2dBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectLog2dBits(0x4000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectLog2dBits(0x7fe0_0000_0000_0000, 0x408f_f800_0000_0000);
}

test "log2 handles zeros infinities and negative domain" {
    try expectLog2hBits(0x0000, 0xfc00);
    try expectLog2hBits(0x8000, 0xfc00);
    try expectLog2hBits(0x7c00, 0x7c00);
    try expectLog2hBits(0xbc00, 0x7e00);
    try expectLog2hBits(0xfc00, 0x7e00);

    try expectLog2fBits(0x0000_0000, 0xff80_0000);
    try expectLog2fBits(0x8000_0000, 0xff80_0000);
    try expectLog2fBits(0x7f80_0000, 0x7f80_0000);
    try expectLog2fBits(0xbf80_0000, 0x7fc0_0000);
    try expectLog2fBits(0xff80_0000, 0x7fc0_0000);

    try expectLog2dBits(0x0000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLog2dBits(0x8000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLog2dBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectLog2dBits(0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectLog2dBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
}

test "log2 quiets NaNs while preserving sign and payload" {
    try expectLog2hBits(0x7c55, 0x7e55);
    try expectLog2hBits(0xfc55, 0xfe55);
    try expectLog2fBits(0x7f81_2345, 0x7fc1_2345);
    try expectLog2fBits(0xff81_2345, 0xffc1_2345);
    try expectLog2dBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLog2dBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "log2 handles subnormal boundaries" {
    try expectLog2hBits(0x0001, 0xce00);
    try expectLog2hBits(0x0400, 0xcb00);
    try expectLog2fBits(0x0000_0001, 0xc315_0000);
    try expectLog2fBits(0x0080_0000, 0xc2fc_0000);
    try expectLog2dBits(0x0000_0000_0000_0001, 0xc090_c800_0000_0000);
    try expectLog2dBits(0x0010_0000_0000_0000, 0xc08f_f000_0000_0000);
}

test "log2d resolves CORE-MATH hard cases" {
    try expectLog2dBits(0x3ff9_0564_eb47_015d, 0x3fe4_a46d_4bb4_f1b7);
    try expectLog2dBits(0x3fff_f556_6506_93a6, 0x3feff09b72649a00);
    try expectLog2dBits(0x0000_12ff_4037_3cbb, 0xc090_1702_598b_4b68);
}

test "log2h visits every finite positive encoding" {
    var raw_bits: u32 = 1;
    while (raw_bits < 0x7c00) : (raw_bits += 1) {
        const input_bits: u16 = @intCast(raw_bits);
        const output_bits: u16 = @bitCast(scalar.log2h(@bitCast(input_bits)));
        try std.testing.expect((output_bits & 0x7c00) != 0x7c00);
    }
}
