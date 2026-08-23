const std = @import("std");
const scalar = @import("scalar.zig");
const root = @import("root.zig");

fn expectH(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u16, @bitCast(scalar.log2p1h(input))));
    try std.testing.expectEqual(expected_bits, @as(u16, @bitCast(root.log2p1(input))));
}

fn expectF(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(scalar.log2p1f(input))));
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(root.log2p1(input))));
}

fn expectD(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(scalar.log2p1d(input))));
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(root.log2p1(input))));
}

test "log2p1 exposes nearest-even H, F, and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.log2p1h) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.log2p1f) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.log2p1d) == fn (f64) f64);
    }
}

test "log2p1 preserves zeros and handles exact powers" {
    try expectH(0x0000, 0x0000);
    try expectH(0x8000, 0x8000);
    try expectH(0x3c00, 0x3c00);
    try expectH(0xb800, 0xbc00);
    try expectH(0x4200, 0x4000);

    try expectF(0x0000_0000, 0x0000_0000);
    try expectF(0x8000_0000, 0x8000_0000);
    try expectF(0x3f80_0000, 0x3f80_0000);
    try expectF(0xbf00_0000, 0xbf80_0000);
    try expectF(0x4040_0000, 0x4000_0000);

    try expectD(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectD(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectD(0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectD(0xbfe0_0000_0000_0000, 0xbff0_0000_0000_0000);
    try expectD(0x4008_0000_0000_0000, 0x4000_0000_0000_0000);
}
test "log2p1 covers pinned CORE-MATH H and F hard cases" {
    try expectH(0x8cf6, 0x8f29);
    try expectH(0x9505, 0x973f);
    try expectH(0xa5cf, 0xa83d);
    try expectH(0x8309, 0x8461);
    try expectH(0x04b3, 0x06c7);
    try expectH(0x0b46, 0x0d3f);
    try expectH(0x25ff, 0x2847);
    try expectH(0x0fb1, 0x118c);
    try expectH(0x1ffb, 0x21bc);

    try expectF(0x4ebd_09e3, 0x41f4_8013);
    try expectF(0xbd6d_142e, 0xbdb0_2a48);
}

test "log2p1 applies stable pole domain infinity and NaN policy" {
    try expectH(0xbc00, 0xfc00);
    try expectH(0xc000, 0x7e00);
    try expectH(0xfc00, 0x7e00);
    try expectH(0x7c00, 0x7c00);
    try expectH(0x7c55, 0x7e55);
    try expectH(0xfe55, 0xfe55);

    try expectF(0xbf80_0000, 0xff80_0000);
    try expectF(0xc000_0000, 0x7fc0_0000);
    try expectF(0xff80_0000, 0x7fc0_0000);
    try expectF(0x7f80_0000, 0x7f80_0000);
    try expectF(0x7f81_2345, 0x7fc1_2345);
    try expectF(0xffc1_2345, 0xffc1_2345);

    try expectD(0xbff0_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectD(0xc000_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectD(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectD(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectD(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectD(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}
