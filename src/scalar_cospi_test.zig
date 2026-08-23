const std = @import("std");
const root = @import("root.zig");
const scalar = @import("scalar.zig");

fn expectCospihBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.cospih(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectCospifBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.cospif(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectCospidBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.cospid(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "cospi exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.cospih) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.cospif) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.cospid) == fn (f64) f64);
    }
}

test "cospi handles signed zero integers and half-integers exactly" {
    try expectCospihBits(0x0000, 0x3c00);
    try expectCospihBits(0x8000, 0x3c00);
    try expectCospihBits(0x3800, 0x0000);
    try expectCospihBits(0xb800, 0x0000);
    try expectCospihBits(0x3c00, 0xbc00);
    try expectCospihBits(0x4000, 0x3c00);
    try expectCospihBits(0x3e00, 0x0000);

    try expectCospifBits(0x0000_0000, 0x3f80_0000);
    try expectCospifBits(0x8000_0000, 0x3f80_0000);
    try expectCospifBits(0x3f00_0000, 0x0000_0000);
    try expectCospifBits(0xbf00_0000, 0x0000_0000);
    try expectCospifBits(0x3f80_0000, 0xbf80_0000);
    try expectCospifBits(0x4000_0000, 0x3f80_0000);
    try expectCospifBits(0x3fc0_0000, 0x0000_0000);

    try expectCospidBits(0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectCospidBits(0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectCospidBits(0x3fe0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectCospidBits(0xbfe0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectCospidBits(0x3ff0_0000_0000_0000, 0xbff0_0000_0000_0000);
    try expectCospidBits(0x4000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectCospidBits(0x3ff8_0000_0000_0000, 0x0000_0000_0000_0000);
}

test "cospi preserves stable NaN bits and canonicalizes infinity" {
    try expectCospihBits(0x7c00, 0x7e00);
    try expectCospihBits(0xfc00, 0x7e00);
    try expectCospihBits(0x7d01, 0x7f01);
    try expectCospihBits(0xfe55, 0xfe55);

    try expectCospifBits(0x7f80_0000, 0x7fc0_0000);
    try expectCospifBits(0xff80_0000, 0x7fc0_0000);
    try expectCospifBits(0x7f80_1234, 0x7fc0_1234);
    try expectCospifBits(0xffc1_2345, 0xffc1_2345);

    try expectCospidBits(0x7ff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectCospidBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectCospidBits(0x7ff0_0000_0000_1234, 0x7ff8_0000_0000_1234);
    try expectCospidBits(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "cospi rounds quarter turns and large exact inputs" {
    try expectCospihBits(0x3400, 0x39a8);
    try expectCospihBits(0x3a00, 0xb9a8);
    try std.testing.expectEqual(@as(f16, -1.0), scalar.cospih(@as(f16, 1023.0)));
    try std.testing.expectEqual(@as(f16, 1.0), scalar.cospih(@as(f16, 1024.0)));

    try expectCospifBits(0x3e80_0000, 0x3f35_04f3);
    try expectCospifBits(0x3f40_0000, 0xbf35_04f3);
    try std.testing.expectEqual(@as(f32, -1.0), scalar.cospif(@as(f32, 8_388_607.0)));
    try std.testing.expectEqual(@as(f32, 1.0), scalar.cospif(@as(f32, 8_388_608.0)));

    try expectCospidBits(0x3fd0_0000_0000_0000, 0x3fe6_a09e_667f_3bcd);
    try expectCospidBits(0x3fe8_0000_0000_0000, 0xbfe6_a09e_667f_3bcd);
    try std.testing.expectEqual(
        @as(f64, -1.0),
        scalar.cospid(@as(f64, 4_503_599_627_370_495.0)),
    );
    try std.testing.expectEqual(
        @as(f64, 1.0),
        scalar.cospid(@as(f64, 4_503_599_627_370_496.0)),
    );
}

test "root cospi dispatches all supported formats" {
    try std.testing.expectEqual(@as(f16, -1.0), root.cospi(@as(f16, 1.0)));
    try std.testing.expectEqual(@as(f32, -1.0), root.cospi(@as(f32, 1.0)));
    try std.testing.expectEqual(@as(f64, -1.0), root.cospi(@as(f64, 1.0)));
}
