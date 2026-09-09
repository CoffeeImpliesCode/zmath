const std = @import("std");
const scalar = @import("scalar.zig");

fn expectPowhBits(base_bits: u16, exponent_bits: u16, expected_bits: u16) !void {
    const base: f16 = @bitCast(base_bits);
    const exponent: f16 = @bitCast(exponent_bits);
    const actual_bits: u16 = @bitCast(scalar.powh(base, exponent));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectPowfBits(base_bits: u32, exponent_bits: u32, expected_bits: u32) !void {
    const base: f32 = @bitCast(base_bits);
    const exponent: f32 = @bitCast(exponent_bits);
    const actual_bits: u32 = @bitCast(scalar.powf(base, exponent));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectPowdBits(base_bits: u64, exponent_bits: u64, expected_bits: u64) !void {
    const base: f64 = @bitCast(base_bits);
    const exponent: f64 = @bitCast(exponent_bits);
    const actual_bits: u64 = @bitCast(scalar.powd(base, exponent));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "pow exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.powh) == fn (f16, f16) f16);
        std.debug.assert(@TypeOf(scalar.powf) == fn (f32, f32) f32);
        std.debug.assert(@TypeOf(scalar.powd) == fn (f64, f64) f64);
    }
}

test "pow computes exact positive and negative integral powers" {
    try expectPowhBits(0x4000, 0x4200, 0x4800);
    try expectPowhBits(0xc000, 0x4200, 0xc800);
    try expectPowhBits(0xc000, 0x4400, 0x4c00);
    try expectPowhBits(0x4400, 0x3800, 0x4000);

    try expectPowfBits(0x4000_0000, 0x4040_0000, 0x4100_0000);
    try expectPowfBits(0xc000_0000, 0x4040_0000, 0xc100_0000);
    try expectPowfBits(0xc000_0000, 0x4080_0000, 0x4180_0000);
    try expectPowfBits(0x4080_0000, 0x3f00_0000, 0x4000_0000);

    try expectPowdBits(0x4000_0000_0000_0000, 0x4008_0000_0000_0000, 0x4020_0000_0000_0000);
    try expectPowdBits(0xc000_0000_0000_0000, 0x4008_0000_0000_0000, 0xc020_0000_0000_0000);
    try expectPowdBits(0xc000_0000_0000_0000, 0x4010_0000_0000_0000, 0x4030_0000_0000_0000);
    try expectPowdBits(0x4010_0000_0000_0000, 0x3fe0_0000_0000_0000, 0x4000_0000_0000_0000);
}

test "pow recognizes odd integers at unit spacing" {
    try expectPowhBits(0xbc00, 0x6401, 0xbc00);
    try expectPowhBits(0x8000, 0x6401, 0x8000);
    try expectPowhBits(0xfc00, 0x6401, 0xfc00);

    try expectPowfBits(0xbf80_0000, 0x4b00_0001, 0xbf80_0000);
    try expectPowfBits(0x8000_0000, 0x4b00_0001, 0x8000_0000);
    try expectPowfBits(0xff80_0000, 0x4b00_0001, 0xff80_0000);

    try expectPowdBits(0xbff0_0000_0000_0000, 0x4330_0000_0000_0001, 0xbff0_0000_0000_0000);
    try expectPowdBits(0x8000_0000_0000_0000, 0x4330_0000_0000_0001, 0x8000_0000_0000_0000);
    try expectPowdBits(0xfff0_0000_0000_0000, 0x4330_0000_0000_0001, 0xfff0_0000_0000_0000);
}

test "pow applies C23 zero pole and signed-zero rules" {
    try expectPowhBits(0x8000, 0x4200, 0x8000);
    try expectPowhBits(0x8000, 0xc200, 0xfc00);
    try expectPowhBits(0x8000, 0x4400, 0x0000);
    try expectPowhBits(0x0000, 0xc400, 0x7c00);

    try expectPowfBits(0x8000_0000, 0x4040_0000, 0x8000_0000);
    try expectPowfBits(0x8000_0000, 0xc040_0000, 0xff80_0000);
    try expectPowfBits(0x8000_0000, 0x4080_0000, 0x0000_0000);
    try expectPowfBits(0x0000_0000, 0xc080_0000, 0x7f80_0000);

    try expectPowdBits(0x8000_0000_0000_0000, 0x4008_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectPowdBits(0x8000_0000_0000_0000, 0xc008_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectPowdBits(0x8000_0000_0000_0000, 0x4010_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectPowdBits(0x0000_0000_0000_0000, 0xc010_0000_0000_0000, 0x7ff0_0000_0000_0000);
}

test "pow applies C23 infinity and domain rules" {
    try expectPowhBits(0xbc00, 0x7c00, 0x3c00);
    try expectPowhBits(0x3800, 0x7c00, 0x0000);
    try expectPowhBits(0x4000, 0xfc00, 0x0000);
    try expectPowhBits(0xc000, 0x3800, 0x7e00);

    try expectPowfBits(0xbf80_0000, 0x7f80_0000, 0x3f80_0000);
    try expectPowfBits(0x3f00_0000, 0x7f80_0000, 0x0000_0000);
    try expectPowfBits(0x4000_0000, 0xff80_0000, 0x0000_0000);
    try expectPowfBits(0xc000_0000, 0x3f00_0000, 0x7fc0_0000);

    try expectPowdBits(0xbff0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectPowdBits(0x3fe0_0000_0000_0000, 0x7ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectPowdBits(0x4000_0000_0000_0000, 0xfff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectPowdBits(0xc000_0000_0000_0000, 0x3fe0_0000_0000_0000, 0x7ff8_0000_0000_0000);
}

test "pow applies first-NaN payload and C23 identity precedence" {
    try expectPowhBits(0xfe55, 0x4000, 0xfe55);
    try expectPowhBits(0x4000, 0x7c55, 0x7e55);
    try expectPowhBits(0x7c55, 0x7e66, 0x7e55);
    try expectPowhBits(0x7e55, 0x0000, 0x3c00);
    try expectPowhBits(0x7e55, 0x8000, 0x3c00);
    try expectPowhBits(0x7c55, 0x0000, 0x7e55);
    try expectPowhBits(0x7c55, 0x8000, 0x7e55);
    try expectPowhBits(0x3c00, 0x7e55, 0x3c00);
    try expectPowhBits(0x3c00, 0x7c55, 0x7e55);

    try expectPowfBits(0xffc1_2345, 0x4000_0000, 0xffc1_2345);
    try expectPowfBits(0x4000_0000, 0x7f81_2345, 0x7fc1_2345);
    try expectPowfBits(0x7f81_2345, 0x7fc2_3456, 0x7fc1_2345);
    try expectPowfBits(0x7fc1_2345, 0x0000_0000, 0x3f80_0000);
    try expectPowfBits(0x7fc1_2345, 0x8000_0000, 0x3f80_0000);
    try expectPowfBits(0x7f81_2345, 0x0000_0000, 0x7fc1_2345);
    try expectPowfBits(0x7f81_2345, 0x8000_0000, 0x7fc1_2345);
    try expectPowfBits(0x3f80_0000, 0x7fc1_2345, 0x3f80_0000);
    try expectPowfBits(0x3f80_0000, 0x7f81_2345, 0x7fc1_2345);

    try expectPowdBits(0xfff8_1234_5678_9abc, 0x4000_0000_0000_0000, 0xfff8_1234_5678_9abc);
    try expectPowdBits(0x4000_0000_0000_0000, 0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectPowdBits(0x7ff0_1234_5678_9abc, 0x7ff8_2345_6789_abcd, 0x7ff8_1234_5678_9abc);
    try expectPowdBits(0x7ff8_1234_5678_9abc, 0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectPowdBits(0x7ff8_1234_5678_9abc, 0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectPowdBits(0x7ff0_1234_5678_9abc, 0x0000_0000_0000_0000, 0x7ff8_1234_5678_9abc);
    try expectPowdBits(0x7ff0_1234_5678_9abc, 0x8000_0000_0000_0000, 0x7ff8_1234_5678_9abc);
    try expectPowdBits(0x3ff0_0000_0000_0000, 0x7ff8_1234_5678_9abc, 0x3ff0_0000_0000_0000);
    try expectPowdBits(0x3ff0_0000_0000_0000, 0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
}

test "pow rounds proven tiny-exponent domains to one" {
    const tiny_positive_f: u32 = 0x2f00_0000;
    const tiny_negative_f: u32 = 0xaf00_0000;
    try expectPowfBits(0x0000_0001, tiny_positive_f, 0x3f80_0000);
    try expectPowfBits(0x7f7f_ffff, tiny_negative_f, 0x3f80_0000);
    try expectPowfBits(0x3f80_0001, tiny_positive_f, 0x3f80_0000);
    try expectPowfBits(0x3f7f_ffff, tiny_negative_f, 0x3f80_0000);

    const tiny_positive_d: u64 = 0x3be0_0000_0000_0000;
    const tiny_negative_d: u64 = 0xbbe0_0000_0000_0000;
    try expectPowdBits(0x0000_0000_0000_0001, tiny_positive_d, 0x3ff0_0000_0000_0000);
    try expectPowdBits(0x7fef_ffff_ffff_ffff, tiny_negative_d, 0x3ff0_0000_0000_0000);
    try expectPowdBits(0x3ff0_0000_0000_0001, tiny_positive_d, 0x3ff0_0000_0000_0000);
    try expectPowdBits(0x3fef_ffff_ffff_ffff, tiny_negative_d, 0x3ff0_0000_0000_0000);
}

test "pow saturates proven huge-exponent domains around one" {
    const huge_positive_f: u32 = 0x5f00_0000;
    const huge_negative_f: u32 = 0xdf00_0000;
    try expectPowfBits(0x3f7f_ffff, huge_positive_f, 0x0000_0000);
    try expectPowfBits(0x3f7f_ffff, huge_negative_f, 0x7f80_0000);
    try expectPowfBits(0x3f80_0001, huge_positive_f, 0x7f80_0000);
    try expectPowfBits(0x3f80_0001, huge_negative_f, 0x0000_0000);
    try expectPowfBits(0xbf80_0001, huge_positive_f, 0x7f80_0000);
    try expectPowfBits(0xbf80_0001, huge_negative_f, 0x0000_0000);

    const huge_positive_d: u64 = 0x47f0_0000_0000_0000;
    const huge_negative_d: u64 = 0xc7f0_0000_0000_0000;
    try expectPowdBits(0x3fef_ffff_ffff_ffff, huge_positive_d, 0x0000_0000_0000_0000);
    try expectPowdBits(0x3fef_ffff_ffff_ffff, huge_negative_d, 0x7ff0_0000_0000_0000);
    try expectPowdBits(0x3ff0_0000_0000_0001, huge_positive_d, 0x7ff0_0000_0000_0000);
    try expectPowdBits(0x3ff0_0000_0000_0001, huge_negative_d, 0x0000_0000_0000_0000);
    try expectPowdBits(0xbff0_0000_0000_0001, huge_positive_d, 0x7ff0_0000_0000_0000);
    try expectPowdBits(0xbff0_0000_0000_0001, huge_negative_d, 0x0000_0000_0000_0000);
}

test "pow keeps domain and NaN precedence before exponent shortcuts" {
    try expectPowfBits(0xc000_0000, 0x2f00_0000, 0x7fc0_0000);
    try expectPowfBits(0x7f81_2345, 0x2f00_0000, 0x7fc1_2345);
    try expectPowfBits(0x8000_0000, 0x5f00_0000, 0x0000_0000);

    try expectPowdBits(0xc000_0000_0000_0000, 0x3be0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectPowdBits(0x7ff0_1234_5678_9abc, 0x3be0_0000_0000_0000, 0x7ff8_1234_5678_9abc);
    try expectPowdBits(0x8000_0000_0000_0000, 0x47f0_0000_0000_0000, 0x0000_0000_0000_0000);
}
