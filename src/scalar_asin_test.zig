const std = @import("std");
const scalar = @import("scalar.zig");

fn expectAsinhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.asinh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAsinfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.asinf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAsindBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.asind(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "asin exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.asinh) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.asinf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.asind) == fn (f64) f64);
    }
}

test "asin preserves signed zero and handles unit endpoints" {
    try expectAsinhBits(0x0000, 0x0000);
    try expectAsinhBits(0x8000, 0x8000);
    try expectAsinhBits(0x3c00, 0x3e48);
    try expectAsinhBits(0xbc00, 0xbe48);

    try expectAsinfBits(0x0000_0000, 0x0000_0000);
    try expectAsinfBits(0x8000_0000, 0x8000_0000);
    try expectAsinfBits(0x3f80_0000, 0x3fc9_0fdb);
    try expectAsinfBits(0xbf80_0000, 0xbfc9_0fdb);

    try expectAsindBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectAsindBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectAsindBits(0x3ff0_0000_0000_0000, 0x3ff9_21fb_5444_2d18);
    try expectAsindBits(0xbff0_0000_0000_0000, 0xbff9_21fb_5444_2d18);
}

test "asin handles finite boundaries and range reduction" {
    try expectAsinhBits(0x0001, 0x0001);
    try expectAsinhBits(0x03ff, 0x03ff);
    try expectAsinhBits(0x0400, 0x0400);
    try expectAsinhBits(0x3400, 0x340b);
    try expectAsinhBits(0x3800, 0x3830);
    try expectAsinhBits(0x3bff, 0x3e28);

    try expectAsinfBits(0x0000_0001, 0x0000_0001);
    try expectAsinfBits(0x007f_ffff, 0x007f_ffff);
    try expectAsinfBits(0x0080_0000, 0x0080_0000);
    try expectAsinfBits(0x3e80_0000, 0x3e81_5f4e);
    try expectAsinfBits(0x3f00_0000, 0x3f06_0a92);
    try expectAsinfBits(0x3f7f_ffff, 0x3fc9_048a);

    try expectAsindBits(0x0000_0000_0000_0001, 0x0000_0000_0000_0001);
    try expectAsindBits(0x000f_ffff_ffff_ffff, 0x000f_ffff_ffff_ffff);
    try expectAsindBits(0x0010_0000_0000_0000, 0x0010_0000_0000_0000);
    try expectAsindBits(0x3fd0_0000_0000_0000, 0x3fd0_2be9_ce0b_87cd);
    try expectAsindBits(0x3fe0_0000_0000_0000, 0x3fe0_c152_382d_7366);
    try expectAsindBits(0x3fef_ffff_ffff_ffff, 0x3ff9_21fb_5044_2d18);
}

test "asin returns stable NaNs for invalid inputs" {
    try expectAsinhBits(0x3c01, 0x7e00);
    try expectAsinhBits(0xbc01, 0x7e00);
    try expectAsinhBits(0x7c00, 0x7e00);
    try expectAsinhBits(0xfc00, 0x7e00);
    try expectAsinhBits(0x7c55, 0x7e55);
    try expectAsinhBits(0xfc55, 0xfe55);

    try expectAsinfBits(0x3f80_0001, 0x7fc0_0000);
    try expectAsinfBits(0xbf80_0001, 0x7fc0_0000);
    try expectAsinfBits(0x7f80_0000, 0x7fc0_0000);
    try expectAsinfBits(0xff80_0000, 0x7fc0_0000);
    try expectAsinfBits(0x7f81_2345, 0x7fc1_2345);
    try expectAsinfBits(0xff81_2345, 0xffc1_2345);

    try expectAsindBits(0x3ff0_0000_0000_0001, 0x7ff8_0000_0000_0000);
    try expectAsindBits(0xbff0_0000_0000_0001, 0x7ff8_0000_0000_0000);
    try expectAsindBits(0x7ff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectAsindBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectAsindBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectAsindBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "asin covers pinned CORE-MATH hard cases" {
    try expectAsinhBits(0x2e98, 0x2e9b);
    try expectAsinhBits(0xae98, 0xae9b);
    try expectAsinhBits(0x3592, 0x35b1);
    try expectAsinhBits(0x36cb, 0x3704);
    try expectAsinhBits(0xb6cb, 0xb704);
    try expectAsinhBits(0xb592, 0xb5b1);

    try expectAsinfBits(0x3f2a_b445, 0x3f3a_dc51);
    try expectAsinfBits(0xbf2a_b445, 0xbf3a_dc51);
    try expectAsinfBits(0x3f08_3a1a, 0x3f0f_a5b2);
    try expectAsinfBits(0xbf08_3a1a, 0xbf0f_a5b2);

    try expectAsindBits(0x3e57_1374_4912_3ef6, 0x3e57_1374_4912_3ef7);
    try expectAsindBits(0xbe57_1374_4912_3ef6, 0xbe57_1374_4912_3ef7);
    try expectAsindBits(0x3ef2_bfff_fffc_233b, 0x3ef2_c000_0000_6ddd);
    try expectAsindBits(0xbef2_bfff_fffc_233b, 0xbef2_c000_0000_6ddd);
    try expectAsindBits(0x3fdd_a4e0_e6c7_17a5, 0x3fde_d25c_5eb8_c916);
    try expectAsindBits(0xbfdd_a4e0_e6c7_17a5, 0xbfde_d25c_5eb8_c916);
}

test "asinh is odd for every finite in-domain binary16 input" {
    var magnitude: u16 = 0;
    while (magnitude <= 0x3c00) : (magnitude += 1) {
        const positive_bits: u16 = @bitCast(scalar.asinh(@bitCast(magnitude)));
        const negative_input = magnitude | 0x8000;
        const negative_bits: u16 = @bitCast(scalar.asinh(@bitCast(negative_input)));
        try std.testing.expectEqual(positive_bits | 0x8000, negative_bits);
    }
}
