const std = @import("std");
const scalar = @import("scalar.zig");

fn expectAcoshBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.acosh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAcosfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.acosf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAcosdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.acosd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "acos exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.acosh) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.acosf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.acosd) == fn (f64) f64);
    }
}

test "acos maps signed zero and unit endpoints exactly" {
    try expectAcoshBits(0x0000, 0x3e48);
    try expectAcoshBits(0x8000, 0x3e48);
    try expectAcoshBits(0x3c00, 0x0000);
    try expectAcoshBits(0xbc00, 0x4248);

    try expectAcosfBits(0x0000_0000, 0x3fc9_0fdb);
    try expectAcosfBits(0x8000_0000, 0x3fc9_0fdb);
    try expectAcosfBits(0x3f80_0000, 0x0000_0000);
    try expectAcosfBits(0xbf80_0000, 0x4049_0fdb);

    try expectAcosdBits(0x0000_0000_0000_0000, 0x3ff9_21fb_5444_2d18);
    try expectAcosdBits(0x8000_0000_0000_0000, 0x3ff9_21fb_5444_2d18);
    try expectAcosdBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectAcosdBits(0xbff0_0000_0000_0000, 0x4009_21fb_5444_2d18);
}

test "acos handles representative range-reduction regions" {
    try expectAcoshBits(0x3400, 0x3d46);
    try expectAcoshBits(0xb400, 0x3f4b);
    try expectAcoshBits(0x3800, 0x3c30);
    try expectAcoshBits(0xb800, 0x4030);

    try expectAcosfBits(0x3e80_0000, 0x3fa8_b807);
    try expectAcosfBits(0xbe80_0000, 0x3fe9_67ae);
    try expectAcosfBits(0x3f00_0000, 0x3f86_0a92);
    try expectAcosfBits(0xbf00_0000, 0x4006_0a92);

    try expectAcosdBits(0x3fd0_0000_0000_0000, 0x3ff5_1700_e0c1_4b25);
    try expectAcosdBits(0xbfd0_0000_0000_0000, 0x3ffd_2cf5_c7c7_0f0c);
    try expectAcosdBits(0x3fe0_0000_0000_0000, 0x3ff0_c152_382d_7366);
    try expectAcosdBits(0xbfe0_0000_0000_0000, 0x4000_c152_382d_7366);
    try expectAcosdBits(0x0000_0000_0000_0001, 0x3ff9_21fb_5444_2d18);
    try expectAcosdBits(0x8000_0000_0000_0001, 0x3ff9_21fb_5444_2d18);
}

test "acos returns stable NaNs for invalid inputs" {
    try expectAcoshBits(0x3c01, 0x7e00);
    try expectAcoshBits(0xbc01, 0x7e00);
    try expectAcoshBits(0x7c00, 0x7e00);
    try expectAcoshBits(0xfc00, 0x7e00);
    try expectAcoshBits(0x7c55, 0x7e55);
    try expectAcoshBits(0xfc55, 0xfe55);

    try expectAcosfBits(0x3f80_0001, 0x7fc0_0000);
    try expectAcosfBits(0xbf80_0001, 0x7fc0_0000);
    try expectAcosfBits(0x7f80_0000, 0x7fc0_0000);
    try expectAcosfBits(0xff80_0000, 0x7fc0_0000);
    try expectAcosfBits(0x7f81_2345, 0x7fc1_2345);
    try expectAcosfBits(0xff81_2345, 0xffc1_2345);

    try expectAcosdBits(0x3ff0_0000_0000_0001, 0x7ff8_0000_0000_0000);
    try expectAcosdBits(0xbff0_0000_0000_0001, 0x7ff8_0000_0000_0000);
    try expectAcosdBits(0x7ff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectAcosdBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectAcosdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectAcosdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "acos covers pinned CORE-MATH hard cases" {
    try expectAcosfBits(0x3288_85a3, 0x3fc9_0fdb);
    try expectAcosfBits(0x3982_6222, 0x3fc9_07b5);

    try expectAcosdBits(0x3fef_ffff_ffff_fdc0, 0x3e98_0000_0000_0024);
    try expectAcosdBits(0x3fb5_3ea6_c725_5e88, 0x3ff7_cdac_b6bb_e707);
    try expectAcosdBits(0x3f4f_d737_be91_4578, 0x3ff9_1e00_6d41_d8d9);
    try expectAcosdBits(0x3fef_ffff_ffff_ff70, 0x3e88_0000_0000_0009);
    try expectAcosdBits(0x3fa3_90e6_939c_d1a6, 0x3ff8_856a_5d32_96a4);
}
