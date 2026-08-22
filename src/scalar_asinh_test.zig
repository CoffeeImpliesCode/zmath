const std = @import("std");
const scalar = @import("scalar.zig");

fn expectAsinhhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.asinhh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAsinhfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.asinhf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAsinhdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.asinhd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "asinh exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.asinhh) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.asinhf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.asinhd) == fn (f64) f64);
    }
}

test "asinh preserves signed zero and rounds representative values" {
    try expectAsinhhBits(0x0000, 0x0000);
    try expectAsinhhBits(0x8000, 0x8000);
    try expectAsinhhBits(0x3800, 0x37b3);
    try expectAsinhhBits(0x3c00, 0x3b0d);
    try expectAsinhhBits(0x4000, 0x3dc6);

    try expectAsinhfBits(0x0000_0000, 0x0000_0000);
    try expectAsinhfBits(0x8000_0000, 0x8000_0000);
    try expectAsinhfBits(0x3f00_0000, 0x3ef6_6165);
    try expectAsinhfBits(0x3f80_0000, 0x3f61_a1b3);
    try expectAsinhfBits(0x4000_0000, 0x3fb8_c90c);

    try expectAsinhdBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectAsinhdBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectAsinhdBits(0x3fe0_0000_0000_0000, 0x3fde_cc2c_aec5_160a);
    try expectAsinhdBits(0x3ff0_0000_0000_0000, 0x3fec_3436_6179_d427);
    try expectAsinhdBits(0x4000_0000_0000_0000, 0x3ff7_1921_8313_d087);
}

test "asinh preserves infinities and quiets NaNs with stable payloads" {
    try expectAsinhhBits(0x7c00, 0x7c00);
    try expectAsinhhBits(0xfc00, 0xfc00);
    try expectAsinhhBits(0x7c55, 0x7e55);
    try expectAsinhhBits(0xfc55, 0xfe55);

    try expectAsinhfBits(0x7f80_0000, 0x7f80_0000);
    try expectAsinhfBits(0xff80_0000, 0xff80_0000);
    try expectAsinhfBits(0x7f81_2345, 0x7fc1_2345);
    try expectAsinhfBits(0xff81_2345, 0xffc1_2345);

    try expectAsinhdBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectAsinhdBits(0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectAsinhdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectAsinhdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "asinh covers pinned CORE-MATH correction cases" {
    try expectAsinhhBits(0x4cb0, 0x4340);
    try expectAsinhhBits(0x5b60, 0x4628);
    try expectAsinhhBits(0x5f88, 0x46df);
    try expectAsinhhBits(0x62c4, 0x4775);
    try expectAsinhdBits(0x3e59_ffff_ffff_ffff, 0x3e59_ffff_ffff_fffe);
    try expectAsinhdBits(0x3e5a_0000_0000_0000, 0x3e59_ffff_ffff_ffff);
    try expectAsinhdBits(0x3e5a_0000_0000_0001, 0x3e5a_0000_0000_0000);

    try expectAsinhdBits(0x3f2f_ffff_ffff_ffff, 0x3f2f_ffff_faaa_aaac);
    try expectAsinhdBits(0x3f30_0000_0000_0000, 0x3f2f_ffff_faaa_aaad);
    try expectAsinhdBits(0x3f30_0000_0000_0001, 0x3f2f_ffff_faaa_aaaf);

    try expectAsinhdBits(0x3f92_ffff_ffff_ffff, 0x3f92_ffb8_902a_9ab5);
    try expectAsinhdBits(0x3f93_0000_0000_0000, 0x3f92_ffb8_902a_9ab6);
    try expectAsinhdBits(0x3f93_0000_0000_0001, 0x3f92_ffb8_902a_9ab7);

    try expectAsinhfBits(0x3e81_5668, 0x3e80_0001);

    try expectAsinhdBits(0x3f32_bf10_2b93_b850, 0x3f32_bf10_2749_baf8);
    try expectAsinhdBits(0x3f93_aad1_e7b1_88fc, 0x3f93_aa82_acdc_a360);

    try expectAsinhdBits(0x3fd0_0f94_7645_0863, 0x3fcf_cb35_067f_343c);
    try expectAsinhdBits(0x4170_fbc6_c02b_1c90, 0x4031_6369_cd53_bb69);
}

test "asinhh is odd for every finite binary16 input" {
    var magnitude: u16 = 0;
    while (magnitude < 0x7c00) : (magnitude += 1) {
        const positive_bits: u16 = @bitCast(scalar.asinhh(@bitCast(magnitude)));
        const negative_bits: u16 = @bitCast(scalar.asinhh(@bitCast(magnitude | 0x8000)));
        try std.testing.expectEqual(positive_bits | 0x8000, negative_bits);
    }
}
