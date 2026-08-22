const std = @import("std");
const scalar = @import("scalar.zig");

fn expectSinhhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u16, @bitCast(scalar.sinhh(input))));
}

fn expectSinhfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u32, @bitCast(scalar.sinhf(input))));
}

fn expectSinhdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    try std.testing.expectEqual(expected_bits, @as(u64, @bitCast(scalar.sinhd(input))));
}

test "hyperbolic sinh exports use nearest-even scalar signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.sinhh) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.sinhf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.sinhd) == fn (f64) f64);
    }
}

test "hyperbolic sinh preserves signed zeros and tiny values" {
    try expectSinhhBits(0x0000, 0x0000);
    try expectSinhhBits(0x8000, 0x8000);
    try expectSinhhBits(0x0001, 0x0001);
    try expectSinhhBits(0x8001, 0x8001);

    try expectSinhfBits(0x0000_0000, 0x0000_0000);
    try expectSinhfBits(0x8000_0000, 0x8000_0000);
    try expectSinhfBits(0x0000_0001, 0x0000_0001);
    try expectSinhfBits(0x8000_0001, 0x8000_0001);

    try expectSinhdBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectSinhdBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectSinhdBits(0x0000_0000_0000_0001, 0x0000_0000_0000_0001);
    try expectSinhdBits(0x8000_0000_0000_0001, 0x8000_0000_0000_0001);
}

test "hyperbolic sinh keeps subnormal identity under nearest-even" {
    try expectSinhhBits(0x03ff, 0x03ff);
    try expectSinhhBits(0x83ff, 0x83ff);
    try expectSinhfBits(0x007f_ffff, 0x007f_ffff);
    try expectSinhfBits(0x807f_ffff, 0x807f_ffff);
    try expectSinhdBits(0x000f_5add_acfe_2191, 0x000f_5add_acfe_2191);
    try expectSinhdBits(0x800f_5add_acfe_2191, 0x800f_5add_acfe_2191);
}

test "hyperbolic sinh preserves infinities and quiets NaNs" {
    try expectSinhhBits(0x7c00, 0x7c00);
    try expectSinhhBits(0xfc00, 0xfc00);
    try expectSinhhBits(0x7d55, 0x7f55);
    try expectSinhhBits(0xfd55, 0xff55);

    try expectSinhfBits(0x7f80_0000, 0x7f80_0000);
    try expectSinhfBits(0xff80_0000, 0xff80_0000);
    try expectSinhfBits(0x7f81_2345, 0x7fc1_2345);
    try expectSinhfBits(0xff81_2345, 0xffc1_2345);

    try expectSinhdBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectSinhdBits(0xfff0_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectSinhdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectSinhdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "hyperbolic sinh rounds representative CORE-MATH values" {
    try expectSinhhBits(0x3000, 0x3003);
    try expectSinhhBits(0x3c00, 0x3cb3);
    try expectSinhhBits(0xbc00, 0xbcb3);
    try expectSinhhBits(0x4500, 0x54a3);
    try expectSinhhBits(0x4900, 0x7161);

    try expectSinhfBits(0x3e00_0000, 0x3e00_5566);
    try expectSinhfBits(0x3f80_0000, 0x3f96_6cfe);
    try expectSinhfBits(0xbf80_0000, 0xbf96_6cfe);
    try expectSinhfBits(0x40a0_0000, 0x4294_680b);
    try expectSinhfBits(0x4120_0000, 0x462c_14ee);

    try expectSinhdBits(0x3fc0_0000_0000_0000, 0x3fc0_0aac_cd00_d2f1);
    try expectSinhdBits(0x3ff0_0000_0000_0000, 0x3ff2_cd9f_c44e_b982);
    try expectSinhdBits(0xbff0_0000_0000_0000, 0xbff2_cd9f_c44e_b982);
    try expectSinhdBits(0x4014_0000_0000_0000, 0x4052_8d01_66f0_7374);
    try expectSinhdBits(0x4024_0000_0000_0000, 0x40c5_829d_ced6_9992);
}

test "hyperbolic sinh rounds finite overflow boundaries" {
    try expectSinhhBits(0x49e4, 0x7bfb);
    try expectSinhhBits(0x49e5, 0x7c00);
    try expectSinhhBits(0xc9e4, 0xfbfb);
    try expectSinhhBits(0xc9e5, 0xfc00);

    try expectSinhfBits(0x42b2_d4fc, 0x7f7f_ffec);
    try expectSinhfBits(0x42b2_d4fd, 0x7f80_0000);
    try expectSinhfBits(0xc2b2_d4fc, 0xff7f_ffec);
    try expectSinhfBits(0xc2b2_d4fd, 0xff80_0000);

    try expectSinhdBits(0x4086_33ce_8fb9_f87d, 0x7fef_ffff_ffff_fd3b);
    try expectSinhdBits(0x4086_33ce_8fb9_f87e, 0x7ff0_0000_0000_0000);
    try expectSinhdBits(0xc086_33ce_8fb9_f87d, 0xffef_ffff_ffff_fd3b);
    try expectSinhdBits(0xc086_33ce_8fb9_f87e, 0xfff0_0000_0000_0000);
}

test "hyperbolic sinhd applies database corrections" {
    try expectSinhdBits(0x3fd3_6430_3e1a_d8f6, 0x3fd3_b07e_0c77_9dda);
    try expectSinhdBits(0xbfd3_6430_3e1a_d8f6, 0xbfd3_b07e_0c77_9dda);
    try expectSinhdBits(0x4084_a869_881f_72ac, 0x7b79_ea75_40a3_d1f9);
}

test "hyperbolic sinhd resolves near-zero high-bit corrections" {
    try expectSinhdBits(0x3fc0_0a0a_44cc_45b0, 0x3fc0_14cb_39ad_e240);
    try expectSinhdBits(0xbfc0_0a0a_44cc_45b0, 0xbfc0_14cb_39ad_e240);
    try expectSinhdBits(0x3fc0_dc68_d5e8_f959, 0x3fc0_e8e7_3dc4_fee4);
    try expectSinhdBits(0xbfc0_dc68_d5e8_f959, 0xbfc0_e8e7_3dc4_fee4);
    try expectSinhdBits(0x3f91_9e03_c96f_0997, 0x3f91_9e3c_be7e_f607);
    try expectSinhdBits(0xbf91_9e03_c96f_0997, 0xbf91_9e3c_be7e_f607);
    try expectSinhdBits(0x3f41_bd15_d167_0050, 0x3f41_bd15_dff0_122b);
    try expectSinhdBits(0xbf41_bd15_d167_0050, 0xbf41_bd15_dff0_122b);
}
