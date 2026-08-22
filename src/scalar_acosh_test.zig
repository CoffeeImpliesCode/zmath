const std = @import("std");
const scalar = @import("scalar.zig");

fn expectAcoshhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.acoshh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAcoshfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.acoshf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectAcoshdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.acoshd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "acosh exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.acoshh) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.acoshf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.acoshd) == fn (f64) f64);
    }
}

test "acosh handles one representative finite values and infinity" {
    try expectAcoshhBits(0x3c00, 0x0000);
    try expectAcoshhBits(0x3e00, 0x3bb3);
    try expectAcoshhBits(0x4000, 0x3d45);
    try expectAcoshhBits(0x7bff, 0x49e4);
    try expectAcoshhBits(0x7c00, 0x7c00);

    try expectAcoshfBits(0x3f80_0000, 0x0000_0000);
    try expectAcoshfBits(0x3fc0_0000, 0x3f76_6165);
    try expectAcoshfBits(0x4000_0000, 0x3fa8_9214);
    try expectAcoshfBits(0x4120_0000, 0x403f_90f7);
    try expectAcoshfBits(0x7f80_0000, 0x7f80_0000);

    try expectAcoshdBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectAcoshdBits(0x3ff8_0000_0000_0000, 0x3fee_cc2c_aec5_160a);
    try expectAcoshdBits(0x4000_0000_0000_0000, 0x3ff5_1242_7198_0435);
    try expectAcoshdBits(0x4024_0000_0000_0000, 0x4007_f21e_d1ce_05d4);
    try expectAcoshdBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
}

test "acosh returns stable NaNs for inputs below one" {
    try expectAcoshhBits(0x0000, 0x7e00);
    try expectAcoshhBits(0x3bff, 0x7e00);
    try expectAcoshhBits(0x8000, 0x7e00);
    try expectAcoshhBits(0xbc00, 0x7e00);
    try expectAcoshhBits(0xfc00, 0x7e00);

    try expectAcoshfBits(0x0000_0000, 0x7fc0_0000);
    try expectAcoshfBits(0x3f7f_ffff, 0x7fc0_0000);
    try expectAcoshfBits(0x8000_0000, 0x7fc0_0000);
    try expectAcoshfBits(0xbf80_0000, 0x7fc0_0000);
    try expectAcoshfBits(0xff80_0000, 0x7fc0_0000);

    try expectAcoshdBits(0x0000_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectAcoshdBits(0x3fef_ffff_ffff_ffff, 0x7ff8_0000_0000_0000);
    try expectAcoshdBits(0x8000_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectAcoshdBits(0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectAcoshdBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
}

test "acosh quiets NaNs while preserving sign and payload" {
    try expectAcoshhBits(0x7c55, 0x7e55);
    try expectAcoshhBits(0xfc55, 0xfe55);
    try expectAcoshfBits(0x7f81_2345, 0x7fc1_2345);
    try expectAcoshfBits(0xff81_2345, 0xffc1_2345);
    try expectAcoshdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectAcoshdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "acosh covers pinned CORE-MATH binary64 hard cases" {
    try expectAcoshdBits(0x3ff5_bff0_41b2_60fe, 0x3fea_6031_cd5f_93ba);
    try expectAcoshdBits(0x3ff9_efdc_a62b_700a, 0x3ff1_04b6_48f1_13a1);
    try expectAcoshdBits(0x3ffa_5bf3_acfd_e4b2, 0x3ff1_5857_20f3_5cd9);
    try expectAcoshdBits(0x4064_5ea1_60dd_c71f, 0x4017_2581_1dcf_6782);
    try expectAcoshdBits(0x4092_a686_e4b5_67ce, 0x401f_1c92_8e7f_1e65);
    try expectAcoshdBits(0x40ec_b62e_ec26_bd78, 0x4027_59a2_ad4c_4d56);
    try expectAcoshdBits(0x40f3_bf80_0964_8dc0, 0x4027_fce9_5ea5_c653);
}

test "acosh covers CORE-MATH range transitions" {
    try expectAcoshfBits(0x3f80_0001, 0x3a00_0000);
    try expectAcoshfBits(0x3f80_2000, 0x3d35_012e);
    try expectAcoshfBits(0x42df_8000, 0x40ad_19bc);
    try expectAcoshfBits(0x4438_8000, 0x40e9_81c4);
    try expectAcoshfBits(0x5015_02f9, 0x41bd_c082);
    try expectAcoshfBits(0x7f7f_c99e, 0x42b2_d48f);

    try expectAcoshdBits(0x3ff0_0000_0000_0001, 0x3e56_a09e_667f_3bcc);
    try expectAcoshdBits(0x3ff1_0000_0000_0000, 0x3fd6_82dd_8ac1_096d);
    try expectAcoshdBits(0x3ff1_e83e_425a_ee63, 0x3fde_f248_83d6_bb59);
    try expectAcoshdBits(0x405b_f000_0000_0000, 0x4015_a337_7f67_b77e);
    try expectAcoshdBits(0x4087_1000_0000_0000, 0x401d_3038_810e_918e);
    try expectAcoshdBits(0x40e0_1000_0000_0000, 0x4026_3041_ffa2_68fa);
    try expectAcoshdBits(0x4190_0000_0000_0000, 0x4032_b708_8723_20e2);
    try expectAcoshdBits(0x4330_0000_0000_0000, 0x4042_5e4f_7b27_37fa);
    try expectAcoshdBits(0x7fe0_0000_0000_0000, 0x4086_2e42_fefa_39ef);
}
