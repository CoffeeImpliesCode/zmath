const std = @import("std");
const root = @import("root.zig");
const scalar = @import("scalar.zig");

fn expectTgammahBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.tgammah(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectTgammafBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.tgammaf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectTgammadBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.tgammad(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "Tgamma family exposes lower-case suffixes" {
    comptime {
        std.debug.assert(@TypeOf(scalar.tgammah(@as(f16, 1.0))) == f16);
        std.debug.assert(@TypeOf(scalar.tgammaf(@as(f32, 1.0))) == f32);
        std.debug.assert(@TypeOf(scalar.tgammad(@as(f64, 1.0))) == f64);
    }
}

test "Tgamma family handles zeros infinities and poles" {
    try expectTgammahBits(0x0000, 0x7c00);
    try expectTgammahBits(0x8000, 0xfc00);
    try expectTgammahBits(0x7c00, 0x7c00);
    try expectTgammahBits(0xfc00, 0x7e00);
    try expectTgammahBits(0xbc00, 0x7e00);
    try expectTgammahBits(0xc000, 0x7e00);

    try expectTgammafBits(0x0000_0000, 0x7f80_0000);
    try expectTgammafBits(0x8000_0000, 0xff80_0000);
    try expectTgammafBits(0x7f80_0000, 0x7f80_0000);
    try expectTgammafBits(0xff80_0000, 0x7fc0_0000);
    try expectTgammafBits(0xbf80_0000, 0x7fc0_0000);
    try expectTgammafBits(0xc000_0000, 0x7fc0_0000);

    try expectTgammadBits(0x0000_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectTgammadBits(0x8000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectTgammadBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectTgammadBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectTgammadBits(0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectTgammadBits(0xc000_0000_0000_0000, 0x7ff8_0000_0000_0000);
}

test "Tgamma wrappers quiet NaNs and preserve payload and sign" {
    try expectTgammahBits(0x7c01, 0x7e01);
    try expectTgammahBits(0xfe55, 0xfe55);
    try expectTgammahBits(0x7d55, 0x7f55);
    try expectTgammafBits(0x7f80_0001, 0x7fc0_0001);
    try expectTgammafBits(0xffc1_2345, 0xffc1_2345);
    try expectTgammafBits(0x7fa1_2345, 0x7fe1_2345);
    try expectTgammadBits(0x7ff0_0000_0000_0001, 0x7ff8_0000_0000_0001);
    try expectTgammadBits(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectTgammadBits(0x7ff1_2345_6789_abcd, 0x7ff9_2345_6789_abcd);
}

test "Tgamma family rounds exact and representative finite values" {
    try expectTgammahBits(0x3c00, 0x3c00);
    try expectTgammahBits(0x4000, 0x3c00);
    try expectTgammahBits(0x4200, 0x4000);
    try expectTgammahBits(0x4400, 0x4600);
    try expectTgammahBits(0x4500, 0x4e00);
    try expectTgammahBits(0x3800, 0x3f17);
    try expectTgammahBits(0xb800, 0xc317);
    try expectTgammahBits(0x3e00, 0x3b17);
    try expectTgammahBits(0xbe00, 0x40ba);
    try expectTgammahBits(0xc100, 0xbb90);

    try expectTgammafBits(0x3f80_0000, 0x3f80_0000);
    try expectTgammafBits(0x4000_0000, 0x3f80_0000);
    try expectTgammafBits(0x4040_0000, 0x4000_0000);
    try expectTgammafBits(0x40a0_0000, 0x41c0_0000);
    try expectTgammafBits(0x3f00_0000, 0x3fe2_dfc5);
    try expectTgammafBits(0xbf00_0000, 0xc062_dfc5);
    try expectTgammafBits(0x3fc0_0000, 0x3f62_dfc5);
    try expectTgammafBits(0xbfc0_0000, 0x4017_3fd8);
    try expectTgammafBits(0xc020_0000, 0xbf71_ffc1);

    try expectTgammadBits(0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectTgammadBits(0x4000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectTgammadBits(0x4008_0000_0000_0000, 0x4000_0000_0000_0000);
    try expectTgammadBits(0x4014_0000_0000_0000, 0x4038_0000_0000_0000);
    try expectTgammadBits(0x3fe0_0000_0000_0000, 0x3ffc_5bf8_91b4_ef6b);
    try expectTgammadBits(0xbfe0_0000_0000_0000, 0xc00c_5bf8_91b4_ef6b);
    try expectTgammadBits(0x3ff8_0000_0000_0000, 0x3fec_5bf8_91b4_ef6b);
    try expectTgammadBits(0xbff8_0000_0000_0000, 0x4002_e7fb_0bcd_f4f2);
    try expectTgammadBits(0xc004_0000_0000_0000, 0xbfee_3ff8_12e3_2183);
}

test "Tgamma family covers overflow underflow and hard cases" {
    try expectTgammahBits(0x0001, 0x7c00);
    try expectTgammahBits(0x8001, 0xfc00);
    try expectTgammahBits(0x4900, 0x7c00);
    try expectTgammahBits(0xcc60, 0x0000);
    try expectTgammahBits(0xcca0, 0x8000);

    try expectTgammafBits(0x0000_0001, 0x7f80_0000);
    try expectTgammafBits(0x8000_0001, 0xff80_0000);
    try expectTgammafBits(0x4210_0000, 0x7f80_0000);
    try expectTgammafBits(0xc22a_0000, 0x8000_0000);
    try expectTgammafBits(0xc22e_0000, 0x0000_0000);
    try expectTgammafBits(0x27de_86a9, 0x5713_4133);
    try expectTgammafBits(0xb63b_efb3, 0xc8ae_5b72);

    try expectTgammadBits(0x0000_0000_0000_0001, 0x7ff0_0000_0000_0000);
    try expectTgammadBits(0x8000_0000_0000_0001, 0xfff0_0000_0000_0000);
    try expectTgammadBits(0xc067_1000_0000_0000, 0x8000_0000_0000_0000);
    try expectTgammadBits(0xc067_3000_0000_0000, 0x0000_0000_0000_0000);
    try expectTgammadBits(0x3fdc_05aa_42cb_27fe, 0x4000_2f80_f15c_9486);
    try expectTgammadBits(0xc013_fc07_c800_57fd, 0xc001_4fd6_b28f_b843);
}

test "root Tgamma dispatch selects every scalar format" {
    try std.testing.expectEqual(
        @as(u16, 0x3f17),
        @as(u16, @bitCast(root.tgamma(@as(f16, 0.5)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3fe2_dfc5),
        @as(u32, @bitCast(root.tgamma(@as(f32, 0.5)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3ffc_5bf8_91b4_ef6b),
        @as(u64, @bitCast(root.tgamma(@as(f64, 0.5)))),
    );
}
