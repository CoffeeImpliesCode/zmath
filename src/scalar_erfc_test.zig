const std = @import("std");
const root = @import("root.zig");
const scalar = @import("scalar.zig");

fn expectErfchBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.erfch(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectErfcfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.erfcf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectErfcdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.erfcd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "Erfc family exposes lower-case suffixes" {
    comptime {
        std.debug.assert(@TypeOf(scalar.erfch(@as(f16, 0.0))) == f16);
        std.debug.assert(@TypeOf(scalar.erfcf(@as(f32, 0.0))) == f32);
        std.debug.assert(@TypeOf(scalar.erfcd(@as(f64, 0.0))) == f64);
    }
}

test "Erfc family handles zeros and infinities" {
    try expectErfchBits(0x0000, 0x3c00);
    try expectErfchBits(0x8000, 0x3c00);
    try expectErfchBits(0x7c00, 0x0000);
    try expectErfchBits(0xfc00, 0x4000);

    try expectErfcfBits(0x0000_0000, 0x3f80_0000);
    try expectErfcfBits(0x8000_0000, 0x3f80_0000);
    try expectErfcfBits(0x7f80_0000, 0x0000_0000);
    try expectErfcfBits(0xff80_0000, 0x4000_0000);

    try expectErfcdBits(0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectErfcdBits(0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectErfcdBits(0x7ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectErfcdBits(0xfff0_0000_0000_0000, 0x4000_0000_0000_0000);
}

test "Erfc family rounds representative finite values" {
    try expectErfchBits(0x3800, 0x37ac);
    try expectErfchBits(0x3c00, 0x3109);
    try expectErfchBits(0xbc00, 0x3f5f);

    try expectErfcfBits(0x3f00_0000, 0x3ef5_810a);
    try expectErfcfBits(0x3f80_0000, 0x3e21_130b);
    try expectErfcfBits(0xbf80_0000, 0x3feb_dd9f);
    try expectErfcfBits(0x4120_0000, 0x0000_0001);

    try expectErfcdBits(0x3fe0_0000_0000_0000, 0x3fde_b021_47ce_245c);
    try expectErfcdBits(0x3ff0_0000_0000_0000, 0x3fc4_2261_62fb_ddd5);
    try expectErfcdBits(0xbff0_0000_0000_0000, 0x3ffd_7bb3_d3a0_8445);
    try expectErfcdBits(0x4000_0000_0000_0000, 0x3f73_28f5_ec35_0e67);
    try expectErfcdBits(0x4010_0000_0000_0000, 0x3e50_8ddd_13bd_35e7);
    try expectErfcdBits(0x4020_0000_0000_0000, 0x39ec_74fc_4121_7dfb);
    try expectErfcdBits(0x403b_0000_0000_0000, 0x0000_0000_0001_9e0f);
}

test "Erfcd preserves exact tiny and saturation cutovers" {
    try expectErfcdBits(0x3c8c_5bf8_91b4_ef6a, 0x3ff0_0000_0000_0000);
    try expectErfcdBits(0xbc9c_5bf8_91b4_ef6a, 0x3ff0_0000_0000_0000);
    try expectErfcdBits(0x403b_39dc_41e4_8bfd, 0x0000_0000_0000_0000);
    try expectErfcdBits(0xc017_744f_8f74_e94b, 0x4000_0000_0000_0000);
}

test "Erfc wrappers quiet NaNs and preserve payload and sign" {
    try expectErfchBits(0x7c01, 0x7e01);
    try expectErfchBits(0xfe55, 0xfe55);
    try expectErfchBits(0x7d55, 0x7f55);
    try expectErfcfBits(0x7f80_0001, 0x7fc0_0001);
    try expectErfcfBits(0xffc1_2345, 0xffc1_2345);
    try expectErfcfBits(0x7fa1_2345, 0x7fe1_2345);
    try expectErfcdBits(0x7ff0_0000_0000_0001, 0x7ff8_0000_0000_0001);
    try expectErfcdBits(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectErfcdBits(0x7ff1_2345_6789_abcd, 0x7ff9_2345_6789_abcd);
}

test "root Erfc dispatch selects every scalar format" {
    try std.testing.expectEqual(@as(u16, 0x3109), @as(u16, @bitCast(root.erfc(@as(f16, 1.0)))));
    try std.testing.expectEqual(@as(u32, 0x3e21_130b), @as(u32, @bitCast(root.erfc(@as(f32, 1.0)))));
    try std.testing.expectEqual(
        @as(u64, 0x3fc4_2261_62fb_ddd5),
        @as(u64, @bitCast(root.erfc(@as(f64, 1.0)))),
    );
}
