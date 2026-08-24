const std = @import("std");
const scalar = @import("scalar.zig");

// Target-format bit comparisons keep exceptional behavior independent of the host.
fn expectTgammaHBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.tgammah(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectTgammaFBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.tgammaf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectTgammaDBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.tgammad(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "tgamma preserves quiet NaN payloads and signs and quiets signaling NaNs" {
    try expectTgammaHBits(0x7e55, 0x7e55);
    try expectTgammaHBits(0xfe55, 0xfe55);
    try expectTgammaHBits(0x7d55, 0x7f55);
    try expectTgammaHBits(0xfd55, 0xff55);

    try expectTgammaFBits(0x7fc1_2345, 0x7fc1_2345);
    try expectTgammaFBits(0xffc1_2345, 0xffc1_2345);
    try expectTgammaFBits(0x7fa1_2345, 0x7fe1_2345);
    try expectTgammaFBits(0xffa1_2345, 0xffe1_2345);

    try expectTgammaDBits(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectTgammaDBits(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectTgammaDBits(0x7ff1_2345_6789_abcd, 0x7ff9_2345_6789_abcd);
    try expectTgammaDBits(0xfff1_2345_6789_abcd, 0xfff9_2345_6789_abcd);
}

test "tgamma handles infinities with target-format bits" {
    try expectTgammaHBits(0x7c00, 0x7c00);
    try expectTgammaHBits(0xfc00, 0x7e00);

    try expectTgammaFBits(0x7f80_0000, 0x7f80_0000);
    try expectTgammaFBits(0xff80_0000, 0x7fc0_0000);

    try expectTgammaDBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectTgammaDBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
}

test "tgamma maps signed zeros to signed infinities" {
    try expectTgammaHBits(0x0000, 0x7c00);
    try expectTgammaHBits(0x8000, 0xfc00);

    try expectTgammaFBits(0x0000_0000, 0x7f80_0000);
    try expectTgammaFBits(0x8000_0000, 0xff80_0000);

    try expectTgammaDBits(0x0000_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectTgammaDBits(0x8000_0000_0000_0000, 0xfff0_0000_0000_0000);
}

test "tgamma maps negative integer poles to canonical quiet NaNs" {
    try expectTgammaHBits(0xbc00, 0x7e00);
    try expectTgammaHBits(0xc000, 0x7e00);

    try expectTgammaFBits(0xbf80_0000, 0x7fc0_0000);
    try expectTgammaFBits(0xc000_0000, 0x7fc0_0000);

    try expectTgammaDBits(0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectTgammaDBits(0xc000_0000_0000_0000, 0x7ff8_0000_0000_0000);
}
