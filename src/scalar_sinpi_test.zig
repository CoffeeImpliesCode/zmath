const std = @import("std");
const scalar = @import("scalar.zig");

fn expectH(input_bits: u16, expected_bits: u16) !void {
    const actual_bits: u16 = @bitCast(scalar.sinpih(@as(f16, @bitCast(input_bits))));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectF(input_bits: u32, expected_bits: u32) !void {
    const actual_bits: u32 = @bitCast(scalar.sinpif(@as(f32, @bitCast(input_bits))));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectD(input_bits: u64, expected_bits: u64) !void {
    const actual_bits: u64 = @bitCast(scalar.sinpid(@as(f64, @bitCast(input_bits))));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "sinpi exports nearest-even binary16, binary32, and binary64 functions" {
    comptime std.debug.assert(@TypeOf(scalar.sinpih) == fn (f16) f16);
    comptime std.debug.assert(@TypeOf(scalar.sinpif) == fn (f32) f32);
    comptime std.debug.assert(@TypeOf(scalar.sinpid) == fn (f64) f64);
}

test "sinpi preserves signed zero and exact half-turn values" {
    try expectH(0x0000, 0x0000);
    try expectH(0x8000, 0x8000);
    try expectH(0x3c00, 0x0000);
    try expectH(0xbc00, 0x8000);
    try expectH(0x3800, 0x3c00);
    try expectH(0xb800, 0xbc00);

    try expectF(0x0000_0000, 0x0000_0000);
    try expectF(0x8000_0000, 0x8000_0000);
    try expectF(0x3f80_0000, 0x0000_0000);
    try expectF(0xbf80_0000, 0x8000_0000);
    try expectF(0x3f00_0000, 0x3f80_0000);
    try expectF(0xbf00_0000, 0xbf80_0000);

    try expectD(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectD(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectD(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectD(0xbff0_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectD(0x3fe0_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectD(0xbfe0_0000_0000_0000, 0xbff0_0000_0000_0000);
}

test "sinpi maps infinities to canonical NaNs and preserves signed NaN payloads" {
    try expectH(0x7c00, 0x7e00);
    try expectH(0xfc00, 0x7e00);
    try expectH(0x7d01, 0x7f01);
    try expectH(0xfd01, 0xff01);
    try expectH(0x7e55, 0x7e55);
    try expectH(0xfe55, 0xfe55);

    try expectF(0x7f80_0000, 0x7fc0_0000);
    try expectF(0xff80_0000, 0x7fc0_0000);
    try expectF(0x7f80_1234, 0x7fc0_1234);
    try expectF(0xff80_1234, 0xffc0_1234);
    try expectF(0x7fc1_2345, 0x7fc1_2345);
    try expectF(0xffc1_2345, 0xffc1_2345);

    try expectD(0x7ff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectD(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectD(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectD(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectD(0x7ff8_2468_ace0_1357, 0x7ff8_2468_ace0_1357);
    try expectD(0xfff8_2468_ace0_1357, 0xfff8_2468_ace0_1357);
}

test "sinpi rounds representative ordinary arguments" {
    try expectH(0x3400, 0x39a8);
    try expectF(0x3e80_0000, 0x3f35_04f3);
    try expectD(0x3fd0_0000_0000_0000, 0x3fe6_a09e_667f_3bcd);
}

test "sinpi handles large exact half-integers with odd symmetry" {
    try expectD(0xc267_e73e_5a0f_b000, 0x3ff0_0000_0000_0000);
}

test "sinpi f64 preserves exact lattice values at the integral cutoff" {
    try expectD(0x432f_ffff_ffff_ffff, 0xbff0_0000_0000_0000);
    try expectD(0xc32f_ffff_ffff_ffff, 0x3ff0_0000_0000_0000);
    try expectD(0x4330_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectD(0xc330_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectD(0x4330_0000_0000_0001, 0x0000_0000_0000_0000);
    try expectD(0xc330_0000_0000_0001, 0x8000_0000_0000_0000);
    try expectD(0x7fef_ffff_ffff_ffff, 0x0000_0000_0000_0000);
    try expectD(0xffef_ffff_ffff_ffff, 0x8000_0000_0000_0000);
}
