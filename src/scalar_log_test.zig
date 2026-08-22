const std = @import("std");
const scalar = @import("scalar.zig");

fn expectLogdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.logd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "logd has the nearest-even binary64 signature" {
    comptime {
        std.debug.assert(@TypeOf(scalar.logd) == fn (f64) f64);
    }
}

test "logd returns exact special values" {
    try expectLogdBits(0x0000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLogdBits(0x8000_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLogdBits(0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectLogdBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectLogdBits(0xbff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectLogdBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
}

test "logd quiets NaNs while preserving sign and payload" {
    try expectLogdBits(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLogdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLogdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "logd handles normal and subnormal boundaries" {
    try expectLogdBits(0x0000_0000_0000_0001, 0xc087_4385_446d_71c3);
    try expectLogdBits(0x0010_0000_0000_0000, 0xc086_232b_dd7a_bcd2);
    try expectLogdBits(0x3fe0_0000_0000_0000, 0xbfe6_2e42_fefa_39ef);
    try expectLogdBits(0x4000_0000_0000_0000, 0x3fe6_2e42_fefa_39ef);
    try expectLogdBits(0x7fef_ffff_ffff_ffff, 0x4086_2e42_fefa_39ef);
}

test "logd rounds values adjacent to one" {
    try expectLogdBits(0x3fef_ffff_ffff_ffff, 0xbca0_0000_0000_0000);
    try expectLogdBits(0x3ff0_0000_0000_0001, 0x3caf_ffff_ffff_ffff);
}

test "logd resolves CORE-MATH hard cases" {
    try expectLogdBits(0x3ffa_6ae5_1423_26b5, 0x3fe0_0bcc_31eb_ded7);
    try expectLogdBits(0x3ffa_6ba9_5b15_20f9, 0x3fe0_0cb9_ef61_45e1);
    try expectLogdBits(0x400f_fcd1_f1e0_c16f, 0x3ff6_2cab_e3b0_07f3);
    try expectLogdBits(0x47c7_1f7c_59ed_e8e0, 0x4055_c0be_a603_6ab8);
}
