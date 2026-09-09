const std = @import("std");
const scalar = @import("scalar.zig");

fn expectTanhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.tanh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectTanfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.tanf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectTandBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.tand(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "Tan scalar functions expose nearest-even signatures" {
    comptime std.debug.assert(@TypeOf(scalar.tanh) == fn (f16) f16);
    comptime std.debug.assert(@TypeOf(scalar.tanf) == fn (f32) f32);
    comptime std.debug.assert(@TypeOf(scalar.tand) == fn (f64) f64);
}

test "Tan scalar functions preserve signed zero and tiny values" {
    try expectTanhBits(0x0000, 0x0000);
    try expectTanhBits(0x8000, 0x8000);
    try expectTanhBits(0x0001, 0x0001);
    try expectTanhBits(0x8001, 0x8001);
    try expectTanfBits(0x0000_0000, 0x0000_0000);
    try expectTanfBits(0x8000_0000, 0x8000_0000);
    try expectTanfBits(0x0000_0001, 0x0000_0001);
    try expectTanfBits(0x8000_0001, 0x8000_0001);
    try expectTandBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectTandBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectTandBits(0x0000_0000_0000_0001, 0x0000_0000_0000_0001);
    try expectTandBits(0x8000_0000_0000_0001, 0x8000_0000_0000_0001);
}

test "Tan scalar functions use stable infinity and NaN results" {
    try expectTanhBits(0x7c00, 0x7e00);
    try expectTanhBits(0xfc00, 0x7e00);
    try expectTanhBits(0x7c55, 0x7e55);
    try expectTanhBits(0xfc55, 0xfe55);
    try expectTanfBits(0x7f80_0000, 0x7fc0_0000);
    try expectTanfBits(0xff80_0000, 0x7fc0_0000);
    try expectTanfBits(0x7f81_2345, 0x7fc1_2345);
    try expectTanfBits(0xff81_2345, 0xffc1_2345);
    try expectTandBits(0x7ff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectTandBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectTandBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectTandBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "tanh handles ordinary, largest, and CORE-MATH exceptional inputs" {
    try expectTanhBits(0x3800, 0x385f);
    try expectTanhBits(0xb800, 0xb85f);
    try expectTanhBits(0x3c00, 0x3e3b);
    try expectTanhBits(0xbc00, 0xbe3b);
    try expectTanhBits(0x7bff, 0xc46c);
    try expectTanhBits(0xfbff, 0x446c);
    try expectTanhBits(0x55ed, 0x3911);
    try expectTanhBits(0xd5ed, 0xb911);
    try expectTanhBits(0x658c, 0x07e7);
    try expectTanhBits(0x6429, 0x05ed);
    try expectTanhBits(0x681e, 0xa070);
    try expectTanhBits(0xe74e, 0xbb7d);
    try expectTanhBits(0xe807, 0xc015);
    try expectTanhBits(0x6f4d, 0xbe1a);
}

test "tanf handles ordinary, huge, and CORE-MATH hard inputs" {
    try expectTanfBits(0x3f00_0000, 0x3f0b_da7b);
    try expectTanfBits(0xbf00_0000, 0xbf0b_da7b);
    try expectTanfBits(0x3f80_0000, 0x3fc7_5923);
    try expectTanfBits(0xbf80_0000, 0xbfc7_5923);
    try expectTanfBits(0x7f7f_ffff, 0xbf1c_9eca);
    try expectTanfBits(0xff7f_ffff, 0x3f1c_9eca);
    try expectTanfBits(0x3f8a_1f62, 0x3fee_fcfb);
    try expectTanfBits(0x4d56_d355, 0x3e74_0182);
    try expectTanfBits(0x57d7_b0ed, 0x3eb0_68e4);
    try expectTanfBits(0x5980_445e, 0x3fe5_0f68);
    try expectTanfBits(0x63fc_86fe, 0x3f2c_bfce);
    try expectTanfBits(0x6a66_2711, 0xc079_9ac2);
    try expectTanfBits(0x6ad3_6709, 0xbf62_b097);
    try expectTanfBits(0x72b5_05bb, 0xbff2_150f);
    try expectTanfBits(0xbf8a_1f62, 0xbfee_fcfb);
    try expectTanfBits(0xcd56_d355, 0xbe74_0182);
    try expectTanfBits(0xd7d7_b0ed, 0xbeb0_68e4);
    try expectTanfBits(0xd980_445e, 0xbfe5_0f68);
    try expectTanfBits(0xe3fc_86fe, 0xbf2c_bfce);
    try expectTanfBits(0xea66_2711, 0x4079_9ac2);
    try expectTanfBits(0xead3_6709, 0x3f62_b097);
    try expectTanfBits(0xf2b5_05bb, 0x3ff2_150f);
}

test "tand handles ordinary, huge, and CORE-MATH accurate inputs" {
    try expectTandBits(0x3fe0_0000_0000_0000, 0x3fe1_7b4f_5bf3_474a);
    try expectTandBits(0xbfe0_0000_0000_0000, 0xbfe1_7b4f_5bf3_474a);
    try expectTandBits(0x3ff0_0000_0000_0000, 0x3ff8_eb24_5cbe_e3a6);
    try expectTandBits(0xbff0_0000_0000_0000, 0xbff8_eb24_5cbe_e3a6);
    try expectTandBits(0x7fe0_0000_0000_0000, 0xbfe5_ce6b_4c0d_02a3);
    try expectTandBits(0x3e9d_ffff_ffff_ff1f, 0x3e9e_0000_0000_0151);
    try expectTandBits(0x3ead_ffff_ffff_fc7c, 0x3eae_0000_0000_0546);
}

test "tand handles reduced zero and pole neighborhoods" {
    try expectTandBits(0x3ff9_21fb_5444_2d18, 0x434d_0296_7c31_cdb5);
    try expectTandBits(0xbff9_21fb_5444_2d18, 0xc34d_0296_7c31_cdb5);
    try expectTandBits(0x4009_21fb_5444_2d18, 0xbca1_a626_3314_5c07);
    try expectTandBits(0xc009_21fb_5444_2d18, 0x3ca1_a626_3314_5c07);
}
