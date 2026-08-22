const std = @import("std");
const scalar = @import("scalar.zig");

fn expectLog1phBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.log1ph(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectLog1pfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.log1pf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectLog1pdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.log1pd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "log1p exposes nearest-even H, F, and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.log1ph) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.log1pf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.log1pd) == fn (f64) f64);
    }
}

test "log1ph handles exact special values" {
    try expectLog1phBits(0x0000, 0x0000);
    try expectLog1phBits(0x8000, 0x8000);
    try expectLog1phBits(0xbc00, 0xfc00);
    try expectLog1phBits(0xc000, 0x7e00);
    try expectLog1phBits(0xfc00, 0x7e00);
    try expectLog1phBits(0x7c00, 0x7c00);
}

test "log1ph quiets NaNs while preserving sign and payload" {
    try expectLog1phBits(0x7e55, 0x7e55);
    try expectLog1phBits(0x7c55, 0x7e55);
    try expectLog1phBits(0xfe55, 0xfe55);
    try expectLog1phBits(0xfc55, 0xfe55);
}

test "log1ph covers tiny ordinary and CORE-MATH hard cases" {
    try expectLog1phBits(0x0001, 0x0001);
    try expectLog1phBits(0x8001, 0x8001);
    try expectLog1phBits(0x3c00, 0x398c);
    try expectLog1phBits(0xb800, 0xb98c);
    try expectLog1phBits(0x9dfd, 0x9e01);
    try expectLog1phBits(0xa7db, 0xa7fb);
    try expectLog1phBits(0x1e03, 0x1dff);
    try expectLog1phBits(0x2f64, 0x2eff);
}

test "log1pf handles exact special values" {
    try expectLog1pfBits(0x0000_0000, 0x0000_0000);
    try expectLog1pfBits(0x8000_0000, 0x8000_0000);
    try expectLog1pfBits(0xbf80_0000, 0xff80_0000);
    try expectLog1pfBits(0xc000_0000, 0x7fc0_0000);
    try expectLog1pfBits(0xff80_0000, 0x7fc0_0000);
    try expectLog1pfBits(0x7f80_0000, 0x7f80_0000);
}

test "log1pf quiets NaNs while preserving sign and payload" {
    try expectLog1pfBits(0x7fc1_2345, 0x7fc1_2345);
    try expectLog1pfBits(0x7f81_2345, 0x7fc1_2345);
    try expectLog1pfBits(0xffc1_2345, 0xffc1_2345);
    try expectLog1pfBits(0xff81_2345, 0xffc1_2345);
}

test "log1pf covers tiny reduction boundaries and ordinary values" {
    try expectLog1pfBits(0x0000_0001, 0x0000_0001);
    try expectLog1pfBits(0x8000_0001, 0x8000_0001);
    try expectLog1pfBits(0x3f80_0000, 0x3f31_7218);
    try expectLog1pfBits(0xbf00_0000, 0xbf31_7218);
    try expectLog1pfBits(0x3c87_ffff, 0x3c86_e228);
    try expectLog1pfBits(0x3c88_0000, 0x3c86_e229);
    try expectLog1pfBits(0x4f00_0000, 0x41ab_e687);
    try expectLog1pfBits(0xbf7f_ffff, 0xc185_1592);
}

test "log1pd handles exact special values" {
    try expectLog1pdBits(0x0000_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectLog1pdBits(0x8000_0000_0000_0000, 0x8000_0000_0000_0000);
    try expectLog1pdBits(0xbff0_0000_0000_0000, 0xfff0_0000_0000_0000);
    try expectLog1pdBits(0xc000_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectLog1pdBits(0xfff0_0000_0000_0000, 0x7ff8_0000_0000_0000);
    try expectLog1pdBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
}

test "log1pd quiets NaNs while preserving sign and payload" {
    try expectLog1pdBits(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLog1pdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectLog1pdBits(0xfff8_1234_5678_9abc, 0xfff8_1234_5678_9abc);
    try expectLog1pdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "log1pd covers tiny reduction boundaries and ordinary values" {
    try expectLog1pdBits(0x0000_0000_0000_0001, 0x0000_0000_0000_0001);
    try expectLog1pdBits(0x8000_0000_0000_0001, 0x8000_0000_0000_0001);
    try expectLog1pdBits(0x3ff0_0000_0000_0000, 0x3fe6_2e42_fefa_39ef);
    try expectLog1pdBits(0xbfe0_0000_0000_0000, 0xbfe6_2e42_fefa_39ef);
    try expectLog1pdBits(0x3faf_ffff_ffff_ffff, 0x3faf_0a30_c011_62a5);
    try expectLog1pdBits(0x3fb0_0000_0000_0000, 0x3faf_0a30_c011_62a6);
    try expectLog1pdBits(0x4340_0000_0000_0000, 0x4042_5e4f_7b27_37fa);
    try expectLog1pdBits(0xbfef_ffff_ffff_ffff, 0xc042_5e4f_7b27_37fa);
}
