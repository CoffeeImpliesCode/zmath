const std = @import("std");
const alias = @import("scalar_logp1_alias.zig");
const log1p_impl = @import("scalar_log1p.zig");
const log1pd_impl = @import("scalar_log1pd.zig");

fn expectLogp1hMatchesLog1p(input_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const expected_bits: u16 = @bitCast(log1p_impl.log1ph(input));
    const actual_bits: u16 = @bitCast(alias.logp1h(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectLogp1fMatchesLog1p(input_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const expected_bits: u32 = @bitCast(log1p_impl.log1pf(input));
    const actual_bits: u32 = @bitCast(alias.logp1f(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectLogp1dMatchesLog1p(input_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const expected_bits: u64 = @bitCast(log1pd_impl.log1pd(input));
    const actual_bits: u64 = @bitCast(alias.logp1d(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "logp1 aliases preserve format-specific signatures" {
    comptime {
        std.debug.assert(@TypeOf(alias.logp1h) == fn (f16) f16);
        std.debug.assert(@TypeOf(alias.logp1f) == fn (f32) f32);
        std.debug.assert(@TypeOf(alias.logp1d) == fn (f64) f64);
    }
}

test "logp1h preserves log1ph bits for zeros finite values infinities and NaNs" {
    const input_bits = [_]u16{
        0x0000,
        0x8000,
        0x0001,
        0x8001,
        0x3800,
        0xb800,
        0x3c00,
        0xbc00,
        0x4000,
        0x7c00,
        0xfc00,
        0x7e55,
        0xfe55,
        0x7c55,
        0xfc55,
    };
    for (input_bits) |bits| try expectLogp1hMatchesLog1p(bits);
}

test "logp1f preserves log1pf bits for zeros finite values infinities and NaNs" {
    const input_bits = [_]u32{
        0x0000_0000,
        0x8000_0000,
        0x0000_0001,
        0x8000_0001,
        0x3f00_0000,
        0xbf00_0000,
        0x3f80_0000,
        0xbf80_0000,
        0x4000_0000,
        0x7f80_0000,
        0xff80_0000,
        0x7fc1_2345,
        0xffc1_2345,
        0x7f81_2345,
        0xff81_2345,
    };
    for (input_bits) |bits| try expectLogp1fMatchesLog1p(bits);
}

test "logp1d preserves log1pd bits for zeros finite values infinities and NaNs" {
    const input_bits = [_]u64{
        0x0000_0000_0000_0000,
        0x8000_0000_0000_0000,
        0x0000_0000_0000_0001,
        0x8000_0000_0000_0001,
        0x3fe0_0000_0000_0000,
        0xbfe0_0000_0000_0000,
        0x3ff0_0000_0000_0000,
        0xbff0_0000_0000_0000,
        0x4000_0000_0000_0000,
        0x7ff0_0000_0000_0000,
        0xfff0_0000_0000_0000,
        0x7ff8_1234_5678_9abc,
        0xfff8_1234_5678_9abc,
        0x7ff0_1234_5678_9abc,
        0xfff0_1234_5678_9abc,
    };
    for (input_bits) |bits| try expectLogp1dMatchesLog1p(bits);
}
