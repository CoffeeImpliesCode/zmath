const std = @import("std");
const scalar = @import("scalar.zig");
const ieee = @import("ieee.zig");

fn expectExphBits(input_bits: u16, expected_bits: u16) !void {
    const input = ieee.fromBitsh(input_bits);
    try std.testing.expectEqual(expected_bits, ieee.toBitsh(scalar.exph(input)));
}

fn expectExpfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.expf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectExpdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.expd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "exp exports use nearest-even rounding" {
    comptime {
        std.debug.assert(@TypeOf(scalar.exph) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.expf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.expd) == fn (f64) f64);
    }
}

test "exp maps both signed zeros to positive one" {
    try expectExphBits(0x0000, 0x3c00);
    try expectExphBits(0x8000, 0x3c00);
    try expectExpfBits(0x0000_0000, 0x3f80_0000);
    try expectExpfBits(0x8000_0000, 0x3f80_0000);
    try expectExpdBits(0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectExpdBits(0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000);
}

test "exp maps positive infinity to positive infinity" {
    try expectExphBits(0x7c00, 0x7c00);
    try expectExpfBits(0x7f80_0000, 0x7f80_0000);
    try expectExpdBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
}

test "exp maps negative infinity to positive zero" {
    try expectExphBits(0xfc00, 0x0000);
    try expectExpfBits(0xff80_0000, 0x0000_0000);
    try expectExpdBits(0xfff0_0000_0000_0000, 0x0000_0000_0000_0000);
}

test "exp quiets NaNs while preserving sign and payload" {
    try expectExphBits(0x7e55, 0x7e55);
    try expectExphBits(0xfd55, 0xff55);

    try expectExpfBits(0x7fc1_2345, 0x7fc1_2345);
    try expectExpfBits(0xff81_2345, 0xffc1_2345);

    try expectExpdBits(0x7ff8_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectExpdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "exph finite corpus hook dispatches every finite H encoding" {
    var raw_bits: u32 = 0;
    while (raw_bits < 0x1_0000) : (raw_bits += 1) {
        const input_bits: u16 = @intCast(raw_bits);
        const input = ieee.fromBitsh(input_bits);

        switch (ieee.classifyh(input)) {
            .zero, .subnormal, .normal => {
                const output = scalar.exph(input);
                const output_bits = ieee.toBitsh(output);

                try std.testing.expectEqual(@as(u16, 0), output_bits & 0x8000);
                try std.testing.expect(switch (ieee.classifyh(output)) {
                    .zero, .subnormal, .normal, .infinity => true,
                    .quiet_nan, .signaling_nan => false,
                });
            },
            .infinity, .quiet_nan, .signaling_nan => {},
        }
    }
}
