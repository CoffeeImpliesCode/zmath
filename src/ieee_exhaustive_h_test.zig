//! Exhaustively verifies every binary16 encoding through the public H bit API.
const std = @import("std");
const ieee = @import("ieee.zig");

test "fromBitsh and toBitsh preserve every H encoding" {
    var raw_bits: u32 = 0;
    while (raw_bits < 0x1_0000) : (raw_bits += 1) {
        const bits: u16 = @intCast(raw_bits);
        const value = ieee.fromBitsh(bits);

        try std.testing.expectEqual(bits, ieee.toBitsh(value));
    }
}

test "toBitsh and fromBitsh preserve every H value encoding" {
    var raw_bits: u32 = 0;
    while (raw_bits < 0x1_0000) : (raw_bits += 1) {
        const bits: u16 = @intCast(raw_bits);
        const value: f16 = @bitCast(bits);
        const round_tripped = ieee.fromBitsh(ieee.toBitsh(value));

        try std.testing.expectEqual(bits, @as(u16, @bitCast(round_tripped)));
    }
}

test "classifyh is stable through every H bit round trip" {
    var raw_bits: u32 = 0;
    while (raw_bits < 0x1_0000) : (raw_bits += 1) {
        const bits: u16 = @intCast(raw_bits);
        const value: f16 = @bitCast(bits);
        const round_tripped = ieee.fromBitsh(ieee.toBitsh(value));

        try std.testing.expectEqual(
            ieee.classifyh(value),
            ieee.classifyh(round_tripped),
        );
    }
}
