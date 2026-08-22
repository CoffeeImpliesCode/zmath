//! Exhaustively verifies every binary16 encoding through the public H bit API.
const std = @import("std");
const ieee = @import("ieee.zig");

test "fromBitsH and toBitsH preserve every H encoding" {
    var raw_bits: u32 = 0;
    while (raw_bits < 0x1_0000) : (raw_bits += 1) {
        const bits: u16 = @intCast(raw_bits);
        const value = ieee.fromBitsH(bits);

        try std.testing.expectEqual(bits, ieee.toBitsH(value));
    }
}

test "toBitsH and fromBitsH preserve every H value encoding" {
    var raw_bits: u32 = 0;
    while (raw_bits < 0x1_0000) : (raw_bits += 1) {
        const bits: u16 = @intCast(raw_bits);
        const value: f16 = @bitCast(bits);
        const round_tripped = ieee.fromBitsH(ieee.toBitsH(value));

        try std.testing.expectEqual(bits, @as(u16, @bitCast(round_tripped)));
    }
}

test "classifyH is stable through every H bit round trip" {
    var raw_bits: u32 = 0;
    while (raw_bits < 0x1_0000) : (raw_bits += 1) {
        const bits: u16 = @intCast(raw_bits);
        const value: f16 = @bitCast(bits);
        const round_tripped = ieee.fromBitsH(ieee.toBitsH(value));

        try std.testing.expectEqual(
            ieee.classifyH(value),
            ieee.classifyH(round_tripped),
        );
    }
}
