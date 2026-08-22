const std = @import("std");
const scalar = @import("scalar.zig");

fn expectExp2hBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.exp2h(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectExp2fBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.exp2f(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectExp2dBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.exp2d(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "exp2 exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.exp2h) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.exp2f) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.exp2d) == fn (f64) f64);
    }
}

test "exp2 maps signed zero and exact integer powers" {
    try expectExp2hBits(0x0000, 0x3c00);
    try expectExp2hBits(0x8000, 0x3c00);
    try expectExp2hBits(0x3c00, 0x4000);
    try expectExp2hBits(0xbc00, 0x3800);

    try expectExp2fBits(0x0000_0000, 0x3f80_0000);
    try expectExp2fBits(0x8000_0000, 0x3f80_0000);
    try expectExp2fBits(0x3f80_0000, 0x4000_0000);
    try expectExp2fBits(0xbf80_0000, 0x3f00_0000);

    try expectExp2dBits(0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectExp2dBits(0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectExp2dBits(0x3ff0_0000_0000_0000, 0x4000_0000_0000_0000);
    try expectExp2dBits(0xbff0_0000_0000_0000, 0x3fe0_0000_0000_0000);
}

test "exp2 handles exact subnormal powers and nearest-even underflow" {
    try expectExp2hBits(0xce00, 0x0001);
    try expectExp2hBits(0xce40, 0x0000);

    try expectExp2fBits(0xc315_0000, 0x0000_0001);
    try expectExp2fBits(0xc316_0000, 0x0000_0000);

    try expectExp2dBits(0xc090_c800_0000_0000, 0x0000_0000_0000_0001);
    try expectExp2dBits(0xc090_cc00_0000_0000, 0x0000_0000_0000_0000);
}

test "exp2 handles finite overflow boundaries" {
    try expectExp2hBits(0x4b80, 0x7800);
    try expectExp2hBits(0x4c00, 0x7c00);

    try expectExp2fBits(0x42fe_0000, 0x7f00_0000);
    try expectExp2fBits(0x4300_0000, 0x7f80_0000);

    try expectExp2dBits(0x408f_f800_0000_0000, 0x7fe0_0000_0000_0000);
    try expectExp2dBits(0x4090_0000_0000_0000, 0x7ff0_0000_0000_0000);
}

test "exp2 handles infinities and quiets NaNs with stable payloads" {
    try expectExp2hBits(0x7c00, 0x7c00);
    try expectExp2hBits(0xfc00, 0x0000);
    try expectExp2hBits(0x7c55, 0x7e55);
    try expectExp2hBits(0xfc55, 0xfe55);

    try expectExp2fBits(0x7f80_0000, 0x7f80_0000);
    try expectExp2fBits(0xff80_0000, 0x0000_0000);
    try expectExp2fBits(0x7f81_2345, 0x7fc1_2345);
    try expectExp2fBits(0xff81_2345, 0xffc1_2345);

    try expectExp2dBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectExp2dBits(0xfff0_0000_0000_0000, 0x0000_0000_0000_0000);
    try expectExp2dBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectExp2dBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "exp2f covers CORE-MATH accurate fallback inputs" {
    try expectExp2fBits(0x3b42_9d37, 0x3f80_4385);
    try expectExp2fBits(0xbcf3_a937, 0x3f7a_c6b1);
    try expectExp2fBits(0xb8d3_d026, 0x3f7f_fb69);
}

test "exp2d covers CORE-MATH database hard cases" {
    try expectExp2dBits(0x3f5e_4596_526b_f94d, 0x3ff0_053f_c2ec_2b53);
    try expectExp2dBits(0x3f98_d040_898b_73f5, 0x3ff0_4560_ec7b_6c8d);
    try expectExp2dBits(0x3ff9_f1a7_d355_cb4f, 0x4008_9d94_8a94_fe17);
    try expectExp2dBits(0xbf24_3c1c_ea9b_d4d9, 0x3fef_ff1f_99c8_5951);
    try expectExp2dBits(0xbf71_fb57_e199_6e26, 0x3fef_e71c_1c5f_e84f);
    try expectExp2dBits(0xbff6_0e58_2caa_34b1, 0x3fd8_9d94_8a94_fe17);
}

test "exp2h finite corpus hook visits every finite encoding" {
    var raw_bits: u32 = 0;
    while (raw_bits < 0x1_0000) : (raw_bits += 1) {
        const input_bits: u16 = @intCast(raw_bits);
        const magnitude = input_bits & 0x7fff;
        if (magnitude >= 0x7c00) continue;

        const output_bits: u16 = @bitCast(scalar.exp2h(@bitCast(input_bits)));
        try std.testing.expectEqual(@as(u16, 0), output_bits & 0x8000);
        try std.testing.expect((output_bits & 0x7c00) != 0x7c00 or (output_bits & 0x03ff) == 0);
    }
}
