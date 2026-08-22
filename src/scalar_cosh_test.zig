const std = @import("std");
const scalar = @import("scalar.zig");

fn expectCoshhBits(input_bits: u16, expected_bits: u16) !void {
    const input: f16 = @bitCast(input_bits);
    const actual_bits: u16 = @bitCast(scalar.coshh(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectCoshfBits(input_bits: u32, expected_bits: u32) !void {
    const input: f32 = @bitCast(input_bits);
    const actual_bits: u32 = @bitCast(scalar.coshf(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

fn expectCoshdBits(input_bits: u64, expected_bits: u64) !void {
    const input: f64 = @bitCast(input_bits);
    const actual_bits: u64 = @bitCast(scalar.coshd(input));
    try std.testing.expectEqual(expected_bits, actual_bits);
}

test "hyperbolic Cosh exports nearest-even H F and D signatures" {
    comptime {
        std.debug.assert(@TypeOf(scalar.coshh) == fn (f16) f16);
        std.debug.assert(@TypeOf(scalar.coshf) == fn (f32) f32);
        std.debug.assert(@TypeOf(scalar.coshd) == fn (f64) f64);
    }
}

test "hyperbolic Cosh maps signed zero and tiny values to one" {
    try expectCoshhBits(0x0000, 0x3c00);
    try expectCoshhBits(0x8000, 0x3c00);
    try expectCoshhBits(0x0001, 0x3c00);

    try expectCoshfBits(0x0000_0000, 0x3f80_0000);
    try expectCoshfBits(0x8000_0000, 0x3f80_0000);
    try expectCoshfBits(0x0000_0001, 0x3f80_0000);

    try expectCoshdBits(0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectCoshdBits(0x8000_0000_0000_0000, 0x3ff0_0000_0000_0000);
    try expectCoshdBits(0x0000_0000_0000_0001, 0x3ff0_0000_0000_0000);
}

test "hyperbolic Cosh rounds representative finite values" {
    try expectCoshhBits(0x3800, 0x3c83); // cosh(0.5)
    try expectCoshhBits(0x3c00, 0x3e2c); // cosh(1)
    try expectCoshhBits(0x4000, 0x4386); // cosh(2)

    try expectCoshfBits(0x3f00_0000, 0x3f90_560c); // cosh(0.5)
    try expectCoshfBits(0x3f80_0000, 0x3fc5_83ab); // cosh(1)
    try expectCoshfBits(0x4000_0000, 0x4070_c7d0); // cosh(2)

    try expectCoshdBits(0x3fe0_0000_0000_0000, 0x3ff2_0ac1_862a_e8d0); // cosh(0.5)
    try expectCoshdBits(0x3ff0_0000_0000_0000, 0x3ff8_b075_51d9_f550); // cosh(1)
    try expectCoshdBits(0x4024_0000_0000_0000, 0x40c5_829d_d053_712d); // cosh(10)
}

test "hyperbolic Cosh handles overflow boundaries" {
    try expectCoshhBits(0x49e4, 0x7bfb);
    try expectCoshhBits(0x49e5, 0x7c00);

    try expectCoshfBits(0x42b2_d4fc, 0x7f7f_ffec);
    try expectCoshfBits(0x42b2_d4fd, 0x7f80_0000);

    try expectCoshdBits(0x4086_33ce_8fb9_f87e, 0x7ff0_0000_0000_0000);
}

test "hyperbolic Cosh maps infinities to positive infinity" {
    try expectCoshhBits(0x7c00, 0x7c00);
    try expectCoshhBits(0xfc00, 0x7c00);
    try expectCoshfBits(0x7f80_0000, 0x7f80_0000);
    try expectCoshfBits(0xff80_0000, 0x7f80_0000);
    try expectCoshdBits(0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
    try expectCoshdBits(0xfff0_0000_0000_0000, 0x7ff0_0000_0000_0000);
}

test "hyperbolic Cosh quiets NaNs while preserving sign and payload" {
    try expectCoshhBits(0x7c55, 0x7e55);
    try expectCoshhBits(0xfc55, 0xfe55);
    try expectCoshfBits(0x7f81_2345, 0x7fc1_2345);
    try expectCoshfBits(0xff81_2345, 0xffc1_2345);
    try expectCoshdBits(0x7ff0_1234_5678_9abc, 0x7ff8_1234_5678_9abc);
    try expectCoshdBits(0xfff0_1234_5678_9abc, 0xfff8_1234_5678_9abc);
}

test "coshd near-zero double-double path matches CORE-MATH" {
    try expectCoshdBits(0x3fb6_0b25_aa47_d147, 0x3ff0_0f31_beb5_db92);
    try expectCoshdBits(0xbfb6_0b25_aa47_d147, 0x3ff0_0f31_beb5_db92);
}

test "coshh is even for every finite binary16 input" {
    var magnitude: u16 = 0;
    while (magnitude < 0x7c00) : (magnitude += 1) {
        const positive_bits: u16 = @bitCast(scalar.coshh(@bitCast(magnitude)));
        const negative_bits: u16 = @bitCast(scalar.coshh(@bitCast(magnitude | 0x8000)));
        try std.testing.expectEqual(positive_bits, negative_bits);
    }
}
