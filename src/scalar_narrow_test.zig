const std = @import("std");
const narrow = @import("scalar_narrow.zig");

fn h(bits: u16) f16 {
    return @bitCast(bits);
}

fn f(bits: u32) f32 {
    return @bitCast(bits);
}

fn d(bits: u64) f64 {
    return @bitCast(bits);
}

fn expectH(expected: u16, actual: f16) !void {
    try std.testing.expectEqual(expected, @as(u16, @bitCast(actual)));
}

fn expectF(expected: u32, actual: f32) !void {
    try std.testing.expectEqual(expected, @as(u32, @bitCast(actual)));
}

test "all twenty-four narrowing symbols return their declared destination format" {
    try expectF(0x3fc0_0000, narrow.fadd(d(0x3ff0_0000_0000_0000), d(0x3fe0_0000_0000_0000)));
    try expectF(0x3f80_0000, narrow.fsub(d(0x3ff8_0000_0000_0000), d(0x3fe0_0000_0000_0000)));
    try expectF(0x4040_0000, narrow.fmul(d(0x3ff8_0000_0000_0000), d(0x4000_0000_0000_0000)));
    try expectF(0x3fc0_0000, narrow.fdiv(d(0x4008_0000_0000_0000), d(0x4000_0000_0000_0000)));
    try expectF(0x4060_0000, narrow.ffma(
        d(0x3ff8_0000_0000_0000),
        d(0x4000_0000_0000_0000),
        d(0x3fe0_0000_0000_0000),
    ));
    try expectF(0x4000_0000, narrow.fsqrt(d(0x4010_0000_0000_0000)));

    try expectH(0x3e00, narrow.f16addf32(f(0x3f80_0000), f(0x3f00_0000)));
    try expectH(0x3c00, narrow.f16subf32(f(0x3fc0_0000), f(0x3f00_0000)));
    try expectH(0x4200, narrow.f16mulf32(f(0x3fc0_0000), f(0x4000_0000)));
    try expectH(0x3e00, narrow.f16divf32(f(0x4040_0000), f(0x4000_0000)));
    try expectH(0x4300, narrow.f16fmaf32(f(0x3fc0_0000), f(0x4000_0000), f(0x3f00_0000)));
    try expectH(0x4000, narrow.f16sqrtf32(f(0x4080_0000)));

    try expectH(0x3e00, narrow.f16addf64(d(0x3ff0_0000_0000_0000), d(0x3fe0_0000_0000_0000)));
    try expectH(0x3c00, narrow.f16subf64(d(0x3ff8_0000_0000_0000), d(0x3fe0_0000_0000_0000)));
    try expectH(0x4200, narrow.f16mulf64(d(0x3ff8_0000_0000_0000), d(0x4000_0000_0000_0000)));
    try expectH(0x3e00, narrow.f16divf64(d(0x4008_0000_0000_0000), d(0x4000_0000_0000_0000)));
    try expectH(0x4300, narrow.f16fmaf64(
        d(0x3ff8_0000_0000_0000),
        d(0x4000_0000_0000_0000),
        d(0x3fe0_0000_0000_0000),
    ));
    try expectH(0x4000, narrow.f16sqrtf64(d(0x4010_0000_0000_0000)));

    try expectF(0x3fc0_0000, narrow.f32addf64(d(0x3ff0_0000_0000_0000), d(0x3fe0_0000_0000_0000)));
    try expectF(0x3f80_0000, narrow.f32subf64(d(0x3ff8_0000_0000_0000), d(0x3fe0_0000_0000_0000)));
    try expectF(0x4040_0000, narrow.f32mulf64(d(0x3ff8_0000_0000_0000), d(0x4000_0000_0000_0000)));
    try expectF(0x3fc0_0000, narrow.f32divf64(d(0x4008_0000_0000_0000), d(0x4000_0000_0000_0000)));
    try expectF(0x4060_0000, narrow.f32fmaf64(
        d(0x3ff8_0000_0000_0000),
        d(0x4000_0000_0000_0000),
        d(0x3fe0_0000_0000_0000),
    ));
    try expectF(0x4000_0000, narrow.f32sqrtf64(d(0x4010_0000_0000_0000)));
}

test "tie-adjacent sums and differences do not round in the source format" {
    const one_f = f(0x3f80_0000);
    try expectH(0x3c00, narrow.f16addf32(one_f, f(0x39ff_ffff)));
    try expectH(0x3c00, narrow.f16addf32(one_f, f(0x3a00_0000)));
    try expectH(0x3c01, narrow.f16addf32(one_f, f(0x3a00_0001)));
    try expectH(0x3c00, narrow.f16subf32(one_f, f(0xb9ff_ffff)));
    try expectH(0x3c00, narrow.f16subf32(one_f, f(0xba00_0000)));
    try expectH(0x3c01, narrow.f16subf32(one_f, f(0xba00_0001)));

    const one_d = d(0x3ff0_0000_0000_0000);
    try expectH(0x3c00, narrow.f16addf64(one_d, d(0x3f3f_ffff_ffff_ffff)));
    try expectH(0x3c00, narrow.f16addf64(one_d, d(0x3f40_0000_0000_0000)));
    try expectH(0x3c01, narrow.f16addf64(one_d, d(0x3f40_0000_0000_0001)));
    try expectH(0x3c00, narrow.f16subf64(one_d, d(0xbf3f_ffff_ffff_ffff)));
    try expectH(0x3c00, narrow.f16subf64(one_d, d(0xbf40_0000_0000_0000)));
    try expectH(0x3c01, narrow.f16subf64(one_d, d(0xbf40_0000_0000_0001)));

    try expectF(0x3f80_0000, narrow.f32addf64(one_d, d(0x3e6f_ffff_ffff_ffff)));
    try expectF(0x3f80_0000, narrow.f32addf64(one_d, d(0x3e70_0000_0000_0000)));
    try expectF(0x3f80_0001, narrow.f32addf64(one_d, d(0x3e70_0000_0000_0001)));
    try expectF(0x3f80_0000, narrow.fsub(one_d, d(0xbe70_0000_0000_0000)));
    try expectF(0x3f80_0001, narrow.fsub(one_d, d(0xbe70_0000_0000_0001)));
    try expectF(0x3f80_0000, narrow.fsub(d(0x3ff0_0000_1000_0000), d(0x1)));
}

test "multiplication retains the exact product past source precision" {
    try expectH(0x3c01, narrow.f16mulf32(f(0x3f80_0fff), f(0x3f80_0001)));
    try expectH(0x3c01, narrow.f16mulf64(
        d(0x3ff0_01ff_ffff_ffff),
        d(0x3ff0_0000_0000_0001),
    ));
    try expectF(0x3f80_0001, narrow.f32mulf64(
        d(0x3ff0_0000_0fff_ffff),
        d(0x3ff0_0000_0000_0001),
    ));
    try expectF(0x3f80_0001, narrow.fmul(
        d(0x3ff0_0000_0fff_ffff),
        d(0x3ff0_0000_0000_0001),
    ));
}

test "division corrects the adjacent destination result from an exact residual" {
    try expectH(0x3c01, narrow.f16divf32(f(0x3f80_1800), f(0x3f80_07ff)));
    try expectH(0x3c01, narrow.f16divf64(
        d(0x3ff0_0200_0000_0800),
        d(0x3ff0_0000_0000_07ff),
    ));
    try expectF(0x3f80_0001, narrow.f32divf64(
        d(0x3ff0_0000_1100_0000),
        d(0x3ff0_0000_00ff_ffff),
    ));
    try expectF(0x3f80_0001, narrow.fdiv(
        d(0x3ff0_0000_1100_0000),
        d(0x3ff0_0000_00ff_ffff),
    ));
}

test "fused multiply-add retains the exact addend past source precision" {
    try expectH(0x3c01, narrow.f16fmaf32(
        f(0x3f80_0001),
        f(0x3f80_0001),
        f(0x39ff_e000),
    ));
    try expectH(0x3c01, narrow.f16fmaf64(
        d(0x3ff0_0000_0000_0001),
        d(0x3ff0_0000_0000_0001),
        d(0x3f3f_ffff_ffff_e000),
    ));
    try expectF(0x3f80_0001, narrow.f32fmaf64(
        d(0x3ff0_0000_0000_0001),
        d(0x3ff0_0000_0000_0001),
        d(0x3e6f_ffff_fc00_0000),
    ));
    try expectF(0x3f80_0001, narrow.ffma(
        d(0x3ff0_0000_0000_0001),
        d(0x3ff0_0000_0000_0001),
        d(0x3e6f_ffff_fc00_0000),
    ));
    try expectF(0x3f80_0000, narrow.ffma(
        d(0x3ff0_0000_1000_0000),
        d(0x3ff0_0000_0000_0000),
        d(0x8000_0000_0000_0001),
    ));
    try expectF(0x3f80_0000, narrow.f32fmaf64(
        d(0x3ff0_0000_1000_0000),
        d(0x3ff0_0000_0000_0000),
        d(0x8000_0000_0000_0001),
    ));
    try expectH(0x3c00, narrow.f16fmaf64(
        d(0x3ff0_0000_1000_0000),
        d(0x3ff0_0000_0000_0000),
        d(0x8000_0000_0000_0001),
    ));
    try expectH(0x3c00, narrow.f16fmaf32(
        f(0x3f80_0001),
        f(0x3f80_0000),
        f(0x8000_0001),
    ));
}

test "square root corrects the adjacent destination result from an exact residual" {
    try expectH(0x3c01, narrow.f16sqrtf32(f(0x3f80_2003)));
    try expectH(0x3c01, narrow.f16sqrtf64(d(0x3ff0_0400_4000_0001)));
    try expectF(0x3f80_0001, narrow.f32sqrtf64(d(0x3ff0_0000_2000_0011)));
    try expectF(0x3f80_0001, narrow.fsqrt(d(0x3ff0_0000_2000_0011)));
}

test "square root rounds the exact midpoint once at both parity choices" {
    // The source-format roots at these midpoints would lose the one-bit side.
    try expectH(0x3c00, narrow.f16sqrtf32(f(0x3f80_2001)));
    try expectH(0x3c00, narrow.f16sqrtf32(f(0x3f80_2002)));
    try expectH(0x3c01, narrow.f16sqrtf32(f(0x3f80_2003)));
    try expectH(0x3c01, narrow.f16sqrtf32(f(0x3f80_6011)));
    try expectH(0x3c02, narrow.f16sqrtf32(f(0x3f80_6012)));
    try expectH(0x3c02, narrow.f16sqrtf32(f(0x3f80_6013)));

    try expectH(0x3c00, narrow.f16sqrtf64(d(0x3ff0_0400_3fff_ffff)));
    try expectH(0x3c00, narrow.f16sqrtf64(d(0x3ff0_0400_4000_0000)));
    try expectH(0x3c01, narrow.f16sqrtf64(d(0x3ff0_0400_4000_0001)));
    try expectH(0x3c01, narrow.f16sqrtf64(d(0x3ff0_0c02_3fff_ffff)));
    try expectH(0x3c02, narrow.f16sqrtf64(d(0x3ff0_0c02_4000_0000)));
    try expectH(0x3c02, narrow.f16sqrtf64(d(0x3ff0_0c02_4000_0001)));

    try expectF(0x3f80_0000, narrow.f32sqrtf64(d(0x3ff0_0000_2000_000f)));
    try expectF(0x3f80_0000, narrow.f32sqrtf64(d(0x3ff0_0000_2000_0010)));
    try expectF(0x3f80_0001, narrow.f32sqrtf64(d(0x3ff0_0000_2000_0011)));
    try expectF(0x3f80_0001, narrow.f32sqrtf64(d(0x3ff0_0000_6000_008f)));
    try expectF(0x3f80_0002, narrow.f32sqrtf64(d(0x3ff0_0000_6000_0090)));
    try expectF(0x3f80_0002, narrow.f32sqrtf64(d(0x3ff0_0000_6000_0091)));
}

test "square root rounds once at zero and the least subnormal" {
    // 2^-50 and 2^-300 are exact squares of the destination midpoints.
    try expectH(0x0000, narrow.f16sqrtf32(f(0x267f_ffff)));
    try expectH(0x0000, narrow.f16sqrtf32(f(0x2680_0000)));
    try expectH(0x0001, narrow.f16sqrtf32(f(0x2680_0001)));
    try expectH(0x0000, narrow.f16sqrtf64(d(0x3ccf_ffff_ffff_ffff)));
    try expectH(0x0000, narrow.f16sqrtf64(d(0x3cd0_0000_0000_0000)));
    try expectH(0x0001, narrow.f16sqrtf64(d(0x3cd0_0000_0000_0001)));
    try expectF(0x0000_0000, narrow.f32sqrtf64(d(0x2d2f_ffff_ffff_ffff)));
    try expectF(0x0000_0000, narrow.f32sqrtf64(d(0x2d30_0000_0000_0000)));
    try expectF(0x0000_0001, narrow.f32sqrtf64(d(0x2d30_0000_0000_0001)));
}

test "square root rounds at the subnormal to normal boundary" {
    try expectH(0x03ff, narrow.f16sqrtf32(f(0x317f_c003)));
    try expectH(0x0400, narrow.f16sqrtf32(f(0x317f_c004)));
    try expectH(0x0400, narrow.f16sqrtf32(f(0x317f_c005)));
    try expectH(0x03ff, narrow.f16sqrtf64(d(0x3e2f_f800_7fff_ffff)));
    try expectH(0x0400, narrow.f16sqrtf64(d(0x3e2f_f800_8000_0000)));
    try expectH(0x0400, narrow.f16sqrtf64(d(0x3e2f_f800_8000_0001)));

    try expectF(0x007f_ffff, narrow.f32sqrtf64(d(0x302f_ffff_c000_001f)));
    try expectF(0x0080_0000, narrow.f32sqrtf64(d(0x302f_ffff_c000_0020)));
    try expectF(0x0080_0000, narrow.f32sqrtf64(d(0x302f_ffff_c000_0021)));
}

test "gradual underflow and overflow ties are rounded in the destination format" {
    const zero_f = f(0x0000_0000);
    try expectH(0x0000, narrow.f16addf32(f(0x3300_0000), zero_f));
    try expectH(0x0001, narrow.f16addf32(f(0x3300_0001), zero_f));
    try expectH(0x7bff, narrow.f16addf32(f(0x477f_efff), zero_f));
    try expectH(0x7c00, narrow.f16addf32(f(0x477f_f000), zero_f));

    const zero_d = d(0x0000_0000_0000_0000);
    try expectH(0x0000, narrow.f16addf64(d(0x3e60_0000_0000_0000), zero_d));
    try expectH(0x0001, narrow.f16addf64(d(0x3e60_0000_0000_0001), zero_d));
    try expectH(0x7bff, narrow.f16addf64(d(0x40ef_fdff_ffff_ffff), zero_d));
    try expectH(0x7c00, narrow.f16addf64(d(0x40ef_fe00_0000_0000), zero_d));

    try expectF(0x0000_0000, narrow.f32addf64(d(0x3690_0000_0000_0000), zero_d));
    try expectF(0x0000_0001, narrow.f32addf64(d(0x3690_0000_0000_0001), zero_d));
    try expectF(0x0000_0001, narrow.f32addf64(d(0x3690_0000_0000_0000), d(0x1)));
    try expectH(0x0001, narrow.f16addf64(d(0x3e60_0000_0000_0000), d(0x1)));
    try expectF(0x7f7f_ffff, narrow.f32addf64(d(0x47ef_ffff_efff_ffff), zero_d));
    try expectF(0x7f80_0000, narrow.f32addf64(d(0x47ef_ffff_f000_0000), zero_d));
    try expectF(0xff80_0000, narrow.fadd(d(0xc7ef_ffff_f000_0000), zero_d));
}

test "signed zero follows the corresponding exact IEC operation" {
    const positive_zero_f = f(0x0000_0000);
    const negative_zero_f = f(0x8000_0000);
    const two_f = f(0x4000_0000);
    try expectH(0x8000, narrow.f16addf32(negative_zero_f, negative_zero_f));
    try expectH(0x0000, narrow.f16addf32(positive_zero_f, negative_zero_f));
    try expectH(0x8000, narrow.f16subf32(negative_zero_f, positive_zero_f));
    try expectH(0x8000, narrow.f16mulf32(negative_zero_f, two_f));
    try expectH(0x8000, narrow.f16divf32(positive_zero_f, f(0xc000_0000)));
    try expectH(0x8000, narrow.f16fmaf32(negative_zero_f, two_f, negative_zero_f));
    try expectH(0x0000, narrow.f16fmaf32(negative_zero_f, two_f, positive_zero_f));
    try expectH(0x8000, narrow.f16sqrtf32(negative_zero_f));

    const positive_zero_d = d(0x0000_0000_0000_0000);
    const negative_zero_d = d(0x8000_0000_0000_0000);
    const two_d = d(0x4000_0000_0000_0000);
    try expectH(0x8000, narrow.f16addf64(negative_zero_d, negative_zero_d));
    try expectH(0x8000, narrow.f16subf64(negative_zero_d, positive_zero_d));
    try expectH(0x8000, narrow.f16mulf64(negative_zero_d, two_d));
    try expectH(0x8000, narrow.f16divf64(positive_zero_d, d(0xc000_0000_0000_0000)));
    try expectH(0x8000, narrow.f16fmaf64(negative_zero_d, two_d, negative_zero_d));
    try expectH(0x8000, narrow.f16sqrtf64(negative_zero_d));

    try expectF(0x8000_0000, narrow.f32addf64(negative_zero_d, negative_zero_d));
    try expectF(0x8000_0000, narrow.f32subf64(negative_zero_d, positive_zero_d));
    try expectF(0x8000_0000, narrow.f32mulf64(negative_zero_d, two_d));
    try expectF(0x8000_0000, narrow.f32divf64(positive_zero_d, d(0xc000_0000_0000_0000)));
    try expectF(0x8000_0000, narrow.f32fmaf64(negative_zero_d, two_d, negative_zero_d));
    try expectF(0x8000_0000, narrow.f32sqrtf64(negative_zero_d));
    try expectF(0x0000_0000, narrow.fsub(d(0x3ff0_0000_0000_0000), d(0x3ff0_0000_0000_0000)));
}

fn expectInvalidH32() !void {
    const zeros = [_]f32{ f(0x0000_0000), f(0x8000_0000) };
    const infinities = [_]f32{ f(0x7f80_0000), f(0xff80_0000) };
    for (infinities) |infinity| {
        for (zeros) |zero| {
            try expectH(0x7e00, narrow.f16mulf32(infinity, zero));
            try expectH(0x7e00, narrow.f16mulf32(zero, infinity));
            try expectH(0x7e00, narrow.f16fmaf32(infinity, zero, f(0x3f80_0000)));
            try expectH(0x7e00, narrow.f16fmaf32(zero, infinity, f(0x3f80_0000)));
        }
        for (infinities) |denominator| {
            try expectH(0x7e00, narrow.f16divf32(infinity, denominator));
        }
    }
    for (zeros) |zero| {
        for (zeros) |denominator| {
            try expectH(0x7e00, narrow.f16divf32(zero, denominator));
        }
    }
    try expectH(0x7e00, narrow.f16addf32(infinities[0], infinities[1]));
    try expectH(0x7e00, narrow.f16addf32(infinities[1], infinities[0]));
    try expectH(0x7e00, narrow.f16subf32(infinities[0], infinities[0]));
    try expectH(0x7e00, narrow.f16subf32(infinities[1], infinities[1]));
    try expectH(0x7e00, narrow.f16fmaf32(infinities[0], f(0x3f80_0000), infinities[1]));
    try expectH(0x7e00, narrow.f16fmaf32(infinities[1], f(0x3f80_0000), infinities[0]));
    try expectH(0x7e00, narrow.f16sqrtf32(f(0xff80_0000)));
    try expectH(0x7e00, narrow.f16sqrtf32(f(0xbf80_0000)));
}

fn expectInvalidH64() !void {
    const zeros = [_]f64{ d(0x0000_0000_0000_0000), d(0x8000_0000_0000_0000) };
    const infinities = [_]f64{ d(0x7ff0_0000_0000_0000), d(0xfff0_0000_0000_0000) };
    for (infinities) |infinity| {
        for (zeros) |zero| {
            try expectH(0x7e00, narrow.f16mulf64(infinity, zero));
            try expectH(0x7e00, narrow.f16mulf64(zero, infinity));
            try expectH(0x7e00, narrow.f16fmaf64(infinity, zero, d(0x3ff0_0000_0000_0000)));
            try expectH(0x7e00, narrow.f16fmaf64(zero, infinity, d(0x3ff0_0000_0000_0000)));
        }
        for (infinities) |denominator| {
            try expectH(0x7e00, narrow.f16divf64(infinity, denominator));
        }
    }
    for (zeros) |zero| {
        for (zeros) |denominator| {
            try expectH(0x7e00, narrow.f16divf64(zero, denominator));
        }
    }
    try expectH(0x7e00, narrow.f16addf64(infinities[0], infinities[1]));
    try expectH(0x7e00, narrow.f16addf64(infinities[1], infinities[0]));
    try expectH(0x7e00, narrow.f16subf64(infinities[0], infinities[0]));
    try expectH(0x7e00, narrow.f16subf64(infinities[1], infinities[1]));
    try expectH(0x7e00, narrow.f16fmaf64(infinities[0], d(0x3ff0_0000_0000_0000), infinities[1]));
    try expectH(0x7e00, narrow.f16fmaf64(infinities[1], d(0x3ff0_0000_0000_0000), infinities[0]));
    try expectH(0x7e00, narrow.f16sqrtf64(d(0xfff0_0000_0000_0000)));
    try expectH(0x7e00, narrow.f16sqrtf64(d(0xbff0_0000_0000_0000)));
}

fn expectInvalidF64() !void {
    const zeros = [_]f64{ d(0x0000_0000_0000_0000), d(0x8000_0000_0000_0000) };
    const infinities = [_]f64{ d(0x7ff0_0000_0000_0000), d(0xfff0_0000_0000_0000) };
    for (infinities) |infinity| {
        for (zeros) |zero| {
            try expectF(0x7fc0_0000, narrow.f32mulf64(infinity, zero));
            try expectF(0x7fc0_0000, narrow.f32mulf64(zero, infinity));
            try expectF(0x7fc0_0000, narrow.fmul(infinity, zero));
            try expectF(0x7fc0_0000, narrow.f32fmaf64(infinity, zero, d(0x3ff0_0000_0000_0000)));
            try expectF(0x7fc0_0000, narrow.f32fmaf64(infinity, zero, d(0x3ff0_0000_0000_0000)));
            try expectF(0x7fc0_0000, narrow.ffma(zero, infinity, d(0x3ff0_0000_0000_0000)));
        }
        for (infinities) |denominator| {
            try expectF(0x7fc0_0000, narrow.f32divf64(infinity, denominator));
            try expectF(0x7fc0_0000, narrow.fdiv(infinity, denominator));
        }
    }
    for (zeros) |zero| {
        for (zeros) |denominator| {
            try expectF(0x7fc0_0000, narrow.f32divf64(zero, denominator));
            try expectF(0x7fc0_0000, narrow.fdiv(zero, denominator));
        }
    }
    try expectF(0x7fc0_0000, narrow.f32addf64(infinities[0], infinities[1]));
    try expectF(0x7fc0_0000, narrow.f32subf64(infinities[0], infinities[0]));
    try expectF(0x7fc0_0000, narrow.fadd(infinities[1], infinities[0]));
    try expectF(0x7fc0_0000, narrow.fsub(infinities[1], infinities[1]));
    try expectF(0x7fc0_0000, narrow.f32fmaf64(infinities[0], d(0x3ff0_0000_0000_0000), infinities[1]));
    try expectF(0x7fc0_0000, narrow.ffma(infinities[1], d(0x3ff0_0000_0000_0000), infinities[0]));
    try expectF(0x7fc0_0000, narrow.f32sqrtf64(d(0xfff0_0000_0000_0000)));
    try expectF(0x7fc0_0000, narrow.fsqrt(d(0xbff0_0000_0000_0000)));
}

test "all infinity and zero invalid combinations construct canonical destination NaNs" {
    try expectInvalidH32();
    try expectInvalidH64();
    try expectInvalidF64();
}

test "NaN conversion quiets and selects the first NaN in argument order" {
    const signaling_f = f(0xffaa_a000);
    const quiet_f = f(0x7fd5_4000);
    try expectH(0xff55, narrow.f16addf32(signaling_f, quiet_f));
    try expectH(0x7eaa, narrow.f16addf32(quiet_f, signaling_f));
    try expectH(0xff55, narrow.f16divf32(f(0x3f80_0000), signaling_f));
    try expectH(0xff55, narrow.f16fmaf32(f(0x3f80_0000), signaling_f, quiet_f));
    try expectH(0x7eaa, narrow.f16fmaf32(f(0x7f80_0000), f(0x0000_0000), quiet_f));
    try expectH(0xff55, narrow.f16sqrtf32(signaling_f));

    const signaling_h_d = d(0xfff5_5400_0000_0000);
    const quiet_h_d = d(0x7ffa_a800_0000_0000);
    try expectH(0xff55, narrow.f16addf64(signaling_h_d, quiet_h_d));
    try expectH(0x7eaa, narrow.f16addf64(quiet_h_d, signaling_h_d));
    try expectH(0xff55, narrow.f16fmaf64(d(0x3ff0_0000_0000_0000), signaling_h_d, quiet_h_d));
    try expectH(0x7eaa, narrow.f16fmaf64(d(0x7ff0_0000_0000_0000), d(0), quiet_h_d));
    try expectH(0xff55, narrow.f16sqrtf64(signaling_h_d));

    const signaling_f_d = d(0xfff0_0246_8000_0000);
    const quiet_f_d = d(0x7ff8_08ac_e000_0000);
    try expectF(0xffc0_1234, narrow.f32addf64(signaling_f_d, quiet_f_d));
    try expectF(0x7fc0_4567, narrow.f32addf64(quiet_f_d, signaling_f_d));
    try expectF(0xffc0_1234, narrow.f32fmaf64(d(0x3ff0_0000_0000_0000), signaling_f_d, quiet_f_d));
    try expectF(0x7fc0_4567, narrow.ffma(d(0x7ff0_0000_0000_0000), d(0), quiet_f_d));
    try expectF(0xffc0_1234, narrow.fsqrt(signaling_f_d));
}

test "NaN conversion is stable across every narrowing operation" {
    const signaling_f = f(0xffaa_a000);
    try expectH(0xff55, narrow.f16addf32(signaling_f, f(0x3f80_0000)));
    try expectH(0xff55, narrow.f16subf32(signaling_f, f(0x3f80_0000)));
    try expectH(0xff55, narrow.f16mulf32(signaling_f, f(0x3f80_0000)));
    try expectH(0xff55, narrow.f16divf32(signaling_f, f(0x3f80_0000)));
    try expectH(0xff55, narrow.f16fmaf32(signaling_f, f(0x3f80_0000), f(0)));
    try expectH(0xff55, narrow.f16sqrtf32(signaling_f));

    const signaling_d = d(0xfff5_5400_0000_0000);
    try expectH(0xff55, narrow.f16addf64(signaling_d, d(0x3ff0_0000_0000_0000)));
    try expectH(0xff55, narrow.f16subf64(signaling_d, d(0x3ff0_0000_0000_0000)));
    try expectH(0xff55, narrow.f16mulf64(signaling_d, d(0x3ff0_0000_0000_0000)));
    try expectH(0xff55, narrow.f16divf64(signaling_d, d(0x3ff0_0000_0000_0000)));
    try expectH(0xff55, narrow.f16fmaf64(signaling_d, d(0x3ff0_0000_0000_0000), d(0)));
    try expectH(0xff55, narrow.f16sqrtf64(signaling_d));

    const signaling_f64 = d(0xfff0_0246_8000_0000);
    try expectF(0xffc0_1234, narrow.f32addf64(signaling_f64, d(0x3ff0_0000_0000_0000)));
    try expectF(0xffc0_1234, narrow.f32subf64(signaling_f64, d(0x3ff0_0000_0000_0000)));
    try expectF(0xffc0_1234, narrow.f32mulf64(signaling_f64, d(0x3ff0_0000_0000_0000)));
    try expectF(0xffc0_1234, narrow.f32divf64(signaling_f64, d(0x3ff0_0000_0000_0000)));
    try expectF(0xffc0_1234, narrow.f32fmaf64(signaling_f64, d(0x3ff0_0000_0000_0000), d(0)));
    try expectF(0xffc0_1234, narrow.f32sqrtf64(signaling_f64));
}
