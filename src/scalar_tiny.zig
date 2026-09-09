const std = @import("std");

pub const Kind = enum { ln2, ln10, inv_ln2, inv_ln10, pi, inv_pi, two_inv_sqrt_pi };
pub const magnitude_limit: u64 = 0x0020_0000_0000_0000;
const sign_mask: u64 = 0x8000_0000_0000_0000;
const significand_max: u64 = magnitude_limit - 1;

const Certificate = struct {
    coefficient: u128,
    lower_n: u64,
    lower_d: u64,
    upper_n: u64,
    upper_d: u64,
};

// Each coefficient is RNE(c * 2^126), checked with directed MPFR-512 bounds.
const certificates = [_]Certificate{
    .{ // ln2
        .coefficient = 0x2c5c_85fd_f473_de6a_f278_ece6_00fc_bdac,
        .lower_n = 3109807662917472,
        .lower_d = 2243252046704767,
        .upper_n = 9617789512135903,
        .upper_d = 6937768616736175,
    },
    .{ // ln10
        .coefficient = 0x935d_8ddd_aaa8_ac16_ea56_d62b_82d3_0a29,
        .lower_n = 41289067906310253,
        .lower_d = 8965807177319597,
        .upper_n = 227480160645689,
        .upper_d = 49396689255444,
    },
    .{ // inv_ln2
        .coefficient = 0x5c55_1d94_ae0b_f85d_df43_ff68_348e_9f44,
        .lower_n = 25507822420239933,
        .lower_d = 8840337596406535,
        .upper_n = 2243252046704767,
        .upper_d = 777451915729368,
    },
    .{ // inv_ln10
        .coefficient = 0x1bcb_7b15_26e5_0e32_a6ab_7555_f5a6_7b86,
        .lower_n = 197586757021776,
        .lower_d = 227480160645689,
        .upper_n = 7730889945933497,
        .upper_d = 8900515972542007,
    },
    .{ // pi
        .coefficient = 0xc90f_daa2_2168_c234_c4c6_628b_80dc_1cd1,
        .lower_n = 42516072084570011,
        .lower_d = 6766643033110660,
        .upper_n = 30246273033735921,
        .upper_d = 4813843863426169,
    },
    .{ // inv_pi
        .coefficient = 0x145f_306d_c9c8_82a5_3f84_eafa_3ea6_9bb8,
        .lower_n = 3769290217798865,
        .lower_d = 5920787228742393,
        .upper_n = 3905598339368982,
        .upper_d = 6134899525417045,
    },
    .{ // 2/sqrt(pi)
        .coefficient = 0x4837_5d41_0a6d_b446_b8ea_453f_b5ff_61a2,
        .lower_n = 7705853849415594,
        .lower_d = 3414567582477941,
        .upper_n = 19180872201022277,
        .upper_d = 8499302699107124,
    },
};

// x = sign * M * 2^-1074 for M < 2^53. The allowed mathematical remainder
// puts the result strictly between M*(C-1)/2^126 and M*(C+1)/2^126 in units
// of the least subnormal. Every output rounding midpoint is a half-integer.
// A midpoint in this interval would give n/M between (C-1)/2^125 and
// (C+1)/2^125. Adjacent Farey fractions bracket that interval below; every
// enclosed rational has denominator at least lower_d + upper_d > M.
comptime {
    std.debug.assert(certificates.len == @typeInfo(Kind).@"enum".fields.len);
    const denominator: u256 = @as(u256, 1) << 125;
    for (certificates) |certificate| {
        std.debug.assert(certificate.coefficient > 1);
        std.debug.assert(certificate.coefficient < std.math.maxInt(u128));
        std.debug.assert(certificate.lower_d > 0);
        std.debug.assert(certificate.upper_d > 0);
        std.debug.assert(certificate.lower_d <= significand_max);
        std.debug.assert(certificate.upper_d <= significand_max);
        std.debug.assert(certificate.lower_d + certificate.upper_d > significand_max);
        std.debug.assert(@as(u128, certificate.upper_n) * certificate.lower_d ==
            @as(u128, certificate.lower_n) * certificate.upper_d + 1);
        std.debug.assert(@as(u256, certificate.lower_n) * denominator <
            @as(u256, certificate.coefficient - 1) * certificate.lower_d);
        std.debug.assert(@as(u256, certificate.coefficient + 1) * certificate.upper_d <
            @as(u256, certificate.upper_n) * denominator);
    }
}

/// Round x*c in the tiny domain without subnormal floating-point arithmetic.
/// Requires |x| < 2^-1021. Also valid for an additional real remainder whose
/// magnitude is at most |x|*2^-128. The certificates exclude all ambiguities.
pub inline fn scale(comptime kind: Kind, value: f64) f64 {
    const bits: u64 = @bitCast(value);
    const magnitude = bits & ~sign_mask;
    std.debug.assert(magnitude < magnitude_limit);
    if (magnitude == 0) return value;

    const significand: u53 = @intCast(magnitude);
    const coefficient = comptime certificates[@intFromEnum(kind)].coefficient;
    const coefficient_low: u64 = @truncate(coefficient);
    const coefficient_high: u64 = @truncate(coefficient >> 64);
    const low_product = @as(u128, significand) * coefficient_low;
    const high_product = @as(u128, significand) * coefficient_high +% (low_product >> 64);
    const low: u64 = @truncate(low_product);
    const whole: u64 = @intCast(high_product >> 62);

    // Above the first normal binade, the result's unit is 2 or 4 subnormals.
    const extra: u2 = if (whole < @as(u64, 1) << 53) 0 else if (whole < @as(u64, 1) << 54) 1 else 2;
    const discarded: u7 = 62 + @as(u7, extra);
    const half: u128 = @as(u128, 1) << (discarded - 1);
    const integer: u64 = @intCast(high_product >> discarded);
    const sticky = (high_product & (half - 1)) != 0 or low != 0;
    const round_up = (high_product & half) != 0 and (sticky or integer & 1 != 0);
    const rounded = integer + @as(u64, @intFromBool(round_up));
    const encoded = rounded + (@as(u64, extra) << 52);
    return @bitCast(encoded | (bits & sign_mask));
}

const Case = struct { kind: Kind, input: u64, expected: u64 };

test "tiny scaling rounds close rational approximations" {
    // MPFR-512 lower and upper bounds agree on each target-format result.
    const cases = [_]Case{
        .{ .kind = .ln2, .input = 0x0007f839d162b07f, .expected = 0x0005862d4288b3b0 },
        .{ .kind = .ln2, .input = 0x0018a5dd21ef35af, .expected = 0x001115aa0fb29a6f },
        .{ .kind = .ln10, .input = 0x001fda5aa7edc4ad, .expected = 0x003256060835340e },
        .{ .kind = .ln10, .input = 0x00002ced101dbc14, .expected = 0x000067722d15991c },
        .{ .kind = .inv_ln2, .input = 0x001f683d7e20db07, .expected = 0x0026a7cead96898f },
        .{ .kind = .inv_ln2, .input = 0x0002c316a14459d8, .expected = 0x0003fc1ce8b1583f },
        .{ .kind = .inv_ln10, .input = 0x0000cee45a2b3239, .expected = 0x000059da203b7828 },
        .{ .kind = .inv_ln10, .input = 0x001f9ef8dcdc6e37, .expected = 0x000dbb9a0a83335c },
        .{ .kind = .pi, .input = 0x00180a39da4f7884, .expected = 0x0032e1848dbb0ff3 },
        .{ .kind = .pi, .input = 0x00111a2a5d7bd479, .expected = 0x002add34cdfd23bc },
        .{ .kind = .inv_pi, .input = 0x001508ecb38f52f9, .expected = 0x0006b213223e5e69 },
        .{ .kind = .inv_pi, .input = 0x0015cba89af1f855, .expected = 0x0006f00f7cd3a40b },
        .{ .kind = .two_inv_sqrt_pi, .input = 0x000c2187f7d0ae75, .expected = 0x000db03774b157d5 },
        .{ .kind = .two_inv_sqrt_pi, .input = 0x001e32121e380334, .expected = 0x00210939b4818ad1 },
    };
    inline for (cases) |case| {
        try std.testing.expectEqual(case.expected, @as(u64, @bitCast(scale(case.kind, @bitCast(case.input)))));
    }
}

test "tiny scaling preserves signs and normal-binade transitions" {
    const cases = [_]Case{
        .{ .kind = .ln2, .input = 0x8000000000000000, .expected = 0x8000000000000000 },
        .{ .kind = .inv_pi, .input = 0x8000000000000001, .expected = 0x8000000000000000 },
        .{ .kind = .pi, .input = 0x0000000000000001, .expected = 0x0000000000000003 },
        .{ .kind = .pi, .input = 0x000a2f9836e4e441, .expected = 0x001fffffffffffff },
        .{ .kind = .pi, .input = 0x000a2f9836e4e442, .expected = 0x0020000000000001 },
        .{ .kind = .inv_ln2, .input = 0x00162e42fefa39ef, .expected = 0x0020000000000000 },
        .{ .kind = .inv_ln2, .input = 0x00162e42fefa39f0, .expected = 0x0020000000000001 },
        .{ .kind = .ln10, .input = 0x001bcb7b1526e50e, .expected = 0x0030000000000000 },
        .{ .kind = .ln10, .input = 0x001bcb7b1526e50f, .expected = 0x0030000000000000 },
    };
    inline for (cases) |case| {
        try std.testing.expectEqual(case.expected, @as(u64, @bitCast(scale(case.kind, @bitCast(case.input)))));
    }
}
