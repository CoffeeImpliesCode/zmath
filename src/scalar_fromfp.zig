const std = @import("std");
const common = @import("scalar_c23_common.zig");
const c23 = @import("scalar_c23_types.zig");

fn BinaryFormat(
    comptime Ieee: type,
    comptime fraction_bits: comptime_int,
    comptime exponent_bits: comptime_int,
    comptime exponent_bias: i32,
) type {
    return struct {
        const Float = Ieee.Float;
        const Bits = Ieee.Bits;
        const Shift = std.math.Log2Int(Bits);

        const fraction_width: i32 = fraction_bits;
        const bias: i32 = exponent_bias;
        const sign_mask: Bits = @as(Bits, 1) << (fraction_bits + exponent_bits);
        const fraction_mask: Bits = (@as(Bits, 1) << fraction_bits) - 1;
        const exponent_mask: Bits =
            ((@as(Bits, 1) << exponent_bits) - 1) << fraction_bits;
        const infinity_bits: Bits = exponent_mask;
        const canonical_nan_bits: Bits =
            infinity_bits | (@as(Bits, 1) << (fraction_bits - 1));

        fn toBits(value: Float) Bits {
            return Ieee.toBits(value);
        }

        fn fromBits(bits: Bits) Float {
            return Ieee.fromBits(bits);
        }

        fn invalidResult() Float {
            return fromBits(canonical_nan_bits);
        }

        fn encodeInteger(magnitude: Bits, negative: bool) Float {
            if (magnitude == 0) {
                return fromBits(if (negative) sign_mask else 0);
            }

            const leading_zeros: u32 = @intCast(@clz(magnitude));
            const top_bit: u32 =
                @as(u32, @intCast(@bitSizeOf(Bits) - 1)) - leading_zeros;
            std.debug.assert(top_bit <= fraction_bits);

            const left_shift: u32 = @as(u32, fraction_bits) - top_bit;
            const shift: Shift = @intCast(left_shift);
            const significand = magnitude << shift;
            const exponent: i32 = @as(i32, @intCast(top_bit)) + bias;
            const exponent_field: Bits = @as(Bits, @intCast(exponent)) << fraction_bits;
            const sign_field: Bits = if (negative) sign_mask else 0;

            return fromBits(sign_field | exponent_field | (significand & fraction_mask));
        }

        fn roundedBits(bits: Bits, round: c23.FromFpRound) Bits {
            const sign = bits & sign_mask != 0;
            const exponent_field: Bits = (bits & exponent_mask) >> fraction_bits;
            const fraction = bits & fraction_mask;

            if (exponent_field == 0) {
                if (fraction == 0) return bits;
            }

            // SAFETY: Both `if` arms below assign `exponent` before its first read.
            var exponent: i32 = undefined;
            // SAFETY: Both `if` arms below assign `significand` before its first read.
            var significand: Bits = undefined;
            if (exponent_field == 0) {
                exponent = 1 - bias;
                significand = fraction;
            } else {
                exponent = @as(i32, @intCast(exponent_field)) - bias;
                significand = (@as(Bits, 1) << fraction_bits) | fraction;
            }

            if (exponent >= fraction_width) return bits;

            var magnitude: Bits = 0;
            if (exponent < 0) {
                const half_significand: Bits = @as(Bits, 1) << fraction_bits;
                const exact_half = exponent == -1 and significand == half_significand;
                const above_half = exponent == -1 and significand > half_significand;
                const increment = switch (round) {
                    .nearest_even => above_half,
                    .nearest_from_zero => exact_half or above_half,
                    .toward_zero => false,
                    .toward_positive => !sign,
                    .toward_negative => sign,
                };
                if (increment) magnitude = 1;
            } else {
                const fractional_bit_count: u32 =
                    @intCast(fraction_width - exponent);
                const shift: Shift = @intCast(fractional_bit_count);
                const fractional_mask: Bits = (@as(Bits, 1) << shift) - 1;
                const remainder = significand & fractional_mask;
                magnitude = significand >> shift;

                const half_shift: Shift = @intCast(fractional_bit_count - 1);
                const half: Bits = @as(Bits, 1) << half_shift;
                const exact_half = remainder == half;
                const above_half = remainder > half;
                const increment = switch (round) {
                    .nearest_even => above_half or (exact_half and (magnitude & 1 != 0)),
                    .nearest_from_zero => exact_half or above_half,
                    .toward_zero => false,
                    .toward_positive => !sign and remainder != 0,
                    .toward_negative => sign and remainder != 0,
                };
                if (increment) magnitude += 1;
            }

            return toBits(encodeInteger(magnitude, sign));
        }

        fn belowPowerOfTwo(bits: Bits, power_exponent: i64, inclusive: bool) bool {
            const magnitude = bits & ~sign_mask;
            if (magnitude == 0) return true;

            const exponent_field: Bits = (magnitude & exponent_mask) >> fraction_bits;
            const exponent: i64 =
                @as(i64, @intCast(exponent_field)) - @as(i64, bias);
            if (exponent < power_exponent) return true;
            if (exponent > power_exponent) return false;
            if (!inclusive) return false;
            return magnitude & fraction_mask == 0;
        }

        fn withinWidth(bits: Bits, round_width: u32, unsigned: bool) bool {
            if (round_width == 0) return false;

            const magnitude = bits & ~sign_mask;
            if (magnitude == 0) return true;
            const negative = bits & sign_mask != 0;

            const width_i64: i64 = @intCast(round_width);
            if (unsigned) {
                if (negative) return false;
                return belowPowerOfTwo(bits, width_i64, false);
            }

            const positive_limit_exponent: i64 = width_i64 - 1;
            if (negative) {
                return belowPowerOfTwo(bits, positive_limit_exponent, true);
            }
            return belowPowerOfTwo(bits, positive_limit_exponent, false);
        }

        fn convert(
            value: Float,
            round: c23.FromFpRound,
            round_width: u32,
            unsigned: bool,
        ) Float {
            const bits = toBits(value);
            const magnitude = bits & ~sign_mask;
            if (round_width == 0) return invalidResult();
            if (magnitude >= infinity_bits) return invalidResult();

            const result_bits = roundedBits(bits, round);
            if (!withinWidth(result_bits, round_width, unsigned)) return invalidResult();
            return fromBits(result_bits);
        }
    };
}

const H = BinaryFormat(common.H, 10, 5, 15);
const F = BinaryFormat(common.F, 23, 8, 127);
const D = BinaryFormat(common.D, 52, 11, 1023);

/// Converts a binary16 value to an integral value using the selected C23 direction.
pub fn fromfph(value: f16, round: c23.FromFpRound, width: u32) f16 {
    return H.convert(value, round, width, false);
}

/// Converts a binary32 value to an integral value using the selected C23 direction.
pub fn fromfpf(value: f32, round: c23.FromFpRound, width: u32) f32 {
    return F.convert(value, round, width, false);
}

/// Converts a binary64 value to an integral value using the selected C23 direction.
pub fn fromfpd(value: f64, round: c23.FromFpRound, width: u32) f64 {
    return D.convert(value, round, width, false);
}

/// Converts a binary16 value to an unsigned integral value using the selected direction.
pub fn ufromfph(value: f16, round: c23.FromFpRound, width: u32) f16 {
    return H.convert(value, round, width, true);
}

/// Converts a binary32 value to an unsigned integral value using the selected direction.
pub fn ufromfpf(value: f32, round: c23.FromFpRound, width: u32) f32 {
    return F.convert(value, round, width, true);
}

/// Converts a binary64 value to an unsigned integral value using the selected direction.
pub fn ufromfpd(value: f64, round: c23.FromFpRound, width: u32) f64 {
    return D.convert(value, round, width, true);
}

/// Converts a binary16 value with the C23 inexact-reporting spelling.
pub fn fromfpxh(value: f16, round: c23.FromFpRound, width: u32) f16 {
    return H.convert(value, round, width, false);
}

/// Converts a binary32 value with the C23 inexact-reporting spelling.
pub fn fromfpxf(value: f32, round: c23.FromFpRound, width: u32) f32 {
    return F.convert(value, round, width, false);
}

/// Converts a binary64 value with the C23 inexact-reporting spelling.
pub fn fromfpxd(value: f64, round: c23.FromFpRound, width: u32) f64 {
    return D.convert(value, round, width, false);
}

/// Converts a binary16 value to an unsigned value with the C23 x spelling.
pub fn ufromfpxh(value: f16, round: c23.FromFpRound, width: u32) f16 {
    return H.convert(value, round, width, true);
}

/// Converts a binary32 value to an unsigned value with the C23 x spelling.
pub fn ufromfpxf(value: f32, round: c23.FromFpRound, width: u32) f32 {
    return F.convert(value, round, width, true);
}

/// Converts a binary64 value to an unsigned value with the C23 x spelling.
pub fn ufromfpxd(value: f64, round: c23.FromFpRound, width: u32) f64 {
    return D.convert(value, round, width, true);
}
