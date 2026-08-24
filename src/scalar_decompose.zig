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

        const fraction_width: u32 = fraction_bits;
        const fraction_width_i32: i32 = fraction_bits;
        const bias: i32 = exponent_bias;
        const sign_mask: Bits = @as(Bits, 1) << (fraction_bits + exponent_bits);
        const fraction_mask: Bits = (@as(Bits, 1) << fraction_bits) - 1;
        const exponent_mask: Bits =
            ((@as(Bits, 1) << exponent_bits) - 1) << fraction_bits;
        const half_exponent_bits: Bits =
            @as(Bits, @intCast(exponent_bias - 1)) << fraction_bits;
        const negative_infinity_bits: Bits = sign_mask | exponent_mask;
        const min_subnormal_exponent: i32 = 1 - exponent_bias - fraction_bits;

        fn toBits(value: Float) Bits {
            return Ieee.toBits(value);
        }

        fn fromBits(bits: Bits) Float {
            return Ieee.fromBits(bits);
        }

        fn classifyBits(bits: Bits) common.Class {
            return Ieee.classifyBits(bits);
        }

        fn quietNaNBits(bits: Bits) Bits {
            return Ieee.quietNaNBits(bits);
        }

        fn signedZero(bits: Bits) Float {
            return fromBits(bits & sign_mask);
        }

        fn highestFractionBit(fraction: Bits) u32 {
            std.debug.assert(fraction != 0);
            std.debug.assert(fraction & ~fraction_mask == 0);

            const leading_zeros: u32 = @intCast(@clz(fraction));
            const top_bit: u32 = @bitSizeOf(Bits) - 1;
            return top_bit - leading_zeros;
        }

        fn normalizedFractionBits(fraction: Bits) Bits {
            const shift: Shift = @intCast(fraction_width - highestFractionBit(fraction));
            return (fraction << shift) & fraction_mask;
        }

        fn normalExponent(bits: Bits) i32 {
            const exponent_field = (bits & exponent_mask) >> fraction_bits;
            return @as(i32, @intCast(exponent_field)) - bias;
        }

        fn subnormalExponent(fraction: Bits) i32 {
            return min_subnormal_exponent + @as(i32, @intCast(highestFractionBit(fraction)));
        }
    };
}

const H = BinaryFormat(common.H, 10, 5, 15);
const F = BinaryFormat(common.F, 23, 8, 127);
const D = BinaryFormat(common.D, 52, 11, 1023);

fn frexpValue(comptime Format: type, value: Format.Float) c23.FrexpResult(Format.Float) {
    const bits = Format.toBits(value);

    return switch (Format.classifyBits(bits)) {
        .zero => .{ .fraction = value, .exponent = 0 },
        .subnormal => blk: {
            const fraction = bits & Format.fraction_mask;
            const result_bits = (bits & Format.sign_mask) |
                Format.half_exponent_bits |
                Format.normalizedFractionBits(fraction);
            break :blk .{
                .fraction = Format.fromBits(result_bits),
                .exponent = Format.subnormalExponent(fraction) + 1,
            };
        },
        .normal => .{
            .fraction = Format.fromBits(
                (bits & Format.sign_mask) |
                    Format.half_exponent_bits |
                    (bits & Format.fraction_mask),
            ),
            .exponent = Format.normalExponent(bits) + 1,
        },
        .infinity => .{ .fraction = value, .exponent = 0 },
        .quiet_nan, .signaling_nan => .{
            .fraction = Format.fromBits(Format.quietNaNBits(bits)),
            .exponent = 0,
        },
    };
}

fn modfValue(comptime Format: type, value: Format.Float) c23.ModfResult(Format.Float) {
    @setFloatMode(.strict);
    const bits = Format.toBits(value);

    return switch (Format.classifyBits(bits)) {
        .zero => .{ .fractional = value, .integral = value },
        .subnormal => .{
            .fractional = value,
            .integral = Format.signedZero(bits),
        },
        .normal => blk: {
            const exponent = Format.normalExponent(bits);
            if (exponent < 0) {
                break :blk .{
                    .fractional = value,
                    .integral = Format.signedZero(bits),
                };
            }
            if (exponent >= Format.fraction_width_i32) {
                break :blk .{
                    .fractional = Format.signedZero(bits),
                    .integral = value,
                };
            }

            const fractional_bit_count: u32 =
                @intCast(Format.fraction_width_i32 - exponent);
            const shift: Format.Shift = @intCast(fractional_bit_count);
            const fractional_mask = (@as(Format.Bits, 1) << shift) - 1;
            if (bits & fractional_mask == 0) {
                break :blk .{
                    .fractional = Format.signedZero(bits),
                    .integral = value,
                };
            }

            const integral = Format.fromBits(bits & ~fractional_mask);
            break :blk .{
                .fractional = value - integral,
                .integral = integral,
            };
        },
        .infinity => .{
            .fractional = Format.signedZero(bits),
            .integral = value,
        },
        .quiet_nan, .signaling_nan => blk: {
            const quiet_nan = Format.fromBits(Format.quietNaNBits(bits));
            break :blk .{ .fractional = quiet_nan, .integral = quiet_nan };
        },
    };
}

fn integerLogb(comptime Format: type, comptime Int: type, value: Format.Float) Int {
    const bits = Format.toBits(value);

    return switch (Format.classifyBits(bits)) {
        .zero => std.math.minInt(Int),
        .subnormal => @intCast(Format.subnormalExponent(bits & Format.fraction_mask)),
        .normal => @intCast(Format.normalExponent(bits)),
        .infinity, .quiet_nan, .signaling_nan => std.math.maxInt(Int),
    };
}

fn logbValue(comptime Format: type, value: Format.Float) Format.Float {
    const bits = Format.toBits(value);

    return switch (Format.classifyBits(bits)) {
        .zero => Format.fromBits(Format.negative_infinity_bits),
        .subnormal => @floatFromInt(
            Format.subnormalExponent(bits & Format.fraction_mask),
        ),
        .normal => @floatFromInt(Format.normalExponent(bits)),
        .infinity => Format.fromBits(Format.exponent_mask),
        .quiet_nan, .signaling_nan => Format.fromBits(Format.quietNaNBits(bits)),
    };
}

/// Decomposes a binary16 value into a signed fraction and a power of two.
pub fn frexph(value: f16) c23.FrexpResult(f16) {
    return frexpValue(H, value);
}

/// Decomposes a binary32 value into a signed fraction and a power of two.
pub fn frexpf(value: f32) c23.FrexpResult(f32) {
    return frexpValue(F, value);
}

/// Decomposes a binary64 value into a signed fraction and a power of two.
pub fn frexpd(value: f64) c23.FrexpResult(f64) {
    return frexpValue(D, value);
}

/// Splits a binary16 value into signed fractional and integral parts.
pub fn modfh(value: f16) c23.ModfResult(f16) {
    return modfValue(H, value);
}

/// Splits a binary32 value into signed fractional and integral parts.
pub fn modff(value: f32) c23.ModfResult(f32) {
    return modfValue(F, value);
}

/// Splits a binary64 value into signed fractional and integral parts.
pub fn modfd(value: f64) c23.ModfResult(f64) {
    return modfValue(D, value);
}

/// Returns the binary exponent of a binary16 value as an i32.
pub fn ilogbh(value: f16) i32 {
    return integerLogb(H, i32, value);
}

/// Returns the binary exponent of a binary32 value as an i32.
pub fn ilogbf(value: f32) i32 {
    return integerLogb(F, i32, value);
}

/// Returns the binary exponent of a binary64 value as an i32.
pub fn ilogbd(value: f64) i32 {
    return integerLogb(D, i32, value);
}

/// Returns the binary exponent of a binary16 value as an i64.
pub fn llogbh(value: f16) i64 {
    return integerLogb(H, i64, value);
}

/// Returns the binary exponent of a binary32 value as an i64.
pub fn llogbf(value: f32) i64 {
    return integerLogb(F, i64, value);
}

/// Returns the binary exponent of a binary64 value as an i64.
pub fn llogbd(value: f64) i64 {
    return integerLogb(D, i64, value);
}

/// Returns the binary exponent of a binary16 value in binary16 format.
pub fn logbh(value: f16) f16 {
    return logbValue(H, value);
}

/// Returns the binary exponent of a binary32 value in binary32 format.
pub fn logbf(value: f32) f32 {
    return logbValue(F, value);
}

/// Returns the binary exponent of a binary64 value in binary64 format.
pub fn logbd(value: f64) f64 {
    return logbValue(D, value);
}
