const std = @import("std");
const common = @import("scalar_c23_common.zig");

const IntegralRound = enum {
    nearest_even,
    nearest_from_zero,
    toward_zero,
    toward_positive,
    toward_negative,
};

const HConfig = struct {
    const Float = f16;
    const Bits = u16;
    const sign_mask: Bits = 0x8000;
    const exponent_mask: Bits = 0x7c00;
    const fraction_mask: Bits = 0x03ff;
    const quiet_nan_bit: Bits = 0x0200;
    const fraction_bits: u7 = 10;
    const exponent_bias: i32 = 15;
    const maximum_exponent: u32 = 0x1f;
};

const FConfig = struct {
    const Float = f32;
    const Bits = u32;
    const sign_mask: Bits = 0x8000_0000;
    const exponent_mask: Bits = 0x7f80_0000;
    const fraction_mask: Bits = 0x007f_ffff;
    const quiet_nan_bit: Bits = 0x0040_0000;
    const fraction_bits: u7 = 23;
    const exponent_bias: i32 = 127;
    const maximum_exponent: u32 = 0xff;
};

const DConfig = struct {
    const Float = f64;
    const Bits = u64;
    const sign_mask: Bits = 0x8000_0000_0000_0000;
    const exponent_mask: Bits = 0x7ff0_0000_0000_0000;
    const fraction_mask: Bits = 0x000f_ffff_ffff_ffff;
    const quiet_nan_bit: Bits = 0x0008_0000_0000_0000;
    const fraction_bits: u7 = 52;
    const exponent_bias: i32 = 1023;
    const maximum_exponent: u32 = 0x7ff;
};

fn roundBits(comptime Config: type, comptime mode: IntegralRound, bits: Config.Bits) Config.Bits {
    @setFloatMode(.strict);

    const magnitude = bits & ~Config.sign_mask;
    if (mode != .nearest_even) {
        if (magnitude >= Config.exponent_mask) {
            if (magnitude > Config.exponent_mask) return bits | Config.quiet_nan_bit;
            return bits;
        }
        if (magnitude == 0) return bits;

        // f16 widens exactly to f32. Nonintegral f16 is below 2^10, so its rounded
        // integer fits exactly in f16.
        const NativeFloat = if (Config == HConfig) f32 else Config.Float;
        const value: NativeFloat = @floatCast(@as(Config.Float, @bitCast(bits)));
        const rounded: NativeFloat = switch (mode) {
            .nearest_even => unreachable,
            .nearest_from_zero => @round(value),
            .toward_zero => @trunc(value),
            .toward_positive => @ceil(value),
            .toward_negative => @floor(value),
        };
        return @bitCast(@as(Config.Float, @floatCast(rounded)));
    }

    // @round breaks ties away from zero; this path keeps nearest-even ties.
    const exponent_field: u32 = @intCast(
        (magnitude & Config.exponent_mask) >> Config.fraction_bits,
    );
    const fraction = magnitude & Config.fraction_mask;

    if (exponent_field == Config.maximum_exponent) {
        if (fraction == 0) return bits;
        return bits | Config.quiet_nan_bit;
    }
    if (magnitude == 0) return bits;

    const unbiased_exponent: i32 = if (exponent_field == 0)
        1 - Config.exponent_bias
    else
        @as(i32, @intCast(exponent_field)) - Config.exponent_bias;
    const fraction_bits_i32: i32 = @intCast(Config.fraction_bits);

    if (unbiased_exponent < 0) {
        const increment = unbiased_exponent == -1 and fraction != 0;
        const one_magnitude: Config.Bits = @as(Config.Bits, @intCast(Config.exponent_bias)) <<
            Config.fraction_bits;
        return (bits & Config.sign_mask) | if (increment) one_magnitude else 0;
    }

    if (unbiased_exponent >= fraction_bits_i32) return bits;

    const shift: std.math.Log2Int(Config.Bits) =
        @intCast(fraction_bits_i32 - unbiased_exponent);
    const discarded_mask: Config.Bits = (@as(Config.Bits, 1) << shift) - 1;
    const discarded = magnitude & discarded_mask;
    if (discarded == 0) return bits;

    const truncated = magnitude & ~discarded_mask;
    const unit: Config.Bits = @as(Config.Bits, 1) << shift;
    const half: Config.Bits = unit >> 1;
    const increment = discarded > half or
        (discarded == half and ((truncated >> shift) & 1) != 0);

    const rounded_magnitude = truncated + if (increment) unit else 0;
    return (bits & Config.sign_mask) | rounded_magnitude;
}

/// Rounds `value` toward positive infinity. The direction is fixed in code and
/// does not depend on the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn ceilh(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, .toward_positive, bits));
}

/// Rounds `value` toward positive infinity. The direction is fixed in code and
/// does not depend on the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn ceilf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, .toward_positive, bits));
}

/// Rounds `value` toward positive infinity. The direction is fixed in code and
/// does not depend on the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn ceild(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, .toward_positive, bits));
}

/// Rounds `value` toward negative infinity. The direction is fixed in code and
/// does not depend on the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn floorh(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, .toward_negative, bits));
}

/// Rounds `value` toward negative infinity. The direction is fixed in code and
/// does not depend on the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn floorf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, .toward_negative, bits));
}

/// Rounds `value` toward negative infinity. The direction is fixed in code and
/// does not depend on the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn floord(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, .toward_negative, bits));
}

/// Rounds `value` toward zero. The direction is fixed in code and does not
/// depend on the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn trunch(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, .toward_zero, bits));
}

/// Rounds `value` toward zero. The direction is fixed in code and does not
/// depend on the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn truncf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, .toward_zero, bits));
}

/// Rounds `value` toward zero. The direction is fixed in code and does not
/// depend on the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn truncd(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, .toward_zero, bits));
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour. The direction is fixed in code and does not depend on the ambient
/// rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn nearbyinth(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, .nearest_even, bits));
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour. The direction is fixed in code and does not depend on the ambient
/// rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn nearbyintf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, .nearest_even, bits));
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour. The direction is fixed in code and does not depend on the ambient
/// rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn nearbyintd(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, .nearest_even, bits));
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour, by calling `nearbyinth`. Unlike the C function of the same name
/// it never consults the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn rinth(value: f16) f16 {
    return nearbyinth(value);
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour, by calling `nearbyintf`. Unlike the C function of the same name
/// it never consults the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn rintf(value: f32) f32 {
    return nearbyintf(value);
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour, by calling `nearbyintd`. Unlike the C function of the same name
/// it never consults the ambient rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn rintd(value: f64) f64 {
    return nearbyintd(value);
}

/// Rounds `value` to the nearest integer, sending a halfway case away from
/// zero. The direction is fixed in code and does not depend on the ambient
/// rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn roundh(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, .nearest_from_zero, bits));
}

/// Rounds `value` to the nearest integer, sending a halfway case away from
/// zero. The direction is fixed in code and does not depend on the ambient
/// rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn roundf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, .nearest_from_zero, bits));
}

/// Rounds `value` to the nearest integer, sending a halfway case away from
/// zero. The direction is fixed in code and does not depend on the ambient
/// rounding mode.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn roundd(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, .nearest_from_zero, bits));
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour, by calling `nearbyinth`. The ambient rounding mode is not read.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn roundevenh(value: f16) f16 {
    return nearbyinth(value);
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour, by calling `nearbyintf`. The ambient rounding mode is not read.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn roundevenf(value: f32) f32 {
    return nearbyintf(value);
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour, by calling `nearbyintd`. The ambient rounding mode is not read.
///
/// An infinity and a zero of either sign are returned unchanged; a NaN is
/// returned with its quiet bit set, keeping its sign and payload.
pub fn roundevend(value: f64) f64 {
    return nearbyintd(value);
}
