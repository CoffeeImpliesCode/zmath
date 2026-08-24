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
    const Bits = u64;
    const sign_mask: Bits = 0x8000_0000_0000_0000;
    const exponent_mask: Bits = 0x7ff0_0000_0000_0000;
    const fraction_mask: Bits = 0x000f_ffff_ffff_ffff;
    const quiet_nan_bit: Bits = 0x0008_0000_0000_0000;
    const fraction_bits: u7 = 52;
    const exponent_bias: i32 = 1023;
    const maximum_exponent: u32 = 0x7ff;
};

fn roundBits(comptime Config: type, bits: Config.Bits, mode: IntegralRound) Config.Bits {
    const magnitude = bits & ~Config.sign_mask;
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
        const increment = switch (mode) {
            .nearest_even => unbiased_exponent == -1 and fraction != 0,
            .nearest_from_zero => unbiased_exponent == -1,
            .toward_zero => false,
            .toward_positive => bits & Config.sign_mask == 0,
            .toward_negative => bits & Config.sign_mask != 0,
        };
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
    const truncated_integer_lsb = (truncated >> shift) & 1;
    const increment = switch (mode) {
        .nearest_even => discarded > half or
            (discarded == half and truncated_integer_lsb != 0),
        .nearest_from_zero => discarded >= half,
        .toward_zero => false,
        .toward_positive => bits & Config.sign_mask == 0,
        .toward_negative => bits & Config.sign_mask != 0,
    };

    const rounded_magnitude = truncated + if (increment) unit else 0;
    return (bits & Config.sign_mask) | rounded_magnitude;
}

pub fn ceilh(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, bits, .toward_positive));
}

pub fn ceilf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, bits, .toward_positive));
}

pub fn ceild(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, bits, .toward_positive));
}

pub fn floorh(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, bits, .toward_negative));
}

pub fn floorf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, bits, .toward_negative));
}

pub fn floord(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, bits, .toward_negative));
}

pub fn trunch(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, bits, .toward_zero));
}

pub fn truncf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, bits, .toward_zero));
}

pub fn truncd(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, bits, .toward_zero));
}

pub fn nearbyinth(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, bits, .nearest_even));
}

pub fn nearbyintf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, bits, .nearest_even));
}

pub fn nearbyintd(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, bits, .nearest_even));
}

pub fn rinth(value: f16) f16 {
    return nearbyinth(value);
}

pub fn rintf(value: f32) f32 {
    return nearbyintf(value);
}

pub fn rintd(value: f64) f64 {
    return nearbyintd(value);
}

pub fn roundh(value: f16) f16 {
    const bits: HConfig.Bits = common.H.toBits(value);
    return common.H.fromBits(roundBits(HConfig, bits, .nearest_from_zero));
}

pub fn roundf(value: f32) f32 {
    const bits: FConfig.Bits = common.F.toBits(value);
    return common.F.fromBits(roundBits(FConfig, bits, .nearest_from_zero));
}

pub fn roundd(value: f64) f64 {
    const bits: DConfig.Bits = common.D.toBits(value);
    return common.D.fromBits(roundBits(DConfig, bits, .nearest_from_zero));
}

pub fn roundevenh(value: f16) f16 {
    return nearbyinth(value);
}

pub fn roundevenf(value: f32) f32 {
    return nearbyintf(value);
}

pub fn roundevend(value: f64) f64 {
    return nearbyintd(value);
}
