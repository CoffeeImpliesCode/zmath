const common = @import("scalar_c23_common.zig");

fn signMask(comptime Format: type) Format.Bits {
    return @as(Format.Bits, 1) << (@bitSizeOf(Format.Bits) - 1);
}

fn positiveZero(comptime Format: type) Format.Float {
    return Format.fromBits(0);
}

fn firstNaN(comptime Format: type, first: Format.Float, second: Format.Float) ?Format.Float {
    return Format.firstNaN2(first, second);
}

fn maximumNonNaN(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    const first_key = Format.finiteOrderKey(Format.toBits(first));
    const second_key = Format.finiteOrderKey(Format.toBits(second));
    if (first_key >= second_key) return first;
    return second;
}

fn minimumNonNaN(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    const first_key = Format.finiteOrderKey(Format.toBits(first));
    const second_key = Format.finiteOrderKey(Format.toBits(second));
    if (first_key <= second_key) return first;
    return second;
}

fn maximumLegacy(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    const first_bits = Format.toBits(first);
    const second_bits = Format.toBits(second);
    const first_nan = Format.isNaNBits(first_bits);
    const second_nan = Format.isNaNBits(second_bits);

    if (first_nan) {
        if (second_nan) return Format.fromBits(Format.quietNaNBits(first_bits));
        return second;
    }
    if (second_nan) return first;
    return maximumNonNaN(Format, first, second);
}

fn minimumLegacy(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    const first_bits = Format.toBits(first);
    const second_bits = Format.toBits(second);
    const first_nan = Format.isNaNBits(first_bits);
    const second_nan = Format.isNaNBits(second_bits);

    if (first_nan) {
        if (second_nan) return Format.fromBits(Format.quietNaNBits(first_bits));
        return second;
    }
    if (second_nan) return first;
    return minimumNonNaN(Format, first, second);
}

fn maximumNumber(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    const first_bits = Format.toBits(first);
    const second_bits = Format.toBits(second);
    const first_nan = Format.isNaNBits(first_bits);
    const second_nan = Format.isNaNBits(second_bits);

    if (first_nan) {
        if (second_nan) return Format.fromBits(Format.quietNaNBits(first_bits));
        return second;
    }
    if (second_nan) return first;
    return maximumNonNaN(Format, first, second);
}

fn minimumNumber(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    const first_bits = Format.toBits(first);
    const second_bits = Format.toBits(second);
    const first_nan = Format.isNaNBits(first_bits);
    const second_nan = Format.isNaNBits(second_bits);

    if (first_nan) {
        if (second_nan) return Format.fromBits(Format.quietNaNBits(first_bits));
        return second;
    }
    if (second_nan) return first;
    return minimumNonNaN(Format, first, second);
}

fn maximumMagnitudeNonNaN(
    comptime Format: type,
    first: Format.Float,
    second: Format.Float,
) Format.Float {
    const magnitude_mask = ~signMask(Format);
    const first_magnitude = Format.toBits(first) & magnitude_mask;
    const second_magnitude = Format.toBits(second) & magnitude_mask;

    if (first_magnitude > second_magnitude) return first;
    if (first_magnitude < second_magnitude) return second;
    return maximumNonNaN(Format, first, second);
}

fn minimumMagnitudeNonNaN(
    comptime Format: type,
    first: Format.Float,
    second: Format.Float,
) Format.Float {
    const magnitude_mask = ~signMask(Format);
    const first_magnitude = Format.toBits(first) & magnitude_mask;
    const second_magnitude = Format.toBits(second) & magnitude_mask;

    if (first_magnitude < second_magnitude) return first;
    if (first_magnitude > second_magnitude) return second;
    return minimumNonNaN(Format, first, second);
}

fn maximumMagnitude(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    if (firstNaN(Format, first, second)) |nan| return nan;
    return maximumMagnitudeNonNaN(Format, first, second);
}

fn minimumMagnitude(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    if (firstNaN(Format, first, second)) |nan| return nan;
    return minimumMagnitudeNonNaN(Format, first, second);
}

fn maximumMagnitudeNumber(
    comptime Format: type,
    first: Format.Float,
    second: Format.Float,
) Format.Float {
    const first_bits = Format.toBits(first);
    const second_bits = Format.toBits(second);
    const first_nan = Format.isNaNBits(first_bits);
    const second_nan = Format.isNaNBits(second_bits);

    if (first_nan) {
        if (second_nan) return Format.fromBits(Format.quietNaNBits(first_bits));
        return second;
    }
    if (second_nan) return first;
    return maximumMagnitudeNonNaN(Format, first, second);
}

fn minimumMagnitudeNumber(
    comptime Format: type,
    first: Format.Float,
    second: Format.Float,
) Format.Float {
    const first_bits = Format.toBits(first);
    const second_bits = Format.toBits(second);
    const first_nan = Format.isNaNBits(first_bits);
    const second_nan = Format.isNaNBits(second_bits);

    if (first_nan) {
        if (second_nan) return Format.fromBits(Format.quietNaNBits(first_bits));
        return second;
    }
    if (second_nan) return first;
    return minimumMagnitudeNonNaN(Format, first, second);
}

fn fdim(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    @setFloatMode(.strict);
    if (firstNaN(Format, first, second)) |nan| return nan;
    if (first > second) return first - second;
    return positiveZero(Format);
}

fn fmax(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    return maximumLegacy(Format, first, second);
}

fn fmin(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    return minimumLegacy(Format, first, second);
}

fn fmaximum(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    if (firstNaN(Format, first, second)) |nan| return nan;
    return maximumNonNaN(Format, first, second);
}

fn fminimum(comptime Format: type, first: Format.Float, second: Format.Float) Format.Float {
    if (firstNaN(Format, first, second)) |nan| return nan;
    return minimumNonNaN(Format, first, second);
}

pub fn fdimh(first: f16, second: f16) f16 {
    return fdim(common.H, first, second);
}

pub fn fdimf(first: f32, second: f32) f32 {
    return fdim(common.F, first, second);
}

pub fn fdimd(first: f64, second: f64) f64 {
    return fdim(common.D, first, second);
}

pub fn fmaxh(first: f16, second: f16) f16 {
    return fmax(common.H, first, second);
}

pub fn fmaxf(first: f32, second: f32) f32 {
    return fmax(common.F, first, second);
}

pub fn fmaxd(first: f64, second: f64) f64 {
    return fmax(common.D, first, second);
}

pub fn fminh(first: f16, second: f16) f16 {
    return fmin(common.H, first, second);
}

pub fn fminf(first: f32, second: f32) f32 {
    return fmin(common.F, first, second);
}

pub fn fmind(first: f64, second: f64) f64 {
    return fmin(common.D, first, second);
}

pub fn fmaximumh(first: f16, second: f16) f16 {
    return fmaximum(common.H, first, second);
}

pub fn fmaximumf(first: f32, second: f32) f32 {
    return fmaximum(common.F, first, second);
}

pub fn fmaximumd(first: f64, second: f64) f64 {
    return fmaximum(common.D, first, second);
}

pub fn fminimumh(first: f16, second: f16) f16 {
    return fminimum(common.H, first, second);
}

pub fn fminimumf(first: f32, second: f32) f32 {
    return fminimum(common.F, first, second);
}

pub fn fminimumd(first: f64, second: f64) f64 {
    return fminimum(common.D, first, second);
}

pub fn fmaximumMagh(first: f16, second: f16) f16 {
    return maximumMagnitude(common.H, first, second);
}

pub fn fmaximumMagf(first: f32, second: f32) f32 {
    return maximumMagnitude(common.F, first, second);
}

pub fn fmaximumMagd(first: f64, second: f64) f64 {
    return maximumMagnitude(common.D, first, second);
}

pub fn fminimumMagh(first: f16, second: f16) f16 {
    return minimumMagnitude(common.H, first, second);
}

pub fn fminimumMagf(first: f32, second: f32) f32 {
    return minimumMagnitude(common.F, first, second);
}

pub fn fminimumMagd(first: f64, second: f64) f64 {
    return minimumMagnitude(common.D, first, second);
}

pub fn fmaximumNumh(first: f16, second: f16) f16 {
    return maximumNumber(common.H, first, second);
}

pub fn fmaximumNumf(first: f32, second: f32) f32 {
    return maximumNumber(common.F, first, second);
}

pub fn fmaximumNumd(first: f64, second: f64) f64 {
    return maximumNumber(common.D, first, second);
}

pub fn fminimumNumh(first: f16, second: f16) f16 {
    return minimumNumber(common.H, first, second);
}

pub fn fminimumNumf(first: f32, second: f32) f32 {
    return minimumNumber(common.F, first, second);
}

pub fn fminimumNumd(first: f64, second: f64) f64 {
    return minimumNumber(common.D, first, second);
}

pub fn fmaximumMagNumh(first: f16, second: f16) f16 {
    return maximumMagnitudeNumber(common.H, first, second);
}

pub fn fmaximumMagNumf(first: f32, second: f32) f32 {
    return maximumMagnitudeNumber(common.F, first, second);
}

pub fn fmaximumMagNumd(first: f64, second: f64) f64 {
    return maximumMagnitudeNumber(common.D, first, second);
}

pub fn fminimumMagNumh(first: f16, second: f16) f16 {
    return minimumMagnitudeNumber(common.H, first, second);
}

pub fn fminimumMagNumf(first: f32, second: f32) f32 {
    return minimumMagnitudeNumber(common.F, first, second);
}

pub fn fminimumMagNumd(first: f64, second: f64) f64 {
    return minimumMagnitudeNumber(common.D, first, second);
}
