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

/// Returns `first - second` when `first > second` and positive zero otherwise,
/// so an equal or negative difference is `+0` and never `-0`.
///
/// A NaN operand is propagated as the first NaN encountered, with its quiet bit
/// set and its sign and payload kept. The returned value is the whole report:
/// there is no `errno`, no floating-point exception flag, and no rounding mode
/// is read or set.
///
/// Callers must keep IEEE round-to-nearest, ties-to-even active.
pub fn fdimh(first: f16, second: f16) f16 {
    return fdim(common.H, first, second);
}

/// The binary32 form of `fdimh`.
pub fn fdimf(first: f32, second: f32) f32 {
    return fdim(common.F, first, second);
}

/// The binary64 form of `fdimh`.
pub fn fdimd(first: f64, second: f64) f64 {
    return fdim(common.D, first, second);
}

/// The larger of two values, treating a NaN operand as missing: a NaN `first`
/// yields `second`, a NaN `second` yields `first`, and two NaNs yield `first`
/// with its quiet bit set. This is the C `fmax` policy, and it is the whole
/// difference from `fmaximumh`, which returns a NaN instead. Otherwise the
/// operand of greater value wins, so a `+0`/`-0` pair gives `+0` and an
/// equal-magnitude opposite-sign pair gives the positive operand.
pub fn fmaxh(first: f16, second: f16) f16 {
    return fmax(common.H, first, second);
}

/// The binary32 form of `fmaxh`; a NaN operand is missing rather than
/// propagated.
pub fn fmaxf(first: f32, second: f32) f32 {
    return fmax(common.F, first, second);
}

/// The binary64 form of `fmaxh`; a NaN operand is missing rather than
/// propagated.
pub fn fmaxd(first: f64, second: f64) f64 {
    return fmax(common.D, first, second);
}

/// The smaller of two values under the same missing-NaN policy as `fmaxh`, so a
/// NaN operand is skipped, two NaNs give `first` with its quiet bit set, and a
/// `+0`/`-0` pair gives `-0`.
pub fn fminh(first: f16, second: f16) f16 {
    return fmin(common.H, first, second);
}

/// The binary32 form of `fminh`; a NaN operand is missing rather than
/// propagated.
pub fn fminf(first: f32, second: f32) f32 {
    return fmin(common.F, first, second);
}

/// The binary64 form of `fminh`; a NaN operand is missing rather than
/// propagated.
pub fn fmind(first: f64, second: f64) f64 {
    return fmin(common.D, first, second);
}

/// The larger of two values, returning a NaN rather than skipping it: the
/// first NaN operand is propagated with its quiet bit set and its sign and
/// payload kept, and the ordered comparison runs only when neither operand is
/// a NaN. That is the whole difference from `fmaxh`. A `+0`/`-0` pair gives
/// `+0` and an equal-magnitude opposite-sign pair gives the positive operand.
pub fn fmaximumh(first: f16, second: f16) f16 {
    return fmaximum(common.H, first, second);
}

/// The binary32 form of `fmaximumh`; a NaN operand is propagated, not skipped.
pub fn fmaximumf(first: f32, second: f32) f32 {
    return fmaximum(common.F, first, second);
}

/// The binary64 form of `fmaximumh`; a NaN operand is propagated, not skipped.
pub fn fmaximumd(first: f64, second: f64) f64 {
    return fmaximum(common.D, first, second);
}

/// The smaller of two values, propagating a NaN as `fmaximumh` does: the first
/// NaN operand is returned with its quiet bit set and its sign and payload
/// kept. A `+0`/`-0` pair gives `-0` and an equal-magnitude opposite-sign pair
/// gives the negative operand.
pub fn fminimumh(first: f16, second: f16) f16 {
    return fminimum(common.H, first, second);
}

/// The binary32 form of `fminimumh`; a NaN operand is propagated, not skipped.
pub fn fminimumf(first: f32, second: f32) f32 {
    return fminimum(common.F, first, second);
}

/// The binary64 form of `fminimumh`; a NaN operand is propagated, not skipped.
pub fn fminimumd(first: f64, second: f64) f64 {
    return fminimum(common.D, first, second);
}

/// Selects the operand of larger magnitude, comparing the two sign-cleared
/// encodings, while propagating a NaN as `fmaximumh` does. An exact magnitude
/// tie — the `+0`/`-0` pair, or equal magnitudes of opposite sign — falls back
/// to the ordered comparison, so the positive operand wins and a tie never
/// settles on `-0` against `+0`. `scalar.zig` re-exports this function as
/// `fmaximum_magh`; the names used in these docs are the in-file ones.
pub fn fmaximumMagh(first: f16, second: f16) f16 {
    return maximumMagnitude(common.H, first, second);
}

/// The binary32 form of `fmaximumMagh`: the larger magnitude wins, an exact
/// magnitude tie falls back to the ordered comparison, and a NaN is propagated.
pub fn fmaximumMagf(first: f32, second: f32) f32 {
    return maximumMagnitude(common.F, first, second);
}

/// The binary64 form of `fmaximumMagh`: the larger magnitude wins, an exact
/// magnitude tie falls back to the ordered comparison, and a NaN is propagated.
pub fn fmaximumMagd(first: f64, second: f64) f64 {
    return maximumMagnitude(common.D, first, second);
}

/// The magnitude counterpart of `fmaximumMagh`: the smaller magnitude wins, an
/// exact magnitude tie falls back to the ordered comparison, and a NaN is
/// propagated.
pub fn fminimumMagh(first: f16, second: f16) f16 {
    return minimumMagnitude(common.H, first, second);
}

/// The binary32 form of `fminimumMagh`: the smaller magnitude wins and a NaN
/// is propagated.
pub fn fminimumMagf(first: f32, second: f32) f32 {
    return minimumMagnitude(common.F, first, second);
}

/// The binary64 form of `fminimumMagh`: the smaller magnitude wins and a NaN
/// is propagated.
pub fn fminimumMagd(first: f64, second: f64) f64 {
    return minimumMagnitude(common.D, first, second);
}

/// The maximum family member spelled `number`. Its body is the same code as the
/// one behind `fmaxh`, so the two agree on every input; unlike `fmaximumh` it
/// does not propagate a NaN.
pub fn fmaximumNumh(first: f16, second: f16) f16 {
    return maximumNumber(common.H, first, second);
}

/// The binary32 form of `fmaximumNumh`; it shares `fmaxf`'s missing-NaN policy.
pub fn fmaximumNumf(first: f32, second: f32) f32 {
    return maximumNumber(common.F, first, second);
}

/// The binary64 form of `fmaximumNumh`; it shares `fmaxd`'s missing-NaN policy.
pub fn fmaximumNumd(first: f64, second: f64) f64 {
    return maximumNumber(common.D, first, second);
}

/// The minimum family member spelled `number`. Its body is the same code as the
/// one behind `fminh`, so the two agree on every input; unlike `fminimumh` it
/// does not propagate a NaN.
pub fn fminimumNumh(first: f16, second: f16) f16 {
    return minimumNumber(common.H, first, second);
}

/// The binary32 form of `fminimumNumh`; it shares `fminf`'s missing-NaN policy.
pub fn fminimumNumf(first: f32, second: f32) f32 {
    return minimumNumber(common.F, first, second);
}

/// The binary64 form of `fminimumNumh`; it shares `fmind`'s missing-NaN policy.
pub fn fminimumNumd(first: f64, second: f64) f64 {
    return minimumNumber(common.D, first, second);
}

/// Combines the magnitude comparison of `fmaximumMagh` with the missing-NaN
/// policy of `fmaxh`: the magnitude tie-break still applies, but a NaN operand
/// is skipped and two NaNs yield `first` quieted, which is the only difference
/// from `fmaximumMagh`.
pub fn fmaximumMagNumh(first: f16, second: f16) f16 {
    return maximumMagnitudeNumber(common.H, first, second);
}

/// The binary32 form of `fmaximumMagNumh`; a NaN operand is skipped rather
/// than propagated.
pub fn fmaximumMagNumf(first: f32, second: f32) f32 {
    return maximumMagnitudeNumber(common.F, first, second);
}

/// The binary64 form of `fmaximumMagNumh`; a NaN operand is skipped rather
/// than propagated.
pub fn fmaximumMagNumd(first: f64, second: f64) f64 {
    return maximumMagnitudeNumber(common.D, first, second);
}

/// Combines the magnitude comparison of `fminimumMagh` with the missing-NaN
/// policy of `fminh`, so the tie-break is the one described on `fmaximumMagh`
/// while a NaN operand is skipped instead of propagated.
pub fn fminimumMagNumh(first: f16, second: f16) f16 {
    return minimumMagnitudeNumber(common.H, first, second);
}

/// The binary32 form of `fminimumMagNumh`; a NaN operand is skipped rather
/// than propagated.
pub fn fminimumMagNumf(first: f32, second: f32) f32 {
    return minimumMagnitudeNumber(common.F, first, second);
}

/// The binary64 form of `fminimumMagNumh`; a NaN operand is skipped rather
/// than propagated.
pub fn fminimumMagNumd(first: f64, second: f64) f64 {
    return minimumMagnitudeNumber(common.D, first, second);
}
