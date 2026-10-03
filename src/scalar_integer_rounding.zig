const std = @import("std");
const c23 = @import("scalar_c23_types.zig");

// C23 integer rounding families with a fixed nearest-even or ties-away policy.
//
// Invalid inputs return an explicit result instead of relying on the floating-point
// environment or trapping during a float-to-integer conversion.

fn invalidResult() c23.IntegerResult {
    return .{ .value = 0, .ok = false };
}

/// The tight `f64` window whose truncation lands in `i64`.
///
/// `0x1p63` is `maxInt(i64) + 1`, the first magnitude that does not fit, so it
/// is the exact threshold. The widest `f64` below it is `0x1.fffffffffffffp62`,
/// which is `2^63 - 1024` and so inside `i64`, and `-0x1p63` is exactly
/// `minInt(i64)`. Both comparisons are false for NaN, so the guard rejects NaN
/// and both infinities together with the out-of-range magnitudes.
fn isFiniteInI64Range(value: f64) bool {
    return value >= -0x1p63 and value < 0x1p63;
}

fn roundNearestEven(comptime T: type, value: T) c23.IntegerResult {
    @setFloatMode(.strict);
    const value_wide: f64 = @floatCast(value);
    if (!isFiniteInI64Range(value_wide)) return invalidResult();

    // The guard above admits only `[-0x1p63, 0x1p63)`, and every `f64` in that
    // window is at most `0x1.fffffffffffffp62`, so the truncation fits `i64`.
    // The clamp is the identity on that window and states the bound here.
    const truncated: i64 = @intFromFloat(std.math.clamp(value_wide, -0x1p63, 0x1.fffffffffffffp62));
    const truncated_value: T = @floatFromInt(truncated);
    const fraction = value - truncated_value;
    const odd = (truncated & 1) != 0;

    if (fraction > 0.5 or (fraction == 0.5 and odd)) {
        if (truncated == std.math.maxInt(i64)) return invalidResult();
        return .{ .value = truncated + 1, .ok = true };
    }
    if (fraction < -0.5 or (fraction == -0.5 and odd)) {
        if (truncated == std.math.minInt(i64)) return invalidResult();
        return .{ .value = truncated - 1, .ok = true };
    }
    return .{ .value = truncated, .ok = true };
}

fn roundTiesAway(comptime T: type, value: T) c23.IntegerResult {
    @setFloatMode(.strict);
    const value_wide: f64 = @floatCast(value);
    if (!isFiniteInI64Range(value_wide)) return invalidResult();

    // At the i64 limits, every representable input is integral, so rounding stays in range.
    // `@round` carries a value between two integer endpoints onto an integer
    // between them, so the rounded value is in `[-0x1p63,
    // 0x1.fffffffffffffp62]` and this clamp is the identity on it.
    return .{
        .value = @intFromFloat(std.math.clamp(@round(value_wide), -0x1p63, 0x1.fffffffffffffp62)),
        .ok = true,
    };
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour. The direction is fixed in code; the conversion never traps, and
/// it neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, and the tie-breaking step stays
/// inside the `i64` range, so `ok` is false only for a NaN, either infinity, or
/// a finite magnitude outside that window. The `value` field is then 0 and must
/// not be read as a converted result. `llrinth` returns the same `i64` width.
pub fn lrinth(value: f16) c23.IntegerResult {
    return roundNearestEven(f16, value);
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour. The direction is fixed in code; the conversion never traps, and
/// it neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, and the tie-breaking step stays
/// inside the `i64` range, so `ok` is false only for a NaN, either infinity, or
/// a finite magnitude outside that window. The `value` field is then 0 and must
/// not be read as a converted result. `llrintf` returns the same `i64` width.
pub fn lrintf(value: f32) c23.IntegerResult {
    return roundNearestEven(f32, value);
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour. The direction is fixed in code; the conversion never traps, and
/// it neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, and the tie-breaking step stays
/// inside the `i64` range, so `ok` is false only for a NaN, either infinity, or
/// a finite magnitude outside that window. The `value` field is then 0 and must
/// not be read as a converted result. `llrintd` returns the same `i64` width.
pub fn lrintd(value: f64) c23.IntegerResult {
    return roundNearestEven(f64, value);
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour. The direction is fixed in code; the conversion never traps, and
/// it neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, and the tie-breaking step stays
/// inside the `i64` range, so `ok` is false only for a NaN, either infinity, or
/// a finite magnitude outside that window. The `value` field is then 0 and must
/// not be read as a converted result. The `ll` prefix names the C23 spelling
/// only; `lrinth` returns the same `i64` width.
pub fn llrinth(value: f16) c23.IntegerResult {
    return roundNearestEven(f16, value);
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour. The direction is fixed in code; the conversion never traps, and
/// it neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, and the tie-breaking step stays
/// inside the `i64` range, so `ok` is false only for a NaN, either infinity, or
/// a finite magnitude outside that window. The `value` field is then 0 and must
/// not be read as a converted result. The `ll` prefix names the C23 spelling
/// only; `lrintf` returns the same `i64` width.
pub fn llrintf(value: f32) c23.IntegerResult {
    return roundNearestEven(f32, value);
}

/// Rounds `value` to the nearest integer, resolving a halfway case to the even
/// neighbour. The direction is fixed in code; the conversion never traps, and
/// it neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, and the tie-breaking step stays
/// inside the `i64` range, so `ok` is false only for a NaN, either infinity, or
/// a finite magnitude outside that window. The `value` field is then 0 and must
/// not be read as a converted result. The `ll` prefix names the C23 spelling
/// only; `lrintd` returns the same `i64` width.
pub fn llrintd(value: f64) c23.IntegerResult {
    return roundNearestEven(f64, value);
}

/// Rounds `value` to the nearest integer, sending a halfway case away from
/// zero. The direction is fixed in code; the conversion never traps, and it
/// neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, so `ok` is false only for a
/// NaN, either infinity, or a finite magnitude outside that window. The `value`
/// field is then 0 and must not be read as a converted result. `llroundh`
/// returns the same `i64` width.
pub fn lroundh(value: f16) c23.IntegerResult {
    return roundTiesAway(f16, value);
}

/// Rounds `value` to the nearest integer, sending a halfway case away from
/// zero. The direction is fixed in code; the conversion never traps, and it
/// neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, so `ok` is false only for a
/// NaN, either infinity, or a finite magnitude outside that window. The `value`
/// field is then 0 and must not be read as a converted result. `llroundf`
/// returns the same `i64` width.
pub fn lroundf(value: f32) c23.IntegerResult {
    return roundTiesAway(f32, value);
}

/// Rounds `value` to the nearest integer, sending a halfway case away from
/// zero. The direction is fixed in code; the conversion never traps, and it
/// neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, so `ok` is false only for a
/// NaN, either infinity, or a finite magnitude outside that window. The `value`
/// field is then 0 and must not be read as a converted result. `llroundd`
/// returns the same `i64` width.
pub fn lroundd(value: f64) c23.IntegerResult {
    return roundTiesAway(f64, value);
}

/// Rounds `value` to the nearest integer, sending a halfway case away from
/// zero. The direction is fixed in code; the conversion never traps, and it
/// neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, so `ok` is false only for a
/// NaN, either infinity, or a finite magnitude outside that window. The `value`
/// field is then 0 and must not be read as a converted result. `llroundh`
/// returns the same `i64` width.
pub fn llroundh(value: f16) c23.IntegerResult {
    return roundTiesAway(f16, value);
}

/// Rounds `value` to the nearest integer, sending a halfway case away from
/// zero. The direction is fixed in code; the conversion never traps, and it
/// neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, so `ok` is false only for a
/// NaN, either infinity, or a finite magnitude outside that window. The `value`
/// field is then 0 and must not be read as a converted result. The `ll` prefix
/// names the C23 spelling only; `lroundf` returns the same `i64` width.
pub fn llroundf(value: f32) c23.IntegerResult {
    return roundTiesAway(f32, value);
}

/// Rounds `value` to the nearest integer, sending a halfway case away from
/// zero. The direction is fixed in code; the conversion never traps, and it
/// neither reads nor changes the floating-point environment.
///
/// Every input in `[-0x1p63, 0x1p63)` converts, so `ok` is false only for a
/// NaN, either infinity, or a finite magnitude outside that window. The `value`
/// field is then 0 and must not be read as a converted result. The `ll` prefix
/// names the C23 spelling only; `lroundd` returns the same `i64` width.
pub fn llroundd(value: f64) c23.IntegerResult {
    return roundTiesAway(f64, value);
}
