const std = @import("std");
const c23 = @import("scalar_c23_types.zig");

// C23 integer rounding families with a fixed nearest-even or ties-away policy.
//
// Invalid inputs return an explicit result instead of relying on the floating-point
// environment or trapping during a float-to-integer conversion.

fn invalidResult() c23.IntegerResult {
    return .{ .value = 0, .ok = false };
}

fn isFiniteInI64Range(value: f64) bool {
    return value >= -0x1p63 and value < 0x1p63;
}

fn roundNearestEven(comptime T: type, value: T) c23.IntegerResult {
    @setFloatMode(.strict);
    const value_wide: f64 = @floatCast(value);
    if (!isFiniteInI64Range(value_wide)) return invalidResult();

    const truncated: i64 = @intFromFloat(value);
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
    return .{ .value = @intFromFloat(@round(value_wide)), .ok = true };
}

pub fn lrinth(value: f16) c23.IntegerResult {
    return roundNearestEven(f16, value);
}

pub fn lrintf(value: f32) c23.IntegerResult {
    return roundNearestEven(f32, value);
}

pub fn lrintd(value: f64) c23.IntegerResult {
    return roundNearestEven(f64, value);
}

pub fn llrinth(value: f16) c23.IntegerResult {
    return roundNearestEven(f16, value);
}

pub fn llrintf(value: f32) c23.IntegerResult {
    return roundNearestEven(f32, value);
}

pub fn llrintd(value: f64) c23.IntegerResult {
    return roundNearestEven(f64, value);
}

pub fn lroundh(value: f16) c23.IntegerResult {
    return roundTiesAway(f16, value);
}

pub fn lroundf(value: f32) c23.IntegerResult {
    return roundTiesAway(f32, value);
}

pub fn lroundd(value: f64) c23.IntegerResult {
    return roundTiesAway(f64, value);
}

pub fn llroundh(value: f16) c23.IntegerResult {
    return roundTiesAway(f16, value);
}

pub fn llroundf(value: f32) c23.IntegerResult {
    return roundTiesAway(f32, value);
}

pub fn llroundd(value: f64) c23.IntegerResult {
    return roundTiesAway(f64, value);
}
