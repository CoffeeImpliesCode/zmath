const std = @import("std");
const scalar = @import("scalar.zig");

const tolerance_ulps: comptime_float = 4.0;

// Use four target-format epsilons as the relative tolerance. For subnormal
// results, use four smallest-subnormal steps as the absolute tolerance.
// These std.math references are an approximate secondary check, not an
// exact-bit oracle.
fn stdlibExp(comptime T: type, input: T) T {
    return std.math.exp(input);
}

fn stdlibLog(comptime T: type, input: T) T {
    if (T == f16) {
        const promoted_input: f32 = @floatCast(input);
        const promoted_result = std.math.log(f32, std.math.e, promoted_input);
        return @floatCast(promoted_result);
    }
    return std.math.log(T, @as(T, std.math.e), input);
}

fn stdlibSin(comptime T: type, input: T) T {
    return std.math.sin(input);
}

fn expectStdlibAgreement(comptime T: type, expected: T, actual: T) !void {
    const relative_tolerance: T = tolerance_ulps * std.math.floatEps(T);
    const absolute_tolerance: T = tolerance_ulps * std.math.floatTrueMin(T);
    const magnitude = @max(@abs(expected), @abs(actual));

    if (magnitude < std.math.floatMin(T)) {
        try std.testing.expectApproxEqAbs(expected, actual, absolute_tolerance);
    } else {
        try std.testing.expectApproxEqRel(expected, actual, relative_tolerance);
    }
}

fn checkExpSamples(comptime T: type, samples: []const T, function: fn (T) T) !void {
    for (samples) |input| {
        try expectStdlibAgreement(T, stdlibExp(T, input), function(input));
    }
}

fn checkLogSamples(comptime T: type, samples: []const T, function: fn (T) T) !void {
    for (samples) |input| {
        try expectStdlibAgreement(T, stdlibLog(T, input), function(input));
    }
}

fn checkSinSamples(comptime T: type, samples: []const T, function: fn (T) T) !void {
    for (samples) |input| {
        try expectStdlibAgreement(T, stdlibSin(T, input), function(input));
    }
}

test "stdlib agreement: exph finite samples" {
    const minimum_subnormal: f16 = @bitCast(@as(u16, 0x0001));
    const samples = [_]f16{
        minimum_subnormal,
        -minimum_subnormal,
        0.001,
        -0.001,
        1.0,
        -1.0,
        10.0,
        -10.0,
    };

    try checkExpSamples(f16, &samples, scalar.exph);
}

test "stdlib agreement: logh finite samples" {
    const minimum_subnormal: f16 = @bitCast(@as(u16, 0x0001));
    const below_one: f16 = @bitCast(@as(u16, 0x3bff));
    const above_one: f16 = @bitCast(@as(u16, 0x3c01));
    const maximum_finite: f16 = @bitCast(@as(u16, 0x7bff));
    const samples = [_]f16{
        minimum_subnormal,
        below_one,
        above_one,
        0.5,
        2.0,
        10.0,
        maximum_finite,
    };

    try checkLogSamples(f16, &samples, scalar.logh);
}

test "stdlib agreement: sinh finite samples" {
    const minimum_subnormal: f16 = @bitCast(@as(u16, 0x0001));
    const maximum_finite: f16 = @bitCast(@as(u16, 0x7bff));
    const samples = [_]f16{
        minimum_subnormal,
        -minimum_subnormal,
        0.001,
        -0.001,
        1.0,
        -1.0,
        std.math.pi / 2.0,
        1000.0,
        maximum_finite,
    };

    try checkSinSamples(f16, &samples, scalar.sinh);
}

test "stdlib agreement: expf finite samples" {
    const minimum_subnormal: f32 = @bitCast(@as(u32, 0x0000_0001));
    const samples = [_]f32{
        minimum_subnormal,
        -minimum_subnormal,
        0.000_001,
        -0.000_001,
        1.0,
        -1.0,
        80.0,
        -80.0,
    };

    try checkExpSamples(f32, &samples, scalar.expf);
}

test "stdlib agreement: logf finite samples" {
    const minimum_subnormal: f32 = @bitCast(@as(u32, 0x0000_0001));
    const below_one: f32 = @bitCast(@as(u32, 0x3f7f_ffff));
    const above_one: f32 = @bitCast(@as(u32, 0x3f80_0001));
    const maximum_finite: f32 = @bitCast(@as(u32, 0x7f7f_ffff));
    const samples = [_]f32{
        minimum_subnormal,
        below_one,
        above_one,
        0.5,
        2.0,
        10.0,
        maximum_finite,
    };

    try checkLogSamples(f32, &samples, scalar.logf);
}

test "stdlib agreement: sinf finite samples" {
    const minimum_subnormal: f32 = @bitCast(@as(u32, 0x0000_0001));
    const maximum_finite: f32 = @bitCast(@as(u32, 0x7f7f_ffff));
    const samples = [_]f32{
        minimum_subnormal,
        -minimum_subnormal,
        0.000_001,
        -0.000_001,
        1.0,
        -1.0,
        std.math.pi / 2.0,
        10_000.25,
        maximum_finite,
    };

    try checkSinSamples(f32, &samples, scalar.sinf);
}

test "stdlib agreement: expd finite samples" {
    const minimum_subnormal: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0001));
    const samples = [_]f64{
        minimum_subnormal,
        -minimum_subnormal,
        0.000_000_000_001,
        -0.000_000_000_001,
        1.0,
        -1.0,
        700.0,
        -700.0,
    };

    try checkExpSamples(f64, &samples, scalar.expd);
}

test "stdlib agreement: logd finite samples" {
    const minimum_subnormal: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0001));
    const below_one: f64 = @bitCast(@as(u64, 0x3fef_ffff_ffff_ffff));
    const above_one: f64 = @bitCast(@as(u64, 0x3ff0_0000_0000_0001));
    const maximum_finite: f64 = @bitCast(@as(u64, 0x7fef_ffff_ffff_ffff));
    const samples = [_]f64{
        minimum_subnormal,
        below_one,
        above_one,
        0.5,
        2.0,
        10.0,
        maximum_finite,
    };

    try checkLogSamples(f64, &samples, scalar.logd);
}

test "stdlib agreement: sind finite samples" {
    const minimum_subnormal: f64 = @bitCast(@as(u64, 0x0000_0000_0000_0001));
    const maximum_finite: f64 = @bitCast(@as(u64, 0x7fef_ffff_ffff_ffff));
    const samples = [_]f64{
        minimum_subnormal,
        -minimum_subnormal,
        0.000_000_000_001,
        -0.000_000_000_001,
        1.0,
        -1.0,
        std.math.pi / 2.0,
        1.0e20,
        maximum_finite,
    };

    try checkSinSamples(f64, &samples, scalar.sind);
}
