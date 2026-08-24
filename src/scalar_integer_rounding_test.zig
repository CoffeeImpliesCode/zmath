const std = @import("std");
const c23 = @import("scalar_c23_types.zig");
const rounding = @import("scalar_integer_rounding.zig");

fn expectValid(expected: i64, actual: c23.IntegerResult) !void {
    try std.testing.expect(actual.ok);
    try std.testing.expectEqual(expected, actual.value);
}

fn expectInvalid(actual: c23.IntegerResult) !void {
    try std.testing.expect(!actual.ok);
    try std.testing.expectEqual(@as(i64, 0), actual.value);
}

fn expectNearestEven(
    comptime T: type,
    comptime operation: fn (T) c23.IntegerResult,
) !void {
    const cases = [_]struct {
        input: T,
        expected: i64,
    }{
        .{ .input = -3.75, .expected = -4 },
        .{ .input = -3.25, .expected = -3 },
        .{ .input = -0.0, .expected = 0 },
        .{ .input = 0.0, .expected = 0 },
        .{ .input = 3.25, .expected = 3 },
        .{ .input = 3.75, .expected = 4 },
        .{ .input = -2.5, .expected = -2 },
        .{ .input = -1.5, .expected = -2 },
        .{ .input = -0.5, .expected = 0 },
        .{ .input = 0.5, .expected = 0 },
        .{ .input = 1.5, .expected = 2 },
        .{ .input = 2.5, .expected = 2 },
    };

    for (cases) |test_case| {
        try expectValid(test_case.expected, operation(test_case.input));
    }
}

fn expectTiesAway(
    comptime T: type,
    comptime operation: fn (T) c23.IntegerResult,
) !void {
    const cases = [_]struct {
        input: T,
        expected: i64,
    }{
        .{ .input = -3.75, .expected = -4 },
        .{ .input = -3.25, .expected = -3 },
        .{ .input = -0.0, .expected = 0 },
        .{ .input = 0.0, .expected = 0 },
        .{ .input = 3.25, .expected = 3 },
        .{ .input = 3.75, .expected = 4 },
        .{ .input = -2.5, .expected = -3 },
        .{ .input = -1.5, .expected = -2 },
        .{ .input = -0.5, .expected = -1 },
        .{ .input = 0.5, .expected = 1 },
        .{ .input = 1.5, .expected = 2 },
        .{ .input = 2.5, .expected = 3 },
    };

    for (cases) |test_case| {
        try expectValid(test_case.expected, operation(test_case.input));
    }
}

fn Operations(comptime T: type) type {
    return struct {
        lrint: *const fn (T) c23.IntegerResult,
        llrint: *const fn (T) c23.IntegerResult,
        lround: *const fn (T) c23.IntegerResult,
        llround: *const fn (T) c23.IntegerResult,
    };
}

fn expectAllOperationsAt(
    comptime T: type,
    operations: Operations(T),
    input: T,
    expected: i64,
) !void {
    try expectValid(expected, operations.lrint(input));
    try expectValid(expected, operations.llrint(input));
    try expectValid(expected, operations.lround(input));
    try expectValid(expected, operations.llround(input));
}

fn expectAllOperationsInvalidAt(
    comptime T: type,
    operations: Operations(T),
    input: T,
) !void {
    try expectInvalid(operations.lrint(input));
    try expectInvalid(operations.llrint(input));
    try expectInvalid(operations.lround(input));
    try expectInvalid(operations.llround(input));
}

test "lrint and llrint use nearest-even for H F and D" {
    try expectNearestEven(f16, rounding.lrinth);
    try expectNearestEven(f16, rounding.llrinth);
    try expectNearestEven(f32, rounding.lrintf);
    try expectNearestEven(f32, rounding.llrintf);
    try expectNearestEven(f64, rounding.lrintd);
    try expectNearestEven(f64, rounding.llrintd);
}

test "lround and llround use ties-away for H F and D" {
    try expectTiesAway(f16, rounding.lroundh);
    try expectTiesAway(f16, rounding.llroundh);
    try expectTiesAway(f32, rounding.lroundf);
    try expectTiesAway(f32, rounding.llroundf);
    try expectTiesAway(f64, rounding.lroundd);
    try expectTiesAway(f64, rounding.llroundd);
}

// Explicit boundary bit patterns keep the valid and invalid sides of i64 conversion visible.
test "integer rounding reports each format's i64 boundary results" {
    const operations_h: Operations(f16) = .{
        .lrint = &rounding.lrinth,
        .llrint = &rounding.llrinth,
        .lround = &rounding.lroundh,
        .llround = &rounding.llroundh,
    };
    const operations_f: Operations(f32) = .{
        .lrint = &rounding.lrintf,
        .llrint = &rounding.llrintf,
        .lround = &rounding.lroundf,
        .llround = &rounding.llroundf,
    };
    const operations_d: Operations(f64) = .{
        .lrint = &rounding.lrintd,
        .llrint = &rounding.llrintd,
        .lround = &rounding.lroundd,
        .llround = &rounding.llroundd,
    };

    const minimum_h: f16 = @bitCast(@as(u16, 0xfbff));
    const maximum_h: f16 = @bitCast(@as(u16, 0x7bff));
    try expectAllOperationsAt(f16, operations_h, minimum_h, -65_504);
    try expectAllOperationsAt(f16, operations_h, maximum_h, 65_504);

    const minimum_f: f32 = @bitCast(@as(u32, 0xdf00_0000));
    const maximum_f: f32 = @bitCast(@as(u32, 0x5eff_ffff));
    try expectAllOperationsAt(f32, operations_f, minimum_f, -9_223_372_036_854_775_808);
    try expectAllOperationsAt(f32, operations_f, maximum_f, 9_223_371_487_098_961_920);

    const minimum_d: f64 = @bitCast(@as(u64, 0xc3e0_0000_0000_0000));
    const maximum_d: f64 = @bitCast(@as(u64, 0x43df_ffff_ffff_ffff));
    try expectAllOperationsAt(f64, operations_d, minimum_d, -9_223_372_036_854_775_808);
    try expectAllOperationsAt(f64, operations_d, maximum_d, 9_223_372_036_854_774_784);

    const nan_h: f16 = @bitCast(@as(u16, 0x7e01));
    const positive_infinity_h: f16 = @bitCast(@as(u16, 0x7c00));
    const negative_infinity_h: f16 = @bitCast(@as(u16, 0xfc00));
    try expectAllOperationsInvalidAt(f16, operations_h, nan_h);
    try expectAllOperationsInvalidAt(f16, operations_h, positive_infinity_h);
    try expectAllOperationsInvalidAt(f16, operations_h, negative_infinity_h);

    const nan_f: f32 = @bitCast(@as(u32, 0x7fc0_0001));
    const positive_infinity_f: f32 = @bitCast(@as(u32, 0x7f80_0000));
    const negative_infinity_f: f32 = @bitCast(@as(u32, 0xff80_0000));
    const positive_out_of_range_f: f32 = @bitCast(@as(u32, 0x5f00_0000));
    const negative_out_of_range_f: f32 = @bitCast(@as(u32, 0xdf00_0001));
    try expectAllOperationsInvalidAt(f32, operations_f, nan_f);
    try expectAllOperationsInvalidAt(f32, operations_f, positive_infinity_f);
    try expectAllOperationsInvalidAt(f32, operations_f, negative_infinity_f);
    try expectAllOperationsInvalidAt(f32, operations_f, positive_out_of_range_f);
    try expectAllOperationsInvalidAt(f32, operations_f, negative_out_of_range_f);

    const nan_d: f64 = @bitCast(@as(u64, 0x7ff8_0000_0000_0001));
    const positive_infinity_d: f64 = @bitCast(@as(u64, 0x7ff0_0000_0000_0000));
    const negative_infinity_d: f64 = @bitCast(@as(u64, 0xfff0_0000_0000_0000));
    const positive_out_of_range_d: f64 = @bitCast(@as(u64, 0x43e0_0000_0000_0000));
    const negative_out_of_range_d: f64 = @bitCast(@as(u64, 0xc3e0_0000_0000_0001));
    try expectAllOperationsInvalidAt(f64, operations_d, nan_d);
    try expectAllOperationsInvalidAt(f64, operations_d, positive_infinity_d);
    try expectAllOperationsInvalidAt(f64, operations_d, negative_infinity_d);
    try expectAllOperationsInvalidAt(f64, operations_d, positive_out_of_range_d);
    try expectAllOperationsInvalidAt(f64, operations_d, negative_out_of_range_d);
}

test "integer rounding preserves H F D result signatures" {
    comptime {
        std.debug.assert(@TypeOf(rounding.lrinth) == fn (f16) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.lrintf) == fn (f32) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.lrintd) == fn (f64) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.llrinth) == fn (f16) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.llrintf) == fn (f32) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.llrintd) == fn (f64) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.lroundh) == fn (f16) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.lroundf) == fn (f32) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.lroundd) == fn (f64) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.llroundh) == fn (f16) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.llroundf) == fn (f32) c23.IntegerResult);
        std.debug.assert(@TypeOf(rounding.llroundd) == fn (f64) c23.IntegerResult);
    }
}
