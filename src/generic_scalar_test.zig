const std = @import("std");

fn expectExp(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.exp(value);
    const expected = if (T == f16)
        zmath.scalar.exph(value)
    else if (T == f32)
        zmath.scalar.expf(value)
    else
        zmath.scalar.expd(value);
    comptime std.debug.assert(@TypeOf(got) == T);
    try std.testing.expectEqual(expected, got);
}

fn expectLog(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.log(value);
    const expected = if (T == f16)
        zmath.scalar.logh(value)
    else if (T == f32)
        zmath.scalar.logf(value)
    else
        zmath.scalar.logd(value);
    comptime std.debug.assert(@TypeOf(got) == T);
    try std.testing.expectEqual(expected, got);
}

fn expectSin(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.sin(value);
    const expected = if (T == f16)
        zmath.scalar.sinh(value)
    else if (T == f32)
        zmath.scalar.sinf(value)
    else
        zmath.scalar.sind(value);
    comptime std.debug.assert(@TypeOf(got) == T);
    try std.testing.expectEqual(expected, got);
}

fn expectCos(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.cos(value);
    const expected = if (T == f16)
        zmath.scalar.cosh(value)
    else if (T == f32)
        zmath.scalar.cosf(value)
    else
        zmath.scalar.cosd(value);
    comptime std.debug.assert(@TypeOf(got) == T);
    try std.testing.expectEqual(expected, got);
}

fn expectTanh(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.tanh(value);
    const expected = if (T == f16)
        zmath.scalar.tanhh(value)
    else if (T == f32)
        zmath.scalar.tanhf(value)
    else
        zmath.scalar.tanhd(value);
    comptime std.debug.assert(@TypeOf(got) == T);
    try std.testing.expectEqual(expected, got);
}

fn expectRsqrt(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.rsqrt(value);
    const expected = if (T == f16)
        zmath.scalar.rsqrth(value)
    else if (T == f32)
        zmath.scalar.rsqrtf(value)
    else
        zmath.scalar.rsqrtd(value);
    comptime std.debug.assert(@TypeOf(got) == T);
    try std.testing.expectEqual(expected, got);
}

pub fn run(comptime zmath: type) !void {
    try expectExp(zmath, f16, @as(f16, 0.5));
    try expectExp(zmath, f32, @as(f32, 0.5));
    try expectExp(zmath, f64, @as(f64, 0.5));

    try expectLog(zmath, f16, @as(f16, 1.5));
    try expectLog(zmath, f32, @as(f32, 1.5));
    try expectLog(zmath, f64, @as(f64, 1.5));

    try expectSin(zmath, f16, @as(f16, 0.5));
    try expectSin(zmath, f32, @as(f32, 0.5));
    try expectSin(zmath, f64, @as(f64, 0.5));

    try expectCos(zmath, f16, @as(f16, 0.5));
    try expectCos(zmath, f32, @as(f32, 0.5));
    try expectCos(zmath, f64, @as(f64, 0.5));
    try expectTanh(zmath, f16, @as(f16, 0.5));
    try expectTanh(zmath, f32, @as(f32, 0.5));
    try expectTanh(zmath, f64, @as(f64, 0.5));
    try expectRsqrt(zmath, f16, @as(f16, 1.5));
    try expectRsqrt(zmath, f32, @as(f32, 1.5));
    try expectRsqrt(zmath, f64, @as(f64, 1.5));
}
