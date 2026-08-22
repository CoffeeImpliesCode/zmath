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

fn expectExp10m1(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.exp10m1(value);
    const expected = if (T == f16)
        zmath.scalar.exp10m1h(value)
    else if (T == f32)
        zmath.scalar.exp10m1f(value)
    else
        zmath.scalar.exp10m1d(value);
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

fn expectErf(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.erf(value);
    const expected = if (T == f16)
        zmath.scalar.erfh(value)
    else if (T == f32)
        zmath.scalar.erff(value)
    else
        zmath.scalar.erfd(value);
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

fn expectSinh(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.sinh(value);
    const expected = if (T == f16)
        zmath.scalar.sinhh(value)
    else if (T == f32)
        zmath.scalar.sinhf(value)
    else
        zmath.scalar.sinhd(value);
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

fn expectCosh(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.cosh(value);
    const expected = if (T == f16)
        zmath.scalar.coshh(value)
    else if (T == f32)
        zmath.scalar.coshf(value)
    else
        zmath.scalar.coshd(value);
    comptime std.debug.assert(@TypeOf(got) == T);
    try std.testing.expectEqual(expected, got);
}
fn expectAcosh(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.acosh(value);
    const expected = if (T == f16)
        zmath.scalar.acoshh(value)
    else if (T == f32)
        zmath.scalar.acoshf(value)
    else
        zmath.scalar.acoshd(value);
    comptime std.debug.assert(@TypeOf(got) == T);
    try std.testing.expectEqual(expected, got);
}
fn expectAcos(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.acos(value);
    const expected = if (T == f16)
        zmath.scalar.acosh(value)
    else if (T == f32)
        zmath.scalar.acosf(value)
    else
        zmath.scalar.acosd(value);
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

fn expectAsinh(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.asinh(value);
    const expected = if (T == f16)
        zmath.scalar.asinhh(value)
    else if (T == f32)
        zmath.scalar.asinhf(value)
    else
        zmath.scalar.asinhd(value);
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

fn expectTan(comptime zmath: type, comptime T: type, value: T) !void {
    const got = zmath.tan(value);
    const expected = if (T == f16)
        zmath.scalar.tanh(value)
    else if (T == f32)
        zmath.scalar.tanf(value)
    else
        zmath.scalar.tand(value);
    comptime std.debug.assert(@TypeOf(got) == T);
    try std.testing.expectEqual(expected, got);
}

pub fn run(comptime zmath: type) !void {
    try expectExp(zmath, f16, @as(f16, 0.5));
    try expectExp(zmath, f32, @as(f32, 0.5));
    try expectExp(zmath, f64, @as(f64, 0.5));

    try expectExp10m1(zmath, f16, @as(f16, 0.5));
    try expectExp10m1(zmath, f32, @as(f32, 0.5));
    try expectExp10m1(zmath, f64, @as(f64, 0.5));
    try expectLog(zmath, f16, @as(f16, 1.5));
    try expectLog(zmath, f32, @as(f32, 1.5));
    try expectLog(zmath, f64, @as(f64, 1.5));

    try expectErf(zmath, f16, @as(f16, 0.5));
    try expectErf(zmath, f32, @as(f32, 0.5));
    try expectErf(zmath, f64, @as(f64, 0.5));

    try expectSin(zmath, f16, @as(f16, 0.5));
    try expectSin(zmath, f32, @as(f32, 0.5));
    try expectSin(zmath, f64, @as(f64, 0.5));
    try expectSinh(zmath, f16, @as(f16, 0.5));
    try expectSinh(zmath, f32, @as(f32, 0.5));
    try expectSinh(zmath, f64, @as(f64, 0.5));

    try expectCos(zmath, f16, @as(f16, 0.5));
    try expectCos(zmath, f32, @as(f32, 0.5));
    try expectCos(zmath, f64, @as(f64, 0.5));
    try expectTan(zmath, f16, @as(f16, 0.5));
    try expectTan(zmath, f32, @as(f32, 0.5));
    try expectTan(zmath, f64, @as(f64, 0.5));

    try expectCosh(zmath, f16, @as(f16, 0.5));
    try expectCosh(zmath, f32, @as(f32, 0.5));
    try expectCosh(zmath, f64, @as(f64, 0.5));

    try expectTanh(zmath, f16, @as(f16, 0.5));
    try expectTanh(zmath, f32, @as(f32, 0.5));
    try expectAcosh(zmath, f16, @as(f16, 1.5));
    try expectAcosh(zmath, f32, @as(f32, 1.5));
    try expectAcosh(zmath, f64, @as(f64, 1.5));
    try expectAcos(zmath, f16, @as(f16, 0.5));
    try expectAcos(zmath, f32, @as(f32, 0.5));
    try expectAcos(zmath, f64, @as(f64, 0.5));

    try expectTanh(zmath, f64, @as(f64, 0.5));
    try expectAsinh(zmath, f16, @as(f16, 0.5));
    try expectAsinh(zmath, f32, @as(f32, 0.5));
    try expectAsinh(zmath, f64, @as(f64, 0.5));
    try expectRsqrt(zmath, f16, @as(f16, 1.5));
    try expectRsqrt(zmath, f32, @as(f32, 1.5));
    try expectRsqrt(zmath, f64, @as(f64, 1.5));
}
