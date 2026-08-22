const std = @import("std");
const ieee = @import("ieee.zig");
const scalar = @import("scalar.zig");

test "IEEE exposes h, f, and d suffix exports" {
    comptime {
        std.debug.assert(@TypeOf(ieee.classifyh) == fn (f16) ieee.Class);
        std.debug.assert(@TypeOf(ieee.classifyf) == fn (f32) ieee.Class);
        std.debug.assert(@TypeOf(ieee.classifyd) == fn (f64) ieee.Class);

        _ = ieee.copySignh;
        _ = ieee.copySignf;
        _ = ieee.copySignd;
        _ = ieee.nextAfterh;
        _ = ieee.nextAfterf;
        _ = ieee.nextAfterd;
        _ = ieee.totalOrderh;
        _ = ieee.totalOrderf;
        _ = ieee.totalOrderd;
    }
}

test "scalar exposes h, f, and d suffix exports" {
    comptime {
        std.debug.assert(@TypeOf(scalar.absh(@as(f16, 0.0))) == f16);
        std.debug.assert(@TypeOf(scalar.absf(@as(f32, 0.0))) == f32);
        std.debug.assert(@TypeOf(scalar.absd(@as(f64, 0.0))) == f64);

        _ = scalar.ulpDistanceh;
        _ = scalar.ulpDistancef;
        _ = scalar.ulpDistanced;

        _ = scalar.rsqrth;
        _ = scalar.rsqrtf;
        _ = scalar.rsqrtd;

        _ = scalar.coshh;
        _ = scalar.coshf;
        _ = scalar.coshd;
        _ = scalar.acoshh;
        _ = scalar.acoshf;
        _ = scalar.acoshd;
        _ = scalar.acosh;
        _ = scalar.acosf;
        _ = scalar.acosd;
        _ = scalar.exp10m1h;
        _ = scalar.exp10m1f;
        _ = scalar.exp10m1d;
    }
}
