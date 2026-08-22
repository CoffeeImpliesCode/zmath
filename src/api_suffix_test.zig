const std = @import("std");
const ieee = @import("ieee.zig");
const scalar = @import("scalar.zig");

test "IEEE exposes H, F, and D suffix exports" {
    comptime {
        std.debug.assert(@TypeOf(ieee.classifyH) == fn (f16) ieee.Class);
        std.debug.assert(@TypeOf(ieee.classifyF) == fn (f32) ieee.Class);
        std.debug.assert(@TypeOf(ieee.classifyD) == fn (f64) ieee.Class);

        _ = ieee.copySignH;
        _ = ieee.copySignF;
        _ = ieee.copySignD;
        _ = ieee.nextAfterH;
        _ = ieee.nextAfterF;
        _ = ieee.nextAfterD;
        _ = ieee.totalOrderH;
        _ = ieee.totalOrderF;
        _ = ieee.totalOrderD;
    }
}

test "scalar exposes H, F, and D suffix exports" {
    comptime {
        std.debug.assert(@TypeOf(scalar.absH(@as(f16, 0.0))) == f16);
        std.debug.assert(@TypeOf(scalar.absF(@as(f32, 0.0))) == f32);
        std.debug.assert(@TypeOf(scalar.absD(@as(f64, 0.0))) == f64);

        _ = scalar.ulpDistanceH;
        _ = scalar.ulpDistanceF;
        _ = scalar.ulpDistanceD;
    }
}
