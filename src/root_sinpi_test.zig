const std = @import("std");
const zmath = @import("root.zig");

test "root exports nearest-even Sinpi family" {
    try std.testing.expectEqual(
        @as(u16, 0x39a8),
        @as(u16, @bitCast(zmath.sinpi(@as(f16, 0.25)))),
    );
    try std.testing.expectEqual(
        @as(u32, 0x3f35_04f3),
        @as(u32, @bitCast(zmath.sinpi(@as(f32, 0.25)))),
    );
    try std.testing.expectEqual(
        @as(u64, 0x3fe6_a09e_667f_3bcd),
        @as(u64, @bitCast(zmath.sinpi(@as(f64, 0.25)))),
    );
}
