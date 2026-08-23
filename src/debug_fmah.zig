const std = @import("std");
const scalar = @import("scalar.zig");

pub fn main() !void {
    const quiet_nan_payload: f16 = @bitCast(@as(u16, 0x7e00 | 0x0123));
    std.debug.print("quiet_nan_payload bits: {x}\n", .{@as(u16, @bitCast(quiet_nan_payload))});
    const result = scalar.fmah(quiet_nan_payload, 1.0, 2.0);
    std.debug.print("result bits: {x}\n", .{@as(u16, @bitCast(result))});
    if (@as(u16, @bitCast(quiet_nan_payload)) == @as(u16, @bitCast(result))) {
        std.debug.print("equal\n", .{});
    } else {
        std.debug.print("not equal\n", .{});
    }
}
