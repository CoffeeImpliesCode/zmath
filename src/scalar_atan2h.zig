// Correctly rounded binary16 atan2 for nearest-even rounding.
// Derived from CORE-MATH binary16/atan2/atan2f16.c.
const atan_impl = @import("scalar_atan.zig");

pub fn atan2h(y: f16, x: f16) f16 {
    @setFloatMode(.strict);
    const yb: u16 = @bitCast(y);
    const xb: u16 = @bitCast(x);
    const ay = yb & 0x7fff;
    const ax = xb & 0x7fff;
    if (ay > 0x7c00) return @bitCast(yb | 0x0200);
    if (ax > 0x7c00) return @bitCast(xb | 0x0200);

    const sy = (yb & 0x8000) != 0;
    const sx = (xb & 0x8000) != 0;
    const pi: f16 = @floatCast(0x1.921fb54442d18p+1);
    const pi2: f16 = @floatCast(0x1.921fb54442d18p+0);
    const pi4_3: f16 = @floatCast(0x1.2d97c7f3321d2p+1);

    if (ay == 0x7c00) {
        if (ax == 0x7c00) return if (sy) -pi4_3 else pi4_3;
        return if (sy) -pi2 else pi2;
    }
    if (ax == 0x7c00) {
        if (sx) return if (sy) -pi else pi;
        return @bitCast(yb & 0x8000);
    }
    if (ay == 0 and ax == 0) {
        if (sx) return if (sy) @bitCast(@as(u16, 0x8000)) else pi;
        return @bitCast(yb & 0x8000);
    }
    if (ay == 0) {
        if (sx) return if (sy) -pi else pi;
        return @bitCast(yb & 0x8000);
    }
    if (ax == 0) return if (sy) -pi2 else pi2;

    const yf: f16 = @bitCast(ay);
    const xf: f16 = @bitCast(ax);
    var angle: f16 = if (ay > ax)
        pi2 - atan_impl.atanh(xf / yf)
    else
        atan_impl.atanh(yf / xf);
    if (sx) angle = pi - angle;
    if (sy) angle = -angle;
    return angle;
}