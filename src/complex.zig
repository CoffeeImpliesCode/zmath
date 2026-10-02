//! Complex arithmetic over `f32` and `f64`.
//!
//! This module owns `Complexf` and `Complexd` and the C99/C23 complex surface
//! over them. Every numerical kernel is a `zmath.scalar` kernel; the module
//! contributes only the complex decompositions, the principal-branch
//! selections, and the overflow-safe reductions a scalar kernel cannot
//! express. No host libm is called, and nothing here reads, sets, or restores
//! the floating-point environment, so a caller must keep round-to-nearest,
//! ties-to-even active.
//!
//! No operation returns `null` or an error union. A C mathematical pole is a
//! value: the IEEE division of the decomposition, which is a signed infinity
//! where the numerator survives and the canonical quiet NaN where it does not.
//!
//! ## Operation table
//!
//! Every function name below is the `f` or `d` spelling of one row. A NaN
//! operand keeps its sign and diagnostic payload and gains the quiet bit; an
//! invalid non-NaN input, and every NaN column marked "canonical", is the
//! positive canonical quiet NaN of the result format: `7fc00000` for `f` and
//! `7ff8000000000000` for `d`.
//!
//! | Family | Result | Branch cut | Principal range | Signed zero | NaN | Infinity | Pole | Overflow | Underflow |
//! | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
//! | `realf` `reald` `imagf` `imagd` | the stored component | none | none | stored exactly | stored exactly | stored exactly | none | none | stored exactly |
//! | `conjf` `conjd` | `Complex` | none | none | imaginary part negated, so `-0` becomes `+0` | the real part is stored; a NaN imaginary part keeps its payload and flips its sign | stored exactly | none | none | none |
//! | `projf` `projd` | `Complex` | none | real axis | imaginary part is always `+0` | first NaN component, quieted, in both parts | `(+inf, +0)` | none | none | none |
//! | `absf` `absd` | magnitude `f32`/`f64` | none | `[+0, +inf]` | `+0` for `+0 + i0` | canonical | `+inf` | none | `+inf` when the true magnitude rounds to infinity | `+0` when the true magnitude rounds to zero |
//! | `argf` `argd` | angle `f32`/`f64` | none | `(-pi, +pi]` | `+0` for `+0 + i0`, `-0` for `+0 - i0`, `+pi` for `-0 + i0`, `-pi` for `-0 - i0` | first NaN component, quieted | `scalar.atan2` quadrant rules | none | none | `scalar.atan2` rounding |
//! | `addf` `addd` `subf` `subd` | `Complex` | none | none | IEEE addition of the components | per component, the first NaN of that component, quieted | IEEE addition, so one `inf` component and one opposing `inf` give canonical | none | `+inf` per IEEE | `+0` per IEEE |
//! | `mulf` `muld` | `Complex` | none | none | IEEE multiplication of the products | both parts, first NaN of the four operands, quieted | IEEE multiplication, so `inf * 0` gives canonical | none | `+inf` per IEEE | `+0` per IEEE |
//! | `divf` `divd` | `Complex` | none | none | IEEE division of the scaled numerator, so a zero `z.re` with a nonzero `z.im` keeps the divisor's zero sign | both parts, first NaN of the four operands, quieted | an infinite divisor gives the signed zeros of the leading term of its reciprocal; an infinite numerator over an infinite divisor and a divisor with two infinite components give canonical | `+0 + i0` divisor gives an infinite real part signed by the naive numerator and a canonical imaginary part; `0/0` gives canonical in both parts | `+inf` per IEEE | `+0` only where the exact quotient underflows, so a subnormal quotient is exact |
//! | `expf` `expd` | `Complex` | none | real part unbounded, imaginary in `(-pi, pi]` | `z.im == 0` returns `e^z.re` with `z.im` verbatim, so `exp(x - i0) = e^x - i0` | canonical in both parts | infinite imaginary part gives canonical in both parts; `+inf` real with nonzero imaginary gives `(+inf, canonical)`; `+inf + i0` gives `(+inf, +0)` and `-inf + i0` gives `(+0, +0)` | none | real part is `e^z.re`, so `+inf` once `z.re` exceeds the range | `+0` real once `z.re` is below the range |
//! | `logf` `logd` | `Complex` | none | real in `(-inf, +inf]`, imaginary in `(-pi, pi]` | `z = +0 + i0` gives `(-inf, +0)` and `z = +0 - i0` gives `(-inf, -0)` | canonical in both parts | infinite real part gives `(+inf, +0)`; infinite imaginary part gives `(+inf, +-pi/2)`; both infinite gives `(+inf, canonical)` | none | `+inf` for an infinite input only | `-inf` real when the magnitude underflows to zero |
//! | `log10f` `log10d` | `Complex` | as `log`, because the imaginary part is that of `log` divided by `ln 10` | real in `(-inf, +inf]`, imaginary in `(-pi/ln 10, pi/ln 10]` | inherits `log` | canonical in both parts | inherits `log` | none | inherits `log` | inherits `log` |
//! | `sqrtf` `sqrtd` | `Complex` | none | real in `[+0, +inf)`, imaginary signed by `z.im` | `+0 + i+-0` returns the input; negative real part gives `+0` real part | canonical in both parts | infinite real part with infinite imaginary part gives `(+inf, +inf)`; `z.re` infinite gives `(+inf, +-inf)`; `z.im` infinite gives `(+inf, +inf)`; `+inf + i+-0` gives `(+inf, +-0)` | none | `+inf` in the part that the input makes unbounded | `+0` in the part that the input makes vanish |
//! | `powf` `powd` | `Complex` | negative real axis for a non-integral real exponent, taken at `arg = +pi` | principal value of `exp(y log x)` | `y = 0 + i+-0` gives `(1, +-0)`; a negative real base with an exact integer exponent gives `+0` imaginary | canonical in both parts | `+0 + i0` base with a negative real or nonzero imaginary exponent gives `(+inf, canonical)`; any other infinite base or exponent gives canonical in both parts | `+0 + i0` base with exponent `+-0` gives `(1, +-0)` | inherits `exp` | inherits `exp` |
//! | `sinf` `sind` | `Complex` | none | unbounded in both parts | `z.im == 0` returns `(sin z.re, z.im)` verbatim | canonical in both parts | infinite imaginary part gives `(canonical, +-inf)`; infinite real part gives canonical in both parts | none | `cosh` and `sinh` overflow drives the scaled parts to `+inf` | a part whose trigonometric factor underflows gives `+0` |
//! | `cosf` `cosd` | `Complex` | none | unbounded in both parts | `z.im == 0` returns `(cos z.re, z.im)` verbatim | canonical in both parts | infinite imaginary part gives `(canonical, -+-inf)`; infinite real part gives canonical in both parts | none | as `sin` | as `sin` |
//! | `tanf` `tand` | `Complex` | none | unbounded in both parts | `z.im == 0` returns `(tan z.re, z.im)` verbatim | canonical in both parts | infinite imaginary part gives `(canonical, +-inf)`; infinite real part gives canonical in both parts | denominator `cos 2a + cosh 2b` vanishes only at `b = 0`, which the `z.im == 0` path serves, so no finite binary32 or binary64 argument reaches it | `cosh 2b` overflow gives a `+-0` real part and `+-1` imaginary part | none |
//! | `asinf` `asind` | `Complex` | real axis outside `[-1, +1]`, approached from the upper half plane | real in `[-pi/2, +pi/2]`, imaginary signed by `z.im` | `z = 0` gives `(+0, -0)` and a real argument inside the unit interval gives `-0` imaginary | canonical in both parts | canonical in both parts | none | none | none |
//! | `acosf` `acosd` | `Complex` | real axis outside `[-1, +1]`, approached from the lower half plane | real in `[0, +pi]`, imaginary signed by `-z.im` | `z = 0` gives `(+pi/2, +0)`, `z = -1` gives `(+pi, +0)` and `z = 1` gives `(+0, +0)` | canonical in both parts | canonical in both parts | none | none | none |
//! | `atanf` `atand` | `Complex` | imaginary axis outside `[-i, +i]` | real in `(-pi/2, +pi/2)` | `z = 0` gives `(+0, +0)` and a real argument gives `+0` imaginary | canonical in both parts | canonical in both parts | none | none | none |
//! | `sinhf` `sinhd` | `Complex` | none | unbounded in both parts | `z.im == 0` returns `(sinh z.re, z.im)` verbatim | canonical in both parts | canonical in both parts | none | `sinh` and `cosh` overflow drives the parts to `+inf` | a part whose trigonometric factor underflows gives `+0` |
//! | `coshf` `coshd` | `Complex` | none | unbounded in both parts | `z.im == 0` returns `(cosh z.re, z.im)` verbatim | canonical in both parts | canonical in both parts | none | as `sinh` | as `sinh` |
//! | `tanhf` `tanhd` | `Complex` | none | real in `(-1, +1)`, imaginary unbounded | `z.im == 0` returns `(tanh z.re, z.im)` verbatim | canonical in both parts | infinite real part gives `(+-1, +0)`; infinite imaginary part gives canonical in both parts | denominator `cosh 2a + cos 2b` vanishes at `a = 0`, `b = +-pi/2`, which gives `(canonical, +-inf)` | `cosh 2a` overflow gives `(+-1, +0)` | none |
//! | `asinhf` `asinhd` | `Complex` | imaginary axis outside `[-i, +i]` | real signed by `z.re`, imaginary in `[-pi/2, +pi/2]` | `z = 0` gives `(0, 0)` | canonical in both parts | canonical in both parts | none | none | none |
//! | `acoshf` `acoshd` | `Complex` | real axis below `+1` | real in `[+0, +inf)`, imaginary in `[0, +pi]` | `z = 1` gives `(+0, +0)` | canonical in both parts | canonical in both parts | none | none | none |
//! | `atanhf` `atanhd` | `Complex` | real axis outside `[-1, +1]` | real in `(-1, +1)`, imaginary in `(-pi/2, +pi/2]` | `z = 0` gives `(0, 0)` | canonical in both parts | canonical in both parts | `z = 1` gives `(+inf, +0)` and `z = -1` gives `(-inf, +0)` | none | none |

const ieee = @import("ieee.zig");
const scalar = @import("scalar.zig");
const scalar_div = @import("scalar_div.zig");

/// A complex value in binary32.
pub const Complexf = struct { re: f32, im: f32 };

/// A complex value in binary64.
pub const Complexd = struct { re: f64, im: f64 };

// ---------------------------------------------------------------------------
// Format bindings
// ---------------------------------------------------------------------------

/// The complex type whose components are `T`.
fn ComplexOf(comptime T: type) type {
    if (T == f32) return Complexf;
    if (T == f64) return Complexd;
    @compileError("zmath.complex supports only f32 and f64");
}

/// The positive infinity of `T`.
fn positiveInfinity(comptime T: type) T {
    if (T == f32) return @bitCast(@as(u32, 0x7f80_0000));
    return @bitCast(@as(u64, 0x7ff0_0000_0000_0000));
}

/// The positive canonical quiet NaN of `T`.
fn canonicalNaN(comptime T: type) T {
    if (T == f32) return @bitCast(@as(u32, 0x7fc0_0000));
    return @bitCast(@as(u64, 0x7ff8_0000_0000_0000));
}

/// Round-to-nearest pi/2 as `T`.
fn halfPiOf(comptime T: type) T {
    if (T == f32) return 1.5707963267948966;
    return 1.5707963267948966;
}

/// Round-to-nearest ln(10) as `T`.
fn lnTenOf(comptime T: type) T {
    if (T == f32) return 2.302585092994046;
    return 2.302585092994046;
}

/// Round-to-nearest sqrt(2) as `T`.
fn sqrtTwoOf(comptime T: type) T {
    if (T == f32) return 1.4142135623730951;
    return 1.4142135623730951;
}

/// The smallest magnitude whose square is still a normal number, that is
/// sqrt(2^-126) for `f32` and sqrt(2^-1022) for `f64`.
fn sqrtMinNormalOf(comptime T: type) T {
    if (T == f32) return 1.0842022e-19;
    return 1.4916681462400413e-154;
}

/// The largest magnitude for which `m * hypot(1, s/m)` with `s <= m` cannot
/// round to infinity, that is sqrt(max / 2) of the format.
fn magnitudeLimitOf(comptime T: type) T {
    if (T == f32) return 1.3043818e19;
    return 9.480751908109176e153;
}

/// The `zmath.scalar` kernels and format limits that a generic body needs.
fn Format(comptime T: type) type {
    return struct {
        const infinity = positiveInfinity(T);
        const nan = canonicalNaN(T);
        const half_pi = halfPiOf(T);
        const ln_10 = lnTenOf(T);
        const sqrt_two = sqrtTwoOf(T);
        const sqrt_min_normal = sqrtMinNormalOf(T);
        const magnitude_limit = magnitudeLimitOf(T);

        inline fn exp(x: T) T {
            if (T == f32) return scalar.expf(x);
            return scalar.expd(x);
        }

        inline fn log(x: T) T {
            if (T == f32) return scalar.logf(x);
            return scalar.logd(x);
        }

        inline fn log10(x: T) T {
            if (T == f32) return scalar.log10f(x);
            return scalar.log10d(x);
        }

        inline fn pow(x: T, y: T) T {
            if (T == f32) return scalar.powf(x, y);
            return scalar.powd(x, y);
        }

        inline fn sqrt(x: T) T {
            if (T == f32) return scalar.sqrtf(x);
            return scalar.sqrtd(x);
        }

        inline fn hypot(x: T, y: T) T {
            if (T == f32) return scalar.hypotf(x, y);
            return scalar.hypotd(x, y);
        }

        inline fn atan2(y: T, x: T) T {
            if (T == f32) return scalar.atan2f(y, x).value;
            return scalar.atan2d(y, x).value;
        }

        inline fn sin(x: T) T {
            if (T == f32) return scalar.sinf(x);
            return scalar.sind(x);
        }

        inline fn cos(x: T) T {
            if (T == f32) return scalar.cosf(x);
            return scalar.cosd(x);
        }

        inline fn tan(x: T) T {
            if (T == f32) return scalar.tanf(x);
            return scalar.tand(x);
        }

        inline fn sinh(x: T) T {
            if (T == f32) return scalar.sinhf(x);
            return scalar.sinhd(x);
        }

        inline fn cosh(x: T) T {
            if (T == f32) return scalar.coshf(x);
            return scalar.coshd(x);
        }

        inline fn tanh(x: T) T {
            if (T == f32) return scalar.tanhf(x);
            return scalar.tanhd(x);
        }

        inline fn ldexp(x: T, exponent: i32) T {
            if (T == f32) return scalar.scalbnf(x, exponent);
            return scalar.scalbnd(x, exponent);
        }

        /// Splits a positive finite magnitude into a fraction in `[1/2, 1)`
        /// and the power of two that the fraction was scaled by.
        inline fn frexp(x: T) scalar.FrexpResult(T) {
            if (T == f32) return scalar.frexpf(x);
            return scalar.frexpd(x);
        }
    };
}

// ---------------------------------------------------------------------------
// Classification and NaN ordering
// ---------------------------------------------------------------------------

inline fn isNaN(comptime T: type, value: T) bool {
    if (T == f32) {
        const bits: u32 = @bitCast(value);
        return (bits & 0x7f80_0000) == 0x7f80_0000 and (bits & 0x007f_ffff) != 0;
    }
    const bits: u64 = @bitCast(value);
    return (bits & 0x7ff0_0000_0000_0000) == 0x7ff0_0000_0000_0000 and
        (bits & 0x000f_ffff_ffff_ffff) != 0;
}

inline fn isInfinity(comptime T: type, value: T) bool {
    if (T == f32) {
        const bits: u32 = @bitCast(value);
        return (bits & 0x7fff_ffff) == 0x7f80_0000;
    }
    const bits: u64 = @bitCast(value);
    return (bits & 0x7fff_ffff_ffff_ffff) == 0x7ff0_0000_0000_0000;
}

/// `value` with the quiet bit set. The caller must pass a NaN of `T`, so the
/// sign and the diagnostic payload survive and only the quiet bit is added.
inline fn quieted(comptime T: type, value: T) T {
    if (T == f32) {
        const bits: u32 = @bitCast(value);
        return @bitCast(bits | 0x0040_0000);
    }
    const bits: u64 = @bitCast(value);
    return @bitCast(bits | 0x0008_0000_0000_0000);
}

/// The first NaN among `a` and `b`, quieted, or `null` when both are numbers.
///
/// The scan fixes the winner, because the IEEE NaN payload a host arithmetic
/// would pick is not specified.
inline fn firstNaN2(comptime T: type, a: T, b: T) ?T {
    if (isNaN(T, a)) return quieted(T, a);
    if (isNaN(T, b)) return quieted(T, b);
    return null;
}

/// The first NaN among four operands, quieted, in the order given.
inline fn firstNaN4(comptime T: type, a: T, b: T, c: T, d: T) ?T {
    if (isNaN(T, a)) return quieted(T, a);
    if (isNaN(T, b)) return quieted(T, b);
    if (isNaN(T, c)) return quieted(T, c);
    if (isNaN(T, d)) return quieted(T, d);
    return null;
}

/// `a * b`, except that an infinite factor paired with a zero factor gives the
/// infinity rather than the IEEE NaN.
///
/// A complex decomposition scales every component by the same factor, so a
/// zero factor that belongs to one component must not destroy the infinity
/// that belongs to the other.
inline fn product(comptime T: type, a: T, b: T) T {
    if (isInfinity(T, a) and b == 0) return a;
    if (isInfinity(T, b) and a == 0) return a * withSign(T, 1, b);
    return a * b;
}

/// The magnitude `value` with the sign of `sign`, from the accepted IEEE
/// kernel rather than a fresh sign-bit shuffle.
inline fn withSign(comptime T: type, value: T, sign: T) T {
    if (T == f32) return ieee.copySignf(value, sign);
    return ieee.copySignd(value, sign);
}

/// `a / b` through the divide that keeps a finite binary64 underflow.
///
/// The binary32 kernel divides natively, and so does every binary64 target
/// with a native FMA. Elsewhere the target lowers the binary64 divide to a
/// software helper that flushes a subnormal quotient to a signed zero, so the
/// reductions below would lose a quotient they are written to survive. Every
/// quotient in this module can reach that range, so they all take the shared
/// exact path, which is the same single divide for a normal, an infinite or a
/// NaN result and so leaves every other result bit alone.
inline fn divide(comptime T: type, a: T, b: T) T {
    if (T == f32) return a / b;
    return scalar_div.div64(a, b);
}

// ---------------------------------------------------------------------------
// Accessors
// ---------------------------------------------------------------------------

/// The real part of `value`.
pub fn realf(value: Complexf) f32 {
    return value.re;
}

/// The real part of `value`.
pub fn reald(value: Complexd) f64 {
    return value.re;
}

/// The imaginary part of `value`.
pub fn imagf(value: Complexf) f32 {
    return value.im;
}

/// The imaginary part of `value`.
pub fn imagd(value: Complexd) f64 {
    return value.im;
}

/// `value` with its imaginary part negated.
pub fn conjf(value: Complexf) Complexf {
    return .{ .re = value.re, .im = -value.im };
}

/// `value` with its imaginary part negated.
pub fn conjd(value: Complexd) Complexd {
    return .{ .re = value.re, .im = -value.im };
}

/// `value` projected onto the real axis.
///
/// An infinite component is the case that separates this from a plain
/// assignment: `inf - i * inf` is a NaN, while the projection of an infinite
/// value is `+inf + i0`.
fn projValue(comptime T: type, value: ComplexOf(T)) ComplexOf(T) {
    if (firstNaN2(T, value.re, value.im)) |nan| return .{ .re = nan, .im = nan };
    if (isInfinity(T, value.re) or isInfinity(T, value.im)) {
        return .{ .re = Format(T).infinity, .im = 0 };
    }
    return .{ .re = value.re, .im = 0 };
}

/// `value` projected onto the real axis.
pub fn projf(value: Complexf) Complexf {
    return projValue(f32, value);
}

/// `value` projected onto the real axis.
pub fn projd(value: Complexd) Complexd {
    return projValue(f64, value);
}

// ---------------------------------------------------------------------------
// Magnitude and argument
// ---------------------------------------------------------------------------

/// `hypot(x, y)` without the intermediate overflow or underflow that a squared
/// form suffers at extreme magnitudes.
inline fn magnitude(comptime T: type, x: T, y: T) T {
    const F = Format(T);
    const a = @abs(x);
    const b = @abs(y);
    const m = @max(a, b);
    if (m == 0) return 0;
    // Below the limit the true magnitude cannot round to infinity, and a
    // subnormal input stays subnormal rather than underflowing to zero, so the
    // kernel's own correctly rounded reduction is both safe and the most
    // accurate answer. Above it the magnitude can round to infinity, so the
    // ratio is reduced to `[0, 1]` first and the scale is applied once.
    // A subnormal ratio is not observable through the kernel, which answers one
    // for any argument this small, so the reduced divide keeps the target's own.
    if (m <= F.magnitude_limit) return F.hypot(a, b);
    return m * F.hypot(1, @min(a, b) / m);
}

fn absValue(comptime T: type, value: ComplexOf(T)) T {
    const F = Format(T);
    // An infinite component wins over a NaN one, because the value is not
    // merely unknown, it is unbounded. A NaN component alone has no magnitude
    // to report, so the result is the canonical quiet NaN.
    if (isInfinity(T, value.re) or isInfinity(T, value.im)) return F.infinity;
    if (isNaN(T, value.re) or isNaN(T, value.im)) return F.nan;
    return magnitude(T, value.re, value.im);
}

/// The magnitude of `value`.
pub fn absf(value: Complexf) f32 {
    return absValue(f32, value);
}

/// The magnitude of `value`.
pub fn absd(value: Complexd) f64 {
    return absValue(f64, value);
}

/// The principal argument of `value` in `(-pi, +pi]`.
///
/// The quadrant, the signed-zero cases, and the pi constant all come from the
/// accepted `atan2` kernel, so this is a routing decision and not a second
/// argument implementation.
fn argValue(comptime T: type, value: ComplexOf(T)) T {
    return Format(T).atan2(value.im, value.re);
}

/// The principal argument of `value` in `(-pi, +pi]`.
pub fn argf(value: Complexf) f32 {
    return argValue(f32, value);
}

/// The principal argument of `value` in `(-pi, +pi]`.
pub fn argd(value: Complexd) f64 {
    return argValue(f64, value);
}

// ---------------------------------------------------------------------------
// Addition, subtraction, multiplication, division
// ---------------------------------------------------------------------------

/// One component of an addition, with a fixed NaN winner.
///
/// Two non-NaN operands can still be invalid, as in `inf + -inf`, and the host
/// indefinite NaN that produces carries no sign the specification can pin, so
/// an invalid result is the positive canonical quiet NaN.
inline fn addTerm(comptime T: type, a: T, b: T) T {
    if (isNaN(T, a)) return quieted(T, a);
    if (isNaN(T, b)) return quieted(T, b);
    const sum = a + b;
    if (isNaN(T, sum)) return Format(T).nan;
    return sum;
}

/// One component of a subtraction, with a fixed NaN winner and the same
/// canonical quiet NaN for an invalid result as `addTerm`.
inline fn subTerm(comptime T: type, a: T, b: T) T {
    if (isNaN(T, a)) return quieted(T, a);
    if (isNaN(T, b)) return quieted(T, b);
    const difference = a - b;
    if (isNaN(T, difference)) return Format(T).nan;
    return difference;
}

/// Componentwise addition without the NaN scan the public entry points do.
inline fn plus(comptime T: type, a: ComplexOf(T), b: ComplexOf(T)) ComplexOf(T) {
    return .{ .re = a.re + b.re, .im = a.im + b.im };
}

/// Componentwise subtraction without the NaN scan the public entry points do.
inline fn minus(comptime T: type, a: ComplexOf(T), b: ComplexOf(T)) ComplexOf(T) {
    return .{ .re = a.re - b.re, .im = a.im - b.im };
}

/// The naive four-multiply product.
///
/// The scaled and Smith products would be more accurate for a cancelling
/// operand pair, but they change the result bits of every case that the naive
/// form already gets right, and there is no oracle here to justify that. An
/// invalid product such as `inf * 0` is reported as the canonical quiet NaN,
/// because the host indefinite NaN that the arithmetic would produce carries
/// no sign the specification can pin.
fn mulValue(comptime T: type, a: ComplexOf(T), b: ComplexOf(T)) ComplexOf(T) {
    if (firstNaN4(T, a.re, a.im, b.re, b.im)) |nan| return .{ .re = nan, .im = nan };
    const re = a.re * b.re - a.im * b.im;
    const im = a.re * b.im + a.im * b.re;
    if (isNaN(T, re) or isNaN(T, im)) {
        const nan = Format(T).nan;
        return .{ .re = nan, .im = nan };
    }
    return .{ .re = re, .im = im };
}

/// `a + b` componentwise.
pub fn addf(a: Complexf, b: Complexf) Complexf {
    return .{
        .re = addTerm(f32, a.re, b.re),
        .im = addTerm(f32, a.im, b.im),
    };
}

/// `a + b` componentwise.
pub fn addd(a: Complexd, b: Complexd) Complexd {
    return .{
        .re = addTerm(f64, a.re, b.re),
        .im = addTerm(f64, a.im, b.im),
    };
}

/// `a - b` componentwise.
pub fn subf(a: Complexf, b: Complexf) Complexf {
    return .{
        .re = subTerm(f32, a.re, b.re),
        .im = subTerm(f32, a.im, b.im),
    };
}

/// `a - b` componentwise.
pub fn subd(a: Complexd, b: Complexd) Complexd {
    return .{
        .re = subTerm(f64, a.re, b.re),
        .im = subTerm(f64, a.im, b.im),
    };
}

/// `a * b` with the naive four-multiply form.
pub fn mulf(a: Complexf, b: Complexf) Complexf {
    return mulValue(f32, a, b);
}

/// `a * b` with the naive four-multiply form.
pub fn muld(a: Complexd, b: Complexd) Complexd {
    return mulValue(f64, a, b);
}

/// The signed zero of `x` divided by a nonzero infinite magnitude, which is
/// the sign of `x` reversed when the divisor is negative. A zero numerator
/// gives `+0` because `0 / y` is positive for every infinite `y`.
inline fn zeroQuotient(comptime T: type, x: T, negate: bool) T {
    if (x == 0) return 0;
    return withSign(T, 0, if (negate) -x else x);
}

/// `a / b` from the conjugate product, with the denominator scaled so that
/// neither its square nor the numerator can overflow or underflow.
///
/// The direct form divides by `b.re^2 + b.im^2`, which overflows for
/// operands past sqrt(max) and underflows to zero for operands below
/// sqrt(min_normal). Outside that window the divisor is divided by its
/// larger-magnitude component first and the scale is applied to the result
/// once, which is the same value with the same range.
fn divValue(comptime T: type, a: ComplexOf(T), b: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN4(T, a.re, a.im, b.re, b.im)) |nan| return .{ .re = nan, .im = nan };

    if (isInfinity(T, b.re) or isInfinity(T, b.im)) {
        // An infinite divisor inverts to a zero magnitude, so the quotient
        // keeps only the signs of the leading term of the reciprocal. An
        // infinite numerator against an infinite divisor, and a divisor whose
        // two components are both infinite, have no direction at all.
        if (isInfinity(T, a.re) or isInfinity(T, a.im)) return .{ .re = F.nan, .im = F.nan };
        if (isInfinity(T, b.re) and isInfinity(T, b.im)) return .{ .re = F.nan, .im = F.nan };
        if (isInfinity(T, b.re)) {
            // 1/b is 1/b.re to a relative error below every binary format, so
            // the quotient is a componentwise division by b.re.
            return .{
                .re = zeroQuotient(T, a.re, b.re < 0),
                .im = zeroQuotient(T, a.im, b.re < 0),
            };
        }
        // 1/b is -i/b.im, which turns the numerator a quarter turn before the
        // division.
        return .{
            .re = zeroQuotient(T, a.im, b.im < 0),
            .im = zeroQuotient(T, a.re, b.im > 0),
        };
    }

    const scale = @max(@abs(b.re), @abs(b.im));
    if (scale == 0) {
        // Zero over zero is undefined, and a nonzero numerator over a zero
        // divisor has infinite magnitude and undefined direction, so only the
        // real part survives.
        if (a.re == 0 and a.im == 0) return .{ .re = F.nan, .im = F.nan };
        if (isInfinity(T, a.re) or isInfinity(T, a.im)) return .{ .re = F.nan, .im = F.nan };
        return .{
            .re = withSign(T, F.infinity, a.re * b.re + a.im * b.im),
            .im = F.nan,
        };
    }

    // `d <= 2 * scale^2`, so the bound is sqrt(max / 2) rather than sqrt(max).
    if (scale >= F.sqrt_min_normal and scale <= F.magnitude_limit) {
        const d = b.re * b.re + b.im * b.im;
        return .{
            .re = divide(T, product(T, a.re, b.re) + product(T, a.im, b.im), d),
            .im = divide(T, product(T, a.im, b.re) - product(T, a.re, b.im), d),
        };
    }

    // With `b = scale * (b.re/scale, b.im/scale)` the conjugate product is
    // `scale^2 * d * (numerator)`, so the scale must come off the quotient as a
    // division. Multiplying by it here would underflow a tiny numerator against
    // a tiny divisor and overflow a large one against a large divisor.
    const divisor_re = divide(T, b.re, scale);
    const divisor_im = divide(T, b.im, scale);
    const d = divisor_re * divisor_re + divisor_im * divisor_im;
    const numerator_re = product(T, a.re, divisor_re) + product(T, a.im, divisor_im);
    const numerator_im = product(T, a.im, divisor_re) - product(T, a.re, divisor_im);
    return .{
        .re = divide(T, divide(T, numerator_re, d), scale),
        .im = divide(T, divide(T, numerator_im, d), scale),
    };
}

/// `a / b`.
pub fn divf(a: Complexf, b: Complexf) Complexf {
    return divValue(f32, a, b);
}

/// `a / b`.
pub fn divd(a: Complexd, b: Complexd) Complexd {
    return divValue(f64, a, b);
}

// ---------------------------------------------------------------------------
// Exponential, logarithm, square root, power
// ---------------------------------------------------------------------------

fn expValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN2(T, z.re, z.im)) |nan| return .{ .re = nan, .im = nan };
    if (isInfinity(T, z.im)) return .{ .re = F.nan, .im = F.nan };
    // exp(x + i0) is exactly e^x, so the imaginary part survives as the input
    // zero and a real overflow does not become a NaN through inf * 0.
    if (z.im == 0) return .{ .re = F.exp(z.re), .im = z.im };
    if (isInfinity(T, z.re)) {
        const scale: T = if (z.re > 0) F.infinity else 0;
        return .{ .re = scale, .im = F.nan };
    }
    const scale = F.exp(z.re);
    return .{
        .re = scale * F.cos(z.im),
        .im = scale * F.sin(z.im),
    };
}

/// `value` on the base-10 scale when `base_ten` is set, and unchanged
/// otherwise. The division keeps `log10(z) = log(z) / ln(10)` to one rounding.
inline fn toBaseTen(comptime T: type, value: T, comptime base_ten: bool) T {
    if (!base_ten) return value;
    return divide(T, value, Format(T).ln_10);
}

/// The answer for an input the polar reduction cannot handle, or `null` when
/// the magnitude and the argument are the answer.
///
/// Every special value of the real part is an infinity or a NaN, which the
/// natural and the base-10 logarithm share, so only the argument is scaled.
fn polarSpecial(comptime T: type, z: ComplexOf(T), comptime base_ten: bool) ?ComplexOf(T) {
    const F = Format(T);
    if (firstNaN2(T, z.re, z.im)) |nan| return .{ .re = nan, .im = nan };
    if (isInfinity(T, z.re) and isInfinity(T, z.im)) {
        return .{ .re = F.infinity, .im = F.nan };
    }
    if (isInfinity(T, z.re)) {
        return .{ .re = F.infinity, .im = toBaseTen(T, withSign(T, 0, z.im), base_ten) };
    }
    if (isInfinity(T, z.im)) {
        return .{ .re = F.infinity, .im = toBaseTen(T, withSign(T, F.half_pi, z.im), base_ten) };
    }
    // The origin has an infinite logarithm and a zero argument angle, not the
    // signed pi that atan2 would return for a negative zero real part.
    if (z.re == 0 and z.im == 0) {
        return .{ .re = -F.infinity, .im = toBaseTen(T, z.im, base_ten) };
    }
    return null;
}

/// `log(z) = ln|z| + i arg(z)`.
fn logValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (polarSpecial(T, z, false)) |special| return special;
    return .{
        .re = F.log(magnitude(T, z.re, z.im)),
        .im = F.atan2(z.im, z.re),
    };
}

/// `log10(z) = log(z) / ln(10)`, with the real part taken from the base-10
/// kernel on the magnitude itself, which is more accurate than dividing the
/// natural logarithm by ln(10).
fn log10Value(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (polarSpecial(T, z, true)) |special| return special;
    return .{
        .re = F.log10(magnitude(T, z.re, z.im)),
        .im = toBaseTen(T, F.atan2(z.im, z.re), true),
    };
}

/// The `sqrtValue` reduction for magnitudes where the direct form over- or
/// underflows: scale by an exact power of two, take the root, then scale back.
///
/// Halving the exponent first and folding the odd factor back in as sqrt(2)
/// keeps both rescalings exact powers of two.
fn sqrtScaled(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    const split = F.frexp(@max(@abs(z.re), @abs(z.im)));
    const x = F.ldexp(z.re, -split.exponent);
    const y = F.ldexp(z.im, -split.exponent);
    const r = F.hypot(x, y);
    var re: T = undefined;
    var im: T = undefined;
    if (x >= 0) {
        // The larger half is the square root of the real part, and the other
        // one comes from the identity 2 * u * v = y, which never cancels.
        re = F.sqrt(0.5 * r + 0.5 * x);
        im = divide(T, y, 2 * re);
    } else {
        im = F.sqrt(0.5 * r - 0.5 * x);
        re = divide(T, @abs(y), 2 * im);
    }
    const half = @divFloor(split.exponent, 2);
    if (split.exponent - 2 * half == 1) {
        re *= F.sqrt_two;
        im *= F.sqrt_two;
    }
    return .{ .re = F.ldexp(re, half), .im = withSign(T, F.ldexp(im, half), z.im) };
}

fn sqrtValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN2(T, z.re, z.im)) |nan| return .{ .re = nan, .im = nan };
    if (isInfinity(T, z.re) and isInfinity(T, z.im)) {
        return .{ .re = F.infinity, .im = F.infinity };
    }
    if (isInfinity(T, z.re)) {
        if (z.im == 0) return .{ .re = F.infinity, .im = z.im };
        return .{ .re = F.infinity, .im = withSign(T, F.infinity, z.im) };
    }
    if (isInfinity(T, z.im)) return .{ .re = F.infinity, .im = F.infinity };
    // The origin is the only value whose real part is a zero, and it keeps
    // the signs of the input.
    if (z.re == 0 and z.im == 0) return z;
    const x = @abs(z.re);
    const y = @abs(z.im);
    const m = @max(x, y);
    // Inside this window both halves of the reduction are normal and their
    // sum is finite, so the direct form is both safe and the most accurate.
    if (m >= F.sqrt_min_normal and m <= F.magnitude_limit) {
        // Both halves use the magnitude of the real part: the real part of the
        // root squares to (r + a) / 2 and the imaginary part to (r - a) / 2,
        // so for a negative a the imaginary half is the larger one.
        const r = F.hypot(x, y);
        if (z.re >= 0) {
            const re = F.sqrt(0.5 * r + 0.5 * z.re);
            return .{ .re = re, .im = divide(T, z.im, 2 * re) };
        }
        const im = F.sqrt(0.5 * r - 0.5 * z.re);
        return .{ .re = divide(T, y, 2 * im), .im = withSign(T, im, z.im) };
    }
    return sqrtScaled(T, z);
}

/// `exp(z) = e^z.re (cos z.im + i sin z.im)`.
pub fn expf(z: Complexf) Complexf {
    return expValue(f32, z);
}

/// `exp(z) = e^z.re (cos z.im + i sin z.im)`.
pub fn expd(z: Complexd) Complexd {
    return expValue(f64, z);
}

/// `log(z) = ln|z| + i arg(z)`.
pub fn logf(z: Complexf) Complexf {
    return logValue(f32, z);
}

/// `log(z) = ln|z| + i arg(z)`.
pub fn logd(z: Complexd) Complexd {
    return logValue(f64, z);
}

/// `log10(z) = log(z) / ln(10)`.
pub fn log10f(z: Complexf) Complexf {
    return log10Value(f32, z);
}

/// `log10(z) = log(z) / ln(10)`.
pub fn log10d(z: Complexd) Complexd {
    return log10Value(f64, z);
}

/// The principal square root, with a nonnegative real part and an imaginary
/// part signed by the input.
pub fn sqrtf(z: Complexf) Complexf {
    return sqrtValue(f32, z);
}

/// The principal square root, with a nonnegative real part and an imaginary
/// part signed by the input.
pub fn sqrtd(z: Complexd) Complexd {
    return sqrtValue(f64, z);
}

/// `pow(x, y)`, the principal value of `exp(y log x)`.
///
/// Three cases bypass that expression because it cannot reproduce them.
/// A zero exponent gives one whatever the base, matching `pow(0, 0) = 1`. A
/// zero base has no logarithm, so a positive real exponent gives zero and any
/// other exponent gives an infinite magnitude with an undefined direction. A
/// negative real base with an exact real integer exponent is real by
/// definition, while `exp(y log x)` would carry the argument reduction error
/// of `pi` in its imaginary part instead of an exact zero.
fn powValue(comptime T: type, x: ComplexOf(T), y: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN4(T, x.re, x.im, y.re, y.im)) |nan| return .{ .re = nan, .im = nan };
    if (y.re == 0 and y.im == 0) return .{ .re = 1, .im = y.im };
    if (x.re == 0 and x.im == 0) {
        if (y.im != 0 or y.re < 0) return .{ .re = F.infinity, .im = F.nan };
        return .{ .re = 0, .im = 0 };
    }
    // An infinite base or exponent leaves the expression without a value, and
    // the host product of an infinity with a zero would report its own
    // indefinite NaN, so the canonical quiet NaN is returned instead.
    if (isInfinity(T, x.re) or isInfinity(T, x.im) or
        isInfinity(T, y.re) or isInfinity(T, y.im))
    {
        return .{ .re = F.nan, .im = F.nan };
    }
    if (x.im == 0 and x.re < 0 and y.im == 0 and y.re == @trunc(y.re)) {
        const scale = F.pow(-x.re, y.re);
        return .{ .re = if (@mod(y.re, 2) == 0) scale else -scale, .im = 0 };
    }
    return expValue(T, mulValue(T, y, logValue(T, x)));
}

/// `pow(x, y)`, the principal value of `exp(y log x)`.
pub fn powf(x: Complexf, y: Complexf) Complexf {
    return powValue(f32, x, y);
}

/// `pow(x, y)`, the principal value of `exp(y log x)`.
pub fn powd(x: Complexd, y: Complexd) Complexd {
    return powValue(f64, x, y);
}

// ---------------------------------------------------------------------------
// Trigonometric families
// ---------------------------------------------------------------------------

fn sinValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN2(T, z.re, z.im)) |nan| return .{ .re = nan, .im = nan };
    if (z.im == 0) return .{ .re = F.sin(z.re), .im = z.im };
    if (isInfinity(T, z.im)) {
        return .{ .re = F.nan, .im = withSign(T, F.infinity, z.im) };
    }
    if (isInfinity(T, z.re)) return .{ .re = F.nan, .im = F.nan };
    // sin(a + ib) = sin a cosh b + i cos a sinh b. The identity is exact, so
    // the two components share one pair of hyperbolic factors.
    const cosh_b = F.cosh(z.im);
    const sinh_b = F.sinh(z.im);
    return .{
        .re = product(T, F.sin(z.re), cosh_b),
        .im = product(T, F.cos(z.re), sinh_b),
    };
}

fn cosValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN2(T, z.re, z.im)) |nan| return .{ .re = nan, .im = nan };
    if (z.im == 0) return .{ .re = F.cos(z.re), .im = z.im };
    if (isInfinity(T, z.im)) {
        return .{ .re = F.nan, .im = -withSign(T, F.infinity, z.im) };
    }
    if (isInfinity(T, z.re)) return .{ .re = F.nan, .im = F.nan };
    // cos(a + ib) = cos a cosh b - i sin a sinh b.
    const cosh_b = F.cosh(z.im);
    const sinh_b = F.sinh(z.im);
    return .{
        .re = product(T, F.cos(z.re), cosh_b),
        .im = -product(T, F.sin(z.re), sinh_b),
    };
}

/// `tan(a + ib) = (sin 2a + i sinh 2b) / (cos 2a + cosh 2b)`.
///
/// The double angles come from the half angles, so no factor of two can
/// overflow an otherwise finite argument. The denominator is a real number at
/// least zero, and it vanishes only when `b` is zero, which the `z.im == 0`
/// path already serves, so the pole is not reachable from a binary32 or
/// binary64 argument.
fn tanValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN2(T, z.re, z.im)) |nan| return .{ .re = nan, .im = nan };
    if (z.im == 0) return .{ .re = F.tan(z.re), .im = z.im };
    if (isInfinity(T, z.im)) {
        return .{ .re = F.nan, .im = withSign(T, F.infinity, z.im) };
    }
    if (isInfinity(T, z.re)) return .{ .re = F.nan, .im = F.nan };
    const sin_a = F.sin(z.re);
    const cos_a = F.cos(z.re);
    const cosh_b = F.cosh(z.im);
    const sinh_b = F.sinh(z.im);
    const numerator_re = 2 * sin_a * cos_a;
    const denominator = cos_a * cos_a - sin_a * sin_a + cosh_b * cosh_b + sinh_b * sinh_b;
    if (isInfinity(T, denominator)) {
        // cosh(2b) overflows long before the bounded quotient
        // sinh(2b) / cosh(2b) = tanh(2b) loses any accuracy, so the overflow
        // case takes the quotient from tanh and the real part from the sign
        // of the vanishing numerator.
        return .{ .re = withSign(T, 0, numerator_re), .im = F.tanh(2 * z.im) };
    }
    return .{
        .re = divide(T, numerator_re, denominator),
        .im = divide(T, 2 * sinh_b * cosh_b, denominator),
    };
}

/// `sin(z)`.
pub fn sinf(z: Complexf) Complexf {
    return sinValue(f32, z);
}

/// `sin(z)`.
pub fn sind(z: Complexd) Complexd {
    return sinValue(f64, z);
}

/// `cos(z)`.
pub fn cosf(z: Complexf) Complexf {
    return cosValue(f32, z);
}

/// `cos(z)`.
pub fn cosd(z: Complexd) Complexd {
    return cosValue(f64, z);
}

/// `tan(z)`.
pub fn tanf(z: Complexf) Complexf {
    return tanValue(f32, z);
}

/// `tan(z)`.
pub fn tand(z: Complexd) Complexd {
    return tanValue(f64, z);
}

// ---------------------------------------------------------------------------
// Hyperbolic families
// ---------------------------------------------------------------------------

fn sinhValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN2(T, z.re, z.im)) |nan| return .{ .re = nan, .im = nan };
    if (z.im == 0) return .{ .re = F.sinh(z.re), .im = z.im };
    if (isInfinity(T, z.re) or isInfinity(T, z.im)) return .{ .re = F.nan, .im = F.nan };
    // sinh(a + ib) = sinh a cos b + i cosh a sin b.
    const cos_b = F.cos(z.im);
    const sin_b = F.sin(z.im);
    return .{
        .re = product(T, F.sinh(z.re), cos_b),
        .im = product(T, F.cosh(z.re), sin_b),
    };
}

fn coshValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN2(T, z.re, z.im)) |nan| return .{ .re = nan, .im = nan };
    if (z.im == 0) return .{ .re = F.cosh(z.re), .im = z.im };
    if (isInfinity(T, z.re) or isInfinity(T, z.im)) return .{ .re = F.nan, .im = F.nan };
    // cosh(a + ib) = cosh a cos b + i sinh a sin b.
    const cos_b = F.cos(z.im);
    const sin_b = F.sin(z.im);
    return .{
        .re = product(T, F.cosh(z.re), cos_b),
        .im = product(T, F.sinh(z.re), sin_b),
    };
}

/// `tanh(a + ib) = (sinh 2a + i sin 2b) / (cosh 2a + cos 2b)`.
///
/// The double angles come from the half angles, so no factor of two can
/// overflow an otherwise finite argument.
fn tanhValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (firstNaN2(T, z.re, z.im)) |nan| return .{ .re = nan, .im = nan };
    if (z.im == 0) return .{ .re = F.tanh(z.re), .im = z.im };
    // The limit along either direction is a real one.
    if (isInfinity(T, z.re)) return .{ .re = withSign(T, 1, z.re), .im = 0 };
    if (isInfinity(T, z.im)) return .{ .re = F.nan, .im = F.nan };
    const cosh_a = F.cosh(z.re);
    const sinh_a = F.sinh(z.re);
    const cos_b = F.cos(z.im);
    const sin_b = F.sin(z.im);
    const denominator = cosh_a * cosh_a + sinh_a * sinh_a + cos_b * cos_b - sin_b * sin_b;
    if (isInfinity(T, denominator)) return .{ .re = withSign(T, 1, z.re), .im = 0 };
    // The denominator is cosh(2a) + cos(2b), which vanishes at a = 0 and
    // b = +-pi/2. There the real numerator is a zero and the imaginary one is
    // not, so the pole is the canonical quiet NaN and a signed infinity rather
    // than a zero over zero.
    if (denominator == 0) {
        return .{
            .re = F.nan,
            .im = withSign(T, F.infinity, 2 * sin_b * cos_b),
        };
    }
    return .{
        .re = divide(T, 2 * sinh_a * cosh_a, denominator),
        .im = divide(T, 2 * sin_b * cos_b, denominator),
    };
}

/// `sinh(z)`.
pub fn sinhf(z: Complexf) Complexf {
    return sinhValue(f32, z);
}

/// `sinh(z)`.
pub fn sinhd(z: Complexd) Complexd {
    return sinhValue(f64, z);
}

/// `cosh(z)`.
pub fn coshf(z: Complexf) Complexf {
    return coshValue(f32, z);
}

/// `cosh(z)`.
pub fn coshd(z: Complexd) Complexd {
    return coshValue(f64, z);
}

/// `tanh(z)`.
pub fn tanhf(z: Complexf) Complexf {
    return tanhValue(f32, z);
}

/// `tanh(z)`.
pub fn tanhd(z: Complexd) Complexd {
    return tanhValue(f64, z);
}

// ---------------------------------------------------------------------------
// Inverse trigonometric and inverse hyperbolic families
// ---------------------------------------------------------------------------

/// `(canonical, canonical)` for an input a C mathematical pole, a branch
/// point at infinity, or a NaN leaves the decomposition without a value the
/// standard pins, and the reported value must not come from the host.
inline fn rejectInfinity(comptime T: type, z: ComplexOf(T)) ?@TypeOf(Format(T).nan) {
    if (firstNaN2(T, z.re, z.im)) |nan| return nan;
    if (isInfinity(T, z.re) or isInfinity(T, z.im)) return Format(T).nan;
    return null;
}

/// `asin(z) = -i log(i z + sqrt(1 - z^2))`.
///
/// The square root is taken as `sqrt(1 - z) sqrt(1 + z)`, which is the
/// principal root of `1 - z^2` on the same branch and keeps every
/// intermediate at the magnitude of `z` instead of its square.
fn asinValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const one: ComplexOf(T) = .{ .re = 1, .im = 0 };
    if (rejectInfinity(T, z)) |nan| return .{ .re = nan, .im = nan };
    const root = mulValue(
        T,
        sqrtValue(T, minus(T, one, z)),
        sqrtValue(T, plus(T, one, z)),
    );
    const w = logValue(T, plus(T, root, .{ .re = -z.im, .im = z.re }));
    return .{ .re = w.im, .im = -w.re };
}

/// `acos(z) = pi/2 - asin(z)`, which is the branch that inverts `cos`.
///
/// The invalid case is rejected before the subtraction, because negating the
/// canonical quiet NaN would report a negative one.
fn acosValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const F = Format(T);
    if (rejectInfinity(T, z)) |nan| return .{ .re = nan, .im = nan };
    const w = asinValue(T, z);
    return .{ .re = F.half_pi - w.re, .im = -w.im };
}

/// `atan(z) = (i/2) (log(1 - i z) - log(1 + i z))`.
fn atanValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const one: ComplexOf(T) = .{ .re = 1, .im = 0 };
    if (rejectInfinity(T, z)) |nan| return .{ .re = nan, .im = nan };
    const iz: ComplexOf(T) = .{ .re = -z.im, .im = z.re };
    const w = minus(
        T,
        logValue(T, minus(T, one, iz)),
        logValue(T, plus(T, one, iz)),
    );
    // The factor i/2 is applied as a naive product, so the origin keeps the
    // positive zero that the real atan(0) is.
    return mulValue(T, .{ .re = 0, .im = 0.5 }, w);
}

/// `asinh(z) = log(z + sqrt(z^2 + 1))`.
///
/// The root is taken as `sqrt(z + i) sqrt(z - i)`, which is the principal
/// root of `z^2 + 1` on the same branch and keeps every intermediate at the
/// magnitude of `z` instead of its square.
fn asinhValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const unit_imag: ComplexOf(T) = .{ .re = 0, .im = 1 };
    if (rejectInfinity(T, z)) |nan| return .{ .re = nan, .im = nan };
    const root = mulValue(
        T,
        sqrtValue(T, plus(T, z, unit_imag)),
        sqrtValue(T, minus(T, z, unit_imag)),
    );
    return logValue(T, plus(T, z, root));
}

/// `acosh(z) = log(z + sqrt(z + 1) sqrt(z - 1))`.
///
/// The two roots are taken separately, which is the branch whose real part is
/// nonnegative: the principal root of `z^2 - 1` would be negative for a real
/// input above `+1`.
fn acoshValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const one: ComplexOf(T) = .{ .re = 1, .im = 0 };
    if (rejectInfinity(T, z)) |nan| return .{ .re = nan, .im = nan };
    const root = mulValue(
        T,
        sqrtValue(T, plus(T, z, one)),
        sqrtValue(T, minus(T, z, one)),
    );
    return logValue(T, plus(T, z, root));
}

/// `atanh(z) = (log(1 + z) - log(1 - z)) / 2`.
fn atanhValue(comptime T: type, z: ComplexOf(T)) ComplexOf(T) {
    const one: ComplexOf(T) = .{ .re = 1, .im = 0 };
    if (rejectInfinity(T, z)) |nan| return .{ .re = nan, .im = nan };
    const w = minus(
        T,
        logValue(T, plus(T, one, z)),
        logValue(T, minus(T, one, z)),
    );
    return .{ .re = 0.5 * w.re, .im = 0.5 * w.im };
}

/// `asin(z)`, with a real part in `[-pi/2, +pi/2]`.
pub fn asinf(z: Complexf) Complexf {
    return asinValue(f32, z);
}

/// `asin(z)`, with a real part in `[-pi/2, +pi/2]`.
pub fn asind(z: Complexd) Complexd {
    return asinValue(f64, z);
}

/// `acos(z)`, with a real part in `[0, +pi]`.
pub fn acosf(z: Complexf) Complexf {
    return acosValue(f32, z);
}

/// `acos(z)`, with a real part in `[0, +pi]`.
pub fn acosd(z: Complexd) Complexd {
    return acosValue(f64, z);
}

/// `atan(z)`, with a real part in `(-pi/2, +pi/2)`.
pub fn atanf(z: Complexf) Complexf {
    return atanValue(f32, z);
}

/// `atan(z)`, with a real part in `(-pi/2, +pi/2)`.
pub fn atand(z: Complexd) Complexd {
    return atanValue(f64, z);
}

/// `asinh(z)`, with an imaginary part in `[-pi/2, +pi/2]`.
pub fn asinhf(z: Complexf) Complexf {
    return asinhValue(f32, z);
}

/// `asinh(z)`, with an imaginary part in `[-pi/2, +pi/2]`.
pub fn asinhd(z: Complexd) Complexd {
    return asinhValue(f64, z);
}

/// `acosh(z)`, with a nonnegative real part.
pub fn acoshf(z: Complexf) Complexf {
    return acoshValue(f32, z);
}

/// `acosh(z)`, with a nonnegative real part.
pub fn acoshd(z: Complexd) Complexd {
    return acoshValue(f64, z);
}

/// `atanh(z)`, with a real part in `(-1, +1)`.
pub fn atanhf(z: Complexf) Complexf {
    return atanhValue(f32, z);
}

/// `atanh(z)`, with a real part in `(-1, +1)`.
pub fn atanhd(z: Complexd) Complexd {
    return atanhValue(f64, z);
}
