//! Shared result types and the `fromfp` rounding-direction enum.
//!
//! Result structs use automatic Zig layout, not a C ABI layout.
//! Do not reinterpret them as C structs. Copying a result retains no resource.
//! No status field requires clearing, and no result needs a cleanup call.

/// Rounding directions used by the C23 `fromfp` and `fromfpx` families.
///
/// These values select the conversion direction. They do not inspect or change
/// the floating-point environment, and a call that takes one neither reads nor
/// sets the process rounding mode.
pub const FromFpRound = enum {
    nearest_even,
    nearest_from_zero,
    toward_zero,
    toward_positive,
    toward_negative,
};

/// A rounded integer value and its validity status.
pub const IntegerResult = struct {
    /// The rounded integer, which is meaningful only when `ok` is true. It is
    /// 0 on every path that sets `ok` to false, where it is a sentinel rather
    /// than a converted value.
    value: i64,
    /// False when the input was a NaN, either infinity, or a finite magnitude
    /// outside the convertible window, so that no conversion traps and no
    /// caller has to reason about a trap.
    ok: bool,
};

/// An atan2 value and whether its result is resolved.
///
/// The binary16 and binary32 kernels report `resolved` true for every input.
/// Only the binary64 kernel can report false.
pub fn Atan2Result(comptime T: type) type {
    return struct {
        /// The angle, in radians in `[-pi, pi]` for `atan2`, or in half
        /// revolutions in `[-1, 1]` for `atan2pi`. When `resolved` is false,
        /// this is the positive canonical quiet NaN and carries no angle.
        value: T,
        /// Whether the kernel settled the result. False means the value is a NaN
        /// sentinel and the call must be retried or the case handled.
        resolved: bool,
    };
}

/// A power result and whether its value is valid.
pub fn PowerResult(comptime T: type) type {
    return struct {
        /// The computed power, compound power or root. It is meaningful only
        /// when `ok` is true. A domain case inside the accepted range returns
        /// `ok` true with the positive canonical quiet NaN here, which is a real
        /// result of the call rather than a failure sentinel.
        value: T,
        /// False only when the integer degree was outside the range `[-64, 64]`,
        /// the range the kernel implements. Then `value` is the positive
        /// canonical quiet NaN and must not be read as a power.
        ok: bool,
    };
}

/// A significand and binary exponent returned by `frexp`.
pub fn FrexpResult(comptime T: type) type {
    return struct {
        /// `value / 2^exponent`, in `[0.5, 1)` in magnitude and carrying the sign
        /// of the input. A zero, either sign, and an infinity come back unchanged
        /// here, and a NaN comes back quieted with its sign and payload kept.
        fraction: T,
        /// The power of two, chosen so `value == fraction * 2^exponent` for
        /// every finite nonzero input. It is 0 for a zero, an infinity and a
        /// NaN, where that identity does not apply.
        exponent: i32,
    };
}

/// Fractional and integral parts returned by `modf`.
pub fn ModfResult(comptime T: type) type {
    return struct {
        /// The fractional part, with the input's sign. It is the whole input
        /// for magnitudes below one, and signed zero for integral inputs or
        /// infinities. A NaN input leaves both fields holding its quieted NaN.
        fractional: T,
        /// The integral part, truncated toward zero, with the input's sign.
        /// It is signed zero for magnitudes below one and the input infinity
        /// for an infinite input.
        integral: T,
    };
}

/// A remainder and the low quotient bits returned by `remquo`.
pub fn RemquoResult(comptime T: type) type {
    return struct {
        /// The IEEE remainder, the same value `remainder` returns, whose
        /// magnitude is at most half of the divisor's. It is the positive
        /// canonical quiet NaN when the operation is invalid.
        value: T,
        /// The low three bits of the integer quotient, signed, so it lies in
        /// `[-7, 7]`. It is 0 on every invalid or NaN path.
        quotient: i32,
    };
}

/// The rounded logarithm of |gamma(x)| and the sign of gamma(x).
pub fn LgammaResult(comptime T: type) type {
    return struct {
        /// The rounded value of `log(abs(Gamma(x)))`. It is `+infinity` at a
        /// finite negative-integer pole and at either zero, and a NaN input
        /// comes back with its quiet bit set and its sign and payload kept.
        value: T,
        /// The sign of `Gamma(x)`, which is `-1` or `1` and never `0`. The `i2`
        /// field type is wider than the values reported: `0` and `-2` are not
        /// produced by any kernel. It is `-1` at a finite negative-integer pole
        /// and at `-0.0`, and `1` at `+0.0`, at either infinity and for a NaN.
        sign: i2,
    };
}

/// A payload operation result and its success status.
pub fn PayloadResult(comptime T: type) type {
    return struct {
        /// The positive quiet NaN or signaling NaN constructed by the selected
        /// setter. When `ok` is false, it is `+0`, a sentinel rather than a NaN
        /// result with a zero payload.
        value: T,
        /// Whether the setter accepted the numeric payload. Quiet setters
        /// accept zero, including `-0`. Signaling setters require nonzero.
        /// Both reject negative nonzero, fractional, non-finite, or out-of-range
        /// inputs. The payload widths are 9, 22, and 51 bits for f16/f32/f64.
        ok: bool,
    };
}

/// A canonicalized value and its success status.
pub fn CanonicalizeResult(comptime T: type) type {
    return struct {
        /// The input with any signaling NaN quieted: the sign and the diagnostic
        /// payload are kept and a quiet NaN is returned bit-for-bit unchanged.
        value: T,
        /// Whether the canonicalization completed. It is true on every path for
        /// all three supported formats, so no input on this surface reports
        /// failure; the field exists to keep the C23 result shape.
        ok: bool,
    };
}
