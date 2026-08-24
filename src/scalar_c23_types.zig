/// Rounding directions used by the C23 `fromfp` and `fromfpx` families.
///
/// These values select the conversion direction. They do not inspect or change
/// the floating-point environment.
pub const FromFpRound = enum {
    nearest_even,
    nearest_from_zero,
    toward_zero,
    toward_positive,
    toward_negative,
};

/// A rounded integer value and its validity status.
pub const IntegerResult = struct {
    value: i64,
    ok: bool,
};

/// An atan2 value and whether its result is resolved.
pub fn Atan2Result(comptime T: type) type {
    return struct {
        value: T,
        resolved: bool,
    };
}

/// A power result and whether its value is valid.
pub fn PowerResult(comptime T: type) type {
    return struct {
        value: T,
        ok: bool,
    };
}

/// A significand and binary exponent returned by `frexp`.
pub fn FrexpResult(comptime T: type) type {
    return struct {
        fraction: T,
        exponent: i32,
    };
}

/// Fractional and integral parts returned by `modf`.
pub fn ModfResult(comptime T: type) type {
    return struct {
        fractional: T,
        integral: T,
    };
}

/// A remainder and the low quotient bits returned by `remquo`.
pub fn RemquoResult(comptime T: type) type {
    return struct {
        value: T,
        quotient: i32,
    };
}

/// The rounded logarithm of |gamma(x)| and the sign of gamma(x).
pub fn LgammaResult(comptime T: type) type {
    return struct {
        value: T,
        sign: i2,
    };
}

/// A payload operation result and its success status.
pub fn PayloadResult(comptime T: type) type {
    return struct {
        value: T,
        ok: bool,
    };
}

/// A canonicalized value and its success status.
pub fn CanonicalizeResult(comptime T: type) type {
    return struct {
        value: T,
        ok: bool,
    };
}
