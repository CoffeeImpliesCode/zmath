// C23 logp1 spellings that directly alias the established log1p kernels.

const log1p_impl = @import("scalar_log1p.zig");
const log1pd_impl = @import("scalar_log1pd.zig");

/// The C23 `logp1` spelling of `log1ph`, which this forwards to unchanged.
///
/// The domain, boundary, signed-zero, infinity and NaN behaviour is the contract
/// of `log1ph`, not a second one.
pub fn logp1h(value: f16) f16 {
    return log1p_impl.log1ph(value);
}

/// The C23 `logp1` spelling of `log1pf`, which this forwards to unchanged.
///
/// The domain, boundary, signed-zero, infinity and NaN behaviour is the contract
/// of `log1pf`, not a second one.
pub fn logp1f(value: f32) f32 {
    return log1p_impl.log1pf(value);
}

/// The C23 `logp1` spelling of `log1pd`, which this forwards to unchanged.
///
/// The domain, boundary, signed-zero, infinity and NaN behaviour is the contract
/// of `log1pd`, not a second one.
pub fn logp1d(value: f64) f64 {
    return log1pd_impl.log1pd(value);
}
