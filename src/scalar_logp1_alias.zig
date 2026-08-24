// C23 logp1 spellings that directly alias the established log1p kernels.

const log1p_impl = @import("scalar_log1p.zig");
const log1pd_impl = @import("scalar_log1pd.zig");

pub fn logp1h(value: f16) f16 {
    return log1p_impl.log1ph(value);
}

pub fn logp1f(value: f32) f32 {
    return log1p_impl.log1pf(value);
}

pub fn logp1d(value: f64) f64 {
    return log1pd_impl.log1pd(value);
}
