const acosh_impl = @import("scalar_acosh.zig");
pub const acoshh = acosh_impl.acoshh;
pub const acoshf = acosh_impl.acoshf;
pub const acoshd = acosh_impl.acoshd;

const asin_impl = @import("scalar_asin.zig");
pub const asinh = asin_impl.asinh;
pub const asinf = asin_impl.asinf;
pub const asind = asin_impl.asind;
const asinh_impl = @import("scalar_asinhh.zig");
pub const asinhh = asinh_impl.asinhh;
const asinhf_impl = @import("scalar_asinhf.zig");
pub const asinhf = asinhf_impl.asinhf;
const asinhd_impl = @import("scalar_asinhd.zig");
pub const asinhd = asinhd_impl.asinhd;

const sinhh_impl = @import("scalar_sinhh.zig");
pub const sinhh = sinhh_impl.sinhh;
const sinhf_impl = @import("scalar_sinhf.zig");
pub const sinhf = sinhf_impl.sinhf;
const sinhd_impl = @import("scalar_sinhd.zig");
pub const sinhd = sinhd_impl.sinhd;

const cbrt_impl = @import("scalar_cbrt.zig");
pub const cbrth = cbrt_impl.cbrth;
pub const cbrtf = cbrt_impl.cbrtf;
pub const cbrtd = cbrt_impl.cbrtd;

const hypot_impl = @import("scalar_hypot.zig");
pub const hypoth = hypot_impl.hypoth;
pub const hypotf = hypot_impl.hypotf;
pub const hypotd = hypot_impl.hypotd;
const atan_impl = @import("scalar_atan.zig");
pub const atanh = atan_impl.atanh;
pub const atanf = atan_impl.atanf;
pub const atand = atan_impl.atand;
const atanh_impl = @import("scalar_atanh.zig");
pub const atanhh = atanh_impl.atanhh;
pub const atanhf = atanh_impl.atanhf;
pub const atanhd = atanh_impl.atanhd;
pub const expm1h = @import("scalar_expm1h.zig").expm1h;
pub const expm1f = @import("scalar_expm1f.zig").expm1f;
pub const expm1d = @import("scalar_expm1d.zig").expm1d;
const logd_impl = @import("logd.zig");
const log1p_impl = @import("scalar_log1p.zig");
const log1pd_impl = @import("scalar_log1pd.zig");

pub const logd = logd_impl.logd;
pub const log1ph = log1p_impl.log1ph;
pub const log1pf = log1p_impl.log1pf;
pub const log1pd = log1pd_impl.log1pd;
pub const sind = @import("scalar_sind.zig").sind;
pub const cosd = @import("scalar_sind.zig").cosd;
const exp2_impl = @import("scalar_exp2.zig");
pub const exp2h = exp2_impl.exp2h;
pub const exp2f = exp2_impl.exp2f;
pub const exp2d = exp2_impl.exp2d;
const rsqrt_impl = @import("scalar_rsqrt.zig");
pub const rsqrth = rsqrt_impl.rsqrth;
pub const rsqrtf = rsqrt_impl.rsqrtf;
pub const rsqrtd = rsqrt_impl.rsqrtd;
pub const tanh = @import("scalar_tanh.zig").tanh;
pub const tanf = @import("scalar_tanf.zig").tanf;
pub const tand = @import("scalar_tand.zig").tand;

const cosh_impl = @import("scalar_cosh.zig");
pub const coshh = cosh_impl.coshh;
pub const coshf = cosh_impl.coshf;
pub const coshd = cosh_impl.coshd;

const tanh_impl = @import("scalar_tanh.zig");
pub const tanhh = tanh_impl.tanhh;
pub const tanhf = tanh_impl.tanhf;
pub const tanhd = tanh_impl.tanhd;

const erfh_impl = @import("scalar_erfh.zig");
pub const erfh = erfh_impl.erfh;
const erff_impl = @import("scalar_erff.zig");
pub const erff = erff_impl.erff;
const erfd_impl = @import("scalar_erfd.zig");
pub const erfd = erfd_impl.erfd;

const sign_mask_h: u16 = 0x8000;
const sign_mask_f: u32 = 0x8000_0000;
const sign_mask_d: u64 = 0x8000_0000_0000_0000;

const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;

const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;

// The binary16 exponential kernel and tables are ported from CORE-MATH
// src/binary16/exp/expf16.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa (MIT license).
const exph_data = @import("scalar_exph_data.zig");
const logf_data = @import("scalar_logf_data.zig");

const logf_quiet_nan_bit: u32 = 0x0040_0000;
const logf_canonical_nan_bits: u32 = 0x7fc0_0000;
const logf_minus_infinity_bits: u32 = 0xff80_0000;
const logf_fast_coefficients = [_]f64{
    0x1.00000006342eap+0,
    -0x1.0001f7fdc3977p-1,
    0x1.554a4e5cae9cfp-2,
};
const logf_accurate_coefficients = [_]f64{
    -0x1p-1,
    0x1.55555555571cap-2,
    -0x1.0000000002d85p-2,
    0x1.9999987d0c963p-3,
    -0x1.555554059a8bbp-3,
    0x1.24aebcf71a38fp-3,
    -0x1.001c73915d758p-3,
};

pub fn exph(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude > infinity_bits_h) {
        const quiet_nan_bit_h: u16 = 0x0200;
        return @bitCast(bits | quiet_nan_bit_h);
    }
    if (magnitude == infinity_bits_h) {
        if ((bits & sign_mask_h) == 0) return value;
        return @bitCast(@as(u16, 0));
    }
    if (magnitude == 0) return @bitCast(@as(u16, 0x3c00));

    // Decompose x into x1 + x2, then use exp(x) = exp(x1) * exp(x2).
    const table1_index = bits >> 5;
    const table2_index = ((bits >> 10) << 5) | (bits & 0x001f);
    const exp_x1: f32 = @bitCast(exph_data.t1[table1_index]);
    const exp_x2: f32 = @bitCast(exph_data.t2[table2_index]);
    return @floatCast(exp_x1 * exp_x2);
}

// Correctly rounded natural logarithm for binary16 values.
//
// Copyright (c) 2025 Maxence Ponsardin and Paul Zimmermann.
//
// This code is derived from CORE-MATH src/binary16/log/logf16.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
const logh_quiet_nan_bit: u16 = 0x0200;
const logh_canonical_nan_bits: u16 = 0x7e00;
const logh_minus_infinity_bits: u16 = 0xfc00;

// For 0 <= i1 < 2^5, T1[i1] approximates log(x1), where the binary16
// encoding of x1 is i1 * 2^10.
const logh_t1 = [32]f32{
    0.0,            -0x1.3687aap+3, -0x1.205966p+3, -0x1.0a2b24p+3,
    -0x1.e7f9c2p+2, -0x1.bb9d3cp+2, -0x1.8f40b6p+2, -0x1.62e43p+2,
    -0x1.3687aap+2, -0x1.0a2b24p+2, -0x1.bb9d3cp+1, -0x1.62e43p+1,
    -0x1.0a2b24p+1, -0x1.62e43p+0,  -0x1.62e43p-1,  0x0p+0,
    0x1.62e43p-1,   0x1.62e43p+0,   0x1.0a2b24p+1,  0x1.62e43p+1,
    0x1.bb9d3cp+1,  0x1.0a2b24p+2,  0x1.3687aap+2,  0x1.62e43p+2,
    0x1.8f40b6p+2,  0x1.bb9d3cp+2,  0x1.e7f9c2p+2,  0x1.0a2b24p+3,
    0x1.205966p+3,  0x1.3687aap+3,  0x1.4cb5ecp+3,  0.0,
};

// For 0 <= i2 < 2^10, T2[i2] approximates log(x2), where 1 <= x2 < 2
// and i2 matches the ten lower bits of the binary16 encoding of x2.
const logh_t2 = [1024]f32{
    0x0p+0,        0x1.ffc00ap-11, 0x1.ff802ap-10, 0x1.7f7048p-9,
    0x1.ff00aap-9, 0x1.3f38a6p-8,  0x1.7ee11ep-8,  0x1.be79c8p-8,
    0x1.fe02a6p-8, 0x1.1ebde2p-7,  0x1.3e7296p-7,  0x1.5e1f7p-7,
    0x1.7dc476p-7, 0x1.9d61aap-7,  0x1.bcf712p-7,  0x1.dc84b2p-7,
    0x1.fc0a8cp-7, 0x1.0dc452p-6,  0x1.1d7f7ep-6,  0x1.2d36cep-6,
    0x1.3cea44p-6, 0x1.4c99ep-6,   0x1.5c45a6p-6,  0x1.6bed94p-6,
    0x1.7b91bp-6,  0x1.8b31fap-6,  0x1.9ace76p-6,  0x1.aa6722p-6,
    0x1.b9fc02p-6, 0x1.c98d18p-6,  0x1.d91a66p-6,  0x1.e8a3eep-6,
    0x1.f829bp-6,  0x1.03d5d8p-5,  0x1.0b94f8p-5,  0x1.135238p-5,
    0x1.1b0d98p-5, 0x1.22c71cp-5,  0x1.2a7ec2p-5,  0x1.32348cp-5,
    0x1.39e87cp-5, 0x1.419a9p-5,   0x1.494accp-5,  0x1.50f93p-5,
    0x1.58a5bap-5, 0x1.60507p-5,   0x1.67f95p-5,   0x1.6fa05ap-5,
    0x1.77459p-5,  0x1.7ee8f2p-5,  0x1.868a84p-5,  0x1.8e2a42p-5,
    0x1.95c83p-5,  0x1.9d645p-5,   0x1.a4feap-5,   0x1.ac9722p-5,
    0x1.b42dd8p-5, 0x1.bbc2cp-5,   0x1.c355dep-5,  0x1.cae73p-5,
    0x1.d276b8p-5, 0x1.da0478p-5,  0x1.e1907p-5,   0x1.e91aa2p-5,
    0x1.f0a30cp-5, 0x1.f829bp-5,   0x1.ffae92p-5,  0x1.0398d6p-4,
    0x1.075984p-4, 0x1.0b194ep-4,  0x1.0ed83ap-4,  0x1.129644p-4,
    0x1.16536ep-4, 0x1.1a0fbap-4,  0x1.1dcb26p-4,  0x1.2185b4p-4,
    0x1.253f62p-4, 0x1.28f834p-4,  0x1.2cb028p-4,  0x1.30674p-4,
    0x1.341d7ap-4, 0x1.37d2d8p-4,  0x1.3b875ap-4,  0x1.3f3bp-4,
    0x1.42edccp-4, 0x1.469fbcp-4,  0x1.4a50d4p-4,  0x1.4e011p-4,
    0x1.51b074p-4, 0x1.555efep-4,  0x1.590cbp-4,   0x1.5cb98ap-4,
    0x1.60658ap-4, 0x1.6410b4p-4,  0x1.67bb08p-4,  0x1.6b6484p-4,
    0x1.6f0d28p-4, 0x1.72b4f8p-4,  0x1.765bf2p-4,  0x1.7a0216p-4,
    0x1.7da766p-4, 0x1.814be2p-4,  0x1.84ef8ap-4,  0x1.88925ep-4,
    0x1.8c345ep-4, 0x1.8fd58ap-4,  0x1.9375e6p-4,  0x1.97156ep-4,
    0x1.9ab424p-4, 0x1.9e520ap-4,  0x1.a1ef1ep-4,  0x1.a58b6p-4,
    0x1.a926d4p-4, 0x1.acc176p-4,  0x1.b05b4ap-4,  0x1.b3f44ep-4,
    0x1.b78c82p-4, 0x1.bb23eap-4,  0x1.beba82p-4,  0x1.c2504cp-4,
    0x1.c5e548p-4, 0x1.c97978p-4,  0x1.cd0cdcp-4,  0x1.d09f72p-4,
    0x1.d4313ep-4, 0x1.d7c23cp-4,  0x1.db527p-4,   0x1.dee1d8p-4,
    0x1.e27076p-4, 0x1.e5fe4ap-4,  0x1.e98b54p-4,  0x1.ed1794p-4,
    0x1.f0a30cp-4, 0x1.f42dbap-4,  0x1.f7b7ap-4,   0x1.fb40bep-4,
    0x1.fec914p-4, 0x1.01285p-3,   0x1.02ebb4p-3,  0x1.04aeb4p-3,
    0x1.067152p-3, 0x1.08338ap-3,  0x1.09f562p-3,  0x1.0bb6d6p-3,
    0x1.0d77e8p-3, 0x1.0f3898p-3,  0x1.10f8e4p-3,  0x1.12b8dp-3,
    0x1.147858p-3, 0x1.16378p-3,   0x1.17f646p-3,  0x1.19b4aap-3,
    0x1.1b72aep-3, 0x1.1d305p-3,   0x1.1eed9p-3,   0x1.20aa72p-3,
    0x1.2266f2p-3, 0x1.242312p-3,  0x1.25dedp-3,   0x1.279a3p-3,
    0x1.29553p-3,  0x1.2b0fdp-3,   0x1.2cca1p-3,   0x1.2e83fp-3,
    0x1.303d72p-3, 0x1.31f694p-3,  0x1.33af58p-3,  0x1.3567bcp-3,
    0x1.371fc2p-3, 0x1.38d76ap-3,  0x1.3a8eb2p-3,  0x1.3c459ep-3,
    0x1.3dfc2cp-3, 0x1.3fb25ap-3,  0x1.41682cp-3,  0x1.431dap-3,
    0x1.44d2b6p-3, 0x1.46877p-3,   0x1.483bccp-3,  0x1.49efccp-3,
    0x1.4ba37p-3,  0x1.4d56b6p-3,  0x1.4f09ap-3,   0x1.50bc2cp-3,
    0x1.526e5ep-3, 0x1.542034p-3,  0x1.55d1aep-3,  0x1.5782ccp-3,
    0x1.59338ep-3, 0x1.5ae3f4p-3,  0x1.5c94p-3,    0x1.5e43b2p-3,
    0x1.5ff308p-3, 0x1.61a202p-3,  0x1.6350a2p-3,  0x1.64fee8p-3,
    0x1.66acd4p-3, 0x1.685a66p-3,  0x1.6a079ep-3,  0x1.6bb47ap-3,
    0x1.6d60fep-3, 0x1.6f0d28p-3,  0x1.70b8fap-3,  0x1.72647p-3,
    0x1.740f9p-3,  0x1.75ba54p-3,  0x1.7764c2p-3,  0x1.790ed4p-3,
    0x1.7ab89p-3,  0x1.7c61f2p-3,  0x1.7e0afep-3,  0x1.7fb3bp-3,
    0x1.815c0ap-3, 0x1.83040cp-3,  0x1.84abb8p-3,  0x1.86530ap-3,
    0x1.87fa06p-3, 0x1.89a0aap-3,  0x1.8b46f8p-3,  0x1.8ceceep-3,
    0x1.8e928ep-3, 0x1.9037d6p-3,  0x1.91dcc8p-3,  0x1.938164p-3,
    0x1.9525aap-3, 0x1.96c99ap-3,  0x1.986d32p-3,  0x1.9a1076p-3,
    0x1.9bb362p-3, 0x1.9d55fap-3,  0x1.9ef83ep-3,  0x1.a09a2ap-3,
    0x1.a23bc2p-3, 0x1.a3dd04p-3,  0x1.a57df2p-3,  0x1.a71e8cp-3,
    0x1.a8bedp-3,  0x1.aa5ecp-3,   0x1.abfe5ap-3,  0x1.ad9da2p-3,
    0x1.af3c94p-3, 0x1.b0db34p-3,  0x1.b2797ep-3,  0x1.b41776p-3,
    0x1.b5b51ap-3, 0x1.b7526ap-3,  0x1.b8ef68p-3,  0x1.ba8c1p-3,
    0x1.bc2868p-3, 0x1.bdc46ap-3,  0x1.bf601cp-3,  0x1.c0fb7ap-3,
    0x1.c29686p-3, 0x1.c4313ep-3,  0x1.c5cba6p-3,  0x1.c765bap-3,
    0x1.c8ff7cp-3, 0x1.ca98eep-3,  0x1.cc320cp-3,  0x1.cdcadap-3,
    0x1.cf6354p-3, 0x1.d0fb8p-3,   0x1.d29358p-3,  0x1.d42aep-3,
    0x1.d5c216p-3, 0x1.d758fcp-3,  0x1.d8ef92p-3,  0x1.da85d6p-3,
    0x1.dc1bcap-3, 0x1.ddb16ep-3,  0x1.df46cp-3,   0x1.e0dbc4p-3,
    0x1.e27076p-3, 0x1.e404dap-3,  0x1.e598eep-3,  0x1.e72cb2p-3,
    0x1.e8c026p-3, 0x1.ea534ap-3,  0x1.ebe62p-3,   0x1.ed78a6p-3,
    0x1.ef0adcp-3, 0x1.f09cc6p-3,  0x1.f22e5ep-3,  0x1.f3bfaap-3,
    0x1.f550a6p-3, 0x1.f6e154p-3,  0x1.f871b2p-3,  0x1.fa01c4p-3,
    0x1.fb9186p-3, 0x1.fd20fcp-3,  0x1.feb024p-3,  0x1.001f7ep-2,
    0x1.00e6c4p-2, 0x1.01ade4p-2,  0x1.0274dcp-2,  0x1.033baep-2,
    0x1.04025ap-2, 0x1.04c8dep-2,  0x1.058f3cp-2,  0x1.065574p-2,
    0x1.071b86p-2, 0x1.07e172p-2,  0x1.08a736p-2,  0x1.096cd6p-2,
    0x1.0a324ep-2, 0x1.0af7ap-2,   0x1.0bbccep-2,  0x1.0c81d4p-2,
    0x1.0d46b6p-2, 0x1.0e0b72p-2,  0x1.0ed006p-2,  0x1.0f9476p-2,
    0x1.1058cp-2,  0x1.111ce4p-2,  0x1.11e0e2p-2,  0x1.12a4bcp-2,
    0x1.13687p-2,  0x1.142bfep-2,  0x1.14ef68p-2,  0x1.15b2acp-2,
    0x1.1675cap-2, 0x1.1738c4p-2,  0x1.17fb98p-2,  0x1.18be48p-2,
    0x1.1980d2p-2, 0x1.1a4338p-2,  0x1.1b057ap-2,  0x1.1bc794p-2,
    0x1.1c898cp-2, 0x1.1d4b5ep-2,  0x1.1e0d0cp-2,  0x1.1ece96p-2,
    0x1.1f8ffap-2, 0x1.20513ap-2,  0x1.211256p-2,  0x1.21d34cp-2,
    0x1.22942p-2,  0x1.2354cep-2,  0x1.241558p-2,  0x1.24d5bep-2,
    0x1.259602p-2, 0x1.26562p-2,   0x1.27161ap-2,  0x1.27d5fp-2,
    0x1.2895a2p-2, 0x1.29553p-2,   0x1.2a149ap-2,  0x1.2ad3ep-2,
    0x1.2b9304p-2, 0x1.2c5204p-2,  0x1.2d10dep-2,  0x1.2dcf96p-2,
    0x1.2e8e2cp-2, 0x1.2f4c9cp-2,  0x1.300aeap-2,  0x1.30c916p-2,
    0x1.31871cp-2, 0x1.3245p-2,    0x1.3302c2p-2,  0x1.33c06p-2,
    0x1.347ddap-2, 0x1.353b32p-2,  0x1.35f866p-2,  0x1.36b578p-2,
    0x1.377266p-2, 0x1.382f32p-2,  0x1.38ebdcp-2,  0x1.39a862p-2,
    0x1.3a64c6p-2, 0x1.3b2106p-2,  0x1.3bdd24p-2,  0x1.3c992p-2,
    0x1.3d54fap-2, 0x1.3e10b2p-2,  0x1.3ecc46p-2,  0x1.3f87b8p-2,
    0x1.404312p-2, 0x1.40fe36p-2,  0x1.41b942p-2,  0x1.42742cp-2,
    0x1.432ef2p-2, 0x1.43e998p-2,  0x1.44a41cp-2,  0x1.455e7cp-2,
    0x1.4618bcp-2, 0x1.46d2dap-2,  0x1.478cd6p-2,  0x1.4846bp-2,
    0x1.490068p-2, 0x1.49b9fep-2,  0x1.4a7374p-2,  0x1.4b2cbep-2,
    0x1.4be5fap-2, 0x1.4c9f0ap-2,  0x1.4d57f8p-2,  0x1.4e10c6p-2,
    0x1.4ec974p-2, 0x1.4f81fep-2,  0x1.503a68p-2,  0x1.50f2bp-2,
    0x1.51aad8p-2, 0x1.5262dep-2,  0x1.531ac4p-2,  0x1.53d288p-2,
    0x1.548a2cp-2, 0x1.5541aep-2,  0x1.55f91p-2,   0x1.56b052p-2,
    0x1.576772p-2, 0x1.581e7p-2,   0x1.58d55p-2,   0x1.598c0ep-2,
    0x1.5a42acp-2, 0x1.5af928p-2,  0x1.5baf84p-2,  0x1.5c65cp-2,
    0x1.5d1bdcp-2, 0x1.5dd1d8p-2,  0x1.5e87b2p-2,  0x1.5f3d6cp-2,
    0x1.5ff308p-2, 0x1.60a882p-2,  0x1.615ddcp-2,  0x1.621316p-2,
    0x1.62c83p-2,  0x1.637d2ap-2,  0x1.643204p-2,  0x1.64e6bep-2,
    0x1.659b58p-2, 0x1.664fd2p-2,  0x1.67042cp-2,  0x1.67b866p-2,
    0x1.686c82p-2, 0x1.69207ep-2,  0x1.69d45ap-2,  0x1.6a8816p-2,
    0x1.6b3bb2p-2, 0x1.6bef3p-2,   0x1.6ca28ep-2,  0x1.6d55ccp-2,
    0x1.6e08eap-2, 0x1.6ebbeap-2,  0x1.6f6ecap-2,  0x1.70218cp-2,
    0x1.70d42ep-2, 0x1.7186b2p-2,  0x1.723914p-2,  0x1.72eb5ap-2,
    0x1.739d8p-2,  0x1.744f86p-2,  0x1.75016ep-2,  0x1.75b338p-2,
    0x1.7664e2p-2, 0x1.77166cp-2,  0x1.77c7dap-2,  0x1.787926p-2,
    0x1.792a56p-2, 0x1.79db66p-2,  0x1.7a8c58p-2,  0x1.7b3d2cp-2,
    0x1.7bedep-2,  0x1.7c9e78p-2,  0x1.7d4eeep-2,  0x1.7dff48p-2,
    0x1.7eaf84p-2, 0x1.7f5fap-2,   0x1.800fap-2,   0x1.80bf8p-2,
    0x1.816f42p-2, 0x1.821ee6p-2,  0x1.82ce6cp-2,  0x1.837dd4p-2,
    0x1.842d1ep-2, 0x1.84dc4ap-2,  0x1.858b58p-2,  0x1.863a48p-2,
    0x1.86e91ap-2, 0x1.8797cep-2,  0x1.884664p-2,  0x1.88f4dep-2,
    0x1.89a338p-2, 0x1.8a5176p-2,  0x1.8aff96p-2,  0x1.8bad98p-2,
    0x1.8c5b7cp-2, 0x1.8d0944p-2,  0x1.8db6eep-2,  0x1.8e647ap-2,
    0x1.8f11e8p-2, 0x1.8fbf3ap-2,  0x1.906c6ep-2,  0x1.911986p-2,
    0x1.91c67ep-2, 0x1.92735cp-2,  0x1.93201ap-2,  0x1.93ccbcp-2,
    0x1.947942p-2, 0x1.9525aap-2,  0x1.95d1f4p-2,  0x1.967e22p-2,
    0x1.972a34p-2, 0x1.97d628p-2,  0x1.9882p-2,    0x1.992dbap-2,
    0x1.99d958p-2, 0x1.9a84dap-2,  0x1.9b303ep-2,  0x1.9bdb86p-2,
    0x1.9c86bp-2,  0x1.9d31bep-2,  0x1.9ddcbp-2,   0x1.9e8786p-2,
    0x1.9f323ep-2, 0x1.9fdcdcp-2,  0x1.a0875cp-2,  0x1.a131bep-2,
    0x1.a1dc06p-2, 0x1.a28632p-2,  0x1.a3304p-2,   0x1.a3da32p-2,
    0x1.a4840ap-2, 0x1.a52dc4p-2,  0x1.a5d762p-2,  0x1.a680e4p-2,
    0x1.a72a4ap-2, 0x1.a7d394p-2,  0x1.a87cc2p-2,  0x1.a925d4p-2,
    0x1.a9cecap-2, 0x1.aa77a4p-2,  0x1.ab2062p-2,  0x1.abc906p-2,
    0x1.ac718cp-2, 0x1.ad19f8p-2,  0x1.adc246p-2,  0x1.ae6a7ap-2,
    0x1.af1294p-2, 0x1.afba9p-2,   0x1.b06272p-2,  0x1.b10a36p-2,
    0x1.b1b1ep-2,  0x1.b2597p-2,   0x1.b300e4p-2,  0x1.b3a83cp-2,
    0x1.b44f78p-2, 0x1.b4f69ap-2,  0x1.b59dap-2,   0x1.b6448ap-2,
    0x1.b6eb5ap-2, 0x1.b7920ep-2,  0x1.b838a8p-2,  0x1.b8df26p-2,
    0x1.b9858ap-2, 0x1.ba2bd2p-2,  0x1.bad1fep-2,  0x1.bb7812p-2,
    0x1.bc1e08p-2, 0x1.bcc3e6p-2,  0x1.bd69a6p-2,  0x1.be0f4ep-2,
    0x1.beb4dap-2, 0x1.bf5a4cp-2,  0x1.bfffa2p-2,  0x1.c0a4dep-2,
    0x1.c14ap-2,   0x1.c1ef06p-2,  0x1.c293f2p-2,  0x1.c338c4p-2,
    0x1.c3dd7ap-2, 0x1.c48218p-2,  0x1.c5269ap-2,  0x1.c5cbp-2,
    0x1.c66f4ep-2, 0x1.c71382p-2,  0x1.c7b79ap-2,  0x1.c85b98p-2,
    0x1.c8ff7cp-2, 0x1.c9a346p-2,  0x1.ca46f6p-2,  0x1.caea8cp-2,
    0x1.cb8e08p-2, 0x1.cc3168p-2,  0x1.ccd4bp-2,   0x1.cd77dep-2,
    0x1.ce1afp-2,  0x1.cebdeap-2,  0x1.cf60cap-2,  0x1.d0039p-2,
    0x1.d0a63ap-2, 0x1.d148ccp-2,  0x1.d1eb44p-2,  0x1.d28da4p-2,
    0x1.d32fe8p-2, 0x1.d3d212p-2,  0x1.d47424p-2,  0x1.d5161cp-2,
    0x1.d5b7fap-2, 0x1.d659bep-2,  0x1.d6fb6ap-2,  0x1.d79cfap-2,
    0x1.d83e72p-2, 0x1.d8dfdp-2,   0x1.d98116p-2,  0x1.da2242p-2,
    0x1.dac354p-2, 0x1.db644cp-2,  0x1.dc052cp-2,  0x1.dca5f4p-2,
    0x1.dd46ap-2,  0x1.dde734p-2,  0x1.de87bp-2,   0x1.df281p-2,
    0x1.dfc85ap-2, 0x1.e0688ap-2,  0x1.e108ap-2,   0x1.e1a89ep-2,
    0x1.e24882p-2, 0x1.e2e84ep-2,  0x1.e388p-2,    0x1.e4279ap-2,
    0x1.e4c71ap-2, 0x1.e56682p-2,  0x1.e605d2p-2,  0x1.e6a508p-2,
    0x1.e74426p-2, 0x1.e7e32cp-2,  0x1.e88218p-2,  0x1.e920ecp-2,
    0x1.e9bfa6p-2, 0x1.ea5e48p-2,  0x1.eafcd2p-2,  0x1.eb9b44p-2,
    0x1.ec399ep-2, 0x1.ecd7dep-2,  0x1.ed7606p-2,  0x1.ee1416p-2,
    0x1.eeb20cp-2, 0x1.ef4fecp-2,  0x1.efedb2p-2,  0x1.f08b6p-2,
    0x1.f128f6p-2, 0x1.f1c674p-2,  0x1.f263dap-2,  0x1.f30126p-2,
    0x1.f39e5cp-2, 0x1.f43b78p-2,  0x1.f4d87ep-2,  0x1.f5756ap-2,
    0x1.f6124p-2,  0x1.f6aefcp-2,  0x1.f74ba2p-2,  0x1.f7e82ep-2,
    0x1.f884a4p-2, 0x1.f921p-2,    0x1.f9bd46p-2,  0x1.fa5974p-2,
    0x1.faf588p-2, 0x1.fb9186p-2,  0x1.fc2d6cp-2,  0x1.fcc93cp-2,
    0x1.fd64f2p-2, 0x1.fe0092p-2,  0x1.fe9c18p-2,  0x1.ff3788p-2,
    0x1.ffd2ep-2,  0x1.00371p-1,   0x1.0084a6p-1,  0x1.00d22ep-1,
    0x1.011facp-1, 0x1.016d1cp-1,  0x1.01ba82p-1,  0x1.0207dcp-1,
    0x1.02552ap-1, 0x1.02a26cp-1,  0x1.02efa4p-1,  0x1.033cdp-1,
    0x1.0389eep-1, 0x1.03d704p-1,  0x1.04240cp-1,  0x1.047108p-1,
    0x1.04bdfap-1, 0x1.050aep-1,   0x1.0557bap-1,  0x1.05a488p-1,
    0x1.05f14cp-1, 0x1.063e04p-1,  0x1.068abp-1,   0x1.06d75p-1,
    0x1.0723e6p-1, 0x1.07707p-1,   0x1.07bceep-1,  0x1.08096p-1,
    0x1.0855c8p-1, 0x1.08a224p-1,  0x1.08ee76p-1,  0x1.093abap-1,
    0x1.0986f4p-1, 0x1.09d324p-1,  0x1.0a1f48p-1,  0x1.0a6b6p-1,
    0x1.0ab76cp-1, 0x1.0b036ep-1,  0x1.0b4f64p-1,  0x1.0b9b4ep-1,
    0x1.0be72ep-1, 0x1.0c3302p-1,  0x1.0c7eccp-1,  0x1.0cca8ap-1,
    0x1.0d163cp-1, 0x1.0d61e4p-1,  0x1.0dad8p-1,   0x1.0df912p-1,
    0x1.0e4498p-1, 0x1.0e9014p-1,  0x1.0edb84p-1,  0x1.0f26e8p-1,
    0x1.0f7242p-1, 0x1.0fbd9p-1,   0x1.1008d4p-1,  0x1.10540cp-1,
    0x1.109f3ap-1, 0x1.10ea5cp-1,  0x1.113574p-1,  0x1.11808p-1,
    0x1.11cb82p-1, 0x1.121678p-1,  0x1.126164p-1,  0x1.12ac44p-1,
    0x1.12f71ap-1, 0x1.1341e4p-1,  0x1.138ca4p-1,  0x1.13d758p-1,
    0x1.142202p-1, 0x1.146ca2p-1,  0x1.14b736p-1,  0x1.1501bep-1,
    0x1.154c3ep-1, 0x1.1596bp-1,   0x1.15e11ap-1,  0x1.162b78p-1,
    0x1.1675cap-1, 0x1.16c014p-1,  0x1.170a5p-1,   0x1.175484p-1,
    0x1.179eacp-1, 0x1.17e8cap-1,  0x1.1832dcp-1,  0x1.187ce4p-1,
    0x1.18c6ep-1,  0x1.1910d4p-1,  0x1.195abcp-1,  0x1.19a498p-1,
    0x1.19ee6cp-1, 0x1.1a3834p-1,  0x1.1a81fp-1,   0x1.1acba4p-1,
    0x1.1b154cp-1, 0x1.1b5ee8p-1,  0x1.1ba87cp-1,  0x1.1bf204p-1,
    0x1.1c3b82p-1, 0x1.1c84f6p-1,  0x1.1cce5ep-1,  0x1.1d17bcp-1,
    0x1.1d611p-1,  0x1.1daa5ap-1,  0x1.1df398p-1,  0x1.1e3cccp-1,
    0x1.1e85f6p-1, 0x1.1ecf16p-1,  0x1.1f182ap-1,  0x1.1f6134p-1,
    0x1.1faa34p-1, 0x1.1ff32ap-1,  0x1.203c16p-1,  0x1.2084f6p-1,
    0x1.20cdcep-1, 0x1.21169ap-1,  0x1.215f5cp-1,  0x1.21a812p-1,
    0x1.21f0cp-1,  0x1.223962p-1,  0x1.2281fcp-1,  0x1.22ca8ap-1,
    0x1.23130ep-1, 0x1.235b88p-1,  0x1.23a3f6p-1,  0x1.23ec5cp-1,
    0x1.2434b6p-1, 0x1.247d08p-1,  0x1.24c54ep-1,  0x1.250d8ap-1,
    0x1.2555bcp-1, 0x1.259de4p-1,  0x1.25e602p-1,  0x1.262e16p-1,
    0x1.26762p-1,  0x1.26be2p-1,   0x1.270614p-1,  0x1.274ep-1,
    0x1.2795e2p-1, 0x1.27ddb8p-1,  0x1.282586p-1,  0x1.286d48p-1,
    0x1.28b5p-1,   0x1.28fcbp-1,   0x1.294454p-1,  0x1.298bfp-1,
    0x1.29d38p-1,  0x1.2a1b06p-1,  0x1.2a6284p-1,  0x1.2aa9f6p-1,
    0x1.2af16p-1,  0x1.2b38bep-1,  0x1.2b8012p-1,  0x1.2bc75ep-1,
    0x1.2c0e9ep-1, 0x1.2c55d6p-1,  0x1.2c9d04p-1,  0x1.2ce426p-1,
    0x1.2d2b4p-1,  0x1.2d725p-1,   0x1.2db956p-1,  0x1.2e0052p-1,
    0x1.2e4744p-1, 0x1.2e8e2cp-1,  0x1.2ed50ap-1,  0x1.2f1bdep-1,
    0x1.2f62aap-1, 0x1.2fa96ap-1,  0x1.2ff022p-1,  0x1.3036dp-1,
    0x1.307d74p-1, 0x1.30c40ep-1,  0x1.310a9ep-1,  0x1.315124p-1,
    0x1.3197ap-1,  0x1.31de14p-1,  0x1.32247ep-1,  0x1.326adep-1,
    0x1.32b134p-1, 0x1.32f78p-1,   0x1.333dc2p-1,  0x1.3383fcp-1,
    0x1.33ca2cp-1, 0x1.341052p-1,  0x1.34566ep-1,  0x1.349c8p-1,
    0x1.34e28ap-1, 0x1.35288ap-1,  0x1.356e8p-1,   0x1.35b46cp-1,
    0x1.35fa4ep-1, 0x1.364028p-1,  0x1.3685f8p-1,  0x1.36cbbep-1,
    0x1.37117cp-1, 0x1.37572ep-1,  0x1.379cd8p-1,  0x1.37e278p-1,
    0x1.38281p-1,  0x1.386d9ep-1,  0x1.38b322p-1,  0x1.38f89cp-1,
    0x1.393e0ep-1, 0x1.398374p-1,  0x1.39c8d4p-1,  0x1.3a0e28p-1,
    0x1.3a5374p-1, 0x1.3a98b6p-1,  0x1.3addeep-1,  0x1.3b231ep-1,
    0x1.3b6844p-1, 0x1.3bad62p-1,  0x1.3bf274p-1,  0x1.3c378p-1,
    0x1.3c7c8p-1,  0x1.3cc178p-1,  0x1.3d0666p-1,  0x1.3d4b4ap-1,
    0x1.3d9026p-1, 0x1.3dd4fap-1,  0x1.3e19c2p-1,  0x1.3e5e82p-1,
    0x1.3ea33ap-1, 0x1.3ee7e6p-1,  0x1.3f2c8cp-1,  0x1.3f7126p-1,
    0x1.3fb5b8p-1, 0x1.3ffa42p-1,  0x1.403ecp-1,   0x1.408338p-1,
    0x1.40c7a4p-1, 0x1.410c08p-1,  0x1.415064p-1,  0x1.4194b6p-1,
    0x1.41d8fep-1, 0x1.421d3ep-1,  0x1.426174p-1,  0x1.42a5a2p-1,
    0x1.42e9c6p-1, 0x1.432de2p-1,  0x1.4371f4p-1,  0x1.43b5fep-1,
    0x1.43f9fep-1, 0x1.443df6p-1,  0x1.4481e4p-1,  0x1.44c5c8p-1,
    0x1.4509a6p-1, 0x1.454d78p-1,  0x1.459142p-1,  0x1.45d504p-1,
    0x1.4618bcp-1, 0x1.465c6cp-1,  0x1.46a012p-1,  0x1.46e3bp-1,
    0x1.472744p-1, 0x1.476adp-1,   0x1.47ae52p-1,  0x1.47f1ccp-1,
    0x1.48353ep-1, 0x1.4878a6p-1,  0x1.48bc04p-1,  0x1.48ff5ap-1,
    0x1.4942a8p-1, 0x1.4985ecp-1,  0x1.49c928p-1,  0x1.4a0c5cp-1,
    0x1.4a4f86p-1, 0x1.4a92a8p-1,  0x1.4ad5cp-1,   0x1.4b18dp-1,
    0x1.4b5bd6p-1, 0x1.4b9ed4p-1,  0x1.4be1cap-1,  0x1.4c24b8p-1,
    0x1.4c679ap-1, 0x1.4caa76p-1,  0x1.4ced48p-1,  0x1.4d3012p-1,
    0x1.4d72d4p-1, 0x1.4db58cp-1,  0x1.4df83cp-1,  0x1.4e3ae2p-1,
    0x1.4e7d82p-1, 0x1.4ec016p-1,  0x1.4f02a4p-1,  0x1.4f4528p-1,
    0x1.4f87a4p-1, 0x1.4fca18p-1,  0x1.500c82p-1,  0x1.504ee4p-1,
    0x1.50913cp-1, 0x1.50d38ep-1,  0x1.5115d6p-1,  0x1.515816p-1,
    0x1.519a4cp-1, 0x1.51dc7ap-1,  0x1.521eap-1,   0x1.5260bep-1,
    0x1.52a2d2p-1, 0x1.52e4dep-1,  0x1.5326e2p-1,  0x1.5368dep-1,
    0x1.53aadp-1,  0x1.53ecbap-1,  0x1.542e9cp-1,  0x1.547076p-1,
    0x1.54b246p-1, 0x1.54f40ep-1,  0x1.5535cep-1,  0x1.557786p-1,
    0x1.55b936p-1, 0x1.55fadcp-1,  0x1.563c7ap-1,  0x1.567e1p-1,
    0x1.56bf9ep-1, 0x1.570122p-1,  0x1.5742ap-1,   0x1.578414p-1,
    0x1.57c58p-1,  0x1.5806e2p-1,  0x1.58483ep-1,  0x1.58899p-1,
    0x1.58cadcp-1, 0x1.590c1ep-1,  0x1.594d58p-1,  0x1.598e8ap-1,
    0x1.59cfb2p-1, 0x1.5a10d4p-1,  0x1.5a51ecp-1,  0x1.5a92fcp-1,
    0x1.5ad404p-1, 0x1.5b1504p-1,  0x1.5b55fcp-1,  0x1.5b96ecp-1,
    0x1.5bd7d4p-1, 0x1.5c18b2p-1,  0x1.5c5988p-1,  0x1.5c9a58p-1,
    0x1.5cdb1ep-1, 0x1.5d1bdcp-1,  0x1.5d5c92p-1,  0x1.5d9d4p-1,
    0x1.5ddde6p-1, 0x1.5e1e82p-1,  0x1.5e5f18p-1,  0x1.5e9fa6p-1,
    0x1.5ee02ap-1, 0x1.5f20a8p-1,  0x1.5f611cp-1,  0x1.5fa18ap-1,
    0x1.5fe1eep-1, 0x1.60224ap-1,  0x1.60629ep-1,  0x1.60a2eap-1,
    0x1.60e32ep-1, 0x1.61236cp-1,  0x1.6163ap-1,   0x1.61a3ccp-1,
    0x1.61e3eep-1, 0x1.62240cp-1,  0x1.62641ep-1,  0x1.62a42cp-1,
};

// For 0 <= i2 < 2^10, T3[i2] approximates log(x2), where x2 is
// subnormal and i2 equals the binary16 encoding of x2.
const logh_t3 = [1024]f32{
    @bitCast(@as(u32, 0xff80_0000)), -0x1.0a2b24p+4, -0x1.fe2804p+3, -0x1.f12e72p+3,
    -0x1.e7f9c2p+3,                  -0x1.e0d5c4p+3, -0x1.db003p+3,  -0x1.d61162p+3,
    -0x1.d1cb7ep+3,                  -0x1.ce069ep+3, -0x1.caa78p+3,  -0x1.c79abap+3,
    -0x1.c4d1ecp+3,                  -0x1.c24238p+3, -0x1.bfe32p+3,  -0x1.bdadfp+3,
    -0x1.bb9d3cp+3,                  -0x1.b9ac98p+3, -0x1.b7d85ap+3, -0x1.b61d7p+3,
    -0x1.b4793ep+3,                  -0x1.b2e98ep+3, -0x1.b16c76p+3, -0x1.b0005p+3,
    -0x1.aea3aap+3,                  -0x1.ad554p+3,  -0x1.ac13f4p+3, -0x1.aadecap+3,
    -0x1.a9b4dcp+3,                  -0x1.a89564p+3, -0x1.a77facp+3, -0x1.a6730ep+3,
    -0x1.a56ef8p+3,                  -0x1.a472e4p+3, -0x1.a37e56p+3, -0x1.a290dep+3,
    -0x1.a1aa18p+3,                  -0x1.a0c9a4p+3, -0x1.9fef2cp+3, -0x1.9f1a62p+3,
    -0x1.9e4afap+3,                  -0x1.9d80b2p+3, -0x1.9cbb4ap+3, -0x1.9bfa88p+3,
    -0x1.9b3e34p+3,                  -0x1.9a861ap+3, -0x1.99d20cp+3, -0x1.9921ep+3,
    -0x1.987566p+3,                  -0x1.97cc7ep+3, -0x1.9726fep+3, -0x1.9684c4p+3,
    -0x1.95e5b2p+3,                  -0x1.9549a6p+3, -0x1.94b086p+3, -0x1.941a36p+3,
    -0x1.93869ap+3,                  -0x1.92f59ap+3, -0x1.926722p+3, -0x1.91db18p+3,
    -0x1.91516ap+3,                  -0x1.90cap+3,   -0x1.9044ccp+3, -0x1.8fc1b8p+3,
    -0x1.8f40b6p+3,                  -0x1.8ec1b4p+3, -0x1.8e44a2p+3, -0x1.8dc97p+3,
    -0x1.8d5012p+3,                  -0x1.8cd87ap+3, -0x1.8c629cp+3, -0x1.8bee68p+3,
    -0x1.8b7bd4p+3,                  -0x1.8b0ad6p+3, -0x1.8a9b62p+3, -0x1.8a2d6cp+3,
    -0x1.89c0eap+3,                  -0x1.8955d4p+3, -0x1.88ec2p+3,  -0x1.8883c4p+3,
    -0x1.881cb8p+3,                  -0x1.87b6f4p+3, -0x1.87527p+3,  -0x1.86ef24p+3,
    -0x1.868d08p+3,                  -0x1.862c14p+3, -0x1.85cc44p+3, -0x1.856d9p+3,
    -0x1.850ffp+3,                   -0x1.84b36p+3,  -0x1.8457d8p+3, -0x1.83fd52p+3,
    -0x1.83a3cap+3,                  -0x1.834b3ap+3, -0x1.82f39cp+3, -0x1.829cecp+3,
    -0x1.824724p+3,                  -0x1.81f24p+3,  -0x1.819e3ap+3, -0x1.814b1p+3,
    -0x1.80f8bap+3,                  -0x1.80a736p+3, -0x1.80568p+3,  -0x1.800694p+3,
    -0x1.7fb76ep+3,                  -0x1.7f690ap+3, -0x1.7f1b64p+3, -0x1.7ece78p+3,
    -0x1.7e8244p+3,                  -0x1.7e36c2p+3, -0x1.7debf2p+3, -0x1.7da1dp+3,
    -0x1.7d5856p+3,                  -0x1.7d0f86p+3, -0x1.7cc758p+3, -0x1.7c7fccp+3,
    -0x1.7c38dep+3,                  -0x1.7bf28ep+3, -0x1.7bacd6p+3, -0x1.7b67b4p+3,
    -0x1.7b2326p+3,                  -0x1.7adf2ap+3, -0x1.7a9bbep+3, -0x1.7a58dep+3,
    -0x1.7a1688p+3,                  -0x1.79d4bcp+3, -0x1.799376p+3, -0x1.7952b4p+3,
    -0x1.791272p+3,                  -0x1.78d2b2p+3, -0x1.78937p+3,  -0x1.7854aap+3,
    -0x1.78165ep+3,                  -0x1.77d88ap+3, -0x1.779b2ep+3, -0x1.775e46p+3,
    -0x1.7721dp+3,                   -0x1.76e5ccp+3, -0x1.76aa38p+3, -0x1.766f12p+3,
    -0x1.763458p+3,                  -0x1.75fa0ap+3, -0x1.75c026p+3, -0x1.7586a8p+3,
    -0x1.754d92p+3,                  -0x1.7514ep+3,  -0x1.74dc94p+3, -0x1.74a4a8p+3,
    -0x1.746d1ep+3,                  -0x1.7435f4p+3, -0x1.73ff28p+3, -0x1.73c8bap+3,
    -0x1.7392a6p+3,                  -0x1.735ceep+3, -0x1.73279p+3,  -0x1.72f28ap+3,
    -0x1.72bddcp+3,                  -0x1.728984p+3, -0x1.72558p+3,  -0x1.7221d2p+3,
    -0x1.71ee74p+3,                  -0x1.71bb6ap+3, -0x1.7188b2p+3, -0x1.715648p+3,
    -0x1.71242cp+3,                  -0x1.70f26p+3,  -0x1.70c0ep+3,  -0x1.708facp+3,
    -0x1.705ec4p+3,                  -0x1.702e26p+3, -0x1.6ffdd2p+3, -0x1.6fcdc6p+3,
    -0x1.6f9e02p+3,                  -0x1.6f6e84p+3, -0x1.6f3f4cp+3, -0x1.6f105ap+3,
    -0x1.6ee1aep+3,                  -0x1.6eb344p+3, -0x1.6e851cp+3, -0x1.6e5738p+3,
    -0x1.6e2994p+3,                  -0x1.6dfc32p+3, -0x1.6dcf0ep+3, -0x1.6da22cp+3,
    -0x1.6d7586p+3,                  -0x1.6d492p+3,  -0x1.6d1cf6p+3, -0x1.6cf10ap+3,
    -0x1.6cc55ap+3,                  -0x1.6c99e4p+3, -0x1.6c6ea8p+3, -0x1.6c43a8p+3,
    -0x1.6c18ep+3,                   -0x1.6bee52p+3, -0x1.6bc3fcp+3, -0x1.6b99dep+3,
    -0x1.6b6ff8p+3,                  -0x1.6b4646p+3, -0x1.6b1cccp+3, -0x1.6af388p+3,
    -0x1.6aca78p+3,                  -0x1.6aa19cp+3, -0x1.6a78f4p+3, -0x1.6a508p+3,
    -0x1.6a283ep+3,                  -0x1.6a002ep+3, -0x1.69d852p+3, -0x1.69b0a6p+3,
    -0x1.69892cp+3,                  -0x1.6961ep+3,  -0x1.693ac6p+3, -0x1.6913dcp+3,
    -0x1.68ed2p+3,                   -0x1.68c694p+3, -0x1.68a034p+3, -0x1.687a04p+3,
    -0x1.6854p+3,                    -0x1.682e2ap+3, -0x1.68088p+3,  -0x1.67e302p+3,
    -0x1.67bdbp+3,                   -0x1.679888p+3, -0x1.67738cp+3, -0x1.674ebap+3,
    -0x1.672a14p+3,                  -0x1.670596p+3, -0x1.66e142p+3, -0x1.66bd18p+3,
    -0x1.669914p+3,                  -0x1.66753ap+3, -0x1.66518ap+3, -0x1.662dfep+3,
    -0x1.660a9cp+3,                  -0x1.65e76p+3,  -0x1.65c44ap+3, -0x1.65a15cp+3,
    -0x1.657e92p+3,                  -0x1.655beep+3, -0x1.65397p+3,  -0x1.651718p+3,
    -0x1.64f4e4p+3,                  -0x1.64d2d4p+3, -0x1.64b0e8p+3, -0x1.648f2p+3,
    -0x1.646d7ap+3,                  -0x1.644bfap+3, -0x1.642a9ap+3, -0x1.64096p+3,
    -0x1.63e846p+3,                  -0x1.63c74ep+3, -0x1.63a67ap+3, -0x1.6385c6p+3,
    -0x1.636532p+3,                  -0x1.6344c2p+3, -0x1.63247p+3,  -0x1.63044p+3,
    -0x1.62e43p+3,                   -0x1.62c44p+3,  -0x1.62a47p+3,  -0x1.6284bep+3,
    -0x1.62652ep+3,                  -0x1.6245bap+3, -0x1.622668p+3, -0x1.620732p+3,
    -0x1.61e81cp+3,                  -0x1.61c922p+3, -0x1.61aa48p+3, -0x1.618b8ap+3,
    -0x1.616ceap+3,                  -0x1.614e68p+3, -0x1.613002p+3, -0x1.6111bap+3,
    -0x1.60f38cp+3,                  -0x1.60d57cp+3, -0x1.60b78ap+3, -0x1.6099b2p+3,
    -0x1.607bf4p+3,                  -0x1.605e54p+3, -0x1.6040dp+3,  -0x1.602364p+3,
    -0x1.600616p+3,                  -0x1.5fe8e2p+3, -0x1.5fcbc8p+3, -0x1.5faec8p+3,
    -0x1.5f91e2p+3,                  -0x1.5f7516p+3, -0x1.5f5866p+3, -0x1.5f3bcep+3,
    -0x1.5f1f5p+3,                   -0x1.5f02eap+3, -0x1.5ee69ep+3, -0x1.5eca6ap+3,
    -0x1.5eae5p+3,                   -0x1.5e924ep+3, -0x1.5e7666p+3, -0x1.5e5a94p+3,
    -0x1.5e3edcp+3,                  -0x1.5e233ap+3, -0x1.5e07bp+3,  -0x1.5dec4p+3,
    -0x1.5dd0e6p+3,                  -0x1.5db5a2p+3, -0x1.5d9a76p+3, -0x1.5d7f62p+3,
    -0x1.5d6464p+3,                  -0x1.5d497cp+3, -0x1.5d2eacp+3, -0x1.5d13f2p+3,
    -0x1.5cf94ep+3,                  -0x1.5cdecp+3,  -0x1.5cc448p+3, -0x1.5ca9e6p+3,
    -0x1.5c8f9ap+3,                  -0x1.5c7562p+3, -0x1.5c5b4p+3,  -0x1.5c4134p+3,
    -0x1.5c273ep+3,                  -0x1.5c0d5cp+3, -0x1.5bf38ep+3, -0x1.5bd9d6p+3,
    -0x1.5bc032p+3,                  -0x1.5ba6a2p+3, -0x1.5b8d28p+3, -0x1.5b73cp+3,
    -0x1.5b5a6ep+3,                  -0x1.5b413p+3,  -0x1.5b2804p+3, -0x1.5b0eeep+3,
    -0x1.5af5eap+3,                  -0x1.5adcfap+3, -0x1.5ac41ep+3, -0x1.5aab54p+3,
    -0x1.5a929ep+3,                  -0x1.5a79fap+3, -0x1.5a616ap+3, -0x1.5a48ecp+3,
    -0x1.5a3082p+3,                  -0x1.5a182ap+3, -0x1.59ffe4p+3, -0x1.59e7bp+3,
    -0x1.59cf8ep+3,                  -0x1.59b78p+3,  -0x1.599f82p+3, -0x1.598798p+3,
    -0x1.596fbep+3,                  -0x1.5957f8p+3, -0x1.594042p+3, -0x1.59289cp+3,
    -0x1.59110ap+3,                  -0x1.58f988p+3, -0x1.58e218p+3, -0x1.58cab8p+3,
    -0x1.58b36ap+3,                  -0x1.589c2cp+3, -0x1.5885p+3,   -0x1.586de4p+3,
    -0x1.5856dap+3,                  -0x1.583fdep+3, -0x1.5828f4p+3, -0x1.58121ap+3,
    -0x1.57fb52p+3,                  -0x1.57e498p+3, -0x1.57cdeep+3, -0x1.57b756p+3,
    -0x1.57a0ccp+3,                  -0x1.578a52p+3, -0x1.5773e8p+3, -0x1.575d8ep+3,
    -0x1.574744p+3,                  -0x1.573108p+3, -0x1.571adep+3, -0x1.5704cp+3,
    -0x1.56eeb4p+3,                  -0x1.56d8b6p+3, -0x1.56c2c8p+3, -0x1.56ace8p+3,
    -0x1.569716p+3,                  -0x1.568154p+3, -0x1.566bap+3,  -0x1.5655fcp+3,
    -0x1.564066p+3,                  -0x1.562adep+3, -0x1.561566p+3, -0x1.55fffap+3,
    -0x1.55ea9ep+3,                  -0x1.55d55p+3,  -0x1.55c01p+3,  -0x1.55aadep+3,
    -0x1.5595bap+3,                  -0x1.5580a4p+3, -0x1.556b9cp+3, -0x1.5556ap+3,
    -0x1.5541b4p+3,                  -0x1.552cd6p+3, -0x1.551804p+3, -0x1.55034p+3,
    -0x1.54ee8ap+3,                  -0x1.54d9ep+3,  -0x1.54c544p+3, -0x1.54b0b6p+3,
    -0x1.549c34p+3,                  -0x1.5487cp+3,  -0x1.547358p+3, -0x1.545efep+3,
    -0x1.544abp+3,                   -0x1.54367p+3,  -0x1.54223cp+3, -0x1.540e16p+3,
    -0x1.53f9fap+3,                  -0x1.53e5eep+3, -0x1.53d1ecp+3, -0x1.53bdf8p+3,
    -0x1.53aa0ep+3,                  -0x1.539632p+3, -0x1.538264p+3, -0x1.536eap+3,
    -0x1.535ae8p+3,                  -0x1.53473ep+3, -0x1.53339ep+3, -0x1.53200ap+3,
    -0x1.530c84p+3,                  -0x1.52f908p+3, -0x1.52e598p+3, -0x1.52d236p+3,
    -0x1.52bedep+3,                  -0x1.52ab9p+3,  -0x1.52985p+3,  -0x1.52851cp+3,
    -0x1.5271f2p+3,                  -0x1.525ed4p+3, -0x1.524bcp+3,  -0x1.5238bap+3,
    -0x1.5225bep+3,                  -0x1.5212ccp+3, -0x1.51ffe6p+3, -0x1.51ed0cp+3,
    -0x1.51da3cp+3,                  -0x1.51c778p+3, -0x1.51b4bep+3, -0x1.51a21p+3,
    -0x1.518f6cp+3,                  -0x1.517cd4p+3, -0x1.516a46p+3, -0x1.5157c2p+3,
    -0x1.51454ap+3,                  -0x1.5132dcp+3, -0x1.512078p+3, -0x1.510e1ep+3,
    -0x1.50fbdp+3,                   -0x1.50e98cp+3, -0x1.50d754p+3, -0x1.50c524p+3,
    -0x1.50b3p+3,                    -0x1.50a0e4p+3, -0x1.508ed4p+3, -0x1.507ccep+3,
    -0x1.506ad2p+3,                  -0x1.5058ep+3,  -0x1.5046f8p+3, -0x1.50351ap+3,
    -0x1.502346p+3,                  -0x1.50117cp+3, -0x1.4fffbcp+3, -0x1.4fee06p+3,
    -0x1.4fdc58p+3,                  -0x1.4fcab6p+3, -0x1.4fb91cp+3, -0x1.4fa78ep+3,
    -0x1.4f9608p+3,                  -0x1.4f848cp+3, -0x1.4f7318p+3, -0x1.4f61aep+3,
    -0x1.4f505p+3,                   -0x1.4f3ef8p+3, -0x1.4f2dacp+3, -0x1.4f1c68p+3,
    -0x1.4f0b2ep+3,                  -0x1.4ef9fcp+3, -0x1.4ee8d4p+3, -0x1.4ed7b6p+3,
    -0x1.4ec6ap+3,                   -0x1.4eb594p+3, -0x1.4ea49p+3,  -0x1.4e9396p+3,
    -0x1.4e82a4p+3,                  -0x1.4e71bcp+3, -0x1.4e60dcp+3, -0x1.4e5006p+3,
    -0x1.4e3f38p+3,                  -0x1.4e2e72p+3, -0x1.4e1db6p+3, -0x1.4e0d02p+3,
    -0x1.4dfc58p+3,                  -0x1.4debb6p+3, -0x1.4ddb1cp+3, -0x1.4dca8cp+3,
    -0x1.4dba02p+3,                  -0x1.4da982p+3, -0x1.4d990cp+3, -0x1.4d889cp+3,
    -0x1.4d7836p+3,                  -0x1.4d67d8p+3, -0x1.4d5782p+3, -0x1.4d4734p+3,
    -0x1.4d36fp+3,                   -0x1.4d26b2p+3, -0x1.4d167ep+3, -0x1.4d0652p+3,
    -0x1.4cf62ep+3,                  -0x1.4ce612p+3, -0x1.4cd5fcp+3, -0x1.4cc5fp+3,
    -0x1.4cb5ecp+3,                  -0x1.4ca5fp+3,  -0x1.4c95fcp+3, -0x1.4c861p+3,
    -0x1.4c762cp+3,                  -0x1.4c665p+3,  -0x1.4c567cp+3, -0x1.4c46bp+3,
    -0x1.4c36eap+3,                  -0x1.4c272ep+3, -0x1.4c1778p+3, -0x1.4c07cap+3,
    -0x1.4bf824p+3,                  -0x1.4be886p+3, -0x1.4bd8eep+3, -0x1.4bc96p+3,
    -0x1.4bb9d8p+3,                  -0x1.4baa58p+3, -0x1.4b9aep+3,  -0x1.4b8b6ep+3,
    -0x1.4b7c04p+3,                  -0x1.4b6ca2p+3, -0x1.4b5d48p+3, -0x1.4b4df4p+3,
    -0x1.4b3ea8p+3,                  -0x1.4b2f62p+3, -0x1.4b2024p+3, -0x1.4b10eep+3,
    -0x1.4b01cp+3,                   -0x1.4af298p+3, -0x1.4ae376p+3, -0x1.4ad45cp+3,
    -0x1.4ac54ap+3,                  -0x1.4ab63ep+3, -0x1.4aa73ap+3, -0x1.4a983cp+3,
    -0x1.4a8946p+3,                  -0x1.4a7a56p+3, -0x1.4a6b6ep+3, -0x1.4a5c8cp+3,
    -0x1.4a4db2p+3,                  -0x1.4a3edep+3, -0x1.4a3012p+3, -0x1.4a214cp+3,
    -0x1.4a128cp+3,                  -0x1.4a03d4p+3, -0x1.49f522p+3, -0x1.49e676p+3,
    -0x1.49d7d2p+3,                  -0x1.49c936p+3, -0x1.49ba9ep+3, -0x1.49ac0ep+3,
    -0x1.499d84p+3,                  -0x1.498f02p+3, -0x1.498084p+3, -0x1.49720ep+3,
    -0x1.4963ap+3,                   -0x1.495536p+3, -0x1.4946d4p+3, -0x1.493878p+3,
    -0x1.492a22p+3,                  -0x1.491bd4p+3, -0x1.490d8ap+3, -0x1.48ff48p+3,
    -0x1.48f10cp+3,                  -0x1.48e2d6p+3, -0x1.48d4a6p+3, -0x1.48c67ep+3,
    -0x1.48b85ap+3,                  -0x1.48aa3ep+3, -0x1.489c28p+3, -0x1.488e18p+3,
    -0x1.48800ep+3,                  -0x1.48720ap+3, -0x1.48640cp+3, -0x1.485614p+3,
    -0x1.484822p+3,                  -0x1.483a36p+3, -0x1.482c52p+3, -0x1.481e72p+3,
    -0x1.481098p+3,                  -0x1.4802c4p+3, -0x1.47f4f8p+3, -0x1.47e73p+3,
    -0x1.47d96ep+3,                  -0x1.47cbb2p+3, -0x1.47bdfcp+3, -0x1.47b04cp+3,
    -0x1.47a2a2p+3,                  -0x1.4794fep+3, -0x1.47876p+3,  -0x1.4779c6p+3,
    -0x1.476c34p+3,                  -0x1.475ea6p+3, -0x1.47511ep+3, -0x1.47439cp+3,
    -0x1.47362p+3,                   -0x1.4728aap+3, -0x1.471b3ap+3, -0x1.470dcep+3,
    -0x1.470068p+3,                  -0x1.46f30ap+3, -0x1.46e5aep+3, -0x1.46d85ap+3,
    -0x1.46cb0ap+3,                  -0x1.46bdcp+3,  -0x1.46b07cp+3, -0x1.46a33ep+3,
    -0x1.469604p+3,                  -0x1.4688d2p+3, -0x1.467ba2p+3, -0x1.466e7ap+3,
    -0x1.466156p+3,                  -0x1.465438p+3, -0x1.46472p+3,  -0x1.463a0cp+3,
    -0x1.462cfep+3,                  -0x1.461ff6p+3, -0x1.4612f2p+3, -0x1.4605f4p+3,
    -0x1.45f8fap+3,                  -0x1.45ec06p+3, -0x1.45df18p+3, -0x1.45d23p+3,
    -0x1.45c54cp+3,                  -0x1.45b86cp+3, -0x1.45ab92p+3, -0x1.459ebep+3,
    -0x1.4591eep+3,                  -0x1.458524p+3, -0x1.45786p+3,  -0x1.456bap+3,
    -0x1.455ee4p+3,                  -0x1.45522ep+3, -0x1.45457ep+3, -0x1.4538d2p+3,
    -0x1.452c2cp+3,                  -0x1.451f8ap+3, -0x1.4512ecp+3, -0x1.450654p+3,
    -0x1.44f9c2p+3,                  -0x1.44ed34p+3, -0x1.44e0aap+3, -0x1.44d426p+3,
    -0x1.44c7a6p+3,                  -0x1.44bb2cp+3, -0x1.44aeb6p+3, -0x1.44a246p+3,
    -0x1.4495dap+3,                  -0x1.448974p+3, -0x1.447d1p+3,  -0x1.4470b4p+3,
    -0x1.44645ap+3,                  -0x1.445806p+3, -0x1.444bb8p+3, -0x1.443f6cp+3,
    -0x1.443326p+3,                  -0x1.4426e6p+3, -0x1.441aaap+3, -0x1.440e72p+3,
    -0x1.44023ep+3,                  -0x1.43f61p+3,  -0x1.43e9e6p+3, -0x1.43ddc2p+3,
    -0x1.43d1ap+3,                   -0x1.43c584p+3, -0x1.43b96ep+3, -0x1.43ad5ap+3,
    -0x1.43a14cp+3,                  -0x1.439542p+3, -0x1.43893cp+3, -0x1.437d3cp+3,
    -0x1.43714p+3,                   -0x1.436548p+3, -0x1.435954p+3, -0x1.434d66p+3,
    -0x1.43417cp+3,                  -0x1.433596p+3, -0x1.4329b4p+3, -0x1.431dd6p+3,
    -0x1.4311fep+3,                  -0x1.43062ap+3, -0x1.42fa5ap+3, -0x1.42ee8ep+3,
    -0x1.42e2c6p+3,                  -0x1.42d704p+3, -0x1.42cb46p+3, -0x1.42bf8ap+3,
    -0x1.42b3d4p+3,                  -0x1.42a822p+3, -0x1.429c76p+3, -0x1.4290ccp+3,
    -0x1.428528p+3,                  -0x1.427986p+3, -0x1.426deap+3, -0x1.426252p+3,
    -0x1.4256bep+3,                  -0x1.424b2ep+3, -0x1.423fa2p+3, -0x1.42341ap+3,
    -0x1.422896p+3,                  -0x1.421d16p+3, -0x1.42119cp+3, -0x1.420624p+3,
    -0x1.41fab2p+3,                  -0x1.41ef42p+3, -0x1.41e3d8p+3, -0x1.41d87p+3,
    -0x1.41cd0ep+3,                  -0x1.41c1bp+3,  -0x1.41b654p+3, -0x1.41aafep+3,
    -0x1.419facp+3,                  -0x1.41945cp+3, -0x1.418912p+3, -0x1.417dccp+3,
    -0x1.417288p+3,                  -0x1.41674ap+3, -0x1.415c1p+3,  -0x1.4150d8p+3,
    -0x1.4145a6p+3,                  -0x1.413a76p+3, -0x1.412f4cp+3, -0x1.412424p+3,
    -0x1.4119p+3,                    -0x1.410de2p+3, -0x1.4102c6p+3, -0x1.40f7aep+3,
    -0x1.40ec9ap+3,                  -0x1.40e18ap+3, -0x1.40d67ep+3, -0x1.40cb76p+3,
    -0x1.40c07p+3,                   -0x1.40b57p+3,  -0x1.40aa72p+3, -0x1.409f7ap+3,
    -0x1.409484p+3,                  -0x1.408992p+3, -0x1.407ea4p+3, -0x1.4073bap+3,
    -0x1.4068d4p+3,                  -0x1.405dfp+3,  -0x1.405312p+3, -0x1.404836p+3,
    -0x1.403d5ep+3,                  -0x1.40328ap+3, -0x1.4027b8p+3, -0x1.401cecp+3,
    -0x1.401222p+3,                  -0x1.40075ep+3, -0x1.3ffc9cp+3, -0x1.3ff1dcp+3,
    -0x1.3fe722p+3,                  -0x1.3fdc6cp+3, -0x1.3fd1b8p+3, -0x1.3fc708p+3,
    -0x1.3fbc5ap+3,                  -0x1.3fb1b2p+3, -0x1.3fa70cp+3, -0x1.3f9c6ap+3,
    -0x1.3f91ccp+3,                  -0x1.3f8732p+3, -0x1.3f7c9ap+3, -0x1.3f7206p+3,
    -0x1.3f6776p+3,                  -0x1.3f5ceap+3, -0x1.3f526p+3,  -0x1.3f47dap+3,
    -0x1.3f3d58p+3,                  -0x1.3f32dap+3, -0x1.3f285ep+3, -0x1.3f1de6p+3,
    -0x1.3f1372p+3,                  -0x1.3f09p+3,   -0x1.3efe92p+3, -0x1.3ef428p+3,
    -0x1.3ee9cp+3,                   -0x1.3edf5cp+3, -0x1.3ed4fcp+3, -0x1.3ecaap+3,
    -0x1.3ec046p+3,                  -0x1.3eb5fp+3,  -0x1.3eab9cp+3, -0x1.3ea14ep+3,
    -0x1.3e9702p+3,                  -0x1.3e8cb8p+3, -0x1.3e8272p+3, -0x1.3e783p+3,
    -0x1.3e6df2p+3,                  -0x1.3e63b6p+3, -0x1.3e597cp+3, -0x1.3e4f48p+3,
    -0x1.3e4516p+3,                  -0x1.3e3ae6p+3, -0x1.3e30bcp+3, -0x1.3e2692p+3,
    -0x1.3e1c6ep+3,                  -0x1.3e124cp+3, -0x1.3e082ep+3, -0x1.3dfe12p+3,
    -0x1.3df3fap+3,                  -0x1.3de9e4p+3, -0x1.3ddfd2p+3, -0x1.3dd5c4p+3,
    -0x1.3dcbb8p+3,                  -0x1.3dc1bp+3,  -0x1.3db7aap+3, -0x1.3dada8p+3,
    -0x1.3da3a8p+3,                  -0x1.3d99acp+3, -0x1.3d8fb4p+3, -0x1.3d85bep+3,
    -0x1.3d7bccp+3,                  -0x1.3d71dcp+3, -0x1.3d67fp+3,  -0x1.3d5e06p+3,
    -0x1.3d542p+3,                   -0x1.3d4a3cp+3, -0x1.3d405cp+3, -0x1.3d368p+3,
    -0x1.3d2ca6p+3,                  -0x1.3d22cep+3, -0x1.3d18fap+3, -0x1.3d0f2ap+3,
    -0x1.3d055ap+3,                  -0x1.3cfb9p+3,  -0x1.3cf1c8p+3, -0x1.3ce802p+3,
    -0x1.3cde4p+3,                   -0x1.3cd482p+3, -0x1.3ccac6p+3, -0x1.3cc10cp+3,
    -0x1.3cb756p+3,                  -0x1.3cada2p+3, -0x1.3ca3f2p+3, -0x1.3c9a44p+3,
    -0x1.3c909ap+3,                  -0x1.3c86f2p+3, -0x1.3c7d4ep+3, -0x1.3c73acp+3,
    -0x1.3c6a0ep+3,                  -0x1.3c6072p+3, -0x1.3c56d8p+3, -0x1.3c4d42p+3,
    -0x1.3c43aep+3,                  -0x1.3c3a1ep+3, -0x1.3c309p+3,  -0x1.3c2706p+3,
    -0x1.3c1d7ep+3,                  -0x1.3c13f8p+3, -0x1.3c0a76p+3, -0x1.3c00f6p+3,
    -0x1.3bf77ap+3,                  -0x1.3beep+3,   -0x1.3be48ap+3, -0x1.3bdb14p+3,
    -0x1.3bd1a4p+3,                  -0x1.3bc834p+3, -0x1.3bbec8p+3, -0x1.3bb56p+3,
    -0x1.3babfap+3,                  -0x1.3ba296p+3, -0x1.3b9934p+3, -0x1.3b8fd6p+3,
    -0x1.3b867cp+3,                  -0x1.3b7d22p+3, -0x1.3b73ccp+3, -0x1.3b6a7ap+3,
    -0x1.3b612ap+3,                  -0x1.3b57dcp+3, -0x1.3b4e9p+3,  -0x1.3b4548p+3,
    -0x1.3b3c02p+3,                  -0x1.3b32cp+3,  -0x1.3b297ep+3, -0x1.3b2042p+3,
    -0x1.3b1706p+3,                  -0x1.3b0dcep+3, -0x1.3b0498p+3, -0x1.3afb66p+3,
    -0x1.3af234p+3,                  -0x1.3ae908p+3, -0x1.3adfdcp+3, -0x1.3ad6b4p+3,
    -0x1.3acd8ep+3,                  -0x1.3ac46ap+3, -0x1.3abb4ap+3, -0x1.3ab22cp+3,
    -0x1.3aa91p+3,                   -0x1.3a9ff8p+3, -0x1.3a96ep+3,  -0x1.3a8dcep+3,
    -0x1.3a84bcp+3,                  -0x1.3a7baep+3, -0x1.3a72a2p+3, -0x1.3a6998p+3,
    -0x1.3a6092p+3,                  -0x1.3a578cp+3, -0x1.3a4e8ap+3, -0x1.3a458cp+3,
    -0x1.3a3c8ep+3,                  -0x1.3a3394p+3, -0x1.3a2a9cp+3, -0x1.3a21a8p+3,
    -0x1.3a18b4p+3,                  -0x1.3a0fc4p+3, -0x1.3a06d8p+3, -0x1.39fdecp+3,
    -0x1.39f504p+3,                  -0x1.39ec1cp+3, -0x1.39e338p+3, -0x1.39da58p+3,
    -0x1.39d178p+3,                  -0x1.39c89cp+3, -0x1.39bfc2p+3, -0x1.39b6eap+3,
    -0x1.39ae16p+3,                  -0x1.39a544p+3, -0x1.399c72p+3, -0x1.3993a6p+3,
    -0x1.398adap+3,                  -0x1.39821p+3,  -0x1.39794ap+3, -0x1.397086p+3,
    -0x1.3967c4p+3,                  -0x1.395f04p+3, -0x1.395648p+3, -0x1.394d8ep+3,
    -0x1.3944d6p+3,                  -0x1.393c2p+3,  -0x1.39336cp+3, -0x1.392abap+3,
    -0x1.39220cp+3,                  -0x1.39196p+3,  -0x1.3910b6p+3, -0x1.39080ep+3,
    -0x1.38ff68p+3,                  -0x1.38f6c6p+3, -0x1.38ee24p+3, -0x1.38e586p+3,
    -0x1.38dceap+3,                  -0x1.38d45p+3,  -0x1.38cbbap+3, -0x1.38c324p+3,
    -0x1.38ba92p+3,                  -0x1.38b2p+3,   -0x1.38a972p+3, -0x1.38a0e6p+3,
    -0x1.38985ep+3,                  -0x1.388fd6p+3, -0x1.38875p+3,  -0x1.387ecep+3,
    -0x1.38764ep+3,                  -0x1.386dcep+3, -0x1.386552p+3, -0x1.385cd8p+3,
    -0x1.385462p+3,                  -0x1.384becp+3, -0x1.384378p+3, -0x1.383b08p+3,
    -0x1.38329ap+3,                  -0x1.382a2cp+3, -0x1.3821c2p+3, -0x1.38195ap+3,
    -0x1.3810f4p+3,                  -0x1.38089p+3,  -0x1.38003p+3,  -0x1.37f7dp+3,
    -0x1.37ef74p+3,                  -0x1.37e718p+3, -0x1.37decp+3,  -0x1.37d66ap+3,
    -0x1.37ce14p+3,                  -0x1.37c5c2p+3, -0x1.37bd72p+3, -0x1.37b524p+3,
    -0x1.37acdap+3,                  -0x1.37a49p+3,  -0x1.379c48p+3, -0x1.379402p+3,
    -0x1.378bcp+3,                   -0x1.37837ep+3, -0x1.377b4p+3,  -0x1.377304p+3,
    -0x1.376ac8p+3,                  -0x1.37629p+3,  -0x1.375a5ap+3, -0x1.375226p+3,
    -0x1.3749f4p+3,                  -0x1.3741c4p+3, -0x1.373994p+3, -0x1.37316ap+3,
    -0x1.37294p+3,                   -0x1.372118p+3, -0x1.3718f2p+3, -0x1.3710cep+3,
    -0x1.3708acp+3,                  -0x1.37008ep+3, -0x1.36f87p+3,  -0x1.36f054p+3,
    -0x1.36e83cp+3,                  -0x1.36e024p+3, -0x1.36d80ep+3, -0x1.36cffcp+3,
    -0x1.36c7eap+3,                  -0x1.36bfdcp+3, -0x1.36b7cep+3, -0x1.36afc4p+3,
    -0x1.36a7bap+3,                  -0x1.369fb2p+3, -0x1.3697aep+3, -0x1.368faap+3,
};

pub fn logh(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude > infinity_bits_h) {
        return @bitCast(bits | logh_quiet_nan_bit);
    }
    if (magnitude == 0) return @bitCast(logh_minus_infinity_bits);
    if ((bits & sign_mask_h) != 0) return @bitCast(logh_canonical_nan_bits);
    if (magnitude == infinity_bits_h) return value;

    const exponent_index = bits >> 10;
    const significand_index = bits & 0x03ff;
    if (exponent_index == 0) return @floatCast(logh_t3[significand_index]);
    return @floatCast(logh_t1[exponent_index] + logh_t2[significand_index]);
}
// Correctly rounded sine for binary16 values.
//
// Copyright (c) 2025 Maxence Ponsardin and Paul Zimmermann.
//
// This code is derived from CORE-MATH src/binary16/sin/sinf16.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
const sinh_data = @import("scalar_sinh_data.zig");

pub fn sinh(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude > infinity_bits_h) {
        const quiet_nan_bit_h: u16 = 0x0200;
        return @bitCast(bits | quiet_nan_bit_h);
    }
    if (magnitude == infinity_bits_h) return @bitCast(@as(u16, 0x7e00));

    // These are the two remaining hard-to-round inputs after optimizing the
    // four lookup tables. Each case also covers the input with the other sign.
    if (magnitude == 0x658c) {
        const result: f32 = if ((bits & sign_mask_h) == 0)
            0x1.f9bd02p-14
        else
            -0x1.f9bd02p-14;
        return @floatCast(result);
    }
    if (magnitude == 0x6429) {
        const result: f32 = if ((bits & sign_mask_h) == 0)
            -0x1.7b4dc2p-14
        else
            0x1.7b4dc2p-14;
        return @floatCast(result);
    }

    // Decompose x into x1 + x2 and evaluate
    // sin(x1) * cos(x2) + cos(x1) * sin(x2).
    const table1_index: usize = bits >> 5;
    const table2_index: usize = ((bits >> 10) << 5) | (bits & 0x001f);
    const sin_x1: f32 = @bitCast(sinh_data.s1[table1_index]);
    const cos_x1: f32 = @bitCast(sinh_data.c1[table1_index]);
    const sin_x2: f32 = @bitCast(sinh_data.s2[table2_index]);
    const cos_x2: f32 = @bitCast(sinh_data.c2[table2_index]);
    const result = @mulAdd(f32, sin_x1, cos_x2, cos_x1 * sin_x2);
    return @floatCast(result);
}
// Correctly rounded cosine for binary16 values.
//
// Copyright (c) 2025 Maxence Ponsardin and Paul Zimmermann.
//
// This code is derived from CORE-MATH src/binary16/cos/cosf16.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
const cosh_data = @import("scalar_cosh_data.zig");

pub fn cosh(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;

    if (magnitude > infinity_bits_h) {
        const quiet_nan_bit_h: u16 = 0x0200;
        return @bitCast(bits | quiet_nan_bit_h);
    }
    if (magnitude == infinity_bits_h) return @bitCast(@as(u16, 0x7e00));
    if (magnitude == 0) return 1.0;

    // Decompose x into x1 + x2 and evaluate
    // cos(x1) * cos(x2) - sin(x1) * sin(x2).
    const table1_index: usize = bits >> 5;
    const table2_index: usize = ((bits >> 10) << 5) | (bits & 0x001f);
    const sin_x1: f32 = @bitCast(cosh_data.s1[table1_index]);
    const cos_x1: f32 = @bitCast(cosh_data.c1[table1_index]);
    const minus_sin_x2: f32 = @bitCast(cosh_data.s2[table2_index]);
    const cos_x2: f32 = @bitCast(cosh_data.c2[table2_index]);
    const result = @mulAdd(f32, cos_x1, cos_x2, sin_x1 * minus_sin_x2);
    return @floatCast(result);
}

// Correctly rounded natural logarithm for binary32 values.
//
// Copyright (c) 2023-2026 Alexei Sibidanov and Paul Zimmermann.
//
// This code is derived from CORE-MATH src/binary32/log/logf.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
pub fn logf(value: f32) f32 {
    @setFloatMode(.strict);

    const input_bits: u32 = @bitCast(value);
    const magnitude = input_bits & magnitude_mask_f;

    if (magnitude == 0) return @bitCast(logf_minus_infinity_bits);
    if (magnitude > infinity_bits_f) {
        return @bitCast(input_bits | logf_quiet_nan_bit);
    }
    if (magnitude == infinity_bits_f) {
        if ((input_bits & sign_mask_f) == 0) return value;
        return @bitCast(logf_canonical_nan_bits);
    }
    if ((input_bits & sign_mask_f) != 0) return @bitCast(logf_canonical_nan_bits);

    var bits = input_bits;
    if (bits < 0x0080_0000) {
        const shift: u32 = @clz(bits) - 8;
        bits <<= @intCast(shift);
        bits -%= shift << 23;
    }
    if (bits == 127 << 23) return 0.0;

    const mantissa = bits & 0x007f_ffff;
    const table_index: usize = @intCast((mantissa + 0x0001_0000) >> 17);
    const normalized_bits: i32 = @bitCast(bits);
    const exponent: i32 = (normalized_bits >> 23) - 127;
    const exponent_wide: f64 = @floatFromInt(exponent);
    const tz_bits: u64 =
        (@as(u64, mantissa) | (@as(u64, 1023) << 23)) << 29;
    const tz: f64 = @bitCast(tz_bits);
    const z = tz * logf_data.tr[table_index] - 1.0;
    const z_squared = z * z;

    var result = ((exponent_wide * 0x1.62e42fefa39efp-1 + logf_data.tl[table_index]) +
        z * logf_fast_coefficients[0]) +
        z_squared * (logf_fast_coefficients[1] + z * logf_fast_coefficients[2]);
    var upper: f32 = @floatCast(result);
    const lower: f32 = @floatCast(result + 0x1.f06p-33);
    if (upper == lower) return upper;

    var accurate = z_squared *
        ((logf_accurate_coefficients[0] + z * logf_accurate_coefficients[1]) +
            z_squared *
                ((logf_accurate_coefficients[2] + z * logf_accurate_coefficients[3]) +
                    z_squared *
                        (logf_accurate_coefficients[4] + z * logf_accurate_coefficients[5] +
                            z_squared * logf_accurate_coefficients[6])));
    const delta_bits: u32 = @bitCast(value - 1.0);
    if ((delta_bits & magnitude_mask_f) < 0x3a80_0000) {
        return @floatCast(z + accurate);
    }

    accurate -= 0x1.0ca86c3898dp-49 * exponent_wide;
    accurate += z;
    accurate += logf_data.tl[table_index] - logf_data.tl[0];
    const exponent_log: f64 = exponent_wide * 0x1.62e42fefa3ap-1;
    result = exponent_log + accurate;
    upper = @floatCast(result);
    const result_bits: u64 = @bitCast(result);
    if ((result_bits & 0x0fff_ffff) == 0) {
        const rounding_remainder = (exponent_log - result) + accurate;
        result += rounding_remainder * 64.0;
        upper = @floatCast(result);
    }
    return upper;
}

pub fn absh(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    return @bitCast(bits & magnitude_mask_h);
}

pub fn absf(value: f32) f32 {
    const bits: u32 = @bitCast(value);
    return @bitCast(bits & magnitude_mask_f);
}

pub fn absd(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    return @bitCast(bits & magnitude_mask_d);
}

pub fn ulpDistanceh(a: f16, b: f16) u16 {
    const a_bits: u16 = @bitCast(a);
    const b_bits: u16 = @bitCast(b);

    if (isNanH(a_bits) or isNanH(b_bits)) return ~@as(u16, 0);

    const a_ordered = orderedBitsH(a_bits);
    const b_ordered = orderedBitsH(b_bits);
    if (a_ordered >= b_ordered) return a_ordered - b_ordered;
    return b_ordered - a_ordered;
}

pub fn ulpDistancef(a: f32, b: f32) u32 {
    const a_bits: u32 = @bitCast(a);
    const b_bits: u32 = @bitCast(b);

    if (isNanF(a_bits) or isNanF(b_bits)) return ~@as(u32, 0);

    const a_ordered = orderedBitsF(a_bits);
    const b_ordered = orderedBitsF(b_bits);
    if (a_ordered >= b_ordered) return a_ordered - b_ordered;
    return b_ordered - a_ordered;
}

pub fn ulpDistanced(a: f64, b: f64) u64 {
    const a_bits: u64 = @bitCast(a);
    const b_bits: u64 = @bitCast(b);

    if (isNanD(a_bits) or isNanD(b_bits)) return ~@as(u64, 0);

    const a_ordered = orderedBitsD(a_bits);
    const b_ordered = orderedBitsD(b_bits);
    if (a_ordered >= b_ordered) return a_ordered - b_ordered;
    return b_ordered - a_ordered;
}

fn orderedBitsH(bits: u16) u16 {
    if ((bits & sign_mask_h) != 0) return (~bits) + 1;
    return bits | sign_mask_h;
}

fn orderedBitsF(bits: u32) u32 {
    if ((bits & sign_mask_f) != 0) return (~bits) + 1;
    return bits | sign_mask_f;
}

fn orderedBitsD(bits: u64) u64 {
    if ((bits & sign_mask_d) != 0) return (~bits) + 1;
    return bits | sign_mask_d;
}

fn isNanH(bits: u16) bool {
    return (bits & magnitude_mask_h) > infinity_bits_h;
}

fn isNanF(bits: u32) bool {
    return (bits & magnitude_mask_f) > infinity_bits_f;
}

fn isNanD(bits: u64) bool {
    return (bits & magnitude_mask_d) > infinity_bits_d;
}

// Correctly rounded natural exponential for binary32 values.
//
// Copyright (c) 2023-2025 Alexei Sibidanov.
//
// This code is derived from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
const expf_quiet_nan_bit: u32 = 0x0040_0000;

const expf_fast_coefficients = [4]f64{
    1.0,
    0x1.62e42fef4c4e7p-1,
    0x1.ebfd1b232f475p-3,
    0x1.c6b19384ecd93p-5,
};

const expf_accurate_coefficients = [6]f64{
    0x1.62e42fefa39efp-1,
    0x1.ebfbdff82c58fp-3,
    0x1.c6b08d702e0edp-5,
    0x1.3b2ab6fb92e5ep-7,
    0x1.5d886e6d54203p-10,
    0x1.430976b8ce6efp-13,
};

const expf_scale_bits = [64]u64{
    0x3ff0000000000000, 0x3ff02c9a3e778061, 0x3ff059b0d3158574,
    0x3ff0874518759bc8, 0x3ff0b5586cf9890f, 0x3ff0e3ec32d3d1a2,
    0x3ff11301d0125b51, 0x3ff1429aaea92de0, 0x3ff172b83c7d517b,
    0x3ff1a35beb6fcb75, 0x3ff1d4873168b9aa, 0x3ff2063b88628cd6,
    0x3ff2387a6e756238, 0x3ff26b4565e27cdd, 0x3ff29e9df51fdee1,
    0x3ff2d285a6e4030b, 0x3ff306fe0a31b715, 0x3ff33c08b26416ff,
    0x3ff371a7373aa9cb, 0x3ff3a7db34e59ff7, 0x3ff3dea64c123422,
    0x3ff4160a21f72e2a, 0x3ff44e086061892d, 0x3ff486a2b5c13cd0,
    0x3ff4bfdad5362a27, 0x3ff4f9b2769d2ca7, 0x3ff5342b569d4f82,
    0x3ff56f4736b527da, 0x3ff5ab07dd485429, 0x3ff5e76f15ad2148,
    0x3ff6247eb03a5585, 0x3ff6623882552225, 0x3ff6a09e667f3bcd,
    0x3ff6dfb23c651a2f, 0x3ff71f75e8ec5f74, 0x3ff75feb564267c9,
    0x3ff7a11473eb0187, 0x3ff7e2f336cf4e62, 0x3ff82589994cce13,
    0x3ff868d99b4492ed, 0x3ff8ace5422aa0db, 0x3ff8f1ae99157736,
    0x3ff93737b0cdc5e5, 0x3ff97d829fde4e50, 0x3ff9c49182a3f090,
    0x3ffa0c667b5de565, 0x3ffa5503b23e255d, 0x3ffa9e6b5579fdbf,
    0x3ffae89f995ad3ad, 0x3ffb33a2b84f15fb, 0x3ffb7f76f2fb5e47,
    0x3ffbcc1e904bc1d2, 0x3ffc199bdd85529c, 0x3ffc67f12e57d14b,
    0x3ffcb720dcef9069, 0x3ffd072d4a07897c, 0x3ffd5818dcfba487,
    0x3ffda9e603db3285, 0x3ffdfc97337b9b5f, 0x3ffe502ee78b3ff6,
    0x3ffea4afa2a490da, 0x3ffefa1bee615a27, 0x3fff50765b6e4540,
    0x3fffa7c1819e90d8,
};

pub fn expf(value: f32) f32 {
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;

    if (magnitude == 0) return 1.0;
    if (magnitude > infinity_bits_f) {
        return @bitCast(bits | expf_quiet_nan_bit);
    }
    if (magnitude == infinity_bits_f) {
        if ((bits & sign_mask_f) == 0) return @bitCast(infinity_bits_f);
        return 0.0;
    }

    return expfFinite(value, bits);
}

fn expfFinite(value: f32, bits: u32) f32 {
    @setFloatMode(.strict);

    const inverse_ln2: f64 = 0x1.71547652b82fep+0;
    const rounding_constant: f64 = 0x1.8p46;
    const value_wide: f64 = value;
    const scaled = inverse_ln2 * value_wide;
    const rounded = scaled + rounding_constant;
    const rounded_bits: u64 = @bitCast(rounded);
    const absolute_bits = bits << 1;

    if (absolute_bits > 0x8562_e42e or absolute_bits < 0x6f93_813e) {
        if (absolute_bits < 0x6f93_813e) {
            return @floatCast(1.0 + value_wide * (1.0 + value_wide * 0.5));
        }
        if (bits > 0xc2ce_8ec0) {
            var result = 0x1p-149 +
                (value_wide + 0x1.9d1d9fccf477p+6) * 0x1.71547652b82edp-150;
            if (result < 0x1p-151) result = 0x1p-151;
            return @floatCast(result);
        }
        if ((bits & sign_mask_f) == 0 and bits > 0x42b1_7217) {
            return @bitCast(infinity_bits_f);
        }
    }

    const correction = rounding_constant - rounded;
    var reduced = scaled + correction;
    const table_index: u6 = @truncate(rounded_bits);
    const exponent_bits = (rounded_bits >> 6) << 52;
    const scale_bits = expf_scale_bits[table_index] +% exponent_bits;
    const scale: f64 = @bitCast(scale_bits);
    var reduced_squared = reduced * reduced;
    var result =
        ((expf_fast_coefficients[0] + reduced * expf_fast_coefficients[1]) +
            reduced_squared *
                (expf_fast_coefficients[2] + reduced * expf_fast_coefficients[3])) *
        scale;
    var upper: f32 = @floatCast(result);
    const lower: f32 = @floatCast(result - result * 1.45e-10);

    if (upper != lower) {
        const inverse_ln2_high: f64 = 0x1.7154765p+0;
        const inverse_ln2_low: f64 = 0x1.5c17f0bbbe88p-31;
        reduced =
            (inverse_ln2_high * value_wide + correction) +
            inverse_ln2_low * value_wide;
        reduced_squared = reduced * reduced;
        const scaled_remainder = scale * reduced;
        result = scale +
            scaled_remainder *
                ((expf_accurate_coefficients[0] +
                    reduced * expf_accurate_coefficients[1]) +
                    reduced_squared *
                        ((expf_accurate_coefficients[2] +
                            reduced * expf_accurate_coefficients[3]) +
                            reduced_squared *
                                (expf_accurate_coefficients[4] +
                                    reduced * expf_accurate_coefficients[5])));
        upper = @floatCast(result);
    }

    return upper;
}
// Correctly rounded binary64 exponential, ported from CORE-MATH revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
// Copyright (c) 2022-2025 Alexei Sibidanov.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

const ExpdDoubleDouble = struct {
    high: f64,
    low: f64,
};

const expd_database: [51]f64 = .{
    0x1.fffffffffffffp-53,  0x1.ba07d73250de7p-14,  0x1.6a4d1af9cc989p-8,   0x1.5a75293a5dcdap-6,
    0x1.42ea46949b3c7p-5,   0x1.7c8bb0cf5d16p-5,    0x1.0948d39a41695p-3,   0x1.a065fefae814fp-3,
    0x1.f6e4c3ced7c72p-3,   0x1.1a0408712e00ap-2,   0x1.bcab27d05abdep-2,   0x1.005ae04256babp-1,
    0x1.273c188aa7b14p+2,   0x1.83d4bcdebb3f4p+2,   0x1.08f51434652c3p+4,   0x1.1d5c2daebe367p+4,
    0x1.c44ce0d716a1ap+4,   0x1.e07e71bfcf06fp+5,   0x1.f7216c4b435c9p+5,   0x1.54cd1fea7663ap+7,
    0x1.d6479eba7c971p+8,   -0x1.664716b68a409p-14, -0x1.a2fefefd580dfp-13, -0x1.ce3f638d0c742p-12,
    -0x1.ceff32831e2c2p-12, -0x1.33accae78b371p-11, -0x1.d792b60084f92p-11, -0x1.7fb235d76cce7p-8,
    -0x1.1ff9b8e8b38bep-7,  -0x1.54511e930898cp-7,  -0x1.5c5ed0ec83666p-6,  -0x1.8c56ff5326197p-6,
    -0x1.a4187f2ca71f9p-6,  -0x1.a8f783d749a8fp-4,  -0x1.bd44fdaed819fp-4,  -0x1.daf693d64fadap-4,
    -0x1.290ea09e36479p-3,  -0x1.8aeb636f3ce35p-3,  -0x1.d3f3799439415p-3,  -0x1.ea16274b0109bp-3,
    -0x1.22e24fa3d5cf9p-1,  -0x1.85068c07fbbf6p-1,  -0x1.bdc7955d1482cp-1,  -0x1.2a9cad9998262p+0,
    -0x1.cc37ef7de7501p+0,  -0x1.02393d5976769p+1,  -0x1.65061daf79a78p+1,  -0x1.e8bdbfcd9144ep+3,
    -0x1.8f80e06f3a04cp+4,  -0x1.59f038076039cp+6,  -0x1.981587ad4542fp+7,
};

// Double-double approximations of 2^(i / 2^6), for 0 <= i < 64.
const expd_table_0: [64]ExpdDoubleDouble = .{ .{ .high = 0x1p+0, .low = 0x0p+0 }, .{ .high = 0x1.02c9a3e778061p+0, .low = -0x1.19083535b085ep-56 }, .{ .high = 0x1.059b0d3158574p+0, .low = 0x1.d73e2a475b466p-55 }, .{ .high = 0x1.0874518759bc8p+0, .low = 0x1.186be4bb285p-57 }, .{ .high = 0x1.0b5586cf9890fp+0, .low = 0x1.8a62e4adc610ap-54 }, .{ .high = 0x1.0e3ec32d3d1a2p+0, .low = 0x1.03a1727c57b52p-59 }, .{ .high = 0x1.11301d0125b51p+0, .low = -0x1.6c51039449b3ap-54 }, .{ .high = 0x1.1429aaea92dep+0, .low = -0x1.32fbf9af1369ep-54 }, .{ .high = 0x1.172b83c7d517bp+0, .low = -0x1.19041b9d78a76p-55 }, .{ .high = 0x1.1a35beb6fcb75p+0, .low = 0x1.e5b4c7b4968e4p-55 }, .{ .high = 0x1.1d4873168b9aap+0, .low = 0x1.e016e00a2643cp-54 }, .{ .high = 0x1.2063b88628cd6p+0, .low = 0x1.dc775814a8494p-55 }, .{ .high = 0x1.2387a6e756238p+0, .low = 0x1.9b07eb6c70572p-54 }, .{ .high = 0x1.26b4565e27cddp+0, .low = 0x1.2bd339940e9dap-55 }, .{ .high = 0x1.29e9df51fdee1p+0, .low = 0x1.612e8afad1256p-55 }, .{ .high = 0x1.2d285a6e4030bp+0, .low = 0x1.0024754db41d4p-54 }, .{ .high = 0x1.306fe0a31b715p+0, .low = 0x1.6f46ad23182e4p-55 }, .{ .high = 0x1.33c08b26416ffp+0, .low = 0x1.32721843659a6p-54 }, .{ .high = 0x1.371a7373aa9cbp+0, .low = -0x1.63aeabf42eae2p-54 }, .{ .high = 0x1.3a7db34e59ff7p+0, .low = -0x1.5e436d661f5e2p-56 }, .{ .high = 0x1.3dea64c123422p+0, .low = 0x1.ada0911f09ebcp-55 }, .{ .high = 0x1.4160a21f72e2ap+0, .low = -0x1.ef3691c309278p-58 }, .{ .high = 0x1.44e086061892dp+0, .low = 0x1.89b7a04ef80dp-59 }, .{ .high = 0x1.486a2b5c13cdp+0, .low = 0x1.3c1a3b69062fp-56 }, .{ .high = 0x1.4bfdad5362a27p+0, .low = 0x1.d4397afec42e2p-56 }, .{ .high = 0x1.4f9b2769d2ca7p+0, .low = -0x1.4b309d25957e4p-54 }, .{ .high = 0x1.5342b569d4f82p+0, .low = -0x1.07abe1db13cacp-55 }, .{ .high = 0x1.56f4736b527dap+0, .low = 0x1.9bb2c011d93acp-54 }, .{ .high = 0x1.5ab07dd485429p+0, .low = 0x1.6324c054647acp-54 }, .{ .high = 0x1.5e76f15ad2148p+0, .low = 0x1.ba6f93080e65ep-54 }, .{ .high = 0x1.6247eb03a5585p+0, .low = -0x1.383c17e40b496p-54 }, .{ .high = 0x1.6623882552225p+0, .low = -0x1.bb60987591c34p-54 }, .{ .high = 0x1.6a09e667f3bcdp+0, .low = -0x1.bdd3413b26456p-54 }, .{ .high = 0x1.6dfb23c651a2fp+0, .low = -0x1.bbe3a683c88aap-57 }, .{ .high = 0x1.71f75e8ec5f74p+0, .low = -0x1.16e4786887a9ap-55 }, .{ .high = 0x1.75feb564267c9p+0, .low = -0x1.0245957316dd4p-54 }, .{ .high = 0x1.7a11473eb0187p+0, .low = -0x1.41577ee04993p-55 }, .{ .high = 0x1.7e2f336cf4e62p+0, .low = 0x1.05d02ba15797ep-56 }, .{ .high = 0x1.82589994cce13p+0, .low = -0x1.d4c1dd41532d8p-54 }, .{ .high = 0x1.868d99b4492edp+0, .low = -0x1.fc6f89bd4f6bap-54 }, .{ .high = 0x1.8ace5422aa0dbp+0, .low = 0x1.6e9f156864b26p-54 }, .{ .high = 0x1.8f1ae99157736p+0, .low = 0x1.5cc13a2e3976cp-55 }, .{ .high = 0x1.93737b0cdc5e5p+0, .low = -0x1.75fc781b57ebcp-57 }, .{ .high = 0x1.97d829fde4e5p+0, .low = -0x1.d185b7c1b85dp-54 }, .{ .high = 0x1.9c49182a3f09p+0, .low = 0x1.c7c46b071f2bep-56 }, .{ .high = 0x1.a0c667b5de565p+0, .low = -0x1.359495d1cd532p-54 }, .{ .high = 0x1.a5503b23e255dp+0, .low = -0x1.d2f6edb8d41e2p-54 }, .{ .high = 0x1.a9e6b5579fdbfp+0, .low = 0x1.0fac90ef7fd32p-54 }, .{ .high = 0x1.ae89f995ad3adp+0, .low = 0x1.7a1cd345dcc82p-54 }, .{ .high = 0x1.b33a2b84f15fbp+0, .low = -0x1.2805e3084d708p-57 }, .{ .high = 0x1.b7f76f2fb5e47p+0, .low = -0x1.5584f7e54ac3ap-56 }, .{ .high = 0x1.bcc1e904bc1d2p+0, .low = 0x1.23dd07a2d9e84p-55 }, .{ .high = 0x1.c199bdd85529cp+0, .low = 0x1.11065895048dep-55 }, .{ .high = 0x1.c67f12e57d14bp+0, .low = 0x1.2884dff483cacp-54 }, .{ .high = 0x1.cb720dcef9069p+0, .low = 0x1.503cbd1e949dcp-56 }, .{ .high = 0x1.d072d4a07897cp+0, .low = -0x1.cbc3743797a9cp-54 }, .{ .high = 0x1.d5818dcfba487p+0, .low = 0x1.2ed02d75b3706p-55 }, .{ .high = 0x1.da9e603db3285p+0, .low = 0x1.c2300696db532p-54 }, .{ .high = 0x1.dfc97337b9b5fp+0, .low = -0x1.1a5cd4f184b5cp-54 }, .{ .high = 0x1.e502ee78b3ff6p+0, .low = 0x1.39e8980a9cc9p-55 }, .{ .high = 0x1.ea4afa2a490dap+0, .low = -0x1.e9c23179c2894p-54 }, .{ .high = 0x1.efa1bee615a27p+0, .low = 0x1.dc7f486a4b6bp-54 }, .{ .high = 0x1.f50765b6e454p+0, .low = 0x1.9d3e12dd8a18ap-54 }, .{ .high = 0x1.fa7c1819e90d8p+0, .low = 0x1.74853f3a5931ep-55 } };

// Double-double approximations of 2^(i / 2^12), for 0 <= i < 64.
const expd_table_1: [64]ExpdDoubleDouble = .{
    .{ .high = 0x1p+0, .low = 0x0p+0 },                               .{ .high = 0x1.000b175effdc7p+0, .low = 0x1.ae8e38c59c72ap-54 },
    .{ .high = 0x1.00162f3904052p+0, .low = -0x1.7b5d0d58ea8f4p-58 }, .{ .high = 0x1.0021478e11ce6p+0, .low = 0x1.4115cb6b16a8ep-54 },
    .{ .high = 0x1.002c605e2e8cfp+0, .low = -0x1.d7c96f201bb2ep-55 }, .{ .high = 0x1.003779a95f959p+0, .low = 0x1.84711d4c35eap-54 },
    .{ .high = 0x1.0042936faa3d8p+0, .low = -0x1.0484245243778p-55 }, .{ .high = 0x1.004dadb113dap+0, .low = -0x1.4b237da2025fap-54 },
    .{ .high = 0x1.0058c86da1c0ap+0, .low = -0x1.5e00e62d6b30ep-56 }, .{ .high = 0x1.0063e3a559473p+0, .low = 0x1.a1d6cedbb948p-54 },
    .{ .high = 0x1.006eff583fc3dp+0, .low = -0x1.4acf197a00142p-54 }, .{ .high = 0x1.007a1b865a8cap+0, .low = -0x1.eaf2ea42391a6p-57 },
    .{ .high = 0x1.0085382faef83p+0, .low = 0x1.da93f90835f76p-56 },  .{ .high = 0x1.00905554425d4p+0, .low = -0x1.6a79084ab093cp-55 },
    .{ .high = 0x1.009b72f41a12bp+0, .low = 0x1.86364f8fbe8f8p-54 },  .{ .high = 0x1.00a6910f3b6fdp+0, .low = -0x1.82e8e14e3110ep-55 },
    .{ .high = 0x1.00b1afa5abcbfp+0, .low = -0x1.4f6b2a7609f72p-55 }, .{ .high = 0x1.00bcceb7707ecp+0, .low = -0x1.e1a258ea8f71ap-56 },
    .{ .high = 0x1.00c7ee448ee02p+0, .low = 0x1.4362ca5bc26f2p-56 },  .{ .high = 0x1.00d30e4d0c483p+0, .low = 0x1.095a56c919d02p-54 },
    .{ .high = 0x1.00de2ed0ee0f5p+0, .low = -0x1.406ac4e81a646p-57 }, .{ .high = 0x1.00e94fd0398ep+0, .low = 0x1.b5a6902767e08p-54 },
    .{ .high = 0x1.00f4714af41d3p+0, .low = -0x1.91b206085932p-54 },  .{ .high = 0x1.00ff93412315cp+0, .low = 0x1.427068ab22306p-55 },
    .{ .high = 0x1.010ab5b2cbd11p+0, .low = 0x1.c1d0660524e08p-54 },  .{ .high = 0x1.0115d89ff3a8bp+0, .low = -0x1.e7bdfb3204be8p-54 },
    .{ .high = 0x1.0120fc089ff63p+0, .low = 0x1.843aa8b9cbbc6p-55 },  .{ .high = 0x1.012c1fecd613bp+0, .low = -0x1.34104ee7edae8p-56 },
    .{ .high = 0x1.0137444c9b5b5p+0, .low = -0x1.2b6aeb6176892p-56 }, .{ .high = 0x1.01426927f5278p+0, .low = 0x1.a8cd33b8a1bb2p-56 },
    .{ .high = 0x1.014d8e7ee8d2fp+0, .low = 0x1.2edc08e5da99ap-56 },  .{ .high = 0x1.0158b4517bb88p+0, .low = 0x1.57ba2dc7e0c72p-55 },
    .{ .high = 0x1.0163da9fb3335p+0, .low = 0x1.b61299ab8cdb8p-54 },  .{ .high = 0x1.016f0169949edp+0, .low = -0x1.90565902c5f44p-54 },
    .{ .high = 0x1.017a28af25567p+0, .low = 0x1.70fc41c5c2d54p-55 },  .{ .high = 0x1.018550706ab62p+0, .low = 0x1.4b9a6e145d76cp-54 },
    .{ .high = 0x1.019078ad6a19fp+0, .low = -0x1.008eff5142bfap-56 }, .{ .high = 0x1.019ba16628de2p+0, .low = -0x1.77669f033c7dep-54 },
    .{ .high = 0x1.01a6ca9aac5f3p+0, .low = -0x1.09bb78eeead0ap-54 }, .{ .high = 0x1.01b1f44af9f9ep+0, .low = 0x1.371231477ece6p-54 },
    .{ .high = 0x1.01bd1e77170b4p+0, .low = 0x1.5e7626621eb5ap-56 },  .{ .high = 0x1.01c8491f08f08p+0, .low = -0x1.bc72b100828a4p-54 },
    .{ .high = 0x1.01d37442d507p+0, .low = -0x1.ce39cbbab8bbep-57 },  .{ .high = 0x1.01de9fe280ac8p+0, .low = 0x1.16996709da2e2p-55 },
    .{ .high = 0x1.01e9cbfe113efp+0, .low = -0x1.c11f5239bf536p-55 }, .{ .high = 0x1.01f4f8958c1c6p+0, .low = 0x1.e1d4eb5edc6b4p-55 },
    .{ .high = 0x1.020025a8f6a35p+0, .low = -0x1.afb99946ee3fp-54 },  .{ .high = 0x1.020b533856324p+0, .low = -0x1.8f06d8a148a32p-54 },
    .{ .high = 0x1.02168143b0281p+0, .low = -0x1.2bf310fc54eb6p-55 }, .{ .high = 0x1.0221afcb09e3ep+0, .low = -0x1.c95a035eb4176p-54 },
    .{ .high = 0x1.022cdece68c4fp+0, .low = -0x1.491793e46834cp-54 }, .{ .high = 0x1.02380e4dd22adp+0, .low = -0x1.3e8d0d9c4909p-56 },
    .{ .high = 0x1.02433e494b755p+0, .low = -0x1.314aa16278aa4p-54 }, .{ .high = 0x1.024e6ec0da046p+0, .low = 0x1.48daf888e965p-55 },
    .{ .high = 0x1.02599fb483385p+0, .low = 0x1.56dc8046821f4p-55 },  .{ .high = 0x1.0264d1244c719p+0, .low = 0x1.45b42356b9d46p-54 },
    .{ .high = 0x1.027003103b10ep+0, .low = -0x1.082ef51b61d7ep-56 }, .{ .high = 0x1.027b357854772p+0, .low = 0x1.2106ed0920a34p-56 },
    .{ .high = 0x1.0286685c9e059p+0, .low = -0x1.fd4cf26ea5d0ep-54 }, .{ .high = 0x1.02919bbd1d1d8p+0, .low = -0x1.09f8775e78084p-54 },
    .{ .high = 0x1.029ccf99d720ap+0, .low = 0x1.64cbba902ca28p-58 },  .{ .high = 0x1.02a803f2d170dp+0, .low = 0x1.4383ef231d206p-54 },
    .{ .high = 0x1.02b338c811703p+0, .low = 0x1.4a47a505b3a46p-54 },  .{ .high = 0x1.02be6e199c811p+0, .low = 0x1.e47120223468p-54 },
};

const expd_accurate_coefficients: [7]ExpdDoubleDouble = .{
    .{ .high = 0x1p+0, .low = 0x0p+0 },
    .{ .high = 0x1p-1, .low = 0x1.712f72ecec2cfp-99 },
    .{ .high = 0x1.5555555555555p-3, .low = 0x1.5555555554d07p-57 },
    .{ .high = 0x1.5555555555555p-5, .low = 0x1.55194d28275dap-59 },
    .{ .high = 0x1.1111111111111p-7, .low = 0x1.12faa0e1c0f7bp-63 },
    .{ .high = 0x1.6c16c16da6973p-10, .low = -0x1.4ba45ab25d2a3p-64 },
    .{ .high = 0x1.a01a019eb7f31p-13, .low = -0x1.9091d845ecd36p-67 },
};

fn expdRoundEvenFinite(value: f64) i64 {
    const truncated: i64 = @intFromFloat(value);
    const truncated_float: f64 = @floatFromInt(truncated);
    const fraction = value - truncated_float;
    const odd = (truncated & 1) != 0;

    if (fraction > 0x1p-1 or (fraction == 0x1p-1 and odd)) return truncated + 1;
    if (fraction < -0x1p-1 or (fraction == -0x1p-1 and odd)) return truncated - 1;
    return truncated;
}

fn expdFastTwoSum(x: f64, y: f64) ExpdDoubleDouble {
    const sum = x + y;
    const rounded_x = sum - x;
    return .{ .high = sum, .low = y - rounded_x };
}

fn expdFastSum(x_high: f64, x_low: f64, y_high: f64, y_low: f64) ExpdDoubleDouble {
    const sum = expdFastTwoSum(x_high, y_high);
    return .{ .high = sum.high, .low = (x_low + y_low) + sum.low };
}

fn expdMultiplyDoubleDouble(
    x_high: f64,
    x_low: f64,
    c_high: f64,
    c_low: f64,
) ExpdDoubleDouble {
    const high = c_high * x_high;
    const low = (c_high * x_low + c_low * x_high) + @mulAdd(f64, c_high, x_high, -high);
    return .{ .high = high, .low = low };
}

fn expdFastFmaDoubleWord(
    a_high: f64,
    a_low: f64,
    b_high: f64,
    b_low: f64,
    c_high: f64,
    c_low: f64,
) ExpdDoubleDouble {
    const high = @mulAdd(f64, a_high, b_high, c_high);
    const correction = c_high - high;
    const product_error = @mulAdd(f64, a_high, b_high, correction);
    const combined_error = product_error + c_low;
    const cross = a_high * b_low + combined_error;
    const low = a_low * b_high + cross;
    return .{ .high = high, .low = low };
}

fn expdEvaluateDoubleDoublePolynomial(
    x_high: f64,
    x_low: f64,
    coefficients: []const ExpdDoubleDouble,
) ExpdDoubleDouble {
    var index = coefficients.len - 1;
    var result = coefficients[index];
    while (index != 0) {
        index -= 1;
        const coefficient = coefficients[index];
        result = expdFastFmaDoubleWord(
            x_high,
            x_low,
            result.high,
            result.low,
            coefficient.high,
            coefficient.low,
        );
    }
    return result;
}

fn expdScaleNormal(value: f64, exponent: i64) f64 {
    const bits: u64 = @bitCast(value);
    const exponent_bits = @as(u64, @bitCast(exponent)) << 52;
    return @bitCast(bits +% exponent_bits);
}

fn expdToDenormal(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    return @bitCast(bits & 0x000f_ffff_ffff_ffff);
}

fn expdDatabaseLookup(value: f64, fallback: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    var lower: i32 = 0;
    var upper: i32 = @intCast(expd_database.len - 1);

    while (lower <= upper) {
        const middle = @divTrunc(lower + upper, 2);
        const middle_index: usize = @intCast(middle);
        const database_bits: u64 = @bitCast(expd_database[middle_index]);
        if (database_bits < value_bits) {
            lower = middle + 1;
        } else if (database_bits == value_bits) {
            const rounding_patterns = [2]u64{ 0x57f5fe2e5bde4075, 0x0000003c1f16b8ed };
            const sign_pattern: u64 = 333811522313371;
            const middle_shift: u6 = @intCast(middle);
            const delta_bits = (((sign_pattern >> middle_shift) & 1) << 63) | 0x3c90_0000_0000_0000;
            const delta: f64 = @bitCast(delta_bits);
            const pattern_index: usize = @intCast(middle >> 5);
            const pattern_shift: u6 = @intCast((middle << 1) & 63);
            const rounding_pattern = (rounding_patterns[pattern_index] >> pattern_shift) & 3;
            const fallback_bits: u64 = @bitCast(fallback);

            var offset: i64 = -1;
            while (offset <= 1) : (offset += 1) {
                const candidate_bits = if (offset < 0)
                    fallback_bits -% 1
                else if (offset > 0)
                    fallback_bits +% 1
                else
                    fallback_bits;
                if ((candidate_bits & 3) == rounding_pattern) {
                    const candidate: f64 = @bitCast(candidate_bits);
                    return candidate + delta;
                }
            }
            break;
        } else {
            upper = middle - 1;
        }
    }
    return fallback;
}

fn expdAccurate(value: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    if (((value_bits >> 52) & 0x7ff) < 0x3c9) return 1.0 + value;

    const reduction_scale = 0x1.71547652b82fep+12;
    const reduction_integer = expdRoundEvenFinite(value * reduction_scale);
    const reduction: f64 = @floatFromInt(reduction_integer);
    const table_0_index: usize = @intCast((reduction_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(reduction_integer & 0x3f);
    const exponent = reduction_integer >> 12;
    const table_0 = expd_table_0[table_0_index];
    const table_1 = expd_table_1[table_1_index];
    const table_product = expdMultiplyDoubleDouble(
        table_0.high,
        table_0.low,
        table_1.high,
        table_1.low,
    );

    const log2_high = 0x1.62e42ffp-13;
    const log2_low = 0x1.718432a1b0e26p-47;
    const log2_tail = 0x1.9ff0342542fc3p-102;
    const reduced = value - log2_high * reduction;
    var reduced_low = log2_low * reduction;
    const reduced_tail = log2_tail * reduction + @mulAdd(f64, log2_low, reduction, -reduced_low);
    const reduced_high = reduced + reduced_low;
    reduced_low = (reduced - reduced_high) + reduced_low + reduced_tail;

    var polynomial = expdEvaluateDoubleDoublePolynomial(
        reduced_high,
        reduced_low,
        &expd_accurate_coefficients,
    );
    polynomial = expdMultiplyDoubleDouble(
        reduced_high,
        reduced_low,
        polynomial.high,
        polynomial.low,
    );

    if (value_bits > 0xc086_232b_dd7a_bcd2) {
        const scale_bits = @as(u64, @intCast(1 - exponent)) << 52;
        const scale: f64 = @bitCast(scale_bits);
        const product = expdMultiplyDoubleDouble(
            polynomial.high,
            polynomial.low,
            table_product.high,
            table_product.low,
        );
        var result = expdFastSum(
            table_product.high,
            table_product.low,
            product.high,
            product.low,
        );
        const scaled_sum = expdFastTwoSum(scale, result.high);
        result.high = scaled_sum.high;
        result.low += scaled_sum.low;
        return expdToDenormal(result.high + result.low);
    }

    var result = ExpdDoubleDouble{ .high = 0.0, .low = 0.0 };
    if (table_product.high == 1.0) {
        const high_sum = expdFastTwoSum(table_product.high, polynomial.high);
        const low_sum = expdFastTwoSum(high_sum.low, polynomial.low);
        result = .{ .high = high_sum.high, .low = low_sum.high };

        var low_bits: u64 = @bitCast(result.low);
        if ((low_bits & 0x000f_ffff_ffff_ffff) == 0) {
            const error_bits: u64 = @bitCast(low_sum.low);
            if (((low_bits ^ error_bits) & sign_mask_d) == 0) {
                low_bits +%= 1;
            } else {
                low_bits -%= 1;
            }
            result.low = @bitCast(low_bits);
        }
    } else {
        const product = expdMultiplyDoubleDouble(
            polynomial.high,
            polynomial.low,
            table_product.high,
            table_product.low,
        );
        result = expdFastSum(
            table_product.high,
            table_product.low,
            product.high,
            product.low,
        );
    }

    result = expdFastTwoSum(result.high, result.low);
    const low_bits: u64 = @bitCast(result.low);
    const ambiguity = (low_bits +% 2) & 0x000f_ffff_ffff_ffff;
    if (ambiguity <= 2) result.high = expdDatabaseLookup(value, result.high);
    return expdScaleNormal(result.high, exponent);
}

pub fn expd(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask_d;

    if (magnitude > infinity_bits_d) {
        return @bitCast(bits | 0x0008_0000_0000_0000);
    }
    if (magnitude == infinity_bits_d) {
        if ((bits & sign_mask_d) != 0) return 0.0;
        return value;
    }
    if (magnitude == 0) return 1.0;

    // exp(x) rounds to one in nearest-even mode for |x| <= 2^-54.
    if (magnitude <= 0x3c90_0000_0000_0000) return 1.0 + value;

    if (magnitude >= 0x4086_2e42_fefa_39f0) {
        if ((bits & sign_mask_d) == 0) return @bitCast(infinity_bits_d);
        if (magnitude >= 0x4087_4910_d52d_3052) return 0.0;
    }

    const reduction_scale = 0x1.71547652b82fep+12;
    const reduction_integer = expdRoundEvenFinite(value * reduction_scale);
    const reduction: f64 = @floatFromInt(reduction_integer);
    const table_0_index: usize = @intCast((reduction_integer >> 6) & 0x3f);
    const table_1_index: usize = @intCast(reduction_integer & 0x3f);
    const exponent = reduction_integer >> 12;
    const table_0 = expd_table_0[table_0_index];
    const table_1 = expd_table_1[table_1_index];
    const table_product = expdMultiplyDoubleDouble(
        table_0.high,
        table_0.low,
        table_1.high,
        table_1.low,
    );

    const log2_high = 0x1.62e42ffp-13;
    const log2_low = 0x1.718432a1b0e26p-47;
    const reduced = (value - log2_high * reduction) + log2_low * reduction;
    const reduced_squared = reduced * reduced;
    const coefficient_0 = 0x1p+0;
    const coefficient_1 = 0x1p-1;
    const coefficient_2 = 0x1.55555557e54ffp-3;
    const coefficient_3 = 0x1.55555553a12f4p-5;
    const polynomial = (coefficient_0 + reduced * coefficient_1) +
        reduced_squared * (coefficient_2 + reduced * coefficient_3);

    var high = table_product.high;
    const table_times_reduced = table_product.high * reduced;
    var low = table_product.low + table_times_reduced * polynomial;
    const error_bound = 1.64e-19;

    if (bits > 0xc086_232b_dd7a_bcd2) {
        const scale_bits = @as(u64, @intCast(1 - exponent)) << 52;
        const scale: f64 = @bitCast(scale_bits);
        const sum = expdFastTwoSum(scale, high);
        high = sum.high;
        low += sum.low;
        const upper = high + (low + error_bound);
        const lower = high + (low - error_bound);
        if (upper != lower) return expdAccurate(value);
        return expdToDenormal(lower);
    }

    const upper = high + (low + error_bound);
    const lower = high + (low - error_bound);
    if (upper != lower) return expdAccurate(value);
    return expdScaleNormal(lower, exponent);
}
// Correctly rounded sine for binary32 values.
//
// Copyright (c) 2022-2025 Alexei Sibidanov.
//
// This code is derived from CORE-MATH src/binary32/sin/sinf.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

const sinf_reduction_pi_inverse: [4]u64 = .{
    0xfe51_63ab_debb_c562,
    0xdb62_9599_3c43_9041,
    0xfc27_57d1_f534_ddc0,
    0xa2f9_836e_4e44_1529,
};

const sinf_b: [4]f64 = .{
    0x1.3bd3cc9be45dcp-6,
    -0x1.03c1f081b0833p-14,
    0x1.55d3c6fc9ac1fp-24,
    -0x1.e1d3ff281b40dp-35,
};

const sinf_a: [4]f64 = .{
    0x1.921fb54442d17p-3,
    -0x1.4abbce6256a39p-10,
    0x1.466bc5a518c16p-19,
    -0x1.32bdc61074ff6p-29,
};

const sinf_table: [32]f64 = .{
    0x0p+0,
    0x1.8f8b83c69a60bp-3,
    0x1.87de2a6aea963p-2,
    0x1.1c73b39ae68c8p-1,
    0x1.6a09e667f3bcdp-1,
    0x1.a9b66290ea1a3p-1,
    0x1.d906bcf328d46p-1,
    0x1.f6297cff75cbp-1,
    0x1p+0,
    0x1.f6297cff75cbp-1,
    0x1.d906bcf328d46p-1,
    0x1.a9b66290ea1a3p-1,
    0x1.6a09e667f3bcdp-1,
    0x1.1c73b39ae68c8p-1,
    0x1.87de2a6aea963p-2,
    0x1.8f8b83c69a60bp-3,
    0x0p+0,
    -0x1.8f8b83c69a60bp-3,
    -0x1.87de2a6aea963p-2,
    -0x1.1c73b39ae68c8p-1,
    -0x1.6a09e667f3bcdp-1,
    -0x1.a9b66290ea1a3p-1,
    -0x1.d906bcf328d46p-1,
    -0x1.f6297cff75cbp-1,
    -0x1p+0,
    -0x1.f6297cff75cbp-1,
    -0x1.d906bcf328d46p-1,
    -0x1.a9b66290ea1a3p-1,
    -0x1.6a09e667f3bcdp-1,
    -0x1.1c73b39ae68c8p-1,
    -0x1.87de2a6aea963p-2,
    -0x1.8f8b83c69a60bp-3,
};

const sinf_hard_cases: [4]struct { arg: u32, rh: f32, rl: f32 } = .{
    .{ .arg = 0x4619_9998, .rh = -0x1.63f4bap-2, .rl = -0x1p-27 },
    .{ .arg = 0x3f3a_dc51, .rh = 0x1.55688ap-1, .rl = -0x1p-26 },
    .{ .arg = 0x3fa7_832a, .rh = 0x1.ee836cp-1, .rl = -0x1p-26 },
    .{ .arg = 0x4116_cbe4, .rh = -0x1.99bc5ap-26, .rl = -0x1p-51 },
};

fn sinf_round_even(value: f64) f64 {
    @setFloatMode(.strict);
    const truncated: i64 = @intFromFloat(value);
    const truncated_float: f64 = @floatFromInt(truncated);
    const fraction = value - truncated_float;

    if (fraction > 0.5) return truncated_float + 1.0;
    if (fraction < -0.5) return truncated_float - 1.0;
    if (fraction == 0.5 and (truncated & 1) != 0) return truncated_float + 1.0;
    if (fraction == -0.5 and (truncated & 1) != 0) return truncated_float - 1.0;
    return truncated_float;
}

fn sinf_rltl(value: f32, quadrant: *i32) f64 {
    @setFloatMode(.strict);
    const value_wide: f64 = value;
    const id_low = -0x1.b1bbead603d8bp-29 * value_wide;
    const id_high = 0x1.45f306ep+2 * value_wide;
    const id = sinf_round_even(id_high);
    const quadrant_bits: u64 = @bitCast(0x1.8p52 + id);
    quadrant.* = @bitCast(@as(u32, @truncate(quadrant_bits)));
    return (id_high - id) + id_low;
}

fn sinf_rltl0(value: f64, quadrant: *i32) f64 {
    @setFloatMode(.strict);
    const id_high = 0x1.45f306dc9c883p+2 * value;
    const id = sinf_round_even(id_high);
    const quadrant_bits: u64 = @bitCast(0x1.8p52 + id);
    quadrant.* = @bitCast(@as(u32, @truncate(quadrant_bits)));
    return id_high - id;
}

fn sinf_rbig(bits: u32, quadrant: *i32) f64 {
    @setFloatMode(.strict);
    const exponent: i32 = @intCast((bits >> 23) & 0xff);
    const mantissa: u64 = @as(u64, bits & 0x007f_ffff) | 0x0080_0000;
    const product_0: u128 = @as(u128, mantissa) * sinf_reduction_pi_inverse[0];
    var product_1: u128 = @as(u128, mantissa) * sinf_reduction_pi_inverse[1];
    product_1 += product_0 >> 64;
    var product_2: u128 = @as(u128, mantissa) * sinf_reduction_pi_inverse[2];
    product_2 += product_1 >> 64;
    var product_3: u128 = @as(u128, mantissa) * sinf_reduction_pi_inverse[3];
    product_3 += product_2 >> 64;

    const product_3_high: u64 = @truncate(product_3 >> 64);
    const product_3_low: u64 = @truncate(product_3);
    const product_2_low: u64 = @truncate(product_2);
    const product_1_low: u64 = @truncate(product_1);
    const shift: i32 = exponent - 147;
    var integer_bits: u64 = undefined;
    var fraction_bits: u64 = undefined;

    if (shift < 64) {
        const left_shift: u6 = @intCast(shift);
        const right_shift: u6 = @intCast(64 - shift);
        integer_bits = (product_3_high << left_shift) | (product_3_low >> right_shift);
        fraction_bits = (product_3_low << left_shift) | (product_2_low >> right_shift);
    } else if (shift == 64) {
        integer_bits = product_3_low;
        fraction_bits = product_2_low;
    } else {
        const left_shift: u6 = @intCast(shift - 64);
        const right_shift: u6 = @intCast(128 - shift);
        integer_bits = (product_3_low << left_shift) | (product_2_low >> right_shift);
        fraction_bits = (product_2_low << left_shift) | (product_1_low >> right_shift);
    }

    var integer_part: i32 = @bitCast(@as(u32, @truncate(integer_bits)));
    const fraction_part: i64 = @bitCast(fraction_bits);
    const sign: i64 = if ((bits >> 31) != 0) -1 else 0;
    const fraction_sign: i64 = fraction_part >> 63;
    integer_part -= @intCast(fraction_sign);
    const reduced = @as(f64, @floatFromInt(fraction_part ^ sign)) * 0x1p-64;
    const integer_sign: i32 = @intCast(sign);
    integer_part = (integer_part ^ integer_sign) - integer_sign;
    quadrant.* = integer_part;
    return reduced;
}

fn sinf_add_sign(value: f32, high: f32, low: f32) f32 {
    @setFloatMode(.strict);
    const sign: f32 = if ((@as(u32, @bitCast(value)) & sign_mask_f) != 0) -1.0 else 1.0;
    return sign * high + sign * low;
}

fn sinf_database(value: f32, reduced: f64) f32 {
    const magnitude = @as(u32, @bitCast(value)) & magnitude_mask_f;
    for (sinf_hard_cases) |entry| {
        if (entry.arg == magnitude) return sinf_add_sign(value, entry.rh, entry.rl);
    }
    return @floatCast(reduced);
}

fn sinf_big(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const absolute_bits = bits << 1;
    if (absolute_bits >= 0xff00_0000) {
        if ((absolute_bits << 8) != 0) return @bitCast(bits | 0x0040_0000);
        return @bitCast(@as(u32, 0x7fc0_0000));
    }

    var quadrant: i32 = 0;
    const reduced = sinf_rbig(bits, &quadrant);
    const reduced_squared = reduced * reduced;
    const reduced_fourth = reduced_squared * reduced_squared;
    const polynomial_a = (sinf_a[0] + reduced_squared * sinf_a[1]) +
        reduced_fourth * (sinf_a[2] + reduced_squared * sinf_a[3]);
    const polynomial_b = (sinf_b[0] + reduced_squared * sinf_b[1]) +
        reduced_fourth * (sinf_b[2] + reduced_squared * sinf_b[3]);
    const sine_table = sinf_table[@intCast(@as(u32, @bitCast(quadrant)) & 31)];
    const cosine_table =
        sinf_table[@intCast((@as(u32, @bitCast(quadrant)) +% 8) & 31)];
    const result = sine_table + reduced * (polynomial_a * cosine_table -
        polynomial_b * (reduced * sine_table));
    return @floatCast(result);
}

pub fn sinf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const absolute_bits = bits << 1;
    var quadrant: i32 = 0;
    var reduced: f64 = 0.0;

    if (absolute_bits > 0x9900_0000 or absolute_bits < 0x7300_0000) {
        if (absolute_bits < 0x7300_0000) {
            if (absolute_bits < 0x6600_0000) {
                if (absolute_bits == 0) return value;
                const absolute_value: f32 = @bitCast(bits & magnitude_mask_f);
                const result = @mulAdd(f32, -value, absolute_value, value);
                return result;
            }
            return (-0x1.555556p-3 * value) * (value * value) + value;
        }
        return sinf_big(value);
    }

    if (absolute_bits < 0x822d_97c8) {
        if (absolute_bits == 0x7e75_b8a2 or absolute_bits == 0x7f4f_0654) {
            return sinf_database(value, 0.0);
        }
        reduced = sinf_rltl0(@as(f64, value), &quadrant);
    } else {
        if (absolute_bits == 0x8c33_3330) return sinf_database(value, 0.0);
        reduced = sinf_rltl(value, &quadrant);
    }

    const reduced_squared = reduced * reduced;
    const reduced_fourth = reduced_squared * reduced_squared;
    const polynomial_a = (sinf_a[0] + reduced_squared * sinf_a[1]) +
        reduced_fourth * (sinf_a[2] + reduced_squared * sinf_a[3]);
    const polynomial_b = (sinf_b[0] + reduced_squared * sinf_b[1]) +
        reduced_fourth * (sinf_b[2] + reduced_squared * sinf_b[3]);
    const sine_table = sinf_table[@intCast(@as(u32, @bitCast(quadrant)) & 31)];
    const cosine_table =
        sinf_table[@intCast((@as(u32, @bitCast(quadrant)) +% 8) & 31)];
    const result = sine_table + polynomial_a * (reduced * cosine_table) -
        polynomial_b * (reduced_squared * sine_table);
    return @floatCast(result);
}
// Correctly rounded cosine for binary32 values.
//
// Copyright (c) 2022-2023 Alexei Sibidanov.
//
// This code is derived from CORE-MATH src/binary32/cos/cosf.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and is used under the MIT License:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

const cosf_reduction_pi_inverse: [4]u64 = .{
    0xfe51_63ab_debb_c562,
    0xdb62_9599_3c43_9041,
    0xfc27_57d1_f534_ddc0,
    0xa2f9_836e_4e44_1529,
};

const cosf_hard_cases: [5]struct { arg: u32, rh: f32, rl: f32 } = .{
    .{ .arg = 0x4096_cbe4, .rh = 0x1.99bc5cp-27, .rl = -0x1p-52 },
    .{ .arg = 0x5922_aa80, .rh = 0x1.115d7ep-1, .rl = -0x1p-26 },
    .{ .arg = 0x5aa4_542c, .rh = 0x1.f48148p-2, .rl = 0x1p-27 },
    .{ .arg = 0x5f18_b878, .rh = 0x1.fe2976p-1, .rl = 0x1p-26 },
    .{ .arg = 0x6115_cb11, .rh = 0x1.f0285ep-1, .rl = -0x1p-26 },
};

const cosf_table: [32]f64 = .{
    0x1p+0,
    0x1.f6297cff75cbp-1,
    0x1.d906bcf328d46p-1,
    0x1.a9b66290ea1a3p-1,
    0x1.6a09e667f3bcdp-1,
    0x1.1c73b39ae68c8p-1,
    0x1.87de2a6aea963p-2,
    0x1.8f8b83c69a60bp-3,
    0x0p+0,
    -0x1.8f8b83c69a60bp-3,
    -0x1.87de2a6aea963p-2,
    -0x1.1c73b39ae68c8p-1,
    -0x1.6a09e667f3bcdp-1,
    -0x1.a9b66290ea1a3p-1,
    -0x1.d906bcf328d46p-1,
    -0x1.f6297cff75cbp-1,
    -0x1p+0,
    -0x1.f6297cff75cbp-1,
    -0x1.d906bcf328d46p-1,
    -0x1.a9b66290ea1a3p-1,
    -0x1.6a09e667f3bcdp-1,
    -0x1.1c73b39ae68c8p-1,
    -0x1.87de2a6aea963p-2,
    -0x1.8f8b83c69a60bp-3,
    0x0p+0,
    0x1.8f8b83c69a60bp-3,
    0x1.87de2a6aea963p-2,
    0x1.1c73b39ae68c8p-1,
    0x1.6a09e667f3bcdp-1,
    0x1.a9b66290ea1a3p-1,
    0x1.d906bcf328d46p-1,
    0x1.f6297cff75cbp-1,
};

const cosf_b: [4]f64 = .{
    0x1.3bd3cc9be45dcp-6,
    -0x1.03c1f081b0833p-14,
    0x1.55d3c6fc9ac1fp-24,
    -0x1.e1d3ff281b40dp-35,
};

const cosf_a: [4]f64 = .{
    0x1.921fb54442d17p-3,
    -0x1.4abbce6256a39p-10,
    0x1.466bc5a518c16p-19,
    -0x1.32bdc61074ff6p-29,
};

fn cosfRltl(value: f32, quadrant: *i32) f64 {
    @setFloatMode(.strict);
    const value_wide: f64 = value;
    const id_low = -0x1.b1bbead603d8bp-29 * value_wide;
    const id_high = 0x1.45f306ep+2 * value_wide;
    const id = sinf_round_even(id_high);
    const quadrant_bits: u64 = @bitCast(0x1.8p52 + id);
    quadrant.* = @bitCast(@as(u32, @truncate(quadrant_bits)));
    return (id_high - id) + id_low;
}

fn cosfRltl0(value: f64, quadrant: *i32) f64 {
    @setFloatMode(.strict);
    const id_high = 0x1.45f306dc9c883p+2 * value;
    const id = sinf_round_even(id_high);
    const quadrant_bits: u64 = @bitCast(0x1.8p52 + id);
    quadrant.* = @bitCast(@as(u32, @truncate(quadrant_bits)));
    return id_high - id;
}

fn cosfRbig(bits: u32, quadrant: *i32) f64 {
    @setFloatMode(.strict);
    const exponent: i32 = @intCast((bits >> 23) & 0xff);
    const mantissa: u64 = @as(u64, bits & 0x007f_ffff) | 0x0080_0000;
    const product_0: u128 = @as(u128, mantissa) * cosf_reduction_pi_inverse[0];
    var product_1: u128 = @as(u128, mantissa) * cosf_reduction_pi_inverse[1];
    product_1 += product_0 >> 64;
    var product_2: u128 = @as(u128, mantissa) * cosf_reduction_pi_inverse[2];
    product_2 += product_1 >> 64;
    var product_3: u128 = @as(u128, mantissa) * cosf_reduction_pi_inverse[3];
    product_3 += product_2 >> 64;

    const product_3_high: u64 = @truncate(product_3 >> 64);
    const product_3_low: u64 = @truncate(product_3);
    const product_2_low: u64 = @truncate(product_2);
    const product_1_low: u64 = @truncate(product_1);
    const shift: i32 = exponent - 147;
    // SAFETY: Every shift branch assigns both values before they are read.
    var integer_bits: u64 = undefined;
    // SAFETY: The same shift branch assigns this value before it is read.
    var fraction_bits: u64 = undefined;

    if (shift < 64) {
        const left_shift: u6 = @intCast(shift);
        const right_shift: u6 = @intCast(64 - shift);
        integer_bits = (product_3_high << left_shift) | (product_3_low >> right_shift);
        fraction_bits = (product_3_low << left_shift) | (product_2_low >> right_shift);
    } else if (shift == 64) {
        integer_bits = product_3_low;
        fraction_bits = product_2_low;
    } else {
        const left_shift: u6 = @intCast(shift - 64);
        const right_shift: u6 = @intCast(128 - shift);
        integer_bits = (product_3_low << left_shift) | (product_2_low >> right_shift);
        fraction_bits = (product_2_low << left_shift) | (product_1_low >> right_shift);
    }

    var integer_part: i32 = @bitCast(@as(u32, @truncate(integer_bits)));
    const fraction_part: i64 = @bitCast(fraction_bits);
    const sign: i64 = if ((bits & sign_mask_f) != 0) -1 else 0;
    const fraction_sign: i64 = fraction_part >> 63;
    integer_part -= @intCast(fraction_sign);
    const reduced = @as(f64, @floatFromInt(fraction_part ^ sign)) * 0x1p-64;
    const integer_sign: i32 = @intCast(sign);
    integer_part = (integer_part ^ integer_sign) - integer_sign;
    quadrant.* = integer_part;
    return reduced;
}

fn cosfDatabase(value: f32, reduced: f64) f32 {
    @setFloatMode(.strict);
    const magnitude = @as(u32, @bitCast(value)) & magnitude_mask_f;
    for (cosf_hard_cases) |entry| {
        if (entry.arg == magnitude) return entry.rh + entry.rl;
    }
    return @floatCast(reduced);
}

fn cosfBig(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const absolute_bits = bits << 1;
    if (absolute_bits >= 0xff00_0000) {
        if ((absolute_bits << 8) != 0) return @bitCast(bits | 0x0040_0000);
        return @bitCast(@as(u32, 0x7fc0_0000));
    }

    var quadrant: i32 = 0;
    const reduced = cosfRbig(bits, &quadrant);
    const reduced_squared = reduced * reduced;
    const reduced_fourth = reduced_squared * reduced_squared;
    const polynomial_a = (cosf_a[0] + reduced_squared * cosf_a[1]) +
        reduced_fourth * (cosf_a[2] + reduced_squared * cosf_a[3]);
    const polynomial_b = (cosf_b[0] + reduced_squared * cosf_b[1]) +
        reduced_fourth * (cosf_b[2] + reduced_squared * cosf_b[3]);
    const cosine_table = cosf_table[@intCast(@as(u32, @bitCast(quadrant)) & 31)];
    const sine_table =
        cosf_table[@intCast((@as(u32, @bitCast(quadrant)) +% 8) & 31)];
    const result = cosine_table + reduced * (polynomial_a * sine_table -
        polynomial_b * (reduced * cosine_table));
    const result_bits: u64 = @bitCast(result);
    const tail = (result_bits + 6) & (~@as(u64, 0) >> 36);
    if (tail <= 12) return cosfDatabase(value, result);
    return @floatCast(result);
}

/// Correctly rounded binary32 cosine. Callers must keep nearest-even rounding
/// active; this function does not inspect or modify the floating-point mode.
pub fn cosf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const absolute_bits = bits << 1;
    if (absolute_bits >= 0xff00_0000) {
        if ((absolute_bits << 8) != 0) return @bitCast(bits | 0x0040_0000);
        return @bitCast(@as(u32, 0x7fc0_0000));
    }

    if (absolute_bits > 0x9900_0000 or absolute_bits < 0x7400_0000) {
        if (absolute_bits < 0x7400_0000) {
            return 0.5 * @mulAdd(f32, value, -value, 2.0);
        }
        return cosfBig(value);
    }

    var quadrant: i32 = 0;
    const reduced: f64 = if (absolute_bits < 0x82a4_1896) blk: {
        if (absolute_bits == 0x812d_97c8) return cosfDatabase(value, 0.0);
        break :blk cosfRltl0(@as(f64, value), &quadrant);
    } else cosfRltl(value, &quadrant);

    const reduced_squared = reduced * reduced;
    const reduced_fourth = reduced_squared * reduced_squared;
    const polynomial_a = (cosf_a[0] + reduced_squared * cosf_a[1]) +
        reduced_fourth * (cosf_a[2] + reduced_squared * cosf_a[3]);
    const polynomial_b = (cosf_b[0] + reduced_squared * cosf_b[1]) +
        reduced_fourth * (cosf_b[2] + reduced_squared * cosf_b[3]);
    const cosine_table = cosf_table[@intCast(@as(u32, @bitCast(quadrant)) & 31)];
    const sine_table =
        cosf_table[@intCast((@as(u32, @bitCast(quadrant)) +% 8) & 31)];
    const result = cosine_table + polynomial_a * (reduced * sine_table) -
        polynomial_b * (reduced_squared * cosine_table);
    return @floatCast(result);
}
