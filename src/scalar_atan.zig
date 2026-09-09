// Correctly rounded arctangent for binary16, binary32, and binary64.
//
// The finite kernels are ported from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa:
//   src/binary16/atan/atanf16.c (Copyright (c) 2025 Paul Zimmermann).
//   src/binary32/atan/atanf.c (Copyright (c) 2022-2025 Alexei Sibidanov).
//   src/binary64/atan/atan.c (Copyright (c) 2023 Alexei Sibidanov).
//
// This code is derived from CORE-MATH and is used under the MIT License:
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

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;

fn copySignH(magnitude: f16, bits: u16) f16 {
    const magnitude_bits: u16 = @bitCast(magnitude);
    return @bitCast((magnitude_bits & ~sign_mask_h) | (bits & sign_mask_h));
}

fn copySignF(magnitude: f32, bits: u32) f32 {
    const magnitude_bits: u32 = @bitCast(magnitude);
    return @bitCast((magnitude_bits & ~sign_mask_f) | (bits & sign_mask_f));
}

fn copySignD(magnitude: f64, bits: u64) f64 {
    const magnitude_bits: u64 = @bitCast(magnitude);
    return @bitCast((magnitude_bits & ~sign_mask_d) | (bits & sign_mask_d));
}

/// Return correctly rounded binary16 atan with caller-selected nearest-even mode.
pub fn atanh(value: f16) f16 {
    @setFloatMode(.strict);

    const input_bits: u16 = @bitCast(value);
    const magnitude_bits = input_bits & magnitude_mask_h;
    if (magnitude_bits > infinity_bits_h) return @bitCast(input_bits | quiet_nan_bit_h);
    if (magnitude_bits == infinity_bits_h) return copySignH(@bitCast(@as(u16, 0x3e48)), input_bits);

    const t: f32 = @floatCast(value);
    const t_bits: u32 = @bitCast(t);
    const absolute_bits = t_bits & magnitude_mask_f;
    const negative = (t_bits & sign_mask_f) != 0;
    const sign: f32 = if (negative) -1.0 else 1.0;
    var reduced: f32 = @bitCast(absolute_bits);
    const do_reduce = absolute_bits > 0x3f80_0000;
    if (do_reduce) reduced = 1.0 / reduced;

    const p0 = [_]f32{ 0x1.fffffcp-1, -0x1.55546cp-2, 0x1.98d0ep-3, -0x1.0c7c54p-3 };
    const p1 = [_]f32{ 0x1.411612p-14, 0x1.ff076cp-1, 0x1.1ee64cp-6, -0x1.a6fc96p-2, 0x1.81dd3ep-3 };
    const p2 = [_]f32{ -0x1.95964cp-8, 0x1.0bad76p+0, -0x1.e652e2p-4, -0x1.e70ab4p-3, 0x1.a5eb88p-4 };
    const p3 = [_]f32{ -0x1.f75b16p-6, 0x1.2d3d1p+0, -0x1.882376p-2, 0x1.d18c96p-13, 0x1.6aa268p-6 };

    var y: f32 = undefined;
    const square = reduced * reduced;
    if (reduced <= 0.25) {
        if (!do_reduce and reduced <= 0x1.d14p-6) {
            if (absolute_bits == 0) return value;
            const signed = sign * reduced;
            return @floatCast(@mulAdd(f32, signed, -0x1p-23, signed));
        }
        if (absolute_bits == 0x3e56_a000) return @floatCast(copySignF(0x1.a72002p-3, t_bits));
        if (absolute_bits == 0x4115_c000) return @floatCast(copySignF(0x1.76dffep+0, t_bits));
        if (absolute_bits == 0x42c3_2000) return @floatCast(copySignF(0x1.8f7ffep+0, t_bits));

        const c5 = @mulAdd(f32, p0[3], square, p0[2]);
        const c1 = @mulAdd(f32, p0[1], square, p0[0]);
        const c1_adjusted = @mulAdd(f32, c5, square * square, c1);
        y = reduced * c1_adjusted;
    } else {
        const p = if (reduced <= 0.5) p1 else if (reduced <= 0.75) p2 else p3;
        const c3 = @mulAdd(f32, p[4], reduced, p[3]);
        const c2 = @mulAdd(f32, c3, reduced, p[2]);
        y = @mulAdd(f32, p[1], reduced, p[0]);
        y = @mulAdd(f32, c2, square, y);
    }

    if (do_reduce) y = 0x1.921fb6p+0 - y;
    if (negative) y = -y;
    return @floatCast(y);
}

/// Return correctly rounded binary32 atan with caller-selected nearest-even mode.
pub fn atanf(value: f32) f32 {
    @setFloatMode(.strict);

    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;
    if (magnitude > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    if (magnitude == infinity_bits_f) return copySignF(@floatCast(0x1.921fb54442d18p+0), bits);

    const exponent = (magnitude >> 23) & 0xff;
    const greater_equal_one = exponent >= 127;
    if (magnitude >= 0x4c70_0518) return copySignF(@floatCast(0x1.921fb54442d18p+0), bits);
    if (exponent < 127 - 13) {
        if (exponent < 127 - 25) {
            if (magnitude == 0) return value;
            const absolute_value: f32 = @bitCast(magnitude);
            return @mulAdd(f32, -value, absolute_value, value);
        }
        const linear_coefficient = -0x1.5555555555555p-2 * value;
        return @mulAdd(f32, linear_coefficient, value * value, value);
    }

    var z: f64 = @floatCast(value);
    if (greater_equal_one) z = 1.0 / z;
    const z2 = z * z;
    const z4 = z2 * z2;
    const z8 = z4 * z4;
    const cn = [_]f64{
        0x1.51eccde075d67p-2, 0x1.a76bb5637f2f2p-1, 0x1.81e0eed20de88p-1,
        0x1.376c8ca67d11dp-2, 0x1.aec7b69202ac6p-5, 0x1.9561899acc73ep-9,
        0x1.bf9fa5b67e6p-16,
    };
    const cd = [_]f64{
        0x1.51eccde075d66p-2,  0x1.dfbdd7b392d28p-1, 0x1p+0,
        0x1.fd22bf0e89b54p-2,  0x1.d91ff8b576282p-4, 0x1.653ea99fc9bbp-7,
        0x1.1e7fcc202340ap-12,
    };
    var cn0 = cn[0] + z2 * cn[1];
    const cn2 = cn[2] + z2 * cn[3];
    var cn4 = cn[4] + z2 * cn[5];
    const cn6 = cn[6];
    cn0 += z4 * cn2;
    cn4 += z4 * cn6;
    cn0 += z8 * cn4;
    cn0 *= z;
    var cd0 = cd[0] + z2 * cd[1];
    const cd2 = cd[2] + z2 * cd[3];
    var cd4 = cd[4] + z2 * cd[5];
    const cd6 = cd[6];
    cd0 += z4 * cd2;
    cd4 += z4 * cd6;
    cd0 += z8 * cd4;
    const result = cn0 / cd0;
    if (!greater_equal_one) return @floatCast(result);

    const pi_over_two_high: f64 = 0x1.9p0;
    const pi_over_two_low: f64 = 0x1.0fdaa22168c23p-7;
    const sign_d: f64 = if ((bits & sign_mask_f) != 0) -1.0 else 1.0;
    const reconstructed = sign_d * pi_over_two_low - result + sign_d * pi_over_two_high;
    return @floatCast(reconstructed);
}

fn fastTwoSum(x: f64, y: f64, low: *f64) f64 {
    const sum = x + y;
    const z = sum - x;
    low.* = y - z;
    return sum;
}

fn addDoubleDouble(x_high: f64, x_low: f64, c_high: f64, c_low: f64, low: *f64) f64 {
    const sum = x_high + c_high;
    const delta = sum - x_high;
    low.* = ((c_high - delta) + (x_high + (delta - sum))) + (x_low + c_low);
    return sum;
}

fn mulDoubleDoubleAcc(x_high: f64, x_low: f64, c_high: f64, c_low: f64, low: *f64) f64 {
    const high_low = c_high * x_low;
    const low_high = c_low * x_high;
    const high_high = c_high * x_high;
    const high_error = @mulAdd(f64, c_high, x_high, -high_high);
    const middle = high_error + (low_high + high_low);
    return fastTwoSum(high_high, middle, low);
}

const atan_table = [_][2]f64{ .{ 0x0p+0, 0x0p+0 }, .{ 0x1.9224e047e368ep-7, 0x1.a3ca6c727c59dp-62 }, .{ 0x1.92346247a91fp-6, 0x1.138b0ef96a186p-64 }, .{ 0x1.2dbaae9a05dbp-5, 0x1.36e7f8a3f5e42p-59 }, .{ 0x1.927278a3b1162p-5, -0x1.ac986efb92662p-64 }, .{ 0x1.f7495ea3f3783p-5, 0x1.06ec8011ee816p-59 }, .{ 0x1.2e239ccff3831p-4, -0x1.858437d431332p-58 }, .{ 0x1.60b9f7597fdecp-4, -0x1.cebd13eb7c513p-60 }, .{ 0x1.936bb8c5b2da2p-4, -0x1.840cac0d81db5p-58 }, .{ 0x1.c63ce377fc802p-4, 0x1.400b0fdaa109ep-58 }, .{ 0x1.f93183a8db9e9p-4, 0x1.0e04e06c86e72p-59 }, .{ 0x1.1626d85a91e7p-3, 0x1.f7ad829163ca7p-59 }, .{ 0x1.2fcac73a6064p-3, -0x1.2680735ce2cd8p-58 }, .{ 0x1.4986a74cf4e57p-3, -0x1.90559690b42e4p-57 }, .{ 0x1.635c990ce0d36p-3, 0x1.91d29110b41aap-58 }, .{ 0x1.7d4ec54fb5968p-3, -0x1.ea90e2718278p-59 }, .{ 0x1.975f5e0553158p-3, -0x1.dc82ac14e3e1cp-61 }, .{ 0x1.b1909efd8b762p-3, -0x1.73a10fd13daafp-58 }, .{ 0x1.cbe4ceb4b4cf2p-3, -0x1.3a7ffbeabda0bp-57 }, .{ 0x1.e65e3f27c9f2ap-3, -0x1.db6627a24d523p-57 }, .{ 0x1.007fa758626aep-2, -0x1.45f97dd3099f6p-57 }, .{ 0x1.0de53475f3b3cp-2, -0x1.6293f68741816p-57 }, .{ 0x1.1b6103d3597e9p-2, -0x1.ab240d40633e9p-57 }, .{ 0x1.28f459ecad74dp-2, -0x1.de34d14e832ep-61 }, .{ 0x1.36a08355c63dcp-2, 0x1.af540d9fb4926p-57 }, .{ 0x1.4466d542bac92p-2, 0x1.da60fdbc82ac4p-57 }, .{ 0x1.5248ae1701b17p-2, -0x1.92a601170138ap-56 }, .{ 0x1.604775fbb27dfp-2, -0x1.7f1fca1d5d15bp-57 }, .{ 0x1.6e649f7d78649p-2, -0x1.4e223ea716c7bp-57 }, .{ 0x1.7ca1a832d0f84p-2, 0x1.b24c824ac51fcp-56 }, .{ 0x1.8b00196b3d022p-2, 0x1.4314cd132ba43p-57 }, .{ 0x1.998188e816bfp-2, -0x1.11f1e0817879ap-56 }, .{ 0x1.a827999fcef32p-2, -0x1.c3dea4dbad538p-57 }, .{ 0x1.b6f3fc8c61e5bp-2, 0x1.60d1b780ee3ebp-57 }, .{ 0x1.c5e87185e67b6p-2, -0x1.ab5edb7dfa545p-59 }, .{ 0x1.d506c82a2c8p-2, -0x1.8e1437048b5bdp-57 }, .{ 0x1.e450e0d273e7ap-2, -0x1.06951c97b050fp-56 }, .{ 0x1.f3c8ad985d9eep-2, -0x1.14af9522ab518p-59 }, .{ 0x1.01b819b5a7cf7p-1, -0x1.aba0d7d97d1f2p-56 }, .{ 0x1.09a4c59bd0d4dp-1, 0x1.095bc4ebc2c42p-59 }, .{ 0x1.11ab7190834ecp-1, 0x1.798826fa27774p-55 }, .{ 0x1.19cd3fe8e405dp-1, 0x1.008f6258fc98fp-55 }, .{ 0x1.220b5ef047825p-1, -0x1.462af7ceb7de6p-58 }, .{ 0x1.2a6709a74f289p-1, -0x1.1184dfd78b472p-56 }, .{ 0x1.32e1889047ffdp-1, 0x1.9141876dc40c5p-56 }, .{ 0x1.3b7c3289ed6f3p-1, 0x1.481c20189726cp-55 }, .{ 0x1.44386db9ce5dbp-1, 0x1.2e851bd025441p-55 }, .{ 0x1.4d17b087b265dp-1, 0x1.13ada9b8bc419p-56 }, .{ 0x1.561b82ab7f99p-1, -0x1.05b4c3c4cbee8p-55 }, .{ 0x1.5f457e4f4812ep-1, -0x1.5619249bd96f1p-55 }, .{ 0x1.6897514751db6p-1, -0x1.b0a0fbcafc671p-57 }, .{ 0x1.7212be621be6dp-1, -0x1.19ff2dc66da45p-55 }, .{ 0x1.7bb99ed2990cfp-1, 0x1.1320449592d92p-55 }, .{ 0x1.858de3b716571p-1, -0x1.1fddcd2f3da8ep-55 }, .{ 0x1.8f9197bf85eebp-1, 0x1.d44a42e35cc97p-57 }, .{ 0x1.99c6e0f634394p-1, -0x1.585a178b4a18dp-56 }, .{ 0x1.a43002ae4285p-1, 0x1.f95a531b3a97p-57 }, .{ 0x1.aecf5f9ba35a6p-1, -0x1.96c2d43ca3392p-60 }, .{ 0x1.b9a77c18c1af2p-1, -0x1.a5bed94b05defp-57 }, .{ 0x1.c4bb009e77983p-1, 0x1.54509d2bff511p-59 }, .{ 0x1.d00cbc7384d2ep-1, -0x1.b4c867cef300cp-57 }, .{ 0x1.db9fa89953fcfp-1, -0x1.ddfac663d6bc6p-62 }, .{ 0x1.e776eafc91706p-1, -0x1.a510683ff7cb6p-56 }, .{ 0x1.f395d9f0e3c92p-1, 0x1.4fdcd8e4e871p-59 }, .{ 0x1p+0, 0x0p+0 }, .{ 0x1.065c900aaf2d8p+0, -0x1.deec7fc9042adp-55 }, .{ 0x1.0ce29d0883c99p+0, -0x1.395ae45e0657dp-55 }, .{ 0x1.139447e6a86eep+0, 0x1.332cf301a97f3p-55 }, .{ 0x1.1a73d55278c4bp+0, -0x1.6cc8c4b78213bp-55 }, .{ 0x1.2183b0c4573ffp+0, 0x1.70a90841da57ap-55 }, .{ 0x1.28c66fdaf8f09p+0, -0x1.ba39bad450eep-57 }, .{ 0x1.303ed61109e2p+0, -0x1.8692946d9f93cp-55 }, .{ 0x1.37efd8d87607ep+0, 0x1.3b711bf765b58p-57 }, .{ 0x1.3fdca42847507p+0, 0x1.c21387985b081p-56 }, .{ 0x1.48089f8bf42ccp+0, -0x1.7ddb19d3d0efcp-55 }, .{ 0x1.507773c537eadp+0, -0x1.f5e354cf971f3p-56 }, .{ 0x1.592d11142fa55p+0, -0x1.00f0ad675330dp-56 }, .{ 0x1.622db63c8ecc2p+0, -0x1.2c93f50ab2c0ep-55 }, .{ 0x1.6b7df862652p+0, 0x1.bec391adc37d5p-56 }, .{ 0x1.7522cbdd428a8p+0, -0x1.9686ddc9ffcf5p-57 }, .{ 0x1.7f218e25a7461p+0, -0x1.8d16529514246p-56 }, .{ 0x1.89801106cc709p+0, -0x1.092f51e9c2803p-55 }, .{ 0x1.9444a7462122ap+0, -0x1.07c06755404c4p-55 }, .{ 0x1.9f7632fa9e871p+0, 0x1.02e0d43abc92bp-55 }, .{ 0x1.ab1c35d8a74eap+0, 0x1.d0184e48af6f7p-58 }, .{ 0x1.b73ee3c3ef16ap+0, 0x1.73be957380bc2p-56 }, .{ 0x1.c3e738086bc0fp+0, -0x1.02b6e26c84462p-56 }, .{ 0x1.d11f0dae40609p+0, 0x1.25c4f3ffa6e1fp-58 }, .{ 0x1.def13b73c1406p+0, -0x1.e302db3c6823fp-58 }, .{ 0x1.ed69b4153a45dp+0, 0x1.3207830326c0ep-56 }, .{ 0x1.fc95abad6cf4ap+0, -0x1.6308cee7927bfp-57 }, .{ 0x1.0641e192ceab3p+1, -0x1.0147ebf0df4c5p-56 }, .{ 0x1.0ea21d716fbf7p+1, -0x1.168533cc41d8bp-56 }, .{ 0x1.17749711a6679p+1, -0x1.52a0b0333e9c5p-57 }, .{ 0x1.20c36c6a7f38ep+1, 0x1.8659eece35395p-57 }, .{ 0x1.2a99f50fd4f4fp+1, 0x1.20fcad18cb36fp-55 }, .{ 0x1.3504f333f9de6p+1, -0x1.52afdbd5a8c74p-56 }, .{ 0x1.4012ce2586a17p+1, -0x1.9747a792907d7p-56 }, .{ 0x1.4bd3d87fe065p+1, 0x1.90c59393b52c8p-56 }, .{ 0x1.585aa4e1530fap+1, 0x1.af6934f13a3a8p-56 }, .{ 0x1.65bc6cc825147p+1, -0x1.8534dcab5ad3ep-59 }, .{ 0x1.74118e4b6a7c8p+1, -0x1.555aa8bfca9a1p-56 }, .{ 0x1.837626d70fdb8p+1, -0x1.56b3fee9ca72bp-58 }, .{ 0x1.940ad30abc792p+1, 0x1.4b3fdd4fdc06cp-58 }, .{ 0x1.a5f59e90600ddp+1, 0x1.285d367c55ddcp-57 }, .{ 0x1.b9633283b6d14p+1, -0x1.8712976f17a16p-59 }, .{ 0x1.ce885653127e7p+1, -0x1.abe8ab65d49fcp-60 }, .{ 0x1.e5a3de972a377p+1, 0x1.cd9be81ad764bp-58 }, .{ 0x1.ff01305ecd8dcp+1, 0x1.742c2922656fap-59 }, .{ 0x1.0d7dc7cff4c9ep+2, -0x1.7c842978bee09p-56 }, .{ 0x1.1d0143e71565fp+2, 0x1.7bc7dea7c3c03p-57 }, .{ 0x1.2e4ff1626b949p+2, 0x1.aefbe25b404e9p-59 }, .{ 0x1.41bfee2424771p+2, -0x1.4bcfaaa95cb2cp-60 }, .{ 0x1.57be4eaa5e11bp+2, 0x1.0fe741e4ec679p-58 }, .{ 0x1.70d751908c1b1p+2, 0x1.fe74a5b0ec709p-58 }, .{ 0x1.8dc25c117782bp+2, 0x1.0ca1c19f710efp-58 }, .{ 0x1.af73f4ca3310fp+2, 0x1.2867b40ba77d6p-58 }, .{ 0x1.d7398d15e70dbp+2, 0x1.0fd4e0d4b1547p-57 }, .{ 0x1.0372fb36b87e2p+3, 0x1.c16c9ecc1621dp-58 }, .{ 0x1.208dbdae055efp+3, 0x1.6b81a36e75e8cp-58 }, .{ 0x1.44e6c595afdccp+3, -0x1.7c22045771848p-58 }, .{ 0x1.7398c57f3f1adp+3, 0x1.970503be105cp-58 }, .{ 0x1.b1d03c03d2f7fp+3, -0x1.f299d010aead2p-60 }, .{ 0x1.046e9fe60a77ep+4, 0x1.d2b61deff33ecp-58 }, .{ 0x1.45affed201b55p+4, 0x1.0e84d9567203ap-64 }, .{ 0x1.b267195b1ffaep+4, -0x1.ad44b44b92653p-64 }, .{ 0x1.45e2455e4aaa7p+5, -0x1.296d577b5e21dp-60 }, .{ 0x1.45eed6854ce99p+6, 0x1.2db53886013cap-63 } };

const atan_index_coefficients = [_][3]u64{ .{ 419, 81, 0 }, .{ 500, 81, 0 }, .{ 582, 163, 0 }, .{ 745, 163, 0 }, .{ 908, 326, 0 }, .{ 1234, 326, 0 }, .{ 1559, 651, 0 }, .{ 2210, 650, 1 }, .{ 2860, 1299, 3 }, .{ 4156, 1293, 4 }, .{ 5444, 2569, 24 }, .{ 7989, 2520, 32 }, .{ 10476, 4917, 168 }, .{ 15224, 4576, 200 }, .{ 19601, 8341, 838 }, .{ 27105, 6648, 731 }, .{ 33036, 10210, 1998 }, .{ 41266, 6292, 1117 }, .{ 46469, 7926, 2048 }, .{ 52375, 4038, 849 }, .{ 55587, 4591, 1291 }, .{ 58906, 2172, 479 }, .{ 60612, 2390, 688 }, .{ 62325, 1107, 247 }, .{ 63192, 1207, 349 }, .{ 64056, 556, 124 }, .{ 64491, 605, 175 }, .{ 64923, 278, 62 }, .{ 65141, 303, 88 }, .{ 65358, 139, 31 }, .{ 65467, 151, 44 } };

fn atanPolyDoubleDouble(x_high: f64, x_low: f64, low: *f64) f64 {
    const coefficients = [_][2]f64{
        .{ -0x1.5555555555555p-2, -0x1.5555555555555p-56 },
        .{ 0x1.999999999999ap-3, -0x1.999999999bcb8p-57 },
        .{ -0x1.2492492492492p-3, -0x1.249242093c016p-57 },
    };
    var high = coefficients[2][0] + low.*;
    var coefficient_low = (coefficients[2][0] - high) + low.* + coefficients[2][1];
    var index: i32 = 2;
    while (index > 0) {
        index -= 1;
        high = mulDoubleDoubleAcc(x_high, x_low, high, coefficient_low, &coefficient_low);
        const temporary_high = high + coefficients[@intCast(index)][0];
        const temporary_low = (coefficients[@intCast(index)][0] - temporary_high) + high;
        high = temporary_high;
        coefficient_low += temporary_low + coefficients[@intCast(index)][1];
    }
    low.* = coefficient_low;
    return high;
}

fn atanRefine(value: f64, approximation: f64) f64 {
    const low_coefficients = [_]f64{
        0x1.c71c71c71c71cp-4, -0x1.745d1745d1265p-4,
        0x1.3b13b115bcbc4p-4, -0x1.1107c41ad3253p-4,
    };
    const value_bits: u64 = @bitCast(value);
    const value_sign: f64 = if ((value_bits & sign_mask_d) == 0) 1.0 else -1.0;
    const phi_bits: u64 = @bitCast(@abs(approximation) * 0x1.45f306dc9c883p6 + 256.5);
    const table_index: usize = @intCast((phi_bits >> 44) & 0xff);

    var h: f64 = undefined;
    var h_low: f64 = undefined;
    if (table_index == 128) {
        h = -1.0 / value;
        h_low = @mulAdd(f64, h, value, 1.0) * h;
    } else {
        const table_value = copySignD(atan_table[table_index][0], value_bits);
        const product = value * table_value;
        const product_low = @mulAdd(f64, value, table_value, -product);
        const difference = value - table_value;
        const sum = 1.0 + product;
        const delta = 1.0 - sum;
        const sum_error = (delta + product) - ((delta + sum) - 1.0) + product_low;
        const reciprocal = 1.0 / sum;
        const reciprocal_low = (@mulAdd(f64, reciprocal, -sum, 1.0) - sum_error * reciprocal) * reciprocal;
        h = reciprocal * difference;
        h_low = @mulAdd(f64, reciprocal, difference, -h) + reciprocal_low * difference;
    }

    var h_squared_low: f64 = undefined;
    const h_squared = mulDoubleDoubleAcc(h, h_low, h, h_low, &h_squared_low);
    const h_fourth = h_squared * h_squared;
    var h_cubed_low: f64 = undefined;
    const h_cubed = mulDoubleDoubleAcc(h, h_low, h_squared, h_squared_low, &h_cubed_low);
    var polynomial_low = h_squared *
        ((low_coefficients[0] + h_squared * low_coefficients[1]) +
            h_fourth * (low_coefficients[2] + h_squared * low_coefficients[3]));
    var polynomial = atanPolyDoubleDouble(h_squared, h_squared_low, &polynomial_low);
    polynomial = mulDoubleDoubleAcc(h_cubed, h_cubed_low, polynomial, polynomial_low, &polynomial_low);

    var high: f64 = undefined;
    var low: f64 = undefined;
    var tail: f64 = undefined;
    if (table_index == 0) {
        high = h;
        low = polynomial;
        tail = polynomial_low;
    } else {
        const delta_f = if (table_index < 128) value_sign * atan_table[table_index][1] else 0.0;
        const signed_index = value_sign * @as(f64, @floatFromInt(table_index));
        high = 0x1.921fb54442dp-7 * signed_index;
        low = 0x1.8469898cc518p-55 * signed_index;
        tail = -0x1.fc8f8cbb5bf8p-104 * signed_index;
        low = addDoubleDouble(low, tail, delta_f, 0.0, &tail);
        low = addDoubleDouble(low, tail, h, h_low, &tail);
        low = addDoubleDouble(low, tail, polynomial, polynomial_low, &tail);
    }

    var second_low: f64 = undefined;
    const first = fastTwoSum(high, low, &second_low);
    const second = fastTwoSum(second_low, tail, &second_low);
    const first_bits: u64 = @bitCast(first);
    var second_bits: u64 = @bitCast(second);
    const first_exponent = (first_bits >> 52) & 0x7ff;
    const second_exponent = (second_bits >> 52) & 0x7ff;
    if (((second_bits +% 1) & (~@as(u64, 0) >> 12)) <= 2 or
        (first_exponent > second_exponent and first_exponent - second_exponent > 103))
    {
        const difficult = [_][3]f64{
            .{ 0x1.0dc89a3b5501p-7, 0x1.0dc70ac228717p-7, 0x1p-61 },
            .{ 0x1.e3fb41d2d226p-8, 0x1.e3f9013a852f8p-8, 0x1p-62 },
            .{ 0x1.7ba49f739829fp-1, 0x1.46ac372243536p-1, 0x1p-109 },
            .{ 0x1.a933fe176b375p-3, 0x1.a33f32ac5ceb5p-3, -0x1p-112 },
            .{ 0x1.bb04a79820063p-8, 0x1.bb02ed5c5e956p-8, -0x1p-115 },
            .{ 0x1.cd30a9499618bp-8, 0x1.cd2eb65f92a46p-8, -0x1p-112 },
            .{ 0x1.f44aa37b8e66bp-7, 0x1.f440b04187c87p-7, -0x1p-112 },
            .{ 0x1.fd2ac95e57ef9p-8, 0x1.fd2829febc03ap-8, -0x1p-112 },
            .{ 0x1.6419079bbf601p-6, 0x1.640aade8f5427p-6, 0x1p-114 },
            .{ 0x1.7ba49f739829fp-1, 0x1.46ac372243536p-1, 0x1p-110 },
            .{ 0x1.d768804487b07p-3, 0x1.cf5676f373ec1p-3, -0x1p-110 },
            .{ 0x1.bb04a79820063p-8, 0x1.bb02ed5c5e956p-8, -0x1p-115 },
        };
        const absolute_value = @abs(value);
        for (difficult) |entry| {
            if (absolute_value == entry[0]) {
                return copySignD(entry[1], value_bits) + value_sign * entry[2];
            }
        }
        if ((second_bits & (~@as(u64, 0) >> 12)) == 0) {
            const third_bits: u64 = @bitCast(second_low);
            if (((third_bits ^ second_bits) >> 63) != 0) {
                second_bits -%= 1;
            } else {
                second_bits +%= 1;
            }
            second_bits = @bitCast(@as(f64, @bitCast(second_bits)));
            return @as(f64, @bitCast(second_bits)) + first;
        }
    }
    return second + first;
}

/// Return correctly rounded binary64 atan with caller-selected nearest-even mode.
pub fn atand(value: f64) f64 {
    @setFloatMode(.strict);

    const bits: u64 = @bitCast(value);
    const absolute_bits = bits & magnitude_mask_d;
    if (absolute_bits > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
    if (absolute_bits == infinity_bits_d) return copySignD(0x1.921fb54442d18p+0, bits);

    var index: i64 = @intCast(absolute_bits >> 51);
    index -= 2030;
    if (absolute_bits < 0x3f7b_21c4_75e6_362a) {
        if (absolute_bits == 0) return value;
        if (absolute_bits < 0x3e40_0000_0000_0000) return @mulAdd(f64, -0x1p-54, value, value);
        const x_squared = value * value;
        const x_cubed = value * x_squared;
        const x_fourth = x_squared * x_squared;
        const coefficients = [_]f64{
            -0x1.5555555555555p-2, 0x1.99999999998c1p-3,
            -0x1.249249176aecp-3,  0x1.c711fd121ae8p-4,
        };
        const polynomial = x_cubed *
            ((coefficients[0] + x_squared * coefficients[1]) +
                x_fourth * (coefficients[2] + x_squared * coefficients[3]));
        const upper = (polynomial + polynomial * 0x4.8p-52) + value;
        const lower = (polynomial - polynomial * 0x2.8p-52) + value;
        if (upper == lower) return upper;
        return atanRefine(value, upper);
    }

    var h: f64 = undefined;
    var high: f64 = undefined;
    var low: f64 = undefined;
    if (absolute_bits > 0x4062_ded8_e34a_9035) {
        high = copySignD(0x1.921fb54442d18p+0, bits);
        low = copySignD(0x1.1a62633145c07p-54, bits);
        if (absolute_bits >= 0x434d_0296_7c31_cdb5) return high + low;
        h = -1.0 / value;
    } else {
        const masked = bits & (~@as(u64, 0) >> 13);
        const upper = masked >> 35;
        const upper_squared = upper * upper >> 16;
        const coeff = atan_index_coefficients[@intCast(index)];
        index = @intCast(((@as(u64, coeff[0]) << 16) + upper * coeff[1] - upper_squared * coeff[2]) >> 25);
        const table_value = copySignD(atan_table[@intCast(index)][0], bits);
        const signed_index = copySignD(@floatFromInt(index), bits);
        const value_sign: f64 = if ((bits & sign_mask_d) == 0) 1.0 else -1.0;
        low = value_sign * atan_table[@intCast(index)][1] +
            0x1.8469898cc517p-55 * signed_index;
        h = (value - table_value) / (1.0 + value * table_value);
        high = 0x1.921fb54442dp-7 * signed_index;
    }
    const h_squared = h * h;
    const h_fourth = h_squared * h_squared;
    const coefficients = [_]f64{ 0x1p+0, -0x1.555555555552bp-2, 0x1.9999999069c2p-3, -0x1.248d2c8444ac6p-3 };
    const polynomial = (coefficients[0] + h_squared * coefficients[1]) +
        h_fourth * (coefficients[2] + h_squared * coefficients[3]);
    low = @mulAdd(f64, h, polynomial, low);
    const error_bound = h * 0x3.fp-52;
    const upper = (low + error_bound) + high;
    const lower = (low - error_bound) + high;
    if (upper == lower) return upper;
    return atanRefine(value, upper);
}
