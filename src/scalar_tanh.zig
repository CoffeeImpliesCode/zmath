// Correctly rounded hyperbolic tangent for binary16, binary32, and binary64 values.
//
// The finite kernels are ported from CORE-MATH revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa:
//   src/binary16/tanh/tanhf16.c (Copyright (c) 2025 Paul Zimmermann).
//   src/binary32/tanh/tanhf.c (Copyright (c) 2022-2025 Alexei Sibidanov).
//   src/binary64/tanh/tanh.c (Copyright (c) 2023-2026 Alexei Sibidanov,
//     Cyprien Peignier, and Paul Zimmermann).
//
// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/). It is used under the MIT License.
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

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;

fn signedOneH(value: f16) f16 {
    const bits: u16 = @bitCast(value);
    return if ((bits & sign_mask_h) == 0) 1.0 else -1.0;
}

fn signedOneF(value: f32) f32 {
    const bits: u32 = @bitCast(value);
    return if ((bits & sign_mask_f) == 0) 1.0 else -1.0;
}

fn quietNaNBits(comptime Bits: type, bits: Bits, quiet_bit: Bits) Bits {
    return bits | quiet_bit;
}

pub fn tanhh(value: f16) f16 {
    @setFloatMode(.strict);
    const input_bits: u16 = @bitCast(value);
    const magnitude = input_bits & magnitude_mask_h;
    if (magnitude > infinity_bits_h) return @bitCast(quietNaNBits(u16, input_bits, 0x0200));
    if (magnitude == infinity_bits_h) return signedOneH(value);
    if (magnitude >= 0x4482) return @floatCast(@as(f32, signedOneH(value)) * 0x1.fffff2p-1);
    if (magnitude <= 0x2744) {
        if (magnitude == 0) return value;
        const value_f: f32 = @floatCast(value);
        return @floatCast(@mulAdd(f32, value_f, -0x1p-12, value_f));
    }
    // The pinned binary32 kernel is exact after binary16 rounding for this range.
    return @floatCast(tanhf(@floatCast(value)));
}

pub fn tanhf(value: f32) f32 {
    @setFloatMode(.strict);
    const input_bits: u32 = @bitCast(value);
    const exponent = (input_bits >> 23) & 0xff;
    const magnitude = input_bits & magnitude_mask_f;
    if (exponent == 0xff) {
        if (magnitude != infinity_bits_f) {
            return @bitCast(quietNaNBits(u32, input_bits, 0x0040_0000));
        }
        return signedOneF(value);
    }
    if (exponent < 115) {
        if (exponent < 102) {
            if (magnitude == 0) return value;
            return @mulAdd(f32, -value, @abs(value), value);
        }
        const value_squared = value * value;
        return @mulAdd(f32, value, -0x1.555556p-2 * value_squared, value);
    }
    if ((magnitude << 1) > (0x4110_2cb3 << 1)) {
        return signedOneF(value) - signedOneF(value) * 0x1p-25;
    }
    const z: f64 = value;
    const z_squared = z * z;
    const z_fourth = z_squared * z_squared;
    const z_eighth = z_fourth * z_fourth;
    const numerator_coefficients = [_]f64{
        0x1p+0,                0x1.30877b8b72d33p-3,  0x1.694aa09ae9e5ep-8,
        0x1.4101377abb729p-14, 0x1.e0392b1db0018p-22, 0x1.2533756e546f7p-30,
        0x1.d62e5abe6ae8ap-41, 0x1.b06be534182dep-54,
    };
    const denominator_coefficients = [_]f64{
        0x1p+0,                0x1.ed99131b0ebeap-2,  0x1.0d27ed6c95a69p-5,
        0x1.7cbdaca0e9fccp-11, 0x1.b4e60b892578ep-18, 0x1.a6f707c5c71abp-26,
        0x1.35a8b6e2cd94cp-35, 0x1.ca8230677aa01p-47,
    };
    var numerator_0 = numerator_coefficients[0] + z_squared * numerator_coefficients[1];
    const numerator_2 = numerator_coefficients[2] + z_squared * numerator_coefficients[3];
    var numerator_4 = numerator_coefficients[4] + z_squared * numerator_coefficients[5];
    const numerator_6 = numerator_coefficients[6] + z_squared * numerator_coefficients[7];
    numerator_0 += z_fourth * numerator_2;
    numerator_4 += z_fourth * numerator_6;
    numerator_0 += z_eighth * numerator_4;
    var denominator_0 = denominator_coefficients[0] + z_squared * denominator_coefficients[1];
    const denominator_2 = denominator_coefficients[2] + z_squared * denominator_coefficients[3];
    var denominator_4 = denominator_coefficients[4] + z_squared * denominator_coefficients[5];
    const denominator_6 = denominator_coefficients[6] + z_squared * denominator_coefficients[7];
    denominator_0 += z_fourth * denominator_2;
    denominator_4 += z_fourth * denominator_6;
    denominator_0 += z_eighth * denominator_4;
    return @floatCast(z * numerator_0 / denominator_0);
}

fn fastTwoSum(x: f64, y: f64, remainder: *f64) f64 {
    const sum = x + y;
    const z = sum - x;
    remainder.* = y - z;
    return sum;
}

fn fastTwoSub(x: f64, y: f64, remainder: *f64) f64 {
    const sum = x - y;
    const z = x - sum;
    remainder.* = z - y;
    return sum;
}

fn mulDdAcc(xh: f64, xl: f64, ch: f64, cl: f64, low: *f64) f64 {
    const product_low_high = xl * ch;
    const product_high_low = xh * cl;
    const product_high = xh * ch;
    var product_high_rest = @mulAdd(f64, xh, ch, -product_high);
    product_high_rest += product_high_low + product_low_high;
    return fastTwoSum(product_high, product_high_rest, low);
}

fn mulDddAcc(xh: f64, xl: f64, c: f64, low: *f64) f64 {
    const product_low = xl * c;
    const product_high = xh * c;
    var product_high_rest = @mulAdd(f64, xh, c, -product_high);
    product_high_rest += product_low;
    return fastTwoSum(product_high, product_high_rest, low);
}

fn polyDd(xh: f64, xl: f64, comptime n: usize, coefficients: *const [n][2]f64, low: *f64) f64 {
    var index: usize = n - 1;
    // SAFETY: fastTwoSum writes sum_low before polyDd reads it.
    var sum_low: f64 = undefined;
    var sum_high = fastTwoSum(coefficients[index][0], low.*, &sum_low);
    var coefficient_low = sum_low + coefficients[index][1];
    while (index > 0) : (index -= 1) {
        sum_high = mulDdAcc(xh, xl, sum_high, coefficient_low, &coefficient_low);
        // SAFETY: fastTwoSum writes term_low before the value is read.
        var term_low: f64 = undefined;
        sum_high = fastTwoSum(coefficients[index - 1][0], sum_high, &term_low);
        coefficient_low = coefficient_low + coefficients[index - 1][1] + term_low;
    }
    low.* = coefficient_low;
    return sum_high;
}

fn asExpAccurate(x: f64, t: f64, th: f64, tl: f64, low: *f64) f64 {
    const coefficients = [_][2]f64{
        .{ 0x1p+0, 0x1.6c16bd194535dp-94 },
        .{ 0x1p-1, -0x1.8259d904fd34fp-93 },
        .{ 0x1.5555555555555p-3, 0x1.53e93e9f26e62p-57 },
    };
    const log_two_high = 0x1.62e42ffp-13;
    const log_two_low = 0x1.718432a1b0e26p-47;
    const log_two_low_low = 0x1.9ff0342542fc3p-102;
    const delta = x - log_two_high * t;
    const delta_low = log_two_low * t;
    const delta_low_low = log_two_low_low * t + @mulAdd(f64, log_two_low, t, -delta_low);
    const delta_high = delta + delta_low;
    const corrected_delta_low = (delta - delta_high) + delta_low + delta_low_low;
    var polynomial_low = delta_high * (0x1.5555555555555p-5 +
        delta_high * (0x1.11111113e93e9p-7 + delta_high * 0x1.6c16c169400a7p-10));
    var result_high = polyDd(
        delta_high,
        corrected_delta_low,
        3,
        &coefficients,
        &polynomial_low,
    );
    result_high = mulDdAcc(
        delta_high,
        corrected_delta_low,
        result_high,
        polynomial_low,
        &polynomial_low,
    );
    result_high = mulDdAcc(th, tl, result_high, polynomial_low, &polynomial_low);
    const z_high = th + result_high;
    const z_low = (th - z_high) + result_high;
    const u_high = z_high + tl;
    const u_low = (z_high - u_high) + tl + z_low;
    const v_high = u_high + polynomial_low;
    low.* = (u_high - v_high) + polynomial_low + u_low;
    return v_high;
}

fn signedOneD(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    return if ((bits & sign_mask_d) == 0) 1.0 else -1.0;
}

const tanhd_t0 = [_][2]f64{
    .{ 0x0p+0, 0x1p+0 },
    .{ -0x1.19083535b085ep-56, 0x1.02c9a3e778061p+0 },
    .{ 0x1.d73e2a475b466p-55, 0x1.059b0d3158574p+0 },
    .{ 0x1.186be4bb285p-57, 0x1.0874518759bc8p+0 },
    .{ 0x1.8a62e4adc610ap-54, 0x1.0b5586cf9890fp+0 },
    .{ 0x1.03a1727c57b52p-59, 0x1.0e3ec32d3d1a2p+0 },
    .{ -0x1.6c51039449b3ap-54, 0x1.11301d0125b51p+0 },
    .{ -0x1.32fbf9af1369ep-54, 0x1.1429aaea92dep+0 },
    .{ -0x1.19041b9d78a76p-55, 0x1.172b83c7d517bp+0 },
    .{ 0x1.e5b4c7b4968e4p-55, 0x1.1a35beb6fcb75p+0 },
    .{ 0x1.e016e00a2643cp-54, 0x1.1d4873168b9aap+0 },
    .{ 0x1.dc775814a8494p-55, 0x1.2063b88628cd6p+0 },
    .{ 0x1.9b07eb6c70572p-54, 0x1.2387a6e756238p+0 },
    .{ 0x1.2bd339940e9dap-55, 0x1.26b4565e27cddp+0 },
    .{ 0x1.612e8afad1256p-55, 0x1.29e9df51fdee1p+0 },
    .{ 0x1.0024754db41d4p-54, 0x1.2d285a6e4030bp+0 },
    .{ 0x1.6f46ad23182e4p-55, 0x1.306fe0a31b715p+0 },
    .{ 0x1.32721843659a6p-54, 0x1.33c08b26416ffp+0 },
    .{ -0x1.63aeabf42eae2p-54, 0x1.371a7373aa9cbp+0 },
    .{ -0x1.5e436d661f5e2p-56, 0x1.3a7db34e59ff7p+0 },
    .{ 0x1.ada0911f09ebcp-55, 0x1.3dea64c123422p+0 },
    .{ -0x1.ef3691c309278p-58, 0x1.4160a21f72e2ap+0 },
    .{ 0x1.89b7a04ef80dp-59, 0x1.44e086061892dp+0 },
    .{ 0x1.3c1a3b69062fp-56, 0x1.486a2b5c13cdp+0 },
    .{ 0x1.d4397afec42e2p-56, 0x1.4bfdad5362a27p+0 },
    .{ -0x1.4b309d25957e4p-54, 0x1.4f9b2769d2ca7p+0 },
    .{ -0x1.07abe1db13cacp-55, 0x1.5342b569d4f82p+0 },
    .{ 0x1.9bb2c011d93acp-54, 0x1.56f4736b527dap+0 },
    .{ 0x1.6324c054647acp-54, 0x1.5ab07dd485429p+0 },
    .{ 0x1.ba6f93080e65ep-54, 0x1.5e76f15ad2148p+0 },
    .{ -0x1.383c17e40b496p-54, 0x1.6247eb03a5585p+0 },
    .{ -0x1.bb60987591c34p-54, 0x1.6623882552225p+0 },
    .{ -0x1.bdd3413b26456p-54, 0x1.6a09e667f3bcdp+0 },
    .{ -0x1.bbe3a683c88aap-57, 0x1.6dfb23c651a2fp+0 },
    .{ -0x1.16e4786887a9ap-55, 0x1.71f75e8ec5f74p+0 },
    .{ -0x1.0245957316dd4p-54, 0x1.75feb564267c9p+0 },
    .{ -0x1.41577ee04993p-55, 0x1.7a11473eb0187p+0 },
    .{ 0x1.05d02ba15797ep-56, 0x1.7e2f336cf4e62p+0 },
    .{ -0x1.d4c1dd41532d8p-54, 0x1.82589994cce13p+0 },
    .{ -0x1.fc6f89bd4f6bap-54, 0x1.868d99b4492edp+0 },
    .{ 0x1.6e9f156864b26p-54, 0x1.8ace5422aa0dbp+0 },
    .{ 0x1.5cc13a2e3976cp-55, 0x1.8f1ae99157736p+0 },
    .{ -0x1.75fc781b57ebcp-57, 0x1.93737b0cdc5e5p+0 },
    .{ -0x1.d185b7c1b85dp-54, 0x1.97d829fde4e5p+0 },
    .{ 0x1.c7c46b071f2bep-56, 0x1.9c49182a3f09p+0 },
    .{ -0x1.359495d1cd532p-54, 0x1.a0c667b5de565p+0 },
    .{ -0x1.d2f6edb8d41e2p-54, 0x1.a5503b23e255dp+0 },
    .{ 0x1.0fac90ef7fd32p-54, 0x1.a9e6b5579fdbfp+0 },
    .{ 0x1.7a1cd345dcc82p-54, 0x1.ae89f995ad3adp+0 },
    .{ -0x1.2805e3084d708p-57, 0x1.b33a2b84f15fbp+0 },
    .{ -0x1.5584f7e54ac3ap-56, 0x1.b7f76f2fb5e47p+0 },
    .{ 0x1.23dd07a2d9e84p-55, 0x1.bcc1e904bc1d2p+0 },
    .{ 0x1.11065895048dep-55, 0x1.c199bdd85529cp+0 },
    .{ 0x1.2884dff483cacp-54, 0x1.c67f12e57d14bp+0 },
    .{ 0x1.503cbd1e949dcp-56, 0x1.cb720dcef9069p+0 },
    .{ -0x1.cbc3743797a9cp-54, 0x1.d072d4a07897cp+0 },
    .{ 0x1.2ed02d75b3706p-55, 0x1.d5818dcfba487p+0 },
    .{ 0x1.c2300696db532p-54, 0x1.da9e603db3285p+0 },
    .{ -0x1.1a5cd4f184b5cp-54, 0x1.dfc97337b9b5fp+0 },
    .{ 0x1.39e8980a9cc9p-55, 0x1.e502ee78b3ff6p+0 },
    .{ -0x1.e9c23179c2894p-54, 0x1.ea4afa2a490dap+0 },
    .{ 0x1.dc7f486a4b6bp-54, 0x1.efa1bee615a27p+0 },
    .{ 0x1.9d3e12dd8a18ap-54, 0x1.f50765b6e454p+0 },
    .{ 0x1.74853f3a5931ep-55, 0x1.fa7c1819e90d8p+0 },
};

const tanhd_t1 = [_][2]f64{
    .{ 0x0p+0, 0x1p+0 },
    .{ 0x1.ae8e38c59c72ap-54, 0x1.000b175effdc7p+0 },
    .{ -0x1.7b5d0d58ea8f4p-58, 0x1.00162f3904052p+0 },
    .{ 0x1.4115cb6b16a8ep-54, 0x1.0021478e11ce6p+0 },
    .{ -0x1.d7c96f201bb2ep-55, 0x1.002c605e2e8cfp+0 },
    .{ 0x1.84711d4c35eap-54, 0x1.003779a95f959p+0 },
    .{ -0x1.0484245243778p-55, 0x1.0042936faa3d8p+0 },
    .{ -0x1.4b237da2025fap-54, 0x1.004dadb113dap+0 },
    .{ -0x1.5e00e62d6b30ep-56, 0x1.0058c86da1c0ap+0 },
    .{ 0x1.a1d6cedbb948p-54, 0x1.0063e3a559473p+0 },
    .{ -0x1.4acf197a00142p-54, 0x1.006eff583fc3dp+0 },
    .{ -0x1.eaf2ea42391a6p-57, 0x1.007a1b865a8cap+0 },
    .{ 0x1.da93f90835f76p-56, 0x1.0085382faef83p+0 },
    .{ -0x1.6a79084ab093cp-55, 0x1.00905554425d4p+0 },
    .{ 0x1.86364f8fbe8f8p-54, 0x1.009b72f41a12bp+0 },
    .{ -0x1.82e8e14e3110ep-55, 0x1.00a6910f3b6fdp+0 },
    .{ -0x1.4f6b2a7609f72p-55, 0x1.00b1afa5abcbfp+0 },
    .{ -0x1.e1a258ea8f71ap-56, 0x1.00bcceb7707ecp+0 },
    .{ 0x1.4362ca5bc26f2p-56, 0x1.00c7ee448ee02p+0 },
    .{ 0x1.095a56c919d02p-54, 0x1.00d30e4d0c483p+0 },
    .{ -0x1.406ac4e81a646p-57, 0x1.00de2ed0ee0f5p+0 },
    .{ 0x1.b5a6902767e08p-54, 0x1.00e94fd0398ep+0 },
    .{ -0x1.91b206085932p-54, 0x1.00f4714af41d3p+0 },
    .{ 0x1.427068ab22306p-55, 0x1.00ff93412315cp+0 },
    .{ 0x1.c1d0660524e08p-54, 0x1.010ab5b2cbd11p+0 },
    .{ -0x1.e7bdfb3204be8p-54, 0x1.0115d89ff3a8bp+0 },
    .{ 0x1.843aa8b9cbbc6p-55, 0x1.0120fc089ff63p+0 },
    .{ -0x1.34104ee7edae8p-56, 0x1.012c1fecd613bp+0 },
    .{ -0x1.2b6aeb6176892p-56, 0x1.0137444c9b5b5p+0 },
    .{ 0x1.a8cd33b8a1bb2p-56, 0x1.01426927f5278p+0 },
    .{ 0x1.2edc08e5da99ap-56, 0x1.014d8e7ee8d2fp+0 },
    .{ 0x1.57ba2dc7e0c72p-55, 0x1.0158b4517bb88p+0 },
    .{ 0x1.b61299ab8cdb8p-54, 0x1.0163da9fb3335p+0 },
    .{ -0x1.90565902c5f44p-54, 0x1.016f0169949edp+0 },
    .{ 0x1.70fc41c5c2d54p-55, 0x1.017a28af25567p+0 },
    .{ 0x1.4b9a6e145d76cp-54, 0x1.018550706ab62p+0 },
    .{ -0x1.008eff5142bfap-56, 0x1.019078ad6a19fp+0 },
    .{ -0x1.77669f033c7dep-54, 0x1.019ba16628de2p+0 },
    .{ -0x1.09bb78eeead0ap-54, 0x1.01a6ca9aac5f3p+0 },
    .{ 0x1.371231477ece6p-54, 0x1.01b1f44af9f9ep+0 },
    .{ 0x1.5e7626621eb5ap-56, 0x1.01bd1e77170b4p+0 },
    .{ -0x1.bc72b100828a4p-54, 0x1.01c8491f08f08p+0 },
    .{ -0x1.ce39cbbab8bbep-57, 0x1.01d37442d507p+0 },
    .{ 0x1.16996709da2e2p-55, 0x1.01de9fe280ac8p+0 },
    .{ -0x1.c11f5239bf536p-55, 0x1.01e9cbfe113efp+0 },
    .{ 0x1.e1d4eb5edc6b4p-55, 0x1.01f4f8958c1c6p+0 },
    .{ -0x1.afb99946ee3fp-54, 0x1.020025a8f6a35p+0 },
    .{ -0x1.8f06d8a148a32p-54, 0x1.020b533856324p+0 },
    .{ -0x1.2bf310fc54eb6p-55, 0x1.02168143b0281p+0 },
    .{ -0x1.c95a035eb4176p-54, 0x1.0221afcb09e3ep+0 },
    .{ -0x1.491793e46834cp-54, 0x1.022cdece68c4fp+0 },
    .{ -0x1.3e8d0d9c4909p-56, 0x1.02380e4dd22adp+0 },
    .{ -0x1.314aa16278aa4p-54, 0x1.02433e494b755p+0 },
    .{ 0x1.48daf888e965p-55, 0x1.024e6ec0da046p+0 },
    .{ 0x1.56dc8046821f4p-55, 0x1.02599fb483385p+0 },
    .{ 0x1.45b42356b9d46p-54, 0x1.0264d1244c719p+0 },
    .{ -0x1.082ef51b61d7ep-56, 0x1.027003103b10ep+0 },
    .{ 0x1.2106ed0920a34p-56, 0x1.027b357854772p+0 },
    .{ -0x1.fd4cf26ea5d0ep-54, 0x1.0286685c9e059p+0 },
    .{ -0x1.09f8775e78084p-54, 0x1.02919bbd1d1d8p+0 },
    .{ 0x1.64cbba902ca28p-58, 0x1.029ccf99d720ap+0 },
    .{ 0x1.4383ef231d206p-54, 0x1.02a803f2d170dp+0 },
    .{ 0x1.4a47a505b3a46p-54, 0x1.02b338c811703p+0 },
    .{ 0x1.e47120223468p-54, 0x1.02be6e199c811p+0 },
};

fn asTanhZero(value: f64) f64 {
    const high_coefficients = [_][2]f64{
        .{ -0x1.5555555555555p-2, -0x1.5555555554cc4p-56 },
        .{ 0x1.1111111111111p-3, 0x1.111110f8c0178p-59 },
        .{ -0x1.ba1ba1ba1ba1cp-5, 0x1.7917c1d676ff5p-59 },
        .{ 0x1.664f4882c10fap-6, -0x1.9d5cb27c0af28p-63 },
        .{ -0x1.226e355e6c23cp-7, -0x1.c9674586913f3p-61 },
        .{ 0x1.d6d3d0e157db3p-9, -0x1.71376fa06ce94p-65 },
        .{ -0x1.7da36452b5e46p-10, -0x1.aba8d51bd9cp-65 },
        .{ 0x1.3558247faa32dp-11, -0x1.e0cfb423aedfdp-65 },
        .{ -0x1.f57d76ea30928p-13, -0x1.c30601213cae9p-67 },
    };
    const low_coefficients = [_]f64{
        0x1.967e0a63ca836p-14,  -0x1.497b99d2a77d1p-15, 0x1.0ae346258cbdep-16,
        -0x1.aade68fb2f076p-18, 0x1.22e609bf8671fp-19,
    };
    const value_squared = value * value;
    const value_squared_low = @mulAdd(f64, value, value, -value_squared);
    var correction = value_squared * (low_coefficients[0] + value_squared *
        (low_coefficients[1] + value_squared * (low_coefficients[2] +
            value_squared * (low_coefficients[3] + value_squared * low_coefficients[4]))));
    var high = polyDd(value_squared, value_squared_low, 9, &high_coefficients, &correction);
    high = mulDddAcc(high, correction, value, &correction);
    high = mulDdAcc(high, correction, value_squared, value_squared_low, &correction);
    // SAFETY: fastTwoSum writes low before the value is read.
    var low: f64 = undefined;
    const result_high = fastTwoSum(value, high, &low);
    high = fastTwoSum(low, correction, &correction);
    var result_bits: u64 = @bitCast(high);
    if ((result_bits & 0x000f_ffff_ffff_ffff) == 0) {
        const correction_bits: u64 = @bitCast(correction);
        if (((correction_bits ^ result_bits) >> 63) != 0) result_bits -|= 1 else result_bits +|= 1;
        high = @bitCast(result_bits);
        if (@abs(value) == 0x1.ac343b179fec4p-3) return asTanhDatabase(value, result_high + high);
    }
    return result_high + high;
}

fn asTanhDatabase(value: f64, fallback: f64) f64 {
    const database = [_][3]f64{
        .{ 0x1.ac343b179fec4p-3, 0x1.a612499c53078p-3, 0x1p-57 },
        .{ 0x1.00764a988bf73p-2, 0x1.f676484c0703bp-3, -0x1p-104 },
        .{ 0x1.17d1e8a63711fp-2, 0x1.110e96a6c2d96p-2, -0x1p-105 },
        .{ 0x1.291c601a05276p-2, 0x1.210b7d0c03743p-2, 0x1p-56 },
        .{ 0x1.36f33d51c264dp-2, 0x1.2dbb7b1c91363p-2, -0x1p-106 },
        .{ 0x1.43eaea23649c3p-2, 0x1.39877ed028641p-2, -0x1p-55 },
        .{ 0x1.d88d7550b2826p-2, 0x1.b9a3637366afdp-2, 0x1p-56 },
        .{ 0x1.e611aa58ab608p-2, 0x1.c493dc899e4a6p-2, -0x1p-55 },
        .{ 0x1.01efe7ac8c15dp-1, 0x1.dc3fe1b524821p-2, -0x1p-104 },
        .{ 0x1.1005ec0bccabbp-1, 0x1.f20b1c8557dedp-2, -0x1p-55 },
        .{ 0x1.33dfeb0fa4bfep-1, 0x1.1372f9ee76e99p-1, 0x1p-55 },
        .{ 0x1.49f24ac5cac35p-1, 0x1.22c495ff06104p-1, -0x1p-104 },
    };
    var first: i32 = 0;
    var last: i32 = database.len - 1;
    while (first <= last) {
        const middle: i32 = @divTrunc(first + last, 2);
        const row = database[@intCast(middle)];
        const absolute = @abs(value);
        if (row[0] < absolute) first = middle + 1 else if (row[0] == absolute) {
            const sign = signedOneD(value);
            return sign * row[1] + sign * row[2];
        } else last = middle - 1;
    }
    return fallback;
}

pub fn tanhd(value: f64) f64 {
    @setFloatMode(.strict);
    const input_bits: u64 = @bitCast(value);
    const magnitude = input_bits & magnitude_mask_d;
    if (magnitude >= 0x4033_0fc1_931f_09ca) {
        if (magnitude > infinity_bits_d) return @bitCast(input_bits | 0x0008_0000_0000_0000);
        const sign = signedOneD(value);
        if (magnitude == infinity_bits_d) return sign;
        return sign * (1.0 - 0x1p-55);
    }
    const scale = -0x1.71547652b82fep+13;
    const rounded = @mulAdd(f64, @as(f64, @bitCast(magnitude)), scale, 0x1.8000004p+25);
    const rounded_bits: u64 = @bitCast(rounded);
    const table_bits = rounded_bits & ~@as(u64, (1 << 27) - 1);
    const table_offset = @as(f64, @bitCast(table_bits)) - 0x1.8p25;
    const index_one: usize = @intCast((rounded_bits >> 27) & 0x3f);
    const index_zero: usize = @intCast((rounded_bits >> 33) & 0x3f);
    const exponent_bits: i64 = @as(i64, @bitCast(rounded_bits << 13)) >> 52;
    const scale_bits: u64 = @as(u64, @intCast(1023 + exponent_bits)) << 52;
    const scale_factor: f64 = @bitCast(scale_bits);
    const ch = [_]f64{ 0x1p+1, 0x1p+1, 0x1.55555557e54ffp+0, 0x1.55555553a12f4p-1 };
    const table_zero_high = tanhd_t0[index_zero][1];
    const table_one_high = tanhd_t1[index_one][1];
    var table_high = table_zero_high * table_one_high;
    // SAFETY: Each non-returning branch assigns table_low before final use.
    var table_low: f64 = undefined;
    if (magnitude < 0x400d_76c8_b439_5810) {
        if (magnitude < 0x3fd0_0000_0000_0000) {
            if (magnitude <= 0x3e4d_12ed_0af1_a27f) {
                if (magnitude == 0) return value;
                return @mulAdd(f64, value, -0x1p-55, value);
            }
            const coefficients = [_]f64{
                -0x1.5555555555555p-2,  0x1.1111111110f33p-3,  -0x1.ba1ba1b9b8ea6p-5,
                0x1.664f4838e0a43p-6,   -0x1.226e17d1bc09bp-7, 0x1.d6c64dfba2565p-9,
                -0x1.7bdd094d327afp-10, 0x1.1535ad0c31d0ep-11,
            };
            const value_squared = value * value;
            const value_cubed = value_squared * value;
            const value_fourth = value_squared * value_squared;
            const value_eighth = value_fourth * value_fourth;
            var polynomial_one = coefficients[4] + value_squared * coefficients[5];
            polynomial_one += value_fourth * (coefficients[6] + value_squared * coefficients[7]);
            var polynomial_zero = coefficients[0] + value_squared * coefficients[1];
            polynomial_zero += value_fourth * (coefficients[2] + value_squared * coefficients[3]);
            polynomial_zero += value_eighth * polynomial_one;
            polynomial_zero *= value_cubed;
            // SAFETY: fastTwoSum writes polynomial_low before the value is read.
            var polynomial_low: f64 = undefined;
            const polynomial_high = fastTwoSum(value, polynomial_zero, &polynomial_low);
            const bound = value_cubed * 0x1.ap-52;
            const lower = polynomial_high + (polynomial_low - bound);
            const upper = polynomial_high + (polynomial_low + bound);
            if (lower == upper) return lower;
            return asTanhZero(value);
        }
        const table_zero_low = tanhd_t0[index_zero][0];
        const table_one_low = tanhd_t1[index_one][0];
        table_low = table_zero_high * table_one_low + table_one_high * table_zero_low +
            @mulAdd(f64, table_zero_high, table_one_high, -table_high);
        table_high *= scale_factor;
        table_low *= scale_factor;
        const log_two_high = -0x1.62e42ffp-14;
        const log_two_low = -0x1.718432a1b0e26p-48;
        const delta = (log_two_high * table_offset - @abs(value)) - log_two_low * table_offset;
        const delta_squared = delta * delta;
        const polynomial = delta * ((ch[0] + delta * ch[1]) +
            delta_squared * (ch[2] + delta * ch[3]));
        var high = table_high;
        var low = table_low + high * polynomial;
        high = fastTwoSum(high, low, &low);
        var product_high = high;
        var product_low = low;
        var reciprocal_high = high;
        var reciprocal_low = low;
        // SAFETY: fastTwoSum writes correction before reciprocal_low reads it.
        var correction: f64 = undefined;
        reciprocal_high = fastTwoSum(1.0, reciprocal_high, &correction);
        reciprocal_low += correction;
        const reciprocal = 1.0 / reciprocal_high;
        reciprocal_low = (reciprocal_low * reciprocal +
            @mulAdd(f64, reciprocal, reciprocal_high, -1.0)) * -reciprocal;
        product_high = mulDdAcc(
            product_high,
            product_low,
            reciprocal,
            reciprocal_low,
            &product_low,
        );
        const bound = high * 0x1.0dp-62;
        high = fastTwoSub(0.5, product_high, &low);
        low -= product_low;
        high *= signedOneD(value) * 2.0;
        low *= signedOneD(value) * 2.0;
        const lower = high + (low - bound);
        const upper = high + (low + bound);
        if (lower == upper) return lower;
    } else {
        const log_two = -0x1.62e42fefa39efp-14;
        const delta = @mulAdd(f64, log_two, table_offset, -@abs(value));
        const delta_squared = delta * delta;
        const polynomial = delta * ((ch[0] + delta * ch[1]) +
            delta_squared * (ch[2] + delta * ch[3]));
        var high = table_high * scale_factor;
        high += (polynomial + (2.0 * 0x1.3p-55 * @abs(value))) * high;
        const bound = high * 0x1.1p-49;
        high = (2.0 * high) / (1.0 + high);
        const sign = signedOneD(value);
        high = sign * high;
        const lower = sign - (high + bound);
        const upper = sign - (high - bound);
        if (lower == upper) return lower;
        const table_zero_low = tanhd_t0[index_zero][0];
        const table_one_low = tanhd_t1[index_one][0];
        table_low = table_zero_high * table_one_low + table_one_high * table_zero_low +
            @mulAdd(f64, table_zero_high, table_one_high, -table_high);
        table_high *= scale_factor;
        table_low *= scale_factor;
    }
    // SAFETY: asExpAccurate writes exponential_low before the value is read.
    var exponential_low: f64 = undefined;
    const exponential_high = asExpAccurate(
        -2.0 * @abs(value),
        table_offset,
        table_high,
        table_low,
        &exponential_low,
    );
    // SAFETY: fastTwoSum writes correction before denominator_low reads it.
    var correction: f64 = undefined;
    var denominator_high = fastTwoSum(1.0, exponential_high, &correction);
    var denominator_low = exponential_low + correction;
    denominator_high = fastTwoSum(denominator_high, denominator_low, &denominator_low);
    const reciprocal = 1.0 / denominator_high;
    const reciprocal_low = (denominator_low * reciprocal +
        @mulAdd(f64, reciprocal, denominator_high, -1.0)) * -reciprocal;
    // SAFETY: mulDdAcc writes product_low before the value is read.
    var product_low: f64 = undefined;
    const product_high = mulDdAcc(
        exponential_high,
        exponential_low,
        reciprocal,
        reciprocal_low,
        &product_low,
    );
    // SAFETY: fastTwoSub writes result_low before the value is read.
    var result_low: f64 = undefined;
    var result_high = fastTwoSub(0.5, product_high, &result_low);
    result_low -= product_low;
    result_high = fastTwoSum(result_high, result_low, &result_low);
    const result = signedOneD(value) * 2.0 * result_high + signedOneD(value) * 2.0 * result_low;
    const low_bits: u64 = @bitCast(result_low);
    if (((low_bits + 32) & 0x000f_ffff_ffff_ffff) < 65) return asTanhDatabase(value, result);
    return result;
}
