// Correctly rounded inverse cosine for binary16, binary32, and binary64.
//
// The finite kernels and tables are ported from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa:
//   src/binary16/acos/acosf16.c (Copyright 2025 Paul Zimmermann).
//   src/binary32/acos/acosf.c (Copyright 2023 Alexei Sibidanov).
//   src/binary64/acos/acos.c (Copyright 2024-2025 Alexei Sibidanov).
//
// The upstream code is distributed under the MIT License:
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

const fma_impl = @import("scalar_fma.zig");

const sign_mask_h: u16 = 0x8000;
const magnitude_mask_h: u16 = 0x7fff;
const infinity_bits_h: u16 = 0x7c00;
const quiet_nan_bit_h: u16 = 0x0200;
const canonical_nan_bits_h: u16 = 0x7e00;

const sign_mask_f: u32 = 0x8000_0000;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const infinity_bits_f: u32 = 0x7f80_0000;
const quiet_nan_bit_f: u32 = 0x0040_0000;
const canonical_nan_bits_f: u32 = 0x7fc0_0000;

const sign_mask_d: u64 = 0x8000_0000_0000_0000;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit_d: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits_d: u64 = 0x7ff8_0000_0000_0000;

fn absD(value: f64) f64 {
    return if (value < 0.0) -value else value;
}

fn copySignD(value: f64, sign: f64) f64 {
    const value_bits: u64 = @bitCast(value);
    const sign_bits: u64 = @bitCast(sign);
    return @bitCast((value_bits & ~sign_mask_d) | (sign_bits & sign_mask_d));
}

fn roundEvenFinite(value: f64) i64 {
    @setFloatMode(.strict);
    // Both reductions and the refinement table supply values in [0, 33).
    const rounded_bits: i64 = @bitCast(value + 0x1.8p52);
    return rounded_bits -% 0x4338_0000_0000_0000;
}

fn fastTwoSum(x: f64, y: f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    const total = x + y;
    const rounded_y = total - x;
    error_value.* = y - rounded_y;
    return total;
}

fn twoSum(x: f64, y: f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    const total = x + y;
    const x_prime = total - y;
    const y_prime = total - x_prime;
    const delta_x = x - x_prime;
    const delta_y = y - y_prime;
    error_value.* = delta_x + delta_y;
    return total;
}

fn fastSum(x_high: f64, x_low: f64, y_high: f64, y_low: f64, error_value: *f64) f64 {
    @setFloatMode(.strict);
    var low: f64 = undefined;
    const high = fastTwoSum(x_high, y_high, &low);
    error_value.* = (x_low + y_low) + low;
    return high;
}

fn sum(x_high: f64, x_low: f64, coefficient_high: f64, coefficient_low: f64, low: *f64) f64 {
    @setFloatMode(.strict);
    var rounded_low: f64 = undefined;
    const high = twoSum(x_high, coefficient_high, &rounded_low);
    low.* = (x_low + coefficient_low) + rounded_low;
    return high;
}

fn mulDd(x_high: f64, x_low: f64, coefficient_high: f64, coefficient_low: f64, low: *f64) f64 {
    @setFloatMode(.strict);
    const high_product = x_high * coefficient_high;
    low.* = (x_high * coefficient_low + x_low * coefficient_high) +
        fma_impl.fma64(x_high, coefficient_high, -high_product);
    return high_product;
}

fn polyDd(
    x_high: f64,
    x_low: f64,
    count: usize,
    coefficients: []const [2]f64,
    low: *f64,
) f64 {
    @setFloatMode(.strict);
    var index = count - 1;
    var correction: f64 = undefined;
    var high = fastTwoSum(coefficients[index][0], low.*, &correction);
    var coefficient_low = coefficients[index][1] + correction;
    while (index > 0) {
        index -= 1;
        high = mulDd(x_high, x_low, high, coefficient_low, &coefficient_low);
        high = fastSum(
            coefficients[index][0],
            coefficients[index][1],
            high,
            coefficient_low,
            &coefficient_low,
        );
    }
    low.* = coefficient_low;
    return high;
}

pub fn acosh(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude = bits & magnitude_mask_h;
    if (magnitude > infinity_bits_h) return @bitCast(bits | quiet_nan_bit_h);
    if (magnitude >= infinity_bits_h or magnitude > 0x3c00) return @bitCast(canonical_nan_bits_h);
    if (magnitude == 0x3c00) {
        if ((bits & sign_mask_h) == 0) return @bitCast(@as(u16, 0));
        return @bitCast(@as(u16, 0x4248));
    }
    if (magnitude == 0) return @bitCast(@as(u16, 0x3e48));

    const input_f: f32 = @floatCast(value);
    const input_bits: u32 = @bitCast(input_f);
    const input_magnitude = input_bits & magnitude_mask_f;
    const negative = (bits & sign_mask_h) != 0;
    const t = if (negative) -input_f else input_f;
    const squared = t * t;
    var result: f32 = undefined;
    if (input_magnitude < 0x3e80_0000) {
        const first = fma_impl.fma(-0x1.25e6b8p-10, t, -0x1.fffb44p-1);
        const third = fma_impl.fma(-0x1.a85b22p-5, t, -0x1.3cc114p-3);
        result = fma_impl.fma(fma_impl.fma(third, squared, first), t, 0x1.921fb4p+0);
        if (negative) {
            result = 0x1.921fb4p+1 - result;
            if (input_magnitude >= 0x3696_0000 and input_magnitude <= 0x369a_0000) {
                result = 0x1.922002p+0;
            }
            if (input_bits == 0xbcf8_0000) result = 0x1.99e002p+0;
            if (input_bits == 0xbcfc_0000) result = 0x1.9a0006p+0;
        }
    } else if (input_magnitude < 0x3f00_0000) {
        const first = fma_impl.fma(-0x1.bd043ap-4, t, -0x1.f515cap-1);
        const third = fma_impl.fma(-0x1.190806p-2, t, 0x1.7e2d5ap-4);
        result = fma_impl.fma(fma_impl.fma(third, squared, first), t, 0x1.91b678p+0);
        if (negative) result = 0x1.921fb4p+1 - result;
    } else {
        const first = fma_impl.fma(0x1.31640cp-4, t, -0x1.ae5c5ep-3);
        const third = fma_impl.fma(0x1.251b5p-8, t, -0x1.98038p-6);
        const polynomial = fma_impl.fma(fma_impl.fma(third, squared, first), t, 0x1.91fa1cp+0);
        result = @sqrt(1.0 - t) * polynomial;
        if (negative) result = 0x1.921fb4p+1 - result;
    }
    return @floatCast(result);
}

fn poly12(value: f64, coefficients: []const f64) f64 {
    @setFloatMode(.strict);
    const squared = value * value;
    const fourth = squared * squared;
    const coefficient_0 = coefficients[0] + value * coefficients[1];
    const coefficient_2 = coefficients[2] + value * coefficients[3];
    const coefficient_4 = coefficients[4] + value * coefficients[5];
    const coefficient_6 = coefficients[6] + value * coefficients[7];
    const coefficient_8 = coefficients[8] + value * coefficients[9];
    const coefficient_10 = coefficients[10] + value * coefficients[11];
    const low_0 = coefficient_0 + coefficient_2 * squared;
    const low_4 = coefficient_4 + coefficient_6 * squared;
    const low_8 = coefficient_8 + squared * coefficient_10;
    return low_0 + fourth * (low_4 + fourth * low_8);
}

pub fn acosf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude = bits & magnitude_mask_f;
    if (magnitude > infinity_bits_f) return @bitCast(bits | quiet_nan_bit_f);
    if (magnitude >= infinity_bits_f or magnitude > 0x3f80_0000) return @bitCast(canonical_nan_bits_f);
    if (magnitude == 0x3f80_0000) {
        if ((bits & sign_mask_f) == 0) return @bitCast(@as(u32, 0));
        return @bitCast(@as(u32, 0x4049_0fdb));
    }
    if (magnitude == 0) return @bitCast(@as(u32, 0x3fc9_0fdb));

    const input_wide: f64 = @floatCast(value);
    const doubled_magnitude = magnitude << 1;
    if (doubled_magnitude < 0x7ec2_a1dc) {
        if (doubled_magnitude < 0x4000_0000) {
            const pi2_high: f32 = 0x1.921fb6p+0;
            const pi2_low: f32 = -0x1.777a5cp-25;
            return pi2_high + pi2_low;
        }
        const coefficients = [_]f64{
            0x1.fffffffd9ccb8p-1,  0x1.5555c94838007p-3,  0x1.32ded4b7c20fap-4,
            0x1.8566df703309ep-5,  -0x1.980c959bec9a3p-6, 0x1.56fbb04998344p-1,
            -0x1.403d8e4c49f52p+2, 0x1.b06c3e9f311eap+4,  -0x1.9ea97c4e2c21fp+6,
            0x1.200b8261cc61bp+8,  -0x1.2274c2799a5c7p+9, 0x1.a558a59cc19d3p+9,
            -0x1.aca4b6a529ffp+9,  0x1.228744703f813p+9,  -0x1.d7dbb0b322228p+7,
            0x1.5c2018c0c0105p+5,
        };
        const squared = input_wide * input_wide;
        const fourth = squared * squared;
        const eighth = fourth * fourth;
        const sixteenth = eighth * eighth;
        const first = (coefficients[0] + squared * coefficients[1]) +
            fourth * (coefficients[2] + squared * coefficients[3]);
        const second = (coefficients[4] + squared * coefficients[5]) +
            fourth * (coefficients[6] + squared * coefficients[7]);
        const third = (coefficients[8] + squared * coefficients[9]) +
            fourth * (coefficients[10] + squared * coefficients[11]);
        const fourth_part = (coefficients[12] + squared * coefficients[13]) +
            fourth * (coefficients[14] + squared * coefficients[15]);
        const result = input_wide * ((first + eighth * second) + sixteenth * (third + eighth * fourth_part));
        const upper: f32 = @floatCast(0x1.921fb54574191p+0 - result);
        const lower: f32 = @floatCast(0x1.921fb543118ap+0 - result);
        if (upper == lower) return upper;
    }
    if (bits == 0x3288_85a3) return 0x1.921fb6p+0 + 0x1p-25;
    if (bits == 0x3982_6222) return 0x1.920f6ap+0 + 0x1p-25;

    if (doubled_magnitude < 0x7e00_0000) {
        const coefficients = [_]f64{
            0x1.555555555529cp-3, 0x1.333333337e0ddp-4,  0x1.6db6db3b4465ep-5,
            0x1.f1c72e13ac306p-6, 0x1.6e89cebe06bc4p-6,  0x1.1c6dcf5289094p-6,
            0x1.c6dbbcc7c6315p-7, 0x1.8f8dc2615e996p-7,  0x1.a5833b7bf15e8p-8,
            0x1.43f44ace1665cp-6, -0x1.0fb17df881c73p-6, 0x1.07520c026b2d6p-5,
        };
        const squared = input_wide * input_wide;
        const correction = (input_wide * squared) * poly12(squared, &coefficients);
        const result = (0x1.921fb54442d18p+0 - input_wide) - correction;
        return @floatCast(result);
    }
    const absolute_value = absD(input_wide);
    const reduced = 1.0 - absolute_value;
    const root = @sqrt(reduced);
    const coefficients = [_]f64{
        0x1.6a09e667f3bcbp+0,  0x1.e2b7dddff2db9p-4,   0x1.b27247ab42dbcp-6,
        0x1.02995cc4e0744p-7,  0x1.5ffb0276ec8eap-9,   0x1.033885a928decp-10,
        0x1.911f2be23f8c7p-12, 0x1.4c3c55d2437fdp-13,  0x1.af477e1d7b461p-15,
        0x1.abd6bdff67dcbp-15, -0x1.1717e86d0fa28p-16, 0x1.6ff526de46023p-16,
    };
    const result = root * poly12(reduced, &coefficients);
    return @floatCast(if ((bits & sign_mask_f) == 0) result else 0x1.921fb54442d18p+1 - result);
}

const acos_coefficients = [33][8]f64{
    .{
        1,                    0,                    0x1.5555555555555p-3, 0x1.33333333333e4p-4,
        0x1.6db6db6d31f82p-5, 0x1.f1c71f6889397p-6, 0x1.6e874b7045b46p-6, 0x1.1f753132271e2p-6,
    },
    .{
        0x1.0055a27e0d033p+0, -0x1.d9ba10494c062p-54, 0x1.57c00cb5d6c4dp-3, 0x1.37881f5649a75p-4,
        0x1.759af49d494ddp-5, 0x1.002e1864dda2ep-5,   0x1.7c2d5d468cdd9p-6, 0x1.292834c025357p-6,
    },
    .{
        0x1.00abe0c129e1ep+0, 0x1.7ceb0ee49d42ap-57, 0x1.5a3385d5c7ba5p-3, 0x1.3bf51056f6637p-4,
        0x1.7dba76b124b37p-5, 0x1.07be4b02e94c4p-5,  0x1.8a6bb92513f01p-6, 0x1.36afd4c615aecp-6,
    },
    .{
        0x1.0102bcffd6acdp+0, -0x1.c2294c65d2e86p-55, 0x1.5caff17351901p-3, 0x1.407abbc04feb3p-4,
        0x1.86179b8005949p-5, 0x1.0f97520bd4e72p-5,   0x1.9950c5c89f3dfp-6, 0x1.44f2344e7b664p-6,
    },
    .{
        0x1.015a397cf0f1cp+0, -0x1.eebd6ccfe3ee3p-55, 0x1.5f3581be7b08bp-3, 0x1.4519ddf1ae531p-4,
        0x1.8eb4b6ee35e92p-5, 0x1.17bc85414cd46p-5,   0x1.a8e5895e3fcf9p-6, 0x1.53fafdc629400p-6,
    },
    .{
        0x1.01b2588811eebp+0, 0x1.7193e5d0a915fp-59, 0x1.61c46a67205d2p-3, 0x1.49d33a6eeae0bp-4,
        0x1.979438563c014p-5, 0x1.20316ae977f05p-5,  0x1.b9339afb53aa4p-6, 0x1.63d6b02c42d0ap-6,
    },
    .{
        0x1.020b1c7df0575p+0, -0x1.dd547e329c1e5p-55, 0x1.645ce0ab901bbp-3, 0x1.4ea79c34fc7a6p-4,
        0x1.a0b8ac08940ecp-5, 0x1.28f9babd0629bp-5,   0x1.ca452cf90a55ep-6, 0x1.7492b016730efp-6,
    },
    .{
        0x1.026487c8c5d71p+0, -0x1.5fd9b68dc3b6ep-54, 0x1.66ff1b67d5d70p-3, 0x1.5397d613373ebp-4,
        0x1.aa24bce3aeb4ap-5, 0x1.3219610b150acp-5,   0x1.dc251825103f1p-6, 0x1.863d5a3932532p-6,
    },
    .{
        0x1.02be9ce0b87cdp+0, 0x1.e5d09da2e0f04p-56, 0x1.69ab5325bc359p-3, 0x1.58a4c3097aab3p-4,
        0x1.b3db3605f46f2p-5, 0x1.3b94821742cabp-5,  0x1.eedee7da72a15p-6, 0x1.98e6179a3e9a0p-6,
    },
    .{
        0x1.03195e4c483f1p+0, -0x1.5db10ad66eacbp-54, 0x1.6c61c22d908f0p-3, 0x1.5dcf46ab9f2cbp-4,
        0x1.bddf049bb1f4dp-5, 0x1.456f7db6ac768p-5,   0x1.013f738bd7bb3p-5, 0x1.ac9d739783d21p-6,
    },
    .{
        0x1.0374cea0c0c9fp+0, -0x1.917bff5241c76p-54, 0x1.6f22a497b2ec0p-3, 0x1.63184d8a79db5p-4,
        0x1.c83339caf946ep-5, 0x1.4faef331019d4p-5,   0x1.0b8917547d678p-5, 0x1.c17533f147e1cp-6,
    },
    .{
        0x1.03d0f082afcc8p+0, -0x1.018bbcddb49ebp-54, 0x1.71ee385efdf06p-3, 0x1.6880cda2d3885p-4,
        0x1.d2db0cbfae54dp-5, 0x1.5a57c56b50c5ep-5,   0x1.16535a40098b2p-5, 0x1.d780730b8ebb8p-6,
    },
    .{
        0x1.042dc6a65ffbfp+0, -0x1.c7ea28dce95d1p-55, 0x1.74c4bd7412f9ep-3, 0x1.6e09c6d2b72bcp-4,
        0x1.ddd9dcda253dep-5, 0x1.656f1f62b5001p-5,   0x1.21a5ae2ac77eep-5, 0x1.eed3bca067f0ep-6,
    },
    .{
        0x1.048b53d05907bp+0, 0x1.634fffed6e2a6p-54, 0x1.77a675d1978bep-3, 0x1.73b4435583415p-4,
        0x1.e9333402ebbf3p-5, 0x1.70fa78fd9f73fp-5,  0x1.2d8804d934fe1p-5, 0x1.03c29691a281cp-5,
    },
    .{
        0x1.04e99ad5e4bcdp+0, -0x1.e97a72fe827e0p-54, 0x1.7a93a5917200cp-3, 0x1.7981584731c05p-4,
        0x1.f4eac9268fae2p-5, 0x1.7cff9c3b19721p-5,   0x1.3a02d9c1e0145p-5, 0x1.10d64a0e56953p-5,
    },
    .{
        0x1.05489e9d99995p+0, 0x1.d177637ec6a2bp-55, 0x1.7d8c930314681p-3, 0x1.7f72262f532e4p-4,
        0x1.0082416e39013p-4, 0x1.8984aac80ddf4p-5,  0x1.471f3caf18eb8p-5, 0x1.1eb1cce6dd570p-5,
    },
    .{
        0x1.05a8621feb16bp+0, -0x1.e5b33b1407c5fp-56, 0x1.809186c2e57ddp-3, 0x1.8587d99442dc8p-4,
        0x1.06c23d1dfcb7fp-4, 0x1.969024036dd22p-5,   0x1.54e6dd4d2af33p-5, 0x1.2d62f439f2a31p-5,
    },
    .{
        0x1.0608e867bff30p+0, 0x1.cbef5d8580027p-55, 0x1.83a2cbd2d8ba2p-3, 0x1.8bc3ab9724c6ep-4,
        0x1.0d377ef1e0c39p-4, 0x1.a428eb7addf84p-5,  0x1.636417bc01ff2p-5, 0x1.3cf8acc7eb2a0p-5,
    },
    .{
        0x1.066a34930ec8dp+0, -0x1.480f445fedad1p-54, 0x1.86c0afb447a74p-3, 0x1.9226e29948d9cp-4,
        0x1.13e44a9b5a3a6p-4, 0x1.b2564fea8b3fep-5,   0x1.72a2023d92458p-5, 0x1.4d8313cec3485p-5,
    },
    .{
        0x1.06cc49d38146cp+0, -0x1.b55394f4fc07bp-55, 0x1.89eb82831feedp-3, 0x1.98b2d2eb9bb23p-4,
        0x1.1acb01e9ab414p-4, 0x1.c12012cbd00c6p-5,   0x1.82ac7c1d15c38p-5, 0x1.5f13925c6edcap-5,
    },
    .{
        0x1.072f2b6f1e601p+0, -0x1.2dcbb05419970p-54, 0x1.8d2397127aebbp-3, 0x1.9f68df88da51dp-4,
        0x1.21ee26a4f62a1p-4, 0x1.d08e707f7ae6fp-5,   0x1.93903dee3feb0p-5, 0x1.71bcfb5c57b59p-5,
    },
    .{
        0x1.0792dcc0fbd20p+0, -0x1.5bf23ee4f9d54p-56, 0x1.9069430ab508ap-3, 0x1.a64a7adb4cd85p-4,
        0x1.29505c8b48349p-4, 0x1.e0aa2921cfa60p-5,   0x1.a55aeb46f4322p-5, 0x1.8593acad3becep-5,
    },
    .{
        0x1.07f76139f761dp+0, 0x1.fa1046481bb82p-54, 0x1.93bcdf091cca6p-3, 0x1.ad59278edc42fp-4,
        0x1.30f46b7261652p-4, 0x1.f17c8a17c843ep-5,  0x1.b81b2619e15b5p-5, 0x1.9aadb395f5ae4p-5,
    },
    .{
        0x1.085cbc61783c1p+0, 0x1.0a6e9efa20176p-54, 0x1.971ec6c1531e4p-3, 0x1.b496797068912p-4,
        0x1.38dd419140184p-4, 0x1.0187bc3357bbbp-4,  0x1.cbe0a3dcafe26p-5, 0x1.b122f4fa499d0p-5,
    },
    .{
        0x1.08c2f1d638e4cp+0, 0x1.b47c159534a3dp-56, 0x1.9a8f592078624p-3, 0x1.bc04165b57ab2p-4,
        0x1.410df5f4bed1dp-4, 0x1.0ab6bdf478c71p-4,  0x1.e0bc44a945c64p-5, 0x1.c90d59bcd5701p-5,
    },
    .{
        0x1.092a054f1a2fcp+0, -0x1.2f657224e9830p-54, 0x1.9e0ef87243a2cp-3, 0x1.c3a3b7366a278p-4,
        0x1.4989cb22e2175p-4, 0x1.1450e5ba7ad39p-4,   0x1.f6c02c8f0ef93p-5, 0x1.e288ffc8d182cp-5,
    },
    .{
        0x1.0991fa9bffbf4p+0, -0x1.ca1140a1abbf4p-58, 0x1.a19e0a8823b80p-3, 0x1.cb772900f9c24p-4,
        0x1.525431f0cbb2ep-4, 0x1.1e5c2d06804e1p-4,   0x1.06ffefa7aa6b8p-4, 0x1.fdb4704dca347p-5,
    },
    .{
        0x1.09fad5a6b68f9p+0, 0x1.aa1f06e92964ep-56, 0x1.a53cf8e28c50ep-3, 0x1.d3804df1de350p-4,
        0x1.5b70cc8fa98dcp-4, 0x1.28def298c979bp-4,  0x1.13482f6347eebp-4, 0x1.0d586de48358cp-4,
    },
    .{
        0x1.0a649a73e61f2p+0, 0x1.74ac0d817e9c7p-55, 0x1.a8ec30dc93891p-3, 0x1.dbc11ea950625p-4,
        0x1.64e371d5616d3p-4, 0x1.33e0023936249p-4,  0x1.204426263066ap-4, 0x1.1cd12e4629723p-4,
    },
    .{
        0x1.0acf4d240ccc4p+0, 0x1.da890f3b40bd3p-54, 0x1.acac23da07797p-3, 0x1.e43bab7741a98p-4,
        0x1.6eb030c631819p-4, 0x1.3f669d2eb516ep-4,  0x1.2e0006ae505aep-4, 0x1.2d58204457c82p-4,
    },
    .{
        0x1.0b3af1f4880bbp+0, 0x1.f450fb78d32bap-56, 0x1.b07d4778263afp-3, 0x1.ecf21db7be0efp-4,
        0x1.78db5465013e4p-4, 0x1.4b7a8376f0996p-4,  0x1.3c88f9f2ef221p-4, 0x1.3f02ad9eb9753p-4,
    },
    .{
        0x1.0ba78d40a9260p+0, -0x1.57b07a441e242p-54, 0x1.b46015c126262p-3, 0x1.f5e6b94713f3dp-4,
        0x1.836967d0afecfp-4, 0x1.5823fdd1707b9p-4,   0x1.4bed355269dc2p-4, 0x1.51e83065121cfp-4,
    },
    .{
        0x1.0c152382d7366p+0, -0x1.ee6913347c2a6p-54, 0x1.b8550d62bfb6ep-3, 0x1.ff1bde0fa3cadp-4,
        0x1.8e5f3ab550989p-4, 0x1.656be8b38ebafp-4,   0x1.5c3c13008a099p-4, 0x1.662225a1b4f77p-4,
    },
};

fn acosRefine(value: f64, approximation: f64) f64 {
    @setFloatMode(.strict);
    const square = value * value;
    const square_error = fma_impl.fma64(value, value, -square);
    var complement_low: f64 = undefined;
    var complement_high = fastTwoSum(1.0, -square, &complement_low);
    complement_low -= square_error;
    complement_high = fastTwoSum(complement_high, complement_low, &complement_low);
    const cosine_high = @sqrt(complement_high);
    const cosine_low = (complement_low - fma_impl.fma64(cosine_high, cosine_high, -complement_high)) *
        (0.5 / cosine_high);
    const table_index_signed = roundEvenFinite(
        absD(approximation - 0x1.921fb54442d18p+0) * 0x1.45f306dc9c883p+4,
    );
    const table_index: usize = @intCast(table_index_signed);
    const sine_table = [33][2]f64{
        .{ 0x0p+0, 0x0p+0 },
        .{ -0x1.912bd0d569a9p-61, 0x1.91f65f10dd814p-5 },
        .{ -0x1.e2718d26ed688p-60, 0x1.917a6bc29b42cp-4 },
        .{ 0x1.13000a89a11ep-58, 0x1.2c8106e8e613ap-3 },
        .{ -0x1.26d19b9ff8d82p-57, 0x1.8f8b83c69a60bp-3 },
        .{ -0x1.42deef11da2c4p-57, 0x1.f19f97b215f1bp-3 },
        .{ -0x1.5d28da2c4612dp-56, 0x1.294062ed59f06p-2 },
        .{ -0x1.efdc0d58cf62p-62, 0x1.58f9a75ab1fddp-2 },
        .{ -0x1.72cedd3d5a61p-57, 0x1.87de2a6aea963p-2 },
        .{ 0x1.5b362cb974183p-57, 0x1.b5d1009e15ccp-2 },
        .{ 0x1.e0d891d3c6841p-58, 0x1.e2b5d3806f63bp-2 },
        .{ -0x1.a5a014347406cp-55, 0x1.073879922ffeep-1 },
        .{ 0x1.b25dd267f66p-55, 0x1.1c73b39ae68c8p-1 },
        .{ -0x1.efcc626f74a6fp-57, 0x1.30ff7fce17035p-1 },
        .{ 0x1.8076a2cfdc6b3p-57, 0x1.44cf325091dd6p-1 },
        .{ -0x1.75720992bfbb2p-55, 0x1.57d69348cecap-1 },
        .{ -0x1.bdd3413b26456p-55, 0x1.6a09e667f3bcdp-1 },
        .{ -0x1.0f537acdf0ad7p-56, 0x1.7b5df226aafafp-1 },
        .{ -0x1.2c5e12ed1336dp-55, 0x1.8bc806b151741p-1 },
        .{ -0x1.30ee286712474p-55, 0x1.9b3e047f38741p-1 },
        .{ 0x1.9f630e8b6dac8p-60, 0x1.a9b66290ea1a3p-1 },
        .{ -0x1.bc69f324e6d61p-55, 0x1.b728345196e3ep-1 },
        .{ -0x1.6e0b1757c8d07p-56, 0x1.c38b2f180bdb1p-1 },
        .{ -0x1.e7b6bb5ab58aep-58, 0x1.ced7af43cc773p-1 },
        .{ 0x1.457e610231ac2p-56, 0x1.d906bcf328d46p-1 },
        .{ -0x1.014c76c126527p-55, 0x1.e212104f686e5p-1 },
        .{ 0x1.760b1e2e3f81ep-55, 0x1.e9f4156c62ddap-1 },
        .{ 0x1.52c7adc6b4989p-56, 0x1.f0a7efb9230d7p-1 },
        .{ 0x1.562172a361fd3p-56, 0x1.f6297cff75cbp-1 },
        .{ -0x1.7a0a8ca13571fp-55, 0x1.fa7557f08a517p-1 },
        .{ -0x1.87df6378811c7p-55, 0x1.fd88da3d12526p-1 },
        .{ -0x1.c57bc2e24aa15p-57, 0x1.ff621e3796d7ep-1 },
        .{ 0x0p+0, 0x1p+0 },
    };
    const cosine_index = 32 - table_index;
    const cosine_table_high = sine_table[cosine_index][1];
    const cosine_table_low = sine_table[cosine_index][0];
    const sine_table_high = sine_table[table_index][1];
    const sine_table_low = sine_table[table_index][0];
    const absolute_value = absD(value);
    const sine_delta_high = absolute_value - sine_table_high;
    const sine_delta_low = -sine_table_low;
    const cosine_delta_high = cosine_high - cosine_table_high;
    const cosine_delta_low = cosine_low - cosine_table_low;
    const magic = 0x1.8p-4;
    const sine_correction = fma_impl.fma64(sine_table_high, cosine_delta_high, magic) - magic;
    const sine_correction_error = fma_impl.fma64(sine_table_high, cosine_delta_high, -sine_correction);
    const cosine_correction = fma_impl.fma64(cosine_table_high, sine_delta_high, magic) - magic;
    const cosine_correction_error = fma_impl.fma64(cosine_table_high, sine_delta_high, -cosine_correction);
    var delta = cosine_correction - sine_correction;
    var delta_error = (cosine_table_high * sine_delta_low + cosine_table_low * sine_delta_high) -
        (sine_table_high * cosine_delta_low + sine_table_low * cosine_delta_high) -
        (sine_correction_error - cosine_correction_error);
    delta = fastTwoSum(delta, delta_error, &delta_error);
    const sign: f64 = if (value < 0.0) -1.0 else 1.0;
    const table_index_float: f64 = @floatFromInt(table_index);
    const table_index_signed_float = table_index_float * sign;
    const coefficients = [5][2]f64{
        .{ 0x1p+0, -0x1.fc2c76456515bp-108 },
        .{ 0x1.5555555555555p-3, 0x1.5555555623513p-57 },
        .{ 0x1.3333333333333p-4, 0x1.9997e3427441bp-59 },
        .{ 0x1.6db6db6db6db7p-5, -0x1.cb95ff08658e6p-62 },
        .{ 0x1.f1c71c71c6d5bp-6, 0x1.b125bccdcc89ep-60 },
    };
    const correction_coefficients = [3]f64{
        0x1.6e8ba2ec8cb69p-6,
        0x1.1c4ea7a15c997p-6,
        0x1.ca8355d39bb67p-7,
    };
    var delta_square_low: f64 = undefined;
    const delta_square_high = mulDd(delta, delta_error, delta, delta_error, &delta_square_low);
    delta *= -sign;
    delta_error *= -sign;
    var polynomial_low = delta_square_high *
        (correction_coefficients[0] + delta_square_high *
            (correction_coefficients[1] + delta_square_high * correction_coefficients[2]));
    var polynomial_high = polyDd(
        delta_square_high,
        delta_square_low,
        coefficients.len,
        &coefficients,
        &polynomial_low,
    );
    polynomial_high = mulDd(delta, delta_error, polynomial_high, polynomial_low, &polynomial_low);
    const jt = 32.0 - table_index_signed_float;
    var phase_high = jt * 0x1.921fb54442dp-5;
    var phase_low = 0x1.8469898cc518p-53 * jt;
    var phase_tail = -0x1.fc8f8cbb5bf6cp-102 * jt;
    phase_low = sum(polynomial_high, polynomial_low, phase_low, phase_tail, &phase_tail);
    phase_high = fastTwoSum(phase_high, phase_low, &phase_low);
    phase_low = fastTwoSum(phase_low, phase_tail, &phase_tail);
    phase_high = fastTwoSum(phase_high, phase_low, &phase_low);
    phase_low = fastTwoSum(phase_low, phase_tail, &phase_tail);
    var phase_low_bits: u64 = @bitCast(phase_low);
    const phase_low_exponent: i64 = @as(i64, @intCast((phase_low_bits >> 52) & 0x7ff)) - 1023;
    var round_exponent = 52 - (107 + phase_low_exponent);
    if (round_exponent < 0) round_exponent = 0;
    if (round_exponent > 52) round_exponent = 52;
    const low_mask = (@as(u64, 1) << 52) - (@as(u64, 1) << @intCast(round_exponent));
    const round_bit: u64 = if (round_exponent == 0)
        (@as(u64, 1) << 63)
    else
        (@as(u64, 1) << @intCast(round_exponent - 1));
    if (((phase_low_bits +% round_bit) & low_mask) == 0) {
        if (value == 0x1.ffffffffffdc0p-1) return 0x1.8000000000024p-22 + 0x1p-76;
        if (value == 0x1.53ea6c7255e88p-4) return 0x1.7cdacb6bbe707p+0 + 0x1p-54;
        if (value == 0x1.fd737be914578p-11) return 0x1.91e006d41d8d8p+0 + 0x1.8p-53;
        if (value == 0x1.fffffffffff70p-1) return 0x1.8000000000009p-23 + 0x1p-77;
        if (value == 0x1.390e6939cd1a6p-5) return 0x1.8856a5d3296a4p+0 - 0x1p-109;
        const phase_tail_bits: u64 = @bitCast(phase_tail);
        if (((phase_tail_bits ^ phase_low_bits) >> 63) != 0) {
            phase_low_bits -%= 1;
        } else {
            phase_low_bits +%= 1;
        }
        phase_low = @bitCast(phase_low_bits);
    }
    return phase_high + phase_low;
}

pub fn acosd(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask_d;
    if (magnitude > infinity_bits_d) return @bitCast(bits | quiet_nan_bit_d);
    if (magnitude >= infinity_bits_d or magnitude > 0x3ff0_0000_0000_0000) {
        return @bitCast(canonical_nan_bits_d);
    }
    if (magnitude == 0x3ff0_0000_0000_0000) {
        if ((bits & sign_mask_d) == 0) return @bitCast(@as(u64, 0));
        return @bitCast(@as(u64, 0x4009_21fb_5444_2d18));
    }
    if (magnitude == 0) return @bitCast(@as(u64, 0x3ff9_21fb_5444_2d18));

    const doubled_magnitude = magnitude << 1;
    var offset_high: f64 = undefined;
    var offset_low: f64 = undefined;
    var reduced: f64 = undefined;
    var root: f64 = undefined;
    var root_low: f64 = undefined;
    var error_bound: f64 = undefined;
    var coefficient_index: usize = undefined;
    if (doubled_magnitude > 0x7fc0_0000_0000_0000) {
        const negative = (bits & sign_mask_d) != 0;
        offset_high = if (negative) 0x1.921fb54442d18p+1 else 0.0;
        offset_low = if (negative) 0x1.1a62633145c07p-53 else 0.0;
        const absolute_value = absD(value);
        reduced = 2.0 - 2.0 * absolute_value;
        coefficient_index = @intCast(roundEvenFinite(reduced * 0x1p5));
        const table_offset: f64 = @floatFromInt(coefficient_index);
        root = copySignD(@sqrt(reduced), value);
        root_low = fma_impl.fma64(root, root, -reduced) * ((-0.5 / reduced) * root);
        reduced = 0.25 * reduced - table_offset * 0x1p-7;
        error_bound = absD(root * reduced) * 0x1.8cp-52 + 0x1p-105;
    } else {
        offset_high = 0x1.921fb54442d18p+0;
        offset_low = 0x1.1a62633145c07p-54;
        if (doubled_magnitude <= 0x7e00_0000_0000_0000) {
            const value_correction: f64 = if (doubled_magnitude <= 0x7919_6767_0d2e_8fe8)
                0.0
            else
                (value * value) * (-0x1.5555555555555p-3 * value);
            var sum_low: f64 = undefined;
            const sum_high = fastTwoSum(offset_high, -value, &sum_low);
            const low = value_correction + (sum_low + offset_low);
            const lower = sum_high + (low - 0x1.34p-79);
            const upper = sum_high + (low + 0x1.34p-79);
            if (lower != upper) return acosRefine(value, lower);
            return lower;
        }
        const square = value * value;
        coefficient_index = @intCast(roundEvenFinite(square * 0x1p7));
        const table_offset: f64 = @floatFromInt(coefficient_index);
        reduced = fma_impl.fma64(value, value, -0x1p-7 * table_offset);
        root = -value;
        root_low = 0.0;
        error_bound = absD(value * reduced) * 0x1.81p-52;
    }
    const coefficient = acos_coefficients[coefficient_index];
    const reduced_squared = reduced * reduced;
    const coefficient_2 = fma_impl.fma64(reduced, coefficient[3], coefficient[2]);
    const coefficient_4 = fma_impl.fma64(reduced, coefficient[5], coefficient[4]);
    const coefficient_6 = fma_impl.fma64(reduced, coefficient[7], coefficient[6]);
    const tail = fma_impl.fma64(reduced_squared, coefficient_6, coefficient_4);
    const polynomial = reduced * fma_impl.fma64(reduced_squared, tail, coefficient_2);
    var high = coefficient[0];
    var low = coefficient[1] + polynomial;
    high = mulDd(root, root_low, high, low, &low);
    high = fastSum(offset_high, offset_low, high, low, &low);
    const lower = high + (low - error_bound);
    const upper = high + (low + error_bound);
    if (lower != upper) return acosRefine(value, lower);
    return lower;
}
