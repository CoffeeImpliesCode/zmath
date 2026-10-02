// Correctly rounded hyperbolic cosine for binary16, binary32, and binary64 values.
//
// The algorithms are derived from CORE-MATH at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and are used under the MIT License.
// Binary16 copyright (c) 2025 Paul Zimmermann.
// Binary32 copyright (c) 2022-2023 Alexei Sibidanov.
// Binary64 copyright (c) 2023-2025 Alexei Sibidanov.

const fma_impl = @import("scalar_fma.zig");

const coshh_data = @import("scalar_coshh_data.zig");

const coshh_magnitude_mask: u16 = 0x7fff;
const coshh_infinity_bits: u16 = 0x7c00;
const coshh_quiet_nan_bit: u16 = 0x0200;

/// Correctly rounded binary16 hyperbolic cosine under caller round-to-nearest-even.
pub fn coshh(value: f16) f16 {
    @setFloatMode(.strict);
    const bits: u16 = @bitCast(value);
    const magnitude = bits & coshh_magnitude_mask;
    if (magnitude > coshh_infinity_bits) return @bitCast(bits | coshh_quiet_nan_bit);
    if (magnitude == coshh_infinity_bits) return @bitCast(coshh_infinity_bits);
    if (magnitude == 0) return 1.0;

    const table1_index: usize = @intCast(magnitude >> 5);
    const table2_index: usize = @intCast(((magnitude >> 10) << 5) | (magnitude & 0x1f));
    const sinh_x1: f32 = @bitCast(coshh_data.s1[table1_index]);
    const cosh_x1: f32 = @bitCast(coshh_data.c1[table1_index]);
    const sinh_x2: f32 = @bitCast(coshh_data.s2[table2_index]);
    const cosh_x2: f32 = @bitCast(coshh_data.c2[table2_index]);
    const result = fma_impl.fma(cosh_x1, cosh_x2, sinh_x1 * sinh_x2);
    return @floatCast(result);
}

const coshf_infinity_bits: u32 = 0x7f80_0000;
const coshf_magnitude_mask: u32 = 0x7fff_ffff;
const coshf_quiet_nan_bit: u32 = 0x0040_0000;

const coshf_tb: [32]u64 = .{
    0x3fe0000000000000, 0x3fe059b0d3158574, 0x3fe0b5586cf9890f, 0x3fe11301d0125b51,
    0x3fe172b83c7d517b, 0x3fe1d4873168b9aa, 0x3fe2387a6e756238, 0x3fe29e9df51fdee1,
    0x3fe306fe0a31b715, 0x3fe371a7373aa9cb, 0x3fe3dea64c123422, 0x3fe44e086061892d,
    0x3fe4bfdad5362a27, 0x3fe5342b569d4f82, 0x3fe5ab07dd485429, 0x3fe6247eb03a5585,
    0x3fe6a09e667f3bcd, 0x3fe71f75e8ec5f74, 0x3fe7a11473eb0187, 0x3fe82589994cce13,
    0x3fe8ace5422aa0db, 0x3fe93737b0cdc5e5, 0x3fe9c49182a3f090, 0x3fea5503b23e255d,
    0x3feae89f995ad3ad, 0x3feb7f76f2fb5e47, 0x3fec199bdd85529c, 0x3fecb720dcef9069,
    0x3fed5818dcfba487, 0x3fedfc97337b9b5f, 0x3feea4afa2a490da, 0x3fef50765b6e4540,
};

fn coshfTableValue(index: i64) f64 {
    const index_bits: u64 = @bitCast(index);
    const exponent_bits: u64 = @bitCast(index >> 5);
    const bits = coshf_tb[@intCast(index_bits & 31)] +% (exponent_bits << 52);
    return @bitCast(bits);
}

/// Correctly rounded binary32 hyperbolic cosine under caller round-to-nearest-even.
pub fn coshf(value: f32) f32 {
    @setFloatMode(.strict);
    const bits: u32 = @bitCast(value);
    const magnitude = bits & coshf_magnitude_mask;
    if (magnitude > coshf_infinity_bits) return @bitCast(bits | coshf_quiet_nan_bit);
    if (magnitude == coshf_infinity_bits) return @bitCast(coshf_infinity_bits);

    const absolute_shifted = bits << 1;
    if (absolute_shifted > 0x8565_a9f8) return @bitCast(coshf_infinity_bits);
    if (absolute_shifted < 0x7c00_0000) {
        if (absolute_shifted < 0x7400_0000) {
            if (absolute_shifted < 0x6600_0000) {
                // The correction is below half an ulp of one in this branch.
                return 1.0;
            }
            return (0.5 * value) * value + 1.0;
        }
        const coefficients = [_]f64{
            0x1.fffffffffffe3p-2,
            0x1.55555555723cfp-5,
            0x1.6c16bee4a5986p-10,
            0x1.a0483fc0328f7p-16,
        };
        const z: f64 = @floatCast(value);
        const z_squared = z * z;
        const z_fourth = z_squared * z_squared;
        return @floatCast(1.0 + z_squared *
            ((coefficients[0] + z_squared * coefficients[1]) +
                z_fourth * (coefficients[2] + z_squared * coefficients[3])));
    }

    const z: f64 = @floatCast(value);
    const iln2 = 0x1.71547652b82fep+5;
    const a = iln2 * z;
    // The finite path bounds a to less than 2^23, so the bias has unit spacing.
    const rounded_bits: u64 = @bitCast(a + 0x1.8p+52);
    const ia: i64 = @bitCast(rounded_bits -% 0x4338_0000_0000_0000);
    const ia_float: f64 = @floatFromInt(ia);
    const h = a - ia_float;
    const h_squared = h * h;
    const jp: i64 = @bitCast(rounded_bits);
    const jm = -jp;
    const sp = coshfTableValue(jp);
    const sm = coshfTableValue(jm);
    const coefficients = [_]f64{
        1.0,
        0x1.62e42fef4c4e7p-6,
        0x1.ebfd1b232f475p-13,
        0x1.c6b19384ecd93p-20,
    };
    const even = coefficients[0] + h_squared * coefficients[2];
    const odd = coefficients[1] + h_squared * coefficients[3];
    var result = sp * (even + h * odd) + sm * (even - h * odd);
    var upper: f32 = @floatCast(result);
    const lower: f32 = @floatCast(result - 1.45e-10 * result);
    if (upper != lower) {
        const accurate_coefficients = [_]f64{
            1.0,
            0x1.62e42fefa39efp-6,
            0x1.ebfbdff82c58fp-13,
            0x1.c6b08d702e0edp-20,
            0x1.3b2ab6fb92e5ep-27,
            0x1.5d886e6d54203p-35,
            0x1.430976b8ce6efp-43,
        };
        const iln2_high = 0x1.7154765p+5;
        const iln2_low = 0x1.5c17f0bbbe88p-26;
        const accurate_h = (iln2_high * z - ia_float) + iln2_low * z;
        const accurate_h_squared = accurate_h * accurate_h;
        const accurate_even = accurate_coefficients[0] + accurate_h_squared * accurate_coefficients[2] +
            (accurate_h_squared * accurate_h_squared) *
                (accurate_coefficients[4] + accurate_h_squared * accurate_coefficients[6]);
        const accurate_odd = accurate_coefficients[1] + accurate_h_squared *
            (accurate_coefficients[3] + accurate_h_squared * accurate_coefficients[5]);
        result = sp * (accurate_even + accurate_h * accurate_odd) +
            sm * (accurate_even - accurate_h * accurate_odd);
        upper = @floatCast(result);
    }
    return upper;
}

const coshd_magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const coshd_infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const coshd_quiet_nan_bit: u64 = 0x0008_0000_0000_0000;

const CoshdDoubleDouble = struct {
    high: f64,
    low: f64,
};

const coshd_t0: [64][2]f64 = .{
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

const coshd_t1: [64][2]f64 = .{
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

const coshd_database: [21][3]f64 = .{
    .{ 0x1.9a5e3cbe1985ep-4, 0x1.01492f72f984bp+0, -0x1p-107 },
    .{ 0x1.52a11832e847dp-3, 0x1.0381e68cac923p+0, 0x1p-104 },
    .{ 0x1.bf0305e2c6c37p-3, 0x1.061f4c39e16f2p+0, 0x1p-107 },
    .{ 0x1.17326ffc09f68p-2, 0x1.099318a43ac8p+0, 0x1p-104 },
    .{ 0x1.3d27bf16d8bdbp-2, 0x1.0c6091056e06ap+0, -0x1p-107 },
    .{ 0x1.03923f2b47c07p-1, 0x1.219c1989e3373p+0, -0x1p-54 },
    .{ 0x1.a6031cd5f93bap-1, 0x1.5bff041b260fep+0, -0x1p-107 },
    .{ 0x1.104b648f113a1p+0, 0x1.9efdca62b700ap+0, -0x1p-109 },
    .{ 0x1.1585720f35cd9p+0, 0x1.a5bf3acfde4b2p+0, 0x1p-105 },
    .{ 0x1.e9cc7ed2e1a7ep+0, 0x1.bb0ff220d8eb5p+1, -0x1p-53 },
    .{ 0x1.43180ea854696p+1, 0x1.91f1122b6b63ap+2, 0x1p-102 },
    .{ 0x1.725811dcf6782p+2, 0x1.45ea160ddc71fp+7, -0x1p-100 },
    .{ 0x1.5afd56f7d565bp+3, 0x1.8ff8e0ccea7cp+14, 0x1p-90 },
    .{ 0x1.759a2ad4c4d56p+3, 0x1.cb62eec26bd78p+15, -0x1p-92 },
    .{ 0x1.7fce95ea5c653p+3, 0x1.3bf8009648dcp+16, 0x1p-88 },
    .{ 0x1.743d5609348acp+4, 0x1.7a87a8bb7fa28p+32, -0x1p-22 },
    .{ 0x1.e07e71bfcf06fp+5, 0x1.91ec4412c344fp+85, 0x1p-24 },
    .{ 0x1.6474c604cc0d7p+6, 0x1.7a8f65ad009bdp+127, -0x1p+20 },
    .{ 0x1.54cd1fea7663ap+7, 0x1.c90810d354618p+244, 0x1p+135 },
    .{ 0x1.2da9e5e6af0bp+8, 0x1.27d6fe867d6f6p+434, 0x1p+329 },
    .{ 0x1.d6479eba7c971p+8, 0x1.62a88613629b6p+677, -0x1p+568 },
};

const coshd_zero_even = [4][2]f64{
    .{ 0x1p-1, -0x1.c7e8db669f624p-111 },
    .{ 0x1.5555555555555p-5, 0x1.5555555556135p-59 },
    .{ 0x1.6c16c16c16c17p-10, -0x1.f49f4a6e838f2p-65 },
    .{ 0x1.a01a01a01a01ap-16, 0x1.a4ffbe15316aap-76 },
};
const coshd_zero_odd = [4]f64{
    0x1.27e4fb7789f5cp-22,
    0x1.1eed8eff9089cp-29,
    0x1.939749ce13dadp-37,
    0x1.ae9891efb6691p-45,
};

fn coshdFastTwoSum(x: f64, y: f64) CoshdDoubleDouble {
    @setFloatMode(.strict);
    const sum = x + y;
    const rounded_x = sum - x;
    return .{ .high = sum, .low = y - rounded_x };
}

fn coshdMuldd(x_high: f64, x_low: f64, c_high: f64, c_low: f64) CoshdDoubleDouble {
    @setFloatMode(.strict);
    const high = c_high * x_high;
    const low = fma_impl.fma64(c_high, x_high, -high) + x_high * c_low + c_high * x_low;
    return .{ .high = high, .low = low };
}

fn coshdPolydd(x_high: f64, x_low: f64, coefficients: []const [2]f64, low: *f64) f64 {
    @setFloatMode(.strict);
    var index: usize = coefficients.len - 1;
    var sum = coshdFastTwoSum(coefficients[index][0], low.*);
    var high = sum.high;
    var coefficient_low = sum.low + coefficients[index][1];
    while (index > 0) : (index -= 1) {
        const product = coshdMuldd(x_high, x_low, high, coefficient_low);
        sum = coshdFastTwoSum(coefficients[index - 1][0], product.high);
        high = sum.high;
        coefficient_low = (product.low + coefficients[index - 1][1]) + sum.low;
    }
    low.* = coefficient_low;
    return high;
}

fn coshdExpAccurate(x: f64, t: f64, t_high: f64, t_low: f64, low: *f64) f64 {
    @setFloatMode(.strict);
    const coefficients = [3][2]f64{
        .{ 0x1p+0, 0x1.6c16bd194535dp-94 },
        .{ 0x1p-1, -0x1.8259d904fd34fp-93 },
        .{ 0x1.5555555555555p-3, 0x1.53e93e9f26e62p-57 },
    };
    const log_two_high = 0x1.62e42ffp-13;
    const log_two_low = 0x1.718432a1b0e26p-47;
    const log_two_tail = 0x1.9ff0342542fc3p-102;
    const delta = x - log_two_high * t;
    var delta_low = log_two_low * t;
    const delta_tail = log_two_tail * t + fma_impl.fma64(log_two_low, t, -delta_low);
    const delta_high = delta + delta_low;
    delta_low = ((delta - delta_high) + delta_low) + delta_tail;
    var polynomial_low = delta_high *
        (0x1.5555555555555p-5 + delta_high *
            (0x1.11111113e93e9p-7 + delta_high * 0x1.6c16c169400a7p-10));
    var polynomial_high = coshdPolydd(delta_high, delta_low, &coefficients, &polynomial_low);
    var product = coshdMuldd(delta_high, delta_low, polynomial_high, polynomial_low);
    polynomial_high = product.high;
    polynomial_low = product.low;
    product = coshdMuldd(t_high, t_low, polynomial_high, polynomial_low);
    const z_high = t_high + product.high;
    const z_low = (t_high - z_high) + product.high;
    const u_high = z_high + t_low;
    const u_low = ((z_high - u_high) + t_low) + z_low;
    const v_high = u_high + product.low;
    low.* = ((u_high - v_high) + product.low) + u_low;
    return v_high;
}

noinline fn coshdDatabase(value: f64, fallback: f64) f64 {
    @setFloatMode(.strict);
    var lower: isize = 0;
    var upper: isize = @intCast(coshd_database.len - 1);
    var middle: isize = @divTrunc(lower + upper, 2);
    const absolute_value = if (value < 0.0) -value else value;
    var iteration: u8 = 0;
    while (iteration < 8 and lower <= upper) : (iteration += 1) {
        const row = coshd_database[@intCast(middle)];
        if (row[0] < absolute_value) {
            lower = middle + 1;
        } else if (row[0] == absolute_value) {
            return row[1] + row[2];
        } else {
            upper = middle - 1;
        }
        middle = @divTrunc(lower + upper, 2);
    }
    return fallback;
}

noinline fn coshdZero(value: f64) f64 {
    @setFloatMode(.strict);
    const square_high = value * value;
    const square_low = fma_impl.fma64(value, value, -square_high);
    var polynomial_low = square_high *
        (coshd_zero_odd[0] + square_high *
            (coshd_zero_odd[1] + square_high *
                (coshd_zero_odd[2] + square_high * coshd_zero_odd[3])));
    var polynomial_high = coshdPolydd(square_high, square_low, coshd_zero_even[0..], &polynomial_low);
    const product = coshdMuldd(polynomial_high, polynomial_low, square_high, square_low);
    polynomial_high = product.high;
    var sum = coshdFastTwoSum(1.0, polynomial_high);
    const result_zero = sum.high;
    sum = coshdFastTwoSum(sum.low, product.low);
    const result_one = sum.high;
    const result_two = sum.low;
    var rounded = result_one;
    var rounded_bits: u64 = @bitCast(rounded);
    if ((rounded_bits & 0x000f_ffff_ffff_ffff) == 0) {
        const correction_bits: u64 = @bitCast(result_two);
        if (((correction_bits ^ rounded_bits) >> 63) != 0) {
            rounded_bits -%= 1;
        } else {
            rounded_bits +%= 1;
        }
        rounded = @bitCast(rounded_bits);
    }
    if ((rounded_bits & 0x000f_ffff_ffff_ffff) == 0x000f_ffff_ffff_ffff) {
        return coshdDatabase(value, result_zero + rounded);
    }
    return result_zero + rounded;
}

noinline fn coshdFinalizeAccurate(value: f64, high_result: f64, low_result: f64) f64 {
    @setFloatMode(.strict);
    const result = coshdFastTwoSum(high_result, low_result);
    const high_bits: u64 = @bitCast(result.high);
    const low_bits: u64 = @bitCast(result.low);
    const high_exponent: i64 = @intCast((high_bits >> 52) & 0x7ff);
    const low_exponent: i64 = @intCast((low_bits >> 52) & 0x7ff);
    const low_magnitude = (low_bits +% 8) & 0x000f_ffff_ffff_ffff;
    const rounded = result.high + result.low;
    if (low_magnitude <= 16 or high_exponent - low_exponent > 103) {
        return coshdDatabase(value, rounded);
    }
    return rounded;
}

noinline fn coshdLargeAccurate(
    value: f64,
    absolute_value: f64,
    reduction: f64,
    scaled_sp: f64,
    table_high_input: f64,
    table_low_input: f64,
) f64 {
    @setFloatMode(.strict);
    var table_low = table_low_input;
    var table_high = coshdExpAccurate(
        absolute_value,
        reduction,
        table_high_input,
        table_low,
        &table_low,
    );
    const normalized = coshdFastTwoSum(table_high, table_low);
    table_high = normalized.high;
    table_low = normalized.low;
    const high_bits: u64 = @bitCast(table_high);
    const low_bits: u64 = @bitCast(table_low);
    const high_exponent: i64 = @intCast((high_bits >> 52) & 0x7ff);
    const low_exponent: i64 = @intCast((low_bits >> 52) & 0x7ff);
    const low_magnitude = (low_bits +% 8) & 0x000f_ffff_ffff_ffff;
    table_high += table_low;
    table_high *= 2.0;
    table_high *= scaled_sp;
    if (low_magnitude <= 16 or high_exponent - low_exponent > 103) {
        return coshdDatabase(value, table_high);
    }
    return table_high;
}

noinline fn coshdMediumAccurate(
    value: f64,
    absolute_value: f64,
    reduction: f64,
    magnitude: u64,
    table_high_input: f64,
    table_low_input: f64,
    negative_high: f64,
    sm: f64,
    table_j0: usize,
    table_j1: usize,
) f64 {
    @setFloatMode(.strict);
    var table_low = table_low_input;
    const table_high = coshdExpAccurate(
        absolute_value,
        reduction,
        table_high_input,
        table_low,
        &table_low,
    );
    var high_result: f64 = undefined;
    var low_result: f64 = undefined;
    if (magnitude > 0x403f_6666_6666_6666) {
        const sum = coshdFastTwoSum(table_high, negative_high);
        high_result = sum.high;
        low_result = ((table_high - high_result) + negative_high) + table_low;
    } else {
        const q0_low = coshd_t0[table_j0][0];
        const q1_low = coshd_t1[table_j1][0];
        const q_product = coshdMuldd(
            coshd_t0[table_j0][1],
            q0_low,
            coshd_t1[table_j1][1],
            q1_low,
        );
        var accurate_negative_low = q_product.low * sm;
        const accurate_negative_high = coshdExpAccurate(
            -absolute_value,
            -reduction,
            q_product.high * sm,
            accurate_negative_low,
            &accurate_negative_low,
        );
        const sum = coshdFastTwoSum(table_high, accurate_negative_high);
        high_result = sum.high;
        low_result = ((table_high - high_result) + accurate_negative_high) +
            accurate_negative_low + table_low;
    }
    return coshdFinalizeAccurate(value, high_result, low_result);
}

noinline fn coshdSmallAccurate(
    value: f64,
    absolute_value: f64,
    reduction: f64,
    table_high_input: f64,
    table_low_input: f64,
    negative_high: f64,
    negative_low_input: f64,
) f64 {
    @setFloatMode(.strict);
    var table_low = table_low_input;
    var negative_low = negative_low_input;
    const table_high = coshdExpAccurate(
        absolute_value,
        reduction,
        table_high_input,
        table_low,
        &table_low,
    );
    const accurate_negative = coshdExpAccurate(
        -absolute_value,
        -reduction,
        negative_high,
        negative_low,
        &negative_low,
    );
    const sum = coshdFastTwoSum(table_high, accurate_negative);
    const high_result = sum.high;
    const low_result =
        ((table_high - high_result) + accurate_negative) + negative_low + table_low;
    return coshdFinalizeAccurate(value, high_result, low_result);
}

/// Correctly rounded binary64 hyperbolic cosine under caller round-to-nearest-even.
pub fn coshd(value: f64) f64 {
    @setFloatMode(.strict);
    const bits: u64 = @bitCast(value);
    const magnitude = bits & coshd_magnitude_mask;
    if (magnitude > 0x4086_33ce_8fb9_f87d) {
        if (magnitude > coshd_infinity_bits) return @bitCast(bits | coshd_quiet_nan_bit);
        return @bitCast(coshd_infinity_bits);
    }

    const absolute_value: f64 = @bitCast(magnitude);
    if (magnitude < 0x3fc0_0000_0000_0000) {
        if (magnitude < 0x3e50_0000_0000_0000) {
            // The correction is below half an ulp of one in this branch.
            return 1.0;
        }
        const coefficients = [5]f64{
            0x1p-1,
            0x1.5555555555554p-5,
            0x1.6c16c16c1d0cp-10,
            0x1.a01a0075066b4p-16,
            0x1.27faff8dcc1c8p-22,
        };
        const square = value * value;
        const fourth = square * square;
        const polynomial = square * ((coefficients[0] + square * coefficients[1]) +
            fourth * ((coefficients[2] + square * coefficients[3]) + fourth * coefficients[4]));
        const error_bound = square * 0x1.84p-51;
        const lower = 1.0 + (polynomial - error_bound);
        const upper = 1.0 + (polynomial + error_bound);
        if (lower == upper) {
            return lower;
        } else {
            @branchHint(.cold);
            return coshdZero(value);
        }
    }

    const scale = 0x1.71547652b82fep+12;
    const reduction_value = fma_impl.fma64(absolute_value, scale, 0x1.8000002p+26);
    const reduction_bits: u64 = @bitCast(reduction_value);
    const rounded_reduction_bits = reduction_bits & ~@as(u64, 0x03ff_ffff);
    const reduction = @as(f64, @bitCast(rounded_reduction_bits)) - 0x1.8p26;
    const il: i64 = @intCast((reduction_bits << 14) >> 40);
    const table_i1: usize = @intCast(@as(u64, @bitCast(il)) & 0x3f);
    const table_i0: usize = @intCast(@as(u64, @bitCast(il >> 6)) & 0x3f);
    const ie: i64 = il >> 12;
    const t0_high = coshd_t0[table_i0][1];
    const t0_low = coshd_t0[table_i0][0];
    const t1_high = coshd_t1[table_i1][1];
    const t1_low = coshd_t1[table_i1][0];
    const table_product = coshdMuldd(t0_high, t0_low, t1_high, t1_low);
    var table_high = table_product.high;
    var table_low = table_product.low;
    const log_two_high = 0x1.62e42ffp-13;
    const log_two_low = 0x1.718432a1b0e26p-47;
    const delta = (absolute_value - log_two_high * reduction) + log_two_low * reduction;
    const delta_square = delta * delta;
    const coefficients = [4]f64{ 1.0, 0.5, 0x1.5555555aaaaaep-3, 0x1.55555551c98cp-5 };
    const positive_polynomial = delta * ((coefficients[0] + delta * coefficients[1]) +
        delta_square * (coefficients[2] + delta * coefficients[3]));

    if (magnitude > 0x4042_5e4f_7b27_37fa) {
        const scaled_sp_bits = (@as(u64, @intCast(1021 + ie))) << 52;
        const scaled_sp: f64 = @bitCast(scaled_sp_bits);
        const error_bound = 0.12e-18 * table_high;
        const lower = table_high + (table_low + table_high * positive_polynomial - error_bound);
        const upper = table_high + (table_low + table_high * positive_polynomial + error_bound);
        if (lower == upper) {
            return (lower * scaled_sp) * 2.0;
        } else {
            @branchHint(.cold);
            return coshdLargeAccurate(
                value,
                absolute_value,
                reduction,
                scaled_sp,
                table_high,
                table_low,
            );
        }
    }

    const sp_bits = (@as(u64, @bitCast(@as(i64, 1022 + ie)))) << 52;
    const sp: f64 = @bitCast(sp_bits);
    const jl = -il;
    const table_j1: usize = @intCast(@as(u64, @bitCast(jl)) & 0x3f);
    const table_j0: usize = @intCast(@as(u64, @bitCast(jl >> 6)) & 0x3f);
    const je: i64 = jl >> 12;
    const sm_bits = (@as(u64, @bitCast(@as(i64, 1022 + je)))) << 52;
    const sm: f64 = @bitCast(sm_bits);
    const minus_delta = -delta;

    if (magnitude > 0x4014_0000_0000_0000) {
        const q0_high = coshd_t0[table_j0][1];
        const q1_high = coshd_t1[table_j1][1];
        const q_high = q0_high * q1_high;
        table_high *= sp;
        table_low *= sp;
        const negative_high = q_high * sm;
        const negative_polynomial = minus_delta *
            ((coefficients[0] + minus_delta * coefficients[1]) +
                delta_square * (coefficients[2] + minus_delta * coefficients[3]));
        const negative_exp = negative_high + negative_high * negative_polynomial;
        const high_result = table_high;
        const low_result = (table_low + negative_exp) + table_high * positive_polynomial;
        const error_bound = 0x1.202p-63 * high_result;
        const lower = high_result + (low_result - error_bound);
        const upper = high_result + (low_result + error_bound);
        if (lower == upper) {
            return lower;
        } else {
            @branchHint(.cold);
            return coshdMediumAccurate(
                value,
                absolute_value,
                reduction,
                magnitude,
                table_high,
                table_low,
                negative_high,
                sm,
                table_j0,
                table_j1,
            );
        }
    } else {
        const q0_high = coshd_t0[table_j0][1];
        const q0_low = coshd_t0[table_j0][0];
        const q1_high = coshd_t1[table_j1][1];
        const q1_low = coshd_t1[table_j1][0];
        const negative_product = coshdMuldd(q0_high, q0_low, q1_high, q1_low);
        table_high *= sp;
        table_low *= sp;
        const negative_high = negative_product.high * sm;
        const negative_low = negative_product.low * sm;
        const negative_polynomial = minus_delta *
            ((coefficients[0] + minus_delta * coefficients[1]) +
                delta_square * (coefficients[2] + minus_delta * coefficients[3]));
        const positive_high = table_high;
        const positive_low = table_low + table_high * positive_polynomial;
        const negative_sum_high = negative_high;
        const negative_sum_low = negative_low + negative_high * negative_polynomial;
        const high_result = positive_high + negative_sum_high;
        const low_result =
            ((positive_high - high_result) + negative_sum_high) + negative_sum_low + positive_low;
        const error_bound = 0x1.c0ap-62 * high_result;
        const lower = high_result + (low_result - error_bound);
        const upper = high_result + (low_result + error_bound);
        if (lower == upper) {
            return lower;
        } else {
            @branchHint(.cold);
            return coshdSmallAccurate(
                value,
                absolute_value,
                reduction,
                table_high,
                table_low,
                negative_high,
                negative_low,
            );
        }
    }
}
