// This file is part of the CORE-MATH project
// (https://core-math.gitlabpages.inria.fr/).
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
const sign_mask: u64 = 0x8000_0000_0000_0000;
const magnitude_mask: u64 = 0x7fff_ffff_ffff_ffff;
const infinity_bits: u64 = 0x7ff0_0000_0000_0000;
const quiet_nan_bit: u64 = 0x0008_0000_0000_0000;
const canonical_nan_bits: u64 = 0x7ff8_0000_0000_0000;
const minus_infinity_bits: u64 = 0xfff0_0000_0000_0000;

const B = [_]struct { c0: u16, c1: i16 }{
    .{ .c0 = 301, .c1 = 27565 },
    .{ .c0 = 7189, .c1 = 24786 },
    .{ .c0 = 13383, .c1 = 22167 },
    .{ .c0 = 18923, .c1 = 19696 },
    .{ .c0 = 23845, .c1 = 17361 },
    .{ .c0 = 28184, .c1 = 15150 },
    .{ .c0 = 31969, .c1 = 13054 },
    .{ .c0 = 35231, .c1 = 11064 },
    .{ .c0 = 37996, .c1 = 9173 },
    .{ .c0 = 40288, .c1 = 7372 },
    .{ .c0 = 42129, .c1 = 5657 },
    .{ .c0 = 43542, .c1 = 4020 },
    .{ .c0 = 44546, .c1 = 2457 },
    .{ .c0 = 45160, .c1 = 962 },
    .{ .c0 = 45399, .c1 = -468 },
    .{ .c0 = 45281, .c1 = -1838 },
    .{ .c0 = 44821, .c1 = -3151 },
    .{ .c0 = 44032, .c1 = -4412 },
    .{ .c0 = 42929, .c1 = -5622 },
    .{ .c0 = 41522, .c1 = -6786 },
    .{ .c0 = 39825, .c1 = -7905 },
    .{ .c0 = 37848, .c1 = -8982 },
    .{ .c0 = 35602, .c1 = -10020 },
    .{ .c0 = 33097, .c1 = -11020 },
    .{ .c0 = 30341, .c1 = -11985 },
    .{ .c0 = 27345, .c1 = -12916 },
    .{ .c0 = 24115, .c1 = -13816 },
    .{ .c0 = 20661, .c1 = -14685 },
    .{ .c0 = 16989, .c1 = -15526 },
    .{ .c0 = 13107, .c1 = -16339 },
    .{ .c0 = 9022, .c1 = -17126 },
    .{ .c0 = 4740, .c1 = -17889 },
};

const r1 = [_]f64{
    0x1.7154800000000p+0,
    0x1.696af49200000p+0,
    0x1.61ab3fb680000p+0,
    0x1.5a18441680000p+0,
    0x1.52adadb480000p+0,
    0x1.4b6b7c9080000p+0,
    0x1.4451b0aa80000p+0,
    0x1.3d5ed8ae00000p+0,
    0x1.3691834680000p+0,
    0x1.2fe9b07400000p+0,
    0x1.2967603680000p+0,
    0x1.2307afe500000p+0,
    0x1.1cca9f7f80000p+0,
    0x1.16b02f0600000p+0,
    0x1.10b85e7880000p+0,
    0x1.0ae04b2e00000p+0,
    0x1.0527f52680000p+0,
    0x1.ff1eb8c400000p-1,
    0x1.f42a1f1800000p-1,
    0x1.e9721d4900000p-1,
    0x1.def6b35700000p-1,
    0x1.d4b21bf000000p-1,
    0x1.caa739bd00000p-1,
    0x1.c0d32a1500000p-1,
    0x1.b735ecf800000p-1,
    0x1.adcc9fbd00000p-1,
    0x1.a497426400000p-1,
    0x1.9b92f24400000p-1,
    0x1.92c2920600000p-1,
    0x1.8a205c5800000p-1,
    0x1.81af33e300000p-1,
    0x1.796c35fe00000p-1,
    0x1.7154800000000p-1,
};

const r2 = [_]f64{
    0x1.0000000000000p+0,
    0x1.ffa7000000000p-1,
    0x1.ff4f000000000p-1,
    0x1.fef6000000000p-1,
    0x1.fe9e000000000p-1,
    0x1.fe45000000000p-1,
    0x1.fded000000000p-1,
    0x1.fd94000000000p-1,
    0x1.fd3c000000000p-1,
    0x1.fce4000000000p-1,
    0x1.fc8c000000000p-1,
    0x1.fc34000000000p-1,
    0x1.fbdc000000000p-1,
    0x1.fb84000000000p-1,
    0x1.fb2c000000000p-1,
    0x1.fad4000000000p-1,
    0x1.fa7c000000000p-1,
    0x1.fa24000000000p-1,
    0x1.f9cd000000000p-1,
    0x1.f975000000000p-1,
    0x1.f91e000000000p-1,
    0x1.f8c6000000000p-1,
    0x1.f86f000000000p-1,
    0x1.f817000000000p-1,
    0x1.f7c0000000000p-1,
    0x1.f769000000000p-1,
    0x1.f711000000000p-1,
    0x1.f6ba000000000p-1,
    0x1.f663000000000p-1,
    0x1.f60c000000000p-1,
    0x1.f5b5000000000p-1,
    0x1.f55e000000000p-1,
    0x1.f507000000000p-1,
};

const l1 = [_][2]f64{
    .{ 0x0.0p+0, 0x0.0p+0 },
    .{ 0x1.1435edc775b51p-26, 0x1.ffe3800000000p-6 },
    .{ 0x1.c8f1cbf9e4073p-26, 0x1.000bc00000000p-4 },
    .{ -0x1.7bf30fa53957bp-26, 0x1.7ff9400000000p-4 },
    .{ 0x1.674d30b6276edp-27, 0x1.fffa600000000p-4 },
    .{ 0x1.99fcf0d796acep-28, 0x1.4000580000000p-3 },
    .{ -0x1.5ffd8b92706d2p-26, 0x1.7ffe500000000p-3 },
    .{ 0x1.bef90bc5a116dp-26, 0x1.bffc300000000p-3 },
    .{ 0x1.fca73b3d53f0dp-26, 0x1.ffffd00000000p-3 },
    .{ -0x1.41024e560e04ep-27, 0x1.2000d00000000p-2 },
    .{ 0x1.484024fad8461p-26, 0x1.3ffcc00000000p-2 },
    .{ 0x1.2d9ab90ba7694p-26, 0x1.5ffe180000000p-2 },
    .{ -0x1.2c998ea30ba7bp-26, 0x1.80013c0000000p-2 },
    .{ -0x1.75fc8682f918ep-27, 0x1.a0024c0000000p-2 },
    .{ -0x1.f02f268a85fb8p-26, 0x1.bffd300000000p-2 },
    .{ -0x1.50030ea7fae4bp-27, 0x1.dffd7c0000000p-2 },
    .{ -0x1.bffb8da5b849dp-28, 0x1.ffff980000000p-2 },
    .{ 0x1.f71993ff95475p-28, 0x1.0fffd40000000p-1 },
    .{ 0x1.88028e67f78fap-29, 0x1.20010c0000000p-1 },
    .{ 0x1.f915f5a0b4e89p-26, 0x1.30018a0000000p-1 },
    .{ 0x1.c2fe288f968f8p-30, 0x1.3fff480000000p-1 },
    .{ 0x1.7375a75ae0837p-26, 0x1.5001300000000p-1 },
    .{ 0x1.91e48be920323p-29, 0x1.6000ea0000000p-1 },
    .{ 0x1.ee7bc0d39a3dbp-29, 0x1.7001120000000p-1 },
    .{ -0x1.899e2ac5f778cp-29, 0x1.7fffb40000000p-1 },
    .{ 0x1.f1f20176130a7p-26, 0x1.8fffae0000000p-1 },
    .{ 0x1.03fc59d34a4f3p-27, 0x1.9fff1e0000000p-1 },
    .{ -0x1.68722010e4653p-32, 0x1.b0012a0000000p-1 },
    .{ -0x1.58783d505a6ecp-26, 0x1.bffeb80000000p-1 },
    .{ 0x1.b212ab9f8d51dp-27, 0x1.d0004e0000000p-1 },
    .{ -0x1.b3ecd767be776p-26, 0x1.dffeb20000000p-1 },
    .{ 0x1.b1ad41f07fc10p-27, 0x1.effd520000000p-1 },
    .{ 0x0.0p+0, 0x1.0000000000000p+0 },
};

const l2 = [_][2]f64{
    .{ 0x0.0p+0, 0x0.0p+0 },
    .{ -0x1.e2b19f9c7b840p-27, 0x1.00e4000000000p-10 },
    .{ -0x1.b1d68137631fep-27, 0x1.ff10000000000p-10 },
    .{ -0x1.aa92227513fc3p-27, 0x1.8026000000000p-9 },
    .{ -0x1.52918e3ab6f5ep-27, 0x1.ff68000000000p-9 },
    .{ 0x1.1a0b9b9010a9cp-27, 0x1.4019000000000p-8 },
    .{ 0x1.09b87f57867ecp-26, 0x1.7fd0000000000p-8 },
    .{ -0x1.3d6b70c673be6p-28, 0x1.c04c000000000p-8 },
    .{ 0x1.76d340a6780abp-27, 0x1.000c800000000p-7 },
    .{ 0x1.e4181e37d9e05p-26, 0x1.1ff8800000000p-7 },
    .{ -0x1.8c181f042b901p-27, 0x1.3fea800000000p-7 },
    .{ 0x1.9890ae7761d66p-27, 0x1.5fe1800000000p-7 },
    .{ -0x1.55977ae613d5fp-26, 0x1.7fde800000000p-7 },
    .{ 0x1.c75f49acf5e56p-27, 0x1.9fe0800000000p-7 },
    .{ -0x1.e1fc84d0d42bcp-28, 0x1.bfe8800000000p-7 },
    .{ -0x1.5ab8a182ed279p-26, 0x1.dff6000000000p-7 },
    .{ -0x1.b389a7d8a21d9p-26, 0x1.0004800000000p-6 },
    .{ -0x1.7376f1a891fa5p-26, 0x1.1010c00000000p-6 },
    .{ -0x1.0915581b87a8ap-28, 0x1.1ff1000000000p-6 },
    .{ 0x1.b7f2eaa894fcap-27, 0x1.3002c00000000p-6 },
    .{ -0x1.f8c41bdd38c23p-26, 0x1.3fe8c00000000p-6 },
    .{ 0x1.04ffb1f3c8215p-28, 0x1.5000000000000p-6 },
    .{ 0x1.a398ccf2b137dp-26, 0x1.5feb400000000p-6 },
    .{ 0x1.373b5bce07f39p-26, 0x1.7008400000000p-6 },
    .{ -0x1.e42040fd9d454p-27, 0x1.7ff9400000000p-6 },
    .{ 0x1.ea00eb5770526p-27, 0x1.8fecc00000000p-6 },
    .{ -0x1.3e2002bb54bb8p-26, 0x1.a012400000000p-6 },
    .{ 0x1.797ec835db8d3p-26, 0x1.b00b400000000p-6 },
    .{ 0x1.bdcd45fedd285p-28, 0x1.c007400000000p-6 },
    .{ -0x1.33b8f09d1a210p-28, 0x1.d006000000000p-6 },
    .{ -0x1.58937a8eea36fp-27, 0x1.e007800000000p-6 },
    .{ -0x1.3e76203ff54cdp-27, 0x1.f00bc00000000p-6 },
    .{ -0x1.6560f19fc3f41p-30, 0x1.0009600000000p-5 },
};

const c = [_]f64{
    -0x1.62e41d56c6400p-2,
    0x1.47fd2632d2d32p-3,
    -0x1.5504497831ba7p-4,
    0x1.7a3314c5bef3cp-5,
};

const t1 = [_]f64{
    0x1.0000000000000p+0,
    0x1.ea4afa0000000p-1,
    0x1.d5818e0000000p-1,
    0x1.c199be0000000p-1,
    0x1.ae89f98000000p-1,
    0x1.9c49180000000p-1,
    0x1.8ace540000000p-1,
    0x1.7a11470000000p-1,
    0x1.6a09e68000000p-1,
    0x1.5ab07e0000000p-1,
    0x1.4bfdad8000000p-1,
    0x1.3dea650000000p-1,
    0x1.306fe08000000p-1,
    0x1.2387a70000000p-1,
    0x1.172b840000000p-1,
    0x1.0b55870000000p-1,
    0x1.0000000000000p-1,
};

const t2 = [_]f64{
    0x1.0000000000000p+0,
    0x1.fe9d968000000p-1,
    0x1.fd3c228000000p-1,
    0x1.fbdba38000000p-1,
    0x1.fa7c180000000p-1,
    0x1.f91d800000000p-1,
    0x1.f7bfdb0000000p-1,
    0x1.f663278000000p-1,
    0x1.f507658000000p-1,
    0x1.f3ac948000000p-1,
    0x1.f252b38000000p-1,
    0x1.f0f9c20000000p-1,
    0x1.efa1bf0000000p-1,
    0x1.ee4aaa0000000p-1,
    0x1.ecf4830000000p-1,
    0x1.eb9f488000000p-1,
};

const t3 = [_]f64{
    0x1.0000000000000p+0,
    0x1.ffe9d20000000p-1,
    0x1.ffd3a58000000p-1,
    0x1.ffbd798000000p-1,
    0x1.ffa74e8000000p-1,
    0x1.ff91248000000p-1,
    0x1.ff7afb8000000p-1,
    0x1.ff64d38000000p-1,
    0x1.ff4eac8000000p-1,
    0x1.ff38868000000p-1,
    0x1.ff22618000000p-1,
    0x1.ff0c3d0000000p-1,
    0x1.fef61a0000000p-1,
    0x1.fedff78000000p-1,
    0x1.fec9d68000000p-1,
    0x1.feb3b60000000p-1,
};

const t4 = [_]f64{
    0x1.0000000000000p+0,
    0x1.fffe9d0000000p-1,
    0x1.fffd3a0000000p-1,
    0x1.fffbd78000000p-1,
    0x1.fffa748000000p-1,
    0x1.fff9118000000p-1,
    0x1.fff7ae8000000p-1,
    0x1.fff64c0000000p-1,
    0x1.fff4e90000000p-1,
    0x1.fff3860000000p-1,
    0x1.fff2238000000p-1,
    0x1.fff0c08000000p-1,
    0x1.ffef5d8000000p-1,
    0x1.ffedfa8000000p-1,
    0x1.ffec980000000p-1,
    0x1.ffeb350000000p-1,
};

const LL = [4][17][3]f64{
    .{
        .{ 0x0.0p+0, 0x0.0p+0, 0x0.0p+0 },
        .{ 0x1.000001fda0000p-4, 0x1.5ed58a7ff2c40p-42, 0x1.512acbb219717p-96 },
        .{ 0x1.fffffda070000p-4, -0x1.b38ae1c2d5400p-47, -0x1.10f11908c5c8dp-93 },
        .{ 0x1.7ffffefb50000p-3, 0x1.785d91bb08320p-43, 0x1.f5ded30ca48c0p-93 },
        .{ 0x1.0000004a60000p-2, 0x1.380abf8fe7eb0p-42, -0x1.ebd3f567df886p-93 },
        .{ 0x1.4000009760000p-2, 0x1.d321385bf0a10p-43, -0x1.c95a3d59b6c33p-94 },
        .{ 0x1.80000081b4000p-2, 0x1.b6e8f25ff2610p-42, -0x1.4072d39d44270p-93 },
        .{ 0x1.c00000f4f4000p-2, 0x1.a5fed334580a0p-43, -0x1.6f8fdabadd77ap-96 },
        .{ 0x1.ffffff9de0000p-2, -0x1.f85a54fb4a600p-43, -0x1.50083897ee638p-96 },
        .{ 0x1.1fffffa35c000p-1, 0x1.795db9d185140p-42, -0x1.594f07d4d693fp-94 },
        .{ 0x1.3fffff9cbc000p-1, 0x1.2fb09a99477a0p-44, 0x1.b2c6248382152p-94 },
        .{ 0x1.5fffff6df2000p-1, -0x1.3ecf24077c3a0p-42, -0x1.382c6b4e2daa6p-93 },
        .{ 0x1.800000552e000p-1, 0x1.7044806f72380p-43, -0x1.6da3eb126c999p-93 },
        .{ 0x1.9fffffc182000p-1, 0x1.5cf031d5d4cf0p-43, 0x1.302c44da79f4bp-93 },
        .{ 0x1.bfffff6b62000p-1, 0x1.48d11d13372a8p-42, 0x1.f045032543c0dp-93 },
        .{ 0x1.dfffff7a42000p-1, -0x1.0cbf3e62fa860p-43, 0x1.1cfaf84227211p-93 },
        .{ 0x1.0000000000000p+0, 0x0.0p+0, 0x0.0p+0 },
    },
    .{
        .{ 0x0.0p+0, 0x0.0p+0, 0x0.0p+0 },
        .{ 0x1.000024a000000p-8, -0x1.799935d030650p-42, -0x1.cacad0b61c964p-95 },
        .{ 0x1.000014a880000p-7, -0x1.86dd926d35ee8p-42, -0x1.277ab54c00e64p-93 },
        .{ 0x1.7ffff7b380000p-7, 0x1.5ffc7ff00e600p-46, 0x1.5870383c1225dp-93 },
        .{ 0x1.000004b940000p-6, -0x1.4aa9da7ea0120p-43, -0x1.47decfe981022p-93 },
        .{ 0x1.4000064380000p-6, -0x1.015d79fd24c70p-42, 0x1.ec26c231ddb4cp-93 },
        .{ 0x1.7ffff8ff80000p-6, -0x1.a62a2b755f048p-42, -0x1.eebe936e9a3f9p-93 },
        .{ 0x1.c000048180000p-6, 0x1.d1a1bbce5f580p-46, -0x1.4dd25e141ad16p-93 },
        .{ 0x1.0000050ee0000p-5, -0x1.b32450688d5e0p-42, 0x1.577a1c80a65b0p-93 },
        .{ 0x1.2000014760000p-5, -0x1.bac8f6c103d40p-44, -0x1.157a2472c9532p-94 },
        .{ 0x1.3fffff2440000p-5, -0x1.fbab0c9d05aa0p-43, -0x1.3680ffabf1a5bp-94 },
        .{ 0x1.5ffffb1ce0000p-5, 0x1.ef9ec2ba49140p-43, -0x1.24850720b73abp-93 },
        .{ 0x1.7ffffd9600000p-5, 0x1.14e3a91abeff0p-42, -0x1.40e3565c114ddp-94 },
        .{ 0x1.a0000321c0000p-5, 0x1.e33c6a0732340p-44, -0x1.7e8c937a981d8p-93 },
        .{ 0x1.bffffc56a0000p-5, -0x1.2076781ff0460p-42, -0x1.96e054ec8f31ep-93 },
        .{ 0x1.dffffdba40000p-5, -0x1.ea8b4476edde0p-43, -0x1.07bfb9fde770fp-97 },
        .{ 0x1.000001fda0000p-4, 0x1.5ed58a7ff2c40p-42, 0x1.512acbb219717p-96 },
    },
    .{
        .{ 0x0.0p+0, 0x0.0p+0, 0x0.0p+0 },
        .{ 0x1.0002866000000p-12, -0x1.329e31412b688p-42, 0x1.279ddc3585863p-93 },
        .{ 0x1.fffed31000000p-12, 0x1.546dcb40b0518p-42, -0x1.6d048f41debecp-93 },
        .{ 0x1.8000388000000p-11, 0x1.fce6789678c88p-42, 0x1.4959dc2fea030p-93 },
        .{ 0x1.0000668400000p-10, 0x1.c5caec0aa6620p-42, -0x1.89d58a2e0d287p-93 },
        .{ 0x1.4000937000000p-10, -0x1.df91bf5db20f8p-42, 0x1.d4fc35005498ap-93 },
        .{ 0x1.8000a2d800000p-10, 0x1.e7464f400a700p-44, -0x1.6365a147ce4a8p-93 },
        .{ 0x1.c00094a400000p-10, -0x1.c55c526e65a68p-42, -0x1.a7ea2cf16d2cep-93 },
        .{ 0x1.0000345400000p-9, 0x1.8e815121136d0p-43, -0x1.68f8f7a924fd6p-94 },
        .{ 0x1.20000f6400000p-9, 0x1.7fb19a44da298p-42, -0x1.a05fc6ed4124dp-93 },
        .{ 0x1.3fffdb7000000p-9, 0x1.bfac3b1abe2e8p-42, -0x1.9ab1f45a1efb1p-95 },
        .{ 0x1.6000516a00000p-9, -0x1.02544d9faf830p-43, 0x1.286c39adffa8cp-93 },
        .{ 0x1.7fffff4200000p-9, -0x1.2878b893da790p-42, 0x1.b55c5e09022b1p-93 },
        .{ 0x1.a00056f400000p-9, -0x1.9d7ccf4019608p-42, 0x1.f4b05553fe340p-93 },
        .{ 0x1.bfffe65200000p-9, -0x1.ee85f47c2b5a0p-42, 0x1.53fa913d4e9a6p-93 },
        .{ 0x1.e0001f7600000p-9, 0x1.2e594fb2af5e0p-43, -0x1.afaec77ca8d86p-99 },
        .{ 0x1.000024a000000p-8, -0x1.799935d030650p-42, -0x1.cacad0b61c964p-95 },
    },
    .{
        .{ 0x0.0p+0, 0x0.0p+0, 0x0.0p+0 },
        .{ 0x1.0014690000000p-16, -0x1.743d8ea40a4d0p-43, -0x1.75086a21c12abp-97 },
        .{ 0x1.0014c18000000p-15, 0x1.8e7fa3cf59080p-43, 0x1.a11bdcdb1ace3p-93 },
        .{ 0x1.7ff17c8000000p-15, 0x1.118dbbcce8470p-42, 0x1.09806ecd19760p-95 },
        .{ 0x1.fffcbb8000000p-15, -0x1.eaa70b2051e38p-42, 0x1.0d287c5aa12b4p-97 },
        .{ 0x1.4004294000000p-14, 0x1.602b04d79b0c8p-42, 0x1.8d2c8ddf31648p-96 },
        .{ 0x1.800a218000000p-14, -0x1.0a23112178630p-42, -0x1.c4652435118f1p-93 },
        .{ 0x1.bff9304000000p-14, -0x1.385dc475c16b0p-43, -0x1.e7288f74bbbdap-96 },
        .{ 0x1.ffff810000000p-14, 0x1.355edcff4e190p-43, 0x1.f8daae6c1db11p-95 },
        .{ 0x1.2002ff2000000p-13, 0x1.3fbb7cd366140p-45, 0x1.93c0e5b4dca92p-93 },
        .{ 0x1.3ffac90000000p-13, -0x1.9bde416de25c0p-44, -0x1.f4df1d74776b3p-93 },
        .{ 0x1.5ffe340000000p-13, -0x1.21735d0955dc0p-42, -0x1.1b0e10765842fp-93 },
        .{ 0x1.8001b52000000p-13, 0x1.0c0949e968810p-43, -0x1.d3771c0b220f3p-94 },
        .{ 0x1.a0054c8000000p-13, 0x1.26936197b1050p-43, 0x1.e59caa334554ap-93 },
        .{ 0x1.bffd6f0000000p-13, 0x1.27097b0b27d40p-43, 0x1.78a37896a0d1fp-93 },
        .{ 0x1.e00132c000000p-13, 0x1.df1e084120280p-44, 0x1.c61f449a4b5d5p-94 },
        .{ 0x1.0002866000000p-12, -0x1.329e31412b688p-42, 0x1.279ddc3585863p-93 },
    },
};

const cy = [_][2]f64{
    .{ 0x1.71547652b82fep+0, 0x1.777d0ffda0d24p-56 },
    .{ -0x1.71547652b82fep-1, -0x1.777d0ffda0d24p-57 },
    .{ 0x1.ec709dc3a03fdp-2, 0x1.d27e96be541e5p-56 },
};

const cl = [_]f64{
    -0x1.71547652b82fep-2,
    0x1.2776c50f1ff14p-2,
    -0x1.ec709dc3eca5dp-3,
};

const DoubleDouble = struct { high: f64, low: f64 };

fn fastTwoSum(x: f64, y: f64, low: *f64) f64 {
    const sum = x + y;
    const rounded_x = sum - x;
    low.* = y - rounded_x;
    return sum;
}

fn addDd(x_high: f64, x_low: f64, y_high: f64, y_low: f64, low: *f64) f64 {
    const sum = x_high + y_high;
    const delta = sum - x_high;
    low.* = ((y_high - delta) + (x_high + (delta - sum))) + (x_low + y_low);
    return sum;
}

fn mulDdAcc(x_high: f64, x_low: f64, y_high: f64, y_low: f64, low: *f64) f64 {
    const high_high = y_high * x_high;
    const high_low = y_high * x_low;
    const low_high = y_low * x_high;
    var product_low = @mulAdd(f64, y_high, x_high, -high_high);
    product_low += low_high + high_low;
    return fastTwoSum(high_high, product_low, low);
}

fn polyDd(x_high: f64, x_low: f64, coefficients: []const [2]f64, low: *f64) f64 {
    var index = coefficients.len - 1;
    var high = coefficients[index][0] + low.*;
    var correction = (coefficients[index][0] - high) + low.*;
    correction += coefficients[index][1];
    while (index != 0) {
        index -= 1;
        high = mulDdAcc(x_high, x_low, high, correction, &correction);
        const sum = high + coefficients[index][0];
        const sum_low = (coefficients[index][0] - sum) + high;
        high = sum;
        correction += sum_low + coefficients[index][1];
    }
    low.* = correction;
    return high;
}

fn log2dRefine(value: f64, approximation: f64) f64 {
    @setFloatMode(.strict);

    var bits: u64 = @bitCast(value);
    const exponent_field: u64 = bits >> 52;
    var exponent: i32 = @as(i32, @intCast(exponent_field)) - 0x3ff;
    if (exponent_field == 0) {
        const leading_zeros: u32 = @clz(bits);
        exponent -= @as(i32, @intCast(leading_zeros)) - 12;
        bits <<= @intCast(leading_zeros - 11);
    }
    bits &= 0x000f_ffff_ffff_ffff;
    bits |= @as(u64, 0x3ff) << 52;
    const normalized: f64 = @bitCast(bits);
    const exponent_float: f64 = @floatFromInt(exponent);

    const reduced = approximation - exponent_float + 0x1.00008p+0;
    const reduced_bits: u64 = @bitCast(reduced);
    const table_index: usize = @intCast((reduced_bits -% (@as(u64, 0x3ff) << 52)) >> 36);
    const table_1_index = table_index >> 12;
    const table_2_index = (table_index >> 8) & 0xf;
    const table_3_index = (table_index >> 4) & 0xf;
    const table_4_index = table_index & 0xf;

    const logarithm = [3]f64{
        LL[0][table_1_index][0] + LL[1][table_2_index][0] +
            (LL[2][table_3_index][0] + LL[3][table_4_index][0]) + exponent_float,
        LL[0][table_1_index][1] + LL[1][table_2_index][1] +
            (LL[2][table_3_index][1] + LL[3][table_4_index][1]),
        LL[0][table_1_index][2] + LL[1][table_2_index][2] +
            (LL[2][table_3_index][2] + LL[3][table_4_index][2]),
    };

    const table_12 = t1[table_1_index] * t2[table_2_index];
    const table_34 = t3[table_3_index] * t4[table_4_index];
    const table_high = table_12 * table_34;
    const table_low = @mulAdd(f64, table_12, table_34, -table_high);
    const delta_high = table_high * normalized;
    const delta_low = @mulAdd(f64, table_high, normalized, -delta_high);
    const product_high = table_low * normalized;
    const product_low = @mulAdd(f64, table_low, normalized, -product_high);

    var normalized_low: f64 = undefined;
    var normalized_high = fastTwoSum(delta_high - 1.0, delta_low, &normalized_low);
    normalized_high = addDd(normalized_high, normalized_low, product_high, product_low, &normalized_low);

    var polynomial_low = normalized_high *
        (cl[0] + normalized_high * (cl[1] + normalized_high * cl[2]));
    var polynomial_high = polyDd(normalized_high, normalized_low, &cy, &polynomial_low);
    polynomial_high = mulDdAcc(normalized_high, normalized_low, polynomial_high, polynomial_low, &polynomial_low);
    polynomial_high = addDd(polynomial_high, polynomial_low, logarithm[1], logarithm[2], &polynomial_low);

    var rounding_low: f64 = undefined;
    const rounding_high = fastTwoSum(logarithm[0], polynomial_high, &rounding_low);
    const first_low = rounding_low;
    const result_high = fastTwoSum(first_low, polynomial_low, &rounding_low);
    var result_bits: u64 = @bitCast(result_high);
    if ((result_bits & 0x000f_ffff_ffff_ffff) == 0) {
        const correction_bits: u64 = @bitCast(rounding_low);
        if (((correction_bits ^ result_bits) >> 63) != 0) {
            result_bits -%= 1;
        } else {
            result_bits +%= 1;
        }
        return @as(f64, @bitCast(result_bits)) + rounding_high;
    }
    return result_high + rounding_high;
}

fn log2dFinite(value: f64, input_bits: u64) f64 {
    @setFloatMode(.strict);

    var bits = input_bits;
    const exponent_field: u64 = bits >> 52;
    var exponent: i32 = @as(i32, @intCast(exponent_field)) - 0x3ff;
    if (exponent_field == 0) {
        const leading_zeros: u32 = @clz(bits);
        exponent -= @as(i32, @intCast(leading_zeros)) - 12;
        bits <<= @intCast(leading_zeros - 11);
    }
    bits &= 0x000f_ffff_ffff_ffff;
    if (bits == 0) return @floatFromInt(exponent);

    const top_index: usize = @intCast(bits >> 47);
    const low_index: u64 = bits & 0x0000_7fff_ffff_ffff;
    const correction = @as(i64, B[top_index].c1) * @as(i64, @intCast(low_index >> 16));
    const table_bits = bits +% (@as(u64, B[top_index].c0) << 33) +%
        @as(u64, @bitCast(correction));
    const table_index: usize = @intCast(table_bits >> 42);
    const table_1_index = table_index >> 5;
    const table_2_index = table_index & 0x1f;
    bits |= @as(u64, 0x3ff) << 52;
    const normalized: f64 = @bitCast(bits);

    const logarithm_2_high = 0x1.71548p+0;
    const logarithm_2_low = -0x1.ad47a2f472159p-22;
    const reciprocal = r1[table_1_index] * r2[table_2_index];
    const product = reciprocal * normalized;
    const product_low = @mulAdd(f64, reciprocal, normalized, -product);
    const delta_high = product - logarithm_2_high;
    const delta = delta_high + product_low;
    const delta_squared = delta * delta;
    const polynomial = delta_squared *
        ((c[0] + delta * c[1]) + delta_squared * (c[2] + delta * c[3]));
    const exponent_float: f64 = @floatFromInt(exponent);
    const logarithm_tail = (l1[table_1_index][1] + l2[table_2_index][1]) + exponent_float;
    const result_high = logarithm_tail + delta_high;
    var result_low = (logarithm_tail - result_high) + delta_high;
    result_low += ((l1[table_1_index][0] + l2[table_2_index][0]) + product_low) +
        delta_high * logarithm_2_low;
    result_low += polynomial;

    const error_bound = 1.6e-22;
    const lower = result_high + (result_low - error_bound);
    const upper = result_high + (result_low + error_bound);
    if (lower == upper) return lower;
    return log2dRefine(value, upper);
}

// Correctly rounded binary64 log2, ported from CORE-MATH
// src/binary64/log2/log2.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa.
// Copyright (c) 2023 Alexei Sibidanov.
//
// This code is derived from CORE-MATH and used under the MIT License.

pub fn log2d(value: f64) f64 {
    const bits: u64 = @bitCast(value);
    const magnitude = bits & magnitude_mask;
    if (magnitude == 0) return @bitCast(minus_infinity_bits);
    if ((bits & sign_mask) != 0) {
        if (magnitude > infinity_bits) return @bitCast(bits | quiet_nan_bit);
        return @bitCast(canonical_nan_bits);
    }
    if (magnitude >= infinity_bits) {
        if (magnitude == infinity_bits) return value;
        return @bitCast(bits | quiet_nan_bit);
    }
    return log2dFinite(value, bits);
}
