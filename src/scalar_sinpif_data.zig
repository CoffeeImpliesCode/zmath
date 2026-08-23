// CORE-MATH binary32 sinpi tables.
//
// Copyright (c) 2022-2025 Alexei Sibidanov.
//
// Derived from CORE-MATH src/binary32/sinpi/sinpif.c at revision
// 3bafed75e3c22399531857d63a3fa0e24bfb42fa and used under the MIT License.

pub const sn = [_]f64{
    0x1.921fb54442d0fp-37, -0x1.4abbce6102b94p-112, 0x1.4669fa3c58463p-189,
};

pub const cn = [_]f64{
    -0x1.3bd3cc9be45cfp-74, 0x1.03c1f08088742p-150, -0x1.55d1e5eff55a5p-228,
};

pub const s = [_]f64{
    0x0p+0,                0x1.91f65f10dd814p-5,  0x1.917a6bc29b42cp-4,  0x1.2c8106e8e613ap-3,
    0x1.8f8b83c69a60bp-3,  0x1.f19f97b215f1bp-3,  0x1.294062ed59f06p-2,  0x1.58f9a75ab1fddp-2,
    0x1.87de2a6aea963p-2,  0x1.b5d1009e15ccp-2,   0x1.e2b5d3806f63bp-2,  0x1.073879922ffeep-1,
    0x1.1c73b39ae68c8p-1,  0x1.30ff7fce17035p-1,  0x1.44cf325091dd6p-1,  0x1.57d69348cecap-1,
    0x1.6a09e667f3bcdp-1,  0x1.7b5df226aafafp-1,  0x1.8bc806b151741p-1,  0x1.9b3e047f38741p-1,
    0x1.a9b66290ea1a3p-1,  0x1.b728345196e3ep-1,  0x1.c38b2f180bdb1p-1,  0x1.ced7af43cc773p-1,
    0x1.d906bcf328d46p-1,  0x1.e212104f686e5p-1,  0x1.e9f4156c62ddap-1,  0x1.f0a7efb9230d7p-1,
    0x1.f6297cff75cbp-1,   0x1.fa7557f08a517p-1,  0x1.fd88da3d12526p-1,  0x1.ff621e3796d7ep-1,
    0x1p+0,                0x1.ff621e3796d7ep-1,  0x1.fd88da3d12526p-1,  0x1.fa7557f08a517p-1,
    0x1.f6297cff75cbp-1,   0x1.f0a7efb9230d7p-1,  0x1.e9f4156c62ddap-1,  0x1.e212104f686e5p-1,
    0x1.d906bcf328d46p-1,  0x1.ced7af43cc773p-1,  0x1.c38b2f180bdb1p-1,  0x1.b728345196e3ep-1,
    0x1.a9b66290ea1a3p-1,  0x1.9b3e047f38741p-1,  0x1.8bc806b151741p-1,  0x1.7b5df226aafafp-1,
    0x1.6a09e667f3bcdp-1,  0x1.57d69348cecap-1,   0x1.44cf325091dd6p-1,  0x1.30ff7fce17035p-1,
    0x1.1c73b39ae68c8p-1,  0x1.073879922ffeep-1,  0x1.e2b5d3806f63bp-2,  0x1.b5d1009e15ccp-2,
    0x1.87de2a6aea963p-2,  0x1.58f9a75ab1fddp-2,  0x1.294062ed59f06p-2,  0x1.f19f97b215f1bp-3,
    0x1.8f8b83c69a60bp-3,  0x1.2c8106e8e613ap-3,  0x1.917a6bc29b42cp-4,  0x1.91f65f10dd814p-5,
    0x0p+0,                -0x1.91f65f10dd814p-5, -0x1.917a6bc29b42cp-4, -0x1.2c8106e8e613ap-3,
    -0x1.8f8b83c69a60bp-3, -0x1.f19f97b215f1bp-3, -0x1.294062ed59f06p-2, -0x1.58f9a75ab1fddp-2,
    -0x1.87de2a6aea963p-2, -0x1.b5d1009e15ccp-2,  -0x1.e2b5d3806f63bp-2, -0x1.073879922ffeep-1,
    -0x1.1c73b39ae68c8p-1, -0x1.30ff7fce17035p-1, -0x1.44cf325091dd6p-1, -0x1.57d69348cecap-1,
    -0x1.6a09e667f3bcdp-1, -0x1.7b5df226aafafp-1, -0x1.8bc806b151741p-1, -0x1.9b3e047f38741p-1,
    -0x1.a9b66290ea1a3p-1, -0x1.b728345196e3ep-1, -0x1.c38b2f180bdb1p-1, -0x1.ced7af43cc773p-1,
    -0x1.d906bcf328d46p-1, -0x1.e212104f686e5p-1, -0x1.e9f4156c62ddap-1, -0x1.f0a7efb9230d7p-1,
    -0x1.f6297cff75cbp-1,  -0x1.fa7557f08a517p-1, -0x1.fd88da3d12526p-1, -0x1.ff621e3796d7ep-1,
    -0x1p+0,               -0x1.ff621e3796d7ep-1, -0x1.fd88da3d12526p-1, -0x1.fa7557f08a517p-1,
    -0x1.f6297cff75cbp-1,  -0x1.f0a7efb9230d7p-1, -0x1.e9f4156c62ddap-1, -0x1.e212104f686e5p-1,
    -0x1.d906bcf328d46p-1, -0x1.ced7af43cc773p-1, -0x1.c38b2f180bdb1p-1, -0x1.b728345196e3ep-1,
    -0x1.a9b66290ea1a3p-1, -0x1.9b3e047f38741p-1, -0x1.8bc806b151741p-1, -0x1.7b5df226aafafp-1,
    -0x1.6a09e667f3bcdp-1, -0x1.57d69348cecap-1,  -0x1.44cf325091dd6p-1, -0x1.30ff7fce17035p-1,
    -0x1.1c73b39ae68c8p-1, -0x1.073879922ffeep-1, -0x1.e2b5d3806f63bp-2, -0x1.b5d1009e15ccp-2,
    -0x1.87de2a6aea963p-2, -0x1.58f9a75ab1fddp-2, -0x1.294062ed59f06p-2, -0x1.f19f97b215f1bp-3,
    -0x1.8f8b83c69a60bp-3, -0x1.2c8106e8e613ap-3, -0x1.917a6bc29b42cp-4, -0x1.91f65f10dd814p-5,
};
