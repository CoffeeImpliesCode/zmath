const sign_mask_h: u16 = 0x8000;
const sign_mask_f: u32 = 0x8000_0000;
const sign_mask_d: u64 = 0x8000_0000_0000_0000;

const magnitude_mask_h: u16 = 0x7fff;
const magnitude_mask_f: u32 = 0x7fff_ffff;
const magnitude_mask_d: u64 = 0x7fff_ffff_ffff_ffff;

const infinity_bits_h: u16 = 0x7c00;
const infinity_bits_f: u32 = 0x7f80_0000;
const infinity_bits_d: u64 = 0x7ff0_0000_0000_0000;

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
