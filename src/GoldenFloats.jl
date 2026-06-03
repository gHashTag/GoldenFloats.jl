# SPDX-License-Identifier: Apache-2.0
# GoldenFloats.jl -- phi-structured floating-point reference implementation in Julia
#
# Companion to:
#   gHashTag/tt-trinity-gamma  (Verilog RTL, Apache-2.0)
#   gHashTag/t27               (canonical .t27 SSOT)
#   gHashTag/zig-golden-float  (Zig reference impl, MIT)
#
# Bit-exact match to:
#   t27/specs/numeric/goldenfloat_family.t27
#   t27/specs/numeric/gf{4,8,12,16,20,24,32,64}.t27
#   t27/specs/numeric/lucas_accumulator.t27
#
# CLAIM STATUS:
#   The GoldenFloat ladder rule e = round((N-1)/phi^2) and the identity
#   phi^2 + phi^-2 = 3 = L_2 are [Verified]. Whether the ladder is
#   numerically superior to posit / OCP-MX / takum is [Open conjecture]
#   (FL-002 in gHashTag/trios-trainer-igla src/ledger.rs).
#
#   This package provides a reference Julia harness so the conjecture
#   can be tested side-by-side against MuFoLAB, libtakum, Posits.jl.

module GoldenFloats

export GFFormat, GOLDEN_FLOAT_FAMILY, get_format, encode, decode,
       LucasAcc, lucas, phi_acc, PHI, PHI_INV, TRINITY,
       GF4, GF8, GF12, GF16, GF20, GF24, GF32, GF64, GF256

const PHI     = 1.6180339887498948482
const PHI_INV = 0.6180339887498948482
const TRINITY = 3.0  # L_2 == phi^2 + phi^-2 == 3

# ------------------------------------------------------------------
# 1. Format descriptor -- mirrors t27 GoldenFloatFormat struct
# ------------------------------------------------------------------

struct GFFormat
    name      :: String
    bits      :: Int
    sign_bits :: Int
    exp_bits  :: Int
    mant_bits :: Int
    bias      :: BigInt
end

# Canonical ladder. Field widths follow t27 goldenfloat_family.t27 verbatim.
# The single rule that produces every rung:
#   exp_bits = round((N-1)/phi^2)
#   mant_bits = N - 1 - exp_bits
#   bias = 2^(exp_bits-1) - 1   (IEEE-style; explicit overrides below match t27)
const GOLDEN_FLOAT_FAMILY = [
    GFFormat("GF4",   4,   1, 1,   2,   0),    # t27 gf4.t27: EXP_BIAS=0
    GFFormat("GF8",   8,   1, 3,   4,   3),    # t27 gf8.t27: EXP_BIAS=3
    GFFormat("GF12",  12,  1, 4,   7,   7),    # 2^(4-1)-1 = 7
    GFFormat("GF16",  16,  1, 6,   9,   31),   # t27 gf16.t27: BIAS=31
    GFFormat("GF20",  20,  1, 7,   12,  63),   # 2^(7-1)-1 = 63
    GFFormat("GF24",  24,  1, 9,   14,  255),  # 2^(9-1)-1 = 255
    GFFormat("GF32",  32,  1, 12,  19,  2047), # t27 gf32.t27: EXP_BIAS=2047
    GFFormat("GF64",  64,  1, 24,  39,  8388607),    # 2^(24-1)-1
    GFFormat("GF256", 256, 1, 97,  158, (BigInt(1)<<96)-1),
]

# Convenience constants
const GF4   = GOLDEN_FLOAT_FAMILY[1]
const GF8   = GOLDEN_FLOAT_FAMILY[2]
const GF12  = GOLDEN_FLOAT_FAMILY[3]
const GF16  = GOLDEN_FLOAT_FAMILY[4]
const GF20  = GOLDEN_FLOAT_FAMILY[5]
const GF24  = GOLDEN_FLOAT_FAMILY[6]
const GF32  = GOLDEN_FLOAT_FAMILY[7]
const GF64  = GOLDEN_FLOAT_FAMILY[8]
const GF256 = GOLDEN_FLOAT_FAMILY[9]

# ------------------------------------------------------------------
# 2. Ladder-rule verification
# ------------------------------------------------------------------

"""
    ladder_rule(N)

The single closed-form rule for the GoldenFloat exponent width:
exp_bits = round((N-1)/phi^2). Returns the exp_bits for an N-bit format.
"""
ladder_rule(N::Int) = round(Int, (N-1)/(PHI*PHI))

"""
    verify_ladder()

Confirms that every entry in GOLDEN_FLOAT_FAMILY satisfies the rule.
Returns true if all 9 widths match; false otherwise.
"""
function verify_ladder()
    all(fmt -> fmt.exp_bits == ladder_rule(fmt.bits) &&
              fmt.mant_bits == fmt.bits - 1 - fmt.exp_bits,
        GOLDEN_FLOAT_FAMILY)
end

# ------------------------------------------------------------------
# 3. Encode / decode -- bit-exact match to t27 gf{N}.t27 modules
# ------------------------------------------------------------------

get_format(name::String) = first(filter(f -> f.name == name, GOLDEN_FLOAT_FAMILY))
get_format(bits::Int)    = first(filter(f -> f.bits == bits, GOLDEN_FLOAT_FAMILY))

"""
    encode(fmt::GFFormat, x::Real) -> UInt

Encode a real number x into the GoldenFloat format `fmt`. Returns an
unsigned bit pattern of width `fmt.bits` (packed in the smallest UInt
that fits). Special cases: x == 0 -> 0; subnormals are flushed to zero.
"""
function encode(fmt::GFFormat, x::Real)
    bits = fmt.bits
    sign = x < 0 ? UInt(1) : UInt(0)
    ax = abs(Float64(x))
    if ax == 0.0
        return zero_pattern(bits)
    end
    e_unbiased = floor(Int, log2(ax))
    # bias is BigInt to accommodate GF256 (2^96-1); narrow to Int for arithmetic
    # on the small-ladder rungs where the value fits in Int64.
    bias_int = fmt.bits <= 64 ? Int(fmt.bias) : error("encode not implemented for GF$(fmt.bits) (bias exceeds Int64)")
    e_biased = e_unbiased + bias_int
    emax = (1 << fmt.exp_bits) - 1
    if e_biased >= emax
        # clamp to max representable (saturate, no Inf in early ladder rungs)
        e_biased = emax - 1
    elseif e_biased <= 0
        # flush-to-zero subnormal
        return zero_pattern(bits)
    end
    # mantissa: round(2^m * (ax / 2^e_unbiased - 1))
    frac = ax / (2.0^e_unbiased) - 1.0
    mant = round(Int, frac * (1 << fmt.mant_bits))
    if mant >= (1 << fmt.mant_bits)
        mant = (1 << fmt.mant_bits) - 1
    end
    raw = (sign << (bits - 1)) | (UInt(e_biased) << fmt.mant_bits) | UInt(mant)
    return pack(bits, raw)
end

"""
    decode(fmt::GFFormat, raw::Integer) -> Float64

Decode a packed GoldenFloat bit pattern back to a Float64. Inverse of `encode`
within the format's representable range.
"""
function decode(fmt::GFFormat, raw::Integer)
    bits = fmt.bits
    mant_mask = (UInt(1) << fmt.mant_bits) - 1
    exp_mask  = (UInt(1) << fmt.exp_bits)  - 1
    sign     = (raw >> (bits - 1)) & 1
    e_biased = (raw >> fmt.mant_bits) & exp_mask
    mant     = raw & mant_mask
    if e_biased == 0 && mant == 0
        return 0.0
    end
    bias_int = fmt.bits <= 64 ? Int(fmt.bias) : error("decode not implemented for GF$(fmt.bits) (bias exceeds Int64)")
    e_unbiased = Int(e_biased) - bias_int
    val = (1.0 + Float64(mant) / Float64(1 << fmt.mant_bits)) * (2.0^e_unbiased)
    return sign == 1 ? -val : val
end

# Internal helpers
zero_pattern(bits::Int) = bits <= 64 ? UInt(0) : BigInt(0)
function pack(bits::Int, raw::Integer)
    bits <= 8  && return UInt8(raw)
    bits <= 16 && return UInt16(raw)
    bits <= 32 && return UInt32(raw)
    bits <= 64 && return UInt64(raw)
    return BigInt(raw)
end

# ------------------------------------------------------------------
# 4. Lucas-EII accumulator -- the [Verified] arithmetic leg
# ------------------------------------------------------------------

"""
    lucas(k)

The k-th Lucas number via the exact integer recurrence L_0=2, L_1=1,
L_k = L_{k-1} + L_{k-2}. Returns an Int (or BigInt for k > 91).
"""
function lucas(k::Integer)
    k == 0 && return 2
    k == 1 && return 1
    if k > 91
        a, b = BigInt(2), BigInt(1)
    else
        a, b = 2, 1
    end
    for _ in 2:k
        a, b = b, a + b
    end
    return b
end

"""
    phi_acc(n)

The accumulator value phi^(2n) + phi^(-2n), computed in Float64.
By the L5 identity this equals the integer Lucas number L_{2n} exactly.
"""
phi_acc(n::Integer) = PHI^(2n) + PHI^(-2n)

"""
    LucasAcc

An integer-exact phi-scaled accumulator: stores the running sum as a
plain BigInt, exploiting phi^(2n) + phi^(-2n) = L_{2n}.

This is the integer-backed-accumulation engineering point: a phi-scaled
partial-sum can carry in unsigned integer storage with no rounding error,
unlike posit / OCP-MX / bf16 floats.
"""
mutable struct LucasAcc
    sum :: BigInt
    LucasAcc() = new(BigInt(0))
end

"""
    push!(acc, n) -- add L_{2n} = phi^(2n)+phi^(-2n) to the running sum, exactly.
"""
function Base.push!(acc::LucasAcc, n::Integer)
    acc.sum += lucas(2n)
    return acc
end

Base.sum(acc::LucasAcc) = acc.sum

end # module
