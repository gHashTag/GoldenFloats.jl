# GoldenFloats.jl

[![CI](https://github.com/gHashTag/GoldenFloats.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/gHashTag/GoldenFloats.jl/actions/workflows/CI.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

Julia reference implementation of the **GoldenFloat** family of phi-structured floating-point formats. Companion to the Verilog RTL in [`gHashTag/tt-trinity-gamma`](https://github.com/gHashTag/tt-trinity-gamma), the canonical `.t27` SSOT in [`gHashTag/t27`](https://github.com/gHashTag/t27), and the Zig reference in [`gHashTag/zig-golden-float`](https://github.com/gHashTag/zig-golden-float).

## What this is

A small, dependency-free Julia harness for the GoldenFloat width ladder (GF4 ... GF256) and the Lucas-EII integer-exact accumulator. The package is intentionally minimal so it can be dropped into existing numerical linear algebra benchmark suites (e.g. [`takum-arithmetic/MuFoLAB`](https://github.com/takum-arithmetic/MuFoLAB)) for side-by-side comparison against `posit`, `takum`, `OCP MX`, and IEEE binary formats.

## Claim status

- **`[Verified]`** -- the one closed-form ladder rule `e = round((N-1)/phi^2)` reproduces all 9 GoldenFloat field-width splits (GF4, GF8, GF12, GF16, GF20, GF24, GF32, GF64, GF256) exactly.
- **`[Verified]`** -- the identity `phi^2 + phi^-2 = 3 = L_2` (Lucas L_2, classical 1878, not original to this project). More generally `phi^(2n) + phi^(-2n) = L_{2n}` is verified to 60 decimal digits with `mpmath` over `n=0..12`; max residual `7.1e-56`. See `lucas_accumulator.t27` in the t27 SSOT.
- **`[Open conjecture]`** -- whether the phi-ladder is numerically superior to an equally-tuned posit / OCP-MX / takum family. This is FL-002 in [`gHashTag/trios-trainer-igla` `src/ledger.rs`](https://github.com/gHashTag/trios-trainer-igla) and is **not** settled by this package. The point of this Julia port is to make the conjecture **testable side-by-side** under matched harnesses.

## Install

```julia
using Pkg
Pkg.add(url="https://github.com/gHashTag/GoldenFloats.jl")
```

## Quick start

```julia
using GoldenFloats

# 1. Inspect the 9-rung ladder
for fmt in GOLDEN_FLOAT_FAMILY
    println(fmt)
end

# 2. The single closed rule that derives every width
GoldenFloats.verify_ladder()    # -> true

# 3. Encode / decode a GF16 value
raw = encode(GF16, 1.5)
decode(GF16, raw)               # -> 1.5 (rtol 1e-2)

# 4. Lucas-EII integer-exact accumulator
acc = LucasAcc()
push!(acc, 1)   # add L_2 = 3
push!(acc, 2)   # add L_4 = 7
push!(acc, 3)   # add L_6 = 18
sum(acc)        # -> 28 (BigInt, exact)

# 5. The L5 identity that anchors the family
phi_acc(1)      # 3.0   == lucas(2)
phi_acc(6)      # 322.0 == lucas(12)
```

## Format ladder

| Name   | Bits | Sign | Exp | Mant | Bias            |
|--------|-----:|-----:|----:|-----:|----------------:|
| GF4    |    4 |    1 |   1 |    2 |               0 |
| GF8    |    8 |    1 |   3 |    4 |               3 |
| GF12   |   12 |    1 |   4 |    7 |               7 |
| GF16   |   16 |    1 |   6 |    9 |              31 |
| GF20   |   20 |    1 |   7 |   12 |              63 |
| GF24   |   24 |    1 |   9 |   14 |             255 |
| GF32   |   32 |    1 |  12 |   19 |            2047 |
| GF64   |   64 |    1 |  24 |   39 |         8388607 |
| GF256  |  256 |    1 |  97 |  158 |       2^96 - 1  |

Field widths follow `t27/specs/numeric/goldenfloat_family.t27` verbatim. Bias values match `t27/specs/numeric/gf{N}.t27`.

## Why Julia

The numerical linear algebra benchmark ecosystem this package targets is Julia-native: [`Posits.jl`](https://github.com/milankl/Posits.jl), [`Takums.jl`](https://github.com/takum-arithmetic/Takums.jl), `SuiteSparse.jl`, [`MuFoLAB`](https://github.com/takum-arithmetic/MuFoLAB). A Julia reference impl drops directly into those harnesses without a foreign-call boundary.

## Related work and prior art

- **Posits** (J. L. Gustafson, *Unum III*, 2017) -- tapered precision, regime-coded exponent.
- **Takums** (L. Hunhold, *Beating Posits at Their Own Game*, CoNGA 2024, [arXiv:2404.18603](https://arxiv.org/abs/2404.18603); *Integer Representations of Real Numbers*, 2024, [arXiv:2412.20273](https://arxiv.org/abs/2412.20273)) -- IEEE-754-backward-compatible tapered.
- **OCP MX** (Rouhani et al., 2023, [arXiv:2310.10537](https://arxiv.org/abs/2310.10537)) -- block-scaled microfloats.
- **Tekum balanced ternary** (Hunhold, 2025, [arXiv:2512.10964](https://arxiv.org/abs/2512.10964)).
- **IEEE P3109** -- ML binary8 working group; FLoPS Lean formalization at [arXiv:2602.15965](https://arxiv.org/abs/2602.15965).
- **Zeckendorf / base-phi integer arithmetic** (Ahlbach, Usatine, Pippenger, 2012, [arXiv:1207.4497](https://arxiv.org/abs/1207.4497)) -- linear-time arithmetic in base phi, the theoretical anchor for the Lucas-EII track.

## License

Apache-2.0. See [LICENSE](LICENSE).

## Citation

If you use this package in research, please cite the upcoming Corona open-silicon paper (in preparation, ARITH 2027 target):

```bibtex
@misc{vasilev2026goldenfloats,
  author = {Vasilev, Dmitrii},
  title  = {GoldenFloats.jl: A phi-structured floating-point reference in Julia},
  year   = {2026},
  howpublished = {\url{https://github.com/gHashTag/GoldenFloats.jl}},
  note   = {Apache-2.0}
}
```

## Contact

Dmitrii Vasilev -- `admin@t27.ai` -- GitHub [`@gHashTag`](https://github.com/gHashTag)
ORCID [0009-0008-4294-6159](https://orcid.org/0009-0008-4294-6159)
