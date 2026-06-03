# SPDX-License-Identifier: Apache-2.0
using Test
using GoldenFloats

@testset "GoldenFloats.jl" begin

    @testset "L5 identity: phi^2 + phi^-2 = 3" begin
        @test isapprox(PHI^2 + PHI^-2, TRINITY; atol=1e-9)
        @test isapprox(PHI * PHI, PHI + 1; atol=1e-9)
        @test lucas(2) == 3
    end

    @testset "Ladder rule e = round((N-1)/phi^2)" begin
        @test GoldenFloats.verify_ladder()
        for fmt in GOLDEN_FLOAT_FAMILY
            @test fmt.sign_bits + fmt.exp_bits + fmt.mant_bits + 1 == fmt.bits + 1
            @test fmt.mant_bits == fmt.bits - 1 - fmt.exp_bits
        end
    end

    @testset "Ladder cardinality and primary" begin
        @test length(GOLDEN_FLOAT_FAMILY) == 9
        @test GF16.bits == 16
        @test GF16.exp_bits == 6
        @test GF16.mant_bits == 9
        @test GF16.bias == 31
    end

    @testset "Lucas recurrence (integer-exact)" begin
        # Hand-verified Lucas seeds
        @test lucas(0)  == 2
        @test lucas(1)  == 1
        @test lucas(2)  == 3
        @test lucas(4)  == 7
        @test lucas(6)  == 18
        @test lucas(8)  == 47
        @test lucas(10) == 123
        @test lucas(12) == 322
    end

    @testset "phi_acc(n) == L_{2n} within f64 tol (F1)" begin
        for n in 0:6
            @test isapprox(phi_acc(n), Float64(lucas(2n)); atol=1e-9)
        end
    end

    @testset "LucasAcc integer-exact integration" begin
        acc = LucasAcc()
        push!(acc, 1)  # adds L_2 = 3
        push!(acc, 2)  # adds L_4 = 7
        push!(acc, 3)  # adds L_6 = 18
        @test sum(acc) == 3 + 7 + 18
        @test sum(acc) == 28
    end

    @testset "GF16 encode/decode round-trip (representable values)" begin
        # Powers of 2 within the GF16 dynamic range
        for x in [1.0, 2.0, 0.5, 4.0, 0.25, 8.0]
            raw = encode(GF16, x)
            y = decode(GF16, raw)
            @test isapprox(y, x; rtol=1e-2)
            raw_neg = encode(GF16, -x)
            y_neg = decode(GF16, raw_neg)
            @test isapprox(y_neg, -x; rtol=1e-2)
        end
        # Zero must round-trip exactly
        @test decode(GF16, encode(GF16, 0.0)) == 0.0
    end

    @testset "GF8 encode/decode round-trip" begin
        for x in [1.0, 2.0, 0.5, 4.0]
            raw = encode(GF8, x)
            y = decode(GF8, raw)
            @test isapprox(y, x; rtol=2e-1)  # coarse precision at 8 bits
        end
    end

    @testset "GF4 ladder presence" begin
        # GF4 is too coarse for meaningful round-trip; just verify metadata
        @test GF4.bits == 4
        @test GF4.exp_bits == 1
        @test GF4.mant_bits == 2
    end

    @testset "Format lookup" begin
        @test get_format("GF16").bits == 16
        @test get_format(32).name == "GF32"
    end

end
