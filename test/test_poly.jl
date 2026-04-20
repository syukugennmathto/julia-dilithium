using Test
using Libdl

include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))

using .Dilithium2.Reduce: reduce32, caddq, freeze
using .Dilithium2.Rounding: power2round, decompose, make_hint, use_hint
using .Dilithium2.Poly: Polynomial, polyeta_pack, polyeta_unpack,
                         polyt1_pack, polyt1_unpack,
                         polyt0_pack, polyt0_unpack,
                         polyz_pack, polyz_unpack
using .Dilithium2.Params: N, Q, D, ETA, GAMMA1, GAMMA2

const C_LIB = joinpath(@__DIR__, "..", "..", "libpqcrystals_dilithium2_ref.dylib")

@testset "reduce32 vs C" begin
    test_vals = Int32[0, 1, -1, Int32(Q-1), Int32(-Q+1),
                      Int32(2*Q), Int32(-2*Q), Int32(100000), Int32(-100000)]
    for v in test_vals
        julia_r = reduce32(v)
        c_r = ccall((:pqcrystals_dilithium2_ref_reduce32, C_LIB), Int32, (Int32,), v)
        if julia_r != c_r
            println("  FAIL reduce32($v): Julia=$julia_r C=$c_r")
        end
        @test julia_r == c_r
    end
    println("✅ reduce32 全て一致")
end

@testset "caddq vs C" begin
    test_vals = Int32[0, 1, -1, Int32(Q-1), Int32(-Q+1), Int32(-100), Int32(100)]
    for v in test_vals
        julia_r = caddq(v)
        c_r = ccall((:pqcrystals_dilithium2_ref_caddq, C_LIB), Int32, (Int32,), v)
        if julia_r != c_r
            println("  FAIL caddq($v): Julia=$julia_r C=$c_r")
        end
        @test julia_r == c_r
    end
    println("✅ caddq 全て一致")
end

@testset "power2round vs C" begin
    test_vals = Int32[0, 1, -1, Int32(Q-1), Int32(1000), Int32(-1000), Int32(Q÷2)]
    for v in test_vals
        julia_a1, julia_a0 = power2round(v)
        c_a0 = Ref(Int32(0))
        c_a1 = ccall((:pqcrystals_dilithium2_ref_power2round, C_LIB), Int32,
                     (Ptr{Int32}, Int32), c_a0, v)
        if julia_a1 != c_a1 || julia_a0 != c_a0[]
            println("  FAIL power2round($v): Julia=($julia_a1,$julia_a0) C=($c_a1,$(c_a0[]))")
        end
        @test julia_a1 == c_a1
        @test julia_a0 == c_a0[]
    end
    println("✅ power2round 全て一致")
end

@testset "decompose vs C" begin
    test_vals = Int32[0, 1, Int32(Q-1), Int32(GAMMA2), Int32(2*GAMMA2),
                      Int32(Q÷2), Int32(100000)]
    for v in test_vals
        julia_a1, julia_a0 = decompose(v)
        c_a0 = Ref(Int32(0))
        c_a1 = ccall((:pqcrystals_dilithium2_ref_decompose, C_LIB), Int32,
                     (Ptr{Int32}, Int32), c_a0, v)
        if julia_a1 != c_a1 || julia_a0 != c_a0[]
            println("  FAIL decompose($v): Julia=($julia_a1,$julia_a0) C=($c_a1,$(c_a0[]))")
        end
        @test julia_a1 == c_a1
        @test julia_a0 == c_a0[]
    end
    println("✅ decompose 全て一致")
end

@testset "make_hint vs C" begin
    test_pairs = [(Int32(0), Int32(0)),
                  (Int32(GAMMA2), Int32(0)),
                  (Int32(GAMMA2+1), Int32(0)),
                  (Int32(Q - GAMMA2), Int32(0)),
                  (Int32(Q - GAMMA2), Int32(1)),
                  (Int32(Q - GAMMA2 + 1), Int32(0))]
    for (a0, a1) in test_pairs
        julia_r = make_hint(a0, a1)
        c_r = ccall((:pqcrystals_dilithium2_ref_make_hint, C_LIB), UInt32,
                    (Int32, Int32), a0, a1)
        if julia_r != Int32(c_r)
            println("  FAIL make_hint($a0,$a1): Julia=$julia_r C=$c_r")
        end
        @test julia_r == Int32(c_r)
    end
    println("✅ make_hint 全て一致")
end

@testset "polyeta pack/unpack roundtrip" begin
    p = Polynomial()
    for i in 1:N
        p.coeffs[i] = Int32(mod(i, 2*ETA+1) - ETA)
    end
    buf = zeros(UInt8, 96)
    polyeta_pack(buf, p)
    p2 = Polynomial()
    polyeta_unpack(p2, buf)
    @test p == p2
    println("✅ polyeta roundtrip 一致")
end

@testset "polyt1 pack/unpack roundtrip" begin
    p = Polynomial()
    for i in 1:N
        p.coeffs[i] = Int32(mod(i, 1024))
    end
    buf = zeros(UInt8, 320)
    polyt1_pack(buf, p)
    p2 = Polynomial()
    polyt1_unpack(p2, buf)
    @test p == p2
    println("✅ polyt1 roundtrip 一致")
end

@testset "polyt0 pack/unpack roundtrip" begin
    p = Polynomial()
    for i in 1:N
        p.coeffs[i] = Int32(mod(i, 8192) - 4095)
    end
    buf = zeros(UInt8, 416)
    polyt0_pack(buf, p)
    p2 = Polynomial()
    polyt0_unpack(p2, buf)
    @test p == p2
    println("✅ polyt0 roundtrip 一致")
end

@testset "polyz pack/unpack roundtrip" begin
    p = Polynomial()
    for i in 1:N
        p.coeffs[i] = Int32(mod(i, 2*GAMMA1) - GAMMA1 + 1)
    end
    buf = zeros(UInt8, 576)
    polyz_pack(buf, p)
    p2 = Polynomial()
    polyz_unpack(p2, buf)
    @test p == p2
    println("✅ polyz roundtrip 一致")
end

println("\n🎉 All poly tests passed!")
