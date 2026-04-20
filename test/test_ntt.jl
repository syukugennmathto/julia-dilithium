using Test
using Libdl

include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))

using .Dilithium2.NTT: ZETAS, ntt!, invntt!, montgomery_reduce
using .Dilithium2.Params: N, Q

const C_LIB = joinpath(@__DIR__, "..", "..", "libpqcrystals_dilithium2_ref.dylib")

@testset "Montgomery Reduce vs C" begin
    test_values = Int64[0, 1, -1, Int64(Q), -Int64(Q),
        Int64(25847) * Int64(100), Int64(-2608894) * Int64(50)]
    for val in test_values
        c_result = ccall((:pqcrystals_dilithium2_ref_montgomery_reduce, C_LIB), Int32, (Int64,), val)
        @test montgomery_reduce(val) == c_result
    end
    println("✅ montgomery_reduce 全て一致")
end

@testset "NTT vs C実装" begin
    a_julia = Int32.(mod.(0:255, 100))
    a_c = copy(a_julia)

    ntt!(a_julia)
    ccall((:pqcrystals_dilithium2_ref_ntt, C_LIB), Cvoid, (Ptr{Int32},), a_c)

    println("Julia NTT[1:5]: ", a_julia[1:5])
    println("C     NTT[1:5]: ", a_c[1:5])
    @test a_julia == a_c
    println("✅ NTT 完全一致")
end

@testset "InvNTT vs C実装" begin
    # まず順NTTで変換したデータを使う
    a_julia = Int32.(mod.(0:255, 100))
    ntt!(a_julia)
    a_c = copy(a_julia)

    invntt!(a_julia)
    ccall((:pqcrystals_dilithium2_ref_invntt_tomont, C_LIB), Cvoid, (Ptr{Int32},), a_c)

    println("Julia invNTT[1:5]: ", a_julia[1:5])
    println("C     invNTT[1:5]: ", a_c[1:5])
    @test a_julia == a_c
    println("✅ InvNTT 完全一致")
end

println("\n🎉 All NTT tests passed!")
