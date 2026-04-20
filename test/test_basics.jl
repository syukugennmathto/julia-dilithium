using Test

include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))

using .Dilithium2

@testset "Params Tests" begin
    @test Dilithium2.Params.N == 256
    @test Dilithium2.Params.Q == 8380417
    @test Dilithium2.Params.CRYPTO_PUBLICKEYBYTES == 1312
end

@testset "SHAKE Tests" begin
    using .Dilithium2.Fips202
    
    input = Vector{UInt8}("test input")
    output256 = shake256(32, input)
    @test length(output256) == 32
    @test output256[1] == 0xe9
end

@testset "Polynomial Tests" begin
    using .Dilithium2.Poly
    
    p = Polynomial()
    @test all(p.coeffs .== 0)
    @test length(p.coeffs) == 256
end

println("\n✅ All basic tests passed!")
