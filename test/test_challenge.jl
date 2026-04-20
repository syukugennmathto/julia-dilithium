using Libdl

include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))
using .Dilithium2.Poly: Polynomial, poly_challenge
using .Dilithium2.Fips202: shake256

const C_LIB = joinpath(@__DIR__, "..", "..", "libpqcrystals_dilithium2_ref.dylib")

# Test 1: 固定シードでpoly_challengeを比較
println("=== Test poly_challenge with different seeds ===")
for test_idx in 1:3
    seed = zeros(UInt8, 32)
    for i in 1:32
        seed[i] = UInt8(test_idx * i)
    end
    
    # Julia
    jp = Polynomial()
    poly_challenge(jp, seed)
    
    # C
    c_coeffs = zeros(Int32, 256)
    c_seed = copy(seed)
    ccall((:pqcrystals_dilithium2_ref_poly_challenge, C_LIB),
          Nothing, (Ptr{Int32}, Ptr{UInt8}),
          c_coeffs, c_seed)
    
    match = (jp.coeffs == c_coeffs)
    println("\nTest $test_idx: ", match ? "✅ MATCH" : "❌ MISMATCH")
    
    if !match
        # 非ゼロ係数の位置を比較
        julia_nonzero = [i for i in 1:256 if jp.coeffs[i] != 0]
        c_nonzero = [i for i in 1:256 if c_coeffs[i] != 0]
        println("  Julia nonzero count: $(length(julia_nonzero))")
        println("  C     nonzero count: $(length(c_nonzero))")
        println("  Julia nonzero positions[1:10]: $julia_nonzero[1:min(10,end)]")
        println("  C     nonzero positions[1:10]: $c_nonzero[1:min(10,end)]")
        
        for i in 1:256
            if jp.coeffs[i] != c_coeffs[i]
                println("  First diff at [$i]: Julia=$(jp.coeffs[i]) C=$(c_coeffs[i])")
                break
            end
        end
    end
end

# Test 2: SHAKE256の整合性を確認
println("\n=== Test SHAKE256 ===")
test_seed = UInt8[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
                  17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
julia_out = shake256(136, test_seed)  # SHAKE256_RATE = 136
println("Julia SHAKE256[1:16]: $(julia_out[1:16])")

c_out = zeros(UInt8, 136)
c_in = copy(test_seed)
ccall((:pqcrystals_fips202_ref_shake256, C_LIB),
      Nothing, (Ptr{UInt8}, UInt64, Ptr{UInt8}, UInt64),
      c_out, 136, c_in, 32)
println("C     SHAKE256[1:16]: $(c_out[1:16])")
println("SHAKE256 match: ", julia_out == c_out)
