using Libdl

include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))
using .Dilithium2.Sign: crypto_sign_keypair
using .Dilithium2.Packing: unpack_pk
using .Dilithium2.Poly: PolyVec

const C_LIB = joinpath(@__DIR__, "..", "..", "libpqcrystals_dilithium2_ref.dylib")

seed = zeros(UInt8, 32)
for i in 1:32
    seed[i] = UInt8(i)
end

pk, sk = crypto_sign_keypair(seed)

println("=== Julia unpack_pk ===")
rho_j = zeros(UInt8, 32)
t1_j = PolyVec(4)
unpack_pk(rho_j, t1_j, pk)
println("rho[1:8] = $(rho_j[1:8])")
println("t1[1][1:8] = $(t1_j.polys[1].coeffs[1:8])")
println("t1[2][1:8] = $(t1_j.polys[2].coeffs[1:8])")

println("\n=== C unpack_pk ===")
c_rho = zeros(UInt8, 32)
c_t1 = zeros(Int32, 4 * 256)  # polyveck = K polynomials
c_pk = copy(pk)
ccall((:pqcrystals_dilithium2_ref_unpack_pk, C_LIB),
      Nothing,
      (Ptr{UInt8}, Ptr{Int32}, Ptr{UInt8}),
      c_rho, c_t1, c_pk)
println("rho[1:8] = $(c_rho[1:8])")
println("t1[1][1:8] = $(c_t1[1:8])")
println("t1[2][1:8] = $(c_t1[257:264])")

println("\n=== Compare ===")
println("rho match: ", rho_j == c_rho)
println("t1[1] match: ", t1_j.polys[1].coeffs == c_t1[1:256])
println("t1[2] match: ", t1_j.polys[2].coeffs == c_t1[257:512])
println("t1[3] match: ", t1_j.polys[3].coeffs == c_t1[513:768])
println("t1[4] match: ", t1_j.polys[4].coeffs == c_t1[769:1024])

if t1_j.polys[1].coeffs != c_t1[1:256]
    println("\n=== t1[1] differences ===")
    for i in 1:256
        if t1_j.polys[1].coeffs[i] != c_t1[i]
            println("  [$i]: Julia=$(t1_j.polys[1].coeffs[i]) C=$(c_t1[i])")
            if i > 10
                println("  ... (more)")
                break
            end
        end
    end
end
