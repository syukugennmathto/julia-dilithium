using Libdl

include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))
using .Dilithium2.Sign: crypto_sign_keypair, crypto_sign_signature, crypto_sign_verify
using .Dilithium2.Packing: unpack_sig
using .Dilithium2.Poly: PolyVec

const C_LIB = joinpath(@__DIR__, "..", "..", "libpqcrystals_dilithium2_ref.dylib")

seed = zeros(UInt8, 32)
for i in 1:32
    seed[i] = UInt8(i)
end

pk, sk = crypto_sign_keypair(seed)
msg = UInt8.("Test" |> collect)
sig = crypto_sign_signature(msg, sk)

println("=== Unpack Julia signature ===")
c_j  = zeros(UInt8, 32)
z_j  = PolyVec(4)
h_j  = PolyVec(4)
ret = unpack_sig(c_j, z_j, h_j, sig)
println("Julia unpack_sig ret = $ret")
println("c[1:8] = $(c_j[1:8])")
println("z[1][1:5] = $(z_j.polys[1].coeffs[1:5])")

# h の非ゼロ数 (global scope)
global h_count = 0
for i in 1:4, j in 1:256
    if h_j.polys[i].coeffs[j] != 0
        global h_count += 1
    end
end
println("h nonzero count = $h_count")

println("\n=== Unpack C signature ===")
c_c  = zeros(UInt8, 32)
z_c  = zeros(Int32, 4 * 256)
h_c  = zeros(Int32, 4 * 256)
c_sig = copy(sig)
ret_c = ccall((:pqcrystals_dilithium2_ref_unpack_sig, C_LIB),
             Int32,
             (Ptr{UInt8}, Ptr{Int32}, Ptr{Int32}, Ptr{UInt8}),
             c_c, z_c, h_c, c_sig)
println("C unpack_sig ret = $ret_c")
println("c[1:8] = $(c_c[1:8])")
println("z[1][1:5] = $(z_c[1:5])")
h_count_c = sum(h_c .!= 0)
println("h nonzero count = $h_count_c")

println("\n=== Compare ===")
println("c match: ", c_j == c_c)
println("z[1] match: ", z_j.polys[1].coeffs == z_c[1:256])
println("z[2] match: ", z_j.polys[2].coeffs == z_c[257:512])
println("z[3] match: ", z_j.polys[3].coeffs == z_c[513:768])
println("z[4] match: ", z_j.polys[4].coeffs == z_c[769:1024])

# h を詳細比較
h_j_flat = Int32[]
for i in 1:4
    append!(h_j_flat, h_j.polys[i].coeffs)
end
h_match = (h_j_flat == h_c)
println("h match: $h_match")

if !h_match
    println("\n=== h detailed comparison ===")
    diff_count = 0
    for i in 1:min(1024, length(h_j_flat))
        if h_j_flat[i] != h_c[i]
            poly_idx = (i-1) ÷ 256
            coeff_idx = (i-1) % 256
            println("  Diff at poly=$poly_idx coeff=$coeff_idx: Julia=$(h_j_flat[i]) C=$(h_c[i])")
            diff_count += 1
            if diff_count > 10  # 最初の10個だけ表示
                println("  ... (more differences)")
                break
            end
        end
    end
end

println("\n=== Try C verify ===")
c_ret = ccall((:pqcrystals_dilithium2_ref_verify, C_LIB),
              Int32,
              (Ptr{UInt8}, UInt64, Ptr{UInt8}, UInt64, Ptr{UInt8}),
              copy(sig), UInt64(length(sig)),
              copy(msg), UInt64(length(msg)),
              copy(pk))
println("C verify result = $c_ret (0=success, -1=fail)")

# Julia verify for comparison
j_ret = crypto_sign_verify(sig, msg, pk)
println("Julia verify result = $j_ret")
