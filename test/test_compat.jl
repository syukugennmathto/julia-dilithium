include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))
using .Dilithium2.Sign: crypto_sign_keypair, crypto_sign_signature, crypto_sign_verify
using .Dilithium2.Packing: pack_pk, unpack_pk, pack_sk, unpack_sk, pack_sig, unpack_sig
using .Dilithium2.Poly: Polynomial, PolyVec

const C_LIB = joinpath(@__DIR__, "..", "..", "libpqcrystals_dilithium2_ref.dylib")

seed = zeros(UInt8, 32)
for i in 1:32
    seed[i] = UInt8(i)
end

println("=== Julia keypair → C verify ===")
pk, sk = crypto_sign_keypair(seed)
msg = UInt8.("Hello Dilithium!" |> collect)

println("\n=== Check pk content BEFORE sign ===")
rho_check = zeros(UInt8, 32)
t1_check = PolyVec(4)
unpack_pk(rho_check, t1_check, pk)
println("pk size = $(length(pk))")
println("t1_check[1].coeffs[1:10] = $(t1_check.polys[1].coeffs[1:10])")
println("t1_check[1].coeffs[100:110] = $(t1_check.polys[1].coeffs[100:110])")
# t1がゼロなら問題！

# もしt1がゼロなら、keypairの中身を見る
println("\n=== Direct t1 from keypair (before packing) ===")
# keypair内部でt1を確認するために、中間値をprintする必要がある
# → sign.jl を一時的に修正するか、ここでkeypairを再実装して確認

sig = crypto_sign_signature(msg, sk)

c_sig = copy(sig)
c_msg = copy(msg)
c_pk  = copy(pk)
ret = ccall((:pqcrystals_dilithium2_ref_verify, C_LIB),
            Int32,
            (Ptr{UInt8}, UInt64, Ptr{UInt8}, UInt64, Ptr{UInt8}),
            c_sig, UInt64(length(sig)),
            c_msg, UInt64(length(msg)),
            c_pk)
println("\nC verify(Julia sig, Julia pk) = $ret  (0=success, -1=fail)")
