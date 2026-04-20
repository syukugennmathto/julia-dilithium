include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))
using .Dilithium2.Sign: crypto_sign_keypair, crypto_sign_signature
using .Dilithium2.Params

seed = zeros(UInt8, 32)
for i in 1:32
    seed[i] = UInt8(i)
end

pk, sk = crypto_sign_keypair(seed)
msg = UInt8.("Test" |> collect)
sig = crypto_sign_signature(msg, sk)

println("Expected sizes from params.jl:")
println("  CRYPTO_PUBLICKEYBYTES  = $CRYPTO_PUBLICKEYBYTES")
println("  CRYPTO_SECRETKEYBYTES  = $CRYPTO_SECRETKEYBYTES")
println("  CRYPTO_BYTES           = $CRYPTO_BYTES")

println("\nActual sizes:")
println("  pk  = $(length(pk))")
println("  sk  = $(length(sk))")
println("  sig = $(length(sig))")

println("\nMatch:")
println("  pk:  ", length(pk) == CRYPTO_PUBLICKEYBYTES)
println("  sk:  ", length(sk) == CRYPTO_SECRETKEYBYTES)
println("  sig: ", length(sig) == CRYPTO_BYTES)

# C側のexpected sizes
println("\nC expected (from params.h):")
println("  CRYPTO_PUBLICKEYBYTES  = 1312")
println("  CRYPTO_SECRETKEYBYTES  = 2544")  # NOT 2560!
println("  CRYPTO_BYTES           = 2420")
