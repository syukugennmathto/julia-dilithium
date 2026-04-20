include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))

using .Dilithium2.Sign: crypto_sign_keypair, crypto_sign_signature, crypto_sign_verify

println("🔑 Keypair生成中...")
seed = zeros(UInt8, 32)  # 固定シード（再現性のため）
for i in 1:32
    seed[i] = UInt8(i)
end

pk, sk = crypto_sign_keypair(seed)
println("  pk: $(length(pk)) bytes")
println("  sk: $(length(sk)) bytes")
println("  pk[1:8] = $(pk[1:8])")
println("  sk[1:8] = $(sk[1:8])")

println("\n✍️  署名中...")
msg = UInt8.("Hello Dilithium!" |> collect)
sig = crypto_sign_signature(msg, sk)
println("  sig: $(length(sig)) bytes")
println("  sig[1:8] = $(sig[1:8])")

println("\n✅ 検証中...")
ret = crypto_sign_verify(sig, msg, pk)
println("  verify result = $ret")
if ret == 0
    println("🎉 署名検証成功!")
else
    println("❌ 署名検証失敗")
end

# 偽メッセージの検証（失敗であるべき）
println("\n🔍 偽メッセージの検証...")
fake_msg = UInt8.("Hello FAKE!" |> collect)
ret2 = crypto_sign_verify(sig, fake_msg, pk)
println("  verify(fake) = $ret2")
if ret2 == -1
    println("✅ 正しく拒否された")
else
    println("❌ 偽メッセージが受け入れられた（バグ）")
end
