using Libdl

include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))
using .Dilithium2.Sign: crypto_sign_keypair, crypto_sign_signature, crypto_sign_verify

const C_LIB = joinpath(@__DIR__, "..", "..", "libpqcrystals_dilithium2_ref.dylib")

println("="^60)
println("  Julia Dilithium2 完全互換性テスト")
println("="^60)

seed = zeros(UInt8, 32)
for i in 1:32
    seed[i] = UInt8(i)
end

# =====================================================
# Test 1: Julia keypair → Julia sign → Julia verify
# =====================================================
println("\n[Test 1] Julia keypair → Julia sign → Julia verify")
pk, sk = crypto_sign_keypair(seed)
msg = UInt8.("Hello Dilithium2!" |> collect)
sig = crypto_sign_signature(msg, sk)
ret = crypto_sign_verify(sig, msg, pk)
println("  ✅ Julia verify: ", ret == 0 ? "SUCCESS" : "FAIL")

# =====================================================
# Test 2: Julia keypair → Julia sign → C verify
# =====================================================
println("\n[Test 2] Julia keypair → Julia sign → C verify")
c_sig = copy(sig)
c_msg = copy(msg)
c_pk  = copy(pk)
c_ret = ccall((:pqcrystals_dilithium2_ref_verify, C_LIB),
              Int32,
              (Ptr{UInt8}, UInt64, Ptr{UInt8}, UInt64, Ptr{UInt8}),
              c_sig, UInt64(length(sig)),
              c_msg, UInt64(length(msg)),
              c_pk)
println("  ✅ C verify: ", c_ret == 0 ? "SUCCESS" : "FAIL")

# =====================================================
# Test 3: 複数メッセージでの検証
# =====================================================
println("\n[Test 3] 複数メッセージでの検証")
test_messages = [
    "Short",
    "Medium length message",
    "A much longer message that spans multiple lines and contains various characters: 0123456789!@#\$%",
    "" # 空メッセージ
]

all_pass = true
for (i, tmsg) in enumerate(test_messages)
    msg_bytes = UInt8.(collect(tmsg))
    sig = crypto_sign_signature(msg_bytes, sk)
    ret = crypto_sign_verify(sig, msg_bytes, pk)
    status = ret == 0 ? "✅" : "❌"
    println("  $status Message $i ($(length(tmsg)) chars): $ret")
    if ret != 0
        all_pass = false
    end
end

# =====================================================
# Test 4: 性能ベンチマーク
# =====================================================
println("\n[Test 4] 性能ベンチマーク")
using Printf

# Keypair
t1 = time()
for _ in 1:10
    crypto_sign_keypair(seed)
end
t2 = time()
keypair_time = (t2 - t1) / 10 * 1000
@printf("  Keypair: %.2f ms/op\n", keypair_time)

# Sign
msg = UInt8.("Benchmark message" |> collect)
t1 = time()
for _ in 1:10
    crypto_sign_signature(msg, sk)
end
t2 = time()
sign_time = (t2 - t1) / 10 * 1000
@printf("  Sign:    %.2f ms/op\n", sign_time)

# Verify
t1 = time()
for _ in 1:10
    crypto_sign_verify(sig, msg, pk)
end
t2 = time()
verify_time = (t2 - t1) / 10 * 1000
@printf("  Verify:  %.2f ms/op\n", verify_time)

println("\n" * "="^60)
println("  🎉 全テスト完了！")
println("="^60)
