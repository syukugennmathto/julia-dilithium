using Libdl
using Printf

include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))
using .Dilithium2.Sign: crypto_sign_keypair, crypto_sign_signature, crypto_sign_verify

const C_LIB = joinpath(@__DIR__, "..", "..", "libpqcrystals_dilithium2_ref.dylib")

println("="^70)
println("  Dilithium2 Performance Comparison: Julia vs C")
println("="^70)

# テスト用データ
seed = zeros(UInt8, 32)
for i in 1:32
    seed[i] = UInt8(i)
end

msg = UInt8.("Performance benchmark message for Dilithium2" |> collect)

# Julia版のkeypair（ウォームアップ＋事前生成）
println("\n[Warmup] Generating Julia keypair...")
pk_j, sk_j = crypto_sign_keypair(seed)
sig_j = crypto_sign_signature(msg, sk_j)
crypto_sign_verify(sig_j, msg, pk_j)

# C版のkeypair（事前生成）
println("[Warmup] Generating C keypair...")
pk_c = zeros(UInt8, 1312)
sk_c = zeros(UInt8, 2544)
ccall((:pqcrystals_dilithium2_ref_keypair, C_LIB),
      Int32, (Ptr{UInt8}, Ptr{UInt8}),
      pk_c, sk_c)

sig_c = zeros(UInt8, 2420)
siglen_c = Ref{UInt64}(0)
ccall((:pqcrystals_dilithium2_ref_signature, C_LIB),
      Int32, (Ptr{UInt8}, Ptr{UInt64}, Ptr{UInt8}, UInt64, Ptr{UInt8}),
      sig_c, siglen_c, msg, UInt64(length(msg)), sk_c)

ccall((:pqcrystals_dilithium2_ref_verify, C_LIB),
      Int32, (Ptr{UInt8}, UInt64, Ptr{UInt8}, UInt64, Ptr{UInt8}),
      sig_c, UInt64(2420), msg, UInt64(length(msg)), pk_c)

println("\nRunning benchmarks (100 iterations each)...\n")

# ==========================================
# Keypair Benchmark
# ==========================================
N = 100

# Julia Keypair
t_start = time()
for _ in 1:N
    crypto_sign_keypair(seed)
end
t_julia_keypair = (time() - t_start) / N * 1000

# C Keypair
t_start = time()
for _ in 1:N
    ccall((:pqcrystals_dilithium2_ref_keypair, C_LIB),
          Int32, (Ptr{UInt8}, Ptr{UInt8}),
          pk_c, sk_c)
end
t_c_keypair = (time() - t_start) / N * 1000

# ==========================================
# Sign Benchmark
# ==========================================

# Julia Sign
t_start = time()
for _ in 1:N
    crypto_sign_signature(msg, sk_j)
end
t_julia_sign = (time() - t_start) / N * 1000

# C Sign
t_start = time()
for _ in 1:N
    ccall((:pqcrystals_dilithium2_ref_signature, C_LIB),
          Int32, (Ptr{UInt8}, Ptr{UInt64}, Ptr{UInt8}, UInt64, Ptr{UInt8}),
          sig_c, siglen_c, msg, UInt64(length(msg)), sk_c)
end
t_c_sign = (time() - t_start) / N * 1000

# ==========================================
# Verify Benchmark
# ==========================================

# Julia Verify
t_start = time()
for _ in 1:N
    crypto_sign_verify(sig_j, msg, pk_j)
end
t_julia_verify = (time() - t_start) / N * 1000

# C Verify
t_start = time()
for _ in 1:N
    ccall((:pqcrystals_dilithium2_ref_verify, C_LIB),
          Int32, (Ptr{UInt8}, UInt64, Ptr{UInt8}, UInt64, Ptr{UInt8}),
          sig_c, UInt64(2420), msg, UInt64(length(msg)), pk_c)
end
t_c_verify = (time() - t_start) / N * 1000

# ==========================================
# Results
# ==========================================

println("┌" * "─"^68 * "┐")
println("│" * " "^20 * "Performance Results (ms/op)" * " "^21 * "│")
println("├" * "─"^68 * "┤")
@printf("│ %-20s │ %12s │ %12s │ %12s │\n", "Operation", "Julia", "C", "Ratio (C/J)")
println("├" * "─"^68 * "┤")
@printf("│ %-20s │ %12.3f │ %12.3f │ %12.2fx │\n", "Keypair", t_julia_keypair, t_c_keypair, t_c_keypair/t_julia_keypair)
@printf("│ %-20s │ %12.3f │ %12.3f │ %12.2fx │\n", "Sign", t_julia_sign, t_c_sign, t_c_sign/t_julia_sign)
@printf("│ %-20s │ %12.3f │ %12.3f │ %12.2fx │\n", "Verify", t_julia_verify, t_c_verify, t_c_verify/t_julia_verify)
println("└" * "─"^68 * "┘")

println("\nInterpretation:")
println("  Ratio < 1.0: Julia is faster")
println("  Ratio = 1.0: Same performance")
println("  Ratio > 1.0: C is faster")

# Overall summary
total_julia = t_julia_keypair + t_julia_sign + t_julia_verify
total_c = t_c_keypair + t_c_sign + t_c_verify

println("\n" * "="^70)
@printf("Total (Keypair + Sign + Verify): Julia=%.3f ms, C=%.3f ms\n", total_julia, total_c)
@printf("Overall speedup: %.2fx (C is %.1f%% of Julia speed)\n", 
        total_julia/total_c, (total_c/total_julia)*100)
println("="^70)
