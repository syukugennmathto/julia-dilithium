"""
多項式ベクトル演算の実装
polyvec.c に対応
"""

module PolyVecOps

using ..Params
using ..Poly: Polynomial
using ..Poly: PolyVec
using ..Poly:
              poly_reduce, poly_caddq, poly_freeze,
              poly_add, poly_sub, poly_shiftl,
              poly_ntt, poly_invntt_tomont, poly_pointwise_montgomery,
              poly_power2round, poly_decompose, poly_make_hint, poly_use_hint,
              poly_chknorm, rej_uniform, rej_eta,
              polyw1_pack, polyz_unpack
using ..Fips202: shake128, shake256

export polyvec_matrix_expand
export polyvecl_uniform_eta, polyveck_uniform_eta
export polyvecl_uniform_gamma1
export polyvecl_ntt, polyveck_ntt
export polyvecl_invntt_tomont, polyveck_invntt_tomont
export polyvec_matrix_pointwise_montgomery
export polyvecl_pointwise_poly_montgomery, polyveck_pointwise_poly_montgomery
export polyveck_reduce, polyvecl_reduce
export polyveck_caddq
export polyveck_add, polyveck_sub, polyvecl_add
export polyveck_shiftl
export polyveck_power2round, polyveck_decompose
export polyveck_make_hint, polyveck_use_hint
export polyvecl_chknorm, polyveck_chknorm
export polyveck_pack_w1

const STREAM128_BLOCKBYTES = 168
const STREAM256_BLOCKBYTES = 136

# =====================================================
# poly_uniform: SHAKE128ベースの一様サンプリング
# C: stream128_init(&state, seed, nonce) → seed||nonce(LE 2byte)をSHAKE128に吸込む
# =====================================================

function poly_uniform(a::Polynomial, seed::Vector{UInt8}, nonce::UInt16)
    POLY_UNIFORM_NBLOCKS = (768 + STREAM128_BLOCKBYTES - 1) ÷ STREAM128_BLOCKBYTES

    # input = seed || nonce(LE)
    input = Vector{UInt8}(undef, length(seed) + 2)
    copyto!(input, seed)
    input[end-1] = UInt8(nonce & 0xFF)
    input[end]   = UInt8((nonce >> 8) & 0xFF)

    buflen = POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES
    buf = shake128(buflen + 3, input)  # 余分に3バイト

    ctr = rej_uniform(a.coeffs, 0, N, buf, buflen)

    while ctr < N
        off = buflen % 3
        # 末尾の余り部分を先頭にコピー + 追加スクイーズ
        # 実際には一回目で全て収まるが念のため
        new_input = copy(input)
        new_buf = shake128(STREAM128_BLOCKBYTES + off + 3, new_input)
        combined = Vector{UInt8}(undef, STREAM128_BLOCKBYTES + off)
        copyto!(combined, 1, buf, buflen - off + 1, off)
        copyto!(combined, off + 1, new_buf, 1, STREAM128_BLOCKBYTES)
        buflen = STREAM128_BLOCKBYTES + off
        ctr += rej_uniform(a.coeffs, ctr, N - ctr, combined, buflen)
        buf = combined
    end
end

# =====================================================
# polyvec_matrix_expand: A = expandA(rho)
# C: nonce = (i << 8) + j  ← 重要！i*L+jではない
# =====================================================

function polyvec_matrix_expand(mat::Vector{PolyVec}, rho::Vector{UInt8})
    for i in 0:(K-1)
        for j in 0:(L-1)
            nonce = UInt16((i << 8) + j)  # C実装と同じ
            poly_uniform(mat[i+1].polys[j+1], rho, nonce)
        end
    end
end

# =====================================================
# Uniform eta sampling
# =====================================================

function poly_uniform_eta(a::Polynomial, seed::Vector{UInt8}, nonce::UInt16)
    POLY_UNIFORM_ETA_NBLOCKS = (136 + STREAM128_BLOCKBYTES - 1) ÷ STREAM128_BLOCKBYTES

    input = Vector{UInt8}(undef, length(seed) + 2)
    copyto!(input, seed)
    input[end-1] = UInt8(nonce & 0xFF)
    input[end]   = UInt8((nonce >> 8) & 0xFF)

    buflen = POLY_UNIFORM_ETA_NBLOCKS * STREAM128_BLOCKBYTES
    buf = shake128(buflen, input)

    ctr = rej_eta(a.coeffs, 0, N, buf, buflen)

    while ctr < N
        buf = shake128(STREAM128_BLOCKBYTES, input)
        ctr += rej_eta(a.coeffs, ctr, N - ctr, buf, STREAM128_BLOCKBYTES)
    end
end

function polyvecl_uniform_eta(v::PolyVec, seed::Vector{UInt8}, nonce::UInt16)
    # C: nonce++ → nonce, nonce+1, nonce+2, ...
    for i in 0:(L-1)
        poly_uniform_eta(v.polys[i+1], seed, UInt16(nonce + i))
    end
end

function polyveck_uniform_eta(v::PolyVec, seed::Vector{UInt8}, nonce::UInt16)
    for i in 0:(K-1)
        poly_uniform_eta(v.polys[i+1], seed, UInt16(nonce + i))
    end
end

# =====================================================
# Uniform gamma1 sampling
# C: poly_uniform_gamma1(&v->vec[i], seed, L*nonce + i);
# =====================================================

function poly_uniform_gamma1(a::Polynomial, seed::Vector{UInt8}, nonce::UInt16)
    POLY_UNIFORM_GAMMA1_NBLOCKS = (576 + STREAM256_BLOCKBYTES - 1) ÷ STREAM256_BLOCKBYTES

    input = Vector{UInt8}(undef, length(seed) + 2)
    copyto!(input, seed)
    input[end-1] = UInt8(nonce & 0xFF)
    input[end]   = UInt8((nonce >> 8) & 0xFF)

    buflen = POLY_UNIFORM_GAMMA1_NBLOCKS * STREAM256_BLOCKBYTES
    buf = shake256(buflen, input)

    polyz_unpack(a, buf)
end

function polyvecl_uniform_gamma1(v::PolyVec, seed::Vector{UInt8}, nonce::UInt16)
    # C: L*nonce + i
    for i in 0:(L-1)
        poly_uniform_gamma1(v.polys[i+1], seed, UInt16(L * nonce + i))
    end
end

# =====================================================
# NTT
# =====================================================

function polyvecl_ntt(v::PolyVec)
    for i in 1:L
        poly_ntt(v.polys[i])
    end
end

function polyveck_ntt(v::PolyVec)
    for i in 1:K
        poly_ntt(v.polys[i])
    end
end

function polyvecl_invntt_tomont(v::PolyVec)
    for i in 1:L
        poly_invntt_tomont(v.polys[i])
    end
end

function polyveck_invntt_tomont(v::PolyVec)
    for i in 1:K
        poly_invntt_tomont(v.polys[i])
    end
end

# =====================================================
# Matrix-vector multiplication
# polyvecl_pointwise_acc_montgomery を内部で使用
# =====================================================

function polyvecl_pointwise_acc_montgomery(w::Polynomial, u::PolyVec, v::PolyVec)
    tmp = Polynomial()
    poly_pointwise_montgomery(w, u.polys[1], v.polys[1])
    for i in 2:L
        poly_pointwise_montgomery(tmp, u.polys[i], v.polys[i])
        poly_add(w, w, tmp)
    end
end

function polyvec_matrix_pointwise_montgomery(t::PolyVec, mat::Vector{PolyVec}, v::PolyVec)
    for i in 1:K
        polyvecl_pointwise_acc_montgomery(t.polys[i], mat[i], v)
    end
end

function polyvecl_pointwise_poly_montgomery(r::PolyVec, a::Polynomial, v::PolyVec)
    for i in 1:L
        poly_pointwise_montgomery(r.polys[i], a, v.polys[i])
    end
end

function polyveck_pointwise_poly_montgomery(r::PolyVec, a::Polynomial, v::PolyVec)
    for i in 1:K
        poly_pointwise_montgomery(r.polys[i], a, v.polys[i])
    end
end

# =====================================================
# Reduce / caddq
# =====================================================

function polyveck_reduce(v::PolyVec)
    for i in 1:K
        poly_reduce(v.polys[i])
    end
end

function polyvecl_reduce(v::PolyVec)
    for i in 1:L
        poly_reduce(v.polys[i])
    end
end

function polyveck_caddq(v::PolyVec)
    for i in 1:K
        poly_caddq(v.polys[i])
    end
end

# =====================================================
# Add / Sub
# =====================================================

function polyveck_add(r::PolyVec, a::PolyVec, b::PolyVec)
    for i in 1:K
        poly_add(r.polys[i], a.polys[i], b.polys[i])
    end
end

function polyvecl_add(r::PolyVec, a::PolyVec, b::PolyVec)
    for i in 1:L
        poly_add(r.polys[i], a.polys[i], b.polys[i])
    end
end

function polyveck_sub(r::PolyVec, a::PolyVec, b::PolyVec)
    for i in 1:K
        poly_sub(r.polys[i], a.polys[i], b.polys[i])
    end
end

# =====================================================
# Shiftl
# =====================================================

function polyveck_shiftl(v::PolyVec)
    for i in 1:K
        poly_shiftl(v.polys[i])
    end
end

# =====================================================
# Rounding
# =====================================================

function polyveck_power2round(t1::PolyVec, t0::PolyVec, t::PolyVec)
    for i in 1:K
        poly_power2round(t1.polys[i], t0.polys[i], t.polys[i])
    end
end

function polyveck_decompose(a1::PolyVec, a0::PolyVec, a::PolyVec)
    for i in 1:K
        poly_decompose(a1.polys[i], a0.polys[i], a.polys[i])
    end
end

function polyveck_make_hint(h::PolyVec, a0::PolyVec, a1::PolyVec)
    s = 0
    for i in 1:K
        s += poly_make_hint(h.polys[i], a0.polys[i], a1.polys[i])
    end
    return s
end

function polyveck_use_hint(b::PolyVec, a::PolyVec, h::PolyVec)
    for i in 1:K
        poly_use_hint(b.polys[i], a.polys[i], h.polys[i])
    end
end

# =====================================================
# Norm check
# =====================================================

function polyvecl_chknorm(v::PolyVec, B::Int32)
    for i in 1:L
        if poly_chknorm(v.polys[i], B) != 0
            return 1
        end
    end
    return 0
end

function polyveck_chknorm(v::PolyVec, B::Int32)
    for i in 1:K
        if poly_chknorm(v.polys[i], B) != 0
            return 1
        end
    end
    return 0
end

# =====================================================
# Pack w1
# =====================================================

function polyveck_pack_w1(r::Vector{UInt8}, v::PolyVec)
    for i in 1:K
        offset = (i-1) * POLYW1_PACKEDBYTES
        tmp = zeros(UInt8, POLYW1_PACKEDBYTES)
        polyw1_pack(tmp, v.polys[i])
        r[offset+1:offset+POLYW1_PACKEDBYTES] = tmp
    end
end

end # module PolyVecOps
