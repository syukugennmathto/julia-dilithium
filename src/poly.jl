"""
多項式演算の実装
poly.c に対応
"""

module Poly

using ..Params
using ..NTT: montgomery_reduce, ntt!, invntt!
using ..Reduce: reduce32, caddq, freeze
using ..Rounding: power2round, decompose, make_hint, use_hint
using ..Fips202: shake128, shake256

export Polynomial, PolyVec
export poly_reduce, poly_caddq, poly_freeze
export poly_add, poly_sub, poly_shiftl
export poly_ntt, poly_invntt_tomont, poly_pointwise_montgomery
export poly_power2round, poly_decompose, poly_make_hint, poly_use_hint
export poly_chknorm
export rej_uniform, rej_eta
export polyeta_pack, polyeta_unpack
export polyt1_pack, polyt1_unpack
export polyt0_pack, polyt0_unpack
export polyz_pack, polyz_unpack
export polyw1_pack

# =====================================================
# データ構造
# =====================================================

mutable struct Polynomial
    coeffs::Vector{Int32}

    function Polynomial()
        new(zeros(Int32, N))
    end

    function Polynomial(coeffs::Vector{Int32})
        @assert length(coeffs) == N
        new(copy(coeffs))
    end
end

struct PolyVec
    polys::Vector{Polynomial}
    k::Int

    function PolyVec(k::Int)
        new([Polynomial() for _ in 1:k], k)
    end
end

Base.copy(p::Polynomial) = Polynomial(copy(p.coeffs))
Base.:(==)(p1::Polynomial, p2::Polynomial) = p1.coeffs == p2.coeffs

# =====================================================
# 基本演算
# =====================================================

function poly_reduce(a::Polynomial)
    for i in 1:N
        a.coeffs[i] = reduce32(a.coeffs[i])
    end
end

function poly_caddq(a::Polynomial)
    for i in 1:N
        a.coeffs[i] = caddq(a.coeffs[i])
    end
end

function poly_freeze(a::Polynomial)
    for i in 1:N
        a.coeffs[i] = freeze(a.coeffs[i])
    end
end

function poly_add(c::Polynomial, a::Polynomial, b::Polynomial)
    for i in 1:N
        c.coeffs[i] = a.coeffs[i] + b.coeffs[i]
    end
end

function poly_sub(c::Polynomial, a::Polynomial, b::Polynomial)
    for i in 1:N
        c.coeffs[i] = a.coeffs[i] - b.coeffs[i]
    end
end

function poly_shiftl(a::Polynomial)
    for i in 1:N
        a.coeffs[i] = a.coeffs[i] << D
    end
end

# =====================================================
# NTT演算
# =====================================================

function poly_ntt(a::Polynomial)
    ntt!(a.coeffs)
end

function poly_invntt_tomont(a::Polynomial)
    invntt!(a.coeffs)
end

function poly_pointwise_montgomery(c::Polynomial, a::Polynomial, b::Polynomial)
    for i in 1:N
        c.coeffs[i] = montgomery_reduce(Int64(a.coeffs[i]) * Int64(b.coeffs[i]))
    end
end

# =====================================================
# Rounding
# =====================================================

function poly_power2round(a1::Polynomial, a0::Polynomial, a::Polynomial)
    for i in 1:N
        a1.coeffs[i], a0.coeffs[i] = power2round(a.coeffs[i])
    end
end

function poly_decompose(a1::Polynomial, a0::Polynomial, a::Polynomial)
    for i in 1:N
        a1.coeffs[i], a0.coeffs[i] = decompose(a.coeffs[i])
    end
end

function poly_make_hint(h::Polynomial, a0::Polynomial, a1::Polynomial)
    s = Int32(0)
    for i in 1:N
        h.coeffs[i] = make_hint(a0.coeffs[i], a1.coeffs[i])
        s += h.coeffs[i]
    end
    return s
end

function poly_use_hint(b::Polynomial, a::Polynomial, h::Polynomial)
    for i in 1:N
        b.coeffs[i] = use_hint(a.coeffs[i], h.coeffs[i])
    end
end

# =====================================================
# Norm check
# =====================================================

function poly_chknorm(a::Polynomial, B::Int32)
    if B > Int32((Q - 1) ÷ 8)
        return 1
    end
    for i in 1:N
        t = a.coeffs[i] >> 31
        t = a.coeffs[i] - (t & (2 * a.coeffs[i]))
        if t >= B
            return 1
        end
    end
    return 0
end

# =====================================================
# Rejection sampling
# =====================================================

"""
rej_uniform: [0, Q-1]の一様乱数サンプリング
C実装と完全に同じ動作
"""
function rej_uniform(a::Vector{Int32}, offset::Int, len::Int, buf::Vector{UInt8}, buflen::Int)
    ctr = 0
    pos = 0

    while ctr < len && pos + 3 <= buflen
        t  = UInt32(buf[pos + 1])
        t |= UInt32(buf[pos + 2]) << 8
        t |= UInt32(buf[pos + 3]) << 16
        t &= UInt32(0x7FFFFF)
        pos += 3

        if t < UInt32(Q)
            a[offset + ctr + 1] = Int32(t)
            ctr += 1
        end
    end

    return ctr
end

"""
rej_eta: [-ETA, ETA]の一様乱数サンプリング (ETA=2)
C実装と完全に同じ動作
"""
function rej_eta(a::Vector{Int32}, offset::Int, len::Int, buf::Vector{UInt8}, buflen::Int)
    ctr = 0
    pos = 0

    while ctr < len && pos < buflen
        t0 = UInt32(buf[pos + 1] & 0x0F)
        t1 = UInt32(buf[pos + 1] >> 4)
        pos += 1

        # ETA == 2
        if t0 < 15
            t0 = t0 - ((205 * t0) >> 10) * 5
            a[offset + ctr + 1] = Int32(2) - Int32(t0)
            ctr += 1
        end
        if t1 < 15 && ctr < len
            t1 = t1 - ((205 * t1) >> 10) * 5
            a[offset + ctr + 1] = Int32(2) - Int32(t1)
            ctr += 1
        end
    end

    return ctr
end

# =====================================================
# Packing (ETA=2, GAMMA1=2^17, GAMMA2=(Q-1)/88)
# =====================================================

"""
polyeta_pack: [-ETA,ETA] の多項式をバイト列にパック (ETA=2)
"""
function polyeta_pack(r::Vector{UInt8}, a::Polynomial)
    for i in 0:(N÷8 - 1)
        t = zeros(UInt8, 8)
        for j in 0:7
            t[j+1] = UInt8(ETA - a.coeffs[8*i + j + 1])
        end

        r[3*i+1]  = (t[1] >> 0) | (t[2] << 3) | (t[3] << 6)
        r[3*i+2]  = (t[3] >> 2) | (t[4] << 1) | (t[5] << 4) | (t[6] << 7)
        r[3*i+3]  = (t[6] >> 1) | (t[7] << 2) | (t[8] << 5)
    end
end

"""
polyeta_unpack: バイト列から [-ETA,ETA] 多項式をアンパック (ETA=2)
"""
function polyeta_unpack(r::Polynomial, a::Vector{UInt8})
    for i in 0:(N÷8 - 1)
        r.coeffs[8*i+1] =  Int32((a[3*i+1] >> 0) & 7)
        r.coeffs[8*i+2] =  Int32((a[3*i+1] >> 3) & 7)
        r.coeffs[8*i+3] =  Int32(((UInt16(a[3*i+1]) >> 6) | (UInt16(a[3*i+2]) << 2)) & 7)
        r.coeffs[8*i+4] =  Int32((a[3*i+2] >> 1) & 7)
        r.coeffs[8*i+5] =  Int32((a[3*i+2] >> 4) & 7)
        r.coeffs[8*i+6] =  Int32(((UInt16(a[3*i+2]) >> 7) | (UInt16(a[3*i+3]) << 1)) & 7)
        r.coeffs[8*i+7] =  Int32((a[3*i+3] >> 2) & 7)
        r.coeffs[8*i+8] =  Int32((a[3*i+3] >> 5) & 7)

        for j in 1:8
            r.coeffs[8*i+j] = Int32(ETA) - r.coeffs[8*i+j]
        end
    end
end

"""
polyt1_pack: t1を10ビットパック
"""
function polyt1_pack(r::Vector{UInt8}, a::Polynomial)
    for i in 0:(N÷4 - 1)
        r[5*i+1] = UInt8((a.coeffs[4*i+1] >> 0) & 0xFF)
        r[5*i+2] = UInt8(((a.coeffs[4*i+1] >> 8) | (a.coeffs[4*i+2] << 2)) & 0xFF)
        r[5*i+3] = UInt8(((a.coeffs[4*i+2] >> 6) | (a.coeffs[4*i+3] << 4)) & 0xFF)
        r[5*i+4] = UInt8(((a.coeffs[4*i+3] >> 4) | (a.coeffs[4*i+4] << 6)) & 0xFF)
        r[5*i+5] = UInt8((a.coeffs[4*i+4] >> 2) & 0xFF)
    end
end

"""
polyt1_unpack: 10ビットパックからt1をアンパック
"""
function polyt1_unpack(r::Polynomial, a::Vector{UInt8})
    for i in 0:(N÷4 - 1)
        r.coeffs[4*i+1] = Int32(((UInt32(a[5*i+1]) >> 0) | (UInt32(a[5*i+2]) << 8)) & 0x3FF)
        r.coeffs[4*i+2] = Int32(((UInt32(a[5*i+2]) >> 2) | (UInt32(a[5*i+3]) << 6)) & 0x3FF)
        r.coeffs[4*i+3] = Int32(((UInt32(a[5*i+3]) >> 4) | (UInt32(a[5*i+4]) << 4)) & 0x3FF)
        r.coeffs[4*i+4] = Int32(((UInt32(a[5*i+4]) >> 6) | (UInt32(a[5*i+5]) << 2)) & 0x3FF)
    end
end

"""
polyt0_pack: t0を13ビットパック
"""
function polyt0_pack(r::Vector{UInt8}, a::Polynomial)
    for i in 0:(N÷8 - 1)
        t = zeros(UInt32, 8)
        for j in 0:7
            t[j+1] = UInt32((1 << (D-1)) - a.coeffs[8*i + j + 1])
        end

        r[13*i+ 1]  =  UInt8(t[1] & 0xFF)
        r[13*i+ 2]  =  UInt8((t[1] >>  8) & 0xFF)
        r[13*i+ 2] |=  UInt8((t[2] <<  5) & 0xFF)
        r[13*i+ 3]  =  UInt8((t[2] >>  3) & 0xFF)
        r[13*i+ 4]  =  UInt8((t[2] >> 11) & 0xFF)
        r[13*i+ 4] |=  UInt8((t[3] <<  2) & 0xFF)
        r[13*i+ 5]  =  UInt8((t[3] >>  6) & 0xFF)
        r[13*i+ 5] |=  UInt8((t[4] <<  7) & 0xFF)
        r[13*i+ 6]  =  UInt8((t[4] >>  1) & 0xFF)
        r[13*i+ 7]  =  UInt8((t[4] >>  9) & 0xFF)
        r[13*i+ 7] |=  UInt8((t[5] <<  4) & 0xFF)
        r[13*i+ 8]  =  UInt8((t[5] >>  4) & 0xFF)
        r[13*i+ 9]  =  UInt8((t[5] >> 12) & 0xFF)
        r[13*i+ 9] |=  UInt8((t[6] <<  1) & 0xFF)
        r[13*i+10]  =  UInt8((t[6] >>  7) & 0xFF)
        r[13*i+10] |=  UInt8((t[7] <<  6) & 0xFF)
        r[13*i+11]  =  UInt8((t[7] >>  2) & 0xFF)
        r[13*i+12]  =  UInt8((t[7] >> 10) & 0xFF)
        r[13*i+12] |=  UInt8((t[8] <<  3) & 0xFF)
        r[13*i+13]  =  UInt8((t[8] >>  5) & 0xFF)
    end
end

"""
polyt0_unpack: 13ビットパックからt0をアンパック
"""
function polyt0_unpack(r::Polynomial, a::Vector{UInt8})
    for i in 0:(N÷8 - 1)
        r.coeffs[8*i+1]  = Int32(UInt32(a[13*i+ 1]))
        r.coeffs[8*i+1] |= Int32(UInt32(a[13*i+ 2]) << 8)
        r.coeffs[8*i+1]  &= Int32(0x1FFF)

        r.coeffs[8*i+2]  = Int32(UInt32(a[13*i+ 2]) >> 5)
        r.coeffs[8*i+2] |= Int32(UInt32(a[13*i+ 3]) << 3)
        r.coeffs[8*i+2] |= Int32(UInt32(a[13*i+ 4]) << 11)
        r.coeffs[8*i+2]  &= Int32(0x1FFF)

        r.coeffs[8*i+3]  = Int32(UInt32(a[13*i+ 4]) >> 2)
        r.coeffs[8*i+3] |= Int32(UInt32(a[13*i+ 5]) << 6)
        r.coeffs[8*i+3]  &= Int32(0x1FFF)

        r.coeffs[8*i+4]  = Int32(UInt32(a[13*i+ 5]) >> 7)
        r.coeffs[8*i+4] |= Int32(UInt32(a[13*i+ 6]) << 1)
        r.coeffs[8*i+4] |= Int32(UInt32(a[13*i+ 7]) << 9)
        r.coeffs[8*i+4]  &= Int32(0x1FFF)

        r.coeffs[8*i+5]  = Int32(UInt32(a[13*i+ 7]) >> 4)
        r.coeffs[8*i+5] |= Int32(UInt32(a[13*i+ 8]) << 4)
        r.coeffs[8*i+5] |= Int32(UInt32(a[13*i+ 9]) << 12)
        r.coeffs[8*i+5]  &= Int32(0x1FFF)

        r.coeffs[8*i+6]  = Int32(UInt32(a[13*i+ 9]) >> 1)
        r.coeffs[8*i+6] |= Int32(UInt32(a[13*i+10]) << 7)
        r.coeffs[8*i+6]  &= Int32(0x1FFF)

        r.coeffs[8*i+7]  = Int32(UInt32(a[13*i+10]) >> 6)
        r.coeffs[8*i+7] |= Int32(UInt32(a[13*i+11]) << 2)
        r.coeffs[8*i+7] |= Int32(UInt32(a[13*i+12]) << 10)
        r.coeffs[8*i+7]  &= Int32(0x1FFF)

        r.coeffs[8*i+8]  = Int32(UInt32(a[13*i+12]) >> 3)
        r.coeffs[8*i+8] |= Int32(UInt32(a[13*i+13]) << 5)
        r.coeffs[8*i+8]  &= Int32(0x1FFF)

        for j in 1:8
            r.coeffs[8*i+j] = Int32(1 << (D-1)) - r.coeffs[8*i+j]
        end
    end
end

"""
polyz_pack: z を18ビットパック (GAMMA1 = 2^17)
"""
function polyz_pack(r::Vector{UInt8}, a::Polynomial)
    for i in 0:(N÷4 - 1)
        t = zeros(UInt32, 4)
        for j in 0:3
            t[j+1] = UInt32(GAMMA1 - a.coeffs[4*i + j + 1])
        end

        r[9*i+1]  = UInt8(t[1] & 0xFF)
        r[9*i+2]  = UInt8((t[1] >> 8) & 0xFF)
        r[9*i+3]  = UInt8((t[1] >> 16) & 0xFF)
        r[9*i+3] |= UInt8((t[2] << 2) & 0xFF)
        r[9*i+4]  = UInt8((t[2] >> 6) & 0xFF)
        r[9*i+5]  = UInt8((t[2] >> 14) & 0xFF)
        r[9*i+5] |= UInt8((t[3] << 4) & 0xFF)
        r[9*i+6]  = UInt8((t[3] >> 4) & 0xFF)
        r[9*i+7]  = UInt8((t[3] >> 12) & 0xFF)
        r[9*i+7] |= UInt8((t[4] << 6) & 0xFF)
        r[9*i+8]  = UInt8((t[4] >> 2) & 0xFF)
        r[9*i+9]  = UInt8((t[4] >> 10) & 0xFF)
    end
end

"""
polyz_unpack: 18ビットパックからz をアンパック (GAMMA1 = 2^17)
"""
function polyz_unpack(r::Polynomial, a::Vector{UInt8})
    for i in 0:(N÷4 - 1)
        r.coeffs[4*i+1]  = Int32(UInt32(a[9*i+1]))
        r.coeffs[4*i+1] |= Int32(UInt32(a[9*i+2]) << 8)
        r.coeffs[4*i+1] |= Int32(UInt32(a[9*i+3]) << 16)
        r.coeffs[4*i+1]  &= Int32(0x3FFFF)

        r.coeffs[4*i+2]  = Int32(UInt32(a[9*i+3]) >> 2)
        r.coeffs[4*i+2] |= Int32(UInt32(a[9*i+4]) << 6)
        r.coeffs[4*i+2] |= Int32(UInt32(a[9*i+5]) << 14)
        r.coeffs[4*i+2]  &= Int32(0x3FFFF)

        r.coeffs[4*i+3]  = Int32(UInt32(a[9*i+5]) >> 4)
        r.coeffs[4*i+3] |= Int32(UInt32(a[9*i+6]) << 4)
        r.coeffs[4*i+3] |= Int32(UInt32(a[9*i+7]) << 12)
        r.coeffs[4*i+3]  &= Int32(0x3FFFF)

        r.coeffs[4*i+4]  = Int32(UInt32(a[9*i+7]) >> 6)
        r.coeffs[4*i+4] |= Int32(UInt32(a[9*i+8]) << 2)
        r.coeffs[4*i+4] |= Int32(UInt32(a[9*i+9]) << 10)
        r.coeffs[4*i+4]  &= Int32(0x3FFFF)

        for j in 1:4
            r.coeffs[4*i+j] = Int32(GAMMA1) - r.coeffs[4*i+j]
        end
    end
end

"""
polyw1_pack: w1を6ビットパック (GAMMA2 = (Q-1)/88)
"""
function polyw1_pack(r::Vector{UInt8}, a::Polynomial)
    for i in 0:(N÷4 - 1)
        r[3*i+1]  = UInt8(a.coeffs[4*i+1] & 0xFF)
        r[3*i+1] |= UInt8((a.coeffs[4*i+2] << 6) & 0xFF)
        r[3*i+2]  = UInt8((a.coeffs[4*i+2] >> 2) & 0xFF)
        r[3*i+2] |= UInt8((a.coeffs[4*i+3] << 4) & 0xFF)
        r[3*i+3]  = UInt8((a.coeffs[4*i+3] >> 4) & 0xFF)
        r[3*i+3] |= UInt8((a.coeffs[4*i+4] << 2) & 0xFF)
    end
end


# =====================================================
# poly_challenge (sign.c で使用)
# SHAKE256(seed) から TAU=39 個の非ゼロ係数を生成
# =====================================================

export poly_challenge

function poly_challenge(c::Polynomial, seed::Vector{UInt8})
    # C実装と同じ: SHAKE256のブロックサイズ = 136
    SHAKE256_RATE = 136

    # SHAKE256(seed) を1ブロック分スクイーズ
    buf = shake256(SHAKE256_RATE, seed)

    # 最初の8バイトから signs を読む (LE)
    signs = UInt64(0)
    for i in 0:7
        signs |= UInt64(buf[i+1]) << (8*i)
    end
    pos = 8  # 次に読むバイト位置 (0-indexed)

    # 全て0にする
    for i in 1:N
        c.coeffs[i] = Int32(0)
    end

    # C: for(i = N-TAU; i < N; ++i)
    # Fisher-Yates の変形: i番目の位置をランダムに選んで交換
    for i in (N - TAU):(N - 1)  # 0-indexed: 217..255
        # b を i 以下の値になるまで棄却サンプリング
        b = 0
        while true
            if pos >= SHAKE256_RATE
                # 追加スクイーズが必要（実際には殺害にはならないが念のため）
                buf = shake256(SHAKE256_RATE, seed)
                pos = 0
            end
            b = Int(buf[pos + 1])
            pos += 1
            if b <= i
                break
            end
        end

        # c[i] = c[b]; c[b] = 1 - 2*(signs & 1)
        # Julia: 1-indexed なので +1
        c.coeffs[i + 1] = c.coeffs[b + 1]
        c.coeffs[b + 1] = Int32(1) - Int32(2) * Int32(signs & UInt64(1))
        signs >>= 1
    end
end
end # module
