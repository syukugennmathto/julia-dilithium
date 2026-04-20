"""
鍵・署名のパッキング/アンパッキング
packing.c に対応
"""

module Packing

using ..Params
using ..Poly: Polynomial, PolyVec,
              polyeta_pack, polyeta_unpack,
              polyt1_pack, polyt1_unpack,
              polyt0_pack, polyt0_unpack,
              polyz_pack, polyz_unpack

export pack_pk, unpack_pk, pack_sk, unpack_sk, pack_sig, unpack_sig

"""
pack_pk: 公開鍵をバイト列にパック
pk = rho || t1
"""
function pack_pk(pk::Vector{UInt8}, rho::Vector{UInt8}, t1::PolyVec)
    pk[1:SEEDBYTES] = rho[1:SEEDBYTES]
    for i in 0:(K-1)
        offset = SEEDBYTES + i * POLYT1_PACKEDBYTES
        tmp = zeros(UInt8, POLYT1_PACKEDBYTES)
        polyt1_pack(tmp, t1.polys[i+1])
        pk[offset+1:offset+POLYT1_PACKEDBYTES] = tmp
    end
end

"""
unpack_pk: バイト列から公開鍵をアンパック
"""
function unpack_pk(rho::Vector{UInt8}, t1::PolyVec, pk::Vector{UInt8})
    rho[1:SEEDBYTES] = pk[1:SEEDBYTES]
    for i in 0:(K-1)
        offset = SEEDBYTES + i * POLYT1_PACKEDBYTES
        polyt1_unpack(t1.polys[i+1], pk[offset+1:offset+POLYT1_PACKEDBYTES])
    end
end

"""
pack_sk: 秘密鍵をバイト列にパック
sk = rho || tr || key || t0 || s1 || s2
"""
function pack_sk(sk::Vector{UInt8}, rho::Vector{UInt8}, tr::Vector{UInt8},
                 key::Vector{UInt8}, t0::PolyVec, s1::PolyVec, s2::PolyVec)
    idx = 1
    # rho
    sk[idx:idx+SEEDBYTES-1] = rho[1:SEEDBYTES]
    idx += SEEDBYTES
    # tr
    sk[idx:idx+CRHBYTES-1] = tr[1:CRHBYTES]
    idx += CRHBYTES
    # key
    sk[idx:idx+SEEDBYTES-1] = key[1:SEEDBYTES]
    idx += SEEDBYTES
    # t0
    for i in 0:(K-1)
        tmp = zeros(UInt8, POLYT0_PACKEDBYTES)
        polyt0_pack(tmp, t0.polys[i+1])
        sk[idx:idx+POLYT0_PACKEDBYTES-1] = tmp
        idx += POLYT0_PACKEDBYTES
    end
    # s1
    for i in 0:(L-1)
        tmp = zeros(UInt8, POLYETA_PACKEDBYTES)
        polyeta_pack(tmp, s1.polys[i+1])
        sk[idx:idx+POLYETA_PACKEDBYTES-1] = tmp
        idx += POLYETA_PACKEDBYTES
    end
    # s2
    for i in 0:(K-1)
        tmp = zeros(UInt8, POLYETA_PACKEDBYTES)
        polyeta_pack(tmp, s2.polys[i+1])
        sk[idx:idx+POLYETA_PACKEDBYTES-1] = tmp
        idx += POLYETA_PACKEDBYTES
    end
end

"""
unpack_sk: バイト列から秘密鍵をアンパック
"""
function unpack_sk(rho::Vector{UInt8}, tr::Vector{UInt8}, key::Vector{UInt8},
                   t0::PolyVec, s1::PolyVec, s2::PolyVec, sk::Vector{UInt8})
    idx = 1
    # rho
    rho[1:SEEDBYTES] = sk[idx:idx+SEEDBYTES-1]
    idx += SEEDBYTES
    # tr
    tr[1:CRHBYTES] = sk[idx:idx+CRHBYTES-1]
    idx += CRHBYTES
    # key
    key[1:SEEDBYTES] = sk[idx:idx+SEEDBYTES-1]
    idx += SEEDBYTES
    # t0
    for i in 0:(K-1)
        polyt0_unpack(t0.polys[i+1], sk[idx:idx+POLYT0_PACKEDBYTES-1])
        idx += POLYT0_PACKEDBYTES
    end
    # s1
    for i in 0:(L-1)
        polyeta_unpack(s1.polys[i+1], sk[idx:idx+POLYETA_PACKEDBYTES-1])
        idx += POLYETA_PACKEDBYTES
    end
    # s2
    for i in 0:(K-1)
        polyeta_unpack(s2.polys[i+1], sk[idx:idx+POLYETA_PACKEDBYTES-1])
        idx += POLYETA_PACKEDBYTES
    end
end

"""
pack_sig: 署名をバイト列にパック
sig = c || z || h
"""
function pack_sig(sig::Vector{UInt8}, c::Vector{UInt8}, z::PolyVec, h::PolyVec)
    idx = 1
    # c (SEEDBYTES)
    sig[idx:idx+SEEDBYTES-1] = c[1:SEEDBYTES]
    idx += SEEDBYTES
    # z
    for i in 0:(L-1)
        tmp = zeros(UInt8, POLYZ_PACKEDBYTES)
        polyz_pack(tmp, z.polys[i+1])
        sig[idx:idx+POLYZ_PACKEDBYTES-1] = tmp
        idx += POLYZ_PACKEDBYTES
    end
    # h のパック: 各多項式の1の位置と、各多項式のオフセット
    # h はバイト列として K バイトの末尾に非ゼロ位置を記録する
    hints = zeros(UInt8, OMEGA + K)
    k = 0
    for i in 0:(K-1)
        for j in 1:N
            if h.polys[i+1].coeffs[j] != 0
                hints[k+1] = UInt8(j - 1)  # 0-indexed
                k += 1
            end
        end
        hints[OMEGA + i + 1] = UInt8(k)
    end
    sig[idx:idx+OMEGA+K-1] = hints
end

"""
unpack_sig: バイト列から署名をアンパック
Returns 0 on success, -1 on failure
"""
function unpack_sig(c::Vector{UInt8}, z::PolyVec, h::PolyVec, sig::Vector{UInt8})
    idx = 1
    # c
    c[1:SEEDBYTES] = sig[idx:idx+SEEDBYTES-1]
    idx += SEEDBYTES
    # z
    for i in 0:(L-1)
        polyz_unpack(z.polys[i+1], sig[idx:idx+POLYZ_PACKEDBYTES-1])
        idx += POLYZ_PACKEDBYTES
    end
    # h のアンパック
    # 全てのhを0にしてから、ヒント位置を読む
    for i in 1:K
        for j in 1:N
            h.polys[i].coeffs[j] = Int32(0)
        end
    end

    # オフセット配列を読む
    offsets = zeros(Int, K + 1)
    offsets[1] = 0
    for i in 1:K
        offsets[i+1] = Int(sig[idx + OMEGA + i - 1])
    end

    # 検証: オフセットが単調増加で範囲内か
    for i in 1:K
        if offsets[i+1] < offsets[i] || offsets[i+1] > OMEGA
            return -1
        end
    end

    # ヒント位置を読んで設定する
    for i in 1:K
        for j in (offsets[i]+1):offsets[i+1]
            pos = Int(sig[idx + j - 1])
            if pos >= N || (j > offsets[i] + 1 && pos <= Int(sig[idx + j - 2]))
                return -1  # 位置が範囲外か単調非増加
            end
            h.polys[i].coeffs[pos + 1] = Int32(1)
        end
    end

    return 0
end

end # module
