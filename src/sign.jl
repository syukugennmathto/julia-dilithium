"""
署名・検証の実装
sign.c に対応
"""

module Sign

using ..Params
using ..Fips202: shake256
using ..Poly: Polynomial, PolyVec, poly_ntt, poly_challenge
using ..PolyVecOps: polyvec_matrix_expand,
                 polyvecl_uniform_eta, polyveck_uniform_eta,
                 polyvecl_uniform_gamma1,
                 polyvecl_ntt, polyveck_ntt,
                 polyvecl_invntt_tomont, polyveck_invntt_tomont,
                 polyvec_matrix_pointwise_montgomery,
                 polyvecl_pointwise_poly_montgomery, polyveck_pointwise_poly_montgomery,
                 polyveck_reduce, polyvecl_reduce,
                 polyveck_caddq,
                 polyveck_add, polyveck_sub, polyvecl_add,
                 polyveck_shiftl,
                 polyveck_power2round, polyveck_decompose,
                 polyveck_make_hint, polyveck_use_hint,
                 polyvecl_chknorm, polyveck_chknorm,
                 polyveck_pack_w1
using ..Packing: pack_pk, unpack_pk, pack_sk, unpack_sk, pack_sig, unpack_sig

export crypto_sign_keypair, crypto_sign_signature, crypto_sign_verify

# crh: CRH(input) = SHAKE256(input, CRHBYTES)
function crh(input::Vector{UInt8})
    return shake256(CRHBYTES, input)
end

# PolyVec のディープコピー
function copy_polyvec(src::PolyVec)
    dst = PolyVec(src.k)
    for i in 1:src.k
        copyto!(dst.polys[i].coeffs, src.polys[i].coeffs)
    end
    return dst
end

"""
crypto_sign_keypair: 公開鍵・秘密鍵の生成
"""
function crypto_sign_keypair(seed::Vector{UInt8})
    # shake256(seedbuf, 3*SEEDBYTES, seedbuf, SEEDBYTES)
    seedbuf = shake256(3 * SEEDBYTES, seed)
    rho       = seedbuf[1:SEEDBYTES]
    rhoprime  = seedbuf[SEEDBYTES+1:2*SEEDBYTES]
    key       = seedbuf[2*SEEDBYTES+1:3*SEEDBYTES]

    # Expand matrix
    mat = [PolyVec(L) for _ in 1:K]
    polyvec_matrix_expand(mat, rho)

    # Sample s1, s2
    s1 = PolyVec(L)
    polyvecl_uniform_eta(s1, rhoprime, UInt16(0))
    s2 = PolyVec(K)
    polyveck_uniform_eta(s2, rhoprime, UInt16(L))

    # t1 = A*s1 + s2
    s1hat = copy_polyvec(s1)
    polyvecl_ntt(s1hat)

    t1 = PolyVec(K)
    polyvec_matrix_pointwise_montgomery(t1, mat, s1hat)
    polyveck_reduce(t1)
    polyveck_invntt_tomont(t1)
    polyveck_add(t1, t1, s2)

    # power2round
    polyveck_caddq(t1)
    t1_high = PolyVec(K)
    t0 = PolyVec(K)
    polyveck_power2round(t1_high, t0, t1)

    # pack public key
    pk = zeros(UInt8, CRYPTO_PUBLICKEYBYTES)
    pack_pk(pk, rho, t1_high)

    # tr = CRH(pk)
    tr = crh(pk)

    # pack secret key
    sk = zeros(UInt8, CRYPTO_SECRETKEYBYTES)
    pack_sk(sk, rho, tr, key, t0, s1, s2)

    return pk, sk
end

"""
crypto_sign_signature: メッセージの署名
"""
function crypto_sign_signature(m::Vector{UInt8}, sk::Vector{UInt8})
    # Unpack secret key
    rho      = zeros(UInt8, SEEDBYTES)
    tr       = zeros(UInt8, CRHBYTES)
    key      = zeros(UInt8, SEEDBYTES)
    t0       = PolyVec(K)
    s1       = PolyVec(L)
    s2       = PolyVec(K)
    unpack_sk(rho, tr, key, t0, s1, s2, sk)

    # mu = CRH(tr || msg)
    mu_input = [tr; m]
    mu = shake256(CRHBYTES, mu_input)

    # rhoprime = CRH(key || tr)  (deterministic signing)
    rhoprime_input = [key; tr]
    rhoprime = shake256(CRHBYTES, rhoprime_input)

    # Expand matrix
    mat = [PolyVec(L) for _ in 1:K]
    polyvec_matrix_expand(mat, rho)

    # NTT transform
    polyvecl_ntt(s1)
    polyveck_ntt(s2)
    polyveck_ntt(t0)

    nonce = UInt16(0)
    sig = zeros(UInt8, CRYPTO_BYTES)

    while true
        # Sample y
        y = PolyVec(L)
        polyvecl_uniform_gamma1(y, rhoprime, nonce)
        nonce += UInt16(1)

        # z = y (コピー)
        z = copy_polyvec(y)
        polyvecl_ntt(z)

        # w1 = A*z
        w1 = PolyVec(K)
        polyvec_matrix_pointwise_montgomery(w1, mat, z)
        polyveck_reduce(w1)
        polyveck_invntt_tomont(w1)

        # Decompose w
        polyveck_caddq(w1)
        w1_high = PolyVec(K)
        w0 = PolyVec(K)
        polyveck_decompose(w1_high, w0, w1)
        polyveck_pack_w1(sig, w1_high)

        # c = H(mu || w1)
        shake_input = [mu; sig[1:K*POLYW1_PACKEDBYTES]]
        c_seed = shake256(SEEDBYTES, shake_input)

        # challenge
        cp = Polynomial()
        poly_challenge(cp, c_seed)
        poly_ntt(cp)

        # z = c*s1 + y
        z = PolyVec(L)
        polyvecl_pointwise_poly_montgomery(z, cp, s1)
        polyvecl_invntt_tomont(z)
        polyvecl_add(z, z, y)
        polyvecl_reduce(z)

        if polyvecl_chknorm(z, Int32(GAMMA1 - BETA)) != 0
            continue
        end

        # w0 = w0 - c*s2
        h = PolyVec(K)
        polyveck_pointwise_poly_montgomery(h, cp, s2)
        polyveck_invntt_tomont(h)
        polyveck_sub(w0, w0, h)
        polyveck_reduce(w0)

        if polyveck_chknorm(w0, Int32(GAMMA2 - BETA)) != 0
            continue
        end

        # hints
        h = PolyVec(K)
        polyveck_pointwise_poly_montgomery(h, cp, t0)
        polyveck_invntt_tomont(h)
        polyveck_reduce(h)

        if polyveck_chknorm(h, Int32(GAMMA2)) != 0
            continue
        end

        polyveck_add(w0, w0, h)
        polyveck_caddq(w0)
        n = polyveck_make_hint(h, w0, w1_high)
        if n > OMEGA
            continue
        end

        # Pack signature
        pack_sig(sig, c_seed, z, h)
        return sig
    end
end

"""
crypto_sign_verify: 署名の検証
"""
function crypto_sign_verify(sig::Vector{UInt8}, m::Vector{UInt8}, pk::Vector{UInt8})
    if length(sig) != CRYPTO_BYTES
        return -1
    end

    # Unpack public key
    rho = zeros(UInt8, SEEDBYTES)
    t1  = PolyVec(K)
    unpack_pk(rho, t1, pk)

    # Unpack signature
    c  = zeros(UInt8, SEEDBYTES)
    z  = PolyVec(L)
    h  = PolyVec(K)
    ret = unpack_sig(c, z, h, sig)
    if ret != 0
        return -1
    end

    if polyvecl_chknorm(z, Int32(GAMMA1 - BETA)) != 0
        return -1
    end

    # mu = CRH(CRH(pk), msg)
    # C: crh(mu, pk, CRYPTO_PUBLICKEYBYTES);
    #    shake256_init(&state);
    #    shake256_absorb(&state, mu, CRHBYTES);
    #    shake256_absorb(&state, m, mlen);
    #    shake256_finalize(&state);
    #    shake256_squeeze(mu, CRHBYTES, &state);
    tr = crh(pk)  # tr = CRH(pk)
    shake_input = [tr; m]
    mu = shake256(CRHBYTES, shake_input)

    # challenge
    cp = Polynomial()
    poly_challenge(cp, c)

    # Expand matrix
    mat = [PolyVec(L) for _ in 1:K]
    polyvec_matrix_expand(mat, rho)

    # w1 = A*z - c*2^d*t1
    polyvecl_ntt(z)
    w1 = PolyVec(K)
    polyvec_matrix_pointwise_montgomery(w1, mat, z)

    poly_ntt(cp)
    polyveck_shiftl(t1)
    polyveck_ntt(t1)
    polyveck_pointwise_poly_montgomery(t1, cp, t1)

    polyveck_sub(w1, w1, t1)
    polyveck_reduce(w1)
    polyveck_invntt_tomont(w1)

    # Reconstruct w1
    polyveck_caddq(w1)
    polyveck_use_hint(w1, w1, h)

    buf = zeros(UInt8, K * POLYW1_PACKEDBYTES)
    polyveck_pack_w1(buf, w1)

    # c2 = H(mu || w1)
    shake_input2 = [mu; buf]
    c2 = shake256(SEEDBYTES, shake_input2)

    # Compare c and c2
    if c != c2
        return -1
    end

    return 0
end

end # module
