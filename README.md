# CRYSTALS-Dilithium2 / ML-DSA (FIPS 204) — Pure Julia Implementation

[![Julia](https://img.shields.io/badge/Julia-1.9+-9558B2?style=flat&logo=julia)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status: Research](https://img.shields.io/badge/Status-Research%20%2F%20Educational-yellow)]()
[![Test Coverage](https://img.shields.io/badge/Coverage-96%25-brightgreen)]()

A pure Julia implementation of **CRYSTALS-Dilithium2**, the post-quantum digital signature scheme standardized as **ML-DSA** under **NIST FIPS 204** (August 13, 2024).

Achieves **bit-exact compatibility** with the official C reference implementation.  
Performance: **within 2.5× of native C** (total ~1.8 ms on commodity hardware).

---

## Background

### Why Post-Quantum Signatures?

Current digital signature schemes (RSA, ECDSA) rely on the hardness of integer factorization and discrete logarithm problems. Shor's algorithm, running on a sufficiently powerful quantum computer, solves both in polynomial time — rendering today's public-key infrastructure obsolete.

The threat is not only future: **"Harvest Now, Decrypt Later"** strategies mean adversaries may already be archiving encrypted traffic for future decryption. Migration to post-quantum cryptography is urgent.

### Why Dilithium / ML-DSA?

NIST's Post-Quantum Cryptography Standardization (2016–2024) selected four schemes:

| Standard | Algorithm | Purpose | Based On |
|---|---|---|---|
| FIPS 203 | ML-KEM | Key Encapsulation | CRYSTALS-Kyber |
| **FIPS 204** | **ML-DSA** | **Digital Signature** | **CRYSTALS-Dilithium** |
| FIPS 205 | SLH-DSA | Digital Signature | SPHINCS+ |
| FIPS 206 (draft) | FN-DSA | Digital Signature | FALCON |

ML-DSA (Dilithium) is NIST's primary recommendation for post-quantum digital signatures, suitable for most applications. Its security is based on the hardness of the **Module Learning With Errors (MLWE)** problem over polynomial rings — a lattice problem for which no efficient quantum algorithm is known.

> **Note:** This implementation targets CRYSTALS-Dilithium2, which corresponds to ML-DSA-44 (NIST Security Level 2). Full ML-DSA conformance (incorporating changes specified in FIPS 204 Appendix D) is planned for a future release.

### Why Julia?

Dilithium's core operations are **matrix-vector multiplications over polynomial rings**:

```
KeyGen:  t = A·s₁ + s₂     (A ∈ R^{k×l}, s₁ ∈ R^l, s₂ ∈ R^k)
Sign:    z = y + c·s₁       (y ∈ R^l, c ∈ R)
```

Julia is designed for exactly this class of computation. Its high-performance array operations, LLVM-backed JIT compilation, and tight BLAS/LAPACK integration make it structurally well-suited for lattice cryptography — despite being conventionally considered unsuitable for low-level cryptographic code. The gap between high-level array strength and low-level memory constraints is precisely what this implementation explores.

---

## Features

- ✅ Bit-exact compatibility with the [official C reference implementation](https://github.com/pq-crystals/dilithium)
- ✅ Key generation, signing, and verification
- ✅ Comprehensive test suite — 96% code coverage
- ✅ C interoperability tests (Julia signatures verified by C code)
- ✅ Benchmark suite included
- ✅ Fully documented implementation pitfalls (see [Implementation Notes](#implementation-notes))

---

## Quick Start

```bash
git clone https://github.com/syukugennmathto/julia-dilithium
cd julia-dilithium
julia -e 'include("src/Dilithium2.jl")'
```

```julia
using .Dilithium2.Sign

# Key generation
seed = rand(UInt8, 32)
pk, sk = crypto_sign_keypair(seed)

# Sign a message
msg = collect(UInt8, "Hello, Post-Quantum World!")
sig = crypto_sign_signature(msg, sk)
println("Signature size: ", length(sig), " bytes")  # => 2420 bytes

# Verify
result = crypto_sign_verify(sig, msg, pk)
println("Verification: ", result == 0 ? "OK" : "FAIL")  # => OK

# Run benchmarks
include("test/benchmark.jl")
```

---

## Performance

Benchmarked on commodity hardware (single-threaded, no hardware acceleration):

| Operation | Julia (ms) | C (ms) | Ratio |
|---|---|---|---|
| KeyGen | 0.413 | 0.184 | 2.24× |
| Sign | 0.907 | 0.317 | 2.86× |
| Verify | 0.476 | 0.232 | 2.05× |
| **Total** | **1.796** | **0.731** | **2.46×** |

**Remaining optimization headroom (~30–50% speedup):**

```julia
# Boundary check elimination
@inbounds for i in 1:N
    result[i] = a[i] + b[i]
end

# SIMD vectorization
@simd for i in 1:N
    c[i] = montgomery_reduce(...)
end
```

Expected post-optimization: **60–70% of native C speed**.

---

## Module Structure

```
src/
├── params.jl       # Constants and parameters
├── fips202.jl      # SHAKE128/256 (Keccak)
├── reduce.jl       # Modular reduction (Montgomery)
├── rounding.jl     # Power2Round, Decompose
├── ntt.jl          # Number Theoretic Transform
├── poly.jl         # Polynomial operations
├── polyvec.jl      # Vector operations
├── packing.jl      # Serialization / byte packing
└── sign.jl         # KeyGen, Sign, Verify
test/
├── test_primitives.jl
├── test_ntt.jl
├── test_packing.jl
├── test_sign.jl
├── test_c_interop.jl   # Verification against C reference
└── benchmark.jl
```

**Key parameters (Dilithium2 / ML-DSA-44):**

| Parameter | Value | Description |
|---|---|---|
| N | 256 | Polynomial degree |
| Q | 8,380,417 | Modulus |
| K, L | 4, 4 | Vector dimensions |
| ETA | 2 | Secret coefficient bound |
| GAMMA1 | 2¹⁷ = 131,072 | Signature coefficient bound |
| GAMMA2 | 95,232 | Low-order rounding |
| Public key | 1,312 bytes | |
| Secret key | 2,544 bytes | |
| Signature | 2,420 bytes | |

---

## Implementation Notes

Six significant challenges were encountered and resolved during development. These are documented here as a reference for future implementers.

### 1. Montgomery Reduction Constant (QINV)

**Symptom:** NTT output completely wrong (`Julia: -3148987`, `C: 287814`).  
**Cause:** Incorrect QINV constant (`20088143` instead of `58728449`).  
**Time:** ~2 hours.  
**Lesson:** Always validate constants against the reference source. Never assume.

```julia
const QINV = 58728449  # NOT 20088143
function montgomery_reduce(a::Int64)
    t = Int32(a * QINV)
    t = Int32((a - Int64(t) * Q) >> 32)
    return t
end
```

### 2. Array Index Offset (0-based vs 1-based)

**Symptom:** NTT output diverges from reference.  
**Cause:** C uses 0-based indexing; Julia uses 1-based. All `zetas[k]` accesses require `ZETAS[k+1]` in Julia.  
**Time:** ~3 hours.

### 3. Integer Overflow Semantics

**Symptom:** Incorrect results in modular arithmetic.  
**Cause:** C silently wraps Int32 overflow; Julia throws or wraps differently by type. Explicit modular arithmetic required throughout.  
**Time:** ~1 hour.

### 4. SubArray vs Vector in Packing Functions

**Symptom:** Public key `t1` field silently written as all-zeros. No error thrown.  
**Cause:** Julia's type system distinguishes `SubArray` and `Vector`. Passing a `SubArray` to packing functions silently fails.  
**Fix:** Use an intermediate `Vector` buffer.  
**Time:** ~4 hours (hardest to diagnose).

```julia
function pack_pk(pk::Vector{UInt8}, rho::Vector{UInt8}, t1::PolyVec)
    for i in 0:(K-1)
        tmp = zeros(UInt8, POLYT1_PACKEDBYTES)  # intermediate buffer
        polyt1_pack(tmp, t1.polys[i+1])
        pk[offset+1:offset+POLYT1_PACKEDBYTES] = tmp
    end
end
```

### 5. CRHBYTES Constant Mismatch

**Symptom:** `CRYPTO_SECRETKEYBYTES = 2560` (Julia) vs `2544` (C). C rejects Julia signatures.  
**Cause:** `CRHBYTES = 64` instead of correct value `48`. 16-byte difference propagates to secret key size.  
**Time:** ~2 hours.

### 6. Stateful SHAKE Emulation

**Symptom:** Rejection sampling loop fails intermittently.  
**Cause:** C uses a stateful Keccak context (squeeze-on-demand). Julia's SHA3 library provides a stateless API. Stateful behavior must be emulated with a buffer.  
**Time:** ~3 hours.

**Total debug time: ~15 hours (30% of total development time of ~50 hours)**  
**Most effective debugging tool: C interoperability tests (caught 65% of bugs)**

---

## Relationship to FIPS 204 (ML-DSA)

This implementation targets **CRYSTALS-Dilithium2** (the competition submission), which differs from the finalized **ML-DSA-44** specified in FIPS 204 (August 13, 2024).

Key differences documented in FIPS 204 Appendix D include:

- Input to `SampleInBall` changed from first 256 bits to the full commitment hash
- `ExpandMask` modified to take output bits from the beginning rather than at an offset
- Hint unpacking algorithm updated with restored malformed-input check

**These changes break interoperability.** An existing Dilithium2 implementation cannot be used as a drop-in ML-DSA-44 replacement. Full ML-DSA-44 conformance is a planned future milestone.

---

## IC Card Implementation Note

A known challenge for embedded / IC card deployment: the NTT constant table size.

ML-DSA uses modulus **Q = 8,380,417**, requiring NTT precomputation tables proportional to Q. This is significantly larger than ML-KEM's modulus (**Q = 3,329**), creating a concrete bottleneck for resource-constrained devices. Optimizing NTT table compression for ML-DSA on IC cards is an open engineering problem. Feedback and contributions targeting this use case are welcome.

---

## Limitations

This implementation is suitable for:
- ✅ Education and research
- ✅ Prototyping and experimentation
- ✅ Security-noncritical applications

This implementation is **not** suitable for:
- ❌ Production cryptographic systems
- ❌ Security-critical applications
- ❌ Embedded / resource-constrained devices (as-is)

**Known limitations:**
- Not constant-time (timing channel vulnerabilities exist)
- No side-channel protections (power analysis, EM analysis)
- Single-threaded
- No hardware acceleration

---

## Roadmap

- [ ] ML-DSA-44 full conformance (FIPS 204 Appendix D changes)
- [ ] ML-DSA-65 and ML-DSA-87 parameter sets
- [ ] `@inbounds` / `@simd` performance optimization
- [ ] Constant-time operations
- [ ] Julia package registry submission
- [ ] Masking countermeasures (side-channel resistance) — research in progress

---

## References

1. Ducas, L., et al. (2021). *CRYSTALS-Dilithium Algorithm Specifications v3.1*. https://pq-crystals.org/dilithium/
2. NIST (2024). *FIPS 204: Module-Lattice-Based Digital Signature Standard (ML-DSA)*. https://doi.org/10.6028/NIST.FIPS.204
3. NIST (2024). *FIPS 203: Module-Lattice-Based Key-Encapsulation Mechanism Standard (ML-KEM)*. https://doi.org/10.6028/NIST.FIPS.203
4. NIST (2024). *FIPS 205: Stateless Hash-Based Digital Signature Standard (SLH-DSA)*. https://doi.org/10.6028/NIST.FIPS.205
5. pq-crystals/dilithium — Official C Reference Implementation. https://github.com/pq-crystals/dilithium
6. Bezanson, J., et al. (2017). *Julia: A Fresh Approach to Numerical Computing*. SIAM Review, 59(1), 65–98.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

This implementation is provided for educational and research purposes.  
**Not recommended for production use without constant-time and side-channel hardening.**

---

*Presented at Quantum Computing EXPO Spring 2026 (Tokyo Big Sight, April 15, 2026) — OSS value confirmed in discussion with TOPPAN Holdings.*
