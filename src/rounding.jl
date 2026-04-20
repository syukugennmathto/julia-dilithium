"""
丸め演算の実装
rounding.c に対応 (GAMMA2 = (Q-1)/88 用)
"""

module Rounding

using ..Params

export power2round, decompose, make_hint, use_hint

"""
power2round: C実装と完全に同じ
"""
function power2round(a::Int32)
    a1 = (a + Int32((1 << (D-1)) - 1)) >> D
    a0 = a - (a1 << D)
    return a1, a0
end

"""
decompose: C実装と完全に同じ (GAMMA2 == (Q-1)/88)
  a1  = (a + 127) >> 7;
  a1  = (a1*11275 + (1 << 23)) >> 24;
  a1 ^= ((43 - a1) >> 31) & a1;
  *a0  = a - a1*2*GAMMA2;
  *a0 -= (((Q-1)/2 - *a0) >> 31) & Q;
"""
function decompose(a::Int32)
    a1 = (a + Int32(127)) >> 7
    a1 = Int32((Int64(a1) * 11275 + (1 << 23)) >> 24)
    # xor: a1 ^= ((43 - a1) >> 31) & a1
    a1 = xor(a1, ((Int32(43) - a1) >> 31) & a1)

    a0 = a - a1 * Int32(2 * GAMMA2)
    a0 -= (((Int32((Q-1)÷2) - a0) >> 31) & Int32(Q))
    return a1, a0
end

"""
make_hint: C実装と完全に同じ
"""
function make_hint(a0::Int32, a1::Int32)
    if a0 <= Int32(GAMMA2) || a0 > Int32(Q - GAMMA2) || (a0 == Int32(Q - GAMMA2) && a1 == Int32(0))
        return Int32(0)
    end
    return Int32(1)
end

"""
use_hint: C実装と完全に同じ (GAMMA2 == (Q-1)/88)
"""
function use_hint(a::Int32, hint::Int32)
    a1, a0 = decompose(a)
    if hint == Int32(0)
        return a1
    end
    if a0 > Int32(0)
        return (a1 == Int32(43)) ? Int32(0) : a1 + Int32(1)
    else
        return (a1 == Int32(0)) ? Int32(43) : a1 - Int32(1)
    end
end

end # module
