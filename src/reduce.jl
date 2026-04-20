"""
モジュラー reduction の実装
reduce.c に対応
"""

module Reduce

using ..Params

export reduce32, caddq, freeze

"""
reduce32: C実装と完全に同じ
  t = (a + (1 << 22)) >> 23;
  t = a - t*Q;
"""
function reduce32(a::Int32)
    t = (a + Int32(1 << 22)) >> 23
    t = a - t * Int32(Q)
    return t
end

"""
caddq: a が負なら Q を加える
  a += (a >> 31) & Q;
"""
function caddq(a::Int32)
    a += (a >> 31) & Int32(Q)
    return a
end

"""
freeze: 標準代表元 [0, Q-1] に縮約
"""
function freeze(a::Int32)
    a = reduce32(a)
    a = caddq(a)
    return a
end

end # module
