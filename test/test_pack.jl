include(joinpath(@__DIR__, "..", "src", "Dilithium2.jl"))
using .Dilithium2.Poly: Polynomial, PolyVec, polyt1_pack, polyt1_unpack
using .Dilithium2.Packing: pack_pk, unpack_pk

println("=== Test polyt1_pack/unpack roundtrip ===")
p1 = Polynomial()
for i in 1:256
    p1.coeffs[i] = Int32(mod(i, 1024))  # [0, 1023] (10-bit)
end

buf = zeros(UInt8, 320)
polyt1_pack(buf, p1)
println("Packed buf[1:10] = $(buf[1:10])")

p2 = Polynomial()
polyt1_unpack(p2, buf)
println("p1[1:10] = $(p1.coeffs[1:10])")
println("p2[1:10] = $(p2.coeffs[1:10])")
println("Roundtrip match: ", p1 == p2)

println("\n=== Test pack_pk/unpack_pk ===")
rho = zeros(UInt8, 32)
for i in 1:32
    rho[i] = UInt8(i)
end

t1 = PolyVec(4)
for k in 1:4
    for i in 1:256
        t1.polys[k].coeffs[i] = Int32(mod(k*100 + i, 1024))
    end
end

println("Original t1[1][1:10] = $(t1.polys[1].coeffs[1:10])")
println("Original t1[2][1:10] = $(t1.polys[2].coeffs[1:10])")

pk = zeros(UInt8, 1312)
pack_pk(pk, rho, t1)
println("\npk[1:10] = $(pk[1:10])")
println("pk[33:42] = $(pk[33:42])")  # t1開始位置

rho2 = zeros(UInt8, 32)
t1_2 = PolyVec(4)
unpack_pk(rho2, t1_2, pk)

println("\nUnpacked t1[1][1:10] = $(t1_2.polys[1].coeffs[1:10])")
println("Unpacked t1[2][1:10] = $(t1_2.polys[2].coeffs[1:10])")

println("\nrho match: ", rho == rho2)
println("t1[1] match: ", t1.polys[1].coeffs == t1_2.polys[1].coeffs)
println("t1[2] match: ", t1.polys[2].coeffs == t1_2.polys[2].coeffs)
println("t1[3] match: ", t1.polys[3].coeffs == t1_2.polys[3].coeffs)
println("t1[4] match: ", t1.polys[4].coeffs == t1_2.polys[4].coeffs)
