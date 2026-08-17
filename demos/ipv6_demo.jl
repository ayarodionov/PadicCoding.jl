# Demo: compress a list of IPv6 addresses with the p-adic arithmetic coder.
#
#   julia --project ipv6_demo.jl
#
# An IPv6 address is 16 bytes. As in ip_demo.jl each byte is a symbol 1..256,
# and we compare:
#   * flat   -- one byte-frequency model shared by all 16 positions
#   * by-pos -- one model per byte position (bytes 0-7, the /64 network
#               prefix, are extremely skewed towards a handful of allocated
#               blocks; bytes 8-15, the interface identifier, range from
#               highly structured -- ::1, low sequential host numbers,
#               EUI-64 -- to fully random privacy-extension addresses)
# against the raw 128-bit baseline.

using PadicCoding
using Random

Random.seed!(1)

# --- synthetic but realistic address list -----------------------------------
# Mimics a mixed server/client log: a handful of common /32 or /48 global
# unicast prefixes (a hosting provider, a CDN, a couple of end-user /48s),
# plus link-local traffic, plus loopback -- with interface IDs that are
# either small sequential host numbers, EUI-64-shaped, or fully random
# (privacy extensions).
function sample_ips6(n)
    prefixes = [
        (UInt8[0x26,0x06,0x47,0x00], 0.20),   # 2606:4700::/32  (CDN-ish)
        (UInt8[0x20,0x01,0x0d,0xb8], 0.15),   # 2001:db8::/32   (docs/example)
        (UInt8[0x24,0x00,0xcb,0x00], 0.15),   # 2400:cb00::/32  (CDN-ish)
        (UInt8[0x2a,0x03,0x20,0x40], 0.10),   # a /32 hosting block
        (UInt8[0xfe,0x80,0x00,0x00], 0.15),   # fe80::/10 link-local
    ]

    function random_iid()
        style = rand()
        if style < 0.15
            # small sequential host number: ...::1 .. ::20
            iid = zeros(UInt8, 8)
            iid[8] = rand(UInt8(1):UInt8(32))
            return iid
        elseif style < 0.40
            # EUI-64 derived: xx:xx:xx:ff:fe:xx:xx:xx
            iid = rand(UInt8, 8)
            iid[4] = 0xff
            iid[5] = 0xfe
            return iid
        else
            # privacy-extension style: fully random
            return rand(UInt8, 8)
        end
    end

    ips = Vector{Vector{UInt8}}(undef, n)
    for i in 1:n
        r = rand()
        acc = 0.0
        net = UInt8[0x00,0x00,0x00,0x00]
        for (p, w) in prefixes
            acc += w
            if r <= acc
                net = p
                break
            end
        end
        if r > acc
            # long tail: random global unicast /32
            net = UInt8[rand(0x20:0x3f), rand(UInt8), rand(UInt8), rand(UInt8)]
        end
        # bytes 4-7: subnet id within the /32 (mostly small / structured)
        subnet = zeros(UInt8, 4)
        if rand() < 0.7
            subnet[4] = rand(UInt8(0):UInt8(16))   # small subnet numbers dominate
        else
            subnet = rand(UInt8, 4)
        end
        ips[i] = vcat(net, subnet, random_iid())
    end
    return ips
end

n = 20_000
ips = sample_ips6(n)

raw_bits = 128 * n
nbytes = 16

# --- flat model: one shared byte model for all 16 positions -----------------
flat_freqs = ones(Int, 256)
for ip in ips, byte in ip
    flat_freqs[Int(byte) + 1] += 1
end
flat_msg = Int[]
for ip in ips, byte in ip
    push!(flat_msg, Int(byte) + 1)
end
flat_entropy = let f = flat_freqs .- 1, tot = sum(f)
    p = filter(>(0), f ./ tot)
    -tot * sum(x -> x * log2(x), p)
end

mflat = PadicModel(flat_freqs; P = 2, N = 40)
flat_encoded = padic_encode(mflat, flat_msg)
flat_bits = length(flat_encoded)
@assert padic_decode(mflat, flat_encoded) == flat_msg

# --- by-position model: one model per byte position --------------------------
pos_freqs = [ones(Int, 256) for _ in 1:nbytes]
for ip in ips, k in 1:nbytes
    pos_freqs[k][Int(ip[k]) + 1] += 1
end
pos_msgs = [Int[] for _ in 1:nbytes]
for ip in ips, k in 1:nbytes
    push!(pos_msgs[k], Int(ip[k]) + 1)
end

pos_entropy_total = 0.0
pos_bits_total = 0
for k in 1:nbytes
    f = pos_freqs[k] .- 1
    tot = sum(f)
    p = filter(>(0), f ./ tot)
    global pos_entropy_total += -tot * sum(x -> x * log2(x), p)

    m = PadicModel(pos_freqs[k]; P = 2, N = 40)
    bits = padic_encode(m, pos_msgs[k])
    @assert padic_decode(m, bits) == pos_msgs[k]
    global pos_bits_total += length(bits)
end

println("addresses: $n  (raw: $(raw_bits) bits = $(raw_bits ÷ 8) bytes)")
println()
println("flat byte model:")
println("  empirical entropy: $(round(flat_entropy, digits=1)) bits total, ",
        "$(round(flat_entropy / n, digits=3)) bits/address")
println("  coded: $flat_bits bits = $(round(flat_bits/8, digits=1)) bytes, ",
        "$(round(flat_bits/n, digits=3)) bits/address, ",
        "ratio vs raw: $(round(raw_bits/flat_bits, digits=2))x")
println()
println("per-position model:")
println("  empirical entropy: $(round(pos_entropy_total, digits=1)) bits total, ",
        "$(round(pos_entropy_total / n, digits=3)) bits/address")
println("  coded: $pos_bits_total bits = $(round(pos_bits_total/8, digits=1)) bytes, ",
        "$(round(pos_bits_total/n, digits=3)) bits/address, ",
        "ratio vs raw: $(round(raw_bits/pos_bits_total, digits=2))x")
println()
println("per-byte-position entropy breakdown (bits, byte 0 = network prefix start):")
for k in 1:nbytes
    f = pos_freqs[k] .- 1
    tot = sum(f)
    p = filter(>(0), f ./ tot)
    e = -sum(x -> x * log2(x), p)
    label = k <= 8 ? "prefix[$k]" : "iid[$(k-8)]"
    println("  byte $k ($label): $(round(e, digits=3)) bits/byte")
end
