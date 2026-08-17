# Demo: compress a list of IPv4 addresses with the p-adic arithmetic coder.
#
#   julia --project ip_demo.jl
#
# Addresses are 4 octets, each treated as a symbol 1..256 (byte value + 1).
# Two static models are compared:
#   * flat   -- one byte-frequency model shared by all 4 octet positions
#               (same approach as demo.jl uses for text)
#   * by-pos -- one model per octet position (position 1 is heavily skewed
#               towards a handful of allocated blocks, position 4 is close
#               to uniform -- a flat model can't see that difference)
# against the 32-bit raw baseline.

using PadicCoding
using Random

Random.seed!(1)

# --- synthetic but realistic address list -----------------------------------
# Mimics a server access log / netflow sample: mostly private-range and a
# handful of common cloud/CDN /16 blocks, plus a long tail of random public IPs.
function sample_ips(n)
    blocks = [
        (10, 0, 0.35),        # 10.0.0.0/8   (private)
        (192, 168, 0.20),     # 192.168.0.0/16 (private)
        (172, 16, 0.10),      # 172.16.0.0/12 (private, approx as /16 here)
        (52, 84, 0.08),       # a cloud/CDN-ish /16
        (104, 16, 0.05),      # another CDN-ish /16
    ]
    ips = Vector{NTuple{4,UInt8}}(undef, n)
    for i in 1:n
        r = rand()
        acc = 0.0
        chosen = nothing
        for (a, b, w) in blocks
            acc += w
            if r <= acc
                chosen = (a, b)
                break
            end
        end
        if chosen === nothing
            # long tail: uniformly random public-looking address
            a = rand(1:223)
            b = rand(0:255)
            ips[i] = (UInt8(a), UInt8(b), rand(UInt8), rand(UInt8))
        else
            a, b = chosen
            ips[i] = (UInt8(a), UInt8(b), rand(UInt8), rand(UInt8))
        end
    end
    return ips
end

n = 20_000
ips = sample_ips(n)

raw_bits = 32 * n

# --- flat model: one shared byte model for all 4 positions ------------------
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

# --- by-position model: one model per octet position -------------------------
pos_freqs = [ones(Int, 256) for _ in 1:4]
for ip in ips, k in 1:4
    pos_freqs[k][Int(ip[k]) + 1] += 1
end
pos_msgs = [Int[] for _ in 1:4]
for ip in ips, k in 1:4
    push!(pos_msgs[k], Int(ip[k]) + 1)
end

pos_entropy_total = 0.0
pos_bits_total = 0
for k in 1:4
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
