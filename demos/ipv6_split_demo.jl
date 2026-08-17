# Demo: "split" IPv6 coding -- treat the /32 network prefix, the subnet id,
# and the interface identifier as three separate streams with separate
# models, instead of one per-byte-position model over all 16 bytes.
#
#   julia --project ipv6_split_demo.jl
#
# Rationale:
#   * prefix (bytes 1-4)  -- in real traffic only a handful of distinct
#     blocks occur. Coding each byte independently (as in ipv6_demo.jl)
#     under-uses the fact that the 4 bytes are *jointly* one of a few fixed
#     tuples: a per-byte model spends bits as if byte2 could vary freely
#     given byte1, when really the whole tuple is drawn from a tiny
#     dictionary. So we code it as ONE dictionary symbol (known block, or
#     ESCAPE + the 4 literal bytes for the long tail).
#   * subnet id (bytes 5-8) -- still has real per-byte structure (mostly
#     small numbers), so keep the per-position model.
#   * interface id (bytes 9-16) -- in our sample this is a mix of
#     structured (EUI-64, sequential host numbers) and fully random
#     (privacy extensions). We compare entropy-coding it (per-position
#     model) against just packing it raw, which is what real systems do
#     once a field is close enough to uniform that modeling isn't worth it.

using PadicCoding
using Random

Random.seed!(1)

const KNOWN_PREFIXES = [
    UInt8[0x26,0x06,0x47,0x00],
    UInt8[0x20,0x01,0x0d,0xb8],
    UInt8[0x24,0x00,0xcb,0x00],
    UInt8[0x2a,0x03,0x20,0x40],
    UInt8[0xfe,0x80,0x00,0x00],
]
const WEIGHTS = [0.20, 0.15, 0.15, 0.10, 0.15]   # remaining 0.25 -> long tail
const ESCAPE = length(KNOWN_PREFIXES) + 1

function random_iid()
    style = rand()
    if style < 0.15
        iid = zeros(UInt8, 8); iid[8] = rand(UInt8(1):UInt8(32)); iid
    elseif style < 0.40
        iid = rand(UInt8, 8); iid[4] = 0xff; iid[5] = 0xfe; iid
    else
        rand(UInt8, 8)
    end
end

function sample_ips6(n)
    prefix_sym = Vector{Int}(undef, n)          # 1..5 known, ESCAPE otherwise
    prefix_lit = Vector{Vector{UInt8}}(undef, n) # literal 4 bytes (always kept)
    subnet     = Vector{Vector{UInt8}}(undef, n)
    iid        = Vector{Vector{UInt8}}(undef, n)

    for i in 1:n
        r = rand()
        acc = 0.0
        sym = ESCAPE
        lit = UInt8[rand(0x20:0x3f), rand(UInt8), rand(UInt8), rand(UInt8)]
        for (k, w) in enumerate(WEIGHTS)
            acc += w
            if r <= acc
                sym = k
                lit = KNOWN_PREFIXES[k]
                break
            end
        end
        prefix_sym[i] = sym
        prefix_lit[i] = lit

        sn = zeros(UInt8, 4)
        if rand() < 0.7
            sn[4] = rand(UInt8(0):UInt8(16))
        else
            sn = rand(UInt8, 4)
        end
        subnet[i] = sn
        iid[i] = random_iid()
    end
    return prefix_sym, prefix_lit, subnet, iid
end

n = 20_000
prefix_sym, prefix_lit, subnet, iid = sample_ips6(n)
raw_bits = 128 * n

# helper: static PadicModel coded size for a column of bytes (1..256 symbols)
function code_bytes(col::Vector{UInt8})
    freqs = ones(Int, 256)
    for b in col; freqs[Int(b) + 1] += 1; end
    msg = Int.(col) .+ 1
    m = PadicModel(freqs; P = 2, N = 40)
    enc = padic_encode(m, msg)
    @assert padic_decode(m, enc) == msg
    return length(enc)
end

# --- 1) prefix: dictionary symbol (+ escaped literal bytes when needed) -----
sym_freqs = ones(Int, ESCAPE)   # 5 known + 1 escape
for s in prefix_sym; sym_freqs[s] += 1; end
msym = PadicModel(sym_freqs; P = 2, N = 24)
sym_enc = padic_encode(msym, prefix_sym)
@assert padic_decode(msym, sym_enc) == prefix_sym
prefix_dict_bits = length(sym_enc)

escape_idx = findall(==(ESCAPE), prefix_sym)
escape_literal_bits = length(escape_idx) * 32   # long tail: just pack raw

# for comparison: coding the prefix with the ipv6_demo.jl per-byte-position
# approach (4 independent byte models)
prefix_perbyte_bits = sum(code_bytes([prefix_lit[i][k] for i in 1:n]) for k in 1:4)

# --- 2) subnet id: per-byte-position model (kept, has real structure) ------
subnet_bits = sum(code_bytes([subnet[i][k] for i in 1:n]) for k in 1:4)

# --- 3) interface id: per-byte-position model vs raw packing ---------------
iid_coded_bits = sum(code_bytes([iid[i][k] for i in 1:n]) for k in 1:8)
iid_raw_bits = 8 * 8 * n  # 8 bytes/address, no coding

println("addresses: $n  (raw: $raw_bits bits = $(raw_bits ÷ 8) bytes)")
println()
println("prefix (bytes 1-4):")
println("  per-byte-position coding : $prefix_perbyte_bits bits ",
        "($(round(prefix_perbyte_bits/n, digits=3)) bits/address)")
println("  dictionary + escape      : $(prefix_dict_bits + escape_literal_bits) bits ",
        "($(round((prefix_dict_bits+escape_literal_bits)/n, digits=3)) bits/address) ",
        "[dict: $prefix_dict_bits, $(length(escape_idx)) escapes x 32b: $escape_literal_bits]")
println()
println("subnet id (bytes 5-8), per-byte-position: $subnet_bits bits ",
        "($(round(subnet_bits/n, digits=3)) bits/address)")
println()
println("interface id (bytes 9-16):")
println("  per-byte-position coding : $iid_coded_bits bits ",
        "($(round(iid_coded_bits/n, digits=3)) bits/address)")
println("  raw packed (no coding)   : $iid_raw_bits bits ",
        "($(round(iid_raw_bits/n, digits=3)) bits/address)")
println()

total_full_percolumn = prefix_perbyte_bits + subnet_bits + iid_coded_bits
total_split_coded    = prefix_dict_bits + escape_literal_bits + subnet_bits + iid_coded_bits
total_split_practical = prefix_dict_bits + escape_literal_bits + subnet_bits + iid_raw_bits

for (label, bits) in (
    ("all per-byte-position (ipv6_demo.jl baseline)", total_full_percolumn),
    ("split: dict prefix + coded subnet + coded iid", total_split_coded),
    ("split: dict prefix + coded subnet + RAW iid",   total_split_practical),
)
    println(rpad(label, 48), ": $bits bits = $(round(bits/8, digits=1)) bytes, ",
            "$(round(bits/n, digits=3)) bits/address, ",
            "ratio vs raw: $(round(raw_bits/bits, digits=2))x")
end
