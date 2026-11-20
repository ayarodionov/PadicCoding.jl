# Experiment: compressed data on a PAM-3 link (3 voltage levels per symbol,
# as in GDDR7, USB4 v2 and 100BASE-T1 Ethernet).
#
#   julia --project demos/pam3_demo.jl [file ...]
#
# A PAM-3 symbol carries at most log2(3) = 1.585 bits.  Binary data has to be
# converted to ternary symbols first; links do it with small block codes such
# as 3 bits -> 2 symbols (9 states for 8 values, 1.5 bits/symbol).  Here the
# compressed output of the P = 2 coder is converted with several block codes
# and with one whole-message base conversion, and compared with the P = 3 coder,
# which writes ternary symbols directly.  Every path is decoded and checked.

using PadicCoding, Printf

# ---- k bits -> m trits block codes (2^k <= 3^m) -----------------------------
function block_encode(bits, k, m)
    bits = vcat(bits, zeros(Int, mod(-length(bits), k)))   # trailing 0s are harmless
    out = Int[]
    for i in 1:k:length(bits)
        v = foldl((a, b) -> 2a + b, bits[i:i+k-1]; init = 0)
        append!(out, reverse(digits(v; base = 3, pad = m)))
    end
    return out
end

function block_decode(trits, k, m)
    out = Int[]
    for i in 1:m:length(trits)
        v = foldl((a, b) -> 3a + b, trits[i:i+m-1]; init = 0)
        append!(out, reverse(digits(v; base = 2, pad = k)))
    end
    return out
end

# ---- whole-message conversion: bits as one big number written in base 3 ----
# T = smallest trit count with 3^T >= 2^L.  The receiver reads back
# Lmax = floor(T log2 3) >= L bits; the extra trailing zeros decode the same.
# Exact integer comparisons, starting from a floating-point estimate.
function trits_for_bits(L)
    T = max(0, ceil(Int, L / log2(3)) - 2)
    p3, p2 = big(3)^T, big(2)^L
    while p3 < p2; p3 *= 3; T += 1; end
    return T
end

function bits_for_trits(T)
    L = max(0, floor(Int, T * log2(3)) - 2)
    p2, p3 = big(2)^(L + 1), big(3)^T
    while p2 <= p3; p2 *= 2; L += 1; end
    return L
end

function bigint_encode(bits)
    L = length(bits)
    T = trits_for_bits(L)
    Lmax = bits_for_trits(T)
    x = parse(BigInt, join(vcat(bits, zeros(Int, Lmax - L))); base = 2)
    return reverse(digits(x; base = 3, pad = T))
end

function bigint_decode(trits)
    T = length(trits)
    Lmax = bits_for_trits(T)
    x = parse(BigInt, join(trits); base = 3)
    return reverse(digits(x; base = 2, pad = Lmax))
end

# ---- experiment --------------------------------------------------------------
function experiment(name, data::Vector{UInt8})
    msg = Int.(data) .+ 1
    n = length(msg)
    freqs = ones(Int, 256)
    for s in msg; freqs[s] += 1; end
    ideal = -sum(log2(freqs[s] / sum(freqs)) for s in msg)

    m2 = PadicModel(freqs; P = 2, N = 40)
    m3 = PadicModel(freqs; P = 3, N = 26)
    bits = padic_encode(m2, msg)
    rows = []

    raw = reduce(vcat, [reverse(digits(b; base = 2, pad = 8)) for b in data])
    push!(rows, ("raw bytes, 3b->2s (no compression)", length(block_encode(raw, 3, 2)), NaN))

    for (k, m) in ((3, 2), (11, 7), (19, 12))
        t = @elapsed sym = block_encode(bits, k, m)
        @assert padic_decode(m2, block_decode(sym, k, m)) == msg
        push!(rows, (@sprintf("P=2 coder, %db->%ds block (%.1f%%)", k, m,
                              100 * k / (m * log2(3))), length(sym), t))
    end

    t = @elapsed sym = bigint_encode(bits)
    @assert padic_decode(m2, bigint_decode(sym)) == msg
    push!(rows, ("P=2 coder, whole-message base conversion", length(sym), t))

    t = @elapsed sym = padic_encode(m3, msg)
    @assert padic_decode(m3, sym) == msg
    push!(rows, ("P=3 coder, direct", length(sym), t))

    floor3 = ceil(Int, ideal / log2(3))
    @printf("\n%s: %d bytes, %.4f bits/byte, P=2 coder output %d bits\n",
            name, n, ideal / n, length(bits))
    @printf("  floor (entropy / log2 3): %d symbols\n", floor3)
    @printf("  %-44s %9s %9s %9s %10s\n", "path", "symbols", "sym/byte", "vs P=3", "convert s")
    ref = rows[end][2]
    for (lbl, s, t) in rows
        @printf("  %-44s %9d %9.4f %+8.2f%% %10s\n", lbl, s, s / n,
                100 * (s / ref - 1), isnan(t) ? "-" : @sprintf("%.3f", t))
    end
end

function corpus()
    root = joinpath(@__DIR__, "..")
    files = vcat(joinpath(root, "README.md"),
                 [joinpath(root, d, f) for d in ("src", "demos")
                  for f in sort(readdir(joinpath(root, d))) if endswith(f, ".jl")])
    out = [("project text + Julia source", reduce(vcat, read.(files)))]
    words = "/usr/share/dict/words"
    isfile(words) && push!(out, ("/usr/share/dict/words, first 150 KB",
                                 read(words)[1:150_000]))
    return out
end

sources = isempty(ARGS) ? corpus() : [(f, read(f)) for f in ARGS]
for (name, data) in sources
    experiment(name, data)
end
println("\nconvert s: time of the binary->ternary step alone, except for the P=3",
        "\ncoder, where it is the whole compression.")
