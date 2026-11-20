# Experiment: DNA storage with the "rotating" ternary code (Goldman et al.,
# Nature 494, 2013).
#
#   julia --project demos/dna_demo.jl [file ...]
#
# DNA synthesis and sequencing are error prone on homopolymers (runs of the
# same base).  Goldman et al. avoid them completely: data is first written as
# trits, and each trit selects one of the 3 bases that differ from the
# previous base.  So the channel is exactly ternary: log2(3) = 1.585 bits per
# nucleotide.  They produce the trits with a ternary Huffman code over bytes.
#
# Here the same static byte model drives
#   (a) a ternary Huffman code (Goldman),
#   (b) the p-adic arithmetic coder with P = 3 (output trits -> rotating code),
#   (c) the P = 2 coder with the plain 2-bits-per-base mapping (A,C,G,T), which
#       is denser but has homopolymers.
# Every path is decoded back from the DNA string and checked.
#
# The second part cuts the data into oligos with a fixed payload of L nt, as in
# real DNA storage (each oligo is synthesized and read separately), and counts
# how many oligos each code needs:
#   (a) ternary Huffman: codewords until the next one does not fit; the rest
#       is padded with a prefix of that codeword, so the decoder stops there;
#   (b) P = 3 coder with EOM, the longest prefix that fits;
#   (c) P = 3 fixed-length coder (padic_fill, no EOM) after a small header
#       with the number of bytes in the oligo.

using PadicCoding, Printf

const BASES = "ACGT"

# ---- rotating code: trits <-> DNA without homopolymers ----------------------
function trits_to_dna(trits)
    prev = 'A'                                   # virtual base before the start
    io = IOBuffer()
    for t in trits
        prev = [b for b in BASES if b != prev][t + 1]
        print(io, prev)
    end
    return String(take!(io))
end

function dna_to_trits(dna)
    prev = 'A'
    out = Int[]
    for b in dna
        push!(out, findfirst(==(b), [c for c in BASES if c != prev]) - 1)
        prev = b
    end
    return out
end

bits_to_dna(bits) = String([BASES[2bits[i] + bits[i+1] + 1]
                            for i in 1:2:length(bits)])
dna_to_bits(dna) = reduce(vcat, [[(i - 1) >> 1, (i - 1) & 1]
                                 for i in (findfirst(==(b), BASES) for b in dna)])

function homopolymers(dna)
    maxrun, run, long = 1, 1, 0
    for i in 2:length(dna)
        run = dna[i] == dna[i-1] ? run + 1 : 1
        run == 4 && (long += 1)                  # runs of length >= 4
        maxrun = max(maxrun, run)
    end
    return maxrun, long
end

# ---- D-ary Huffman code with codewords --------------------------------------
function huffman_code(p::Vector{Float64}, D::Int)
    n = length(p)
    pad = (D - 1 - (n - 1) % (D - 1)) % (D - 1)  # dummy leaves: full D-ary tree
    prob = vcat(p, zeros(pad))
    kids = [Int[] for _ in eachindex(prob)]
    live = collect(eachindex(prob))
    while length(live) > 1
        sort!(live; by = i -> prob[i])
        k = live[1:D]
        live = live[D+1:end]
        push!(prob, sum(prob[k])); push!(kids, k); push!(live, length(prob))
    end
    code = [Int[] for _ in 1:n]
    stack = [(live[1], Int[])]
    while !isempty(stack)
        v, pre = pop!(stack)
        if isempty(kids[v])
            v <= n && (code[v] = pre)
        else
            for (d, c) in enumerate(kids[v])
                push!(stack, (c, vcat(pre, d - 1)))
            end
        end
    end
    return code
end

function huffman_decode(code, digits)
    table = Dict(c => s for (s, c) in enumerate(code))
    out, cur = Int[], Int[]
    for d in digits
        push!(cur, d)
        s = get(table, cur, 0)
        s > 0 && (push!(out, s); empty!(cur))
    end
    return out
end

# ---- experiment --------------------------------------------------------------
function experiment(name, data::Vector{UInt8})
    msg = Int.(data) .+ 1
    n = length(msg)
    freqs = ones(Int, 256)
    for s in msg; freqs[s] += 1; end
    p = freqs ./ sum(freqs)
    ideal = -sum(log2(p[s]) for s in msg)        # bits, static model

    rows = []

    # (a) Goldman: ternary Huffman -> rotating code
    code = huffman_code(p, 3)
    trits = reduce(vcat, code[msg])
    dna = trits_to_dna(trits)
    @assert huffman_decode(code, dna_to_trits(dna)) == msg
    push!(rows, ("ternary Huffman + rotating", dna))

    # (b) p-adic P = 3 -> rotating code
    m3 = PadicModel(freqs; P = 3, N = 26)
    dna = trits_to_dna(padic_encode(m3, msg))
    @assert padic_decode(m3, dna_to_trits(dna)) == msg
    push!(rows, ("p-adic P=3 + rotating", dna))

    # (c) p-adic P = 2 -> 2 bits per base (homopolymers allowed)
    m2 = PadicModel(freqs; P = 2, N = 40)
    bits = padic_encode(m2, msg)
    isodd(length(bits)) && push!(bits, 0)        # trailing 0s decode the same
    dna = bits_to_dna(bits)
    @assert padic_decode(m2, dna_to_bits(dna)) == msg
    push!(rows, ("p-adic P=2 + 2 bits/base", dna))

    @printf("\n%s: %d bytes, static-model entropy %.4f bits/byte\n", name, n, ideal / n)
    @printf("  floor with no homopolymers: %d nt (%.4f nt/byte)\n",
            ceil(Int, ideal / log2(3)), ideal / log2(3) / n)
    @printf("  %-28s %9s %8s %9s %8s %9s\n",
            "method", "nt", "nt/byte", "vs P=3", "max run", "runs>=4")
    ref = length(rows[2][2])
    for (lbl, d) in rows
        mr, lr = homopolymers(d)
        @printf("  %-28s %9d %8.4f %+8.2f%% %8d %9d\n",
                lbl, length(d), length(d) / n, 100 * (length(d) / ref - 1), mr, lr)
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

# ---- oligos with a fixed payload of L nt -----------------------------------
function huffman_oligos(code, msg, L)
    oligos = Vector{Vector{Int}}()
    i = 1
    while i <= length(msg)
        o = Int[]
        while i <= length(msg) && length(o) + length(code[msg[i]]) <= L
            append!(o, code[msg[i]]); i += 1
        end
        i <= length(msg) && append!(o, code[msg[i]][1:L-length(o)])  # incomplete codeword
        push!(oligos, o)
    end
    return oligos
end

# Longest prefix of msg[i:end] whose code fits in L digits; enc(k) codes msg[i:i+k-1].
function longest_fit(enc, L, nleft, guess)
    k = clamp(guess, 0, nleft)
    while k > 0 && length(enc(k)) > L; k -= 1; end
    while k < nleft && length(enc(k + 1)) <= L; k += 1; end
    return k
end

function eom_oligos(m, msg, L, avg)
    oligos, counts, i = Vector{Vector{Int}}(), Int[], 1
    while i <= length(msg)
        enc(k) = padic_encode(m, msg[i:i+k-1])
        k = longest_fit(enc, L, length(msg) - i + 1, floor(Int, L / avg) - 6)
        k == 0 && error("L too small")
        push!(oligos, enc(k)); push!(counts, k); i += k
    end
    return oligos, counts
end

# Header: h trits with k - kbase.  h and kbase are global parameters.
function fill_oligos(m, msg, L, h)
    oligos, counts, i = Vector{Vector{Int}}(), Int[], 1
    while i <= length(msg)
        k, d = padic_fill(m, msg, L - h; start = i)
        k == 0 && error("L too small")
        push!(oligos, d); push!(counts, k); i += k
    end
    return oligos, counts
end

function oligo_experiment(name, data::Vector{UInt8}, Ls)
    msg = Int.(data) .+ 1
    n = length(msg)
    freqs = ones(Int, 256)
    for s in msg; freqs[s] += 1; end
    p = freqs ./ sum(freqs)
    m = PadicModel(freqs; P = 3, N = 26)
    code = huffman_code(p, 3)
    avg = -sum(log(3, p[s]) for s in msg) / n          # trits per byte

    @printf("\n%s: oligos with a fixed payload of L nt (whole file %d bytes)\n", name, n)
    @printf("  %5s  %-34s %8s %10s %9s\n", "L", "method", "oligos", "bytes/oligo", "vs Huffman")
    for L in Ls
        rows = []

        ho = huffman_oligos(code, msg, L)
        @assert reduce(vcat, [huffman_decode(code, o) for o in ho]) == msg
        push!(rows, ("ternary Huffman, padded", length(ho)))

        eo, ec = eom_oligos(m, msg, L, avg)
        @assert reduce(vcat, [padic_decode(m, o) for o in eo]) == msg
        push!(rows, ("P=3 with EOM", length(eo)))

        h = 0                                          # smallest header that fits
        fo, fc = fill_oligos(m, msg, L, h)
        while 3^h < maximum(fc) - minimum(fc) + 1
            h += 1
            fo, fc = fill_oligos(m, msg, L, h)
        end
        kbase = minimum(fc)
        trits = [vcat(reverse(digits(k - kbase; base = 3, pad = h)),
                      d, zeros(Int, L - h - length(d))) for (k, d) in zip(fc, fo)]
        dec = Int[]
        for t in trits
            k = foldl((a, b) -> 3a + b, t[1:h]; init = 0) + kbase
            append!(dec, padic_decode_fixed(m, t[h+1:end], k))
            @assert dna_to_trits(trits_to_dna(t)) == t
        end
        @assert dec == msg
        push!(rows, (@sprintf("P=3 fixed length + %d-trit header", h), length(fo)))

        for (lbl, cnt) in rows
            @printf("  %5d  %-34s %8d %10.2f %+8.2f%%\n", L, lbl, cnt, n / cnt,
                    100 * (cnt / rows[1][2] - 1))
        end
    end
end

sources = isempty(ARGS) ? corpus() : [(f, read(f)) for f in ARGS]
for (name, data) in sources
    experiment(name, data)
end
for (name, data) in sources
    oligo_experiment(name, data, (100, 150, 200, 300))
end
