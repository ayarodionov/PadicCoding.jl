# Benchmark of the p-adic arithmetic coder.
#
#   julia --project benchmark.jl [file]
#
# Measures encoding and decoding throughput of the revised (PR + AR) coder for
# P = 2, 3, 5 and of the basic PR-only coder for P = 2, over messages drawn
# from a skewed 256-symbol distribution (or the bytes of the given file).
# Self-contained: plain warmup + best-of-k timing, no external packages.

using PadicCoding
using Random
using Printf

"Best-of-`k` wall time of `f()` in seconds (after one warmup call)."
function bench(f; k::Int = 5)
    f()                                   # warmup / compile
    best = Inf
    for _ in 1:k
        t = @elapsed f()
        best = min(best, t)
    end
    return best
end

function make_message(len::Int, rng)
    # geometric-ish skewed byte distribution, entropy well below 8 bits
    w = [1.0 / (1 + i ÷ 8) for i in 0:255]
    cw = cumsum(w) ./ sum(w)
    return [searchsortedfirst(cw, rand(rng)) for _ in 1:len]
end

function report(name, len, tenc, tdec, nbits, P)
    bytes = len / 1e6
    @printf("%-18s len %8d | encode %8.3f ms %6.2f MB/s | decode %8.3f ms %6.2f MB/s | %5.3f bits/sym\n",
            name, len, 1e3 * tenc, bytes / tenc, 1e3 * tdec, bytes / tdec,
            nbits * log2(P) / len)
end

rng = MersenneTwister(42)

msg, source = if isempty(ARGS)
    make_message(1_000_000, rng), "random skewed bytes"
else
    Int.(read(ARGS[1])) .+ 1, ARGS[1]
end

freqs = ones(Int, 256)
for a in msg
    freqs[a] += 1
end

println("source: $source, $(length(msg)) symbols")
println()

for P in (2, 3, 5)
    N = P == 2 ? 40 : P == 3 ? 26 : 18
    m = PadicModel(freqs; P = P, N = N)
    for len in (10_000, 100_000, length(msg))
        sub = view(msg, 1:len)
        bits = padic_encode(m, sub)
        @assert padic_decode(m, bits) == sub
        tenc = bench(() -> padic_encode(m, sub))
        tdec = bench(() -> padic_decode(m, bits))
        report("PR+AR  P=$P N=$N", len, tenc, tdec, length(bits), P)
    end
    println()
end

# Basic PR-only coder for comparison (P = 2).  A large N keeps the interval
# wide enough that integer truncation rarely hurts.
let m = PadicModel(freqs; P = 2, N = 40)
    for len in (10_000, 100_000, length(msg))
        sub = view(msg, 1:len)
        bits = padic_encode_pr(m, sub)
        @assert padic_decode_pr(m, bits) == sub
        tenc = bench(() -> padic_encode_pr(m, sub))
        tdec = bench(() -> padic_decode_pr(m, bits))
        report("PR-only P=2 N=40", len, tenc, tdec, length(bits), 2)
    end
end
println()

# Optimized word-level P = 2 coder (PadicCoding2.jl), bit-identical output.
let m = PadicModel2(freqs; N = 40)
    for len in (10_000, 100_000, length(msg))
        sub = view(msg, 1:len)
        bits = padic_encode2(m, sub)
        @assert padic_decode2(m, bits) == sub
        tenc = bench(() -> padic_encode2(m, sub))
        tdec = bench(() -> padic_decode2(m, bits))
        report("fast   P=2 N=40", len, tenc, tdec, length(bits), 2)
    end
end
