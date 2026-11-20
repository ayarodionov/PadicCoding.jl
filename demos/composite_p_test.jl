# Experiment: does the p-adic coder need P to be prime?
#
#   julia --project demos/composite_p_test.jl
#
# The paper states P prime (so that Q_p is a field), but the coder only uses
# base-P digits: hat() reverses them, and ordp(x - y, P) counts the trailing
# base-P digits shared by x and y, which is meaningful for any base P >= 2.
# This script stress-tests padic_encode / padic_decode for prime and composite
# P, with tight grid levels N that force frequent AR rescaling, and counts how
# often AR actually fires so the test is known to cover it.
#
# It also tests the width bound.  The paper (p. 25) says PR + AR keep the
# interval at least 2/P^2 - 1/P^N wide, so total frequency <= 2 P^(N-2) - 1.
# That is too optimistic: an interval that straddles a level-1 point but has
# l[1] < P-1 or r[1] > 0 is not AR-rescaled, and it can be just over 1/P^2
# wide (for P = 2: [0.24, 0.51)).  The bound that holds is
# width >= P^(N-2) + 1 grid units, so total <= P^(N-2) + 1, which is what
# the PadicModel constructor now checks.  Each trial is run at the smallest N
# allowed by each bound; models at the paper's limit are built with the plain
# struct constructor, which skips the check.

using PadicCoding, Random, Printf

# Count AR rescalings: wrap the AR predicate (same body as in src/).
const AR_HITS = Ref(0)
@eval PadicCoding function ar_predicate(lh::I128, rh::I128, P::Int)
    hit = pdigit(rh, 0, P) - pdigit(lh, 0, P) == 1 &&
          pdigit(lh, 1, P) == P - 1 &&
          pdigit(rh, 1, P) == 0
    hit && (Main.AR_HITS[] += 1)
    return hit
end

# Largest N with P^N <= 2^90, so (width * cumulative frequency) fits Int128.
nmax(P) = floor(Int, 90 / log2(P))
# Smallest N allowed for a total frequency by the paper's bound and by the
# corrected one.
paper_ok(P, N, total) = 2 * big(P)^(N - 2) - 1 >= total
safe_ok(P, N, total)  = big(P)^(N - 2) + 1 >= total
function nmin(P, total, ok)
    N = 2
    while !ok(P, N, total); N += 1; end
    return N
end

# PadicModel without the width check (to run at the paper's limit).
function unchecked_model(f, P, N)
    cum = cumsum(vcat(0, f, 1))
    return PadicModel(P, N, PadicCoding.I128(P)^N, cum, cum[end])
end

function random_freqs(rng)
    k = rand(rng, (1, 2, 3, 5, 17, 256))
    kind = rand(rng, 1:3)
    kind == 1 && return rand(rng, 1:10, k)                 # near uniform
    kind == 2 && return [2^rand(rng, 0:16) for _ in 1:k]   # strongly skewed
    f = ones(Int, k); f[rand(rng, 1:k)] = 10^6; return f   # one dominant symbol
end

function random_msg(rng, freqs)
    n = rand(rng, (0, 1, 2, 10, 300, 2000))
    kind = rand(rng, 1:4)
    kind == 1 && return fill(argmax(freqs), n)             # most likely symbol
    kind == 2 && return fill(argmin(freqs), n)             # least likely symbol
    kind == 3 && return rand(rng, 1:length(freqs), n)      # ignores the model
    c = cumsum(freqs)                                      # drawn from the model
    return [searchsortedfirst(c, rand(rng, 1:c[end])) for _ in 1:n]
end

# One trial: true if the message decodes back.
function roundtrip(m, msg)
    out = try padic_encode(m, msg) catch; return false, 0 end
    dec = try padic_decode(m, out) catch; return false, 0 end
    return dec == msg, length(out)
end

ideal_bits(m, msg) = -sum((log2((m.cum[a+1] - m.cum[a]) / m.total) for a in msg); init = 0.0) -
                     log2(1 / m.total)                    # + EOM

Ps = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 27, 30, 64, 100, 210, 256, 1000, 65536]
TRIALS = 600

@printf("%7s %4s %6s %12s %12s %12s %10s %11s\n", "P", "Nmax", "trials",
        "fail paperN", "fail safeN", "fail roomy", "AR fires", "overhead d")
for P in Ps
    rng = Xoshiro(P)
    fails = [0, 0, 0]
    ar = 0
    over = Float64[]
    for _ in 1:TRIALS
        f = random_freqs(rng)
        total = sum(f) + 1
        np, ns = nmin(P, total, paper_ok), nmin(P, total, safe_ok)
        ns > nmax(P) && continue                           # does not fit Int128
        msg = random_msg(rng, f)
        for (i, N) in enumerate((np, ns, rand(rng, ns:nmax(P))))
            m = i == 1 ? unchecked_model(f, P, N) : PadicModel(f; P = P, N = N)
            AR_HITS[] = 0
            ok, len = roundtrip(m, msg)
            ar += AR_HITS[]
            ok || (fails[i] += 1)
            ok && i == 3 && push!(over, len - ideal_bits(m, msg) / log2(P))
        end
    end
    @printf("%7d %4d %6d %12d %12d %12d %10d %11.2f\n",
            P, nmax(P), TRIALS, fails..., ar, sum(over) / length(over))
end
println("\nfail paperN / safeN: failures at the smallest N allowed by the paper's / the corrected bound;",
        "\nfail roomy: random N above that.  overhead d: mean output length minus the ideal",
        "\nlength (roomy N), in P-digits.")

# Is P = 4 (or 16) just P = 2 with its bits grouped?  Compare the streams.
println("\nSame message, same grid R = 2^40:  P = 2^k output unpacked to bits vs P = 2 output")
rng = Xoshiro(7)
for (P, N) in ((4, 20), (16, 10), (1024, 4))
    same = 0
    for _ in 1:200
        f = random_freqs(rng)
        sum(f) + 1 <= 2^(40 - 2 * round(Int, log2(P))) + 1 || continue
        msg = random_msg(rng, f)
        k = round(Int, log2(P))
        b2 = padic_encode(PadicModel(f; P = 2, N = 40), msg)
        bp = padic_encode(PadicModel(f; P = P, N = N), msg)
        unpacked = reduce(vcat, [reverse(digits(d; base = 2, pad = k)) for d in bp]; init = Int[])
        trim(v) = v[1:something(findlast(!=(0), v), 0)]
        same += trim(unpacked) == trim(b2)
    end
    @printf("  P = %4d: identical (up to trailing zeros) in %d of 200 messages\n", P, same)
end

# Throughput on text, prime vs composite P.
text = read(joinpath(@__DIR__, "..", "README.md"))
msg = Int.(text) .+ 1
f = ones(Int, 256); for s in msg; f[s] += 1; end
println("\nThroughput on README.md ($(length(msg)) bytes), N = nmax(P):")
for P in (2, 3, 4, 5, 16, 256, 65536)
    m = PadicModel(f; P = P, N = nmax(P))
    padic_encode(m, msg[1:10])
    te = @elapsed out = padic_encode(m, msg)
    td = @elapsed dec = padic_decode(m, out)
    @assert dec == msg
    @printf("  P = %6d  N = %3d  enc %6.2f MB/s  dec %6.2f MB/s  %.4f bits/byte\n",
            P, nmax(P), length(msg) / te / 1e6, length(msg) / td / 1e6,
            length(out) * log2(P) / length(msg))
end
