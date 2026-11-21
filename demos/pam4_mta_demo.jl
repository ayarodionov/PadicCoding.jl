# Experiment: PAM-4 without full swings (0 <-> 3), as targeted by the MTA
# ("maximum transition avoidance") coding of GDDR6X.
#
#   julia --project demos/pam4_mta_demo.jl
#
# A PAM-4 symbol carries 2 bits, but with the transitions 0 -> 3 and 3 -> 0
# forbidden, a level is followed by 3 choices (after 0 or 3) or 4 (after 1 or
# 2).  The number of valid sequences grows like lambda^n, where lambda is the
# largest eigenvalue of the transition matrix, so the capacity is
# log2(lambda) = log2((3 + sqrt 17) / 2) = 1.8325 bits per symbol.
#
# Fixed frames of n symbols are compared:
#   * block codes: the largest k such that the frame has 2^k valid sequences,
#     either valid after any previous level (first symbol 1 or 2) or with one
#     codebook per previous level.  The second is the best any code can do for
#     the frame, but for large n it needs a table of 2^k codewords;
#   * the p-adic coder as a distribution matcher: the *decoder* maps k input
#     bits to n levels with a Markov model that only allows legal transitions
#     (maxentropic probabilities, one model per previous level), and the
#     *encoder* recovers the bits.  This works when the final interval of the
#     n levels contains only one k-bit point, which is checked for every input
#     (all inputs for k <= 14, sampled above).

using PadicCoding, Random, Printf, LinearAlgebra

const LEVELS = 0:3
allowed(p) = [q for q in LEVELS if abs(q - p) != 3]       # no 0 <-> 3
const A = [abs(p - q) != 3 ? 1.0 : 0.0 for p in LEVELS, q in LEVELS]

# ---- capacity and the maxentropic Markov chain --------------------------------
function perron(A)
    F = eigen(Symmetric(A))
    i = argmax(F.values)
    return F.values[i], abs.(F.vectors[:, i])
end
const LAMBDA, V = perron(A)
const CAP = log2(LAMBDA)

# P(q | p) = A[p, q] v[q] / (lambda v[p])
maxent(p) = [V[q+1] / (LAMBDA * V[p+1]) for q in allowed(p)]

# ---- block-code limits ---------------------------------------------------------
# count[p+1] = number of valid n-sequences after previous level p.
function counts_after(n)
    c = ones(BigInt, 4)                                   # length 0
    for _ in 1:n
        c = [sum(c[q+1] for q in allowed(p)) for p in LEVELS]   # c_n(p) = sum c_{n-1}(q)
    end
    return c
end
# Sequences of length n with first symbol 1 or 2 (valid after any level).
free_count(n) = n == 0 ? big(1) : sum(counts_after(n - 1)[q+1] for q in (1, 2))
block_bits_free(n) = ndigits(free_count(n); base = 2) - 1             # floor(log2)
block_bits_state(n) = ndigits(minimum(counts_after(n)); base = 2) - 1

# ---- the matcher ---------------------------------------------------------------
const SCALE = 2^20
const MODELS = [PadicModel(max.(1, round.(Int, SCALE .* maxent(p))); P = 2, N = 40,
                           eom = false) for p in LEVELS]

# Per-step model: the state is the previous level, found from p0 and the
# symbols so far (indexes into allowed(previous level)).
function model_fn(p0)
    return (t, prefix) -> begin
        p = p0
        for a in prefix
            p = allowed(p)[a]
        end
        MODELS[p+1]
    end
end

function to_levels(p0, idx)
    lv, p = Int[], p0
    for a in idx
        p = allowed(p)[a]; push!(lv, p)
    end
    return lv
end

trim(v) = v[1:something(findlast(!=(0), v), 0)]

# Send k bits u as n levels after previous level p0; true if they come back.
function roundtrip(u, n, p0)
    idx = padic_decode_fixed(model_fn(p0), u, n)
    back = trim(padic_encode_fixed(model_fn(p0), idx))
    return length(back) <= length(u) && back == trim(u)
end

function inputs(k, rng)
    k <= 14 && return ([digits(x; base = 2, pad = k) for x in 0:2^k-1], "all")
    return (vcat([rand(rng, 0:1, k) for _ in 1:500],
                 [zeros(Int, k), ones(Int, k), [isodd(i) ? 1 : 0 for i in 1:k]]), "sampled")
end

function matcher_bits(n, rng)
    k = block_bits_state(n)                               # nothing can do better
    while k > 0
        us, how = inputs(k, rng)
        all(roundtrip(u, n, p0) for u in us for p0 in LEVELS) && return k, how
        k -= 1
    end
    return 0, "-"
end

# ---- report -------------------------------------------------------------------
rng = Xoshiro(4)
@printf("capacity log2((3 + sqrt 17)/2) = %.4f bits/symbol (unconstrained PAM-4: 2)\n", CAP)
println("\nbits per frame of n symbols, no 0<->3 transitions anywhere (also across frames)")
@printf("  %4s %7s %10s %11s %8s %9s %9s %9s %9s\n", "n", "n·cap", "block,free",
        "block,state", "p-adic", "checked", "p-adic/n", "vs state", "vs 7b4s")
for n in (4, 8, 16, 32, 64, 128)
    k, how = matcher_bits(n, rng)
    kb = block_bits_state(n)
    @printf("  %4d %7.1f %10d %11d %8d %9s %9.4f %8.2f%% %8.2f%%\n", n, n * CAP,
            block_bits_free(n), kb, k, how, k / n, 100 * (k / kb - 1),
            100 * (k / (1.75n) - 1))
end
println("\nvs 7b4s: against 7 bits per 4 symbols (1.75 bits/symbol), which GDDR6X's MTA",
        "\nappears to use.  For n = 4 that equals block,state; a code valid after any",
        "\nlevel (block,free) fits only 6 bits, so 7b4s needs the previous level or",
        "\nextra signalling (the MTA patents use the DBI line).")

# A long stream, frames of n = 64 with the state carried over: no full swings.
n = 64
k, _ = matcher_bits(n, rng)
bits = rand(rng, 0:1, 2000k)
levels, p0 = Int[], 1
for f in 1:2000
    u = bits[(f-1)*k+1:f*k]
    idx = padic_decode_fixed(model_fn(p0), u, n)
    back = trim(padic_encode_fixed(model_fn(p0), idx))
    @assert back == trim(u)
    lv = to_levels(p0, idx)
    append!(levels, lv)
    global p0 = lv[end]
end
swings(lv) = count(i -> abs(lv[i] - lv[i-1]) == 3, 2:length(lv))
gray = [(0, 1, 3, 2)[2bits[i] + bits[i+1] + 1] for i in 1:2:length(bits)]
@printf("\nstream of %d bits: p-adic matcher %d symbols, %d full swings;", length(bits),
        length(levels), swings(levels))
@printf("\n  plain PAM-4 (2 bits/symbol, Gray) %d symbols, %d full swings\n",
        length(gray), swings(gray))
