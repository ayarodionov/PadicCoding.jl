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
#
# The last part is energy shaping.  A transition p -> q costs
# E = ALPHA (q - p)^2 + BETA (3 - q): switching energy grows with the square
# of the voltage step, and the optional static term models a termination that
# draws more current at lower levels (BETA = 0 by default).  This is a simple
# model, not measured hardware.  Weighting the transitions by 2^(-mu E) gives,
# for each mu, the Markov chain with the most bits per symbol for its mean
# energy; the matcher uses it as its model.  It is compared with block codes
# that send k bits in 4 symbols using, for each previous level, the 2^k
# lowest-energy codewords.

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
# symbols so far (indexes into allowed(previous level)).  The coder calls it
# once per symbol in order, so the level is updated one step at a time.
function model_fn(p0, models = MODELS)
    t0, p = 0, p0                    # p = level after the first t0 symbols
    return (t, prefix) -> begin
        if t - 1 == t0 + 1
            p = allowed(p)[prefix[end]]; t0 += 1
        elseif t - 1 != t0
            t0, p = 0, p0
            for a in prefix
                p = allowed(p)[a]; t0 += 1
            end
        end
        models[p+1]
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
function roundtrip(u, n, p0, models = MODELS)
    idx = padic_decode_fixed(model_fn(p0, models), u, n)
    back = trim(padic_encode_fixed(model_fn(p0, models), idx))
    return length(back) <= length(u) && back == trim(u)
end

function inputs(k, rng)
    k <= 14 && return ([digits(x; base = 2, pad = k) for x in 0:2^k-1], "all")
    return (vcat([rand(rng, 0:1, k) for _ in 1:500],
                 [zeros(Int, k), ones(Int, k), [isodd(i) ? 1 : 0 for i in 1:k]]), "sampled")
end

function matcher_bits(n, rng; models = MODELS, kstart = block_bits_state(n))
    k = min(kstart, block_bits_state(n))                  # nothing can do better
    while k > 0
        us, how = inputs(k, rng)
        all(roundtrip(u, n, p0, models) for u in us for p0 in LEVELS) && return k, how
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

# ---- energy shaping -------------------------------------------------------------
const ALPHA = 1.0               # switching energy per (level step)^2
const BETA = 0.0                # static energy per level below the top one
energy(p, q) = ALPHA * (q - p)^2 + BETA * (3 - q)

# The maxentropic chain for weights B[p, q] = 2^(-mu E): P(q | p) = B v_q / (lambda v_p).
function shaped_chain(mu)
    B = [abs(p - q) != 3 ? 2.0^(-mu * energy(p, q)) : 0.0 for p in LEVELS, q in LEVELS]
    F = eigen(B)
    i = argmax(real.(F.values))
    lam, v = real(F.values[i]), abs.(real.(F.vectors[:, i]))
    T = [B[p+1, q+1] * v[q+1] / (lam * v[p+1]) for p in LEVELS, q in LEVELS]
    G = eigen(Matrix(T'))                                 # stationary distribution
    pi = abs.(real.(G.vectors[:, argmax(real.(G.values))])); pi ./= sum(pi)
    H = -sum(pi[p+1] * T[p+1, q+1] * log2(T[p+1, q+1]) for p in LEVELS, q in LEVELS
             if T[p+1, q+1] > 0)
    E = sum(pi[p+1] * T[p+1, q+1] * energy(p, q) for p in LEVELS, q in LEVELS)
    return T, H, E
end

shaped_models(T) = [PadicModel(max.(1, round.(Int, SCALE .* [T[p+1, q+1] for q in allowed(p)]));
                               P = 2, N = 40, eom = false) for p in LEVELS]

seq_energy(p0, lv) = sum(energy(i == 1 ? p0 : lv[i-1], lv[i]) for i in eachindex(lv))

# Largest k (bits per frame) that every tested input round-trips with; binary
# search, assuming that a k that works also works for smaller k.
function fixed_bits(n, models, kmax, rng)
    lo, hi = 0, kmax                       # lo works, hi + 1 is not tested
    while lo < hi
        mid = (lo + hi + 1) ÷ 2
        us, _ = inputs(mid, rng)
        if all(roundtrip(u, n, p0, models) for u in us for p0 in LEVELS)
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

# Fixed rate: frames of n symbols, k bits each, state carried across frames.
function matcher_stream(models, n, k, frames, rng)
    E, p0 = 0.0, 1
    for _ in 1:frames
        u = rand(rng, 0:1, k)
        idx = padic_decode_fixed(model_fn(p0, models), u, n)
        @assert trim(padic_encode_fixed(model_fn(p0, models), idx)) == trim(u)
        lv = to_levels(p0, idx)
        E += seq_energy(p0, lv)
        p0 = lv[end]
    end
    return k / n, E / (frames * n), E / (frames * k)
end

# Digits the encoder has committed after coding idx (no termination): every
# point of the final interval starts with them.
function committed(model_at, idx)
    e = PadicCoding.Encoder()
    for t in eachindex(idx)
        PadicCoding.encode_step!(e, model_at(t, view(idx, 1:t-1)), idx[t])
    end
    return e.out
end

# Variable rate: one continuous stream.  The decoder turns all input bits into
# levels; the receiver's encoder commits digits as levels arrive, and the
# stream needs the first S levels after which it has committed every input bit.
function variable_stream(models, H, nbits, rng)    # H: bits per level of the chain
    bits = rand(rng, 0:1, nbits + 256)   # more data follows the measured part
    S = ceil(Int, 1.1 * nbits / H) + 200
    idx = padic_decode_fixed(model_fn(1, models), bits, S)
    at = model_fn(1, models)
    e = PadicCoding.Encoder()
    used = 0
    for t in 1:S
        PadicCoding.encode_step!(e, at(t, view(idx, 1:t-1)), idx[t])
        length(e.out) >= nbits && (used = t; break)
    end
    used == 0 && error("not enough levels")
    @assert e.out[1:nbits] == bits[1:nbits]
    lv = to_levels(1, idx[1:used])
    E = seq_energy(1, lv)
    return nbits / used, E / used, E / nbits
end

# Block code: k bits -> 4 symbols, the 2^k lowest-energy valid codewords per
# previous level.
function block_codebooks(k)
    books = Vector{Vector{Vector{Int}}}()
    for p in LEVELS
        words = Vector{Vector{Int}}()
        for x in 0:255
            w = digits(x; base = 4, pad = 4)
            all(abs(w[i] - (i == 1 ? p : w[i-1])) != 3 for i in 1:4) && push!(words, w)
        end
        sort!(words; by = w -> (seq_energy(p, w), w))
        length(words) >= 2^k || error("$k bits do not fit")
        push!(books, words[1:2^k])
    end
    return books
end

function block_stream(k, frames, rng)
    books = block_codebooks(k)
    E, p0 = 0.0, 1
    for _ in 1:frames
        w = books[p0+1][rand(rng, 1:2^k)]
        E += seq_energy(p0, w)
        p0 = w[end]
    end
    return k / 4, E / (4frames), E / (k * frames)
end

function energy_report(rng)
    println("\n==== Energy shaping: E(p -> q) = $(ALPHA) (q - p)^2 + $(BETA) (3 - q) ====")
    gray = rand(rng, 0:3, 100_000)
    Eg = seq_energy(1, gray)
    @printf("  plain PAM-4 (full swings allowed): 2.0000 bits/symbol, E/symbol %.3f, E/bit %.3f\n",
            Eg / length(gray), Eg / (2length(gray)))
    n = 128
    @printf("\n  %-34s %9s %9s %8s   %s\n", "code", "bits/sym", "E/symbol", "E/bit",
            "theory (n -> inf): bits/sym, E/bit")
    for k in (7, 6, 5)
        @printf("  %-34s %9.4f %9.3f %8.3f\n", "block $(k)b4s, low-energy codebooks",
                block_stream(k, 50_000, rng)...)
    end
    for mu in (0.0, 0.25, 0.5, 1.0, 2.0)
        T, H, E = shaped_chain(mu)
        models = shaped_models(T)
        k = fixed_bits(n, models, floor(Int, n * H) + 2, rng)
        @printf("  %-34s %9.4f %9.3f %8.3f   %.4f, %.3f\n",
                "p-adic fixed, mu = $(mu), n = $n", matcher_stream(models, n, k, 500, rng)..., H, E / H)
        @printf("  %-34s %9.4f %9.3f %8.3f\n",
                "p-adic variable (one stream)", variable_stream(models, H, 60_000, rng)...)
    end
    println("\nfixed: the same k bits in every frame (what a link needs; k from sampled",
            "\ninputs above 14 bits).  variable: one continuous stream, as many bits per",
            "\nlevel as the chain carries on average, so the input must be buffered.")
end

energy_report(rng)
