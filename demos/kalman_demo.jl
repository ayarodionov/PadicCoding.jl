# Demo: Kalman-filter predictive models + the adaptive p-adic coder.
#
#   julia --project demos/kalman_demo.jl
#
# A smooth numeric signal (a synthetic 1-D GPS-like position track, 24-bit
# integer samples) is coded sample by sample. Each sample becomes a residual
# symbol around a predicted centre c:
#
#   x - c in -K .. K   -> one of 2K+1 residual bins
#   otherwise          -> ESC, followed by the raw value as 3 uniform bytes
#
# and the per-step distribution over those symbols comes from one of
#
#   * static delta   -- two-pass histogram of x[t] - x[t-1] (table not counted)
#   * adaptive delta -- one-pass adaptive frequency counter of x[t] - x[t-1]
#   * fixed Kalman   -- constant-velocity Kalman filter; its predictive
#                       Gaussian N(mu, S) is discretised into the bins. The
#                       best (q, r) of the bank below, chosen in hindsight.
#   * Kalman mixture -- a bank of Kalman filters with different process /
#                       measurement noise, mixed by their recent predictive
#                       likelihood (Bayesian model averaging with forgetting),
#                       so the model adapts when the signal changes regime.
#
# The coder is `padic_encode2(model_at, msg)`: `model_at(t, prefix)` builds a
# fresh PadicModel2 for every step. The decoder runs an identical copy of the
# filter on the samples decoded so far, so no model parameters are sent.
# Encoder and decoder must compute bit-identical frequencies; that holds here
# (same code, same machine). Across platforms, use fixed-point filter
# arithmetic.

using PadicCoding
using Random
using Printf

const K     = 255             # residual bins -K .. K
const NBINS = 2K + 1
const ESC   = NBINS + 1       # escape symbol
const SCALE = 2^22            # probability -> integer frequency scale
const NGRID = 40
const XBITS = 24              # raw sample width
const XMAX  = 2^XBITS - 1
const UNIFORM_BYTE = PadicModel2(ones(Int, 256); N = NGRID)

# ---------------------------------------------------------------------------
# synthetic signal: regimes of motion and measurement quality, plus outliers
# ---------------------------------------------------------------------------

# (name, acceleration sigma, measurement sigma), units per sample
const REGIMES = [("parked",       0.0,  3.0),
                 ("cruise",       0.3,  3.0),
                 ("maneuver",     4.0,  3.0),
                 ("urban canyon", 0.3, 40.0)]

function make_track(n, rng)
    x = Vector{Int}(undef, n)
    regime = Vector{Int}(undef, n)
    p, v = 2.0^23, 20.0
    t = 1
    while t <= n
        k = rand(rng, 1:length(REGIMES))
        _, sa, sm = REGIMES[k]
        k == 1 && (v = 0.0)
        k != 1 && abs(v) < 1 && (v = 30 * randn(rng))
        for _ in 1:min(rand(rng, 500:2000), n - t + 1)
            v += sa * randn(rng)
            p += v
            z = p + sm * randn(rng)
            rand(rng) < 0.003 &&                         # multipath glitch
                (z += rand(rng, (-1, 1)) * rand(rng, 300:3000))
            x[t] = clamp(round(Int, z), 0, XMAX)
            regime[t] = k
            t += 1
        end
    end
    return x, regime
end

# ---------------------------------------------------------------------------
# discretising predictive Gaussians into residual-bin frequencies
# ---------------------------------------------------------------------------

"erfc with fractional error < 1.2e-7 (Numerical Recipes `erfcc`)."
function erfc_nr(x::Float64)
    z = abs(x)
    t = 1 / (1 + z / 2)
    a = t * exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 +
        t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 +
        t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
    return x >= 0 ? a : 2 - a
end
normcdf(z) = erfc_nr(-z / sqrt(2)) / 2

"""
Frequencies of the 2K+1 residual bins + ESC around centre `c` for the mixture
`sum(w[i] * N(mu[i], s2[i]))`. Every symbol keeps frequency >= 1.
"""
function mixture_freqs(c, w, mu, s2)
    prob = zeros(NBINS)
    for i in eachindex(w)
        s = sqrt(s2[i])
        lo = normcdf((c - K - 0.5 - mu[i]) / s)
        for j in 1:NBINS
            hi = normcdf((c - K + j - 0.5 - mu[i]) / s)
            prob[j] += w[i] * (hi - lo)
            lo = hi
        end
    end
    tail = max(0.0, 1 - sum(prob))
    f = Vector{Int}(undef, NBINS + 1)
    for j in 1:NBINS
        f[j] = 1 + floor(Int, prob[j] * SCALE)
    end
    f[ESC] = 1 + floor(Int, tail * SCALE)
    return f
end

# ---------------------------------------------------------------------------
# predictors: `predict(p) -> (centre, freqs)` and `observe!(p, x)`
# ---------------------------------------------------------------------------

"Adaptive frequency counter over x[t] - x[t-1] (centre = previous sample)."
mutable struct DeltaCounter
    last::Int
    counts::Vector{Int}
end
DeltaCounter() = DeltaCounter(2^23, ones(Int, NBINS + 1))
predict(p::DeltaCounter) = p.last, p.counts
function observe!(p::DeltaCounter, x)
    e = x - p.last
    p.counts[abs(e) <= K ? e + K + 1 : ESC] += 32
    if sum(p.counts) > 2^16                    # halve: forget old statistics
        p.counts .= (p.counts .+ 1) .>> 1
    end
    p.last = x
end

"Static two-pass delta histogram (centre = previous sample)."
mutable struct StaticDelta
    last::Int
    freqs::Vector{Int}
end
predict(p::StaticDelta) = p.last, p.freqs
observe!(p::StaticDelta, x) = (p.last = x)

function static_delta(x)
    f = ones(Int, NBINS + 1)
    last = 2^23
    for xi in x
        e = xi - last
        f[abs(e) <= K ? e + K + 1 : ESC] += 1
        last = xi
    end
    return () -> StaticDelta(2^23, copy(f))
end

"Constant-velocity Kalman filter, state (position, velocity), dt = 1."
mutable struct Kalman
    q::Float64                       # acceleration (process) noise variance
    r::Float64                       # measurement noise variance
    p::Float64; v::Float64           # predicted state
    P11::Float64; P12::Float64; P22::Float64   # predicted covariance
    skipped::Int
end
Kalman(q, r) = Kalman(q, r, 2.0^23, 0.0, 1e12, 0.0, 1e4, 0)

predictive(k::Kalman) = k.p, k.P11 + k.r       # mean, variance of next sample

function update!(k::Kalman, z)
    mu, S = predictive(k)
    nu = z - mu
    if nu^2 > 64S && k.skipped < 3             # gate isolated outliers
        k.skipped += 1
    else
        k.skipped = 0
        K1, K2 = k.P11 / S, k.P12 / S
        k.p += K1 * nu
        k.v += K2 * nu
        k.P22 -= K2 * k.P12
        k.P12 *= 1 - K1
        k.P11 *= 1 - K1
    end
    # time update to the next sample
    k.p += k.v
    k.P11 += 2k.P12 + k.P22 + k.q / 4
    k.P12 += k.P22 + k.q / 2
    k.P22 += k.q
end

"Bank of Kalman filters mixed by forgetting-weighted predictive likelihood."
mutable struct KalmanMix
    filters::Vector{Kalman}
    logw::Vector{Float64}
    forget::Float64
end
KalmanMix(params; forget = 0.9) =
    KalmanMix([Kalman(q, r) for (q, r) in params], zeros(length(params)), forget)

function weights(m::KalmanMix)
    w = exp.(m.logw .- maximum(m.logw))
    return w ./ sum(w)
end

function predict(m::KalmanMix)
    w = weights(m)
    mus = [predictive(k)[1] for k in m.filters]
    s2s = [predictive(k)[2] for k in m.filters]
    c = clamp(round(Int, sum(w .* mus)), 0, XMAX)
    return c, mixture_freqs(c, w, mus, s2s)
end

function observe!(m::KalmanMix, x)
    for (i, k) in enumerate(m.filters)
        mu, S = predictive(k)
        m.logw[i] = m.forget * m.logw[i] - (x - mu)^2 / 2S - log(S) / 2
        update!(k, x)
    end
end

# ---------------------------------------------------------------------------
# samples <-> symbols, and the per-step model function for the coder
# ---------------------------------------------------------------------------

"Map samples to symbols with predictor `pred`; also return ideal bit costs."
function tokenize(pred, x)
    syms = Int[]
    cost = zeros(length(x))
    for (t, xi) in enumerate(x)
        c, f = predict(pred)
        e = xi - c
        s = abs(e) <= K ? e + K + 1 : ESC
        push!(syms, s)
        cost[t] = -log2(f[s] / sum(f))
        if s == ESC
            for sh in (16, 8, 0)
                push!(syms, ((xi >> sh) & 0xff) + 1)
            end
            cost[t] += 24
        end
        observe!(pred, xi)
    end
    return syms, cost
end

"""
Stateful `model_at(t, prefix)` for `padic_encode2` / `padic_decode2`: replays
the previous symbol into the predictor, then returns the model for symbol t.
Decoded samples are collected in `xs`.
"""
function model_fn(pred)
    phase = 0                      # 0: residual symbol, 1..3: raw byte
    c, acc = 0, 0
    xs = Int[]
    f = function (t, prefix)
        if t > 1
            s = prefix[end]
            if phase == 0
                if s == ESC
                    phase, acc = 1, 0
                else
                    push!(xs, c + s - K - 1)
                    observe!(pred, xs[end])
                end
            else
                acc = (acc << 8) | (s - 1)
                if phase == 3
                    push!(xs, acc)
                    observe!(pred, acc)
                    phase = 0
                else
                    phase += 1
                end
            end
        end
        phase == 0 || return UNIFORM_BYTE
        c, freqs = predict(pred)
        return PadicModel2(freqs; N = NGRID)
    end
    return f, xs
end

"Encode `x` with fresh predictors from `mk()`, verify the roundtrip."
function run(mk, x)
    syms, cost = tokenize(mk(), x)
    enc, _ = model_fn(mk())
    bits = padic_encode2(enc, syms)
    dec, xs = model_fn(mk())
    @assert padic_decode2(dec, bits) == syms
    @assert xs == x
    return length(bits), cost
end

# ---------------------------------------------------------------------------

rng = MersenneTwister(2026)
n = 20_000
x, regime = make_track(n, rng)

bank = [(q, r) for q in (0.0, 0.1, 16.0) for r in (9.0, 1600.0)]

results = Tuple{String,Int,Vector{Float64}}[]
push!(results, ("static delta (2-pass)", run(static_delta(x), x)...))
push!(results, ("adaptive delta counter", run(DeltaCounter, x)...))

best = nothing
for (q, r) in bank
    b, cost = run(() -> KalmanMix([(q, r)]), x)
    if best === nothing || b < best[2]
        global best = ("fixed Kalman q=$q r=$r", b, cost)
    end
end
push!(results, best)
push!(results, ("Kalman mixture (6 filters)", run(() -> KalmanMix(bank), x)...))

raw = XBITS * n
println("samples: $n  (raw: $raw bits = $(XBITS) bits/sample)")
println("regimes: ", join(["$(REGIMES[k][1]) $(count(==(k), regime))"
                           for k in eachindex(REGIMES)], ", "))
println()
@printf("%-28s %10s %12s %10s\n", "model", "bits", "bits/sample", "vs raw")
for (name, b, _) in results
    @printf("%-28s %10d %12.3f %9.2fx\n", name, b, b / n, raw / b)
end
println()
println("ideal cost per regime (bits/sample):")
@printf("%-28s", "")
for (name, _, _) in REGIMES
    @printf(" %13s", name)
end
println()
for (name, _, cost) in results
    @printf("%-28s", name)
    for k in eachindex(REGIMES)
        @printf(" %13.3f", sum(cost[regime .== k]) / count(==(k), regime))
    end
    println()
end
