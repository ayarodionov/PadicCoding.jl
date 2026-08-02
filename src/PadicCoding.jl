#
# P-adic arithmetic coding
#
# Julia implementation of the algorithms from:
#   Anatoly Rodionov, Sergey Volkov, "P-adic arithmetic coding", 2007.
#
# A message is mapped to a semi interval [l, r) of [0, 1) whose edges lie on an
# equidistant grid G(P^N); edge points are kept as integer indexes in the ring
# mod P^N.  The IP ("hat", ^) transformation maps a grid index to the p-adic
# number of the path leading to it on the coding P-tree (it just reverses the
# N base-P digits).  Whenever the paths to the two edges of the current message
# interval get a common prefix, that prefix is pushed to the output and the
# interval is rescaled (PR rescaling).  The revised algorithm adds AR rescaling
# (based on the usual Archimedean distance) for intervals straddling a level-1
# grid point, which bounds the interval width from below and makes the coder
# work with fixed precision integers for arbitrary messages.
#
# Two encoder/decoder pairs are provided:
#   * padic_encode_pr / padic_decode_pr  -- basic algorithm, PR rescaling only
#     (pseudo code on pp. 17-19 of the paper);
#   * padic_encode / padic_decode        -- revised algorithm, PR + AR
#     (pseudo code on pp. 25-27 of the paper).
#
module PadicCoding

export PadicModel, eom,
    padic_encode, padic_decode,
    padic_encode_pr, padic_decode_pr,
    PadicModel2, padic_encode2, padic_decode2

const I128 = Int128

# ---------------------------------------------------------------------------
# p-adic helpers
# ---------------------------------------------------------------------------

"""
    hat(x, P, N)

The IP (Index-Path) transformation, operator `^` of the paper: reverse the
N base-P digits of `x`.  Maps a grid index to its path presented as a p-adic
integer number and back (it is idempotent: `hat(hat(x)) == x`).
"""
function hat(x::I128, P::Int, N::Int)
    y = zero(I128)
    for _ in 1:N
        y = y * P + x % P
        x ÷= P
    end
    return y
end

"p-adic valuation `ord_P(x)`: max k such that P^k divides x (x ≠ 0)."
function ordp(x::I128, P::Int)
    x = abs(x)
    k = 0
    while x % P == 0
        x ÷= P
        k += 1
    end
    return k
end

"`com_{P,N}(x, y)`: length of the common path of two level-N paths."
compath(x::I128, y::I128, P::Int, N::Int) = x == y ? N : ordp(x - y, P)

"i-th base-P digit of x (digit 0 is the least significant): operator `x[i]`."
pdigit(x::I128, i::Int, P::Int) = Int((x ÷ I128(P)^i) % P)

"`cut(x, 1, 1)`: remove the path digit at position 1 (used by AR rescaling)."
cut1(x::I128, P::Int) = x % P + P * (x ÷ (I128(P) * P))

"""
    rescale(x, n, P, N)

`lift(res(x^, n), n)^`: drop the first n path digits of index `x` and return
to index representation at level N (PR rescaling of one interval edge).
"""
rescale(x::I128, n::Int, P::Int, N::Int) =
    hat(hat(x, P, N) ÷ I128(P)^n, P, N)

# ---------------------------------------------------------------------------
# Model: static frequencies + EOM, subdivides [l, r) on the ring mod P^N
# ---------------------------------------------------------------------------

"""
    PadicModel(freqs; P = 2, N)

Static model over the alphabet `1 .. length(freqs)` with the given positive
symbol frequencies; a special EOM symbol (index `length(freqs) + 1`) with
frequency 1 is appended.  `P` must be a prime, `N` is the grid level: all
interval edges are indexes of the grid G(P^N).

The revised (PR + AR) algorithm guarantees the current interval is never
shorter than 2/P^2 - 1/P^N, so the total frequency must not exceed
2*P^(N-2) - 1.
"""
struct PadicModel
    P::Int
    N::Int
    R::I128           # P^N, number of grid intervals
    cum::Vector{Int}  # cumulative frequencies, cum[a] .. cum[a+1] for symbol a
    total::Int
end

function PadicModel(freqs::AbstractVector{<:Integer}; P::Int = 2, N::Int)
    isempty(freqs) && throw(ArgumentError("alphabet must not be empty"))
    all(>(0), freqs) || throw(ArgumentError("frequencies must be positive"))
    cum = cumsum(vcat(0, Vector{Int}(freqs), 1))  # EOM gets weight 1
    total = cum[end]
    R = I128(P)^N
    total <= 2 * (R ÷ (I128(P) * P)) - 1 ||
        throw(ArgumentError("grid level N too small for total frequency $total"))
    return PadicModel(P, N, R, cum, total)
end

"Index of the End-Of-Message symbol of the model."
eom(m::PadicModel) = length(m.cum) - 1

"`M.code(a, l, r)`: new subinterval for symbol `a` (r == 0 stands for P^N)."
function mcode(m::PadicModel, a::Int, l::I128, r::I128)
    R = r == 0 ? m.R : r
    w = R - l
    w >= m.total ||
        error("message interval too narrow (width $w < total $(m.total)); " *
              "increase N or use the PR+AR coder")
    lnew = l + (w * m.cum[a]) ÷ m.total
    rnew = l + (w * m.cum[a + 1]) ÷ m.total
    return lnew, (rnew == m.R ? zero(I128) : rnew)
end

"`M.decode(g, l, r)`: symbol whose subinterval contains grid point `g`."
function mdecode(m::PadicModel, g::I128, l::I128, r::I128)
    R = r == 0 ? m.R : r
    w = R - l
    off = g - l
    0 <= off < w || error("decoder point outside of message interval")
    lo, hi = 1, length(m.cum) - 1
    while lo < hi                       # binary search on cumulative bounds
        mid = (lo + hi + 1) >> 1
        if (w * m.cum[mid]) ÷ m.total <= off
            lo = mid
        else
            hi = mid - 1
        end
    end
    lnew, rnew = mcode(m, lo, l, r)
    return lnew, rnew, lo
end

# ---------------------------------------------------------------------------
# Input / output of P-bits (symbols 0 .. P-1)
# ---------------------------------------------------------------------------

mutable struct BitReader
    bits::Vector{Int}
    pos::Int
end
BitReader(bits) = BitReader(bits, 1)

"Next P-bit from the stream; an exhausted stream yields the dropped trailing 0s."
getb!(io::BitReader) =
    io.pos <= length(io.bits) ? (b = io.bits[io.pos]; io.pos += 1; b) : 0

"`(I.getB(n) • P^T_n)^`: read n P-bits and build the index contribution."
function getpath!(io::BitReader, n::Int, P::Int)
    x = zero(I128)
    for _ in 1:n
        x = x * P + getb!(io)
    end
    return x
end

"`O.pushB(ext(x, n))`: push the first n path digits of x (m0 first)."
function pushdigits!(out::Vector{Int}, x::I128, n::Int, P::Int)
    for i in 0:n-1
        push!(out, pdigit(x, i, P))
    end
end

"`O.pushB(p, n)`: push P-bit p n times."
function pushrep!(out::Vector{Int}, p::Int, n::Int)
    for _ in 1:n
        push!(out, p)
    end
end

# ---------------------------------------------------------------------------
# selectPoint: the point of [l, r) with the minimal path
# ---------------------------------------------------------------------------

"""
    select_path(l, r, P, N, R)

`selectPoint(l, r) = min(x^ : l <= x < r)` returned in path representation:
the p-adic number of the shortest path into the final interval.  Minimizing
the path as an integer means minimizing the index digits lexicographically
from the least significant one up, which is done greedily digit by digit.
"""
function select_path(l::I128, r::I128, P::Int, N::Int, R::I128)
    lo = l
    hi = r == 0 ? R : r
    u = zeros(Int, N)                # index digits, least significant first
    for k in 1:N
        for d in 0:P-1
            a0 = lo + mod(d - lo, P) # first a >= lo with a % P == d
            if a0 < hi
                u[k] = d
                lo = a0 ÷ P
                hi = cld(hi - d, P)
                break
            end
        end
    end
    q = zero(I128)                   # path number: q = sum u[k] * P^(N-k)
    for k in 1:N
        q = q * P + u[k]
    end
    return q
end

# ---------------------------------------------------------------------------
# Basic algorithm: PR rescaling only (paper pp. 17-19)
# ---------------------------------------------------------------------------

"""
    padic_encode_pr(m, msg) -> Vector{Int}

Basic p-adic encoder (PR rescaling only).  `msg` contains symbols
`1 .. length(freqs)`.  Returns the stream of P-bits.  May fail on messages
whose interval straddles a level-1 grid point for too long -- use
[`padic_encode`](@ref) for the robust variant.
"""
function padic_encode_pr(m::PadicModel, msg::AbstractVector{<:Integer})
    P, N = m.P, m.N
    out = Int[]
    l = r = zero(I128)
    for a in msg
        l, r = mcode(m, Int(a), l, r)
        n = compath(hat(l, P, N), hat(mod(r - 1, m.R), P, N), P, N)
        if n > 0
            pushdigits!(out, hat(l, P, N), n, P)
            l = rescale(l, n, P, N)
            r = rescale(r, n, P, N)
        end
    end
    l, r = mcode(m, eom(m), l, r)
    q = select_path(l, r, P, N, m.R)
    pushdigits!(out, q, ndigits(q, base = P), P)  # lnz: keep no trailing zeros
    return out
end

"""
    padic_decode_pr(m, bits) -> Vector{Int}

Decoder matching [`padic_encode_pr`](@ref).
"""
function padic_decode_pr(m::PadicModel, bits::Vector{Int})
    P, N = m.P, m.N
    io = BitReader(bits)
    out = Int[]
    l = r = zero(I128)
    g = getpath!(io, N, P)
    while true
        l, r, a = mdecode(m, g, l, r)
        a == eom(m) && break
        push!(out, a)
        n = compath(hat(l, P, N), hat(mod(r - 1, m.R), P, N), P, N)
        if n > 0
            l = rescale(l, n, P, N)
            r = rescale(r, n, P, N)
            g = rescale(g, n, P, N)
            g += getpath!(io, n, P)
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Revised algorithm: PR + AR rescaling (paper pp. 25-27)
# ---------------------------------------------------------------------------

# AR? predicate: the interval sits in the smallest neighbourhood of a level-1
# grid point: r^[0] - l^[0] == 1, l^[1] == P-1, r^[1] == 0.
function ar_predicate(lh::I128, rh::I128, P::Int)
    return pdigit(rh, 0, P) - pdigit(lh, 0, P) == 1 &&
           pdigit(lh, 1, P) == P - 1 &&
           pdigit(rh, 1, P) == 0
end

# toLeft?(n, l) = l[0] >= n : the stable point lies to the left of l.
toleft(sp::Int, lh::I128, P::Int) = pdigit(lh, 0, P) >= sp
# toRight?(n, r) = r[0] < n OR r == n : the stable point lies to the right of r
# (the path {n, 0, ..., 0} is the p-adic number n).
toright(sp::Int, rh::I128, P::Int) = pdigit(rh, 0, P) < sp || rh == sp

"""
    padic_encode(m, msg) -> Vector{Int}

Revised p-adic encoder with PR and AR rescaling.  `msg` contains symbols
`1 .. length(freqs)`; the EOM symbol is appended automatically.  Returns the
stream of P-bits (0 .. P-1); for P = 2 these are plain bits.
"""
function padic_encode(m::PadicModel, msg::AbstractVector{<:Integer})
    P, N = m.P, m.N
    out = Int[]
    l = r = zero(I128)
    sp = 0    # stable point of AR (level-1 grid point index, 0 < sp < P)
    spn = 0   # number of times AR rescaling was applied

    # Flush a pending AR transformation once the interval leaves the stable
    # point: push the common path it determines and rescale by one digit.
    function flush_sp!()
        spn == 0 && return
        lh = hat(l, P, N)
        rh = hat(r, P, N)
        if toleft(sp, lh, P)
            push!(out, sp)
            pushrep!(out, 0, spn)
        elseif toright(sp, rh, P)
            push!(out, sp - 1)
            pushrep!(out, P - 1, spn)
        else
            return
        end
        l = rescale(l, 1, P, N)
        r = rescale(r, 1, P, N)
        sp = 0
        spn = 0
    end

    for a in msg
        l, r = mcode(m, Int(a), l, r)
        flush_sp!()
        if spn == 0                                    # PR rescaling
            lh = hat(l, P, N)
            n = compath(lh, hat(mod(r - 1, m.R), P, N), P, N)
            if n > 0
                pushdigits!(out, lh, n, P)
                l = rescale(l, n, P, N)
                r = rescale(r, n, P, N)
            end
        end
        while true                                     # AR rescaling
            lh = hat(l, P, N)
            rh = hat(r, P, N)
            ar_predicate(lh, rh, P) || break
            sp == 0 && (sp = pdigit(rh, 0, P))
            spn += 1
            l = hat(cut1(lh, P), P, N)                 # lift(cut(l^,1,1),1)^
            r = hat(cut1(rh, P), P, N)
        end
    end

    l, r = mcode(m, eom(m), l, r)
    flush_sp!()
    if spn == 0
        q = select_path(l, r, P, N, m.R)
        pushdigits!(out, q, ndigits(q, base = P), P)
    else
        push!(out, sp)  # we already have a point of level 1
    end
    return out
end

"""
    padic_decode(m, bits) -> Vector{Int}

Decoder matching [`padic_encode`](@ref).
"""
function padic_decode(m::PadicModel, bits::Vector{Int})
    P, N = m.P, m.N
    io = BitReader(bits)
    out = Int[]
    l = r = zero(I128)
    sp = 0
    spn = 0
    g = getpath!(io, N, P)
    while true
        l, r, a = mdecode(m, g, l, r)
        a == eom(m) && break
        push!(out, a)

        if spn != 0
            lh = hat(l, P, N)
            rh = hat(r, P, N)
            if toleft(sp, lh, P) || toright(sp, rh, P)
                l = rescale(l, 1, P, N)
                r = rescale(r, 1, P, N)
                g = rescale(g, 1, P, N)
                g += getpath!(io, 1, P)
                sp = 0
                spn = 0
            end
        end

        if spn == 0                                    # PR rescaling
            n = compath(hat(l, P, N), hat(mod(r - 1, m.R), P, N), P, N)
            if n > 0
                l = rescale(l, n, P, N)
                r = rescale(r, n, P, N)
                g = rescale(g, n, P, N)
                g += getpath!(io, n, P)
            end
        end

        while true                                     # AR rescaling
            lh = hat(l, P, N)
            rh = hat(r, P, N)
            ar_predicate(lh, rh, P) || break
            sp == 0 && (sp = pdigit(rh, 0, P))
            spn += 1
            l = hat(cut1(lh, P), P, N)
            r = hat(cut1(rh, P), P, N)
            g = hat(cut1(hat(g, P, N), P), P, N)
            g += getpath!(io, 1, P)
        end
    end
    return out
end

include("PadicCoding2.jl")
using .PadicCoding2: PadicModel2, padic_encode2, padic_decode2

end # module
