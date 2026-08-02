#
# Optimized p-adic arithmetic coder specialized to P = 2.
#
# Same revised (PR + AR) algorithm as `padic_encode` / `padic_decode` in
# PadicCoding.jl and bit-identical output, but with the p-adic machinery
# lowered to word operations on UInt64:
#
#   * the IP (`^`) transformation is a hardware bit reversal,
#   * ord_2 / com is `trailing_zeros` of an XOR,
#   * the ring mod 2^N is plain masking (R = 1 << N),
#   * the AR stable point is always 1, so only the counter `spn` is tracked,
#   * the output is a packed BitVector.
#
# Requires 2 <= N <= 62; the total model frequency must be < 2^(N-1).
#
module PadicCoding2

export PadicModel2, padic_encode2, padic_decode2

"""
    PadicModel2(freqs; N)

Static binary (P = 2) model over the alphabet `1 .. length(freqs)` with an
appended EOM symbol of frequency 1, on the grid G(2^N).
"""
struct PadicModel2
    N::Int
    R::UInt64            # 2^N
    mask::UInt64         # R - 1
    shift::Int           # 64 - N, for bit reversal
    cum::Vector{UInt64}
    total::UInt64
end

function PadicModel2(freqs::AbstractVector{<:Integer}; N::Int)
    isempty(freqs) && throw(ArgumentError("alphabet must not be empty"))
    all(>(0), freqs) || throw(ArgumentError("frequencies must be positive"))
    2 <= N <= 62 || throw(ArgumentError("N must be in 2 .. 62"))
    cum = cumsum(vcat(UInt64(0), UInt64.(freqs), UInt64(1)))
    total = cum[end]
    total <= (UInt64(1) << (N - 1)) - 1 ||
        throw(ArgumentError("grid level N too small for total frequency $total"))
    R = UInt64(1) << N
    return PadicModel2(N, R, R - 1, 64 - N, cum, total)
end

eom2(m::PadicModel2) = length(m.cum) - 1

# ---------------------------------------------------------------------------
# word-level p-adic primitives
# ---------------------------------------------------------------------------

"IP transformation: reverse the N low bits of x."
@inline hat(m::PadicModel2, x::UInt64) = bitreverse(x) >> m.shift

"com_{2,N}: length of the common path of two level-N paths."
@inline compath(x::UInt64, y::UInt64, N::Int) =
    x == y ? N : trailing_zeros(x ⊻ y)

"cut(x, 1, 1): remove the path bit at position 1."
@inline cut1(x::UInt64) = (x & 0x1) | ((x >> 2) << 1)

"lift(res(x^, n), n)^: PR rescaling of one interval edge."
@inline resc(m::PadicModel2, x::UInt64, n::Int) = hat(m, hat(m, x) >> n)

# AR?: l^ ends in bits (l1 l0) = 10, r^ ends in (r1 r0) = 01.
@inline ar2(lh::UInt64, rh::UInt64) = (lh & 0x3 == 0x2) && (rh & 0x3 == 0x1)

# ---------------------------------------------------------------------------
# model on the ring mod 2^N (r == 0 stands for 2^N)
# ---------------------------------------------------------------------------

@inline function mcode(m::PadicModel2, a::Int, l::UInt64, r::UInt64)
    R = r == 0 ? m.R : r
    w = R - l
    w >= m.total ||
        error("message interval too narrow (width $w < total $(m.total)); increase N")
    lnew = l + UInt64((UInt128(w) * m.cum[a]) ÷ m.total)
    rnew = l + UInt64((UInt128(w) * m.cum[a + 1]) ÷ m.total)
    return lnew, rnew & m.mask
end

function mdecode(m::PadicModel2, g::UInt64, l::UInt64, r::UInt64)
    R = r == 0 ? m.R : r
    w = R - l
    off = g - l
    off < w || error("decoder point outside of message interval")
    lo, hi = 1, length(m.cum) - 1
    while lo < hi
        mid = (lo + hi + 1) >> 1
        if UInt64((UInt128(w) * m.cum[mid]) ÷ m.total) <= off
            lo = mid
        else
            hi = mid - 1
        end
    end
    lnew, rnew = mcode(m, lo, l, r)
    return lnew, rnew, lo
end

# ---------------------------------------------------------------------------
# bit i/o
# ---------------------------------------------------------------------------

"Push the n low bits of x, least significant first (O.pushB(ext(x, n)))."
@inline function pushlow!(out::BitVector, x::UInt64, n::Int)
    for i in 0:n-1
        push!(out, (x >> i) & 0x1 == 0x1)
    end
end

mutable struct BitReader2{V<:AbstractVector}
    bits::V
    pos::Int
end

@inline getb!(io::BitReader2) =
    io.pos <= length(io.bits) ? (b = Bool(io.bits[io.pos]); io.pos += 1;
                                 UInt64(b)) : UInt64(0)

"Read n bits and build the index contribution ((I.getB(n) • P^T_n)^)."
@inline function getpath!(io::BitReader2, n::Int)
    x = UInt64(0)
    for _ in 1:n
        x = (x << 1) | getb!(io)
    end
    return x
end

# ---------------------------------------------------------------------------
# selectPoint
# ---------------------------------------------------------------------------

"Minimal path (as p-adic number) into [l, r), built greedily bit by bit."
function select_path2(m::PadicModel2, l::UInt64, r::UInt64)
    lo = l
    hi = r == 0 ? m.R : r
    q = UInt64(0)
    for _ in 1:m.N
        a0 = lo + (lo & 0x1)          # first even index >= lo
        d = UInt64(0)
        if a0 >= hi                   # no even index in range: take odd
            d = UInt64(1)
            a0 = lo | 0x1
        end
        q = (q << 1) | d
        lo = a0 >> 1
        hi = (hi - d + 1) >> 1        # cld(hi - d, 2)
    end
    return q
end

# ---------------------------------------------------------------------------
# encoder / decoder (revised algorithm, PR + AR)
# ---------------------------------------------------------------------------

# Flush a pending AR transformation; for P = 2 the stable point is always 1.
@inline function flush_sp!(out::BitVector, m::PadicModel2,
                           l::UInt64, r::UInt64, spn::Int)
    spn == 0 && return l, r, spn
    lh = hat(m, l)
    rh = hat(m, r)
    if lh & 0x1 == 0x1                        # toLeft?(1, l^)
        push!(out, true)
        for _ in 1:spn
            push!(out, false)
        end
    elseif rh & 0x1 == 0x0 || rh == 1         # toRight?(1, r^)
        push!(out, false)
        for _ in 1:spn
            push!(out, true)
        end
    else
        return l, r, spn
    end
    return resc(m, l, 1), resc(m, r, 1), 0
end

"""
    padic_encode2(m, msg) -> BitVector

Optimized P = 2 encoder; bit-identical to `PadicCoding.padic_encode` with the
same frequencies and N.
"""
function padic_encode2(m::PadicModel2, msg::AbstractVector{<:Integer})
    out = BitVector()
    l = r = UInt64(0)
    spn = 0
    for a in msg
        l, r = mcode(m, Int(a), l, r)
        l, r, spn = flush_sp!(out, m, l, r, spn)
        if spn == 0                                    # PR rescaling
            lh = hat(m, l)
            n = compath(lh, hat(m, (r - 1) & m.mask), m.N)
            if n > 0
                pushlow!(out, lh, n)
                l = resc(m, l, n)
                r = resc(m, r, n)
            end
        end
        while true                                     # AR rescaling
            lh = hat(m, l)
            rh = hat(m, r)
            ar2(lh, rh) || break
            spn += 1
            l = hat(m, cut1(lh))
            r = hat(m, cut1(rh))
        end
    end
    l, r = mcode(m, eom2(m), l, r)
    l, r, spn = flush_sp!(out, m, l, r, spn)
    if spn == 0
        q = select_path2(m, l, r)
        pushlow!(out, q, 64 - leading_zeros(q))        # drop trailing zeros
    else
        push!(out, true)                               # point of level 1
    end
    return out
end

"""
    padic_decode2(m, bits) -> Vector{Int}

Decoder matching [`padic_encode2`](@ref); accepts any 0/1 vector.
"""
function padic_decode2(m::PadicModel2, bits::AbstractVector)
    io = BitReader2(bits, 1)
    out = Int[]
    l = r = UInt64(0)
    spn = 0
    g = getpath!(io, m.N)
    while true
        l, r, a = mdecode(m, g, l, r)
        a == eom2(m) && break
        push!(out, a)

        if spn != 0
            lh = hat(m, l)
            rh = hat(m, r)
            if lh & 0x1 == 0x1 || rh & 0x1 == 0x0 || rh == 1
                l = resc(m, l, 1)
                r = resc(m, r, 1)
                g = resc(m, g, 1)
                g += getpath!(io, 1)
                spn = 0
            end
        end

        if spn == 0                                    # PR rescaling
            n = compath(hat(m, l), hat(m, (r - 1) & m.mask), m.N)
            if n > 0
                l = resc(m, l, n)
                r = resc(m, r, n)
                g = resc(m, g, n)
                g += getpath!(io, n)
            end
        end

        while true                                     # AR rescaling
            lh = hat(m, l)
            rh = hat(m, r)
            ar2(lh, rh) || break
            spn += 1
            l = hat(m, cut1(lh))
            r = hat(m, cut1(rh))
            g = hat(m, cut1(hat(m, g)))
            g += getb!(io)
        end
    end
    return out
end

end # module
