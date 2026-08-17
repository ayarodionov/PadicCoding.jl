// Optimized p-adic arithmetic coder specialized to P = 2.
//
// Direct C++ port of src/PadicCoding2.jl (see that file for the derivation
// and the paper "P-adic arithmetic coding" by Rodionov & Volkov). Same
// revised (PR + AR) algorithm, bit-identical output to the Julia coder for
// the same frequencies and N. All p-adic machinery is lowered to word
// operations on uint64_t:
//
//   * the IP ('^') transformation is a 64-bit bit reversal,
//   * ord_2 / com is trailing_zeros of an XOR,
//   * the ring mod 2^N is plain masking (R = 1 << N),
//   * the AR stable point is always 1, so only the counter spn is tracked.
//
// Requires 2 <= N <= 62; the total model frequency must be < 2^(N-1).

#pragma once

#include <cstdint>
#include <vector>
#include <stdexcept>

namespace padic2 {

using u64 = uint64_t;
using u128 = __uint128_t;

// ---------------------------------------------------------------------------
// bit reversal (branch-free swap-mask trick)
// ---------------------------------------------------------------------------

inline u64 bitreverse64(u64 x) {
    x = ((x & 0x5555555555555555ULL) << 1) | ((x >> 1) & 0x5555555555555555ULL);
    x = ((x & 0x3333333333333333ULL) << 2) | ((x >> 2) & 0x3333333333333333ULL);
    x = ((x & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    x = ((x & 0x00FF00FF00FF00FFULL) << 8) | ((x >> 8) & 0x00FF00FF00FF00FFULL);
    x = ((x & 0x0000FFFF0000FFFFULL) << 16) | ((x >> 16) & 0x0000FFFF0000FFFFULL);
    x = (x << 32) | (x >> 32);
    return x;
}

// ---------------------------------------------------------------------------
// model
// ---------------------------------------------------------------------------

struct PadicModel2 {
    int N;
    u64 R;       // 2^N
    u64 mask;    // R - 1
    int shift;   // 64 - N, for bit reversal
    std::vector<u64> cum;
    u64 total;

    PadicModel2(const std::vector<u64>& freqs, int N_) : N(N_) {
        if (freqs.empty()) throw std::invalid_argument("alphabet must not be empty");
        if (N < 2 || N > 62) throw std::invalid_argument("N must be in 2..62");
        cum.reserve(freqs.size() + 2);
        cum.push_back(0);
        for (u64 f : freqs) {
            if (f == 0) throw std::invalid_argument("frequencies must be positive");
            cum.push_back(cum.back() + f);
        }
        cum.push_back(cum.back() + 1);  // EOM symbol, frequency 1
        total = cum.back();
        if (total > (u64(1) << (N - 1)) - 1)
            throw std::invalid_argument("grid level N too small for total frequency");
        R = u64(1) << N;
        mask = R - 1;
        shift = 64 - N;
    }

    int eom() const { return static_cast<int>(cum.size()) - 1; }
};

// ---------------------------------------------------------------------------
// word-level p-adic primitives
// ---------------------------------------------------------------------------

inline u64 hat(const PadicModel2& m, u64 x) { return bitreverse64(x) >> m.shift; }

inline int compath(u64 x, u64 y, int N) {
    if (x == y) return N;
    return __builtin_ctzll(x ^ y);
}

inline u64 cut1(u64 x) { return (x & 0x1) | ((x >> 2) << 1); }

inline u64 resc(const PadicModel2& m, u64 x, int n) { return hat(m, hat(m, x) >> n); }

inline bool ar2(u64 lh, u64 rh) { return (lh & 0x3) == 0x2 && (rh & 0x3) == 0x1; }

// ---------------------------------------------------------------------------
// model on the ring mod 2^N (r == 0 stands for 2^N)
// ---------------------------------------------------------------------------

// m.cum[i] is the right cumulative boundary of symbol i (m.cum[0] == 0), so
// symbol a's interval is [m.cum[a-1], m.cum[a]).
inline void mcode(const PadicModel2& m, int a, u64 l, u64 r, u64& lnew, u64& rnew) {
    u64 R = (r == 0) ? m.R : r;
    u64 w = R - l;
    lnew = l + static_cast<u64>((u128(w) * m.cum[a - 1]) / m.total);
    u64 rn = l + static_cast<u64>((u128(w) * m.cum[a]) / m.total);
    rnew = rn & m.mask;
}

inline int mdecode(const PadicModel2& m, u64 g, u64 l, u64 r, u64& lnew, u64& rnew) {
    u64 R = (r == 0) ? m.R : r;
    u64 w = R - l;
    u64 off = g - l;
    int lo = 1, hi = static_cast<int>(m.cum.size()) - 1;
    while (lo < hi) {
        int mid = (lo + hi + 1) >> 1;
        if (static_cast<u64>((u128(w) * m.cum[mid - 1]) / m.total) <= off) lo = mid;
        else hi = mid - 1;
    }
    mcode(m, lo, l, r, lnew, rnew);
    return lo;
}

// ---------------------------------------------------------------------------
// bit i/o: packed into a vector<uint8_t>, one bit per byte's LSB packed later
// via std::vector<bool>-like bitset; we use a simple growable bit buffer.
// ---------------------------------------------------------------------------

struct BitOut {
    std::vector<u64> words;
    size_t nbits = 0;

    inline void reserve_bits(size_t n) { words.reserve((n >> 6) + 2); }

    inline void push(bool b) {
        size_t w = nbits >> 6, o = nbits & 63;
        if (o == 0) words.push_back(0);
        if (b) words[w] |= (u64(1) << o);
        ++nbits;
    }
    inline void push_low(u64 x, int n) {
        for (int i = 0; i < n; ++i) push((x >> i) & 1);
    }
    inline void push_rep(bool b, int n) {
        for (int i = 0; i < n; ++i) push(b);
    }
};

struct BitIn {
    const std::vector<u64>& words;
    size_t nbits;
    size_t pos = 0;

    BitIn(const std::vector<u64>& w, size_t n) : words(w), nbits(n) {}

    inline u64 get() {
        if (pos >= nbits) return 0;
        size_t w = pos >> 6, o = pos & 63;
        ++pos;
        return (words[w] >> o) & 1;
    }
    inline u64 get_path(int n) {
        u64 x = 0;
        for (int i = 0; i < n; ++i) x = (x << 1) | get();
        return x;
    }
};

// ---------------------------------------------------------------------------
// selectPoint: minimal path into [l, r), built greedily bit by bit
// ---------------------------------------------------------------------------

inline u64 select_path2(const PadicModel2& m, u64 l, u64 r) {
    u64 lo = l;
    u64 hi = (r == 0) ? m.R : r;
    u64 q = 0;
    for (int i = 0; i < m.N; ++i) {
        u64 a0 = lo + (lo & 0x1);
        u64 d = 0;
        if (a0 >= hi) {
            d = 1;
            a0 = lo | 0x1;
        }
        q = (q << 1) | d;
        lo = a0 >> 1;
        hi = (hi - d + 1) >> 1;
    }
    return q;
}

// ---------------------------------------------------------------------------
// encoder / decoder (revised algorithm, PR + AR)
// ---------------------------------------------------------------------------

inline void flush_sp(BitOut& out, const PadicModel2& m, u64& l, u64& r, int& spn) {
    if (spn == 0) return;
    u64 lh = hat(m, l), rh = hat(m, r);
    if ((lh & 0x1) == 0x1) {                 // toLeft?(1, l^)
        out.push(true);
        out.push_rep(false, spn);
    } else if ((rh & 0x1) == 0x0 || rh == 1) { // toRight?(1, r^)
        out.push(false);
        out.push_rep(true, spn);
    } else {
        return;
    }
    l = resc(m, l, 1);
    r = resc(m, r, 1);
    spn = 0;
}

template <typename Msg>
BitOut padic_encode2(const PadicModel2& m, const Msg& msg) {
    BitOut out;
    out.reserve_bits(msg.size() * 12 + 64);  // rough headroom, grows if needed
    u64 l = 0, r = 0;
    int spn = 0;
    for (int a : msg) {
        mcode(m, a, l, r, l, r);
        flush_sp(out, m, l, r, spn);
        if (spn == 0) {
            u64 lh = hat(m, l);
            int n = compath(lh, hat(m, (r - 1) & m.mask), m.N);
            if (n > 0) {
                out.push_low(lh, n);
                l = resc(m, l, n);
                r = resc(m, r, n);
            }
        }
        while (true) {
            u64 lh = hat(m, l), rh = hat(m, r);
            if (!ar2(lh, rh)) break;
            ++spn;
            l = hat(m, cut1(lh));
            r = hat(m, cut1(rh));
        }
    }
    mcode(m, m.eom(), l, r, l, r);
    flush_sp(out, m, l, r, spn);
    if (spn == 0) {
        u64 q = select_path2(m, l, r);
        int n = (q == 0) ? 0 : (64 - __builtin_clzll(q));
        out.push_low(q, n);
    } else {
        out.push(true);
    }
    return out;
}

inline std::vector<int> padic_decode2(const PadicModel2& m, const BitOut& bits) {
    BitIn io(bits.words, bits.nbits);
    std::vector<int> out;
    u64 l = 0, r = 0;
    int spn = 0;
    u64 g = io.get_path(m.N);
    while (true) {
        u64 lnew, rnew;
        int a = mdecode(m, g, l, r, lnew, rnew);
        l = lnew; r = rnew;
        if (a == m.eom()) break;
        out.push_back(a);

        if (spn != 0) {
            u64 lh = hat(m, l), rh = hat(m, r);
            if ((lh & 0x1) == 0x1 || (rh & 0x1) == 0x0 || rh == 1) {
                l = resc(m, l, 1);
                r = resc(m, r, 1);
                g = resc(m, g, 1);
                g += io.get_path(1);
                spn = 0;
            }
        }
        if (spn == 0) {
            int n = compath(hat(m, l), hat(m, (r - 1) & m.mask), m.N);
            if (n > 0) {
                l = resc(m, l, n);
                r = resc(m, r, n);
                g = resc(m, g, n);
                g += io.get_path(n);
            }
        }
        while (true) {
            u64 lh = hat(m, l), rh = hat(m, r);
            if (!ar2(lh, rh)) break;
            ++spn;
            l = hat(m, cut1(lh));
            r = hat(m, cut1(rh));
            g = hat(m, cut1(hat(m, g)));
            g += io.get();
        }
    }
    return out;
}

} // namespace padic2
