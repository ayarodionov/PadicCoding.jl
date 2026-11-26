# PadicCoding.jl demos and experiments

See the [README](../README.md) for an overview, usage and the API. Run the demos
from a clone with `julia --project demos/<name>.jl`.

## Kalman-filter demo

`demos/kalman_demo.jl` codes a synthetic 1-D GPS-like track (20 000 24-bit
samples switching between parked / cruise / maneuver / noisy "urban canyon"
regimes, with rare multipath outliers). Each sample is a residual bin around a
predicted centre (or an escape + raw value), with the per-step distribution
from the adaptive coder:

| Model | bits/sample | vs raw |
|---|---|---|
| static delta histogram (two-pass, table not counted) | 8.07 | 2.97x |
| adaptive delta counter | 7.41 | 3.24x |
| best single Kalman filter (hindsight-tuned) | 7.29 | 3.29x |
| mixture of 6 Kalman filters, likelihood-weighted with forgetting | **5.72** | **4.19x** |

The Kalman predictive Gaussian `N(mu, S)` supplies both the centre and a
calibrated width; mixing filters with different noise levels lets the model
follow regime changes, which no single fixed filter can.

## Binary vs. ternary demo

`demos/ternary_demo.jl` asks whether a ternary machine (digits 0, 1, 2)
compresses better because 3 is closer to `e` than 2. It codes four sources in
both bases and converts the output to bits (1 trit = log2 3 = 1.585 bits):

| Source | arith P=2 | arith P=3 | Huffman D=2 | Huffman D=3 | cost P=2 / P=3 |
|---|---|---|---|---|---|
| uniform {1/3, 1/3, 1/3} | 1.5858 | 1.5858 | 1.6663 | **1.5850** | 95148 / 90045 |
| dyadic {1/2, 1/4, 1/4} | 1.4961 | 1.4962 | **1.4955** | 1.5850 | 89768 / 84960 |
| uniform bytes | 8.0004 | 8.0006 | **8.0000** | 8.0504 | 480024 / 454302 |
| text bytes | 4.0961 | 4.0961 | 4.1392 | 4.1932 | 245764 / 232590 |

Rates are bits/symbol for 30 000 symbols. Information does not depend on the
base: the arithmetic coder reaches the entropy for both `P`. A symbol-by-symbol
Huffman code is better only in the base that matches the probabilities. The
`e` argument is *radix economy*: if a digit cell costs in proportion to its
number of states, storing information costs `b / ln b` per nat, which is
smallest at `b = e`. Ternary is 5.36% cheaper than binary in that model (the
last column, digits × radix), but the compressed output is no smaller.

## DNA storage demo

DNA synthesis and sequencing make many errors on homopolymers (runs of the
same base). Goldman et al. (*Nature* 494, 2013) avoid them completely with a
*rotating code*: data is written as trits, and each trit selects one of the 3
bases that differ from the previous base. The channel is then exactly
ternary, `log2 3 = 1.585` bits per nucleotide (nt). They produce the trits
with a ternary Huffman code over bytes. `demos/dna_demo.jl` uses the same
static byte model for three paths and decodes each one back from the DNA
string:

| Path (150 KB of `/usr/share/dict/words`) | nt | nt/byte | vs. P = 3 | longest run | runs ≥ 4 |
|---|---|---|---|---|---|
| floor: entropy / log2 3 | 393 043 | 2.6203 | | | |
| ternary Huffman + rotating code (Goldman) | 402 700 | 2.6847 | +2.45% | 1 | 0 |
| **P = 3 coder + rotating code** | **393 055** | **2.6204** | 0 | 1 | 0 |
| P = 2 coder, 2 bits per base | 311 488 | 2.0766 | −20.75% | 10 | 3 712 |

- The P = 3 coder writes trits directly and ends 12 nt above the floor. The
  Huffman code spends a whole number of trits on each byte and loses 2.45%
  here (1.6% on the project's own text and source).
- Writing 2 bits per base is 21% denser, but it gives thousands of runs of 4
  or more equal bases. The rotating code never repeats a base.

### Short oligos

In a real system the stream is cut into short *oligos* (oligonucleotides):
single-stranded DNA molecules made by chemical synthesis, base by base.
Errors and cost grow with length, so synthesis gives reliable strands of only
about 100–300 nt. The oligos are made in parallel, kept as an unordered pool
and read back by sequencing, so each one also carries an index and often
redundancy. Goldman et al. used 117-nt oligos: a 100-nt segment of the coded
stream plus index and parity, with consecutive segments overlapping by 75 nt
(each base stored four times).

An arithmetic code makes one substitution corrupt the rest of its block, so
each oligo is coded on its own. With `padic_encode`, every oligo pays for its
EOM symbol (frequency 1 of about 150 000, so `log3 150000 ≈ 11` trits). The
fixed-length coder (`padic_fill`) drops the EOM symbol and instead stores the
number of bytes in a 3–5-trit header. The demo fills oligos with a payload of
`L` nt (index and redundancy not counted) and counts how many each code needs:

| L (nt) | ternary Huffman, padded | P = 3 with EOM | P = 3 fixed length + header |
|---|---|---|---|
| 100 | 4 066 | 4 434 (+9.05%) | 4 073 (+0.17%) |
| 150 | 2 702 | 2 836 (+4.96%) | 2 683 (**−0.70%**) |
| 200 | 2 024 | 2 083 (+2.92%) | 2 011 (**−0.64%**) |
| 300 | 1 347 | 1 362 (+1.11%) | 1 335 (**−0.89%**) |

(150 KB of `/usr/share/dict/words`; the Huffman oligos are padded with a
prefix of the next codeword, so the decoder stops there without a header.)

- Fixed length removes the EOM penalty, which was the whole problem with
  short oligos.
- Against Huffman the result is only about even. The arithmetic coder's
  1.6–2.5% advantage on long streams is mostly used up by the header and the
  1–2 termination trits per oligo: it wins by less than 1% from 150 nt up on
  this text, and loses 0.07–2% on the project's own text and source.

## PAM-3 demo

A PAM-3 link sends one of 3 voltage levels per symbol, at most
`log2 3 = 1.585` bits. Binary data is converted to symbols with small block
codes. `demos/pam3_demo.jl` converts the output of the P = 2 coder with
several block codes and with one whole-message base conversion, and compares
it with the P = 3 coder, which writes symbols directly. Every path is decoded
and checked:

| Path (150 KB of `/usr/share/dict/words`) | symbols | vs. P = 3 |
|---|---|---|
| raw bytes, 3 bits → 2 symbols (no compression) | 800 000 | +103.5% |
| P = 2 coder, 3 bits → 2 symbols (94.6% efficient) | 415 318 | +5.66% |
| P = 2 coder, 11 bits → 7 symbols (99.1%) | 396 445 | +0.86% |
| P = 2 coder, 19 bits → 12 symbols (99.9%) | 393 468 | +0.11% |
| P = 2 coder, whole message converted to base 3 | 393 055 | 0.00% |
| **P = 3 coder, direct** | **393 055** | 0 |

The floor is 393 043 symbols; the project's own text gives the same
percentages.

- Against 3 bits → 2 symbols, the P = 3 coder saves the expected 5.7%
  (`1 − 1.5 / log2 3`).
- Larger block codes recover almost all of it: 19 bits → 12 symbols is only
  0.11% behind. Converting the whole message as one big integer matches the
  P = 3 coder exactly, but it needs the whole message first and big-integer
  arithmetic.
- So the P = 3 coder is not better at compression here. Its advantage is that
  it needs no separate conversion stage: it is a streaming, fixed-precision
  coder that already speaks the channel's alphabet.

### Fixed frames

Links send fixed frames, so the demo also uses the fixed-length coder. For
incompressible data (what links usually carry), the P = 3 coder with a
uniform model over bits works as a base converter of any size. The table
gives the bits that always fit in a frame of `F` symbols:

| F | 3 bits → 2 symbols | 11 → 7 (+ 3 → 2 for the rest) | P = 3 coder | `floor(F log2 3)` |
|---|---|---|---|---|
| 7 | 9 | 11 | 11 | 11 |
| 112 | 168 | 176 | 177 | 177 |
| 161 | 240 | 253 | 255 | 255 |
| 1000 | 1500 | 1571 | 1584 | 1584 |

For `F ≤ 7` all inputs are checked; for larger frames the coder's numbers
come from 2000 random inputs plus constant and alternating ones, so they are
not a guarantee. For compressible text in frames, with a byte-count header in
each frame, the P = 3 coder fits 0.7–1.0% more bytes than the P = 2 coder
followed by 11 → 7, and about 5% more than with 3 → 2.

So on a link the fixed-length coder gains at most about 1% over an 11 → 7
code, and it costs error resilience: one wrong symbol corrupts the rest of the
frame, where 11 → 7 confines the damage to 11 bits. It also needs arithmetic
for every symbol, where a block code is a small lookup table.

## Composite P

The paper takes `P` prime, but the coder only uses base-`P` digits: `hat`
reverses them and `ordp(x - y, P)` counts the trailing digits two indexes
share, which is meaningful for any base `P ≥ 2`. `demos/composite_p_test.jl`
runs 600 random trials for each of 19 values of `P` (2, 3, 4, 5, 6, 7, 8, 9,
10, 12, 16, 27, 30, 64, 100, 210, 256, 1000, 65536), with grid levels that
force hundreds of thousands of AR rescalings:

- every message round-trips for every `P`, prime or not, at the corrected
  width bound (and fails for 0–14% of trials at the paper's bound, for prime
  and composite `P` alike);
- the rate is the same for every `P` (5.1404–5.1407 bits/byte on this README);
- a larger `P` is faster, since each output digit carries more bits: in the
  generic coder (the test's timing run) P = 256 encodes about 15 times faster
  than P = 2. A byte-wise
  range coder is in effect the case `P = 256`;
- `P = 4` is not just `P = 2` with paired bits: rescaling happens at
  different points, and the streams agree for only 40–64% of messages.

## PAM-4 without full swings

A PAM-4 symbol carries 2 bits, but jumps between the outer levels (0 ↔ 3) cause
crosstalk and power noise, so GDDR6X uses MTA ("maximum transition
avoidance") coding to forbid them. That makes PAM-4 a constrained channel,
like DNA without homopolymers, except that the number of choices is not fixed:
3 after level 0 or 3, 4 after level 1 or 2. So no single base `P` fits. The
capacity is `log2((3 + √17)/2) = 1.8325` bits per symbol.

`demos/pam4_mta_demo.jl` runs the p-adic coder backwards, as a
*distribution matcher*: `padic_decode_fixed` with per-step models (one per
previous level, maxentropic probabilities, legal transitions only, no EOM)
turns `k` input bits into `n` levels, and `padic_encode_fixed` recovers the
bits. It works when the final interval of the `n` levels holds only one
`k`-bit point; the demo checks that for every input (all inputs up to
`k = 14`, sampled above). Bits per frame of `n` symbols:

| n | block code valid after any level | block code per previous level (best possible) | p-adic matcher | bits/symbol | vs. 7 bits / 4 symbols |
|---|---|---|---|---|---|
| 4 | 6 | 7 | 7 | 1.750 | 0 |
| 8 | 13 | 14 | 14 | 1.750 | 0 |
| 16 | 28 | 29 | 28 | 1.750 | 0 |
| 32 | 57 | 58 | 58 | 1.813 | +3.6% |
| 64 | 116 | 117 | 116 | 1.813 | +3.6% |
| 128 | 233 | 234 | 234 | 1.828 | +4.5% |

- The matcher reaches the best possible frame code, or comes within one bit,
  with fixed-precision arithmetic and no tables. A block code for `n = 128`
  would need one codebook of `2^234` entries per previous level.
- A stream of 232 000 bits in frames of 64 symbols, with the previous level
  carried across frames, has no full swings; plain Gray-coded PAM-4 has
  14 567.
- GDDR6X's MTA appears to send 7 bits per 4 symbols (1.75 bits/symbol), the
  best possible for `n = 4`. A code valid after any level fits only 6 bits
  there; the MTA patents use the DBI line as extra signalling.
- Long frames buy up to 4.5% more throughput, at the usual price: one wrong
  symbol corrupts the rest of the frame, and the frame adds latency. The
  checks above `k = 14` are sampled, not proven.

### Energy shaping

The demo can also bias the model toward cheap transitions. With a simple
energy model, a transition `p → q` costs `(q − p)^2` (switching energy grows
with the square of the voltage step; an optional static per-level term is off
by default). Weighting the transitions by `2^(−μ E)` gives, for each `μ`, the
Markov chain with the most bits per symbol for its mean energy. The baselines
are block codes that send `k` bits in 4 symbols using, for each previous
level, the `2^k` lowest-energy codewords:

| code | bits/symbol | energy/bit |
|---|---|---|
| plain PAM-4, full swings allowed | 2.000 | 1.246 |
| block 7b4s, low-energy codebooks | 1.750 | 0.721 |
| block 6b4s, low-energy codebooks | 1.500 | 0.527 |
| block 5b4s, low-energy codebooks | 1.250 | 0.408 |
| p-adic, one stream, μ = 0.25 | 1.790 | 0.651 |
| p-adic, one stream, μ = 0.5 | 1.681 | 0.516 |
| p-adic, one stream, μ = 1 | 1.457 | 0.382 |
| p-adic, fixed 128-symbol frames, μ = 0.25 | 1.500 | 0.767 |
| p-adic, fixed 128-symbol frames, μ = 1 | 0.906 | 0.617 |

- As one continuous stream, the shaped matcher reaches the theoretical
  trade-off (e.g. 1.7896 bits/symbol vs. 1.7868 in theory at μ = 0.25) and
  beats the low-energy block codes on both counts: more bits per symbol *and*
  less energy per bit (1.79 vs. 1.75 bits/symbol at 0.651 vs. 0.721 per bit).
  But the number of bits per frame then varies, so the input must be
  buffered.
- With a fixed number of bits per frame, which a link needs, it is *worse*
  than the block codes. Every input must fit, so the rate is set by the most
  probable frame (for example one that stays on a single level), not by the
  average.
- So energy shaping with this coder suits a buffered, variable-rate setting.
  For fixed-rate links, low-energy block codes (possibly generated from the
  same model) are the better choice.

Here, unlike plain compression, the coder brings something a fixed base
cannot: the Markov model handles a constraint whose number of choices changes
from symbol to symbol. `P` itself does not matter (`P = 2` is used).

## Is P ≠ 2 worth it?

For compression ratio, no. Information does not depend on the base, and
every gain measured above is against a weak baseline (3 → 2 packing, Huffman
on long streams). Against a good binary-to-ternary conversion the gain is
below 1%, and in realistic settings (short oligos, link frames with error
handling) about zero. `P ≠ 2` does buy:

1. convenience when the medium is natively `P`-ary (DNA with the rotating
   code, PAM-3): one streaming, fixed-precision coder, no conversion stage;
2. speed, from fewer output steps for large `P`;
3. a unified view: the algorithm works for any base, which ties p-adic
   coding to ordinary range coding.

The most useful result of these experiments was the corrected width bound,
which applies to every `P`, including `P = 2`.

## Performance

`benchmark.jl` on an Apple-silicon Mac, Julia 1.12, 10^6 symbols of a skewed
256-symbol source (~7.15 bits/symbol entropy; output rate matches it for
every coder):

| Coder | Encode | Decode |
|---|---|---|
| `PadicCoding`, PR+AR, P=2, N=40 | 0.2 MB/s | 0.2 MB/s |
| `PadicCoding`, PR+AR, P=3, N=26 | 0.4 MB/s | 0.3 MB/s |
| `PadicCoding`, PR+AR, P=5, N=18 | 0.6 MB/s | 0.5 MB/s |
| `PadicCoding2`, P=2, N=40 | **23 MB/s** | **8.7 MB/s** |

`demos/cpp/padic2.hpp` is a direct C++ port of `PadicCoding2` (same PR+AR
algorithm lowered to `uint64_t` word ops, bit-identical output). Built with
`g++ -O3`, it runs the same 10^6-symbol benchmark at roughly 25 MB/s encode /
9 MB/s decode — modestly faster than the Julia specialization, since Julia's
JIT already compiles this kind of tight bit-manipulation loop to comparably
efficient native code:

```
g++ -O3 -std=c++17 demos/cpp/padic2_bench.cpp -o padic2_bench
./padic2_bench
```
