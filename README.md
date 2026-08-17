# PadicCoding.jl

[![CI](https://github.com/ayarodionov/PadicCoding.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ayarodionov/PadicCoding.jl/actions/workflows/CI.yml)

A Julia package implementing the compression algorithm from

> A. Rodionov and S. Volkov, "p-Adic arithmetic coding", *Contemporary
> Mathematics*, 508, 201–213, 2010.
>
> A. Rodionov and S. Volkov, "p-Adic arithmetic coding", 29 pp., 2007,
> <http://arxiv.org/abs/0704.0834v1>
> (`p-adic_arithm_coding.pdf` in this directory).

The coder is a generalization of integer arithmetic coding. A message is
mapped to a semi interval of `[0, 1)` whose edges lie on an equidistant grid
`G(P^N)` (`P` prime); edge points are kept as integer indexes in the ring
mod `P^N`. The *IP* transformation (operator `^` of the paper) maps a grid
index to the p-adic integer number of the path leading to it on the coding
P-tree — it simply reverses the `N` base-P digits. Whenever the paths to the
two interval edges get a common prefix (their p-adic distance drops below 1),
that prefix is pushed to the output and the interval is rescaled (*PR*
rescaling). The revised algorithm adds *AR* rescaling, based on the ordinary
Archimedean distance, for intervals straddling a level-1 grid point; together
the two rescalings keep the interval width above `2/P^2 - 1/P^N`, so the coder
works with fixed-precision integers on arbitrarily long messages. For `P = 2`
the algorithm is an integer variant of classic arithmetic coding; for special
models it reproduces Huffman and Golomb-Rice codes exactly.

## Files

| File | Contents |
|---|---|
| `src/PadicCoding.jl` | Reference implementation for any prime `P`: module `PadicCoding` |
| `src/PadicCoding2.jl` | Optimized `P = 2` specialization: submodule `PadicCoding2` |
| `test/runtests.jl` | Test suite (paper examples, roundtrips, cross-checks) |
| `demos/benchmark.jl` | Throughput benchmark of all coders |
| `demos/demo.jl` | Compresses a text file and prints the rate vs. entropy |
| `demos/ip_demo.jl` | IPv4 address list compression (flat vs. per-octet-position model) |
| `demos/ipv6_demo.jl` | IPv6 address list compression (flat vs. per-byte-position model) |
| `demos/ipv6_split_demo.jl` | IPv6 compression with a dictionary prefix model + per-field coding |
| `p-adic_arithm_coding.pdf` | The paper |

No external packages are required (only `Test`, `Random`, `Printf` from the
standard library).

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/ayarodionov/PadicCoding.jl")
```

or, for local development from a clone:

```julia
Pkg.develop(path = "path/to/PadicCoding.jl")
```

## Usage

```julia
using PadicCoding

# Static model: alphabet 1..4 with weights 1/2, 1/4, 1/8, 1/8
# (an EOM symbol with frequency 1 is appended automatically).
m = PadicModel([4, 2, 1, 1]; P = 2, N = 32)

msg  = [1, 2, 1, 1, 3, 4, 1]          # symbols are 1-based indexes
bits = padic_encode(m, msg)           # Vector{Int} of P-bits (0 .. P-1)
@assert padic_decode(m, bits) == msg
```

The optimized binary coder has the same interface and produces bit-identical
output for the same frequencies and `N` (exported by `using PadicCoding` as
well):

```julia
m2   = PadicModel2([4, 2, 1, 1]; N = 32)
bits = padic_encode2(m2, msg)         # packed BitVector
@assert padic_decode2(m2, bits) == msg
```

Run the tests, demo and benchmark from a clone with

```
julia --project -e 'using Pkg; Pkg.test()'
julia --project demos/demo.jl [file]
julia --project demos/benchmark.jl [file]
julia --project demos/ip_demo.jl
julia --project demos/ipv6_demo.jl
julia --project demos/ipv6_split_demo.jl
```

## API

### `PadicCoding` (any prime `P`)

- `PadicModel(freqs; P = 2, N)` — static frequency model over the alphabet
  `1 .. length(freqs)`; appends an EOM symbol with frequency 1. The total
  frequency must not exceed `2*P^(N-2) - 1` (the width guarantee of the
  revised algorithm). Internally all indexes are `Int128`.
- `padic_encode(m, msg)` / `padic_decode(m, bits)` — revised algorithm with
  PR + AR rescaling (pseudocode on pp. 25–27 of the paper). This is the pair
  to use.
- `padic_encode_pr(m, msg)` / `padic_decode_pr(m, bits)` — basic algorithm
  with PR rescaling only (pp. 17–19), kept for study and comparison.
- `eom(m)` — index of the EOM symbol.

### `PadicCoding2` (optimized, `P = 2` only)

- `PadicModel2(freqs; N)` with `2 ≤ N ≤ 62`.
- `padic_encode2(m, msg) -> BitVector`, `padic_decode2(m, bits)`.

The specialization lowers the p-adic machinery to word operations on
`UInt64`: `^` is a hardware `bitreverse`, `ord_2`/`com` is `trailing_zeros`
of an XOR, `cut`/`res`/`lift` are shifts and masks, the ring mod `2^N` is
plain masking, and the AR stable point is hardwired to 1.

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

## Implementation notes

- **`selectPoint`** (the shortest-path point of the final interval) is found
  in `O(N·P)` instead of scanning the interval: minimizing the path as a
  p-adic integer is equivalent to minimizing the index digits
  lexicographically from the least significant one up, which is done greedily
  digit by digit.
- **Trailing zeros** of the final point are dropped, as the paper prescribes
  (`lnz`); the decoder's bit reader pads an exhausted stream with zeros.
- **A deviation from the literal pseudocode:** in the encoder's AR flush the
  `toLeft?`/`toRight?` branches are `if`/`elseif`, not two independent `if`s.
  After the first branch fires and resets `sp = 0`, the literal second `if`
  can trigger spuriously when `r = 0` (since `r^ == sp` becomes `0 == 0`) and
  would emit an invalid `sp - 1 = -1` digit.
- **PR-only coder and straddling intervals:** the paper notes that a message
  interval straddling a level-1 point defeats PR rescaling (nothing is ever
  pushed). In exact arithmetic the interval would shrink forever; with integer
  truncation the straddle eventually breaks and the PR-only coder survives
  with a small compression loss. It can still fail (`message interval too
  narrow`) for unlucky model/message combinations — the PR+AR coder cannot.
- Both `padic_encode` and `padic_encode2` implement the same algorithm and
  their outputs are bit-for-bit identical; the test suite cross-decodes
  streams between the two implementations.
