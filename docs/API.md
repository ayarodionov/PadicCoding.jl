# PadicCoding.jl API

See the [README](../README.md) for an overview, usage examples and the demos.

## `PadicCoding` (any `P ≥ 2`)

- `PadicModel(freqs; P = 2, N, eom = true)` — static frequency model over the
  alphabet `1 .. length(freqs)`; appends an EOM symbol with frequency 1
  (`eom = false`: no EOM symbol, for the fixed-length functions only). `P` may be
  any integer `≥ 2`, prime or not. The total frequency must not exceed
  `P^(N-2) + 1` (the width guarantee of the revised algorithm, see the
  [implementation notes](../README.md#implementation-notes)). Internally all indexes are `Int128`, so
  `P^N * total` must stay below `2^126`.
- `padic_encode(m, msg)` / `padic_decode(m, bits)` — revised algorithm with
  PR + AR rescaling (pseudocode on pp. 25–27 of the [paper](p-adic_arithm_coding.pdf)). This is the pair
  to use.
- `padic_encode_fixed(m, msg)` / `padic_decode_fixed(m, bits, n)` — the same
  coder for messages whose length `n` the decoder knows: no EOM symbol is
  coded, which saves about `log_P(total)` digits per message (about 11 trits
  for a byte model built from 150 KB of text). Useful for many short blocks.
- `padic_encode_fixed(model_at, msg)` / `padic_decode_fixed(model_at, bits, n)`
  — per-step models: symbol `t` uses `model_at(t, prefix)`, called once per
  symbol in order; all models must share `P` and `N`.
- `padic_fill(m, msg, L; start = 1) -> (k, digits)` — fixed code length: codes
  the longest prefix `msg[start : start+k-1]` that fits in `L` digits (for
  example one DNA oligo, a short synthetic DNA strand of about 100–300
  nucleotides, or one link frame); padding `digits` with zeros to `L`
  does not change the result. Decode with `padic_decode_fixed(m, digits, k)`,
  so `k` must be sent too, e.g. in a small header.
- `padic_encode_pr(m, msg)` / `padic_decode_pr(m, bits)` — basic algorithm
  with PR rescaling only (pp. 17–19 of the paper), kept for study and comparison.
- `eom(m)` — index of the EOM symbol.

## `PadicCoding2` (optimized, `P = 2` only)

- `PadicModel2(freqs; N)` with `2 ≤ N ≤ 62` and total frequency at most
  `2^(N-2) + 1`.
- `padic_encode2(m, msg) -> BitVector`, `padic_decode2(m, bits)`.
- `padic_encode2(model_at, msg)`, `padic_decode2(model_at, bits)` — adaptive
  coding: symbol `t` uses the model `model_at(t, prefix)`, where `prefix` holds
  the symbols before `t` (the EOM uses `t = length(msg) + 1`). `model_at` is
  called once per step, in order, so it may be a stateful closure (an adaptive
  counter, a Kalman filter, a neural predictor, ...); the decoder must use a
  fresh, identical copy of it. All models must share `N`.

The specialization lowers the p-adic machinery to word operations on
`UInt64`: `^` is a hardware `bitreverse`, `ord_2`/`com` is `trailing_zeros`
of an XOR, `cut`/`res`/`lift` are shifts and masks, the ring mod `2^N` is
plain masking, and the AR stable point is hardwired to 1.
