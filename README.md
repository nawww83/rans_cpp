# rans_cpp

Byte-wise range ANS codec implementation in C++.

A compact and fast entropy coder designed for high-throughput compression and decompression workloads. It emphasizes efficient integer arithmetic, low-overhead symbol handling, and a tunable speed/compression trade-off.

## Why rans_cpp?

- High throughput on modern CPUs
- Byte-oriented rANS design
- Optimized for 32-bit arithmetic and compact frequency tables
- Fast four-lane symbol counting
- Selectable frequency quantization: from exact (arithmetic-coding-like) to Huffman-like
- Self-describing stream: the decoder does not need to be told the encoder mode

## Features

- Quad interleaved rANS
- 16-bit flushing and 32-bit arithmetic
- Direct symbol-indexed encoder tables (no per-block sort)
- Compact bit-packed frequency table
- Automatic scale selection via parameter `L` (power of two, `256 <= L <= 65536`)
- Free-remainder rANS transformation using only `div` (no `mod` required)
- Four-lane symbol counter that breaks read-modify-write dependencies
- Four frequency-normalization modes with specialized encoder paths

## Quantization modes

The encoder can normalize symbol frequencies in four ways, trading compression
ratio for speed:

| Mode | Frequencies | Encoder division | Notes |
|------|-------------|------------------|-------|
| `Precise` | arbitrary integers | hardware `div` | Closest to arithmetic coding; best CR, slowest |
| `Pow2` | `2^a` | shift `x >> a` | Equivalent to Huffman code lengths; fastest, worst CR |
| `Pow23` | `2^a * 3^b` | shift + reciprocal multiply | Intermediate |
| `Pow235` | `2^a * 3^b * 5^c` | shift + reciprocal multiply | Closest quantized mode to `Precise` |

All modes keep the invariant `sum(freq) == L` and `freq >= 1` for every symbol
present in the input, where `L` is a power of two. In the code the modes are
`Precise`, `Pow2`, `Pow23` and `Pow235`; the benchmark prints the last two as
`Pow{2,3}` and `Pow{2,3,5}`, and this document uses both spellings.

### How the quantized frequencies are built

- `Precise` — proportional normalization with a minimum frequency of 1.
- `Pow2` — length-limited Huffman via **package-merge** with maximum code
  length `log2(L)`, so `f = 2^a`. This is the optimal prefix-code solution
  under the Kraft equality `sum(2^-len) == 1`.
- `Pow23` / `Pow235` — two-pass **water-filling**:
  1. For every symbol take the largest smooth value `<= t_i`, where
     `t_i = count_i * L / N` (the lower bracket).
  2. Distribute the remaining `L - sum(f)` one smooth step at a time, always
     choosing the symbol with the largest marginal gain per unit of frequency.
     A final sum-preserving pair replacement fixes the last 0–1 units.
  This keeps `sum(f) == L` exactly and gives near-`Precise` redundancy at a
  higher encoding speed.

### Encoder division paths

| Mode | Hot-loop division |
|------|-------------------|
| `Precise` | `x / f` (hardware `div`) |
| `Pow2` | `x >> a` |
| `Pow23` / `Pow235` | `x >> a`, then multiplication by `ceil(2^48 / odd)` and `>> 48`; exact because `y < odd * Bflush`, so no correction step is needed |

## Stream header

The header stores the frequency table compactly. The mode is written into the
stream, so the caller does not pass it to the decoder.

| Field | Size | Meaning |
|-------|------|---------|
| `M - 1` | 1 byte | number of distinct symbols |
| `mode \| logL` | 1 byte | mode (2 bits) and `logL` (5 bits), `L = 1 << logL` |
| symbol records | bit-packed | see below |
| `x, y, z, w` | 4 × u32 | final rANS states |
| `in_size` | u32 | decompressed size |

Each symbol record is bit-packed (LSB-first), with 8 bits for the symbol and:

| Mode | Payload | Payload bits | Total bits/symbol |
|------|---------|-------------:|------------------:|
| `Precise` | `f - 1` | 16 | 24 |
| `Pow2` | `a` | 5 | 13 |
| `Pow23` | `a`, `b` | 9 | 17 |
| `Pow235` | `a`, `b`, `c` | 12 | 20 |

`Total bits/symbol` includes the 8-bit symbol. Compared with storing `f - 1`
for every mode, this saves up to 11 bits per symbol (352 bytes for a
256-symbol alphabet). With `mode == Precise` the stream layout is byte-for-byte
identical to the older 16-bit frequency table.

## Benchmark

`main.cpp` benchmarks all four modes on identical data.

### Data distribution

Each byte is generated as

```text
k ~ Geometric(p),   byte = (k mod 256) XOR (inversion ? 0xFF : 0x00)
```

where `k` is the number of failures before the first success,
`P(k) = p (1-p)^k`. Keeping the low 8 bits is exact, not a lossy shortcut:
the folded tail is a geometric series, so

```text
P(byte = b) = p (1-p)^b / (1 - (1-p)^256),   b = 0..255,
```

i.e. the geometric distribution **truncated to the 256-symbol alphabet** and
renormalized. For small `p` (`p << 1/256`) this is almost uniform
(`H ~ 8` bits); for `p -> 1` it concentrates on `0x00` (`H -> 0`). The sweep
`p = max(0.0001, q/128)`, `q = 0..128`, covers the whole range.

The `inversion` flag XORs the bytes with `0xFF`; it mirrors the alphabet
(`0x00 <-> 0xFF`) and checks that the codec does not depend on symbol order.

`N` is drawn uniformly from `[1, N_max]` (`N_max = 8 MB` by default), so the
sweep also exercises header overhead on small inputs. For every `p` the same
buffer is encoded and decoded by all four modes, and each round-trip is
verified.

### Fixed-p runs

The main table aggregates the 129 `p` values (min/median). Because a single
median hides the behavior at the extremes, `main.cpp` also prints a
`Fixed-p benchmarks` table (for `N >= N_max / 4`) with

```text
p = 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.6, 0.9
```

and the same `N` and inversion flag. It shows, for example, that `Pow2` stays
within ~0.7 percentage points of `Precise` up to `p = 0.3` but loses badly at
`p = 0.9` (very skewed data), while `Pow23`/`Pow235` track `Precise` across the
whole range.

### Redundancy

For each buffer the benchmark computes the normalized histogram `f[256]`
(sum `L`, `L` a power of two) and its entropy

```text
H = -sum_i (f[i]/L) * log2(f[i]/L)      [bits per byte]
```

The ideal compression ratio is `ideal_CR = 8 / H`, the actual one is
`CR = N / compressed_size`, and

```text
redundancy = (ideal_CR / CR - 1) * 100   [%]
```

`compressed_size` includes the rANS header (packed frequencies, `L`, original
size) and the final encoder state, so small-input runs measure the header cost
as well as the coding cost. Runs with `H <= 0.001` are reported as 0%
(entropy is too small for a meaningful ratio).

Reported per mode:

- `CR min` / `CR med` — compression ratio `N / compressed_size`
- `Red% min` / `Red% med` — redundancy against the ideal entropy,
  `(ideal_CR / CR - 1) * 100`
- `Enc MB/s` — encode throughput
- `Dec MB/s` — decode throughput

Throughput is sampled for `N >= N_max / 4`. `N_max` is set in `main.cpp`
(8 MB for the slow test; smaller values are provided as comments).

## Build

The code builds with GCC, Clang and MSVC (C++20). On MSVC the project adds
`/utf-8`, `/Zc:__cplusplus` and `/permissive-` automatically.

### CMake (recommended)

GCC/Clang:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Visual Studio:

```bash
cmake -S . -B build -G "Visual Studio 18 2026" -A x64
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The benchmark is built as `rans`, the tests as `rans_tests`.
Enable CPU-specific tuning with `-DRANS_NATIVE=ON`
(`-march=native` on GCC/Clang, `/arch:AVX2` on MSVC).

### Direct compilation

GCC/Clang:

```bash
g++ main.cpp -O3 -std=c++20 -o rans
```

MSVC:

```bat
cl /std:c++20 /EHsc /utf-8 /Zc:__cplusplus /permissive- /O2 main.cpp /Fe:rans.exe
```

## Usage

```bash
./rans
```

The project is intended as a C++ implementation of byte-wise range ANS
compression and decompression routines, suitable for benchmarking and
experimentation.

## Tests

The tests are registered with CTest, so the CMake workflow above is enough:

```bash
ctest --test-dir build --output-on-failure
```

They can also be built directly:

```bash
g++ test_quantizer.cpp -O2 -std=c++20 -o test_quantizer
./test_quantizer
```

Asserts are intentionally kept enabled for `rans_tests` even in Release builds,
because the test suite is assert-based.

The test covers:

- exactness of the reciprocal-division path on every smooth divisor,
- `sum(f) == L` and membership in the allowed set for each mode,
- round-trip encode/decode on random and degenerate inputs (1 symbol,
  2 symbols, tiny buffers, full 256-symbol alphabet).

## Performance

Absolute throughput depends on CPU, data distribution and frequency scaling.
Representative medians from the distribution sweep on Linux (Intel Core
i7-12700, GCC 11.4, pinned to one core with `taskset -c 3`, inputs 2-8 MB):

| Mode | Redundancy | Enc MB/s | Dec MB/s |
|------|-----------:|---------:|---------:|
| `Precise` | ~0.15 % | ~410 | ~499 |
| `Pow235` | ~0.16 % | ~462 | ~498 |
| `Pow23` | ~0.20 % | ~461 | ~500 |
| `Pow2` | ~1.79 % | ~472 | ~512 |

`Pow2` is fastest but least compact. `Precise` and `Pow235` are essentially
tied on redundancy, but `Pow235` encodes ~13% faster because it replaces
hardware division with a shift and an exact magic reciprocal. `Pow23` is
slightly less compact and about as fast. Decompression is essentially
mode-independent, since every mode decodes with the same multiplication-based
rANS step: the four modes land within ~3% of each other.

The encoder precomputes `fr * Bflush`, `L - fr`, and the 48-bit magic
reciprocal, so the non-power-of-two modes need neither a division nor a
correction step in the hot loop. The remaining limit is the dependent rANS
state update and renormalization, not the division alone: per `p` the
water-filling modes encode 6-29% faster than `Precise` (`Pow2` up to 53%), and
12-15% faster on the sweep medians — not the 2-3x that removing `div` from the
hot loop suggests on paper. The gain is largest where the alphabet is small and
the stream is skewed (`p` around 0.001 or 0.9), and smallest around
`p = 0.03-0.1`.

The fixed-`p` runs show where `Pow2` breaks. Up to `p = 0.3` it stays within
0.7 percentage points of `Precise` (`p = 0.3`: CR 2.702 against 2.720), but at
`p = 0.9` one symbol dominates and integer code lengths cannot express a
frequency ratio of thousands to one: CR 7.200 against 15.326, i.e. 1.111
bits/symbol where the entropy is 0.521. `Pow23` and `Pow235` track `Precise`
across the whole range (CR 15.323 at `p = 0.9`). At the other extreme, near
`p = 0.0001` all modes compress the same (CR 1.000), but encode at 635
(`Precise`), 775 (`Pow23`), 771 (`Pow235`) and 895 (`Pow2`) MB/s.

Two changes account for most of the encode gain in this version. Measured by
rebuilding the previous revision and running both builds back to back on the
same core with the same inputs (2 MiB), they are worth 1.04-1.12x at
`p <= 0.1`, 1.11-1.18x at `p = 0.3`, 1.30-1.39x at `p = 0.6` and 1.58-1.78x
at `p = 0.9`.

First, byte counting uses four independent histograms instead of one: at
`p = 0.9` a single histogram serializes four read-modify-write operations on
the same element, while the four-lane version runs them in parallel. The
counting pass itself becomes 3.3x faster at `p = 0.6` and 3.8x at `p = 0.9`,
which alone is +25-30% encode at `p = 0.6` and +54-62% at `p = 0.9`, +8-10% at
`p = 0.3`, and +1-4% at low `p`.
Second, the encoder tables are indexed directly by symbol, so the hot loop no
longer loads `symbol_to_idx[symbol]` and the per-block frequency sort is gone:
another +2-10% across the sweep (median ~+8%). The header is now written in
symbol order; compressed size differs from the sorted order only by a few
bytes per block.

Decoding is unchanged and remains latency-bound. The 64 KiB `mLut` does not fit
in L1, but it is L2-resident and four independent streams hide most of that
latency: clamping `Lmax` to 32 KiB or 16 KiB makes the table L1-resident and
speeds decode up by 3-9% on 2 MiB blocks (439 -> 461 and 678 -> 738 MB/s at
`p = 0.1` and `p = 0.9`), so the cache limit is real, but shrinking `L` is the
wrong trade in general: `choose_L()` deliberately picks `L > N / min_count` so
that no observed symbol is rounded to zero, and a clamp would force the
dynamic-range fixup on short blocks with a wide alphabet.

Compiling with `-DRANS_NATIVE=ON` (`-march=native`) does not pay off on this
CPU: best-of-25 decode came out at 424-430 MB/s against 432-447 MB/s for the
generic build, while encode moved by at most ~4% either way. The hot loops are
scalar integer code — a byte LUT lookup plus dependent multiply and shift
chains — so there is little for AVX2 to vectorize. The option is kept for
machines where the compiler does find something.

Fixed per-block costs matter for small buffers. At 64 KiB and `p = 0.1` the
histogram takes 15 µs, the quantizer 13 µs (`Pow2`), 34 µs (`Pow23`) or 21 µs
(`Pow235`), the `Precise` path counts and normalizes in 16 µs, encode takes
138-172 µs and decode ~125 µs. On such blocks the `Pow23` quantizer alone is
20% of encode time, and `Pow23`/`Pow235` end up *slower* to encode than
`Precise` (172 and 160 µs against ~150): water-filling costs 21-34 µs on top of
the histogram, while the `Precise` normalization adds only about a microsecond
to the same count. At 2 MiB the same quantizers take 18-46 µs against ~5 ms of
encode (under 1%), while the histogram is ~490 µs (~9%, and it is already
included in the throughput table). The decoder also fills up to 64 KiB of
`mLut` per block (~1-2 µs), which is negligible for large blocks but visible
for tiny ones.

## Choosing a mode

The mode is stored in the stream, so a block can be encoded with whichever mode
fits it and the decoder does not need to be told. Redundancy figures below are
the sweep medians from the table above (the definition is in the Redundancy
section); the per-symbol costs come from separate runs at `N` = 64 KiB.

| Situation | Mode | Why |
|-----------|------|-----|
| Default, unknown data | `Pow235` | `Precise`-level ratio, ~+13% encode, 20 header bits/symbol |
| Very small blocks, ratio is all that matters | `Pow23` | 3 header bits/symbol less than `Pow235` |
| Near-uniform data, or encode speed is everything | `Pow2` | Fastest encoder, 13 header bits/symbol |
| Large blocks, ratio is all that matters | `Precise` | Best ratio when `L` is large |

In more detail:

- **`Pow235` is the sensible default.** Its median redundancy (0.16%) is within
  a hundredth of a percentage point of `Precise` (0.15%), and it encodes about
  13% faster on the sweep and 21-29% faster at the extremes. Its quantizer is
  also the cheapest of the two water-filling modes (21 µs against 34 µs at
  64 KiB), and its header is 20 bits per symbol.
- **`Pow23` buys exactly one thing: a shorter header.** Three bits per symbol
  less means 96 bytes per block at a full 256-symbol alphabet. Everything else
  is slightly worse — redundancy 0.20% against 0.16%, and the most expensive
  quantizer — so it only wins when the header dominates, i.e. on very small
  blocks, and only if encode time does not matter there. Otherwise prefer
  `Pow235`.
- **`Pow2` is not a general-purpose choice.** It is a Huffman code: its
  frequencies are integer powers of two, so it cannot express a dominant
  symbol. Up to `p = 0.3` (entropy above ~2.9 bits/byte) it stays within 0.7
  percentage points of `Precise` and encodes 10-53% faster; at `p = 0.9` it
  compresses to CR 7.2 instead of 15.3. Use it for data that is close to
  uniform (already-compressed or encrypted-looking input, byte planes of
  numeric arrays) or when the encoder is the bottleneck and the ratio is not.
- **`Precise` is not automatically the best ratio.** When `L` is small — short
  blocks, or a wide alphabet with very rare symbols — its cumulative rounding is
  worse than water-filling. At `N` = 64 KiB, `p = 0.003` the codec picks
  `L` = 512, and there `Precise` costs +0.042 bit/symbol against +0.020 for
  `Pow23`; given `L` = 4096 on the same data the two swap places (+0.0005
  against +0.005). For large blocks with a large `L` it is the most accurate
  mode and wins.
- **Decoding does not constrain the choice**: all four modes decode within ~3%
  of each other, because the decoder only sees the frequencies, not the mode.

Rules of thumb:

- Start with `Pow235`; there is rarely a reason to switch to `Precise` unless
  the blocks are large and every 0.01% of ratio counts.
- If a single byte value takes more than about half of the block, do not use
  `Pow2`: that is exactly the case where a power-of-two frequency grid fails.
- Below ~64 KiB the water-filling (21-34 µs) is a visible share of encode time
  and `Precise`'s rounding starts to hurt, so `Pow235` is the best compromise;
  `Pow23` only if the block is so small that 96 header bytes really matter.

## Notes

This implementation is optimized for speed and is designed for experimentation
and performance-oriented use cases. The code focuses on efficient rANS
execution while keeping the design compact and approachable.
