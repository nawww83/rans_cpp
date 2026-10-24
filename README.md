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
- Division-free rANS hot loop: shift plus an exact 48-bit reciprocal multiply (no `div` or `mod`)
- Four-lane symbol counter that breaks read-modify-write dependencies
- Four frequency-normalization modes sharing two specialized encoder paths

## Quantization modes

The encoder can normalize symbol frequencies in four ways, trading compression
ratio for speed:

| Mode | Frequencies | Encoder division | Notes |
|------|-------------|------------------|-------|
| `Precise` | arbitrary integers | shift + reciprocal multiply | Closest to arithmetic coding; best CR on large blocks |
| `Pow2` | `2^a` | shift `x >> a` | Equivalent to Huffman code lengths; fastest, worst CR |
| `Pow23` | `2^a * 3^b` | shift + reciprocal multiply | Intermediate |
| `Pow235` | `2^a * 3^b * 5^c` | shift + reciprocal multiply | Closest quantized mode to `Precise`; shorter header |

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

Every frequency is `f = odd * 2^tz` with `odd` odd, so the same exact path
serves all non-`Pow2` modes:

| Mode | Hot-loop division |
|------|-------------------|
| `Precise` / `Pow23` / `Pow235` | `x >> tz`, then multiplication by `ceil(2^48 / odd)` and `>> 48`; exact because `y < odd * Bflush`, so no correction step is needed |
| `Pow2` | `x >> a` (here `odd == 1`, so the reciprocal is the identity; the multiply is skipped) |

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

- exactness of the reciprocal-division path for every odd divisor below 2^16,
  over the full renormalization range `y < odd * Bflush`,
- `sum(f) == L` and membership in the allowed set for each mode,
- round-trip encode/decode on random and degenerate inputs (1 symbol,
  2 symbols, tiny buffers, full 256-symbol alphabet).

## Performance

Absolute throughput depends on CPU, data distribution and frequency scaling.
The numbers below were measured on an Intel Core i7-8565U (Whiskey Lake, 4C/8T,
15 W), GCC 13.3, pinned to one core with `taskset -c 2`. This machine has no
usable PMU and no fixed CPU frequency, so wall-clock runs vary by about ±5%;
the figures are medians over three interleaved rounds. Inputs are 2 MiB (the
same 129-distribution sweep as `main.cpp`, with `N` pinned to 2 MiB instead of
drawn randomly), and each row aggregates that sweep:

| Mode | Redundancy | Enc MB/s | Dec MB/s |
|------|-----------:|---------:|---------:|
| `Precise` | ~0.12 % | ~256 | ~331 |
| `Pow235` | ~0.12 % | ~259 | ~332 |
| `Pow23` | ~0.15 % | ~259 | ~334 |
| `Pow2` | ~1.78 % | ~264 | ~338 |

All non-`Pow2` modes now encode at essentially the same speed. The encoder
uses one exact shift-plus-reciprocal path for any frequency `f = odd * 2^tz`,
so `Precise` no longer pays for hardware division: on the sweep it moved from
~197 to ~256 MB/s (+30%), and from ~270 to ~371 MB/s at `p = 0.9` (+37%).
`Pow2` remains the fastest encoder but is the least compact. Decompression is
mode-independent within noise, since every mode decodes with the same
multiplication-based rANS step.

The remaining limit is the dependent rANS state update and renormalization,
not the division: even with the reciprocal path, encode is far from the 2-3x
that removing `div` suggests on paper. The reciprocal gain is largest where the
stream is skewed (`p` around 0.9) and smallest around `p = 0.03-0.1`.

At 64 KiB the medians are 0.66% (`Precise`), 0.58% (`Pow235`), 0.57%
(`Pow23`) and 2.07% (`Pow2`) redundancy; encode lands at 246-249 MB/s for all
modes there, so the mode choice on small blocks is about ratio and header size,
not speed.

The fixed-`p` runs show where `Pow2` breaks. Up to `p = 0.3` it stays within
0.7 percentage points of `Precise` (`p = 0.3`: CR 2.702 against 2.720), but at
`p = 0.9` one symbol dominates and integer code lengths cannot express a
frequency ratio of thousands to one: CR 7.200 against 15.326, i.e. 1.111
bits/symbol where the entropy is 0.521. `Pow23` and `Pow235` track `Precise`
across the whole range (CR 15.323 at `p = 0.9`). At the other extreme, near
`p = 0.0001`, all modes compress the same (CR 1.000); on this CPU they encode
at roughly 350-420 MB/s there, and the spread is within run-to-run noise.

This revision adds two changes. First, `Precise` uses the 48-bit reciprocal
described above; measured by rebuilding the previous revision and running both
builds back to back on the same core (2 MiB), that is +30% on the sweep and
+37% at `p = 0.9`. Second, `SmoothSet` caches `log2(v)` and water-filling
keeps the current smooth index per symbol instead of searching the 128 KiB LUT
on every step, making the `Pow23`/`Pow235` quantizer 4-5x faster (530 -> 105 µs
at `p = 0.01`, 81 -> 19 µs at `p = 0.1`). That is under 1% of encode at 2 MiB,
but visible on 64 KiB blocks.

The previous revision's gains are unchanged: measured the same way, four-lane
counting and symbol-indexed tables are worth 1.04-1.12x at `p <= 0.1`,
1.11-1.18x at `p = 0.3`, 1.30-1.39x at `p = 0.6` and 1.58-1.78x at `p = 0.9`.

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

Decoding is latency-bound and its throughput is unchanged within noise.
Making the `DecoderInfo` fields `u32` removes zero-extends (about 9% fewer
instructions in a callgrind run) and fixes the `f == 65536` truncation, but
wall-clock decode stays inside the ±5% run-to-run spread. The 64 KiB `mLut`
does not fit in L1, but it is L2-resident and four independent streams hide
most of that latency: clamping `Lmax` to 32 KiB or 16 KiB makes the table
L1-resident and speeds decode up by a few percent on 2 MiB blocks, so the
cache limit is real. Shrinking `L` is the wrong trade in general:
`choose_L()` deliberately picks `L > N / min_count` so that no observed
symbol is rounded to zero, and a clamp would force the dynamic-range fixup on
short blocks with a wide alphabet.

Compiling with `-DRANS_NATIVE=ON` (`-march=native`) does not pay off on this
CPU: best-of-25 decode came out at 424-430 MB/s against 432-447 MB/s for the
generic build, while encode moved by at most ~4% either way. The hot loops are
scalar integer code — a byte LUT lookup plus dependent multiply and shift
chains — so there is little for AVX2 to vectorize. The option is kept for
machines where the compiler does find something.

Fixed per-block costs matter for small buffers. At 64 KiB and `p = 0.1` the
histogram takes ~15 µs, the `Pow23`/`Pow235` quantizer ~10-12 µs, and encode
takes ~260 µs, so the quantizer is a few percent of encode; the `Precise` path
counts and normalizes in about the same time as the histogram. The worst case
for water-filling is near-uniform data with a small `L`: at `p = 0.01` the
quantizer still takes ~105 µs even on 64 KiB blocks, roughly a third of encode,
while `Precise` and `Pow2` do not pay that cost. At 2 MiB the quantizer is
under 1% of encode, while the histogram is ~490 µs (~6%, and it is already
included in the throughput table). The decoder also fills up to 64 KiB of
`mLut` per block (~1-2 µs), which is negligible for large blocks but visible
for tiny ones.

## Choosing a mode

The mode is stored in the stream, so a block can be encoded with whichever mode
fits it and the decoder does not need to be told. Redundancy figures below are
the sweep medians from the table above (2 MiB); the 64 KiB figures are quoted
where they change the recommendation. The definition of redundancy is in the
Redundancy section.

| Situation | Mode | Why |
|-----------|------|-----|
| Large blocks (`L` large) | `Precise` | Best ratio; same encode speed as the water-filling modes |
| Small blocks, ratio matters | `Pow235` | Better ratio than `Precise` when `L` is small, 4 header bits/symbol less |
| Very small blocks, header dominates | `Pow23` | 3 header bits/symbol less than `Pow235` |
| Near-uniform data, or encode speed is everything | `Pow2` | Fastest encoder, 13 header bits/symbol |

In more detail:

- **`Precise` is the default for large blocks.** At 2 MiB its median
  redundancy is 0.12%, tied with `Pow235`, and since it now uses the same
  reciprocal path it encodes at the same speed (256 against 259 MB/s, within
  noise). Its header is 24 bits/symbol, four more than `Pow235`, which is
  negligible on large blocks.
- **`Pow235` is the better choice when `L` is small.** At 64 KiB its median
  redundancy is 0.58% against 0.66% for `Precise`: water-filling handles the
  rounding better when `L` is small, and its header is 4 bits/symbol shorter.
  Its quantizer costs ~10-12 µs at `p = 0.1` on a 64 KiB block.
- **`Pow23` buys exactly one thing: a shorter header.** Three bits per symbol
  less means 96 bytes per block at a full 256-symbol alphabet. Everything else
  is slightly worse or tied — 0.57% redundancy against 0.58% for `Pow235` at
  64 KiB — so it only wins when the header dominates, i.e. on very small
  blocks. Otherwise prefer `Pow235`.
- **`Pow2` is not a general-purpose choice.** It is a Huffman code: its
  frequencies are integer powers of two, so it cannot express a dominant
  symbol. Up to `p = 0.3` (entropy above ~2.9 bits/byte) it stays within 0.7
  percentage points of `Precise`, and its encode speed is now essentially the
  same as the other modes (264 against 256 MB/s on the sweep; it is faster only
  on near-uniform data); at `p = 0.9` it compresses to CR 7.2 instead of 15.3.
  Use it for data that is close to uniform (already-compressed or
  encrypted-looking input, byte planes of numeric arrays), or when the 13-bit
  header matters on tiny blocks.
- **`Precise` is not automatically the best ratio.** When `L` is small — short
  blocks, or a wide alphabet with very rare symbols — its cumulative rounding is
  worse than water-filling. At `N` = 64 KiB, `p = 0.003` the codec picks
  `L` = 512, and there `Precise` costs +0.042 bit/symbol against +0.020 for
  `Pow23`; given `L` = 4096 on the same data the two swap places (+0.0005
  against +0.005). For large blocks with a large `L` it is the most accurate
  mode and wins.
- **Decoding does not constrain the choice**: all four modes decode within a few
  percent of each other, because the decoder only sees the frequencies, not the
  mode.

Rules of thumb:

- Start with `Precise` for large blocks and `Pow235` up to ~64 KiB. `Pow23`
  only if the block is so small that 96 header bytes really matter.
- If a single byte value takes more than about half of the block, do not use
  `Pow2`: that is exactly the case where a power-of-two frequency grid fails.
- Water-filling is cheap on skewed data (~1 µs) but not on near-uniform data
  with a small `L` (~105 µs at `p = 0.01`, even on 64 KiB); `Precise` and
  `Pow2` avoid that cost.

## Notes

This implementation is optimized for speed and is designed for experimentation
and performance-oriented use cases. The code focuses on efficient rANS
execution while keeping the design compact and approachable.
