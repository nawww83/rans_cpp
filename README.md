# rans_cpp

Byte-wise range ANS codec implementation in C++.

A compact and fast entropy coder designed for high-throughput compression and decompression workloads. It emphasizes efficient integer arithmetic, low-overhead symbol handling, and a tunable speed/compression trade-off.

## Why rans_cpp?

- High throughput on modern CPUs
- Byte-oriented rANS design
- Optimized for 32-bit arithmetic and compact frequency tables
- Fast symbol counting and frequency sorting
- Selectable frequency quantization: from exact (arithmetic-coding-like) to Huffman-like
- Self-describing stream: the decoder does not need to be told the encoder mode

## Features

- Quad interleaved rANS
- 16-bit flushing and 32-bit arithmetic
- Symbol frequency sorting for cache-friendly decoding
- Compact bit-packed frequency table
- Automatic scale selection via parameter `L` (power of two, `256 <= L <= 65536`)
- Free-remainder rANS transformation using only `div` (no `mod` required)
- Fast symbol counter
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
present in the input, where `L` is a power of two.

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
| `Pow23` / `Pow235` | `x >> a`, then reciprocal multiplication by `1 / (3^b)` or `1 / (3^b * 5^c)` with an exact correction step |

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
Representative medians from the distribution sweep on a WSL laptop
(Intel Core i7 8th Gen, GCC 13, input near 4 MB):

| Mode | Redundancy | Enc MB/s | Dec MB/s |
|------|-----------:|---------:|---------:|
| `Precise` | ~0.16 % | ~160 | ~295 |
| `Pow235` | ~0.16 % | ~180 | ~295 |
| `Pow23` | ~0.20 % | ~180 | ~295 |
| `Pow2` | ~1.8 % | ~205 | ~300 |

`Pow2` is fastest but least compact. `Precise` and `Pow235` are essentially
tied on redundancy, but `Pow235` encodes faster because it replaces division
with a shift and a reciprocal multiply. `Pow23` is slightly less compact and
about as fast as `Pow235`. Decompression is essentially mode-independent,
since every mode decodes with the same multiplication-based rANS step.

## Notes

This implementation is optimized for speed and is designed for experimentation
and performance-oriented use cases. The code focuses on efficient rANS
execution while keeping the design compact and approachable.
