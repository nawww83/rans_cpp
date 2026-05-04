# rans_cpp
Byte-wise range ANS codec C++ implementation

# Features
- Quad interleave mode
- 16-bit flushing, 32-bit arithmetic
- Symbol frequency sorting
- 16-bit per frequency in the frequency table
- Auto scale L selection
- Free-remainder rANS transformation: only "div" is used without "mod" function
- Fast symbol counter

# Build
g++ main.cpp -O3 -std=c++20 -o rans

# Run
./rans

# Average performance
- Compression: 350 MB/s per core
- Decompression: 510 MB/s per core

Under the conditions:
- The input size is in (2...8) MB
- CPU 12th Gen Intel© Core™ i7-12700 (~4500 MHz)
- Compiler GCC 11.4.1, Linux
