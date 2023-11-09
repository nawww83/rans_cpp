# rans_cpp
Byte-wise range ANS codec C++ implementation

# Features
- Two-stream mode: even/odd input symbols are encoded by its own stream
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
- Compression: 160 MB/s
- Decompression: 225 MB/s

On conditions:
- The input size is more than 2MB
- CPU Intel i7-8565U CPU 4.2GHz
- Compiler GCC 11.4
- Ubuntu 22.04 WSL