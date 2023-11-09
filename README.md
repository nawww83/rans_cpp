# rans_cpp
Range ANS codec C++ implemetation

# Features
- Two-stream mode: even/odd input bytes are encoded by its own stream
- 16-bit flushing, 32-bit arithmetic
- Frequencies sorting
- 16-bit per frequency in the frequency table
- Auto scale L selection
- Free-remainder exact rANS transformation: only "div" is used without "mod" function
- Fast byte counter

# Build
g++ main.cpp -O3 -std=c++20 -o rans

# Run
./rans

# Average performance
- Compression: 160MB/s
- Decompression: 225 MB/s
On conditions:
- The input size is more than 2MB
- CPU Intel i7-8565U CPU 4.2GHz
- Compiler GCC 11.4
- Ubuntu 22.04 WSL