# LUNA
### LUNA: Quasi-Optimally Succinct Designated-Verifier Zero-Knowledge Arguments from Lattices

This software implements the newly proposed vector encryption scheme Half-GSW (HGSW) of the paper [LUNA: Quasi-Optimally Succinct Designated-Verifier Zero-Knowledge Arguments from Lattices](https://eprint.iacr.org/2022/1690.pdf). The construction is called "**LUNA**". The presented HGSW implementation is written in C++ and compiled using gcc 11.4.0, with the C++17 standard. The implementation runs for Ng = 2^10, 2^12, and 2^14 in the given order and outputs times for the components Setup, Encrypt, Add, and Decrypt for each Ng value. At the end of each run, it also does a validation check to ensure every operation is done correctly. Note that it is better to run the implementation with a PC having at least 32GB RAM since the step Ng = 2^14 needs up to ~20GB of memory to store all the computed variables, and the use of swap memory can create inaccuracies in the runtimes of the components. The prerequisites for running the implementation are as follows:

  - A Linux-based OS
  - A CPU with AVX support (for the pseudorandom generator with AES-256 CTR mode to work)
  - Cmake installation (minimum version 3.5.1)
  - gcc/g++ installation (minimum version 11.4.0; a lower version supporting intrinsic `__uint128_t` will probably work, this is a recommendation)
  - PALISADE Homomorphic Addition Software Library installation (refer to *PALISADE_installation.md*)
  - OpenMP installation (for PALISADE to work)

Please download the files to your home folder. The implementation can be built and run from a Linux Terminal with the given chain of commands:

`cd ~/HGSW/build && cmake .. && cd ~/HGSW/build && make && cd ~/HGSW/build/hgsw/bin/ && ./HGSW_test`
