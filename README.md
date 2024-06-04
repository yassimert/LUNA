# LUNA
### LUNA: Quasi-Optimally Succinct Designated-Verifier Zero-Knowledge Arguments from Lattices

This software implements the newly proposed vector encryption scheme Half-GSW (HGSW) of the paper [LUNA: Quasi-Optimally Succinct Designated-Verifier Zero-Knowledge Arguments from Lattices](https://eprint.iacr.org/2022/1690.pdf). The construction is called "**LUNA**". The presented HGSW implementation is written in C++ and compiled using gcc 11.4.0, with the C++17 standard. The implementation runs for Ng ~= 2<sup>10</sup>, 2<sup>12</sup>, and 2<sup>14</sup> in the given order and outputs times for the components Setup, Encrypt, Add, and Decrypt for each Ng value. At the end of each run, it also does a validation check to ensure every operation is done correctly. Note that it is better to run the implementation with a PC with at least 32GB of RAM if Ng ~= 2<sup>16</sup> is used, as this requires up to ~20GB of memory to store all the computed variables, and the use of swap memory can create inaccuracies in the runtimes of the components. 

#### Prerequisites

  - A Linux-based OS 
  - A CPU with AVX support (for the pseudorandom generator with AES-256 CTR mode to work)
  - Cmake installation (minimum version 3.5.1)
  - gcc/g++ installation (minimum version 11.4.0; a lower version supporting intrinsic `__uint128_t` will probably work, this is a recommendation)
  - PALISADE Homomorphic Addition Software Library installation (refer to *PALISADE_installation.md*)

#### Implementation File Structure

```bash
├── hgsw                        # folder to keep files related to HGSW
│   ├── hgsw.hpp                # HGSW Encryption Scheme
│   ├── hgsw_algos.hpp          # high-level algorithms used in HGSW (Setup, Encrypt, Add, Decrypt)
│   ├── hgsw_funcs.hpp          # low-level functions called by the high-level algorithms
│   ├── params.hpp              # parameters of HGSW
│   ├── test_hgsw.hpp           # driver program
│   ├── CMakeLists.txt          # cmake settings file to include subdirectories/files to build
│   └── rng                     # folder to keep files related to uniform random integer sampler
│       ├── aes256ctr.h         # SIMD AES256CTR implementation, used in the sampler
│       ├── rng.h               # file to generate random bytes, used in the sampler
│       └── sampler.h           # uniform random integer sampler mod p and mod qbar
├── build                       # folder to keep build files
├── CMakeLists.txt              # cmake settings file to include necessary build flags and directories
├── LICENSE                     # MIT License
├── PALISADE_Installation.md    # instructions to install PALISADE library
├── luna-install.sh             # script to install LUNA
└── README.md
```

#### Instructions to Run the Implementation

If you want to install LUNA using the script *luna-install.sh*, first you need to give executable rights to it with:

`sudo chmod +x ./luna-install.sh`

Then, simply run LUNA with:

`./luna-install.sh`

in your home folder. This file needs a sudo password to run.

If you want to manually install LUNA, after installing the dependencies, please download the files to your home folder. The implementation can be built and run from a Linux Terminal with the given chain of commands:

`cd ~/LUNA/HGSW/build && cmake .. && cd ~/LUNA/HGSW/build && make && cd ~/LUNA/HGSW/build/hgsw/bin/ && ./HGSW_test`
