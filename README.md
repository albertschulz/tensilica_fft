# Accelarated FFT Implementation using Tensilica TIE Extension

This project contains:
- reference C implementation of DIT and DIF algorithm (fft and inverse fft)
- accelerated implementation using TIE extensions (fused operations, SIMD, FLIX)

## Performance Results for DIT algorithm
The following results apply for compilation with optimization level -O2 and feedback optimization.

N=8:      243
N=256:    9533
N=1024:   41573

Note:
This task was done for the HW/SW Codesign Lab at TU Dresden.
