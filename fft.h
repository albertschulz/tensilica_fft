#ifndef FFT_H
#define FFT_H

#include "fixed.h"

/*
 *	fix_fft_dit() - perform fast Fourier transform using dit algorithm
 *
 *  if n>0 FFT is done, if n<0 inverse FFT is done
 *	fr[n],fi[n] are real,imaginary arrays, INPUT AND RESULT.
 *	size of data = 2^m
 *  set inverse to 0=dft, 1=idft
 */
int fix_fft_dit(fixed *fr, fixed *fi, int m, int inverse);

/*
 *	fix_fft_dif() - perform fast Fourier transform using dif algorithm
 *
 *  if n>0 FFT is done, if n<0 inverse FFT is done
 *	fr[n],fi[n] are real,imaginary arrays, INPUT AND RESULT.
 *	size of data = 2^m
 *  set inverse to 0=dft, 1=idft
 */
int fix_fft_dif(fixed fr[], fixed fi[], int m, int inverse);

#endif //FFT_H


