#ifndef FFT_H
#define FFT_H

#include "fixed.h"

/*
 *	fix_fft() - perform fast Fourier transform.
 *
 *  if n>0 FFT is done, if n<0 inverse FFT is done
 *	fr[n],fi[n] are real,imaginary arrays, INPUT AND RESULT.
 *	size of data = 2^m
 *  set inverse to 0=dft, 1=idft
 */
int fix_fft(fixed *fr, fixed *fi, int m, int inverse);

#endif //FFT_H


