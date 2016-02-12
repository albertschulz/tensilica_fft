#ifndef TIE_DEFINES_H
#define TIE_DEFINES_H


/*
 * Enable these Settings to speed up FFT calculation
 */

// Speeds up bit reversal from 454 to 112 cycles
#define FFT_TIE_FAST_BIT_REVERSAL	1

// Speeds up calculation of one butterfly node
#define FFT_TIE_BUTTERFLY_CALC 		1

// Speeds up twiddle factor calculation
#define FFT_TIE_CALC_TWIDDLE_FACTORS 1

#endif // TIE_DEFINES_H
