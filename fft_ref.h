#ifndef FFT_REF_H
#define FFT_REF_H

#include "fixed.h"

int fix_fft_dit_ref(fixed *fr, fixed *fi, int m, int inverse);
int fix_fft_dif_ref(fixed fr[], fixed fi[], int m, int inverse);

#endif //FFT_REF_H
