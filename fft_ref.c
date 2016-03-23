#include "fft_ref.h"
#include "sinwave.h"

int fix_fft_dif_ref(fixed fr[], fixed fi[], int m, int inverse)
{
    int mr,nn,i,j,l,k,istep, n, scale, shift;

    fixed qr,qi;		//even input
    fixed tr,ti;		//odd input
    fixed wr,wi;		//twiddle factor

    //number of input data
    n = 1<<m;

    if(n > N_WAVE) return -1;

    mr = 0;
    nn = n - 1;
    scale = 0;

    l = n>>1;
    k = LOG2_N_WAVE-m;
    while(l > 0)
    {
        if(inverse)
        {
            /* variable scaling, depending upon data */
            shift = 0;
            for(i=0; i<n; ++i)
            {
                j = fr[i];
                if(j < 0) j = -j;

                m = fi[i];
                if(m < 0) m = -m;

                if(j > 16383 || m > 16383)
                {
                    shift = 1;
                    break;
                }
            }
            if(shift) ++scale;
        }
        else
        {
            /* fixed scaling, for proper normalization -
               there will be log2(n) passes, so this
               results in an overall factor of 1/n,
               distributed to maximize arithmetic accuracy. */
            shift = 1;
        }

        /* it may not be obvious, but the shift will be performed
           on each data point exactly once, during this pass. */
        istep = l << 1;		//step width of current butterfly
        for(m=0; m<l; ++m)
        {
            j = m << k;
            /* 0 <= j < N_WAVE/2 */
            wr =  Sinewave[j+N_WAVE/4];
            wi = -Sinewave[j];

            if(inverse) wi = -wi;
            if(shift)
            {
                wr >>= 1;
                wi >>= 1;
            }

            for(i=m; i<n; i+=istep)
            {
                j = i + l;
                tr = fr[i] - fr[j];
                ti = fi[j] - fi[i];

                qr = fr[i] + fr[j];
                qi = fi[i] + fi[j];

                if(shift)
                {
                        qr >>= 1;
                        qi >>= 1;
                }

                fr[i] = qr;
                fi[i] = qi;
                fr[j] = fix_mpy(wr,tr) + fix_mpy(wi,ti);
                fi[j] = fix_mpy(wi,tr) - fix_mpy(wr,ti);
            }
        }

        ++k;
        l >>= 1;
    }
    
    /* decimation in frequency - re-order data */
    for(m=1; m<=nn; ++m) {
        l = n;
        do{
        	l >>= 1;
        }while(mr+l > nn);
        mr = (mr & (l-1)) + l;

        if(mr <= m) continue;
        tr = fr[m];
        fr[m] = fr[mr];
        fr[mr] = tr;

        ti = fi[m];
        fi[m] = fi[mr];
        fi[mr] = ti;
    }

    return scale;
}

int fix_fft_dit_ref(fixed fr[], fixed fi[], int m, int inverse)
{
    int mr,nn,i,j,l,k,istep, n, scale, shift;
    
    fixed qr,qi;		//even input
    fixed tr,ti;		//odd input
    fixed wr,wi;		//twiddle factor
    
    //number of input data
    n = 1<<m;

    if(n > N_WAVE) return -1;

    mr = 0;
    nn = n - 1;
    scale = 0;

    /* decimation in time - re-order data */
    for(m=1; m<=nn; ++m) {
        l = n;
        do{
        	l >>= 1;
        }while(mr+l > nn);
        mr = (mr & (l-1)) + l;

        if(mr <= m) continue;
        tr = fr[m];
        fr[m] = fr[mr];
        fr[mr] = tr;
        
        ti = fi[m];
        fi[m] = fi[mr];
        fi[mr] = ti;
    }
	
    
    l = 1;
    k = LOG2_N_WAVE-1;
    while(l < n)
    {
        if(inverse)
        {
            /* variable scaling, depending upon data */
            shift = 0;
            for(i=0; i<n; ++i)
            {
                j = fr[i];
                if(j < 0) j = -j;
                
                m = fi[i];
                if(m < 0) m = -m;
                
                if(j > 16383 || m > 16383)
                {
                    shift = 1;
                    break;
                }
            }
            if(shift) ++scale;
        }
        else
        {
            /* fixed scaling, for proper normalization -
               there will be log2(n) passes, so this
               results in an overall factor of 1/n,
               distributed to maximize arithmetic accuracy. */
            shift = 1;
        }
        
        /* it may not be obvious, but the shift will be performed
           on each data point exactly once, during this pass. */
        istep = l << 1;		//step width of current butterfly
        for(m=0; m<l; ++m)
        {
            j = m << k;
            /* 0 <= j < N_WAVE/2 */
            wr =  Sinewave[j+N_WAVE/4];
            wi = -Sinewave[j];
            
            if(inverse) wi = -wi;
            if(shift)
            {
                wr >>= 1;
                wi >>= 1;
            }
                        
            for(i=m; i<n; i+=istep)
            {
            	
                j = i + l;
                tr = fix_mpy(wr,fr[j]) - fix_mpy(wi,fi[j]);
                ti = fix_mpy(wr,fi[j]) + fix_mpy(wi,fr[j]);
                
                qr = fr[i];
                qi = fi[i];
                
                if(shift)
                {
                        qr >>= 1;
                        qi >>= 1;
                }
                
                fr[j] = qr - tr;
                fi[j] = qi - ti;
                fr[i] = qr + tr;
                fi[i] = qi + ti;
            }
        }
        --k;
        l = istep;
    }

    return scale;
}
