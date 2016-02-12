#include "fft.h"
#include "sinwave.h"
#include <xtensa/tie/xt_booleans.h>
#include <xtensa/tie/fft.h>
#include "tie_defines.h"

int fix_fft(fixed fr[], fixed fi[], int m, int inverse)
{
    int mr,nn,i,j,l,k,istep, n, scale, shift;
    
    fixed qr,qi;		//even input
    fixed tr,ti;		//odd input
    fixed wr,wi;		//twiddle factor
    
    //number of input data (n = 2^m), m is the number of stages
    n = 1 << m;
    
    int size = m;

    if(n > N_WAVE) return -1;

    mr = 0;
    nn = n - 1;
    scale = 0;

    /* decimation in time - re-order data */
    for(m=1; m<=nn; ++m) 
    {
    	
#if FFT_TIE_FAST_BIT_REVERSAL 
    	mr = FFT_REVERSE_BITS(m, size);
#else
    	l = n;
    	do
    	{
    	    l >>= 1;
    	} while(mr+l > nn);
    	    	        
    	mr = (mr & (l-1)) + l;
#endif

        if(mr <= m) continue;
        
        // swap contents of memory (real part)
        tr = fr[m];
        fr[m] = fr[mr];
        fr[mr] = tr;
        
        // swap contents of memory (imaginary part)
        ti = fi[m];
        fi[m] = fi[mr];
        fi[mr] = ti;
    }
	    
    l = 1;
    k = LOG2_N_WAVE-1;
    while(l < n) // 1 run per stage
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
            shift = 1; // evaluate putting this out of loop
        }
        
        /* it may not be obvious, but the shift will be performed
           on each data point exactly once, during this pass. */
        istep = l << 1;		//step width of current butterfly
        
        // for each twiddle factor all butterfly nodes are computed (in inner for loop)
        for(m=0; m<l; ++m)
        {
            j = m << k; // j = m * (2^k)
            /* 0 <= j < N_WAVE/2 */
            
            // Calculate twiddle factor for this stage and stepwidth
            wr =  Sinewave[j+N_WAVE/4];
            wi = -Sinewave[j];
            
            if(inverse) wi = -wi;
            if(shift)
            {
            	// simply, scaling of twiddle factors to ensure maximum arithmetic precision
                wr >>= 1;
                wi >>= 1;
            }
            
#if FFT_TIE_BUTTERFLY_CALC
            int twiddle = ((wr << 16) | (wi & 0xffff));
#endif
            
            // all butterfly compute node executions with one specific twiddle factor
            for(i=m; i<n; i = i+istep)
            {
                j = i + l;
                
                //// Implementation of FFT compute node (see task Fig.2)
                
                
				// load even value (complex)
				qr = fr[i];
				qi = fi[i];
                
#if FFT_TIE_BUTTERFLY_CALC
				
				#define aligned_by_16 __attribute__ ((aligned(16)))
				fixed out_data[4] aligned_by_16;
				VR *p_out;
				
				int q = ((qr << 16) | (qi & 0xffff));
				int f = ((fr[j] << 16) | (fi[j]) & 0xffff);
                
                VR output = FFT_CALC_BUTTERFLY(q, f, twiddle, shift);
                
                p_out = (VR *)out_data;
                *p_out = output;
                
                fr[j] = out_data[3];
                fi[j] = out_data[2];
                fr[i] = out_data[1];
                fi[i] = out_data[0];
#else

                fixed d = fix_mpy(wr, fr[j]);
                fixed e = fix_mpy(wi, fi[j]);
                fixed g = fix_mpy(wi, fr[j]);
                fixed f = fix_mpy(wr, fi[j]);
                
                tr = d - e; // complex mult (real)
                ti = f + g; // complex mult (imag)
                
				if (shift) {
					// simply, scaling of even samples to match result of multiplication with scaled twiddle factor (scaling of twiddle factor before the loop)
					qr >>= 1;
					qi >>= 1;
				}

				// Summation in upper and lower FFT compute nodes (see task Fig.2)
				
				fr[j] = qr - tr;
				fi[j] = qi - ti;
				fr[i] = qr + tr;
				fi[i] = qi + ti;
#endif
            }
        }
        --k;
        l = istep;
    }

    return scale;
}
