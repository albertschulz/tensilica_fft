#include "fft.h"
#include "sinwave.h"
#include <xtensa/tie/xt_booleans.h>
#include <xtensa/tie/fft.h>

#define aligned_by_16 __attribute__ ((aligned(16)))

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
       	mr = FFT_REVERSE_BITS(m, size);
    	
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
            
            fixed twiddle_out[2] aligned_by_16;
            int *twiddle_out_ptr = (int *)twiddle_out;
            
            *twiddle_out_ptr = FFT_CALC_TWIDDLE_FACTOR(j, inverse, shift);


            int twiddle = *twiddle_out_ptr;
            
            if (istep == 2) {
            	
            	i = m;
            	
                j = i + l;
                
                //// Implementation of FFT compute node (see task Fig.2)
                
				// load even value (complex)
				fixed qr1 = fr[i];
				fixed qi1 = fi[i];
				fixed qr2 = fr[i+2];
				fixed qi2 = fi[i+2];
				fixed qr3 = fr[i+4];
				fixed qi3 = fi[i+4];
				fixed qr4 = fr[i+6];
				fixed qi4 = fi[i+6];
				
				fixed out_data1[4] aligned_by_16;
				fixed out_data2[4] aligned_by_16;
				fixed out_data3[4] aligned_by_16;
				fixed out_data4[4] aligned_by_16;
				
				vec4x16 *p_out1 = (vec4x16 *)out_data1;
				vec4x16 *p_out2 = (vec4x16 *)out_data2;
				vec4x16 *p_out3 = (vec4x16 *)out_data3;
				vec4x16 *p_out4 = (vec4x16 *)out_data4;
				
				int q1 = ((qr1 << 16) | (qi1 & 0xffff));
				int q2 = ((qr2 << 16) | (qi2 & 0xffff));
				int q3 = ((qr3 << 16) | (qi3 & 0xffff));
				int q4 = ((qr4 << 16) | (qi4 & 0xffff));
				
				int f1 = ((fr[j] << 16) | (fi[j]) & 0xffff);
				int f2 = ((fr[j+2] << 16) | (fi[j+2]) & 0xffff);
				int f3 = ((fr[j+4] << 16) | (fi[j+4]) & 0xffff);
				int f4 = ((fr[j+6] << 16) | (fi[j+6]) & 0xffff);
                
				*p_out1 = FFT_CALC_BUTTERFLY(q1, f1, twiddle, shift);
				*p_out2 = FFT_CALC_BUTTERFLY(q2, f2, twiddle, shift);
				*p_out3 = FFT_CALC_BUTTERFLY(q3, f3, twiddle, shift);
				*p_out4 = FFT_CALC_BUTTERFLY(q4, f4, twiddle, shift);
                                
				fr[j]   = out_data1[3];
				fr[j+2] = out_data2[3];
				fr[j+4] = out_data3[3];
				fr[j+6] = out_data4[3];
                
				fi[j]   = out_data1[2];
				fi[j+2] = out_data2[2];
				fi[j+4] = out_data3[2];
				fi[j+6] = out_data4[2];
                
                fr[i]   = out_data1[1];
                fr[i+2] = out_data2[1];
                fr[i+4] = out_data3[1];
                fr[i+6] = out_data4[1];
                
                fi[i]   = out_data1[0];
                fi[i+2] = out_data2[0];
                fi[i+4] = out_data3[0];
                fi[i+6] = out_data4[0];
            	
            }
            else {
	            // all butterfly compute node executions with one specific twiddle factor
	            for(i=m; i<n; i = i+istep)
	            {
	                j = i + l;
	                
	                //// Implementation of FFT compute node (see task Fig.2)
	                
					// load even value (complex)
					qr = fr[i];
					qi = fi[i];
					
					fixed out_data[4] aligned_by_16;
					vec4x16 *p_out = (vec4x16 *)out_data;
					
					int q = ((qr << 16) | (qi & 0xffff));
					int f = ((fr[j] << 16) | (fi[j]) & 0xffff);
	                
					*p_out = FFT_CALC_BUTTERFLY(q, f, twiddle, shift);
	                                
	                fr[j] = out_data[3];
	                fi[j] = out_data[2];
	                fr[i] = out_data[1];
	                fi[i] = out_data[0];
	            }
            }
        }
        --k;
        l = istep;
    }

    return scale;
}



//
// original FFT implementation
// and auxiliary modules
//
// for reference use only
//


/*
        fix_mpy() - fixed-point multiplication
*/
fixed fix_mpy_ref(fixed a, fixed b)
{
    FIX_MPY(a,a,b);
    return a;
}


int fix_fft_ref(fixed fr[], fixed fi[], int m, int inverse)
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
                tr = fix_mpy_ref(wr,fr[j]) - fix_mpy(wi,fi[j]);
                ti = fix_mpy_ref(wr,fi[j]) + fix_mpy(wi,fr[j]);
                
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




