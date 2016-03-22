#include "fft.h"
#include "sinwave.h"
#include <xtensa/tie/xt_booleans.h>
#include <xtensa/tie/fft.h>
#include <stdio.h>

#define aligned_by_16 __attribute__ ((aligned(16)))
#define aligned_by_8 __attribute__ ((aligned(8)))
#define aligned_by_4 __attribute__ ((aligned(4)))
#define fixed_complex int

#define SHUFFLE 0
#define REVERSE_SHUFFLE 1

#define UPPER 1
#define LOWER 0

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
        
        // swap contents of memory (real part & imaginary part)
        tr = fr[m];
        ti = fi[m];
        
        fr[m] = fr[mr];
        fi[m] = fi[mr];
        
        fr[mr] = tr;
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
            for(i=0; i<n; i=i+8)
            {
				vect8x16 j_vector = FFT_SIMD_LOAD(fr, i);
				vect8x16 m_vector = FFT_SIMD_LOAD(fi, i);
            	
            	if (FFT_CHECK_SHIFT_CONDITION(j_vector, m_vector))
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

        istep = l << 1; // (2*L)		//step width of current butterfly
        
        // Handling for first 3 Stages
        if (istep == 2)
        {
	        for (i=0; i<n; i = i+8)
	        {
	        	//
	        	// 1. Stage
	        	// 
	        	k = LOG2_N_WAVE-1;
	                   
	        	
        		FFT_SIMD_LOAD_REAL(fr, i);
        		FFT_SIMD_LOAD_IMAG(fi, i);
				
        		// Butterfly (+ shuffle integrated) Berechnung aus States mit gleichen Twiddle Faktoren		
        		FFT_CALC_TWIDDLE_FACTORx4_TO_STATES(k, inverse, shift);
        		FFT_CALC_4_BUTTERFLIES_FROM_STATES(shift);

				//
				// 2. Stage
				//
				
	        	// 2 Twiddle Faktoren berechnen
				--k;
				
				// Butterfly (+ shuffle integrated) Berechnung aus States mit unterschiedlichen Twiddle Faktoren
				FFT_CALC_TWIDDLE_FACTORx4_TO_STATES(k, inverse, shift);
				FFT_CALC_4_BUTTERFLIES_FROM_STATES_2(shift);
                
				//
		        // 3. Stage
				// 
		        
				--k;
	        	
				// Butterfly Berechnung aus States mit 4 unterschiedlichen Twiddle Faktoren
				FFT_CALC_TWIDDLE_FACTORx4_TO_STATES(k, inverse, shift);
				FFT_CALC_4_BUTTERFLIES_FROM_STATES_4(shift);
                
                // Werte speichern und shuffeln
				FFT_SIMD_SHUFFLE_STORE_REAL(fr, i, REVERSE_SHUFFLE);
				FFT_SIMD_SHUFFLE_STORE_IMAG(fi, i, REVERSE_SHUFFLE);
	        }
	        
	        // Für nachfolgende Berechnungen Schrittweite auf 8 erhöhen
	        istep = 8;
        }
        else { // Stages greater than 3
        	
        	for (i=0; i<n; i = i+istep)
        	{
	        	for (m = i; m<l+i; m+=4)
	        	{
	        		// Load Values
					vect8x16 even_odd_r;
					vect8x16 even_odd_i;
	
					// Green
					FFT_SIMD_LOAD_EXTENDED(fr, m, even_odd_r, UPPER);
					FFT_SIMD_LOAD_EXTENDED(fi, m, even_odd_i, UPPER);
					FFT_SIMD_LOAD_EXTENDED(fr, m+l, even_odd_r, LOWER);
					FFT_SIMD_LOAD_EXTENDED(fi, m+l, even_odd_i, LOWER);
	
					// Calculate twiddle factors
					
					int j1 = m<<k;
					int j2 = (m+1)<<k;
					int j3 = (m+2)<<k;
					int j4 = (m+3)<<k;
					
					fixed_complex tw1 = FFT_CALC_TWIDDLE_FACTOR(j1, inverse, shift);
					fixed_complex tw2 = FFT_CALC_TWIDDLE_FACTOR(j2, inverse, shift);
					fixed_complex tw3 = FFT_CALC_TWIDDLE_FACTOR(j3, inverse, shift);
					fixed_complex tw4 = FFT_CALC_TWIDDLE_FACTOR(j4, inverse, shift);
					
					fixed_complex twiddle_vect[] = {tw4, tw3, tw2, tw1};
					
					// Do the actual calculation
					vect8x16 evenodd_r_out = FFT_CALC_4_BUTTERFLIES_REAL(even_odd_r, even_odd_i, *(vect8x16*)twiddle_vect, shift);
					vect8x16 evenodd_i_out = FFT_CALC_4_BUTTERFLIES_IMAG(even_odd_r, even_odd_i, *(vect8x16*)twiddle_vect, shift);
					
					evenodd_r_out = FFT_SHUFFLE(evenodd_r_out);
					evenodd_i_out = FFT_SHUFFLE(evenodd_i_out);
					
					// Store Values
					FFT_SIMD_STORE_EXTENDED(fr, m, evenodd_r_out, UPPER);
					FFT_SIMD_STORE_EXTENDED(fi, m, evenodd_i_out, UPPER);
					FFT_SIMD_STORE_EXTENDED(fr, m+l, evenodd_r_out, LOWER);
					FFT_SIMD_STORE_EXTENDED(fi, m+l, evenodd_i_out, LOWER);
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




