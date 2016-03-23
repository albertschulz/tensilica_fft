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

int fix_fft_dit(fixed fr[], fixed fi[], int m, int inverse)
{
    int mr,nn,i,j,l,k, n, scale, shift;
    
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
        
        // Handling for first 3 Stages
        if (l == 1)
        {
        	register vect8x16 real, imag;
        	
	        for (i=0; i<n; i = i+8)
	        {
	        	//
	        	// 1. Stage
	        	// 
	        	k = LOG2_N_WAVE-1;
	                   
	        	
        		real = FFT_SIMD_LOAD(fr, i);
        		imag = FFT_SIMD_LOAD(fi, i);
				
        		// Butterfly (+ shuffle integrated) Berechnung aus States mit gleichen Twiddle Faktoren		
        		FFT_CALC_TWIDDLE_FACTORx4_TO_STATES(k, inverse, shift);
        		FFT_CALC_4_BUTTERFLIES_FROM_STATES(shift, real, imag);

				//
				// 2. Stage
				//
				
	        	// 2 Twiddle Faktoren berechnen
				--k;
				
				// Butterfly (+ shuffle integrated) Berechnung aus States mit unterschiedlichen Twiddle Faktoren
				FFT_CALC_TWIDDLE_FACTORx4_TO_STATES(k, inverse, shift);
				FFT_CALC_4_BUTTERFLIES_FROM_STATES_2(shift, real, imag);
                
				//
		        // 3. Stage
				// 
		        
				--k;
	        	
				// Butterfly Berechnung aus States mit 4 unterschiedlichen Twiddle Faktoren
				FFT_CALC_TWIDDLE_FACTORx4_TO_STATES(k, inverse, shift);
				FFT_CALC_4_BUTTERFLIES_FROM_STATES_4(shift, real, imag);
                
                // Werte speichern und shuffeln
				FFT_SIMD_STORE_SHUFFLED(fr, i, real, REVERSE_SHUFFLE);
				FFT_SIMD_STORE_SHUFFLED(fi, i, imag, REVERSE_SHUFFLE);
	        }
	        
	        // F�r nachfolgende Berechnungen Schrittweite auf 8 erh�hen
	        l = 4;
        }
        else { // Stages greater than 3
        	
        	register vect8x16 real_1, imag_1, real_2, imag_2;
        	
        	WUR_REG_K(k);
        	
        	for (i=0; i<n; i = i+2*l)
        	{
	        	for (m = i; m<l+i; m+=8)
	        	{
	        		// Load Values
					FFT_SIMD_LOAD_INTERLEAVED(fr, m, real_1, real_2, UPPER);
					FFT_SIMD_LOAD_INTERLEAVED(fi, m, imag_1, imag_2, UPPER);
					
					FFT_SIMD_LOAD_INTERLEAVED(fr, m+l, real_1, real_2, LOWER);
					FFT_SIMD_LOAD_INTERLEAVED(fi, m+l, imag_1, imag_2, LOWER);
	
					// Calculate twiddle factors
					vect8x16 twiddles1 = FFT_SIMD_CALC_TWIDDLE_FACTORS(m, inverse, shift);
					vect8x16 twiddles2 = FFT_SIMD_CALC_TWIDDLE_FACTORS(m+4, inverse, shift);
					
					// Do the actual calculation
					FFT_CALC_4_BUTTERFLIES_AND_SHUFFLE_DIT(real_1, imag_1, twiddles1, shift);
					FFT_CALC_4_BUTTERFLIES_AND_SHUFFLE_DIT(real_2, imag_2, twiddles2, shift);
					
					// Store Values
					FFT_SIMD_STORE_INTERLEAVED(fr, m, real_1, real_2, UPPER);
					FFT_SIMD_STORE_INTERLEAVED(fi, m, imag_1, imag_2, UPPER);
					
					FFT_SIMD_STORE_INTERLEAVED(fr, m+l, real_1, real_2, LOWER);
					FFT_SIMD_STORE_INTERLEAVED(fi, m+l, imag_1, imag_2, LOWER);
	        	}
        	}
        }
        --k;
        l <<= 1;
    }

    return scale;
}

int fix_fft_dif(fixed fr[], fixed fi[], int m, int inverse)
{
    int mr,nn,i,j,l,k, n, scale, shift;
    
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
	    
    l = n>>1;
    k = LOG2_N_WAVE-m;
    
    while(l > 0) // 1 run per stage
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
        
        // Handling for first 3 Stages
        if (l == 4)
        {
	        for (i=0; i<n; i = i+8)
	        {
	        	register vect8x16 real, imag;
	        	
	        	k = 7;
	        	//
	        	// Thirst Last Stage
	        	// 
	        	
	        	real = FFT_SIMD_LOAD_SHUFFLED(fr, i, SHUFFLE);
	        	imag = FFT_SIMD_LOAD_SHUFFLED(fi, i, SHUFFLE);
				

				// Butterfly Berechnung aus States mit 4 unterschiedlichen Twiddle Faktoren
				FFT_CALC_TWIDDLE_FACTORx4_TO_STATES(k, inverse, shift);
				FFT_CALC_4_BUTTERFLIES_FROM_STATES_4_SHUFFLE(shift, real, imag);

				//
				// Second Last Stage
				//
				
	        	// 2 Twiddle Faktoren berechnen
				++k;
				
				// Butterfly (+ shuffle integrated) Berechnung aus States mit unterschiedlichen Twiddle Faktoren
				FFT_CALC_TWIDDLE_FACTORx4_TO_STATES(k, inverse, shift);
				FFT_CALC_4_BUTTERFLIES_FROM_STATES_2_SHUFFLE(shift, real, imag);
                
				//
		        // Last Stage
				// 
				
				++k;
        		// Butterfly (+ shuffle integrated) Berechnung aus States mit gleichen Twiddle Faktoren		
        		FFT_CALC_TWIDDLE_FACTORx4_TO_STATES(k, inverse, shift);
        		FFT_CALC_4_BUTTERFLIES_FROM_STATES_NOSHUFFLE(shift, real, imag);

                // Werte speichern und shuffeln
        		FFT_SIMD_STORE(fr, i, real);
        		FFT_SIMD_STORE(fi, i, imag);
	        }
	        
	        // Return from the while loop, since all calculations are done after this stage
	        break;
        }
        else { // Really first stages ...
        	
        	WUR_REG_K(k);
        	        	
        	for (i=0; i<n; i = i+(2*l))
        	{
	        	for (m = i; m<l+i; m+=8)
	        	{
	        		vect8x16 real_1;
					vect8x16 imag_1;
					vect8x16 real_2;
					vect8x16 imag_2;
	        							
	        		// Load Values
					FFT_SIMD_LOAD_INTERLEAVED(fr, m, real_1, real_2, UPPER);
					FFT_SIMD_LOAD_INTERLEAVED(fi, m, imag_1, imag_2, UPPER);
					
					FFT_SIMD_LOAD_INTERLEAVED(fr, m+l, real_1, real_2, LOWER);
					FFT_SIMD_LOAD_INTERLEAVED(fi, m+l, imag_1, imag_2, LOWER);
	
					// Calculate twiddle factors
					vect8x16 twiddles1 = FFT_SIMD_CALC_TWIDDLE_FACTORS(m, inverse, shift);
					vect8x16 twiddles2 = FFT_SIMD_CALC_TWIDDLE_FACTORS(m+4, inverse, shift);
					
					// Do the actual calculation
					FFT_CALC_4_BUTTERFLIES_AND_SHUFFLE_DIF(real_1, imag_1, twiddles1, shift);
					FFT_CALC_4_BUTTERFLIES_AND_SHUFFLE_DIF(real_2, imag_2, twiddles2, shift);
					
					// Store Values
					FFT_SIMD_STORE_INTERLEAVED(fr, m, real_1, real_2, UPPER);
					FFT_SIMD_STORE_INTERLEAVED(fi, m, imag_1, imag_2, UPPER);
					
					FFT_SIMD_STORE_INTERLEAVED(fr, m+l, real_1, real_2, LOWER);
					FFT_SIMD_STORE_INTERLEAVED(fi, m+l, imag_1, imag_2, LOWER);
	        	}
        	}
        }
        ++k;
        l >>= 1;
    }
    
    
    
       /* decimation in frequency - re-order data */
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
    

    return scale;
}
