#include "fft.h"
#include "sinwave.h"
#include <xtensa/tie/xt_booleans.h>
#include <xtensa/tie/fft.h>
#define COMPILER_BARRIER() asm volatile("" ::: "memory")

#define aligned_by_16 __attribute__ ((aligned(16)))
#define aligned_by_8 __attribute__ ((aligned(8)))
#define aligned_by_4 __attribute__ ((aligned(4)))
#define fixed_complex int

#define SHUFFLE 0
#define REVERSE_SHUFFLE 1

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
 
	        // for each twiddle factor all butterfly nodes are computed (in inner for loop)
	        for(m=0; m<l; m = m+1)
	        {
	            j = m << k; // j = m * (2^k)
	            /* 0 <= j < N_WAVE/2 */
	            
	            // Calculate twiddle factor for this stage and stepwidth
	            fixed_complex twiddle = FFT_CALC_TWIDDLE_FACTOR(j, inverse, shift); // TODO: use m,k instead of j
	
	            // all butterfly compute node executions with one specific twiddle factor
	            for(i=m; i<n; i = CALC_I(i, l))
	            {	                
	                //// Implementation of FFT compute node (see task Fig.2)
	                
	                // Load fr_addr, fi_addr into special register
	                SET_FR_FI_ADDR(fr, fi, i);
	                
	                LOAD_INTO_REAL_REG(0, 0);
	                LOAD_INTO_IMAG_REG(0, 0);
	                COMPILER_BARRIER();
	                LOAD_INTO_REAL_REG(1, l);
	                LOAD_INTO_IMAG_REG(1, l);
	                COMPILER_BARRIER();
	                LOAD_INTO_REAL_REG(2, l);
	                LOAD_INTO_IMAG_REG(2, l);
	                COMPILER_BARRIER();
	                LOAD_INTO_REAL_REG(3, l);
	                LOAD_INTO_IMAG_REG(3, l);
	                COMPILER_BARRIER();
					LOAD_INTO_REAL_REG(4, l);
					LOAD_INTO_IMAG_REG(4, l);
					COMPILER_BARRIER();
					LOAD_INTO_REAL_REG(5, l);
					LOAD_INTO_IMAG_REG(5, l);
					COMPILER_BARRIER();
					LOAD_INTO_REAL_REG(6, l);
					LOAD_INTO_IMAG_REG(6, l);
					COMPILER_BARRIER();
					LOAD_INTO_REAL_REG(7, l);
					LOAD_INTO_IMAG_REG(7, l);
					COMPILER_BARRIER();
	                
	                FFT_CALC_4_BUTTERFLIES_FROM_STATES_4_WITH_TWIDDLE(twiddle, shift);

	                STORE_FROM_REAL_REG(0, 0);
	                STORE_FROM_IMAG_REG(0, 0);
	                COMPILER_BARRIER();
	                STORE_FROM_REAL_REG(1, l);
	                STORE_FROM_IMAG_REG(1, l);
	                COMPILER_BARRIER();
	                STORE_FROM_REAL_REG(2, l);
	                STORE_FROM_IMAG_REG(2, l);
	                COMPILER_BARRIER();
	                STORE_FROM_REAL_REG(3, l);
	                STORE_FROM_IMAG_REG(3, l);
	                COMPILER_BARRIER();
	                STORE_FROM_REAL_REG(4, l);
	                STORE_FROM_IMAG_REG(4, l);
	                COMPILER_BARRIER();
	                STORE_FROM_REAL_REG(5, l);
	                STORE_FROM_IMAG_REG(5, l);
	                COMPILER_BARRIER();
	                STORE_FROM_REAL_REG(6, l);
	                STORE_FROM_IMAG_REG(6, l);
	                COMPILER_BARRIER();
	                STORE_FROM_REAL_REG(7, l);
	                STORE_FROM_IMAG_REG(7, l);
	                COMPILER_BARRIER();
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




