#include "fft.h"
#include "sinwave.h"
#include <xtensa/tie/xt_booleans.h>
#include <xtensa/tie/fft.h>

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

        istep = l << 1; // (2^L)		//step width of current butterfly
        
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
	            for(i=m; i<n; i = CALC_I(i, istep)) // Use Fused Multiply Add for i calculation
	            {
	                j = i + l;
	                
	                //// Implementation of FFT compute node (see task Fig.2)
	                
	                // Load fr_addr into special register
	                SET_FR_I_ADDR(fr, i);
	                
					// Load Real Values
	                fixed qr1 = fr[i];
	                fixed qr2 = fr[i + istep];
	                fixed qr3 = fr[i + istep*2];
	                fixed qr4 = fr[i + istep*3];
	                fixed fr1 = fr[i +         + l];
	                fixed fr2 = fr[i + istep   + l];
	                fixed fr3 = fr[i + istep*2 + l];
	                fixed fr4 = fr[i + istep*3 + l];
	                
	                // even, odd, even, odd... (i, j, i, j, ...)
	                LOAD_INTO_REAL_REG(0);
	                LOAD_INTO_REAL_REG(l);
	                LOAD_INTO_REAL_REG(istep);
	                LOAD_INTO_REAL_REG(istep+l);
					LOAD_INTO_REAL_REG(istep*2);
					LOAD_INTO_REAL_REG(istep*2+l);
					LOAD_INTO_REAL_REG(istep*3);
					LOAD_INTO_REAL_REG(istep*3+l);
					
	                // Load Complex Values
	                fixed qi1 = fi[i];
	                fixed qi2 = fi[i + istep];
	                fixed qi3 = fi[i + istep*2];
	                fixed qi4 = fi[i + istep*3];
	                fixed fi1 = fi[i +         + l];
	                fixed fi2 = fi[i + istep   + l];
	                fixed fi3 = fi[i + istep*2 + l];
	                fixed fi4 = fi[i + istep*3 + l];
					
	                fixed out_data1[4] aligned_by_16;
	                fixed out_data2[4] aligned_by_16;
	                fixed out_data3[4] aligned_by_16;
	                fixed out_data4[4] aligned_by_16;
	                
	                vec4x16 *p_out1 = (vec4x16 *)out_data1;
	                vec4x16 *p_out2 = (vec4x16 *)out_data2;
	                vec4x16 *p_out3 = (vec4x16 *)out_data3;
	                vec4x16 *p_out4 = (vec4x16 *)out_data4;
					
	                fixed q_complex1[2] = {qi1, qr1};
	                fixed q_complex2[2] = {qi2, qr2};
	                fixed q_complex3[2] = {qi3, qr3};
	                fixed q_complex4[2] = {qi4, qr4};
	                
	                fixed f_complex1[2] = {fi1, fr1};
	                fixed f_complex2[2] = {fi2, fr2};
	                fixed f_complex3[2] = {fi3, fr3};
	                fixed f_complex4[2] = {fi4, fr4};
	                
	                // TODO: Combine to one butterfly calculation
	                *p_out1 = FFT_CALC_BUTTERFLY(*(int*)q_complex1, *(int*)f_complex1, twiddle, shift);
	                *p_out2 = FFT_CALC_BUTTERFLY(*(int*)q_complex2, *(int*)f_complex2, twiddle, shift);
	                *p_out3 = FFT_CALC_BUTTERFLY(*(int*)q_complex3, *(int*)f_complex3, twiddle, shift);
	                *p_out4 = FFT_CALC_BUTTERFLY(*(int*)q_complex4, *(int*)f_complex4, twiddle, shift);
	                
					// Real Values
	                fr[i] 				= out_data1[3];
	                fr[i+istep] 		= out_data2[3];
	                fr[i+istep*2] 		= out_data3[3];
	                fr[i+istep*3]		= out_data4[3];
	                fr[i + 		   + l] = out_data1[1];
	                fr[i + istep   + l] = out_data2[1];
	                fr[i + istep*2 + l] = out_data3[1];
	                fr[i + istep*3 + l] = out_data4[1];
	                
	                // Imag Values
	                fi[i] 				= out_data1[2];
	                fi[i + istep] 		= out_data2[2];
	                fi[i + istep*2] 	= out_data3[2];
	                fi[i + istep*3] 	= out_data4[2];
	                fi[i + 		   + l] = out_data1[0];
	                fi[i + istep   + l] = out_data2[0];
					fi[i + istep*2 + l] = out_data3[0];
					fi[i + istep*3 + l] = out_data4[0];
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




