#include "fft.h"
#include "sinwave.h"
#include <xtensa/tie/xt_booleans.h>
#include <xtensa/tie/fft.h>

#define aligned_by_16 __attribute__ ((aligned(16)))
#define aligned_by_8 __attribute__ ((aligned(8)))
#define aligned_by_4 __attribute__ ((aligned(4)))
#define fixed_complex int

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
        istep = l << 1; // (2^L)		//step width of current butterfly
        
        
        
        if (istep == 2)
        {
        	int m = 0;
        	
        	j = m << k; // j = m * (2^k)
                        
        	fixed_complex twiddle = FFT_CALC_TWIDDLE_FACTOR(j, inverse, shift);
            
        	for (i=m; i<n; i = i+8)
        	{
        		j = i + l;
        		                
                //// Implementation of FFT compute node (see task Fig.2)

				LOAD_REAL(fr, i);
				LOAD_IMAG(fi, i);
				
				FFT_CALC_4_BUTTERFLIES_FROM_STATES(twiddle, shift);

                STORE_REAL(fr, i);
                STORE_IMAG(fi, i);
        	}
        }
        if (istep == 4) 
        {
        	// Werte umsortieren
        	
        	// Alles f�r 2 Werte f�r m: 0, 1
        	
        	// 2 Twiddle Faktoren berechnen
        	int j1 = 0 << k; // m=0
        	int j2 = 1 << k; // m=1
        	                        
        	fixed_complex tw1 = FFT_CALC_TWIDDLE_FACTOR(j1, inverse, shift);
        	fixed_complex tw2 = FFT_CALC_TWIDDLE_FACTOR(j2, inverse, shift);
        	        	
        	// Butterfly Berechnung mit 2 verschiedenen Twiddle Faktoren
        	for (i=0; i<n; i = i+8)
        	{
                //// Implementation of FFT compute node (see task Fig.2)

				fixed qr1 = fr[i];
				fixed qr2 = fr[i+1];
				fixed fr1 = fr[i+2];
				fixed fr2 = fr[i+3];
				fixed qr3 = fr[i+4];
				fixed qr4 = fr[i+5];
				fixed fr3 = fr[i+6];
				fixed fr4 = fr[i+7];
				
				fixed qi1 = fi[i];
				fixed qi2 = fi[i+1];
				fixed fi1 = fi[i+2];
				fixed fi2 = fi[i+3];
				fixed qi3 = fi[i+4];
				fixed qi4 = fi[i+5];
				fixed fi3 = fi[i+6];
				fixed fi4 = fi[i+7];
                				
				fixed input_q_1[4] aligned_by_8 = {qi1, qr1, qi3, qr3};
				fixed input_f_1[4] aligned_by_8 = {fi1, fr1, fi3, fr3};
				fixed input_q_2[4] aligned_by_8 = {qi2, qr2, qi4, qr4};
				fixed input_f_2[4] aligned_by_8 = {fi2, fr2, fi4, fr4};
				
				fixed out_data5[8] aligned_by_16;
				fixed out_data6[8] aligned_by_16;
				vect8x16 *p_out5 = (vect8x16 *)out_data5;
				vect8x16 *p_out6 = (vect8x16 *)out_data6;
				
				vec4x16 *input_q_1_vector = (vec4x16 *)input_q_1;
				vec4x16 *input_f_1_vector = (vec4x16 *)input_f_1;
				vec4x16 *input_q_2_vector = (vec4x16 *)input_q_2;
				vec4x16 *input_f_2_vector = (vec4x16 *)input_f_2;
				
				// Alle Daten aus Speicher laden, und shufflen
				FFT_SHUFFLE_READ_REAL(fr, i);
				FFT_SHUFFLE_READ_IMAG(fi, i);
				
				
				*p_out5 = FFT_CALC_2_BUTTERFLIES(*input_q_1_vector, *input_f_1_vector, tw1, shift);
				*p_out6 = FFT_CALC_2_BUTTERFLIES(*input_q_2_vector, *input_f_2_vector, tw2, shift);
				
				// TODO: Butterfly Berechnung aus States mit unterschiedlichen Twiddle Faktoren
				
				int twiddle_array[2] aligned_by_4 = {tw2, tw1};
				vec4x16 *twiddle_array_ptr = (vec4x16 *)twiddle_array;
				
				FFT_CALC_4_BUTTERFLIES_FROM_STATES_2(*twiddle_array_ptr, shift);
                              
				fr[i] 	= out_data5[1];
				fr[i+1] = out_data6[1];
                fr[i+2] = out_data5[3];
                fr[i+3] = out_data6[3];
                fr[i+4] = out_data5[5];
                fr[i+5] = out_data6[5];
                fr[i+6] = out_data5[7];
                fr[i+7] = out_data6[7];
                
                fi[i] 	= out_data5[0];
                fi[i+1] = out_data6[0];
                fi[i+2] = out_data5[2];
                fi[i+3] = out_data6[2];
                fi[i+4] = out_data5[4];
                fi[i+5] = out_data6[4];
                fi[i+6] = out_data5[6];
                fi[i+7] = out_data6[6];
                
                // Werte speichern
        	}
        }
        else if (istep == 8)
        {
        	// Werte umsortieren
        	
        	// 4 Twiddle Faktoren berechnen
        	
        	// Butterfly Berechnung mit 4 versch. Twiddle Faktoren
        	
        	// Werte speichern
        }
        
        // Werte in RAM ablegen (evtl. au�erhalb der Schleife)
        
        
        // for each twiddle factor all butterfly nodes are computed (in inner for loop)
        for(m=0; m<l; ++m)
        {
            j = m << k; // j = m * (2^k)
            /* 0 <= j < N_WAVE/2 */
            
            // Calculate twiddle factor for this stage and stepwidth
            fixed_complex twiddle = FFT_CALC_TWIDDLE_FACTOR(j, inverse, shift);

            if (istep == 2 || istep == 4) 
            {}
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




