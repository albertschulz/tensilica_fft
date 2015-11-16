#include "fft.h"
#include "sinwave.h"

int fix_fft(fixed fr[], fixed fi[], int numberOfStages, int inverse)
{
    int mr,nn,i,j,l,k,istep, numberOfInputData, scale, shift;
    
    fixed qr,qi;		//even input
    fixed tr,ti;		//odd input
    fixed wr,wi;		//twiddle factor
    
    //number of input data (n = 2^m), m is the number of stages
    numberOfInputData = 1 << numberOfStages;

    if(numberOfInputData > N_WAVE) return -1;

    mr = 0;
    nn = numberOfInputData - 1;
    scale = 0;

    /* decimation in time - re-order data */
    for(numberOfStages=1; numberOfStages<=nn; ++numberOfStages) 
    {
        l = numberOfInputData;
        do
        {
        	l >>= 1;
        } while(mr+l > nn);
        
        mr = (mr & (l-1)) + l;

        if(mr <= numberOfStages) continue;
        
        tr = fr[numberOfStages];
        fr[numberOfStages] = fr[mr];
        fr[mr] = tr;
        
        ti = fi[numberOfStages];
        fi[numberOfStages] = fi[mr];
        fi[mr] = ti;
    }
	    
    l = 1;
    k = LOG2_N_WAVE-1;
    while(l < numberOfInputData)
    {
        if(inverse)
        {
            /* variable scaling, depending upon data */
            shift = 0;
            for(i=0; i<numberOfInputData; ++i)
            {
                j = fr[i];
                if(j < 0) j = -j;
                
                numberOfStages = fi[i];
                if(numberOfStages < 0) numberOfStages = -numberOfStages;
                
                if(j > 16383 || numberOfStages > 16383)
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
        for(numberOfStages=0; numberOfStages<l; ++numberOfStages)
        {
            j = numberOfStages << k;
            /* 0 <= j < N_WAVE/2 */
            wr =  Sinewave[j+N_WAVE/4];
            wi = -Sinewave[j];
            
            if(inverse) wi = -wi;
            if(shift)
            {
                wr >>= 1;
                wi >>= 1;
            }
            for(i=numberOfStages; i<numberOfInputData; i+=istep)
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
