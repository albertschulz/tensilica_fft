#include        <stdio.h>
#include        <math.h>
#include		"fixed.h"
#include 		"fft.h"

#define M       10
#define NumberOfPoints       (1<<M) // 2^M

fixed real[NumberOfPoints] __attribute__ ((section(".dram0.data")));
fixed imag[NumberOfPoints] __attribute__ ((section(".dram1.data")));

fixed  real_ref[NumberOfPoints], imag_ref[NumberOfPoints];


// modules for compating results: squaredError() and generateInputDataRef()
int squaredError(fixed array1[], fixed array2[], int n)
{
	int i;
	fixed sum = 0;
	for (i=0; i < n; i++)
	{
		sum += (array1[i]-array2[i])*(array1[i]-array2[i]);
	}
	return sum;
}

void generateInputData(fixed real_array[], fixed imag_array[])
{	
	int i;
    for(i=0; i < NumberOfPoints; i++)
    {
    	real_array[i] = 1000*cos(i*2*3.1415926535/NumberOfPoints);
    	imag_array[i] = 0;
    }
}

/*
 * Prints the content of both arrays
 */
void printData(fixed array1[], fixed array2[], int n)
{
#if M <= 4
    int i = 0;
    for (i=0; i<n; i++)
    {
    	printf("%d: %d, %d\n", i, array1[i], array2[i]);
    }
#else
    printf("Notice: Will not print data since too much sampling points\n");
#endif
}

int main()
{
	generateInputData(real, imag);
	generateInputData(real_ref, imag_ref);
	
    printf("\nInput Data\n");
	printData(real, imag, NumberOfPoints);
	
    //FFT
	
	// Default implementation needs 2446 cycles
    fix_fft(real, imag, M, 0);
    fix_fft_ref(real_ref, imag_ref, M, 0);
    
    printf("\nFFT\n");
    printData(real, imag, NumberOfPoints);
    
    //IFFT
    
    // Default implementation needs 3470 cycles 
    // (this has been calculated by decrease the number 
    // of cycles for fft+ifft by the number of cycles 
    // when running the fft only. Only cycles for fix_fft() 
    // function calls are considered)
    fix_fft(real, imag, M, 1);
    printf("\nFFT (reference)\n");
    printData(real_ref, imag_ref, NumberOfPoints);
    
    fix_fft_ref(real_ref, imag_ref, M, 1);
    
    
    printf("\nIFFT (reference)\n");
    printData(real_ref, imag_ref, NumberOfPoints);
    
    printf("\nIFFT\n");
    printData(real, imag, NumberOfPoints);
    
    fixed error_real, error_imag;
	error_real = squaredError(real_ref, real, NumberOfPoints);
	error_imag = squaredError(imag_ref, imag, NumberOfPoints);
	printf("\nSquared Error (Real, Imag): (%d, %d)", error_real, error_imag);
	return 0;
}
