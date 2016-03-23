#include        <stdio.h>
#include        <math.h>
#include		"fixed.h"
#include 		"fft.h"
#include 		"fft_ref.h"

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
void printMultipleArrays(char *title1, char *title2, fixed array1[], fixed array2[], fixed array3[], fixed array4[], int n)
{
#if M <= 4
	printf("%-16s %s\n", title1, title2);
    int i = 0;
    for (i=0; i<n; i++)
    {
    	printf("%2d: %5d, %5d | %5d, %5d\n", i, array1[i], array2[i], array3[i], array4[i]);
    }
#else
    printf("Notice: Will not print data since too much sampling points\n");
#endif
}

void printInputData(fixed array1[], fixed array2[], int n)
{
#if M <= 4
	printf("\nInput Data\n");
    int i = 0;
    for (i=0; i<n; i++)
    {
    	printf("%2d: %5d, %5d\n", i, array1[i], array2[i]);
    }
#else
    printf("Notice: Will not print data since too much sampling points\n");
#endif
}


int main()
{
	generateInputData(real, imag);
	generateInputData(real_ref, imag_ref);
	
    printInputData(real, imag, NumberOfPoints);
	int scale = 0;
	int scale_ref = 0;

    //FFT
    scale 		= fix_fft_dit(real, imag, M, 0);
    scale_ref 	= fix_fft_dit_ref(real_ref, imag_ref, M, 0);
    printMultipleArrays("FFT_DIT", "FFT_DIT (reference)", real, imag, real_ref, imag_ref, NumberOfPoints);
    printf("Scales FFT (Ref, Custom): (%d, %d)\n", scale_ref, scale);
    
    //IFFT
    scale 		= fix_fft_dit(real, imag, M, 1);
    scale_ref 	= fix_fft_dit_ref(real_ref, imag_ref, M, 1);
    printMultipleArrays("IFFT_DIT", "IFFT_DIT (reference)", real, imag, real_ref, imag_ref, NumberOfPoints);
    printf("Scales IFFT (Ref, Custom): (%d, %d)\n", scale_ref, scale);
    
    fixed error_real, error_imag;
	error_real = squaredError(real_ref, real, NumberOfPoints);
	error_imag = squaredError(imag_ref, imag, NumberOfPoints);
	printf("Squared Error (Real, Imag): (%d, %d)\n", error_real, error_imag);
	
	printf("-----------------------------------\n");
	
	generateInputData(real, imag);
	generateInputData(real_ref, imag_ref);
	
    printInputData(real, imag, NumberOfPoints);
	
    //FFT
    scale = fix_fft_dif(real, imag, M, 0);
    scale_ref = fix_fft_dif_ref(real_ref, imag_ref, M, 0);
    printMultipleArrays("FFT_DIF", "FFT_DIF (reference)", real, imag, real_ref, imag_ref, NumberOfPoints);
    printf("Scales FFT (Ref, Custom): (%d, %d)\n", scale_ref, scale);
    
    //IFFT
    scale = fix_fft_dif(real, imag, M, 1);
    scale_ref = fix_fft_dif_ref(real_ref, imag_ref, M, 1);
    printMultipleArrays("IFFT_DIF", "IFFT_DIF (reference)", real, imag, real_ref, imag_ref, NumberOfPoints);
    printf("Scales IFFT (Ref, Custom): (%d, %d)\n", scale_ref, scale);
    
	error_real = squaredError(real_ref, real, NumberOfPoints);
	error_imag = squaredError(imag_ref, imag, NumberOfPoints);
	printf("Squared Error (Real, Imag): (%d, %d)\n", error_real, error_imag);
	
	return 0;
}
