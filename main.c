#include        <stdio.h>
#include        <math.h>
#include		"fixed.h"
#include 		"fft.h"

#define M       3
#define NumberOfPoints       (1<<M)

/*
 * select TIE version:
 * TODO: add different tie versions to have a speedup
 */
#define FFT_TIE	0

fixed real[NumberOfPoints], imag[NumberOfPoints];

void generateInputData() 
{
	//const float pi = 3.1415926535;
	
	int i;
    for(i=0; i < NumberOfPoints; i++)
    {
        real[i] = 1000*cos(i*2*3.1415926535/NumberOfPoints);
        imag[i] = 0;
    }
}

/*
 * Prints the content of both arrays
 */
void printData(fixed array1[], fixed array2[], int n)
{
    int i = 0;
    for (i=0; i<n; i++)
    {
    	printf("%d: %d, %d\n", i, array1[i], array2[i]);
    }
}

int main()
{
	generateInputData();
	
    printf("\nInput Data\n");
	printData(real, imag, NumberOfPoints);
	
    //FFT
    fix_fft(real, imag, M, 0);
    
    printf("\nFFT\n");
    printData(real, imag, NumberOfPoints);
    
    //IFFT
    fix_fft(real, imag, M, 1);
    
    printf("\nIFFT\n");
    printData(real, imag, NumberOfPoints);
	
	return 0;
}
