/**
 * NTTW FFT Test Program
 * \file fft_test.c
 * \brief FFT Test Program for the FRTW C Library.
 *
 * This file is part of FRTW Library.
 *
 * FRTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FRTW is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with FRTW.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2013
*/
#include <stdio.h>

#include <nttw/global.h>
#include <nttw/timing.h>

#include "array_complex.h"
#include "fourier.h"

int main(int argc, char *argv[])
{
    fftw_complex *data, *result, *inversion;
    size_t j, N = 0; //Use size_t instead of int for portability
    unsigned long long duration1 = 0, duration2 = 0;
    int rank = 1, sizes[1];

    printf(">| FFT Test\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(data));
    if(argc != 2)
    {
        printf(">| %s <n>\n",argv[0]);
        printf(">| Computes the FFT of size n, where n is dyadic (16,32,etc.).\n");
        return EXIT_FAILURE;
    }
    N = atoi(argv[1]);
    sizes[0] = (int)N;

    //-----------------------------------------
    printf(">| Allocating memory\n");
    data = fftw_array(N);
    result = fftw_array(N);
    inversion = fftw_array(N);

    printf(">| Data: \n");
    for(j = 0; j < N; j ++)
    {
        data[j][0] = 1.0 + j;
        printf(FORMAT_CSVOUTPUT_FLOAT_STR,data[j][0]);
    }
    printf("\n");

    //-----------------------------------------
    ///FFT
    printf(">| FFT TEST ----------- \n");

    START_TIMER;

    ///FNTT
    fft(rank,sizes,data,result);

    STOP_TIMER;
    duration1 = MICROSECONDS_ELAPSED;

    printf(">| FFT (Real): \n");
    for(j = 0; j < N; j ++)
        printf(FORMAT_CSVOUTPUT_FLOAT_STR,result[j][0]);
    printf("\n\n");
    printf(">| FFT (Imag): \n");
    for(j = 0; j < N; j ++)
        printf(FORMAT_CSVOUTPUT_FLOAT_STR,result[j][1]);
    printf("\n\n");

    RESTART_TIMER;

    ///iFFT
    ifft(rank,sizes,result,inversion);

    STOP_TIMER;
    duration2 = MICROSECONDS_ELAPSED;

    ///Norm
    for(j = 0; j < sizes[0]; j ++) ///Norm
        inversion[j][0] /= N;

    printf(">| iFFT: \n");
    for(j = 0; j < N; j ++)
        printf(FORMAT_CSVOUTPUT_FLOAT_STR,inversion[j][0]);
    printf("\n\n");

    printf(">| FFT Took %llu usecs.\n", duration1);
    printf(">| iFFT Took %llu usecs.\n", duration2);

    //-----------------------------------------
    printf(">| Deallocating memory\n");
    fftw_free(data);
    fftw_free(result);
    fftw_free(inversion);
    printf(">| Operation Complete\n");
    return EXIT_SUCCESS;
}
