/**
 * NTTW Convolution Test Program
 * \file fntt_test.c
 * \brief Convolution Test Program for the NTTW C Library.
 *
 * This file is part of NTTW Library.
 *
 * NTTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DGV is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with NTTW.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2013
*/
#include <stdio.h>

#include "nttw/array.h"
#include "nttw/timing.h"
#include "nttw/number.h"

int main(int argc, char *argv[])
{
    nttw_integer *data, *kernel, *result;
    nttw_big_integer *dataBig, *kernelBig, *resultBig;
    size_t j, N = 0; //Use size_t instead of int for portability
    unsigned long long duration1 = 0, duration2 = 0, duration3 = 0;

    printf(">| Convolution Test\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(data));
    if(argc != 2)
    {
        printf(">| %s <n>\n",argv[0]);
        printf(">| Computes the convolution of size n, where n is dyadic (16,32,etc.).\n");
        return EXIT_FAILURE;
    }
    N = atoi(argv[1]);

    //-----------------------------------------
    printf(">| Allocating memory\n");
    data = array_1D(N);
    kernel = array_1D(N);
    result = array_1D(N);

    //init
    init_1D(data, N, 0);
    init_1D(kernel, N, 0);
    init_1D(result, N, 0);

    printf(">| Data: \n");
    data[0] = 2;
    data[1] = 1;
    for(j = 0; j < N; j ++)
    {
        //data[j] = ( 1 + (nttw_integer)j )%MODULUS;
        printf(FORMAT_CSVOUTPUT_STR,data[j]);
    }
    printf("\n");

    printf(">| Kernel: \n");
    kernel[1] = 2;
    kernel[2] = 1;
    for(j = 0; j < N; j ++)
    {
        //kernel[j] = ( 1 + (nttw_integer)j )%MODULUS;
        printf(FORMAT_CSVOUTPUT_STR,kernel[j]);
    }
    printf("\n");

    //-----------------------------------------
    ///FNTT
    ///MODULUS and PROOT and defined in number.h
    printf(">| FNTT ----------- \n");

    START_TIMER;

    ///FNTT
    fntt(data,N,PROOT,NTTW_FORWARD);
    fntt(kernel,N,PROOT,NTTW_FORWARD);

    STOP_TIMER;
	  duration1 = MICROSECONDS_ELAPSED;

    printf(">| FNTT Data: \n");
    for(j = 0; j < N; j ++)
        printf(FORMAT_CSVOUTPUT_STR,data[j]);
    printf("\n\n");
    printf(">| FNTT Kernel: \n");
    for(j = 0; j < N; j ++)
        printf(FORMAT_CSVOUTPUT_STR,kernel[j]);
    printf("\n\n");

    RESTART_TIMER;

    ///Convolve
    for(j = 0; j < N; j ++)
        result[j] = ((nttw_big_integer)data[j]*kernel[j])%MODULUS;

    STOP_TIMER;
	  duration2 = MICROSECONDS_ELAPSED;

	  printf(">| Convolution (FNTT Space): \n");
    for(j = 0; j < N; j ++)
        printf(FORMAT_CSVOUTPUT_STR,result[j]);
    printf("\n\n");

    RESTART_TIMER;

    ///iFNTT
    fntt(result,N,PROOT,NTTW_INVERSE);

    STOP_TIMER;
	  duration3 = MICROSECONDS_ELAPSED;

    ///Norm
    ntt_norm(result,N);

    printf(">| Convolution iFNTT: \n");
    for(j = 0; j < N; j ++)
        printf(FORMAT_CSVOUTPUT_STR,result[j]);
    printf("\n\n");

    printf(">| FNTTs Took %llu usecs.\n", duration1);
    printf(">| Convolution Took %llu usecs.\n", duration2);
    printf(">| iFNTT Took %llu usecs.\n", duration3);

    //-----------------------------------------
    printf(">| Deallocating memory\n");
    free_array(data);
    free_array(kernel);
    free_array(result);
    free_array(dataBig);
    free_array(kernelBig);
    free_array(resultBig);
    printf(">| Operation Complete\n");
    return EXIT_SUCCESS;
}
