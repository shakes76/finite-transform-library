/**
 * NTTW FNTT 2D Test Program
 * \file fntt_2D_test.c
 * \brief FNTT 2D Test Program for the NTTW C Library.
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
 * \author Shekhar S. Chandra, 2008-9
*/
#include <stdio.h>

#include "nttw/array.h"
#include "nttw/number.h"
#include "nttw/image.h"

const char *filename = "fntt_2D.pgm";

int main(int argc, char *argv[])
{
    nttw_integer *data, *result;
    size_t j, k, N = 0, size; //Use size_t instead of int for portability
    int outputToStd = TRUE;

    printf(">| FNTT 2D Test\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(data));
    if(argc != 2)
    {
        printf(">| %s <n>\n",argv[0]);
        printf(">| Computes the FNTT of size n, where n is dyadic (16,32,etc.).\n");
        return EXIT_FAILURE;
    }
    N = atoi(argv[1]);
    size = N*N;

    if(N > 64)
        outputToStd = FALSE;

    //-----------------------------------------
    printf(">| Allocating memory\n");
    data = array_1D(size);
    result = array_1D(size);

    if(outputToStd)
        printf(">| Data: \n");
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
        {
            data[j*N+k] = ( 1 + (nttw_integer)j*N + (nttw_integer)k )%MODULUS;
            if(outputToStd)
                printf(FORMAT_CSVOUTPUT_STR,data[j*N+k]);
        }
        if(outputToStd)
            printf("\n");
    }

    //-----------------------------------------
    ///FNTT
    ///MODULUS and PROOT and defined in number.h
    printf(">| FNTT 2D TEST ----------- \n");
    fntt_2D(data,result,N,NTTW_FORWARD);

if(outputToStd)
{
    printf(">| FNTT: \n");
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
            printf(FORMAT_CSVOUTPUT_STR,result[j*N+k]);
        printf("\n");
    }
}

    ///iFNTT
    init_1D(data,size,0);
    fntt_2D(result,data,N,NTTW_INVERSE);

    ///Norm
    ///Data normed with functions

if(outputToStd)
{
    printf(">| iFNTT: \n");
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
            printf(FORMAT_CSVOUTPUT_STR,data[j*N+k]);
        printf("\n");
    }
}
else
{
    printf(">| Writing PGM: %s\n",filename);
    writePGM(data,N,N,255,filename,FALSE);
}

    //-----------------------------------------
    printf(">| Deallocating memory\n");
    free_array(data);
    free_array(result);
    printf(">| Operation Complete\n");
    return EXIT_SUCCESS;
}
