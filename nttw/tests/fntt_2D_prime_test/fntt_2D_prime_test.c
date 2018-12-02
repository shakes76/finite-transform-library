/**
 * NTTW FNTT 2D Prime length Test Program
 * \file fntt_2D_prime_test.c
 * \brief FNTT 2D Prime length Test Program for the NTTW C Library.
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
#include "nttw/array.h"
#include "nttw/image.h"
#include "nttw/prime.h"
#include "nttw/number.h"

///Globals
const int totalRows = 1000;
const int totalCols = 10;
const char *filename = "primes.csv";

int main(int argc, char *argv[])
{
    nttw_integer *data, *result;
    size_t j, k, N = 0, size; //Use size_t instead of int for portability
    nttw_integer *primes, root, modulus, proot;
    const int normTransform = TRUE;

    printf(">| FNTT 2D (Prime length) Test\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(data));
    if(argc != 2)
    {
        printf(">| %s <n>\n",argv[0]);
        printf(">| Computes the FNTT of size n, where n is prime (13,17,etc.).\n");
        return EXIT_FAILURE;
    }
    N = atoi(argv[1]);
    size = N*N;

    //-----------------------------------------
    printf(">| Allocating memory\n");
    data = array_1D(size);
    result = array_1D(size);

    printf(">| Data: \n");
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
        {
            data[j*N+k] = 1 + (nttw_integer)j*N + (nttw_integer)k;
            printf("%u, ",data[j*N+k]);
        }
        printf("\n");
    }

    //-----------------------------------------
    ///Determine roots for the prime
    ///Read primes list
    if( !readCSV(&primes,totalRows,totalCols,filename) )
    {
        printf(">| Could not open primes file.\n");
        exit(EXIT_FAILURE);
    }
    ///Compute roots
    root = findFirstPrimitiveRoot(N);
    modulus = findAlternatePrime(primes,totalRows*totalCols,N);
    proot = findFirstPrimitiveRoot(modulus);
    printf(">| Root: %u, Prime': %u, Primitive Root: %u.\n",root,modulus,proot);

    //-----------------------------------------
    ///FNTT
    ///MODULUS and PROOT and defined in number.h
    printf(">| FNTT 2D TEST ----------- \n");
    fntt_2D_prime(data,result,N,root,modulus,proot,NTTW_FORWARD,normTransform);

    printf(">| FNTT: \n");
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
            printf("%u, ",result[j*N+k]);
        printf("\n");
    }

    ///iFNTT
    fntt_2D_prime(result,data,N,root,modulus,proot,NTTW_INVERSE,normTransform);

    ///Norm
    ///Data normed with functions

    printf(">| iFNTT: \n");
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
            printf("%u, ",data[j*N+k]);
        printf("\n");
    }

    //-----------------------------------------
    printf(">| Deallocating memory\n");
    free_array(primes);
    free_array(data);
    free_array(result);
    printf(">| Operation Complete\n");
    return EXIT_SUCCESS;
}
