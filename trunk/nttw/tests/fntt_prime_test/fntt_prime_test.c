/**
 * NTTW FNTT Prime Test Program
 * \file fntt_prime_test.c
 * \brief FNTT Prime Test Program for the NTTW C Library.
 *
 * This file is part of NTTW Library.
 *
 * NTTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DGV is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NTTW.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2008-9
*/
#include "array.h"
#include "prime.h"
#include "image.h"
#include "number32.h"

const int totalRows = 1000;
const int totalCols = 10;
const char *filename = "primes.csv";

int main(int argc, char *argv[])
{
    nttw_integer *data, *result, *inversion, *primes;
    size_t j, N = 0; //Use size_t instead of int for portability
    nttw_integer proot, root, modulus;
    int normTransform = TRUE;

    printf(">| FNTT Prime Test\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Library Integer Size of %u bits\n",sizeof(nttw_integer)*8);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(data));
    if(argc != 2)
    {
        printf(">| %s <n>\n",argv[0]);
        printf(">| Computes the FNTT of size n, where n is prime (13,17,etc.).\n");
        return EXIT_FAILURE;
    }
    N = atoi(argv[1]);

    //-----------------------------------------
    ///Find primitive root
    root = findFirstPrimitiveRoot(N);
    printf(">| Length N: %u has root of %i.\n",N,root);

    if(!readCSV(&primes,totalRows,totalCols,filename))
    {
        printf(">| Could not open primes list. Exiting.\n");
        exit(EXIT_FAILURE);
    }
    printf(">| Prime List loaded.\n");
    modulus = findAlternatePrime(primes,totalRows*totalCols,N);
    printf(">| Found alternate prime.\n");
    proot = findFirstPrimitiveRoot(modulus);
    printf(">| Found primitive root.\n");
    printf(">| Prime': %i, Root: %i\n",modulus,proot);

    //-----------------------------------------
    printf(">| Allocating memory\n");
    data = array_1D(N);
    result = array_1D(N);
    inversion = array_1D(N);

    printf(">| Data: \n");
    for(j = 0; j < N; j ++)
    {
        data[j] = ( 1 + (nttw_integer)j )%modulus;
        printf(FORMAT_CSVOUTPUT_STR,data[j]);
    }
    printf("\n");

    //-----------------------------------------
    ///FNTT
    ///MODULUS and PROOT and defined in number.h
    printf(">| FNTT TEST ----------- \n");
    fntt_prime(data,result,N,root,modulus,proot,NTTW_FORWARD,normTransform);

    printf(">| FNTT: \n");
    for(j = 0; j < N; j ++)
        printf(FORMAT_CSVOUTPUT_STR,result[j]);
    printf("\n\n");

    ///iFNTT
    fntt_prime(result,inversion,N,root,modulus,proot,NTTW_INVERSE,normTransform);

    printf(">| iFNTT: \n");
    for(j = 0; j < N; j ++)
        printf(FORMAT_CSVOUTPUT_STR,inversion[j]);
    printf("\n\n");

    //-----------------------------------------
    printf(">| Deallocating memory\n");
    free_array(primes);
    free_array(data);
    free_array(result);
    free_array(inversion);
    printf(">| Operation Complete\n");
    return EXIT_SUCCESS;
}
