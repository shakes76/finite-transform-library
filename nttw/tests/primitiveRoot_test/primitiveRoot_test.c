/**
 * NTTW Prime' Test Program
 * \file primitiveRoot_test.c
 * \brief Prime' Test Program for the NTTW C Library.
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
#include "nttw/image.h"
#include "nttw/prime.h"

const size_t totalRows = 1000; //1000
const size_t totalCols = 10;
const char *filename = "primes.csv";
const char *primeDash_filename = "primeDash.csv";

int main(int argc, char *argv[])
{
    nttw_integer *primes;
    const size_t totalPrimes = totalRows*totalCols, upto = 2;
    nttw_integer proot, modulus;

    printf(">| Prime' Test\n");
    printf(">| Copyright Shekhar Chandra, 2013\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Library Integer Size of %u bits\n",sizeof(nttw_integer)*8);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(primes));
  
    if(argc != 2)
    {
        printf(">| %s <modulus>\n",argv[0]);
        printf(">| Computes the primitive root of modulus.\n");
        return EXIT_FAILURE;
    }
    modulus = atol(argv[1]);

    //-----------------------------------------
    ///Load primes list
    if(!readCSV(&primes,totalRows,totalCols,filename))
    {
        fprintf(stderr,">| Could not open primes list. Exiting.\n");
        exit(EXIT_FAILURE);
    }
    printf(">| Primes List loaded.\n");

    //-----------------------------------------
    ///Determine prime' values
    printf(">| Modulus: %lu\n",modulus);
    ///Find primitive root
    proot = findFirstPrimitiveRoot_Full(primes,totalPrimes,modulus);
    printf(">| Primitive Root: %u\n",proot);

    //-----------------------------------------
    printf(">| Deallocating memory\n");
    free_array(primes);
    printf(">| Operation Complete\n");
    return EXIT_SUCCESS;
}
