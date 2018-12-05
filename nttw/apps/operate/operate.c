/**
 * Operate (in the number field)
 * Computes the arithmetic of two images in the number field
 * Outputs image result.
 * \file operate.c
 * \brief Image arithmetic Program for the NTTW C Library.
 *
 * This file is part of NTTW Library.
 *
 * NTTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NTTW is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with NTTW.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2010-13
*/
//NTTW
#include "nttw/timing.h"
#include "nttw/prime.h"
#include "nttw/number.h"
#include "nttw/image.h"

///Global
const int totalRows = 1000;
const int totalCols = 10;
const char *filename = "primes.csv";

int main(int argc, char *argv[])
{
    nttw_integer *image, *operand, *result, *primes, root = 0, modulus = 0, proot = 0;
    int j, k, rows1, cols1, rows2, cols2, binaryFile = FALSE, dyadicTransform = TRUE;
    unsigned long long duration = 0;

    printf(">| Image Operation Program.\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %lu bits\n",BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 4)
    {
        printf(">| Usage: %s <filename 1> <filename 2> <output>\n",argv[0]);
        printf(">| filename 1 & 2 are loaded and operated on to produce output.\n");
        printf(">| Addition only supported atm.\n");
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }

    //--------------------------------------------
    ///Load Image
    if(!readPGM(&image,&rows1,&cols1,argv[1],binaryFile))
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }
    if(!readPGM(&operand,&rows2,&cols2,argv[2],binaryFile))
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }
    if(rows1 % 2 == 1) ///Odd #rows then assume prime-length
    {
        dyadicTransform = FALSE;
        ///Read primes list
        if( !readCSV(&primes,totalRows,totalCols,filename) )
        {
            printf(">| Could not open primes file.\n");
            exit(EXIT_FAILURE);
        }
        ///Compute roots
        root = findFirstPrimitiveRoot(rows1);
        modulus = findAlternatePrime(primes,totalRows*totalCols,rows1);
        proot = findFirstPrimitiveRoot(modulus);
        printf(">| Root: %u, Prime': %u, Primitive Root: %u.\n",root,modulus,proot);
    }

    ///Output parameters to be used
    result = array_1D(rows1*cols1);

    //--------------------------------------------
    ///Fast FRT
    fprintf(stderr,">| Operating... ");
    START_TIMER;

    if(dyadicTransform)
    {
        for(j = 0; j < rows1; j ++)
            for(k = 0; k < cols1; k ++)
                result[j*rows1+k] = MODADD(image[j*rows1+k], operand[j*rows1+k], MODULUS);
    }
    else
    {
        for(j = 0; j < rows1; j ++)
            for(k = 0; k < cols1; k ++)
                result[j*rows1+k] = MODADD(image[j*rows1+k], operand[j*rows1+k], modulus);
    }

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);

    ///Save Result
    writePGM(result,rows1,cols1,255,argv[3],binaryFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(image);
    free_array(operand);
    free_array(result);
    if(!dyadicTransform)
        free_array(primes);

    return EXIT_SUCCESS;
}
