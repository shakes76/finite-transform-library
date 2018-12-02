/**
 * NTTW PGM Read/Write Test Program
 * \file pgm_test.c
 * \brief PGM Read/Write Test Program for the NTTW C Library.
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

const int binary = FALSE;
const char *filename = "pgm_test.pgm";

int main(int argc, char *argv[])
{
    nttw_integer *data, *result;
    size_t j, k, N = 0, size; //Use size_t instead of int for portability
    int readRows = 0, readCols = 0;

    printf(">| PGM Read/Write Test\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(data));
    if(argc != 2)
    {
        printf(">| %s <n>\n",argv[0]);
        printf(">| Allocates, Prints and Deallocates array of size n\n");
        return EXIT_FAILURE;
    }
    N = atoi(argv[1]);
    size = N*N;

    //-----------------------------------------
    printf(">| Allocating memory\n");
    data = array_1D(size);

    printf(">| Data: \n");
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
        {
            data[j*N+k] = 1 + (nttw_integer)j*N + (nttw_integer)k;
            printf(FORMAT_CSVOUTPUT_STR,data[j*N+k]);
        }
        printf("\n");
    }

    printf(">| Writing data\n");
    writePGM(data,N,N,N*N,filename,binary);
    printf(">| Reading  data\n");
    readPGM(&result,&readRows,&readCols,filename,binary);

    printf(">| Result: \n");
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
            printf(FORMAT_CSVOUTPUT_STR,result[j*N+k]);
        printf("\n");
    }

    //-----------------------------------------
    printf(">| Deallocating memory\n");
    free_array(data);
    free_array(result);
    printf(">| Operation Complete\n");
    return EXIT_SUCCESS;
}
