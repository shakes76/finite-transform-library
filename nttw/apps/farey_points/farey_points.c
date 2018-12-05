/**
 * Farey Points within an NxN lattice space
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

#include "nttw/prime.h"
#include "nttw/timing.h"
#include "nttw/image.h"

int main(int argc, char *argv[])
{
	nttw_integer *lattice;
	size_t j, k, N = 0, size; //Use size_t instead of int for portability

	printf(">| Farey Points\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(lattice));
    if(argc != 3)
    {
        printf(">| %s <N> <output>\n",argv[0]);
        printf(">| Computes the Farey Points within Lattice of size NxN.\n");
        printf(">| Output will be a PGM file. Indices are zero-indexed.\n");
        return EXIT_FAILURE;
    }
    N = atoi(argv[1]);
    size = N*N;

    lattice = array_1D(size);
    init_1D(lattice,size,0);

	START_TIMER;
	for(j = 1; j < N; j ++)
        for(k = 1; k < N; k ++)
            if(gcd(j,k) == 1)
                lattice[j*N+k] = 1;
    STOP_TIMER;

	printf(">| Took %llu usecs.\n",MICROSECONDS_ELAPSED);

	writePGM(lattice,N,N,1,argv[2],FALSE);

    free_array(lattice);
	return 0;
}
