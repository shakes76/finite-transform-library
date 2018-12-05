/**
 * FRTW Ghosts Library
 * \file ghosts.c
 * \brief Finite Ghosts Source for the FRTW C Library.
 *
 * This source provides all the operations related to Finite Ghosts.
 * This includes convolution and rows-difference algorithms.
 *
 * This file is part of the FRTW Library.
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
 * along with FRTW. If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2008-10
*/
#include <nttw/image.h>
#include <nttw/prime.h>
#include <nttw/number.h>
//FRTW Library
#include "ghosts.h"
#include "radon.h"

size_t determineMissingAngles(const size_t N, nttw_big_integer *finiteAngles, nttw_integer *perpFlags, const size_t noOfProjs, nttw_big_integer **missingAngles, nttw_big_integer **missingAngleLookup)
{
    ///Dec Other variables
    size_t n, count, noOfGhosts, size;
    int dyadicSize = TRUE;

    if(N % 2 == 1) ///Assume prime if odd
        dyadicSize = FALSE;

    if(dyadicSize)
        size = N+N/2;
    else
        size = N+1;

    noOfGhosts = size-noOfProjs;

    ///Allocate
    *missingAngles = array_1D_big(noOfGhosts);
    *missingAngleLookup = array_1D_big(size);

    init_1D_big(*missingAngles, noOfGhosts, 0);
    init_1D_big(*missingAngleLookup, size, TRUE);

    for(n = 0; n < noOfProjs; n ++) ///Build angle lookup
    {
        if(perpFlags[n])
            (*missingAngleLookup)[ finiteAngles[n]+N ] = FALSE;
        else if(finiteAngles[n] > 0)
            (*missingAngleLookup)[ N-finiteAngles[n] ] = FALSE; ///Flip m's for reverse shifts (N- term)
        else
            (*missingAngleLookup)[ finiteAngles[n] ] = FALSE;
    }

    for(n = 0, count = 0; n < size; n ++)
    {
        if((*missingAngleLookup)[n])
        {
            (*missingAngles)[count] = n;
            count ++;
        }
    }

    return noOfGhosts;
}

nttw_integer* ghostEigenvalues(const size_t N, nttw_integer *root, nttw_integer *modulus, nttw_integer *proot)
{
    ///Constants
    const size_t size = N*N;
    const int normTransform = TRUE;

    ///Dec Ghost arrays
    nttw_integer *operators, *eigenvalues, *opRowPtr, *eigenRowPtr;

    ///Dec Other variables
    size_t j, index;

    ///Determine prime length parameters
    findPrimeLengthParameters(N, root, modulus, proot); //!< Reads file!

    ///Allocate arrays
    operators = array_1D(size);
    eigenvalues = array_1D(size);

    ///Determine Eigenvalues
    init_1D(operators,size,0); //Set to zero
    for(j = 0; j < N-1; j ++) //!< Compute 0 <= m < p operators
    {
        opRowPtr = &operators[j*N];
        eigenRowPtr = &eigenvalues[j*N];

        ///Setup Dilation Bases
        opRowPtr[0] += 1;
        opRowPtr[(j+1)] += *modulus-1;

        ///Find Eigenvalues of Operator
        fntt_prime(opRowPtr,eigenRowPtr,N,*root,*modulus,*proot,NTTW_FORWARD,normTransform);
    }
    ///Compute m = p operator
    index = N-1;
    opRowPtr = &operators[index*N];
    opRowPtr[0] += 1;
    opRowPtr[1] += *modulus-1;
    eigenRowPtr = &eigenvalues[index*N];
    ///Find Eigenvalues of m = p Operator
    fntt_prime(opRowPtr,eigenRowPtr,N,*root,*modulus,*proot,NTTW_FORWARD,normTransform);

    ///Deallocate
    free_array(operators);

    return eigenvalues;
}

nttw_integer* convolveGhosts_1D(const size_t N, nttw_integer* eigenvalues, const nttw_integer root, const nttw_integer modulus, const nttw_integer proot, nttw_big_integer *missing, nttw_big_integer *missingLookup, const size_t noOfGhosts)
{
    ///Constants
    const size_t size = N*N;
    const int normTransform = TRUE;

    ///Dec Ghost Arrays
    nttw_integer *deltaFunc, *deltaFuncNTT, *ghosts, *ghostsImage, *ghostRowPtr, *eigenRowPtr;

    ///Dec Other variables
    size_t j, k, l, index;

    ///Allocate arrays
    deltaFunc = array_1D(N);
    deltaFuncNTT = array_1D(N);
    ghosts = array_1D(size+N);
    ghostsImage = array_1D(size);

    ///Construct delta function
    init_1D(deltaFunc, N, 0);
    deltaFunc[0] += 1;

    ///FNTT of delta function
    fntt_prime(deltaFunc,deltaFuncNTT,N,root,modulus,proot,NTTW_FORWARD,normTransform);;

    ///Convolve ghosts in FRT space
    init_1D(ghosts,size+N,0);

    for(j = 0; j < N; j ++) //!< For each known projection (0 <= m < p)
    {
        if(missingLookup[j]) //!<So if ghost projection skip
            continue;

        ghostRowPtr = &ghosts[j*N];
        for(l = 0; l < N; l ++) //!<Copy NTT of Delta function to known projection
            ghostRowPtr[l] = deltaFuncNTT[l];

        for(k = 0; k < noOfGhosts; k ++) //!<For each ghost (0 <= m < p)
        {
            index = missing[k]-j-1; //!< Determine 1D shift
            index += N;
            index %= N; //!< Make shift positive
            eigenRowPtr = &eigenvalues[index*N];
            for(l = 0; l < N; l ++) //!<For each bin
                ghostRowPtr[l] = (ghostRowPtr[l] * (nttw_big_integer)eigenRowPtr[l])%modulus;
        }
    }
    ///Known m = p projection
    index = N;
    ghostRowPtr = &ghosts[index*N];
    for(l = 0; l < N; l ++) //!<Copy NTT of Delta function to m = p projection
        ghostRowPtr[l] = deltaFuncNTT[l];
    index = N-1;
    eigenRowPtr = &eigenvalues[index*N]; //!< All ghosts have same m = p 1D ghost operator
    for(k = 0; k < noOfGhosts; k ++) //!<For each ghost
    {
        for(l = 0; l < N; l ++) //!<For each bin
            ghostRowPtr[l] = (ghostRowPtr[l] * (nttw_big_integer)eigenRowPtr[l])%modulus;
    }

    inrt2(ghosts,ghostsImage,N,root,modulus,proot,TRUE); ///iNRT

    ///Deallocate
    free_array(deltaFunc);
    free_array(deltaFuncNTT);
    free_array(ghosts);

    return ghostsImage;
}

nttw_integer* deghost(const size_t N, nttw_integer* frtSpace, nttw_integer* ghost, const size_t noOfGhosts)
{
    ///Constants
    const size_t size = N*N;
    const int normTransform = FALSE;

    ///Dec Ghost Arrays
    nttw_integer *nttOfFRT, *eigenvalues, *imageEigenvalues, *deghostedImage, *ghostRowPtr, *eigenRowPtr;

    ///Dec Row arrays
    nttw_integer *deghostRow, *conRow, *subConRow;

    ///Dec Other variables
    size_t j, k;
    nttw_integer root, modulus, proot;
    nttw_big_integer inv, Isum = 0;
    int row; //Needs to be signed

    ///Allocate
    nttOfFRT = array_1D(size+N);
    eigenvalues = array_1D(size);
    imageEigenvalues = array_1D(size);
    deghostedImage = array_1D(size);
    conRow = array_1D(N);
    subConRow = array_1D(N);

    ///Determine prime length parameters
    findPrimeLengthParameters(N, &root, &modulus, &proot); //!< Reads file!

    ///Compute Multiplicative inverse of N mod p'
    inv = minverse_big(N, modulus);

    ///Ensure FRT space is modulus and take NTTs
    for(j = 0; j < N+1; j ++)
    {
        for(k = 0; k < N; k ++)
            frtSpace[j*N + k] %= modulus;

        ghostRowPtr = &frtSpace[j*N];
        eigenRowPtr = &nttOfFRT[j*N];

        fntt_prime(ghostRowPtr,eigenRowPtr,N,root,modulus,proot,NTTW_FORWARD,FALSE);
    }
    printf("Isum: %lu\n", nttOfFRT[0]); //First term in NTT of projection is DC

    ///Invert FRT space to get ghosted image
    Isum = inrt2(nttOfFRT,deghostedImage,N,root,modulus,proot,FALSE); ///iNRT
//    printf("Isum: %lu\n", Isum);
//    printf("Image of FRT Space:\n");
    for(j = 0; j < N; j ++) //Norm by p
    {
        for(k = 0; k < N; k ++)
        {
            ///Undo Isum subtraction in inrt2 routine
            deghostedImage[j*N + k] = MODADD(deghostedImage[j*N + k],Isum,modulus);
            ///Norm by p
            deghostedImage[j*N + k] *= inv;
            deghostedImage[j*N + k] %= modulus;
            ///Minus Isum
            deghostedImage[j*N + k] = MODSUB(deghostedImage[j*N + k],Isum,modulus);
            ///Norm by p
//            deghostedImage[j*N + k] *= inv;
//            deghostedImage[j*N + k] %= modulus;
//            printf("%lu, ", deghostedImage[j*N + k]);
        }
//        printf("\n");
    }

    ///Determine ghost row eigenvalues
//    printf("Ghost Space:\n");
    for(j = 0; j < N; j ++)
    {
        ghostRowPtr = &ghost[j*N];
        eigenRowPtr = &eigenvalues[j*N];

//        for(k = 0; k < N; k ++)
//            printf("%lu, ", ghostRowPtr[k]);
//        printf("\n");

        ///Find Eigenvalues of Operator
        fntt_prime(ghostRowPtr,eigenRowPtr,N,root,modulus,proot,NTTW_FORWARD,normTransform);
    }

    ///Determine Image row eigenvalues
//    printf("Image Space:\n");
    for(j = 0; j < N; j ++)
    {
        ghostRowPtr = &deghostedImage[j*N];
        eigenRowPtr = &imageEigenvalues[j*N];

//        for(k = 0; k < N; k ++)
//            printf("%lu, ", ghostRowPtr[k]);
//        printf("\n");

        ///Find Eigenvalues of Operator
        fntt_prime(ghostRowPtr,eigenRowPtr,N,root,modulus,proot,NTTW_FORWARD,normTransform);
    }

    ///Deconvolve Ghosts
    init_1D(deghostedImage, size, 0);
    init_1D(conRow, N, 0);
    init_1D(subConRow, N, 0);

    for(row = N-noOfGhosts-1; row >= 0; row --) //!<Row to deghost
    {
        deghostRow = &deghostedImage[row*N];

        for(j = 0; j < noOfGhosts+1; j ++) ///Deghost current row
        {
            ///Deconvolve
            for(k = 0; k < N; k ++)
                conRow[k] = ( imageEigenvalues[(j+row)*N + k] * (nttw_big_integer)eigenvalues[j*N + k] )%modulus; //!< Need big int here?

            ///Accumulate result (using sub-circs)
            if(j > 0)
            {
                for(k = 0; k < N; k ++)
                    subConRow[k] = MODADD(subConRow[k],conRow[k],modulus);
            }
            else
            {
                for(k = 0; k < N; k ++)
                    subConRow[k] = conRow[k];
            }
        }

        fntt_prime(subConRow,deghostRow,N,root,modulus,proot,NTTW_INVERSE,TRUE);

//        printf("deghostRow %lu: \n", row);
        for(k = 0; k < N; k ++)
        {
            deghostRow[k] *= inv;
            deghostRow[k] %= modulus;
//            printf("%lu, ", deghostRow[k]);
        }
//        printf("\n");

        ///Backsubstitute result
        for(k = 0; k < N; k ++)
            imageEigenvalues[row*N + k] = MODSUB(imageEigenvalues[row*N + k], subConRow[k], modulus);
    }

    ///Deallocate
    free_array(nttOfFRT);
    free_array(eigenvalues);
    free_array(imageEigenvalues);
    free_array(conRow);
    free_array(subConRow);

    return deghostedImage;
}
