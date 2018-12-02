/**
 * FRTW Complex Valued Array Library
 * \file array_complex.c
 * \brief Complex Array Source for the FRTW C Library.
 *
 * This wraps all the FFTW malloc and free functions for producing arrays for 1D/2D.
 * Functions for copying and differencing arrays are included also.
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
 * along with FRTW. If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2008-9
*/
#include "array_complex.h"

//Complex Arrays
fftw_complex* fftw_array(const size_t size)
{
    fftw_complex *data = NULL;

    data = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));

    return data;
}

void fftw_init_1D(fftw_complex *data, const size_t size, const nttw_real value)
{
    size_t j;

    for (j = 0; j < size; j ++)
    {
        data[j][0] = value;
        data[j][1] = value;
    }
}

void array_to_fftw_array(nttw_integer *source, fftw_complex *dest, const int size)
{
    int j;

    for (j = 0; j < size; j++)
        dest[j][0] = (double) source[j]; //[0] is the real part
}

void arraySigned_to_fftw_array(long *source, fftw_complex *dest, const int size)
{
    int j;

    for (j = 0; j < size; j++)
        dest[j][0] = (double) source[j]; //[0] is the real part
}

void fftw_array_to_array(fftw_complex *source, nttw_integer *dest, const int size)
{
    int j;

    for (j = 0; j < size; j++)
        dest[j] = (nttw_integer) source[j][0]; //[0] is the real part
}

void fftw_array_to_arraySigned(fftw_complex *source, long *dest, const int size)
{
    int j;

    for (j = 0; j < size; j++)
    {
        if(source[j][0] < 0)
            dest[j] = (long) (source[j][0] - 0.5); //[0] is the real part
        else
            dest[j] = (long) (source[j][0] + 0.5); //[0] is the real part
    }
}

//2D Ptr Operations
nttw_integer** ptrArray(const size_t size)
{
    nttw_integer **data;

    data = (nttw_integer **)malloc(size*sizeof(nttw_integer *));

    return data;
}

void free_ptrArray(nttw_integer **data)
{
    free(data);
}

void free_ptrArray_All(nttw_integer **data, const size_t size)
{
    size_t j;

    for(j = 0; j < size; j ++)
        free(data[j]);

    free_ptrArray(data);
}
