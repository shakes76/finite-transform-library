/**
 * NTTW Array Module
 * \file array.c
 * \brief Array Source/Object for the NTTW C Library.
 *
 * This file implements the functions that wrap all the malloc and free functions for producing arrays for 1D/2D with aligned memory.
 * Functions for initializing arrays are also provided.
 *
 * This file is part of NTTW Library.
 *
 * NTTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NTTW is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NTTW. If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2008-9
*/
#include "array.h"

void array(void **pointer, const size_t nbytes)
{
    *pointer = malloc(nbytes);

    if(*pointer == NULL)
    {
        //fprintf(stderr,"ERROR: Failed Memory Allocation.\n");
        exit(EXIT_FAILURE);
    }
}

nttw_integer* array_1D(const size_t size)
{
    nttw_integer *data = NULL;

    data = (nttw_integer *)malloc(size*sizeof(nttw_integer));

    if(data == NULL)
    {
        //fprintf(stderr,"ERROR: Failed Memory Allocation.\n");
        exit(EXIT_FAILURE);
    }

    return data;
}

long* arraySigned_1D(const size_t size)
{
    long *data = NULL;

    data = (long *)malloc(size*sizeof(long));

    if(data == NULL)
    {
        //fprintf(stderr,"ERROR: Failed Memory Allocation.\n");
        exit(EXIT_FAILURE);
    }

    return data;
}

unsigned char* arrayUChar_1D(const size_t size)
{
    unsigned char *data = NULL;

    data = (unsigned char *)malloc(size*sizeof(unsigned char));

    if(data == NULL)
    {
        //fprintf(stderr,"ERROR: Failed Memory Allocation.\n");
        exit(EXIT_FAILURE);
    }

    return data;
}

nttw_big_integer* array_1D_big(const size_t size)
{
    nttw_big_integer *data = NULL;

    data = (nttw_big_integer *)malloc(size*sizeof(nttw_big_integer));

    if(data == NULL)
    {
        //fprintf(stderr,"ERROR: Failed Memory Allocation.\n");
        exit(EXIT_FAILURE);
    }

    return data;
}

void init_1D(nttw_integer *data, const size_t size, const nttw_integer value)
{
    size_t j;

    for (j = 0; j < size; j ++)
        data[j] = value;
}

void initSigned_1D(long *data, const size_t size, const long value)
{
    size_t j;

    for (j = 0; j < size; j ++)
        data[j] = value;
}

void initUChar_1D(unsigned char *data, const size_t size, const unsigned char value)
{
    size_t j;

    for (j = 0; j < size; j ++)
        data[j] = value;
}

void init_1D_big(nttw_big_integer *data, const size_t size, const nttw_big_integer value)
{
    size_t j;

    for (j = 0; j < size; j ++)
        data[j] = value;
}

void free_array(void *pointer)
{
    free(pointer);
}
