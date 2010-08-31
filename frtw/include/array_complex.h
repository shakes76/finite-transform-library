/**
 * FRTW Complex Valued Array Library
 * \file array_complex.h
 * \brief Complex Array Header/Object for the FRTW C Library.
 *
 * This wraps all the FFTW malloc and free functions for producing arrays for 1D/2D.
 * Functions for copying and differencing arrays are included also.
 *
 * This file is part of FRTW Library.
 *
 * FRTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FRTW is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FRTW. If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2008-10
*/
#ifndef ARRAY_COMPLEX_H_INCLUDED
#define ARRAY_COMPLEX_H_INCLUDED

#include <stdlib.h>
#include <fftw3.h>
#include <nttw/global.h>
/**
* \mainpage FRTW C Library API
*
* \section intro_sec Introduction
*
* The various Discrete Radon Transform Algorithms Optimised for C and performance. Array, Simple Timing and Imaging (PGM and CSV) modules are also provided via the NTTW library.
* This library requires that the sister NTTW library be installed.
*
* \section features Features
*
* Features an Array module, a micro-second Timing module an PGM Imaging module via NTTW.  Fast Fourier Transforms module is provided via the external FFTW library.
* The Discrete Radon Transform and its fast implementation are provided for dyadic and prime lengths. The Mojette Transform and the Fast Mojette algorithm is also implemented (see "Fast Mojette Transform for Discrete Tomography" by Chandra et al.).
*/

/**
 * \defgroup FFTW_Arrays FFTW (Complex) Array Objects
*/
//@{
/**
 * \fn fftw_array(const size_t size)
 * \brief Forms a FFTW complex array of size
 * \return Pointer to the first element of the array.
 *
 * Constructs a FFTW complex array.
 * The fftw_malloc function is used to allocate the memory and garbage collection
 * is NOT automated, please use fftw_free() from <fftw3.h> to deallocate the memory.
*/
NTTW_DLL_SYM fftw_complex* fftw_array(const size_t size);
/**
 * \fn fftw_init_1D(fftw_complex *data, const size_t size, const nttw_real value)
 * \brief Initializes FFTW array (real and imaginary parts) to value
*/
NTTW_DLL_SYM void fftw_init_1D(fftw_complex *data, const size_t size, const nttw_real value);

/**
 * \fn array_to_fftw_array(nttw_integer *source, fftw_complex *dest, const int size)
 * \brief Copies the array to the real FFTW array components
 *
 * IMPORTANT: Only the Real part of the FFTW is copied!
*/
NTTW_DLL_SYM void array_to_fftw_array(nttw_integer *source, fftw_complex *dest, const int size);
/**
 * \fn arraySigned_to_fftw_array(long *source, fftw_complex *dest, const int size)
 * \brief Copies the signed array to the real FFTW array components
 *
 * IMPORTANT: Only the Real part of the FFTW is copied!
*/
NTTW_DLL_SYM void arraySigned_to_fftw_array(long *source, fftw_complex *dest, const int size);
/**
 * \fn fftw_array_to_array(fftw_complex *source, nttw_integer *dest, const int size)
 * \brief Copies the real FFTW array components to array
 *
 * IMPORTANT: Only the Real part of the FFTW is copied!
*/
NTTW_DLL_SYM void fftw_array_to_array(fftw_complex *source, nttw_integer *dest, const int size);
/**
 * \fn fftw_array_to_arraySigned(fftw_complex *source, long *dest, const int size)
 * \brief Copies the real FFTW array components to signed array
 *
 * IMPORTANT: Only the Real part of the FFTW is copied!
*/
NTTW_DLL_SYM void fftw_array_to_arraySigned(fftw_complex *source, long *dest, const int size);
//@}

/**
 * \defgroup Two_Dimensional_Pointer_Arrays Two Dimensional (2D) Pointer Array Objects
 */
//@{
/**
 * \fn ptrArray(const size_t size)
 * \brief Constructs a pointer array of type nttw_integer* which is accessible from
 * the zeroth index.
*/
NTTW_DLL_SYM nttw_integer** ptrArray(const size_t size);
/**
 * \fn free_ptrArray(nttw_integer **data)
 * \brief Deletes or frees the pointer array that would have been created by ptrArray()
 *
 * Must be used with ptrArray() that allocates a pointer array. The pointers themselves
 * are NOT deallocated. Use free() from <stdlib.h> to free them.
*/
NTTW_DLL_SYM void free_ptrArray(nttw_integer **data);
/**
 * \fn free_ptrArray_All(nttw_integer **data, const size_t size)
 * \brief Deletes or frees the pointer array and subsequent 1-D arrays that each of the pointers
 * are pointing too that would have been created by ptrArray().
 *
 * Must be used with ptrArray() that allocates a pointer array. The pointers in the array
 * are deallocated.
*/
NTTW_DLL_SYM void free_ptrArray_All(nttw_integer **data, const size_t size);
//@}
#endif // ARRAY_COMPLEX_H_INCLUDED
