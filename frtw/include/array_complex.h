/**
 * FRTW Complex Valued Array Library
 * \file array_complex.h
 * \brief Complex Array Header/Object for the FRTW C Library.
 *
 * This wraps all the FFTW malloc and free functions for producing arrays for 1D/2D.
 * Functions for copying and differencing arrays are included also.
 *
 * \author Shekhar S. Chandra, 2008-9
*/
#ifndef ARRAY_COMPLEX_H_INCLUDED
#define ARRAY_COMPLEX_H_INCLUDED

#include <fftw3.h>
#include <nttw/global.h>
/**
* \mainpage FRTW C Library API
*
* \section intro_sec Introduction
*
* The various Discrete Radon Transform Algorithms Optimised for C and performance. Array, Timing and Imaging modules are also provided via the NTTW library. 
* This library requires that the sister NTTW library be installed
*
* \section features Features
*
* Features an Array module, a micro-second Timing module an Imaging module via NTTW.  Fast Fourier Transforms module is provided via the external FFTW library. 
* The Discrete Radon Transform and its fast implementation are provided for dyadic and prime lengths.
*/

/**
 * \defgroup FFTW_Arrays FFTW Array Objects
*/
//@{
/**
 * \fn fftw_array(const int size)
 * \brief Forms a FFTW complex array of size
 * \return Pointer to the first element of the array.
 *
 * Constructs a FFTW complex array.
 * The fftw_malloc function is used to allocate the memory and garbage collection
 * is NOT automated, please use fftw_free() from <fftw3.h> to deallocate the memory.
*/
NTTW_DLL_SYM fftw_complex* fftw_array(const size_t size);
/**
 * \fn fftw_init_1D(fftw_complex *data, const int size, const nttw_real value)
 * \brief Initializes FFTW array (real and imaginary parts) to value
*/
NTTW_DLL_SYM void fftw_init_1D(fftw_complex *data, const size_t size, const nttw_real value);

/**
 * \fn array_to_fftw_array(nttw_integer *source, fftw_complex *dest, const int size)
 * \brief Copies the real FFTW array components to array
 *
 * IMPORTANT: Only the Real part of the FFTW is copied!
*/
NTTW_DLL_SYM void array_to_fftw_array(nttw_integer *source, fftw_complex *dest, const int size);
/**
 * \fn fftw_array_to_array(fftw_complex *source, nttw_integer *dest, const int size)
 * \brief Copies the real FFTW array components to array
 *
 * IMPORTANT: Only the Real part of the FFTW is copied!
*/
NTTW_DLL_SYM void fftw_array_to_array(fftw_complex *source, nttw_integer *dest, const int size);
//@}
#endif // ARRAY_COMPLEX_H_INCLUDED
