/**
 * NTTW Array Module
 * \file array.h
 * \brief Array Header/Object for the NTTW C Library.
 *
 * This file defines the functions that wrap all the malloc and free functions for producing arrays for 1D/2D with aligned memory.
 * Functions for initializing arrays are also provided.
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
 * along with NTTW. If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2008-9
*/
#ifndef ARRAY_H_INCLUDED
#define ARRAY_H_INCLUDED

#include <stdlib.h>

#include "global.h"

/**
 * \defgroup Arrays Arrays allocated as One Dimensional (1D) Array Objects
 * \brief Arrays for arbitrary dimensions. Done by constructing 1D arrays. Handles memory alignment also
 *
 * The main function is array_1D() and free array by using free_array(). These functions handle memory alignment independently of platform.
 * \todo Add memory alignment code.
*/
//@{
/**
 * \fn array(void **pointer, const size_t nbytes)
 * \brief Forms an array of size
 * \param pointer Pointer to the array to be allocated. It's a Double pointer because it is passed by reference.
 * \param nbytes The bytes to be allocated. Use the sizeof() function here.
 *
 * Constructs a array for 1D reference. Reference into array as j*columns+k, where j is the current row and k the current column.
 * The malloc function is used to allocate the memory and garbage collection
 * is NOT automated, please use free_array() to deallocate the memory.
*/
NTTW_DLL_SYM void array(void **pointer, const size_t nbytes);

/**
 * \fn array_1D(const size_t size)
 * \brief Forms an 1D array of size
 * \return Pointer to the first element of the array.
 *
 * Constructs a 1-dimensional array of type nttw_integer (defined in global.h).
 * Allocates 16-bit integers in 32-bit mode or 32-bit integers in 64-bit mode.
 * The malloc function is used to allocate the memory and garbage collection
 * is NOT automated, please use free_array() to deallocate the memory.
*/
NTTW_DLL_SYM nttw_integer* array_1D(const size_t size);

/**
 * \fn arraySigned_1D(const size_t size)
 * \brief Forms an 1D signed array of size
 * \return Pointer to the first element of the array.
 *
 * Constructs a 1-dimensional array of type long.
 * Allocates signed integers of maximum machine word size.
 * The malloc function is used to allocate the memory and garbage collection
 * is NOT automated, please use free_array() to deallocate the memory.
*/
NTTW_DLL_SYM long* arraySigned_1D(const size_t size);

/**
 * \fn arrayUChar_1D(const size_t size)
 * \brief Forms an 1D unsigned character array of size
 * \return Pointer to the first element of the array.
 *
 * Constructs a 1-dimensional array of type unsigned char.
 * The malloc function is used to allocate the memory and garbage collection
 * is NOT automated, please use free_array() to deallocate the memory.
*/
NTTW_DLL_SYM unsigned char* arrayUChar_1D(const size_t size);

/**
 * \fn array_1D_big(const size_t size)
 * \brief Forms an 1D array of size
 * \return Pointer to the first element of the array.
 *
 * Constructs a 1-dimensional array of type nttw_integer (defined in global.h).
 * Allocats 32-bit integers in 32-bit mode or 64-bit integers in 64-bit mode.
 * The malloc function is used to allocate the memory and garbage collection
 * is NOT automated, please use free_array() to deallocate the memory.
*/
NTTW_DLL_SYM nttw_big_integer* array_1D_big(const size_t size);

/**
 * \fn init_1D(nttw_integer *data, const size_t size, const nttw_integer value)
 * \brief Initializes the 1D array to value provided.
*/
NTTW_DLL_SYM void init_1D(nttw_integer *data, const size_t size, const nttw_integer value);

/**
 * \fn initSigned_1D(long *data, const size_t size, const long value)
 * \brief Initializes a signed 1D array to value provided.
*/
NTTW_DLL_SYM void initSigned_1D(long *data, const size_t size, const long value);

/**
 * \fn initUChar_1D(unsigned char *data, const size_t size, const unsigned char value)
 * \brief Initializes a unsigned character 1D array to value provided.
*/
NTTW_DLL_SYM void initUChar_1D(unsigned char *data, const size_t size, const unsigned char value);

/**
 * \fn init_1D_big(nttw_big_integer *data, const size_t size, const nttw_big_integer value)
 * \brief Initializes the 1D array to value provided.
*/
NTTW_DLL_SYM void init_1D_big(nttw_big_integer *data, const size_t size, const nttw_big_integer value);

/**
 * \fn free_array(void *pointer)
 * \brief Frees the Byte Aligned array allocated by array() or array_1D() correctly.
*/
NTTW_DLL_SYM void free_array(void *pointer);
//@}
#endif // ARRAY_H_INCLUDED
