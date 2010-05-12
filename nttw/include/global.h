/**
 * NTTW Library Aliases
 * \file global.h
 * \brief Aliases/Globals Header/Object for the NTTW C Library.
 *
 * This header provides all the relevant macros, alaises and constants used in the library.
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
 * \author Shekhar S. Chandra, 2009
*/
#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED

/**
* \mainpage NTTW C Library API
*
* \section intro_sec Introduction
*
* The various Number Theoretic Transform (NTT) Algorithms Optimised for C and performance. Array, Timing and Imaging modules are also provided.
*
* \section features Features
*
* Features an Array module, a micro-second Timing module an Imaging module and an Number Theoretic module for dyadic and prime lengths.
* There are 1D and 2D versions of the NTTs.
*/

#include <limits.h>

#define NTTW_VERSION       0x000001
#define NTTW_VERSION_STR   "1.0.0"

#if defined(WIN32)

    #if defined(NTTW_DLL)
        /**
        Used for DLL generation purposes (Windows Specific) Import/Export.
        Templates classes cannot be imported hence its own variable.
        The export command is
        */
        #if defined(NTTW_MAKEDLL)     // create a NTTW DLL library
            #define NTTW_DLL_SYM  __declspec(dllexport)
        #else                        // use a NTTW DLL library
            #define NTTW_DLL_SYM  __declspec(dllimport)
        #endif
    #endif // NTTW_DLL
#endif // WIN32

/**
    \def NTTW_DLL_SYM
    \brief DLL Function Symbol for Windows. It is empty for other OSes.
*/
#ifndef NTTW_DLL_SYM
    #define NTTW_DLL_SYM
#endif

/**
    \def CACHE_ALIGN
    \brief Align the integers to match SSE register sizes. UNUSED.
*/
#ifndef CACHE_ALIGN
    #define CACHE_ALIGN
#endif

/**
    \def TRUE
    \brief Define for True Boolean
*/
#ifndef TRUE
    #define TRUE 1
#endif
/**
    \def FALSE
    \brief Define for False Boolean
*/
#ifndef FALSE
    #define FALSE 0
#endif

/**
    \def NTTW_FORWARD
    \brief Define for Forward transform of NTTs
*/
#ifndef NTTW_FORWARD
    #define NTTW_FORWARD -1
#endif
/**
    \def NTTW_INVERSE
    \brief Define for Inverse transform of NTTs
*/
#ifndef NTTW_INVERSE
    #define NTTW_INVERSE 1
#endif

/**
    \def ALIGNOF
    \brief Macro for determining the alignment of a variable. Cross platform.
*/
#if defined WIN32
    #define ALIGNOF(type) __alignof(type)
#else
    #define ALIGNOF(type) __alignof__(type)
#endif

/**
    \def MACHINE_MAX
    \brief Define for the machine architecture type
*/
#define MACHINE_MAX LONG_MAX
/**
    \def ROUND
    \brief Macro for rounding an floating point number.
*/
#define ROUND(x) ( (x)>=0 ? (long)((x)+0.5) : (long)((x)-0.5) )
/**
    \def IS64BITMODE
    \brief Macro for determining if 64-bit mode. Returns Boolean.
*/
#define IS64BITMODE ( MACHINE_MAX>INT_MAX ? 1 : 0 )
/**
    \def SWAP
    \brief Macro for swapping the contents of two variables.
*/
#define SWAP(a, b) {a ^= b; b ^= a; a ^= b;}

/**
    \def BITS
    \brief Define for the number of bits of the local machine.
*/
#ifndef BITS
    #define BITS 8*sizeof(long) //Machine bits
#endif

/**
    \def NTTW_64
    \brief Define building a 64-bit version of the library.
*/
#if defined NTTW_64 //If 64-bit build required use max of 64-bit integers
    #ifndef NTTW_INT_MAX
        #define NTTW_INT_MAX UINT_MAX
    #endif
    #ifndef NTTW_INT_MIN
        #define NTTW_INT_MIN 0
    #endif

    #ifndef NTTW_BIGINT_MAX
        #define NTTW_BIGINT_MAX ULONG_MAX //ULLONG_MAX
    #endif
    #ifndef NTTW_BIGINT_MIN
        #define NTTW_BIGINT_MIN 0
    #endif

    #ifndef FORMAT_INPUT_STR
        #define FORMAT_INPUT_STR "%u"
    #endif
    #ifndef FORMAT_OUTPUT_STR
        #define FORMAT_OUTPUT_STR "%u "
    #endif
    #ifndef FORMAT_CSVOUTPUT_STR
        #define FORMAT_CSVOUTPUT_STR "%u,"
    #endif

    #ifndef FORMAT_BIG_OUTPUT_STR
        #define FORMAT_BIG_OUTPUT_STR "%llu "
    #endif
    #ifndef FORMAT_BIG_CSVOUTPUT_STR
        #define FORMAT_BIG_CSVOUTPUT_STR "%llu,"
    #endif

    typedef unsigned nttw_integer;
    typedef double nttw_real;
    typedef unsigned long long nttw_big_integer;
#else //Else 32-bit build required use max of 32-bit integers
    #ifndef NTTW_INT_MAX
        #define NTTW_INT_MAX USHRT_MAX
    #endif
    #ifndef NTTW_INT_MIN
        #define NTTW_INT_MIN 0
    #endif

    #ifndef NTTW_BIGINT_MAX
        #define NTTW_BIGINT_MAX UINT_MAX
    #endif
    #ifndef NTTW_BIGINT_MIN
        #define NTTW_BIGINT_MIN 0
    #endif

    #ifndef FORMAT_INPUT_STR
        #define FORMAT_INPUT_STR "%hu"
    #endif
    #ifndef FORMAT_OUTPUT_STR
        #define FORMAT_OUTPUT_STR "%hu "
    #endif
    #ifndef FORMAT_CSVOUTPUT_STR
        #define FORMAT_CSVOUTPUT_STR "%hu,"
    #endif

    #ifndef FORMAT_BIG_OUTPUT_STR
        #define FORMAT_BIG_OUTPUT_STR "%u "
    #endif
    #ifndef FORMAT_BIG_CSVOUTPUT_STR
        #define FORMAT_BIG_CSVOUTPUT_STR "%u,"
    #endif

    typedef unsigned short nttw_integer;
    typedef float nttw_real;
    typedef unsigned nttw_big_integer;
#endif

#endif // GLOBAL_H_INCLUDED
