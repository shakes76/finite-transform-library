/**
 * NTTW NTT 32-bit Module
 * \file number32.h
 * \brief NTT Header/Object for the NTTW C Library.
 *
 * This file defines the functions for the Generic 32-bit Modulus NTTs.
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
#ifndef MOD32_H_INCLUDED
#define MOD32_H_INCLUDED

#include <stdlib.h>

#include "global.h"

/**
// Modulus (For Reasonable N)
#define Modulus 2113929217ULL //63*2^25+1
#define PrimitiveRoot 5
#define Modulus 4293918721ULL //2^32-2^20+1
#define PrimitiveRoot 19
// Modulus (For Reasonable N)
#define Modulus 1811939329U //
#define PrimitiveRoot 13
// Modulus (For Reasonable N)
#define MODULUS 2147473409U //
#define PROOT 3
#define Modulus 18446742974197923841ULL //2^64-2^40+1
#define PrimitiveRoot 19
#define MODULUS 65537U //upto 2^16 (p=2^16+1)
#define PROOT 3
#define MODULUS 40961 //upto 2^13 (p=5*2^13+1)
#define PROOT 3
*/

#if defined NTTW_64
    #define MODULUS 2147473409U //
    #define PROOT 3
#else
    #define MODULUS 40961 //upto 2^13 (p=5*2^13+1)
    #define PROOT 3
#endif

#define MODSUB(a, b, m) ( ((a) < (b)) ? ((a)-(b)) + m : ((a)-(b)) )
#define MODADD(a ,b, m) ( ((a+b) >= m) ? (a+b) - m : (a+b) )
#define MODADD2(a ,b, m) ( ((a+b) < a || (a+b) >= m) ? (a+b) - m : (a+b) )
#define MODMUL(a, b, m) ( ((a*b) >= (m)) ? (unsigned long long)((a*b) - (m)) : (unsigned long long)(a*b) )

/**
 * \defgroup Number_Theoretic_Transforms Number Theoretic Transforms (NTTs) and Helpers
 * \brief The Number Theoretic Transform functions.
 *
 * The main function for NTTs is fntt() and it requires that the sequence be dyadic in length.
 * Two dimensional transform is also assumed to be dyadic and can be done via fntt_2D().
 */
//@{
/**
* \fn rearrange(nttw_integer *data, const size_t n)
* \brief Bit-reversal of Data of size n. n should be dyadic.
*/
NTTW_DLL_SYM void rearrange(nttw_integer *data, const size_t n);
NTTW_DLL_SYM void rearrange_big(nttw_big_integer *data, const size_t n);

/**
 * \fn bigshr(nttw_integer *d, nttw_integer *s, const size_t n)
 * \brief Shifts s right one bit to d, returns carry bit.
* \author Mikko Tommila et al., APFloat Library 2005.
*/
NTTW_DLL_SYM int bigshr(nttw_integer *d, nttw_integer *s, const size_t n);
NTTW_DLL_SYM int bigshr_big(nttw_big_integer *d, nttw_big_integer *s, const size_t n);

/**
 * \fn pow_mod(nttw_integer base, nttw_integer exp)
 * \brief Computers the power of the base modulo the defined MODULUS.
 * \author Mikko Tommila et al., APFloat Library 2005.
*/
NTTW_DLL_SYM nttw_integer pow_mod(nttw_integer base, nttw_integer exp);

/**
 * \fn pow_mod2(nttw_integer base, nttw_integer exp, nttw_integer modulus)
 * \brief Computers the power of the base modulo the given modulus.
*/
NTTW_DLL_SYM nttw_integer pow_mod2(nttw_integer base, nttw_integer exp, nttw_integer modulus);

/**
 * \fn pow_mod2_long(nttw_integer base, nttw_integer exp, nttw_integer modulus)
 * \brief Computers the power of the base modulo the given modulus. Force unsigned long long's in calculations.
*/
NTTW_DLL_SYM nttw_integer pow_mod2_long(nttw_integer base, nttw_integer exp, nttw_integer modulus);

/**
 * \fn pow_mod_long(nttw_integer base, nttw_integer exp)
 * \brief Computers the power of the base modulo the given modulus. Big Integer version.
*/
NTTW_DLL_SYM nttw_integer pow_mod_long(nttw_integer base, nttw_integer exp);

/**
 * \fn pow_mod_long(nttw_integer base, nttw_integer exp)
 * \brief Computers the power of the base modulo the given modulus. Big Integer version.
*/
NTTW_DLL_SYM nttw_integer pow_mod_direct(nttw_integer base, nttw_integer exp, nttw_integer modulus);

/**
 * \fn ntt(nttw_integer *data, const size_t nn, const nttw_integer pr, const int isign)
 * \brief Computes the 1D Number Theoretic Transform (NTT) using a direct algorithm (slow).
 *
 * Computes the 1D Number Theoretic Transform (NTT). The result is NOT normalised within the function.
 * The transform is done out-of-place, allocating the output, you must free. Use the ntt_norm() function for normalisation.
*/
NTTW_DLL_SYM nttw_integer* ntt(nttw_integer *data, const size_t nn, const nttw_integer pr, const int isign);

/**
 * \fn fntt(nttw_integer *data, const size_t nn, const nttw_integer pr, const int isign)
 * \brief Computes the 1D Fast Number Theoretic Transform (FNTT) using the Cooley-Tukey algorithm.
 *
 * Computes the 1D Fast Number Theoretic Transform (FNTT). The result is NOT normalised within the function.
 * The transform is done inplace, destroying the input. Use the ntt_norm() function for normalisation.
*/
NTTW_DLL_SYM void fntt(nttw_integer *data, const size_t nn, const nttw_integer pr, const int isign);

/**
 * \fn fntt_2D(nttw_integer *data, nttw_integer *result, const size_t nn, const int isign)
 * \brief Computes the 2D Fast Number Theoretic Transform (FNTT) using the Cooley-Tukey algorithm.
 *
 * Computes the 2D Fast Number Theoretic Transform (FNTT). The result is normalised within the function.
 * The transform is done inplace, destroying the input but copied to result. The data is assumed to be a square array.
*/
NTTW_DLL_SYM void fntt_2D(nttw_integer *data, nttw_integer *result, const size_t nn, const int isign);

/**
 * \fn fntt_prime(const nttw_integer *inData, nttw_integer *outData, const size_t p, const nttw_integer root, const nttw_integer primeDash, const nttw_integer proot, int isign, const int norm)
 * \brief Computes the Prime Length Fast Number Theoretic Transform (FNTT) using Rader's algorithm.
 *
 * Computes the Prime Length Fast Number Theoretic Transform (FNTT). The result is normalised given the Boolean.
 * The transform is copied to result (outData).
*/
NTTW_DLL_SYM void fntt_prime(const nttw_integer *inData, nttw_integer *outData, const size_t p, const nttw_integer root, const nttw_integer primeDash, const nttw_integer proot, int isign, const int norm);

/**
 * \fn fntt_2D_prime(nttw_integer *data, nttw_integer *result, const size_t n, const nttw_integer root, const nttw_integer primeDash, const nttw_integer proot, const int isign, const int norm)
 * \brief Computes the 2D Prime Length Fast Number Theoretic Transform (FNTT) using Rader's algorithm.
 *
 * Computes the 2D Prime Length Fast Number Theoretic Transform (FNTT). The result is normalised given the Boolean.
 * The transform is copied to result (outData). The data is assumed to be a square array.
*/
NTTW_DLL_SYM void fntt_2D_prime(nttw_integer *data, nttw_integer *result, const size_t n, const nttw_integer root, const nttw_integer primeDash, const nttw_integer proot, const int isign, const int norm);

/**
 * \fn padData_Rader(nttw_integer *data, const size_t p, const size_t newSize, const nttw_integer primeDash)
 * \brief Pads the data to the nearest highly composite length as needed by Rader's algorithm.
*/
NTTW_DLL_SYM nttw_integer* padData_Rader(nttw_integer *data, const size_t p, const size_t newSize, const nttw_integer primeDash);

/**
 * \fn padTransformMatrix_Rader(nttw_integer *transData, const size_t p, const size_t newSize)
 * \brief Pads the transform matrix to the nearest highly composite length as needed by Rader's algorithm.
*/
NTTW_DLL_SYM nttw_integer* padTransformMatrix_Rader(nttw_integer *transData, const size_t p, const size_t newSize);

/**
 * \fn extractResult_Rader(nttw_integer *paddedData, nttw_integer *data, const size_t p)
 * \brief Extracts the data for the result from the padded result.
*/
NTTW_DLL_SYM void extractResult_Rader(nttw_integer *paddedData, nttw_integer *data, const size_t p);

/**
 * \fn ntt_norm(nttw_integer *data, const size_t n)
 * \brief Normalises data of length n in the number field of MODULUS.
*/
NTTW_DLL_SYM void ntt_norm(nttw_integer *data, const size_t n);
//@}
#endif //MOD32_H
