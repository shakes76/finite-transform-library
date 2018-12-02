/**
 * NTTW Prime/Number Theory Module
 * \file prime.h
 * \brief Prime/Number Theory Header/Object for the NTTW C Library.
 *
 * This file defines the functions for Prime numbers and other Number Theory related methods.
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
#ifndef PRIME_H_INCLUDED
#define PRIME_H_INCLUDED

#include <stdlib.h>

#include "global.h"

/**
 * \defgroup NT Number Theory Routines
 * \brief Routines for basic number theory used for NTTs.
 */
//@{
/**
    \fn loadPrimes_Small(nttw_integer **data, size_t *size)
    \brief Loads the first 100 primes in data. Memory is allocated and data is passed by reference.
*/
NTTW_DLL_SYM void loadPrimes_Small(nttw_integer **data, size_t *size);

/**
 * \fn isPrimeHighlyComposite(const size_t p, size_t *power)
 * \brief Is the p-1 a dyadic number. So the Fermat prime 257 will return true.
*/
NTTW_DLL_SYM size_t isPrimeHighlyComposite(const size_t p, size_t *power);

/**
 * \fn newHighlyCompositeSize(const size_t p)
 * \brief Return the closest dyadic number (16,32,etc.) greater than 2*p-4. Used for Rader's FFT algorithm.
*/
NTTW_DLL_SYM size_t newHighlyCompositeSize(const size_t p);

/**
    \fn findClosestPrime(const nttw_integer *primeList, const size_t size, const nttw_integer number, size_t *nthPrime, int *isPrime)
    \brief Finds the closest prime which is larger than the number in the primes list provided.

    \param nthPrime - Stores the nth prime (index) which is closest. e.g. it may be the 21st prime which is closest.
    \param isPrime - Will be set to TRUE if number is a prime.
*/
NTTW_DLL_SYM nttw_integer findClosestPrime(const nttw_integer *primeList, const size_t size, const nttw_integer number, size_t *nthPrime, int *isPrime);

/**
 * \fn findAlternatePrime(const nttw_integer *primeList, const size_t size, const nttw_integer number)
 * \brief Finds an alternate prime p' of the form p' = k*p+1, where p is prime.
*/
NTTW_DLL_SYM nttw_integer findAlternatePrime(const nttw_integer *primeList, const size_t size, const nttw_integer number);

/**
    \fn findFactors(const nttw_integer *primeList, const size_t size, const nttw_integer number, nttw_integer **factors, size_t *noOfFactors)
    \brief Finds the unique prime factors of the number using the primes list provided.

    \param factors - Is a 1D array passed by reference which will contain the found factors.
*/
NTTW_DLL_SYM void findFactors(const nttw_integer *primeList, const size_t size, const nttw_integer number, nttw_integer **factors, size_t *noOfFactors);

/**
 * \fn findFirstPrimitiveRoot(const nttw_integer prime)
 * \brief Finds the smallest primitive root of the prime.
*/
NTTW_DLL_SYM nttw_integer findFirstPrimitiveRoot(const nttw_integer prime);
/**
 * \fn findFirstPrimitiveRoot_Full(const nttw_integer *primes, const size_t noOfPrimes, const nttw_integer modulus)
 * \brief Finds the smallest primitive root of the prime using the list of primes provided.
*/
NTTW_DLL_SYM nttw_integer findFirstPrimitiveRoot_Full(const nttw_integer *primes, const size_t noOfPrimes, const nttw_integer modulus);

/**
 * \fn findPrimeLengthParameters(const nttw_integer prime, nttw_integer *root, nttw_integer *modulus, nttw_integer *proot)
 * \brief Finds the parameters need for prime length NTTs.
*/
NTTW_DLL_SYM void findPrimeLengthParameters(const nttw_integer prime, nttw_integer *root, nttw_integer *modulus, nttw_integer *proot);

/**
* \fn euclidean(nttw_big_integer a, nttw_big_integer b, nttw_big_integer *x, nttw_big_integer *y)
* \brief The Extended Euclidean Algorithm for solving congruences. C Version.
* Returns d = gcd(a,b) and x, y so that ax + by = d.
* When a and b are coprime, then x is the modular multiplicative inverse of a modulo b.
*/
NTTW_DLL_SYM nttw_big_integer euclidean(nttw_big_integer a, nttw_big_integer b, nttw_big_integer *x, nttw_big_integer *y);
/**
* \fn euclidean_long(long a, long b, long *x, long *y)
* \brief The Extended Euclidean Algorithm for solving congruences signed version. C Version.
* Returns d = gcd(a,b) and x, y so that ax + by = d.
* When a and b are coprime, then x is the modular multiplicative inverse of a modulo b.
*/
NTTW_DLL_SYM long euclidean_long(long a, long b, long *x, long *y);

/**
 * \fn gcd(nttw_big_integer u, nttw_big_integer v)
 * \brief Finds the GCD using Stein's Binary Algorithm. Code from Wikipedia.Org
*/
NTTW_DLL_SYM nttw_big_integer gcd(nttw_big_integer u, nttw_big_integer v);

/**
 * \fn modpow(nttw_integer base, nttw_integer exponent, nttw_integer modulus)
 * \brief Computes the modular exponentiation of the base to the exponent modulo modulus. Code from Wikipedia.Org
 *
 * This code is not as fast as the powmod() function in number32.h file. This function is provided to keep the Prime Module
 * self contained.
*/
NTTW_DLL_SYM nttw_integer modpow(nttw_integer base, nttw_integer exponent, nttw_integer modulus);

/**
 * \fn minverse(const nttw_integer num, const nttw_integer mod)
 * \brief Computes the Multiplicative inverse using the Extended Euclidean Algorithm.
*/
NTTW_DLL_SYM nttw_big_integer minverse(const nttw_integer num, const nttw_integer mod);
/**
 * \fn minverse_big(const nttw_big_integer num, const nttw_big_integer mod)
 * \brief Computes the Multiplicative inverse using the Extended Euclidean Algorithm. Big Integer Version
*/
NTTW_DLL_SYM nttw_big_integer minverse_big(const nttw_big_integer num, const nttw_big_integer mod);
/**
 * \fn minverse_signed(const long num, const long mod)
 * \brief Computes the Multiplicative inverse using the Extended Euclidean Algorithm. Signed Version
*/
NTTW_DLL_SYM long minverse_signed(const long num, const long mod);
//@}
#endif // PRIME_H_INCLUDED
