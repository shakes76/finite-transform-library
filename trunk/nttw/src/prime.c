/**
 * NTTW Prime/Number Theory Module
 * \file prime.c
 * \brief Prime/Number Theory Source/Object for the NTTW C Library.
 *
 * This file implements the functions for Prime numbers and other Number Theory related methods.
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
#include <stdio.h>

#include "array.h" //Order important for gcc 4.3.3
#include "prime.h"

void loadPrimes_Small(nttw_integer **data, size_t *size)
{
    size_t j;
    const size_t noOfPrimes = 100;
    short primes[] =
    {   2,      3,      5,      7,      11,     13,     17,     19,     23,     29,
        31,     37,     41,     43,     47,     53,     59,     61,     67,     71,
        73,     79,     83,     89,     97,    101,    103,    107,    109,    113,
        127,   131,    137,    139,    149,    151,    157,    163,    167,    173,
        179,   181,    191,    193,    197,    199,    211,    223,    227,    229,
        233,   239,    241,    251,    257,    263,    269,    271,    277,    281,
        283,   293,    307,    311,    313,    317,    331,    337,    347,    349,
        353,   359,    367,    373,    379,    383,    389,    397,    401,    409,
        419,   421,    431,    433,    439,    443,    449,    457,    461,    463,
        467,   479,    487,    491,    499,    503,    509,    521,    523,    541};

    *data = array_1D(noOfPrimes);
    for(j = 0; j < noOfPrimes; j ++)
        (*data)[j] = primes[j];
    *size = noOfPrimes;
}

size_t isPrimeHighlyComposite(const size_t p, size_t *power)
{
    ///Check is p-1 is dyadic. This is normally a Fermat prime.
    size_t count = 0, comp = p-1;

    while(comp % 2 == 0) //Is divisible by 2?
    {
        comp /= 2; //Compiler should optimize this to left shift here
        count ++;
    }
    *power = count;

    if(comp != 1)
        return FALSE;
    else
        return TRUE;
}

size_t newHighlyCompositeSize(const size_t p)
{
    size_t comp = 2;

    while(comp < p || comp <= 2*p-4) //Ensure largest dyadic size for Rader
        comp *= 2; //Compiler should optimize this to right shift here

    return comp;
}

nttw_integer findClosestPrime(const nttw_integer *primeList, const size_t size, const nttw_integer number, size_t *nthPrime, int *isPrime)
{
    size_t j;
    nttw_integer closePrime = 0;

    *isPrime = FALSE;
    for(j = 0; j < size; j ++)
    {
        if(number == primeList[j])
        {
            closePrime = (nttw_integer)primeList[j];
            *nthPrime = j;
            *isPrime = TRUE;
            break;
        }
        else if(number > primeList[j])
            continue;
        else
        {
            closePrime = (nttw_integer)primeList[j];
            *nthPrime = j;
            break;
        }
    }

    return closePrime;
}

nttw_integer findAlternatePrime(const nttw_integer *primeList, const size_t size, const nttw_integer number)
{
    size_t j, k;
    int found = FALSE;
    nttw_integer tmpPrime = 2, altPrime = 0;

    for(j = 1; j < size; j ++)
    {
        tmpPrime = (nttw_integer) j*number+1;
        for(k = 0; k < size; k ++)
        {
            if(tmpPrime == primeList[k])
            {
                altPrime = tmpPrime;
                found = TRUE;
                break;
            }
            else if(tmpPrime < primeList[k])
                break; //Not found
        }
        if(found)
            break;
    }

    return altPrime;
}

void findFactors(const nttw_integer *primeList, const size_t size, const nttw_integer number, nttw_integer **factors, size_t *noOfFactors)
{
    size_t j, approxNoOfFactors = 1;
    int isPrime = FALSE;
    nttw_integer closePrime, *tmpFactors;

    closePrime = findClosestPrime(primeList,size,number,&approxNoOfFactors,&isPrime);

    if(isPrime)
    {
        *factors = array_1D(1);
        (*factors)[0] = number;
        *noOfFactors = 1;
    }
    else
    {
        tmpFactors = array_1D(approxNoOfFactors);

        *noOfFactors = 0;
        for(j = 0; j < size; j ++)
        {
            if(number % primeList[j] == 0)
            {
                tmpFactors[*noOfFactors] = (nttw_integer)primeList[j];
                (*noOfFactors) ++;
            }

            if(primeList[j] > number/2) //Only need to find factors upto half of number
                break;
        }

        ///Copy into proper sized factors array
        *factors = array_1D(*noOfFactors);

        for(j = 0; j < *noOfFactors; j ++)
            (*factors)[j] = tmpFactors[j];

        free_array(tmpFactors);
    }
}

nttw_integer findFirstPrimitiveRoot(const nttw_integer modulus)
{
    size_t j, k, noOfFactors, noOfPrimes, count = 0;
    nttw_integer root = 0, power = modulus-1;
    nttw_integer *primes, *factors;

    loadPrimes_Small(&primes,&noOfPrimes);

    findFactors(primes,noOfPrimes,power,&factors,&noOfFactors);

    for(k = 0; k < noOfPrimes; k ++)
    {
        count = 0;

        if(primes[k] >= modulus) //Only need todo numbers < modulus
            break;

        for(j = 0; j < noOfFactors; j ++)
        {
            if( modpow(primes[k],power/factors[j],modulus) == 1 ) //Cant be root
                break;
            else
                count ++;
        }

        if(count == noOfFactors) //Found root
        {
            root = (nttw_integer)primes[k];
            break;
        }
    }

    free_array(primes);
    free_array(factors);
    return root;
}

nttw_integer findFirstPrimitiveRoot_Full(const nttw_integer *primes, const size_t noOfPrimes, const nttw_integer modulus)
{
    size_t j, k, noOfFactors, count = 0;
    nttw_integer root = 0, power = modulus-1;
    nttw_integer *factors;

    if(noOfPrimes < 1)
        return 0;

    findFactors(primes,noOfPrimes,power,&factors,&noOfFactors);

    for(k = 0; k < noOfPrimes; k ++)
    {
        count = 0;

        if(primes[k] >= modulus) //Only need todo numbers < modulus
            break;

        for(j = 0; j < noOfFactors; j ++)
        {
            if( modpow(primes[k],power/factors[j],modulus) == 1 ) //Cant be root
                break;
            else
                count ++;
        }

        if(count == noOfFactors) //Found root
        {
            root = (nttw_integer)primes[k];
            break;
        }
    }

    free_array(factors);
    return root;
}

nttw_big_integer euclidean(nttw_big_integer a, nttw_big_integer b, nttw_big_integer *x, nttw_big_integer *y)
{
    nttw_big_integer q,r,x_1,x_2,y_1,y_2;

    if( b == 0 ){
        *x = 1;  *y = 0;
        return a;
    }
    x_1 = 0;   x_2 = 1;
    y_1 = 1;   y_2 = 0;
    while( b != 0 ){
        q = a/b;
        r = a - q*b;
        *x = x_2 - q*x_1;
        *y = y_2 - q*y_1;
        a = b;
        b = r;
        x_2 = x_1;  x_1 = *x;
        y_2 = y_1;  y_1 = *y;
    }
    *x = x_2;  *y = y_2;

    return a;
}

long euclidean_long(long a, long b, long *x, long *y)
{
    long q,r,x_1,x_2,y_1,y_2;

    if( b == 0 ){
        *x = 1;  *y = 0;
        return a;
    }
    x_1 = 0;   x_2 = 1;
    y_1 = 1;   y_2 = 0;
    while( b != 0 ){
        q = a/b;
        r = a - q*b;
        *x = x_2 - q*x_1;
        *y = y_2 - q*y_1;
        a = b;
        b = r;
        x_2 = x_1;  x_1 = *x;
        y_2 = y_1;  y_1 = *y;
    }
    *x = x_2;  *y = y_2;

    return a;
}

nttw_integer modpow(nttw_integer base, nttw_integer exponent, nttw_integer modulus)
{
    nttw_integer result = 1;

    while (exponent > 0)
    {
        if ((exponent & 1) == 1)
        {
            // multiply in this bit's contribution while using modulus to keep result small
            result = (result * base) % modulus;
        }
        // move to the next bit of the exponent, square (and mod) the base accordingly
        exponent >>= 1;
        base = (base * base) % modulus;
    }

    return result;
}

nttw_big_integer minverse(const nttw_integer num, const nttw_integer mod)
{
    nttw_big_integer inv = 0, d, y;

    d = euclidean((nttw_big_integer)num,mod,&inv,&y); //Multi Inv of p-1
    inv = (inv + mod)%mod; //Ensure x is positive

    return inv;
}

nttw_big_integer minverse_big(const nttw_big_integer num, const nttw_big_integer mod)
{
    nttw_big_integer inv = 0, d, y;

    d = euclidean(num,mod,&inv,&y); //Multi Inv of p-1
    inv = (inv + mod)%mod; //Ensure x is positive

    return inv;
}

long minverse_signed(const long num, const long mod)
{
    long inv = 0, d, y;

    d = euclidean_long((nttw_big_integer)num,mod,&inv,&y); //Multi Inv of p-1

    return inv;
}
