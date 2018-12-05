/**
 * NTTW NTT 32-bit Module
 * \file number32.c
 * \brief NTT Source for the NTTW C Library.
 *
 * This file implements the functions for the Generic 32-bit Modulus NTTs.
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
#include <stdio.h>

#include "nttw/number.h"
#include "nttw/array.h"
#include "nttw/prime.h"

void rearrange(nttw_integer *data, const size_t n)
{
    size_t target = 0, position = 0, mask = n;

    ///For all of input signal
    for (position = 0; position < n; ++position)
    {
        ///Ignore swapped entries
        if (target > position)
        {
            ///Swap
            SWAP(data[target],data[position]);
        }

        ///Bit mask
        mask = n;
        ///While bit is set
        while (target & (mask >>= 1))
            ///Drop bit
            target &= ~mask;
        ///The current bit is 0 - set it
        target |= mask;
    }
}

void rearrange_big(nttw_big_integer *data, const size_t n)
{
    size_t target = 0, position = 0, mask = n;

    ///For all of input signal
    for (position = 0; position < n; ++position)
    {
        ///Ignore swapped entries
        if (target > position)
        {
            ///Swap
            SWAP(data[target],data[position]);
        }

        ///Bit mask
        mask = n;
        ///While bit is set
        while (target & (mask >>= 1))
            ///Drop bit
            target &= ~mask;
        ///The current bit is 0 - set it
        target |= mask;
    }
}

// Shifts s right one bit to d, returns carry bit
int bigshr (nttw_integer *d, nttw_integer *s, const size_t n)
{
    size_t t;
    nttw_integer tmp;
    int b = 0;

    if (!n) return 0;

    d += n;
    s += n;

    for (t = 0; t < n; t++)
    {
        d--;
        s--;

        tmp = (*s >> 1) + (b ? 0x80000000 : 0);
        b = *s & 1;
        *d = tmp;                               // Works also if d = s
    }

    return b;
}

int bigshr_big(nttw_big_integer *d, nttw_big_integer *s, const size_t n)
{
    size_t t;
    nttw_big_integer tmp;
    int b = 0;

    if (!n) return 0;

    d += n;
    s += n;

    for (t = 0; t < n; t++)
    {
        d--;
        s--;

        tmp = (*s >> 1) + (b ? 0x80000000 : 0);
        b = *s & 1;
        *d = tmp;                               // Works also if d = s
    }

    return b;
}

// Uses standard O(log exp) scheme
// For extreme portability this quite slow, so avoid whenever possible
nttw_integer pow_mod(nttw_integer base, nttw_integer exp)
{
    nttw_integer r;

    if (exp == 0) return 1;

    while (!bigshr (&exp, &exp, 1))
        base = ((nttw_big_integer)base * base)%MODULUS;

    r = base;

    while (exp > 0)
    {
        base = ((nttw_big_integer)base * base)%MODULUS;
        if (bigshr (&exp, &exp, 1))
            r = ((nttw_big_integer)r * base)%MODULUS;
    }

    return r;
}

nttw_integer pow_mod2(nttw_integer base, nttw_integer exp, nttw_integer modulus)
{
    nttw_integer r;

    if (exp == 0) return 1;

    while (!bigshr (&exp, &exp, 1))
        base = ((nttw_big_integer)base * base)%modulus;

    r = base;

    while (exp > 0)
    {
        base = ((nttw_big_integer)base * base)%modulus;
        if (bigshr (&exp, &exp, 1))
            r = ((nttw_big_integer)r * base)%modulus;
    }

    return r;
}

nttw_integer pow_mod2_long(nttw_integer base, nttw_integer exp, nttw_integer modulus)
{
    nttw_integer r;

    if (exp == 0) return 1;

    while (!bigshr (&exp, &exp, 1))
        base = ((unsigned long long)base * base)%modulus;

    r = base;

    while (exp > 0)
    {
        base = ((unsigned long long)base * base)%modulus;
        if (bigshr (&exp, &exp, 1))
            r = ((unsigned long long)r * base)%modulus;
    }

    return r;
}

nttw_integer pow_mod_long(nttw_integer base, nttw_integer exp)
{
    nttw_integer r;

    if (exp == 0) return 1;

    while (!bigshr (&exp, &exp, 1))
        base = ((unsigned long long)base * base)%MODULUS;

    r = base;

    while (exp > 0)
    {
        base = ((unsigned long long)base * base)%MODULUS;

        if (bigshr (&exp, &exp, 1))
            r = ((unsigned long long)r * base)%MODULUS;
    }

    return r;
}

nttw_integer pow_mod_direct(nttw_integer base, nttw_integer exp, nttw_integer modulus)
{
  nttw_integer r = 1;
  size_t e = 0;

  if (exp == 0) return r;

  for(e = 0; e < exp; e ++)
      r = ( r * (nttw_big_integer)base )%modulus;

  return r;
}

///Direct NTT
nttw_integer* ntt(nttw_integer *data, const size_t nn, const nttw_integer pr, const int isign)
{
    size_t m, i, j;
    nttw_integer w, wtemp, *result;

    if (isign > 0)
//        w = pow_mod(pr, MODULUS - 1 - (MODULUS - 1) / (nttw_integer)nn);
        w = pow_mod_direct(pr, MODULUS - 1 - (MODULUS - 1) / (nttw_integer)nn, MODULUS);
    else
//        w = pow_mod(pr, (MODULUS - 1) / (nttw_integer)nn);
        w = pow_mod_direct(pr, (MODULUS - 1) / (nttw_integer)nn, MODULUS);

    result = array_1D(nn);
    init_1D(result, nn, 0);
    for(j = 0; j < nn; j ++)
    {
        for(i = 0; i < nn; i ++)
        {
            //double width for multiply, keep type consistent
//            wtemp = ( data[i]*(nttw_big_integer)pow_mod(w, (nttw_integer)(j*i)) )%MODULUS; //fast
            wtemp = ( data[i]*(nttw_big_integer)pow_mod_direct(w, (nttw_integer)(j*i), MODULUS) )%MODULUS; //slow
            result[j] = MODADD(result[j], wtemp, MODULUS);
        }
    }

    return result;
}

///Cooley-Tukey NTT
void fntt(nttw_integer *data, const size_t nn, const nttw_integer pr, const int isign)
{
    size_t m, i, j, istep, mmax;
    nttw_integer w, wt, wr, wtemp;

    if (isign > 0)
        w = pow_mod(pr, MODULUS - 1 - (MODULUS - 1) / (nttw_integer)nn);
    else
        w = pow_mod(pr, (MODULUS - 1) / (nttw_integer)nn);

    rearrange(data, nn);

    mmax = 1;
    while (nn > mmax)
    {
        istep = mmax << 1;
        wr = wt = pow_mod(w, (nttw_integer)nn / (nttw_integer)istep);

        // Optimize first step when wr = 1
        for (i = 0; i < nn; i += istep)
        {
            j = i + mmax;
            wtemp = data[j];
            data[j] = MODSUB(data[i], wtemp, MODULUS);
            data[i] = MODADD(data[i], wtemp, MODULUS);
        }

        for (m = 1; m < mmax; m++)
        {
            for (i = m; i < nn; i += istep)
            {
                j = i + mmax;
                wtemp = (wr * (nttw_big_integer)data[j])%MODULUS;
                data[j] = MODSUB(data[i], wtemp, MODULUS);
                data[i] = MODADD(data[i], wtemp, MODULUS);
            }
            wr = (wr * (nttw_big_integer)wt)%MODULUS;
        }
        mmax = istep;
    }
}

///Cooley-Tukey
///data is corrupted in the process (for performance)
void fntt_2D(nttw_integer *data, nttw_integer *result, const size_t nn, const int isign)
{
    ///Pointers used where-ever possible, norm when copying
    size_t j, k;
    nttw_integer *ptrResult;
    nttw_big_integer inv, d, y;

    d = euclidean((nttw_big_integer)nn,MODULUS,&inv,&y); //Multi Inv of p-1
    inv = (inv + MODULUS)%MODULUS; //Ensure x is positive

    ///NT Rows
    for (j = 0; j < nn; j ++)
    {
        ptrResult = &data[j*nn];

        fntt(ptrResult, nn, PROOT, isign);
    }

    ptrResult = (nttw_integer *)malloc(nn*sizeof(nttw_integer));

    ///NT Columns
    for (j = 0; j < nn; j++)
    {
        for (k = 0; k < nn; k ++)
            ptrResult[k] = (data[k*nn+j] * inv)%MODULUS; ///Stops modulo overrun, div by N early

        fntt(ptrResult, nn, PROOT, isign);

        for (k = 0; k < nn; k ++) ///Inverse so Copy and Norm
            result[k*nn+j] = ptrResult[k];
    }

    free(ptrResult);
}

void fntt_prime(const nttw_integer *inData, nttw_integer *outData, const size_t p, const nttw_integer root, const nttw_integer primeDash, const nttw_integer proot, int isign, const int norm)
{
    nttw_integer sum, w, *orderedData, *orderedDataPadded, *W, *WPadded;
    size_t j, modulus = p, convolveSize = p-1, dyad;
    unsigned notHighlyComposite = FALSE;
    nttw_big_integer d, x1, x2, y;
    size_t index = 1, *indices;

    if (isign > 0) //Inverse Transform
        w = pow_mod2(proot, primeDash-1 - (primeDash-1)/(nttw_integer)p, primeDash);
    else
        w = pow_mod2(proot, (primeDash-1)/(nttw_integer)p, primeDash);

    ///Indices
    indices = (size_t *)malloc(sizeof(size_t)*convolveSize);
    indices[0] = index;
    for (j = 0; j < convolveSize-1; j ++) //1 to N-1
    {
        index *= root;
        index %= modulus;
        indices[convolveSize-1-j] = index;
    } //Works

    ///Reorder data and W matrix
    orderedData = array_1D(convolveSize);
    W = array_1D(convolveSize);
    index = 1;
    for (j = 0; j < convolveSize; j ++)
    {
        orderedData[j] = inData[ indices[j] ]%primeDash;
        W[j] = pow_mod2_long(w, (nttw_integer)index, primeDash); ///\todo Check casting here
        index *= root;
        index %= modulus;
    }

    ///Ensure convolution will be highly composite
    ///If not then pad sequences
    if ( !isPrimeHighlyComposite(p,&dyad) )
    {
        notHighlyComposite = TRUE;
        dyad = newHighlyCompositeSize(p);
        orderedDataPadded = padData_Rader(orderedData,p,dyad,primeDash);
        WPadded = padTransformMatrix_Rader(W,p,dyad);
        convolveSize = dyad;
    }
    else
    {
        orderedDataPadded = orderedData;
        WPadded = W;
    }

    ///Convolution
    ///Inv's
    d = euclidean((nttw_big_integer)convolveSize,MODULUS,&x1,&y); //Multi Inv of p-1
    x1 = (x1 + MODULUS)%MODULUS; //Ensure x is positive

    ///Data Eigenvalues
    isign *= -1;
    fntt(orderedDataPadded,convolveSize,PROOT,isign);

    ///W Eigenvalues
    fntt(WPadded,convolveSize,PROOT,isign);

    ///Convolve
    for (j = 0; j < convolveSize; j ++)
    {
        orderedDataPadded[j] = (orderedDataPadded[j] * (nttw_big_integer)x1)%MODULUS; //Divide by p-1
        orderedDataPadded[j] = (orderedDataPadded[j] * (nttw_big_integer)WPadded[j])%MODULUS;
    }

    ///Convolution from Convolution Eigenvalues
    isign *= -1;
    fntt(orderedDataPadded,convolveSize,PROOT,isign);

    ///Norm NTTs and add zeroth element
    if (isign > 0 && norm)
    {
        d = euclidean((nttw_big_integer)p,primeDash,&x2,&y); //Multi Inv of p
        x2 = (x2 + primeDash)%primeDash; //Ensure x is positive
    }
    convolveSize = p-1;
    for (j = 0; j < convolveSize; j ++)
    {
        orderedDataPadded[j] += inData[0]; //Plus zeroth element
        orderedDataPadded[j] %= primeDash;
    }

    ///Compute DC
    sum = 0;
    for (j = 0; j < p; j ++)
    {
        sum += inData[j];
        sum %= primeDash;
    }

    ///Reorder Result
    if (isign > 0 && norm) //Inverse
        sum = (sum * (nttw_big_integer)x2)%primeDash; //DC Divide by p
    outData[0] = sum; //DC
    if (notHighlyComposite) ///Copy padded result
        extractResult_Rader(orderedDataPadded,orderedData,p);

    for (j = 0; j < convolveSize; j ++)
    {
        index = (convolveSize-j)%convolveSize;
        index = indices[index];
        if (isign > 0 && norm) ///Norm
            orderedData[j] = (orderedData[j] * (nttw_big_integer)x2)%primeDash;
        outData[index] = orderedData[j];
    }

    free_array(orderedDataPadded);
    free_array(WPadded);
    if (notHighlyComposite)
    {
        free_array(orderedData);
        free_array(W);
    }
    free(indices);
}

void fntt_2D_prime(nttw_integer *data, nttw_integer *result, const size_t n, const nttw_integer root, const nttw_integer primeDash, const nttw_integer proot, const int isign, const int norm)
{
    size_t j, k;
    nttw_integer *ptrData, *ptrResult;
    const int normTransform = FALSE;
    nttw_integer inv = 1;

    if (isign > 0 && norm)
        inv = pow_mod2((nttw_integer)n*(nttw_integer)n,primeDash-2,primeDash);

    ///NT Rows
    for (j = 0; j < n; j ++)
    {
        ptrData = &data[j*n];
        ptrResult = &result[j*n];

        fntt_prime(ptrData, ptrResult, n, root, primeDash, proot, isign, normTransform);
    }

    ptrData = array_1D(n);
    ptrResult = array_1D(n);

    ///NT Columns
    for (j = 0; j < n; j++)
    {
        for (k = 0; k < n; k ++)
            ptrData[k] = result[k*n+j];

        fntt_prime(ptrData, ptrResult, n, root, primeDash, proot, isign, normTransform);

        if (isign > 0 && norm) //Inverse
            for (k = 0; k < n; k ++) ///Inverse so Copy and Norm
                result[k*n+j] = (ptrResult[k] * (nttw_big_integer)inv)%primeDash; //Div by nn
        else
            for (k = 0; k < n; k ++) ///Copy
                result[k*n+j] = ptrResult[k];
    }

    free_array(ptrData);
    free_array(ptrResult);
}

nttw_integer* padData_Rader(nttw_integer *data, const size_t p, const size_t newSize, const nttw_integer primeDash)
{
    const size_t pad = newSize - p + 1;
    size_t j;
    nttw_integer *newData;

    newData = array_1D(newSize);
    init_1D(newData,newSize,0); //Necessary

    newData[0] = data[0]%primeDash;
    for (j = 1; j < p-1; j ++) ///Copy data leaving N'- p + 1 zeros between zeroth and first element
        newData[j + pad] = data[j]%primeDash;

    return newData;
}

nttw_integer* padTransformMatrix_Rader(nttw_integer *transData, const size_t p, const size_t newSize)
{
    size_t j;
    nttw_integer *newData;

    newData = array_1D(newSize);

    for (j = 0; j < newSize; j ++) ///Copy data leaving N'- p + 1 zeros between zeroth and first element
        newData[j] = transData[j%(p-1)];

    return newData;
}

void extractResult_Rader(nttw_integer *paddedData, nttw_integer *data, const size_t p)
{
    size_t j;

    for (j = 0; j < p-1; j ++) ///Copy data from the first N-1 elements
        data[j] = paddedData[j];
}

void ntt_norm(nttw_integer *data, const size_t n)
{
    size_t j;
    nttw_big_integer inv = 0, d, y;

    d = euclidean((nttw_big_integer)n,MODULUS,&inv,&y); //Multi Inv of p-1
    inv = (inv + MODULUS)%MODULUS; //Ensure x is positive

    for (j = 0; j < n; j++)
        data[j] = (data[j] * inv)%MODULUS;
}
