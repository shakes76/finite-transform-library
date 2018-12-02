/**
 * FRTW Radon Library
 * \file radon.c
 * \brief Radon Transforms Source for the FRTW C Library.
 *
 * This header provides all the discrete operations related to the Radon Transform.
 * It includes the Discrete Radon Transform (FRT) as well as other discretized
 * versions of the Radon Transform.
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
#include <nttw/image.h>
#include <nttw/prime.h>
#include <nttw/number.h>

#include "radon.h"
#include "array_complex.h"

nttw_integer getIsum(nttw_integer *radon, const size_t n)
{
    size_t j;
    nttw_integer Isum = 0;

    for(j = 0; j < n; j ++)
        Isum += radon[j]; //Sum first projection

    return Isum;
}

double getIsum_double(double *radon, const size_t n)
{
    size_t j;
    double Isum = 0;

    for(j = 0; j < n; j ++)
        Isum += radon[j]; //Sum first projection

    return Isum;
}

nttw_integer getIsum_Integer(nttw_integer *radon, const size_t n, const nttw_integer modulus)
{
    size_t j;
    nttw_integer Isum = 0;

    for(j = 0; j < n; j ++)
    {
        Isum += radon[j]; //Sum first projection
        Isum %= modulus;
    }

    return Isum;
}

void drt(const nttw_integer *field, nttw_integer *bins, const int p)
{
    int t = 0, y = 0, m = 0, x = 0, index, row;
    init_1D(bins,(p+1)*p,0); ///Init

    ///Get Projections where m denotes projection number
    for (m = 0, y = 0; m < p; m ++) ///Projection angle
    {
        for (x = 0; x < p; x ++) ///Row
        {
            row = x*p; ///Compute Row
            for (t = 0; t < p; t ++) ///Translate
            {
                index = row + (t + y)%p; ///Compute Column (in 1D)
                bins[m*p+t] += field[index];
            }
            y += m; //!< Next pixel
        }
    }

    for (t = 0; t < p; t ++) ///Columns
        for (y = 0; y < p; y ++) ///Translate
            bins[p*p+t] += field[t*p+y];
}

void drt_dyadic(const nttw_integer *field, nttw_integer *bins, const int n)
{
    int t = 0, y = 0, m = 0, s = 0, x = 0, index, row;
    init_1D(bins,(n+n/2)*n,0); ///Init

    ///Get Projections where m denotes projection number
    for (m = 0, y = 0; m < n; m ++) ///Projection angle
    {
        for (x = 0; x < n; x ++) ///Row
        {
            row = x*n; ///Compute Row
            for (t = 0; t < n; t ++) ///Translate
            {
                index = row + (t + y)%n; ///Compute Column (in 1D)
                bins[m*n+t] += field[index];
            }
            y += m; //!< Next pixel
        }
    }

    for (s = 0, x = 0; s < n/2; s ++) ///Perp Projection Angle
    {
        for (y = 0; y < n; y ++) ///Column
        {
            for (t = 0; t < n; t ++) ///Translate
            {
                index = ( (x+t)%n )*n + y; ///Compute Row (in 1D)
                bins[(s+n)*n+t] += field[index];
            }
            x += 2*s; //!< Next pixel
        }
    }
}

void drt_blockcopy(nttw_integer *datain, nttw_integer *dataout, const int p)
{
/// variables /
    int Rrow,Rcol,Irow,Icol=p*(p+1),wrap;
    nttw_integer *tempi, *tempo;

/// Initialise /
    tempo = dataout;
    while (Icol--)
        *tempo++ = 0;

/// transform m = 0 /
    tempi = datain;
    tempo = dataout;
    for (Irow=0;Irow<p;Irow++)
    {
        for (Icol=0;Icol<p;Icol++)
            *tempo += *tempi++;
        tempo++;
    }
/// transform (1 < m < p-1)
    tempi = datain;
    for (Irow=0;Irow<p;Irow++)
    {
        tempo = dataout+p;
        wrap = p;
        for (Rrow=1;Rrow<=p;Rrow++)
        {
            wrap -= Irow;
            if (wrap<0)
            {
                wrap+=p;
                tempi -= p;
            }
            tempi += Irow;
            for (Rcol=0;Rcol<wrap;Rcol++)
                *tempo++ += *tempi++;
            tempi -= p;
            for (Rcol=wrap;Rcol<p;Rcol++)
                *tempo++ += *tempi++;
        }
        if (Irow == 0)
            tempi += p;
    }
}

void frt(nttw_integer *field, nttw_integer *bins, const int p)
{
    const int size = p*p;
    int m, j, fftRank = 2, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;

    sizes[0] = p;
    sizes[1] = p;
    complexField = fftw_array(size);
    fftResult = fftw_array(size);
    slice = fftw_array(p);
    projComplex = fftw_array(p);

    ///FFT image
    array_to_fftw_array(field,complexField,size);
    fft(fftRank,sizes,complexField,fftResult);

    ///Extract slices
    fftRank = 1;
    for(m = 0; m < p+1; m ++)
    {
        getSlice(m,fftResult,slice,p);
        ifft(fftRank,sizes,slice,projComplex);
        //Copy and norm
        for (j = 0; j < p; j ++)
            bins[m*p+j] = (nttw_integer) ( projComplex[j][0]/(double)p + 0.5 ); //Round because of epsilons
            //bins[m*p+j] = (nttw_integer) ( projComplex[j][0] + 0.5 ); //Round because of epsilons
    }

    ///Clean up
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);
}

void frt_double(double *field, double *bins, const int p)
{
    const int size = p*p;
    int m, j, fftRank = 2, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;

    sizes[0] = p;
    sizes[1] = p;
    complexField = fftw_array(size);
    fftResult = fftw_array(size);
    slice = fftw_array(p);
    projComplex = fftw_array(p);

    ///FFT image
    for(j = 0; j < size; j ++) ///Copy data to real part
        complexField[j][0] = field[j]; //[0] is the real part
    fft(fftRank,sizes,complexField,fftResult);

    ///Extract slices
    fftRank = 1;
    for(m = 0; m < p+1; m ++)
    {
        getSlice(m,fftResult,slice,p);
        ifft(fftRank,sizes,slice,projComplex);
        //Copy and norm
        for (j = 0; j < p; j ++)
            bins[m*p+j] = projComplex[j][0]/(double)p; //Round because of epsilons
            //bins[m*p+j] = projComplex[j][0]; //Round because of epsilons
    }

    ///Clean up
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);
}

void frt_dyadic(nttw_integer *field, nttw_integer *bins, const int n)
{
    const int size = n*n;
    int m, j, fftRank = 2, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;

    sizes[0] = n;
    sizes[1] = n;
    complexField = fftw_array(size);
    fftResult = fftw_array(size);
    slice = fftw_array(n);
    projComplex = fftw_array(n);

    ///FFT image
    array_to_fftw_array(field,complexField,size);
    fft(fftRank,sizes,complexField,fftResult);

    ///Extract slices
    fftRank = 1;
    for(m = 0; m < n; m ++)
    {
        getSlice(m,fftResult,slice,n);

        for(j = 0; j < n; j ++)
        {
            slice[j][0] /= n;
            slice[j][1] /= n;
        }

        ifft(fftRank,sizes,slice,projComplex);

        //Copy and norm
        for (j = 0; j < n; j ++)
            //bins[m*n+j] = (nttw_integer) ( projComplex[j][0]/(double)n + 0.5 ); //Round because of epsilons
            bins[m*n+j] = (nttw_integer) ( projComplex[j][0] + 0.5 ); //Round because of epsilons
    }

    for(m = 0; m < n/2; m ++)
    {
        getSlice_Perp(m,fftResult,slice,n);

        for(j = 0; j < n; j ++)
        {
            slice[j][0] /= n;
            slice[j][1] /= n;
        }

        ifft(fftRank,sizes,slice,projComplex);
        //Copy and norm
        for (j = 0; j < n; j ++)
            //bins[(m+n)*n+j] = (nttw_integer) ( projComplex[j][0]/(double)n + 0.5 ); //Round because of epsilons
            bins[(m+n)*n+j] = (nttw_integer) ( projComplex[j][0] + 0.5 ); //Round because of epsilons
    }

    ///Clean up
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);
}

void frt_dyadic_double(double *field, double *bins, const int n)
{
    const int size = n*n;
    int m, j, fftRank = 2, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;

    sizes[0] = n;
    sizes[1] = n;
    complexField = fftw_array(size);
    fftResult = fftw_array(size);
    slice = fftw_array(n);
    projComplex = fftw_array(n);

    ///FFT image
    for(j = 0; j < size; j ++) ///Copy data to real part
        complexField[j][0] = field[j]; //[0] is the real part
    fft(fftRank,sizes,complexField,fftResult);

    ///Extract slices
    fftRank = 1;
    for(m = 0; m < n; m ++)
    {
        getSlice(m,fftResult,slice,n);

        for(j = 0; j < n; j ++)
        {
            slice[j][0] /= n;
            slice[j][1] /= n;
        }

        ifft(fftRank,sizes,slice,projComplex);

        //Copy and norm
        for (j = 0; j < n; j ++)
            //bins[m*n+j] = projComplex[j][0]/(double)n; //Round because of epsilons
            bins[m*n+j] = projComplex[j][0]; //Round because of epsilons
    }

    for(m = 0; m < n/2; m ++)
    {
        getSlice_Perp(m,fftResult,slice,n);

        for(j = 0; j < n; j ++)
        {
            slice[j][0] /= n;
            slice[j][1] /= n;
        }

        ifft(fftRank,sizes,slice,projComplex);
        //Copy and norm
        for (j = 0; j < n; j ++)
            //bins[(m+n)*n+j] = projComplex[j][0]/(double)n; //Round because of epsilons
            bins[(m+n)*n+j] = projComplex[j][0]; //Round because of epsilons
    }

    ///Clean up
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);
}

void frt_dyadic_signed(long *field, long *bins, const int n)
{
    const int size = n*n;
    int m, j, fftRank = 2, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;

    sizes[0] = n;
    sizes[1] = n;
    complexField = fftw_array(size);
    fftResult = fftw_array(size);
    slice = fftw_array(n);
    projComplex = fftw_array(n);

    ///FFT image
    arraySigned_to_fftw_array(field,complexField,size);
    fft(fftRank,sizes,complexField,fftResult);

    ///Extract slices
    fftRank = 1;
    for(m = 0; m < n; m ++)
    {
        getSlice(m,fftResult,slice,n);

        for(j = 0; j < n; j ++)
        {
            slice[j][0] /= n;
            slice[j][1] /= n;
        }

        ifft(fftRank,sizes,slice,projComplex);

        //Copy and norm
        for (j = 0; j < n; j ++)
        {
            if(projComplex[j][0] < 0)
                bins[m*n+j] = (long) ( projComplex[j][0] - 0.5 ); //Round because of epsilons
            else
                bins[m*n+j] = (long) ( projComplex[j][0] + 0.5 ); //Round because of epsilons
        }
    }

    for(m = 0; m < n/2; m ++)
    {
        getSlice_Perp(m,fftResult,slice,n);

        for(j = 0; j < n; j ++)
        {
            slice[j][0] /= n;
            slice[j][1] /= n;
        }

        ifft(fftRank,sizes,slice,projComplex);
        //Copy and norm
        for (j = 0; j < n; j ++)
        {
            if(projComplex[j][0] < 0)
                bins[m*n+j] = (long) ( projComplex[j][0] - 0.5 ); //Round because of epsilons
            else
                bins[m*n+j] = (long) ( projComplex[j][0] + 0.5 ); //Round because of epsilons
        }
    }

    ///Clean up
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);
}

nttw_integer idrt(nttw_integer *bins, nttw_integer *result, const int p)
{
    nttw_integer Isum = 0;
    const int translate_x = 1;
    int translate_y = 0, y = 0, j = 0, m = 0, k = 0;
    init_1D(result,p*p,0);

    ///Calculate the Isum from the first row
    for (k = 0; k < p; k ++)
        Isum += bins[k];

    ///Get Projections where j denotes projection number
    for (j = 0, y = 0; j < p; j ++)
    {
        for (m = 0; m < p; m++)
        {
            y += p-j;
            k = (m + translate_x)%p;
            for (translate_y = 0; translate_y < p; translate_y ++)
                result[j*p+translate_y] += bins[k*p+(translate_y + y)%p];
        }
    }

    ///Get the last projection of rows
    for (j = 0; j < p; j ++)
        for (k = 0; k < p; k++)
        {
            result[j*p+k] += bins[p*p+j];
            result[j*p+k] = (result[j*p+k] - Isum)/p;
        }

    return Isum;
}

nttw_integer idrt_blockcopy(nttw_integer *datain, nttw_integer *dataout, const int p)
{
/// variables
    int Icol = p*p, Irow, Rcol, Rrow, m, wrap;
    nttw_integer *tempi, *tempo, total = 0;

/// Initialise
    tempo = dataout;
    while (Icol--)
        *tempo++ = 0;

/// inverse transform m = 0
    tempo = dataout;
    tempi = datain;
    for (Irow=0;Irow<p;Irow++)
    {
        for (Icol=0;Icol<p;Icol++)
            *tempo++ = *tempi;
        total += *tempi++;
    }

/// inverse transform (1 < m < p-1)
    for (Rrow=1;Rrow<=p;Rrow++)
    {
        tempo = dataout;
        wrap = p;
        m = p-Rrow;
        for (Irow=0;Irow<p;Irow++)
        {
            if (wrap<0)
            {
                wrap+=p;
                tempi-=p;
            }
            for (Rcol=0;Rcol<wrap;Rcol++)
                *tempo++ += *tempi++;
            tempi-=p;
            for (Rcol=wrap;Rcol<p;Rcol++)
                *tempo++ += *tempi++;
            wrap-=m;
            tempi+=m;
        }
    }
/// subtract total and divide by p
    tempo = dataout;
    Icol = p*p;
    while (Icol--)
        *tempo++ = (*tempo - total)/p; //!< \todo Leak Identified !!! Not accessed but pointer ends up pointing outside.

    return total;
}

nttw_integer ifrt(nttw_integer *bins, nttw_integer *result, const int p, const int norm)
{
    const int size = p*p;
    int m, j, fftRank = 1, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;
    nttw_integer Isum = 0;

    Isum = getIsum(bins,p);
    sizes[0] = p;
    sizes[1] = p;
    complexField = fftw_array(size);
    fftResult = fftw_array(size);
    fftw_init_1D(fftResult,size,0);
    slice = fftw_array(p);
    projComplex = fftw_array(p);

    ///FFT projections
    for(m = 0; m < p+1; m ++)
    {
        for (j = 0; j < p; j ++)
        {
            projComplex[j][0] = (double) bins[m*p+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        setSlice(m,fftResult,slice,p);
    }

    ///Inverse FFT to get result
    fftRank = 2;
    ifft(fftRank,sizes,fftResult,complexField);
    fftw_array_to_array(complexField,result,size);
    if(norm)
    {
        for(m = 0; m < p; m ++)
            for (j = 0; j < p; j ++)
                result[m*p+j] = (nttw_integer) (result[m*p+j]/(double)p - Isum + 0.5)/p; //Avoid truncation
                //result[m][j] = (nttw_integer) (result[m][j]/(double)size + 0.5); //Avoid truncation
    }
    /*else
    {
        for(m = 0; m < p; m ++)
            for (j = 0; j < p; j ++)
                result[m][j] = (nttw_integer) (result[m][j]/(double)p - Isum + 0.5); //Avoid truncation
    }*/

    ///Cleanup
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);

    return Isum;
}

double ifrt_double(double *bins, double *result, const int p, const int norm)
{
    const int size = p*p;
    int m, j, fftRank = 1, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;
    double Isum = 0;

    Isum = getIsum_double(bins,p);
    sizes[0] = p;
    sizes[1] = p;
    complexField = fftw_array(size);
    fftResult = fftw_array(size);
    fftw_init_1D(fftResult,size,0);
    slice = fftw_array(p);
    projComplex = fftw_array(p);

    ///FFT projections
    for(m = 0; m < p+1; m ++)
    {
        for (j = 0; j < p; j ++)
        {
            projComplex[j][0] = bins[m*p+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        setSlice(m,fftResult,slice,p);
    }

    ///Inverse FFT to get result
    fftRank = 2;
    ifft(fftRank,sizes,fftResult,complexField);
    for(j = 0; j < size; j ++) ///Copy data from real part
        result[j] = complexField[j][0]; //[0] is the real part
    if(norm)
    {
        for(m = 0; m < p; m ++)
            for (j = 0; j < p; j ++)
                result[m*p+j] = ( result[m*p+j]/(double)p - Isum )/p; //Avoid truncation
                //result[m][j] = (result[m][j]/(double)size); //Avoid truncation
    }
    /*else
    {
        for(m = 0; m < p; m ++)
            for (j = 0; j < p; j ++)
                result[m][j] = (result[m][j]/(double)p - Isum); //Avoid truncation
    }*/

    ///Cleanup
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);

    return Isum;
}

nttw_integer ifrt_signed(nttw_integer *bins, long *result, const int p, const int norm)
{
    const int size = p*p;
    int m, j, fftRank = 1, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;
    long Isum = 0;

    Isum = getIsum(bins,p);
    sizes[0] = p;
    sizes[1] = p;
    complexField = fftw_array(size);
    fftResult = fftw_array(size);
    fftw_init_1D(fftResult,size,0);
    slice = fftw_array(p);
    projComplex = fftw_array(p);

    ///FFT projections
    for(m = 0; m < p+1; m ++)
    {
        for (j = 0; j < p; j ++)
        {
            projComplex[j][0] = (double) bins[m*p+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        setSlice(m,fftResult,slice,p);
    }

    ///Inverse FFT to get result
    fftRank = 2;
    ifft(fftRank,sizes,fftResult,complexField);
    fftw_array_to_arraySigned(complexField,result,size);
    if(norm)
    {
        for(m = 0; m < p; m ++)
            for (j = 0; j < p; j ++)
            {
                if(result[m*p+j] < 0)
                    result[m*p+j] = (long) (result[m*p+j]/(double)p - Isum - 0.5)/p; //Avoid truncation
                else
                    result[m*p+j] = (long) (result[m*p+j]/(double)p - Isum + 0.5)/p; //Avoid truncation
            }
    }
    else
    {
        for(m = 0; m < p; m ++)
            for (j = 0; j < p; j ++)
            {
                if(result[m*p+j] < 0)
                    result[m*p+j] = (long) (result[m*p+j]/(double)p - Isum - 0.5); //Avoid truncation
                else
                    result[m*p+j] = (long) (result[m*p+j]/(double)p - Isum + 0.5); //Avoid truncation
            }
    }

    ///Cleanup
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);

    return Isum;
}

void ifrt2(fftw_complex *bins, fftw_complex *result, fftw_complex *Isum, const int p, const int norm)
{
    const int size = p*p;
    int m, j, fftRank = 1, sizes[2];
    fftw_complex *fftResult, *slice;

    *Isum[0] = 0.0;
    *Isum[1] = 0.0;
    for(j = 0; j < p; j ++)
    {
        *Isum[0] += bins[j][0];
        *Isum[1] += bins[j][1];
    }

    sizes[0] = p;
    sizes[1] = p;
    fftResult = fftw_array(size);
    fftw_init_1D(fftResult,size,0);
    slice = fftw_array(p);

    ///FFT projections
    for(m = 0; m < p+1; m ++)
    {
        for (j = 0; j < p; j ++)
        {
            slice[j][0] = bins[m*p+j][0];
            slice[j][1] = bins[m*p+j][1];
        }
        //fft(fftRank,sizes,projComplex,slice);
        setSlice(m,fftResult,slice,p);
    }

    ///Inverse FFT to get result
    fftRank = 2;
    ifft(fftRank,sizes,fftResult,result);
    for(m = 0; m < p; m ++)
        for (j = 0; j < p; j ++)
        {
            result[m*p+j][0] = result[m*p+j][0]/(double)p*p; //Avoid truncation
            result[m*p+j][1] = result[m*p+j][1]/(double)p*p; //Avoid truncation
        }
    if(norm)
    {
        for(m = 0; m < p; m ++)
            for (j = 0; j < p; j ++)
            {
                result[m*p+j][0] = result[m*p+j][0]/(double)p; //Avoid truncation
                result[m*p+j][1] = result[m*p+j][1]/(double)p; //Avoid truncation
            }
    }
    /*else
    {
        for(m = 0; m < p; m ++)
            for (j = 0; j < p; j ++)
                result[m][j] = (nttw_integer) (result[m][j]/(double)p - Isum + 0.5); //Avoid truncation
    }*/

    ///Cleanup
    fftw_free(fftResult);
    fftw_free(slice);
}

nttw_integer ifrt_dyadic(nttw_integer *bins, nttw_integer *result, const int n, const int norm)
{
    const int size = n*n;
    int m, j, fftRank = 1, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;
    nttw_big_integer *oversampling;
    nttw_integer Isum = 0;

    Isum = getIsum(bins,n);
    oversampling = array_1D_big(n);
    sizes[0] = n;
    sizes[1] = n;
    complexField = fftw_array(size);
    fftw_init_1D(complexField,size,0);
    fftResult = fftw_array(size);
    fftw_init_1D(fftResult,size,0);
    slice = fftw_array(n);
    projComplex = fftw_array(n);

    dyadic_1D_filter(oversampling,n);

    ///FFT projections
    for(m = 0; m < n; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = (double) bins[m*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }
        setSlice(m,fftResult,slice,n);
    }

    for(m = 0; m < n/2; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = (double) bins[(m+n)*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }
        setSlice_Perp(m,fftResult,slice,n);
    }

    ///Correct for oversampling
    //dyadic_oversample(oversampling,n); ///Get oversample filter (no need to orient for Fourier)
    //filter_oversampling(fftResult,oversampling,n);

    ///Inverse FFT to get result
    fftRank = 2;
    ifft(fftRank,sizes,fftResult,complexField);
    fftw_array_to_array(complexField,result,size);
    if(norm)
    {
        for(m = 0; m < n; m ++)
            for (j = 0; j < n; j ++)
                //result[m][j] = (nttw_integer) (result[m][j]/(double)n - Isum + 0.5)/n; //Avoid truncation
                result[m*n+j] = (nttw_integer) (result[m*n+j]/(double)n + 0.5); //Avoid truncation
    }

    ///Cleanup
    free_array(oversampling);
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);

    return Isum;
}

double ifrt_dyadic_double(double *bins, double *result, const int n, const int norm)
{
    const int size = n*n;
    int m, j, fftRank = 1, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;
    nttw_big_integer *oversampling;
    double Isum = 0;

    Isum = getIsum_double(bins,n);
    oversampling = array_1D_big(n);
    sizes[0] = n;
    sizes[1] = n;
    complexField = fftw_array(size);
    fftw_init_1D(complexField,size,0);
    fftResult = fftw_array(size);
    fftw_init_1D(fftResult,size,0);
    slice = fftw_array(n);
    projComplex = fftw_array(n);

    dyadic_1D_filter(oversampling,n);

    ///FFT projections
    for(m = 0; m < n; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = bins[m*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }
        setSlice(m,fftResult,slice,n);
    }

    for(m = 0; m < n/2; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = bins[(m+n)*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }
        setSlice_Perp(m,fftResult,slice,n);
    }

    ///Correct for oversampling
    //dyadic_oversample(oversampling,n); ///Get oversample filter (no need to orient for Fourier)
    //filter_oversampling(fftResult,oversampling,n);

    ///Inverse FFT to get result
    fftRank = 2;
    ifft(fftRank,sizes,fftResult,complexField);
    for(j = 0; j < size; j ++) ///Copy data from real part
        result[j] = complexField[j][0]; //[0] is the real part
    if(norm)
    {
        for(m = 0; m < n; m ++)
            for (j = 0; j < n; j ++)
                //result[m][j] = (result[m][j]/(double)n - Isum)/n; //Avoid truncation
                result[m*n+j] = result[m*n+j]/(double)n; //Avoid truncation
    }

    ///Cleanup
    free_array(oversampling);
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);

    return Isum;
}

nttw_integer ifrt_dyadic_signed(nttw_integer *bins, long *result, const int n, const int norm)
{
    const int size = n*n;
    int m, j, fftRank = 1, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;
    nttw_big_integer *oversampling;
    nttw_integer Isum = 0;

    Isum = getIsum(bins,n);
    oversampling = array_1D_big(n);
    sizes[0] = n;
    sizes[1] = n;
    complexField = fftw_array(size);
    fftw_init_1D(complexField,size,0);
    fftResult = fftw_array(size);
    fftw_init_1D(fftResult,size,0);
    slice = fftw_array(n);
    projComplex = fftw_array(n);

    dyadic_1D_filter(oversampling,n);

    ///FFT projections
    for(m = 0; m < n; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = (double) bins[m*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }
        setSlice(m,fftResult,slice,n);
    }

    for(m = 0; m < n/2; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = (double) bins[(m+n)*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }
        setSlice_Perp(m,fftResult,slice,n);
    }

    ///Correct for oversampling
    //dyadic_oversample(oversampling,n); ///Get oversample filter (no need to orient for Fourier)
    //filter_oversampling(fftResult,oversampling,n);

    ///Inverse FFT to get result
    fftRank = 2;
    ifft(fftRank,sizes,fftResult,complexField);
    fftw_array_to_arraySigned(complexField,result,size);
    if(norm)
    {
        for(m = 0; m < n; m ++)
            for (j = 0; j < n; j ++)
            {
                if(result[m*n+j] < 0)
                    result[m*n+j] = (long) (result[m*n+j]/(double)n - 0.5); //Avoid truncation
                else
                    result[m*n+j] = (long) (result[m*n+j]/(double)n + 0.5); //Avoid truncation
            }
    }

    ///Cleanup
    free_array(oversampling);
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);

    return Isum;
}

nttw_integer ifrt_dyadic_signed2(long *bins, long *result, const int n, const int norm)
{
    const int size = n*n;
    int m, j, fftRank = 1, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;
    nttw_big_integer *oversampling;
    nttw_integer Isum = 0;

    //Isum = getIsum(bins,n);
    oversampling = array_1D_big(n);
    sizes[0] = n;
    sizes[1] = n;
    complexField = fftw_array(size);
    fftw_init_1D(complexField,size,0);
    fftResult = fftw_array(size);
    fftw_init_1D(fftResult,size,0);
    slice = fftw_array(n);
    projComplex = fftw_array(n);

    dyadic_1D_filter(oversampling,n);

    ///FFT projections
    for(m = 0; m < n; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = (double) bins[m*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }
        setSlice(m,fftResult,slice,n);
    }

    for(m = 0; m < n/2; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = (double) bins[(m+n)*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }
        setSlice_Perp(m,fftResult,slice,n);
    }

    ///Correct for oversampling
    //dyadic_oversample(oversampling,n); ///Get oversample filter (no need to orient for Fourier)
    //filter_oversampling(fftResult,oversampling,n);

    ///Inverse FFT to get result
    fftRank = 2;
    ifft(fftRank,sizes,fftResult,complexField);
    fftw_array_to_arraySigned(complexField,result,size);
    if(norm)
    {
        for(m = 0; m < n; m ++)
            for (j = 0; j < n; j ++)
            {
                if(result[m*n+j] < 0)
                    result[m*n+j] = (long) (result[m*n+j]/(double)n - 0.5); //Avoid truncation
                else
                    result[m*n+j] = (long) (result[m*n+j]/(double)n + 0.5); //Avoid truncation
            }
    }

    ///Cleanup
    free_array(oversampling);
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);

    return Isum;
}

nttw_integer ifrt_dyadic_unsample(nttw_integer *bins, long *result, nttw_big_integer *sampling, const int n, const int norm)
{
    const int size = n*n;
    int m, j, fftRank = 1, sizes[2];
    fftw_complex *complexField, *fftResult, *slice, *projComplex;
    nttw_big_integer *oversampling;
    nttw_integer Isum = 0;

    Isum = getIsum(bins,n);
    oversampling = array_1D_big(size);
    sizes[0] = n;
    sizes[1] = n;
    complexField = fftw_array(size);
    fftw_init_1D(complexField,size,0);
    fftResult = fftw_array(size);
    fftw_init_1D(fftResult,size,0);
    slice = fftw_array(n);
    projComplex = fftw_array(n);

    //dyadic_1D_filter(oversampling,n);

    ///FFT projections
    for(m = 0; m < n; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = (double) bins[m*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        /*for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }*/
        setSlice(m,fftResult,slice,n);
    }

    for(m = 0; m < n/2; m ++)
    {
        for (j = 0; j < n; j ++)
        {
            projComplex[j][0] = (double) bins[(m+n)*n+j];
            projComplex[j][1] = 0.0;
        }
        fft(fftRank,sizes,projComplex,slice);
        /*for (j = 0; j < n; j ++)
        {
            slice[j][0] /= oversampling[j]*n;
            slice[j][1] /= oversampling[j]*n;
        }*/
        setSlice_Perp(m,fftResult,slice,n);
    }
    printf("Unfiltered FFT Space:\n");
    for(j = 0; j < n; j ++)
    {
        for(m = 0; m < n; m ++)
            printf("%f, ",fftResult[j*n+m][0]);
        printf("\n");
    }

    ///Correct for oversampling
    dyadic_oversample(oversampling,n); ///Get oversample filter (no need to orient for Fourier)
    printf("Filter:\n");
    for(j = 0; j < n; j ++)
    {
        for(m = 0; m < n; m ++)
            printf(FORMAT_BIG_CSVOUTPUT_STR,oversampling[j*n+m]);
        printf("\n");
    }
    ///Unsample
    for(j = 0; j < n; j ++)
        for(m = 0; m < n; m ++)
        {
//            if(oversampling[j*n+m] > 1)
//                oversampling[j*n+m] -= sampling[j*n+m];
            if(sampling[j*n+m] > 0)
                oversampling[j*n+m] = 0;
        }
    printf("Filter Unsampled:\n");
    for(j = 0; j < n; j ++)
    {
        for(m = 0; m < n; m ++)
            printf(FORMAT_BIG_CSVOUTPUT_STR,oversampling[j*n+m]);
        printf("\n");
    }
    ///Filter
    filter_oversampling(fftResult,oversampling,n);
    printf("FFT Space:\n");
    for(j = 0; j < n; j ++)
    {
        for(m = 0; m < n; m ++)
            printf("%f, ",fftResult[j*n+m][0]);
        printf("\n");
    }

    ///Inverse FFT to get result
    fftRank = 2;
    ifft(fftRank,sizes,fftResult,complexField);
    printf("Inverted FFT Space:\n");
    for(j = 0; j < n; j ++)
    {
        for(m = 0; m < n; m ++)
            printf("%f, ",complexField[j*n+m][0]);
        printf("\n");
    }
    fftw_array_to_arraySigned(complexField,result,size);
    if(norm)
    {
        for(m = 0; m < n; m ++)
            for (j = 0; j < n; j ++)
                //result[m][j] = (nttw_integer) (result[m][j]/(double)n - Isum + 0.5)/n; //Avoid truncation
                if(result[m*n+j] < 0)
                    result[m*n+j] = (long) (result[m*n+j]/(double)(n*n) - 0.5); //Avoid truncation
                else
                    result[m*n+j] = (long) (result[m*n+j]/(double)(n*n) + 0.5); //Avoid truncation
    }
    else
    {
        for(m = 0; m < n; m ++)
            for (j = 0; j < n; j ++)
                //result[m][j] = (nttw_integer) (result[m][j]/(double)n - Isum + 0.5)/n; //Avoid truncation
                if(result[m*n+j] < 0)
                    result[m*n+j] = (long) (result[m*n+j]/(double)n - 0.5); //Avoid truncation
                else
                    result[m*n+j] = (long) (result[m*n+j]/(double)n + 0.5); //Avoid truncation
    }

    ///Cleanup
    free_array(oversampling);
    fftw_free(complexField);
    fftw_free(fftResult);
    fftw_free(slice);
    fftw_free(projComplex);

    return Isum;
}

void nrt(nttw_integer *field, nttw_integer *bins, const size_t p)
{
    const size_t size = p*p;
    size_t m, j;
    const int totalRows = 1000;
    const int totalCols = 10;
    const char *filename = "primes.csv";
    nttw_integer *nttResult, *slice, *ptrData, *primes;
    nttw_integer root, modulus, proot;

    nttResult = array_1D(size);
    ptrData = array_1D(p);
    slice = array_1D(p);

    init_1D(nttResult,size,0);

    ///Read primes list
    if( !readCSV(&primes,totalRows,totalCols,filename) )
    {
        printf(">| Could not open primes file.\n");
        exit(EXIT_FAILURE);
    }
    ///Compute roots
    root = findFirstPrimitiveRoot(p);
    modulus = findAlternatePrime(primes,totalRows*totalCols,p);
    proot = findFirstPrimitiveRoot(modulus);
    printf(">| Root: %i, Prime': %i, Primitive Root: %i\n",root,modulus,proot);

    ///FNTT image
    fntt_2D_prime(field,nttResult,p,root,modulus,proot,NTTW_FORWARD,TRUE); //Normalises only the NTT not the FRT

    ///Set slices
    for(m = 0; m < p+1; m ++)
    {
        getSlice_Integer(m,nttResult,slice,p);

        fntt_prime(slice,ptrData,p,root,modulus,proot,NTTW_INVERSE,TRUE);

        for (j = 0; j < p; j ++)
            bins[m*p+j] = ptrData[j];
    }

    ///Clean up
    free_array(primes);
    free_array(nttResult);
    free_array(slice);
    free_array(ptrData);
}

void nrt2(nttw_integer *field, nttw_integer *nttOfBins, const size_t p, const nttw_integer root, const nttw_integer modulus, const nttw_integer proot)
{
    const size_t size = p*p;
    size_t m, j;
    nttw_integer *nttResult, *slice;

    nttResult = array_1D(size);
    slice = array_1D(p);

    ///FNTT image
    fntt_2D_prime(field,nttResult,p,root,modulus,proot,NTTW_FORWARD,TRUE); //Forward doesn't Normalise
    //printf("FNTT DC: %lu\n", nttResult[0]);

    ///Set slices
    for(m = 0; m < p+1; m ++)
    {
        getSlice_Integer(m,nttResult,slice,p);

        for(j = 0; j < p; j ++)
            nttOfBins[m*p+j] = slice[j];
    }

    ///Clean up
    free_array(nttResult);
    free_array(slice);
}

void nrt_dyadic(nttw_integer *field, nttw_integer *bins, const size_t n)
{
    const size_t size = n*n;
    size_t m, j;
    nttw_integer *nttResult, *slice;
    nttw_big_integer inv, d, y;

    nttResult = array_1D(size);
    slice = array_1D(n);

    init_1D(nttResult,size,0);

    d = euclidean((nttw_big_integer)n,MODULUS,&inv,&y); //Multi Inv of p-1
    inv = (inv + MODULUS)%MODULUS; //Ensure x is positive

    ///FNTT image
    fntt_2D(field,nttResult,n,NTTW_FORWARD); //Normalises only the NTT not the FRT

    ///Set slices
    for(m = 0; m < n; m ++)
    {
        getSlice_Integer(m,nttResult,slice,n);

        fntt(slice,n,PROOT,NTTW_INVERSE);

        for (j = 0; j < n; j ++)
            bins[m*n+j] = slice[j];
    }

    for(m = 0; m < n/2; m ++)
    {
        getSlice_Perp_Integer(m,nttResult,slice,n);

        fntt(slice,n,PROOT,NTTW_INVERSE);

        for (j = 0; j < n; j ++)
            bins[(m+n)*n+j] = slice[j];
    }

    ///Clean up
    free_array(nttResult);
    free_array(slice);
}

nttw_integer inrt(nttw_integer *bins, nttw_integer *result, const int p, const int norm)
{
    const int size = p*p;
    const int totalRows = 1000;
    const int totalCols = 10;
    const char *filename = "primes.csv";
    int m, j;
    const int normTransform = FALSE;
    nttw_integer *nttResult, *slice, *primes, *ptrData, Isum = 0;
    nttw_integer root, modulus, proot;
    nttw_big_integer d, x, y;

    nttResult = array_1D(size);
    slice = array_1D(p);

    init_1D(nttResult,size,0);

    ///Read primes list
    if( !readCSV(&primes,totalRows,totalCols,filename) )
    {
        printf(">| Could not open primes file.\n");
        exit(EXIT_FAILURE);
    }
    ///Compute roots
    root = findFirstPrimitiveRoot(p);
    modulus = findAlternatePrime(primes,totalRows*totalCols,p);
    proot = findFirstPrimitiveRoot(modulus);
    Isum = getIsum(bins,p);
    Isum %= modulus;
    //printf(">| Root: %i, Prime': %i, Primitive Root: %i, Isum (mod): %i\n",root,modulus,proot,Isum);

    ///Set slices (0 <= m < p)
    for(m = 0; m < p; m ++)
    {
        ptrData = &bins[m*p];
        fntt_prime(ptrData,slice,p,root,modulus,proot,NTTW_FORWARD,normTransform);
        setSlice_Integer(m,nttResult,slice,p,modulus);
    }
    ///m = p slice
    ptrData = &bins[p*p];
    fntt_prime(ptrData,slice,p,root,modulus,proot,NTTW_FORWARD,normTransform);
    setSlice_Integer(p,nttResult,slice,p,modulus);

    d = euclidean(p,modulus,&x,&y); ///Inverse of DC
    x = (x + modulus)%modulus; //Make positive
    printf("Modulus Inv: %lu, Isum: %li\n",x,Isum);

    ///iFNTT 2D image
    fntt_2D_prime(nttResult,result,p,root,modulus,proot,NTTW_INVERSE,FALSE);

    ///Norm image
    for(m = 0; m < p; m ++)
        for(j = 0; j < p; j ++)
        {
            result[m*p+j] = (result[m*p+j] * (nttw_big_integer)x)%modulus; ///Norm
            result[m*p+j] = MODSUB(result[m*p+j],Isum,modulus);
            if(norm)
                result[m*p+j] = (result[m*p+j] * (nttw_big_integer)x)%modulus; ///Norm
        }

    ///Clean up
    free_array(nttResult);
    free_array(slice);
    free_array(primes);

    return Isum;
}

nttw_integer inrt2(nttw_integer *bins, nttw_integer *result, const int p, const nttw_integer root, const nttw_integer modulus, const nttw_integer proot, const int norm)
{
    const int size = p*p;
    int m, n;
    nttw_integer *nttResult, *slice, Isum = 0;

    Isum = bins[0]; //First term in NTT of projection is DC
    //printf("Isum: %lu\n", Isum);

    nttResult = array_1D(size);

    init_1D(nttResult,size,0);

    ///Set slices (0 <= m < p)
    for(m = 0; m < p; m ++)
    {
        slice = &bins[m*p];
        setSlice_Integer(m,nttResult,slice,p,modulus);
    }
    ///m = p slice
    slice = &bins[p*p];
    setSlice_Integer(p,nttResult,slice,p,modulus);

    ///Norm image

    ///iFNTT 2D image
    fntt_2D_prime(nttResult,result,p,root,modulus,proot,NTTW_INVERSE,norm);
    //printf("iFNTT DC: %lu\n", nttResult[0]);

    for(m = 0; m < p; m ++)
        for(n = 0; n < p; n ++)
            result[m*p+n] = MODSUB(result[m*p+n],Isum,modulus);

    ///Clean up
    free_array(nttResult);

    return Isum;
}

nttw_integer inrt_dyadic(nttw_integer *bins, nttw_integer *result, const size_t n, const int norm)
{
    const size_t size = n*n;
    size_t m, j;
    nttw_integer *nttResult, *slice, Isum = 0;
    nttw_big_integer inv, d, y;
    nttw_big_integer *oversampling;

    nttResult = array_1D(size);
    oversampling = array_1D_big(n);

    init_1D(nttResult,size,0);

    Isum = getIsum(bins,n);
    dyadic_1D_filter_Integer(oversampling,n,MODULUS);

    if(norm)
    {
        d = euclidean((nttw_big_integer)n,MODULUS,&inv,&y); //Multi Inv of p-1
        inv = (inv + MODULUS)%MODULUS; //Ensure x is positive
        for(j = 0; j < n; j ++)
            oversampling[j] = (inv * oversampling[j])%MODULUS;
    }

    ///Set slices
    for(m = 0; m < n; m ++)
    {
        slice = &bins[m*n];

        fntt(slice,n,PROOT,NTTW_FORWARD);

        for(j = 0; j < n; j ++)
            slice[j] = (slice[j] * oversampling[j])%MODULUS;

        setSlice_Integer(m,nttResult,slice,n,MODULUS);
    }

    for(m = 0; m < n/2; m ++)
    {
        slice = &bins[(m+n)*n];

        fntt(slice,n,PROOT,NTTW_FORWARD);

        for(j = 0; j < n; j ++)
            slice[j] = (slice[j] * oversampling[j])%MODULUS;

        setSlice_Perp_Integer(m,nttResult,slice,n,MODULUS);
    }

    ///Correct for oversampling
    //dyadic_oversample(oversampling,n); ///Get oversample filter
    //filter_oversampling_Integer(nttResult,oversampling,n,MODULUS); //Works

    ///FNTT image
    fntt_2D(nttResult,result,n,NTTW_INVERSE); //Normalises only the NTT not the FRT

    ///Clean up
    free_array(nttResult);
    free_array(oversampling);

    return Isum;
}

void getSlice(int m, fftw_complex *data, fftw_complex *slice, const int size)
{
    int k, index = 0;

    if(m < size && m >= 0)
    {
        slice[0][0] = data[0][0];
        slice[0][1] = data[0][1];
        for(k = 1; k < size; k ++)
        {
            index = ( ( size-(k*m)%size )%size )*size + k;
            slice[k][0] = data[index][0];
            slice[k][1] = data[index][1];
        }
    }
    else if(m == size)
    {
        for(k = 0; k < size; k ++)
        {
            index = k*size;
            slice[k][0] = data[index][0];
            slice[k][1] = data[index][1];
        }
    }
}

void getSlice_Integer(const size_t m, const nttw_integer *data, nttw_integer *slice, const size_t size)
{
    size_t k, index = 0;

    if(m < size)
    {
        slice[0] = data[0];
        for(k = 1; k < size; k ++)
        {
            index = ( ( size-(k*m)%size )%size )*size + k;
            slice[k] = data[index];
        }
    }
    else if(m == size)
    {
        for(k = 0; k < size; k ++)
        {
            index = k*size;
            slice[k] = data[index];
        }
    }
}

void getSlice_Perp(int s, fftw_complex *data, fftw_complex *slice, const int size)
{
    int k, index = 0;

    if(s < size/2 && s >= 0)
    {
        slice[0][0] = data[0][0];
        slice[0][1] = data[0][1];
        for(k = 1; k < size; k ++)
        {
            index = ((size-k*2*s)%size + size)%size + k*size;
            slice[k][0] = data[index][0];
            slice[k][1] = data[index][1];
        }
    }
}

void getSlice_Perp_Integer(const size_t s, const nttw_integer *data, nttw_integer *slice, const size_t size)
{
    size_t k, index = 0;

    if(s < size/2)
    {
        slice[0] = data[0];
        for(k = 1; k < size; k ++)
        {
            index = ((size-k*2*s)%size + size)%size + k*size;
            slice[k] = data[index];
        }
    }
}

void setSlice(int m, fftw_complex *data, fftw_complex *slice, const int size)
{
    int k, index = 0;

    if(m < size && m >= 0)
    {
        data[0][0] += slice[0][0];
        data[0][1] += slice[0][1];
        for(k = 1; k < size; k ++)
        {
            index = ( ( size-(k*m)%size )%size )*size + k;
            data[index][0] += slice[k][0];
            data[index][1] += slice[k][1];
        }
    }
    else if(m == size)
    {
        for(k = 0; k < size; k ++)
        {
            index = k*size;
            data[index][0] += slice[k][0];
            data[index][1] += slice[k][1];
        }
    }
}

void setSlice_Integer(const size_t m, nttw_integer *data, const nttw_integer *slice, const size_t size, const nttw_integer modulus)
{
    size_t k, index = 0;

//    if(m < size && m >= 0)
//    {
//        data[0] += slice[0];
//        data[0] = MODADD2(slice[0], data[0], modulus); //data[0] += slice[0]
//        for(k = 1; k < size; k ++)
//        {
//            index = ( ( size-(k*m)%size )%size )*size + k;
//            data[index] += slice[k];
//            data[index] = MODADD2(slice[k], data[index], modulus); //data[index] += slice[k]
//        }
//    }
//    else if(m == size)
//    {
//        for(k = 0; k < size; k ++)
//        {
//            index = k*size;
//            data[index] += slice[k];
//            data[index] = MODADD2(slice[k], data[index], modulus); //data[index] += slice[k]
//        }
//    }

    if(m < size)
    {
        //data[0] += slice[0];
        //data[0] %= modulus;
        data[0] = MODADD(data[0], slice[0], modulus); //data[0] += slice[0]
        for(k = 1; k < size; k ++)
        {
            index = ( ( size-(k*m)%size )%size )*size + k;
            //data[index] += slice[k];
            //data[index] %= modulus;
            data[index] = MODADD(data[index], slice[k], modulus); //data[index] += slice[k]
        }
    }
    else if(m == size)
    {
        for(k = 0; k < size; k ++)
        {
            index = k*size;
            //data[index] += slice[k];
            //data[index] %= modulus;
            data[index] = MODADD(data[index], slice[k], modulus); //data[index] += slice[k]
        }
    }
}

void setSlice_Perp(int s, fftw_complex *data, fftw_complex *slice, const int size)
{
    int k, index = 0;

    if(s < size/2 && s >= 0)
    {
        data[0][0] += slice[0][0];
        data[0][1] += slice[0][1];
        for(k = 1; k < size; k ++)
        {
            index = ((size-k*2*s)%size + size)%size + k*size;
            data[index][0] += slice[k][0];
            data[index][1] += slice[k][1];
        }
    }
}

void setSlice_Perp_Integer(const size_t s, nttw_integer *data, const nttw_integer *slice, const size_t size, const nttw_integer modulus)
{
    size_t k, index = 0;

//    if(s < size/2 && s >= 0)
//    {
//        data[0] += slice[0];
//        data[0] = MODADD2(slice[0], data[0], modulus); //data[0] += slice[0]
//        for(k = 1; k < size; k ++)
//        {
//            index = ((size-k*2*s)%size + size)%size + k*size;
//            data[index] += slice[k];
//            data[index] = MODADD2(slice[k], data[index], modulus); //data[index] += slice[k]
//        }
//    }

    if(s < size/2)
    {
        //data[0] += slice[0];
        //data[0] %= modulus;
        data[0] = MODADD(data[0], slice[0], modulus); //data[0] += slice[0]
        for(k = 1; k < size; k ++)
        {
            index = ((size-k*2*s)%size + size)%size + k*size;
            //data[index] += slice[k];
            //data[index] %= modulus;
            data[index] = MODADD(data[index], slice[k], modulus); //data[index] += slice[k]
        }
    }
}

void dyadic_oversample(nttw_big_integer *data, const int n)
{
    int j, k;
    nttw_big_integer *gcd_table;
    nttw_big_integer x, y;

    gcd_table = array_1D_big(n);

    gcd_table[0] = n + n/2;
    gcd_table[1] = 1;
    for(j = 2; j < n; j ++)
        gcd_table[j] = euclidean(j,n,&x,&y);

    for(j = 0; j < n; j ++)
        for(k = 0; k < n; k ++)
        {
            if(gcd_table[j] < gcd_table[k])
                data[j*n+k] = gcd_table[j];
            else
                data[j*n+k] = gcd_table[k];
        }

    free_array(gcd_table);
}

void dyadic_1D_filter(nttw_big_integer *data, const int n)
{
    int j;
    nttw_big_integer x, y;

    data[0] = n + n/2;
    data[1] = 1;
    for(j = 2; j < n; j ++)
        data[j] = euclidean(j,n,&x,&y);
}

void dyadic_1D_filter_Integer(nttw_big_integer *data, const size_t n, const nttw_integer modulus)
{
    size_t j;
    nttw_big_integer d, x, y;
    nttw_big_integer *gcd_table;

    gcd_table = array_1D_big(n);

    gcd_table[0] = (nttw_big_integer)(n + n/2);
    gcd_table[1] = 1;
    for(j = 2; j < n; j ++)
        gcd_table[j] = euclidean((nttw_big_integer)j,(nttw_big_integer)n,&x,&y);

    d = euclidean(gcd_table[0],modulus,&x,&y); ///Inverse of DC
    x = (x + modulus)%modulus; //Make positive
    data[0] = x;
    for(j = 1; j < n; j ++) ///Inverses
    {
        d = euclidean(gcd_table[j],modulus,&x,&y);
        x = (x + modulus)%modulus; //Make positive
        data[j] = x;
    }

    free_array(gcd_table);
}

void filter_oversampling(fftw_complex *fftSpace, nttw_big_integer *samples, const int n)
{
    int j, k;

    for(j = 0; j < n; j ++)
        for(k = 0; k < n; k ++)
        {
            if(samples[j*n+k] != 0)
            {
                fftSpace[j*n+k][0] /= (double) samples[j*n+k]; ///Correct Coefficients
                fftSpace[j*n+k][1] /= (double) samples[j*n+k]; ///Correct Coefficients
            }
            else
            {
                fftSpace[j*n+k][0] = 0;
                fftSpace[j*n+k][1] = 0;
            }
        }
}

void filter_oversampling_Integer(nttw_integer *data, nttw_integer *filter, const int n, const nttw_integer modulus)
{
    int j, k, count = 1;
    nttw_integer *tmpData, power;
    nttw_big_integer *gcd_table, d, x, y;

    gcd_table = array_1D_big(n);
    tmpData = array_1D(n*n);

    for(j = 0; j < n; j ++)
        for(k = 0; k < n; k ++)
            tmpData[j*n+k] = data[j*n+k];

    d = euclidean(filter[0],modulus,&x,&y); ///Inverse of DC
    x = (x + modulus)%modulus; //Make positive
    gcd_table[0] = x;
    for(j = 1; j < n; j ++) ///Inverses
    {
        power = 1ULL << j; //2^j
        d = euclidean(power,modulus,&x,&y);
        x = (x + modulus)%modulus; //Make positive
        gcd_table[j] = x;
    }

    for(j = 0; j < n; j ++)
        for(k = 0; k < n; k ++)
        {
            count = 1;
            if(j != 0 || k !=0) ///Ignore DC
            {
                if(filter[j*n+k] > 1) ///Ignore samples that are not multiple
                {
                    power = filter[j*n+k] >> 1ULL;
                    while(power > 1 && count < n) ///Determine power
                    {
                        count ++;
                        power >>= 1ULL;
                    }
                    tmpData[j*n+k] = (tmpData[j*n+k] * (nttw_big_integer)gcd_table[count])%modulus; ///Correct Coefficients
                }
            }
        }
    tmpData[0] = (tmpData[0] * (nttw_big_integer)gcd_table[0])%modulus; ///Correct DC

    for(j = 0; j < n; j ++)
        for(k = 0; k < n; k ++)
            data[j*n+k] = tmpData[j*n+k];

    free_array(gcd_table);
    free_array(tmpData);
}
