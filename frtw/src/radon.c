#include <nttw/image.h>
#include <nttw/prime.h>
#include <nttw/number32.h>

#include "radon.h"
#include "array_complex.h"

nttw_integer getIsum(nttw_integer *radon, const size_t n)
{
    size_t j;
    nttw_integer Isum = 0;

    for(j = 0; j < n; j ++)
        Isum += radon[j];

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

    for (y = 0; y < p; y ++) ///Columns
        for (t = 0; t < p; t ++) ///Translate
            bins[p*p+t] += field[y*p+t];
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
            slice[j][0] /= n*n;
            slice[j][1] /= n*n;
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
            slice[j][0] /= n*n;
            slice[j][1] /= n*n;
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
            slice[j][0] /= oversampling[j];
            slice[j][1] /= oversampling[j];
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
            slice[j][0] /= oversampling[j];
            slice[j][1] /= oversampling[j];
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

void filter_oversampling(fftw_complex *fftSpace, nttw_integer *samples, const int n)
{
    int j, k;

    for(j = 0; j < n; j ++)
        for(k = 0; k < n; k ++)
        {
            fftSpace[j*n+k][0] /= (double) samples[j*n+k]; ///Correct Coefficients
            fftSpace[j*n+k][1] /= (double) samples[j*n+k]; ///Correct Coefficients
        }
}
