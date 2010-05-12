#include "array_complex.h"

//Complex Arrays
fftw_complex* fftw_array(const size_t size)
{
    fftw_complex *data = NULL;

    data = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));

    return data;
}

void fftw_init_1D(fftw_complex *data, const size_t size, const nttw_real value)
{
    size_t j;

    for (j = 0; j < size; j ++)
    {
        data[j][0] = value;
        data[j][1] = value;
    }
}

void array_to_fftw_array(nttw_integer *source, fftw_complex *dest, const int size)
{
    int j;

    for (j = 0; j < size; j++)
        dest[j][0] = (double) source[j]; //[0] is the real part
}

void fftw_array_to_array(fftw_complex *source, nttw_integer *dest, const int size)
{
    int j;

    for (j = 0; j < size; j++)
        dest[j] = (nttw_integer) source[j][0]; //[0] is the real part
}
