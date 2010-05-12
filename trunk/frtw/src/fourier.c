#include "fourier.h"

void fft(int rank, const int *sizes, fftw_complex *data, fftw_complex *fftOfData)
{
    fftw_plan p;

    p = fftw_plan_dft(rank, sizes, data, fftOfData, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    fftw_destroy_plan(p);
}

void ifft(int rank, const int *sizes, fftw_complex *fftOfData, fftw_complex *data)
{
    fftw_plan p;

    p = fftw_plan_dft(rank, sizes, fftOfData, data, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(p);

    fftw_destroy_plan(p);
}

void fft_norm(fftw_complex *data, const size_t n)
{
    size_t j;

    for(j = 0; j < n; j ++)
    {
        data[j][0] /= n;
        data[j][1] /= n;
    }
}
