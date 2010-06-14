/**
 * FRTW FFT Library which uses the External FFTW library (which is not associated with FRTW)
 * \file fourier.c
 * \brief FFT Source for the FRTW C Library.
 *
 * This wraps all the FFTW FFT 1D/2D members.
 *
 * This file is part of FRTW Library.
 *
 * FRTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FRTW is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FRTW. If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2008-9
*/
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
