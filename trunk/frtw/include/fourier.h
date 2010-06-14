/**
 * FRTW FFT Library which uses the External FFTW library (which is not associated with FRTW)
 * \file fourier.h
 * \brief FFT Header/Object for the FRTW C Library.
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
#include <fftw3.h>
#include <nttw/array.h>

/**
 * \defgroup DFT_Members Discrete Fourier Transforms and Helpers
 */
//@{
/*! \fn fft(int rank, const int *sizes, fftw_complex *data, fftw_complex *fftOfData)
    \brief The Complex valued Fast Fourier Transform (FFT). Uses the FFTW C library.
*/
NTTW_DLL_SYM void fft(int rank, const int *sizes, fftw_complex *data, fftw_complex *fftOfData);

/*! \fn fft(int rank, const int *sizes, fftw_complex *fftOfData, fftw_complex *data)
    \brief The Complex valued Inverse Fast Fourier Transform (iFFT). Uses the FFTW C library.
*/
NTTW_DLL_SYM void ifft(int rank, const int *sizes, fftw_complex *fftOfData, fftw_complex *data);

/**
 * \fn fft_norm(fftw_complex *data, const size_t n)
 * \brief Normalises data of length n in the complex field.
*/
NTTW_DLL_SYM void fft_norm(fftw_complex *data, const size_t n);
//@}
