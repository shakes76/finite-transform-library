/**
 * FRTW FFT Library which uses the External FFTW library
 * \file fourier.h
 * \brief FFT Header/Object for the FRTW C Library.
 *
 * This wraps all the FFTW FFT 1D/2D members.
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
