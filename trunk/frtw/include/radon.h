/**
 * FRTW Radon Library
 * \file radon.h
 * \brief Radon Transforms Header/Object for the FRTW C Library.
 *
 * This header provides all the discrete operations related to the Radon Transform.
 * It includes the Discrete Radon Transform (FRT) as well as other discretized
 * versions of the Radon Transform.
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
#ifndef RADON_H_INCLUDED
#define RADON_H_INCLUDED

#include <nttw/array.h>

#include "fourier.h"

/**
 * \defgroup Radon_Transform Discrete Radon Transforms
 */
//@{
/**
 * \fn getIsum(nttw_integer *radon, const size_t n)
 * \brief Returns the Isum of the FRT space.
*/
NTTW_DLL_SYM nttw_integer getIsum(nttw_integer *radon, const size_t n);
/**
 * \fn getIsum_double(double *radon, const size_t n)
 * \brief Returns the Isum of the FRT space as a double.
*/
NTTW_DLL_SYM double getIsum_double(double *radon, const size_t n);
/**
 * \fn getIsum_Integer(nttw_integer *radon, const size_t n, const nttw_integer modulus)
 * \brief Returns the Isum (modulo the given modulus) of the FRT space.
*/
NTTW_DLL_SYM nttw_integer getIsum_Integer(nttw_integer *radon, const size_t n, const nttw_integer modulus);

/**
 * \fn drt(const nttw_integer *field, nttw_integer *bins, const int p)
 * \brief
 * The Discrete Radon Transform (DRT) of Matus & Flusser, 1993 for 2D arrays. O(n^3).
 * \author Shekhar S. Chandra, 2009
 * \param field Contains the data to be transformed. Must be p by p array.
 * \param bins It is assumed that this is already a pointer that has
 * memory is allocated for each bin. The function only does internal
 * initializing with no allocation. It must have size of p+1 by p.
 * \param p The prime which is the size of the square array.
*/
NTTW_DLL_SYM void drt(const nttw_integer *field, nttw_integer *bins, const int p);
/**
 * \fn drt_dyadic(const nttw_integer *field, nttw_integer *bins, const int n)
 * \brief
 * The Dyadic Discrete Radon Transform (DRT) of Hsung et al, 1996 for 2D arrays. O(n^3).
 * \author Shekhar S. Chandra, 2009
 * \param field Contains the data to be transformed. Must be n by n array.
 * \param bins It is assumed that this is already a pointer that has
 * memory is allocated for each bin. The function only does internal
 * initializing with no allocation. It must have size of n+n/2 by n.
 * \param n The dyadic length (32,64,...512,1024,2048,...) which is the size of the square array.
*/
NTTW_DLL_SYM void drt_dyadic(const nttw_integer *field, nttw_integer *bins, const int n);
/**
 * \fn drt_blockcopy(nttw_integer *datain, nttw_integer *dataout, const int p)
 * \brief
 * The Discrete Radon Transform (FRT) of Matus & Flusser, 1993 for 2D arrays.
 * This function uses the intelligent block copy method developed by
 * Imants Svalbe.
 * \author Andrew Kingston & Imants Svalbe, 2005.
 * 
 * Note that the definition of the m=p projection is non-standard compared to other routines 
 * in the libraries.
 * \param datain Contains the data to be transformed. Must be p by p array in
 * 1D form.
 * \param dataout It is assumed that this is already a pointer that has
 * memory is allocated for each pixel. The function only does internal
 * initializing with no allocation. It must have size of p+1 by p.
 * \param p The prime which is the size of the square array.
*/
NTTW_DLL_SYM void drt_blockcopy(nttw_integer *datain, nttw_integer *dataout, const int p);

/**
 * \fn frt(nttw_integer *field, nttw_integer *bins, const int p)
 * \brief
 * The Fast Radon Transform (FRT) of Matus & Flusser, 1993 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT). This is for prime lengths.
 * \author Shekhar S. Chandra, 2008
 * \param field Contains the data to be transformed. Must be p by p array.
 * \param bins It is assumed that this is already a pointer that has
 * memory is allocated for each bin. The function only does internal
 * initializing with no allocation. It must have size of p+1 by p.
 * \param p The prime which is the size of the square array.
*/
NTTW_DLL_SYM void frt(nttw_integer *field, nttw_integer *bins, const int p);
/**
 * \fn frt_double(double *field, double *bins, const int p)
 * \brief
 * The Fast Radon Transform (FRT) of Matus & Flusser, 1993 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT). This is for prime lengths.

 * This function is the double version. There is minimal rounding and casting between the FFT and the image.
 * \author Shekhar S. Chandra, 2009-10
 * \param field Contains the data to be transformed. Must be p by p array.
 * \param bins It is assumed that this is already a pointer that has
 * memory is allocated for each bin. The function only does internal
 * initializing with no allocation. It must have size of p+1 by p.
 * \param p The prime which is the size of the square array.
*/
NTTW_DLL_SYM void frt_double(double *field, double *bins, const int p);
/**
 * \fn frt_dyadic(nttw_integer *field, nttw_integer *bins, const int n)
 * \brief
 * The Dyadic Fast Radon Transform (FRT) of Hsung et al, 1996 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).
 * \author Shekhar S. Chandra, 2009
 * \param field Contains the data to be transformed. Must be n by n array.
 * \param bins It is assumed that this is already a pointer that has
 * memory is allocated for each bin. The function only does internal
 * initializing with no allocation. It must have size of n+n/2 by n.
 * \param n The dyadic size which is the size of the square array.
*/
NTTW_DLL_SYM void frt_dyadic(nttw_integer *field, nttw_integer *bins, const int n);
/**
 * \fn frt_dyadic_double(double *field, double *bins, const int n)
 * \brief
 * The Dyadic Fast Radon Transform (FRT) of Hsung et al, 1996 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).

 * This function is the double version. There is minimal rounding and casting between the FFT and the image.
 * \author Shekhar S. Chandra, 2009
 * \param field Contains the data to be transformed. Must be n by n array.
 * \param bins It is assumed that this is already a pointer that has
 * memory is allocated for each bin. The function only does internal
 * initializing with no allocation. It must have size of n+n/2 by n.
 * \param n The dyadic size which is the size of the square array.
*/
NTTW_DLL_SYM void frt_dyadic_double(double *field, double *bins, const int n);
/**
 * \fn frt_dyadic_signed(long *field, long *bins, const int n)
 * \brief
 * The Dyadic Fast Radon Transform (FRT) of Hsung et al, 1996 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).
 * \author Shekhar S. Chandra, 2009
 * \param field Contains the signed data to be transformed. Must be n by n array.
 * \param bins It is assumed that this is already a pointer that has
 * memory is allocated for each bin. The function only does internal
 * initializing with no allocation. It must have size of n+n/2 by n.
 * \param n The dyadic size which is the size of the square array.
*/
NTTW_DLL_SYM void frt_dyadic_signed(long *field, long *bins, const int n);

/**
 * \fn idrt(nttw_integer *bins, nttw_integer *result, const int p)
 * \brief
 * The Inverse Discrete Radon Transform (iDRT) of Matus & Flusser, 1993 for 2D arrays. O(n^3).
 * \author Shekhar S. Chandra, 2008
 * \param bins Contains the data to be inverse transformed and should be of
 * size p+1 by p.
 * \param result It is assumed that this is already a pointer that has
 * memory is allocated for each pixel. The function only does internal
 * initializing with no allocation. It must have size of p by p.
 * \param p The prime which is the size of the square array.
*/
NTTW_DLL_SYM nttw_integer idrt(nttw_integer *bins, nttw_integer *result, const int p);
/**
 * \fn idrt_blockcopy(nttw_integer *datain, nttw_integer *dataout, const int p)
 * \brief
 * The Inverse Discrete Radon Transform (iDRT) of Matus & Flusser, 1993 for 2D arrays.
 * This function uses the intelligent block copy method developed by
 * Imants Svalbe. Note m=0 and m=p projections are interchanged in comparision to other DRT functions!
 * \author Andrew Kingston & Imants Svalbe, 2005.
 *
 * Note that the definition of the m=p projection is non-standard compared to other routines 
 * in the libraries.
 * \param datain Contains the data to be inverse transformed. Must be p+1 by p
 * array in 1D form.
 * \param dataout It is assumed that this is already a pointer that has
 * memory is allocated for each pixel. The function only does internal
 * initializing with no allocation. It must have size of p by p.
 * \param p The prime which is the size of the square array.
*/
NTTW_DLL_SYM nttw_integer idrt_blockcopy(nttw_integer *datain, nttw_integer *dataout, const int p);

/**
 * \fn ifrt(nttw_integer *bins, nttw_integer *result, const int p, const int norm)
 * \brief
 * The Inverse Fast Radon Transform (iFRT) of Matus & Flusser, 1993 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).
 * \author Shekhar S. Chandra, 2008
 * \param bins Contains the data to be inverse transformed and should be of
 * size p+1 by p.
 * \param result It is assumed that this is already a pointer that has
 * memory is allocated for each pixel. The function only does internal
 * initializing with no allocation. It must have size of p by p.
 * \param p The prime which is the size of the square array.
 * \param norm Normalise the result?
*/
NTTW_DLL_SYM nttw_integer ifrt(nttw_integer *bins, nttw_integer *result, const int p, const int norm);
/**
 * \fn ifrt_double(double *bins, double *result, const int p, const int norm)
 * \brief
 * The Inverse Fast Radon Transform (iFRT) of Matus & Flusser, 1993 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).

 * This function is the double version. There is minimal rounding and casting between the FFT and the image.
 * \author Shekhar S. Chandra, 2009-10
 * \param bins Contains the data to be inverse transformed and should be of
 * size p+1 by p.
 * \param result It is assumed that this is already a pointer that has
 * memory is allocated for each pixel. The function only does internal
 * initializing with no allocation. It must have size of p by p.
 * \param p The prime which is the size of the square array.
 * \param norm Normalise the result?
*/
NTTW_DLL_SYM double ifrt_double(double *bins, double *result, const int p, const int norm);
/**
 * \fn ifrt_signed(nttw_integer *bins, long *result, const int p, const int norm)
 * \brief
 * The Inverse Fast Radon Transform (iFRT) of Matus & Flusser, 1993 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).

 * This is the signed version, ideally suited for signed data such as noise.
 * \author Shekhar S. Chandra, 2009-10
 * \param bins Contains the data to be inverse transformed and should be of
 * size p+1 by p.
 * \param result It is assumed that this is already a pointer that has
 * memory is allocated for each pixel. The function only does internal
 * initializing with no allocation. It must have size of p by p.
 * \param p The prime which is the size of the square array.
 * \param norm Normalise the result?
*/
NTTW_DLL_SYM nttw_integer ifrt_signed(nttw_integer *bins, long *result, const int p, const int norm);
/**
 * \fn ifrt_dyadic(nttw_integer *bins, nttw_integer *result, const int n, const int norm)
 * \brief
 * The Inverse Dyadic Fast Radon Transform (FRT) of Hsung et al, 1996 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).
 * \author Shekhar S. Chandra, 2009
 * \param bins Contains the data to be inverse transformed. Must be n by n+n/2 array.
 * \param result The recontructed image. It must have size of n by n.
 * \param n The dyadic size which is the size of the square array.
*/
NTTW_DLL_SYM nttw_integer ifrt_dyadic(nttw_integer *bins, nttw_integer *result, const int n, const int norm);
/**
 * \fn ifrt_dyadic_double(double *bins, double *result, const int n, const int norm)
 * \brief
 * The Inverse Dyadic Fast Radon Transform (FRT) of Hsung et al, 1996 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).

 * This function is the double version. There is minimal rounding and casting between the FFT and the image.
 * \author Shekhar S. Chandra, 2009-10
 * \param bins Contains the data to be inverse transformed. Must be n by n+n/2 array.
 * \param result The recontructed image. It must have size of n by n.
 * \param n The dyadic size which is the size of the square array.
*/
NTTW_DLL_SYM double ifrt_dyadic_double(double *bins, double *result, const int n, const int norm);
/**
 * \fn ifrt_dyadic_signed(nttw_integer *bins, long *result, const int n, const int norm)
 * \brief
 * The Inverse Dyadic Fast Radon Transform (FRT) of Hsung et al, 1996 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).
 * \author Shekhar S. Chandra, 2009
 * \param bins Contains the data to be inverse transformed. Must be n by n+n/2 array.
 * \param result The signed recontructed image. It must have size of n by n.
 * \param n The dyadic size which is the size of the square array.
*/
NTTW_DLL_SYM nttw_integer ifrt_dyadic_signed(nttw_integer *bins, long *result, const int n, const int norm);
/**
 * \fn ifrt_dyadic_signed2(long *bins, long *result, const int n, const int norm)
 * \brief
 * The Inverse Dyadic Fast Radon Transform (FRT) of Hsung et al, 1996 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT). Complete signed array version.
 * \author Shekhar S. Chandra, 2009
 * \param bins Contains the signed data to be inverse transformed. Must be n by n+n/2 array.
 * \param result The signed recontructed image. It must have size of n by n.
 * \param n The dyadic size which is the size of the square array.
*/
NTTW_DLL_SYM nttw_integer ifrt_dyadic_signed2(long *bins, long *result, const int n, const int norm);

/**
 * \fn getSlice(int m, fftw_complex *data, fftw_complex *slice, const int size)
 * \brief Extracts the slice at slope m within the FFT space given by data.
*/
NTTW_DLL_SYM void getSlice(int m, fftw_complex *data, fftw_complex *slice, const int size);
/**
 * \fn getSlice_Perp(int s, fftw_complex *data, fftw_complex *slice, const int size)
 * \brief Extracts the slice at (perp) slope s within the FFT space given by data.
*/
NTTW_DLL_SYM void getSlice_Perp(int s, fftw_complex *data, fftw_complex *slice, const int size);

/**
 * \fn setSlice(int m, fftw_complex *data, fftw_complex *slice, const int size)
 * \brief Sets the slice at slope m within the FFT space given by data.
*/
NTTW_DLL_SYM void setSlice(int m, fftw_complex *data, fftw_complex *slice, const int size);
/**
 * \fn setSlice_Perp(int s, fftw_complex *data, fftw_complex *slice, const int size)
 * \brief Sets the slice at (perp) slope s within the FFT space given by data.
*/
NTTW_DLL_SYM void setSlice_Perp(int s, fftw_complex *data, fftw_complex *slice, const int size);

/**
 * \fn dyadic_oversample(nttw_big_integer *data, const int n)
 * \brief Constructs the oversample filter required to make FFT space isotropic for dyadic sizes.
*/
NTTW_DLL_SYM void dyadic_oversample(nttw_big_integer *data, const int n);
/**
 * \fn dyadic_1D_filter(nttw_big_integer *data, const int n)
 * \brief Constructs the 1D oversample filter required to make FFT space isotropic for dyadic sizes.
*/
NTTW_DLL_SYM void dyadic_1D_filter(nttw_big_integer *data, const int n);
/**
 * \fn filter_oversampling(fftw_complex *fftSpace, nttw_big_integer *samples, const int n)
 * \brief Filters the oversampling within the FFT space for dyadic sizes given by the filter.

 * Use the dyadic_oversample() function to construct the filter.
*/
NTTW_DLL_SYM void filter_oversampling(fftw_complex *fftSpace, nttw_big_integer *samples, const int n);
//@}
#endif // RADON_H_INCLUDED
