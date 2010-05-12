/**
 * FRTW Radon Library
 * \file radon.h
 * \brief Radon Transforms Header/Object for the FRTW C Library.
 *
 * This header provides all the discrete operations related to the Radon Transform.
 * It includes the Discrete Radon Transform (FRT) as well as other discretized
 * versions of the Radon Transform.
 *
 * \author Shekhar S. Chandra, 2009
*/
#ifndef RADON_H_INCLUDED
#define RADON_H_INCLUDED

#include <nttw/array.h>

#include "fourier.h"

/**
 * \defgroup Radon_Transform Discrete Radon Transforms
 */
//@{
NTTW_DLL_SYM nttw_integer getIsum(nttw_integer *radon, const size_t n);
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
 * \fn ifrt_dyadic(nttw_integer *field, nttw_integer *bins, const int n)
 * \brief
 * The Inverse Dyadic Fast Radon Transform (FRT) of Hsung et al, 1996 for 2D arrays. Fourier Slice Theorem version.
 * Complexity is linear logarithmic in time (same as a 2D FFT).
 * \author Shekhar S. Chandra, 2009
 * \param bins Contains the data to be inverse transformed. Must be n by n+n/2 array.
 * \param result The recontructed image. It must have size of n by n.
 * \param n The dyadic size which is the size of the square array.
*/
NTTW_DLL_SYM nttw_integer ifrt_dyadic(nttw_integer *bins, nttw_integer *result, const int n, const int norm);

NTTW_DLL_SYM void getSlice(int m, fftw_complex *data, fftw_complex *slice, const int size);
NTTW_DLL_SYM void getSlice_Perp(int s, fftw_complex *data, fftw_complex *slice, const int size);
NTTW_DLL_SYM void setSlice(int m, fftw_complex *data, fftw_complex *slice, const int size);
NTTW_DLL_SYM void setSlice_Perp(int s, fftw_complex *data, fftw_complex *slice, const int size);
NTTW_DLL_SYM void dyadic_oversample(nttw_big_integer *data, const int n);
NTTW_DLL_SYM void dyadic_1D_filter(nttw_big_integer *data, const int n);
NTTW_DLL_SYM void filter_oversampling(fftw_complex *fftSpace, nttw_integer *samples, const int n);

//@}
#endif // RADON_H_INCLUDED
