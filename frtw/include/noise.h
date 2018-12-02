/**
 * FRTW Noise Library
 * \file radon.h
 * \brief Random Number Generators (RNGs) for Noise Header/Object for the FRTW C Library.
 *
 * This header provides all the RNGs for producing noise of various types.
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
 * \author Steve Park & Dave Geyer, 1998
*/
/**
Copyright Information (python-cvxopt.copyright):
RNGS Random Number Generation - Multiple Streams (Sep. 22, 1998) by
	Steve Park & Dave Geyer.  See www.cs.wm.edu/~va/software/park/park.html.
	Software is in the public domain as was found out in an email correspondence
	with Virginia Torczon <va@cs.wm.edu> and Lawrence M. Leemis
	<leemis@MATH.WM.EDU> (Professor Steve Park passed away):

	Dear Soeren,

	Glad to have been of help.  And I'm glad to learn from Larry that the
	software is, indeed, in the public domain.

	Virginia
*/
#ifndef NOISE_H_INCLUDED
#define NOISE_H_INCLUDED

#include <stdlib.h>
#include <math.h>
#include <nttw/global.h>

#define RNG_MODULUS    2147483647 /* DON'T CHANGE THIS VALUE                  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define STREAMS    256        /* # of streams, DON'T CHANGE THIS VALUE    */
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */
#define RMSE(x) ( sqrt(x) )
#define PSNR(x) ( 10.0*log10(255*255/(x)) )

static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */

/**
 * \fn Random(void)
 * \brief Returns 1 with probability p or 0 with probability 1 - p. Use n > 0 and 0.0 < p < 1.0
*/
NTTW_DLL_SYM double Random(void);
/**
 * \fn Exponential(double m)
 * \brief Returns an exponentially distributed positive real number. Use m > 0.0
*/
NTTW_DLL_SYM double Exponential(double m);
/**
 * \fn Poisson(double m)
 * \brief Returns a Poisson distributed non-negative integer. Use m > 0.0
*/
NTTW_DLL_SYM long Poisson(double m);
/**
 * \fn Normal(double m, double s)
 * \brief Returns a normal (Gaussian) distributed real number. Use m > 0.0 and s > 0.0
 *
 * Uses a very accurate approximation of the normal idf due to Odeh & Evans,
 * J. Applied Statistics, 1974, vol 23, pp 96-97.
*/
NTTW_DLL_SYM double Normal(double m, double s);
/**
 * \fn mse(long *data1, long *data2, const size_t rows, const size_t cols)
 * \brief Computes and returns the Mean Squared Error of two arrays.
*/
NTTW_DLL_SYM double mse(long *data1, long *data2, const size_t rows, const size_t cols);

#endif // NOISE_H_INCLUDED
