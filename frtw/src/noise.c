/**
 * FRTW Noise Library
 * \file noise.c
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
 * \author Shekhar S. Chandra, 2008-9
*/
#include <math.h>

#include "noise.h"

double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed
 * between 0.0 and 1.0.
 * ----------------------------------------------------------------
 */
{
  const long Q = RNG_MODULUS / MULTIPLIER;
  const long R = RNG_MODULUS % MULTIPLIER;
        long t;

  t = MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
  if (t > 0)
    seed[stream] = t;
  else
    seed[stream] = t + RNG_MODULUS;
  return ((double) seed[stream] / RNG_MODULUS);
}

double Exponential(double m)
/* =========================================================
 * Returns an exponentially distributed positive real number.
 * NOTE: use m > 0.0
 * =========================================================
 */
{
  return (-m * log(1.0 - Random()));
}

long Poisson(double m)
/* ==================================================
 * Returns a Poisson distributed non-negative integer.
 * NOTE: use m > 0
 * ==================================================
 */
{
  double t = 0.0;
  long   x = 0;

  while (t < m) {
    t += Exponential(1.0);
    x++;
  }
  return (x - 1);
}

double Normal(double m, double s)
/* ========================================================================
 * Returns a normal (Gaussian) distributed real number.
 * NOTE: use s > 0.0
 *
 * Uses a very accurate approximation of the normal idf due to Odeh & Evans,
 * J. Applied Statistics, 1974, vol 23, pp 96-97.
 * ========================================================================
 */
{
  const double p0 = 0.322232431088;     const double q0 = 0.099348462606;
  const double p1 = 1.0;                const double q1 = 0.588581570495;
  const double p2 = 0.342242088547;     const double q2 = 0.531103462366;
  const double p3 = 0.204231210245e-1;  const double q3 = 0.103537752850;
  const double p4 = 0.453642210148e-4;  const double q4 = 0.385607006340e-2;
  double u, t, p, q, z;

  u   = Random();
  if (u < 0.5)
    t = sqrt(-2.0 * log(u));
  else
    t = sqrt(-2.0 * log(1.0 - u));
  p   = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4)));
  q   = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));
  if (u < 0.5)
    z = (p / q) - t;
  else
    z = t - (p / q);
  return (m + s * z);
}

double mse(long *data1, long *data2, const size_t rows, const size_t cols)
{
    double MSE = 0.0;
    long tmpValue = 0;
    size_t j, k;
    const size_t size = rows*cols;

    for(j = 0; j < rows; j ++)
        for(k = 0; k < cols; k ++)
        {
            tmpValue = data1[j*cols+k] - data2[j*cols+k];
            MSE += tmpValue*tmpValue;
        }
    MSE /= size;

    //RMSE = sqrt(MSE);
    //PSNR = 10.0*log10(255*255/MSE);
    return MSE;
}
