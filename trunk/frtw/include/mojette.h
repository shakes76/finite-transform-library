/**
 * FRTW Mojette Library
 * \file mojette.h
 * \brief Mojette Transforms Header/Object for the FRTW C Library.
 *
 * This header provides all the operations related to the Mojette Transform.
 * It includes Farey angle related routines as well as other discretized
 * versions of the Mojette Transform.
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
#ifndef MOJETTE_H_INCLUDED
#define MOJETTE_H_INCLUDED

#include <math.h>

#include "vector.h"

/**
 * \defgroup Mojette_Transform Mojette Transforms & Related
 */
//@{
/**
 * \typedef struct mojetteProjection
 * \brief Struct representing a Mojette Projection
*/
typedef struct
{
    size_t binMax;
    nttw_integer *bins;
} mojetteProjection;

/**
 * \fn mojetteArray_1D(vector *angles, const size_t mu, const size_t rows, const size_t cols)
 * \brief Allocates a Mojette Projection set (array) for the angle set provided.
*/
NTTW_DLL_SYM mojetteProjection* mojetteArray_1D(vector *angles, const size_t mu, const size_t rows, const size_t cols);
/**
 * \fn free_mojetteArray(mojetteProjection *set, const size_t mu)
 * \brief Deallocates a Mojette Projection set (array).
*/
NTTW_DLL_SYM void free_mojetteArray(mojetteProjection *set, const size_t mu);
/**
 * \fn writeMojetteSet(mojetteProjection *set, const size_t mu, const char *filename)
 * \brief Writes a Mojette Projection set (array) to file.
*/
NTTW_DLL_SYM int writeMojetteSet(mojetteProjection *set, const size_t mu, const char *filename);
/**
 * \fn outputMojetteSet(mojetteProjection *set, const size_t mu)
 * \brief Writes a Mojette Projection set (array) to stdout.
*/
NTTW_DLL_SYM int outputMojetteSet(mojetteProjection *set, const size_t mu);
/**
 * \fn readMojetteSet(mojetteProjection **set, size_t *mu, const char *filename)
 * \brief Reads a Mojette Projection set (array) from file.
*/
NTTW_DLL_SYM int readMojetteSet(mojetteProjection **set, size_t *mu, const char *filename);

/**
 * \fn isKatzCriterion(const vector *projectionSet, const int projectionsNumber, const int rows, const int cols)
 * \brief Uses the set of projections to determine if the Katz Criterion is met.
 *
 * The Katz criterion signifies whether there is sufficient information in the
 * projections to reconstruct the data from the projections.
 * Katz Criterion is met for a \f$ N\times N \f$ when
 * \f[
 * N < 1+\max\left(\sum_{i=0}^{\mu-1}|p_i|,\sum_{i=0}^{\mu-1}|q_i|\right),
 * \f]
 * where \f$ \mu \f$ is the total number of projections and \f$ p \f$ and \f$ q \f$ are
 * part of the projection vector \f$ \vec{\rho} = p + qi \f$ .
 * \return A Boolean representing if the Katz Criterion is met.
*/
NTTW_DLL_SYM int isKatzCriterion(const vector *projectionSet, const int projectionsNumber, const int rows, const int cols);

/**
 * \fn totalFareyAngles(const int n)
 * \brief Returns the upper bound for the number of Farey angles less than n.
 * Limit is approximately 3n^2/pi^2 = 0.3034 n^2. Its to be used with the compact routines.
*/
NTTW_DLL_SYM int totalFareyAngles(const int n);
/**
 * \fn nextFareyAngle(const int n, const vector angle1, const vector angle2, vector *result)
 * \brief Computes the next Farey fraction/angle from the two fractions provided.
 *
 * The computation is done using the standard method.
*/
NTTW_DLL_SYM void nextFareyAngle(const int n, const vector angle1, const vector angle2, vector *result);
/**
 * \fn nextFareyAngle_Compact(const int n, const vector angle1, const vector angle2, vector *result)
 * \brief Computes the next Farey fraction/angle from the two fractions provided that has minimal l1 norm.
 *
 * The computation is done using the standard method with a l1 norm constraint.
*/
NTTW_DLL_SYM void nextFareyAngle_Compact(const int n, const vector angle1, const vector angle2, vector *result);
/**
 * \fn allFareyAngles(const int n)
 * \brief Computes all the Farey fractions of order n and outputs to stdout.
*/
NTTW_DLL_SYM void allFareyAngles(const int n);
/**
 * \fn allFareyAngles_Compact(const int n)
 * \brief Computes all the Farey fractions of order n that has minimal l1 norm and outputs to stdout.
*/
NTTW_DLL_SYM void allFareyAngles_Compact(const int n);

/**
 * \fn finiteAngle(const vector fareyAngle, const nttw_big_integer modulus, nttw_integer *perp)
 * \brief Computes the corresponding finite angle (m or s) for the Farey fraction for the given modulus.
 *
 * The mapping was developed by Chandra et al (arXiv:1006.1965v1 [physics.med-ph]).
*/
NTTW_DLL_SYM nttw_big_integer finiteAngle(const vector fareyAngle, const nttw_big_integer modulus, nttw_integer *perp);
/**
 * \fn markFiniteAngle(const nttw_big_integer value, const nttw_integer perp, nttw_integer *mLookUpTable, nttw_integer *sLookUpTable, nttw_integer *rejected)
 * \brief Marks the lookup table if the finite angle is unique in the given finite set.
 *
 * The mapping was developed by Chandra et al (arXiv:1006.1965v1 [physics.med-ph]).
*/
NTTW_DLL_SYM void markFiniteAngle(const nttw_big_integer value, const nttw_integer perp, nttw_integer *mLookUpTable, nttw_integer *sLookUpTable, nttw_integer *rejected);
/**
 * \fn computeAngle(const size_t N, long p, long q, nttw_integer *mLookUp, nttw_integer *sLookUp, vector *farey, nttw_big_integer *finite, nttw_integer *perps, nttw_big_integer *angleCount)
 * \brief Determines the finite angle for a given Farey fraction and rejects it if it is not unique. Marks the lookup if unique.
 *
 * The mapping was developed by Chandra et al (arXiv:1006.1965v1 [physics.med-ph]).
*/
NTTW_DLL_SYM void computeAngle(const size_t N, long p, long q, nttw_integer *mLookUp, nttw_integer *sLookUp, vector *farey, nttw_big_integer *finite, nttw_integer *perps, nttw_big_integer *angleCount);

/**
 * \fn fmt_angleSet(const size_t N, const size_t P, const size_t Q, vector *farey, nttw_big_integer *finite, nttw_integer *perps)
 * \brief Determines the finite angle set for the FMT using the first encountered Farey fractions during generation of the l1 norm Farey set of order n.
 *
 * The angle set was developed by Chandra et al (arXiv:1006.1965v1 [physics.med-ph]).
*/
NTTW_DLL_SYM nttw_big_integer fmt_angleSet(const size_t N, const size_t P, const size_t Q, vector *farey, nttw_big_integer *finite, nttw_integer *perps);
/**
 * \fn fmt_angleSet_Simple(const size_t N, vector *farey, nttw_big_integer *finite, nttw_integer *perps)
 * \brief Determines the finite angle set for the FMT using the simple [1,w] set.
 *
 * The angle set was developed by Chandra et al (arXiv:1006.1965v1 [physics.med-ph]).
*/
NTTW_DLL_SYM nttw_big_integer fmt_angleSet_Simple(const size_t N, vector *farey, nttw_big_integer *finite, nttw_integer *perps);
/**
 * \fn fmt_angleSet_L1(const size_t N, const size_t P, const size_t Q, vector *farey, nttw_big_integer *finite, nttw_integer *perps)
 * \brief Determines the l1 norm finite angle set for the FMT by generating the Farey set for order n and sorting in ascending l1 norm. The shortest fractions are chosen to form the set.
 *
 * The angle set was developed by Chandra et al (arXiv:1006.1965v1 [physics.med-ph]).
*/
NTTW_DLL_SYM nttw_big_integer fmt_angleSet_L1(const size_t N, const size_t P, const size_t Q, vector *farey, nttw_big_integer *finite, nttw_integer *perps);

/**
 * \fn mt(const vector *angles, const nttw_integer *data, mojetteProjection *set, const size_t mu, const size_t rows, const size_t cols)
 * \brief The Mojette Transform (MT) of the data taken at the angles provided. The set contains the projection data.
*/
NTTW_DLL_SYM void mt(const vector *angles, const nttw_integer *data, mojetteProjection *set, const size_t mu, const size_t rows, const size_t cols);
/**
 * \fn fmt(const size_t P, const size_t Q, const nttw_integer *image, const vector *farey, const nttw_big_integer *finite, nttw_integer *perps, const size_t N)
 * \brief The Fast Mojette Transform (FMT) of the data taken at the angles provided. The FMT is simulated by taking the MT and then mapping it to FRT space.
 *
 * The FMT was developed by Chandra et al (arXiv:1006.1965v1 [physics.med-ph]).
*/
NTTW_DLL_SYM nttw_integer* fmt(const size_t P, const size_t Q, const nttw_integer *image, const vector *farey, const nttw_big_integer *finite, nttw_integer *perps, const size_t N);
/**
 * \fn fmt_noise(const size_t P, const size_t Q, const nttw_integer *image, const vector *farey, const nttw_big_integer *finite, nttw_integer *perps, const size_t N, const float SNR)
 * \brief The Fast Mojette Transform (FMT) of the data taken at the angles provided with noise at SNR. The FMT is simulated by taking the MT and then mapping it to FRT space.
 *
 * The FMT was developed by Chandra et al (arXiv:1006.1965v1 [physics.med-ph]).
*/
NTTW_DLL_SYM nttw_integer* fmt_noise(const size_t P, const size_t Q, const nttw_integer *image, const vector *farey, const nttw_big_integer *finite, nttw_integer *perps, const size_t N, const float SNR);

/**
 * \fn mt2frt(const vector *farey, const nttw_big_integer *finite, const nttw_integer *perps, const mojetteProjection *set, const size_t mu, const size_t rows, nttw_integer *frtSpace, const size_t N)
 * \brief Mojette projections are mapped to the FRT space using the finite mapping developed by Chandra et al.
*/
NTTW_DLL_SYM void mt2frt(const vector *farey, const nttw_big_integer *finite, const nttw_integer *perps, const mojetteProjection *set, const size_t mu, const size_t rows, nttw_integer *frtSpace, const size_t N);
/**
 * \fn mt2frt_noise(const vector *farey, const nttw_big_integer *finite, const nttw_integer *perps, const mojetteProjection *set, const size_t mu, const size_t rows, nttw_integer *frtSpace, const size_t N, const float SNR)
 * \brief Mojette projections are mapped to the FRT space with noise using the finite mapping developed by Chandra et al.
*/
NTTW_DLL_SYM void mt2frt_noise(const vector *farey, const nttw_big_integer *finite, const nttw_integer *perps, const mojetteProjection *set, const size_t mu, const size_t rows, nttw_integer *frtSpace, const size_t N, const float SNR);
//@}
#endif // MOJETTE_H_INCLUDED
