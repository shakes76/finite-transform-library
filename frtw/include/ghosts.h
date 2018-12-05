/**
 * FRTW Ghosts Library
 * \file ghosts.h
 * \brief Finite Ghosts Header/Object for the FRTW C Library.
 *
 * This header provides all the operations related to Finite Ghosts.
 * This includes convolution and rows-difference algorithms.
 *
 * This file is part of the FRTW Library.
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
 * \author Shekhar S. Chandra, 2008-10
*/
#ifndef GHOSTS_H_INCLUDED
#define GHOSTS_H_INCLUDED

//C NTTW Library
#include <nttw/global.h>

/**
 * \defgroup Ghost_Convolution Finite Ghost Convolution Module
 */
//@{
/**
    \fn determineMissingAngles(const size_t N, nttw_big_integer *finiteAngles, nttw_integer *perpFlags, const size_t noOfProjs, nttw_big_integer **missingAngles, nttw_big_integer **missingAngleLookup)
    \brief Determines what are the missing finite angle values and creates a lookup table also.
*/
NTTW_DLL_SYM size_t determineMissingAngles(const size_t N, nttw_big_integer *finiteAngles, nttw_integer *perpFlags, const size_t noOfProjs, nttw_big_integer **missingAngles, nttw_big_integer **missingAngleLookup);

/**
    \fn ghostEigenvalues(const size_t N, nttw_integer *root, nttw_integer *modulus, nttw_integer *proot)
    \brief Returns the eigenvalues of the ghost operators for a given N.
*/
NTTW_DLL_SYM nttw_integer* ghostEigenvalues(const size_t N, nttw_integer *root, nttw_integer *modulus, nttw_integer *proot);

/**
    \fn convolveGhosts_1D(const size_t N, nttw_integer* eigenvalues, const nttw_integer root, const nttw_integer modulus, const nttw_integer proot, nttw_big_integer *missing, nttw_big_integer *missingLookup, const size_t noOfGhosts)
    \brief Returns the ghost constructed by the 1D convolution approach out of ghost operator eigenvalues.
*/
NTTW_DLL_SYM nttw_integer* convolveGhosts_1D(const size_t N, nttw_integer* eigenvalues, const nttw_integer root, const nttw_integer modulus, const nttw_integer proot, nttw_big_integer *missing, nttw_big_integer *missingLookup, const size_t noOfGhosts);

/**
    \fn
    \brief Returns the deghosted image using the deconvolution approach.
*/
NTTW_DLL_SYM nttw_integer* deghost(const size_t N, nttw_integer* frtSpace, nttw_integer* ghost, const size_t noOfGhosts);

//@}

#endif // GHOSTS_H_INCLUDED
