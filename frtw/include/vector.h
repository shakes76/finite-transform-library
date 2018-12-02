/**
 * FRTW Vector Library
 * \file vector.h
 * \brief Multi-Dimensional Vector Header/Object for the FRTW C Library.
 *
 * This object provides various Vector operations/methods. This wraps the array object
 * to provide vector functionality.
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
#ifndef VECTOR_H_INCLUDED
#define VECTOR_H_INCLUDED

#include <nttw/array.h>

/**
 * \typedef vector
 * \brief This designates the data type that is used to represent a mathematical vector.
*/
typedef long* vector;

/**
 * \defgroup Vector Vector Structures and Operations.
*/
//@{
/**
 * \fn vector vector_2D()
 * \brief Constructs a 2 dimensional vector and returns it.
 *
 * Use the free() operation from <stdlib.h> to delete it.
*/
NTTW_DLL_SYM vector vector_2D();
/**
 * \fn vector vector_2D_Set(const long x, const long y)
 * \brief Constructs a 2 dimensional vector and returns it.
 *
 * Use the free() operation from <stdlib.h> to delete it.
*/
NTTW_DLL_SYM vector vector_2D_Set(const long x, const long y);
/**
 * \fn vector vectorArray_2D(const size_t size)
 * \brief Constructs a 2 dimensional vector and returns it.
 *
 * Use the free() operation from <stdlib.h> to delete it.
*/
NTTW_DLL_SYM vector* vectorArray_2D(const size_t size);
/**
 * \fn getComponent(const vector inVector, const int j)
 * \brief Returns the \f$ x_j \f$ coordinate/component of the vector.
*/
NTTW_DLL_SYM long getComponent(const vector inVector, const int j);
/**
 * \fn setComponent(vector inVector, const int j, const long value)
 * \brief Sets component j of inVector to value provided.
*/
NTTW_DLL_SYM void setComponent(vector inVector, const int j, const long value);
/**
 * \fn getX(const vector inVector)
 * \brief Returns the x component of the vector inVector.
*/
NTTW_DLL_SYM long getX(const vector inVector);
/**
 * \fn setX(vector inVector, const long value)
 * \brief Sets the x component of the vector inVector.
*/
NTTW_DLL_SYM void setX(vector inVector, const long value);
/**
 * \fn getY(const vector inVector)
 * \brief Returns the y component of the vector inVector.
*/
NTTW_DLL_SYM long getY(const vector inVector);
/**
 * \fn setY(vector inVector, const long value)
 * \brief Sets the y component of the vector inVector.
*/
NTTW_DLL_SYM void setY(vector inVector, const long value);
/**
 * \fn setXY(vector inVector, const long xValue, const long yValue)
 * \brief Sets the x and y components of the vector inVector.
*/
NTTW_DLL_SYM void setXY(vector inVector, const long xValue, const long yValue);
/**
 * \fn free_vector(vector vec)
 * \brief Free/Deallocates memory of the vector vec.
*/
NTTW_DLL_SYM void free_vector(vector vec);
/**
 * \fn free_vectorArray(vector *matrix, const size_t size)
 * \brief Free/Deallocates memory used by vector array provided.
*/
NTTW_DLL_SYM void free_vectorArray(vector *matrix, const size_t size);
//@}
#endif // VECTOR_H_INCLUDED
