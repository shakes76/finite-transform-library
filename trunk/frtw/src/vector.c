/**
 * FRTW Vector Library
 * \file vector.c
 * \brief Multi-Dimensional Vector Source for the FRTW C Library.
 *
 * This object provides various Vector operations/methods. This wraps the array object
 * to provide vector functionality.
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
#include "vector.h"

vector vector_2D()
{
    vector tmp = arraySigned_1D(2);

    setXY(tmp,0,0);

    return tmp;
}

vector vector_2D_Set(const long x, const long y)
{
    vector tmp = arraySigned_1D(2);

    setXY(tmp,x,y);

    return tmp;
}

vector* vectorArray_2D(const size_t size)
{
    int j;
    vector *matrix;

    matrix = (vector*)malloc(size*sizeof(vector));

    for(j = 0; j < size; j ++)
        matrix[j] = vector_2D();

    return matrix;
}

long getComponent(const vector inVector, const int j)
{
    return inVector[j];
}

void setComponent(vector inVector, const int j, const long value)
{
    inVector[j] = value;
}

long getX(const vector inVector)
{
    return inVector[0];
}

void setX(vector inVector, const long value)
{
    inVector[0] = value;
}

long getY(const vector inVector)
{
    return inVector[1];
}

void setY(vector inVector, const long value)
{
    inVector[1] = value;
}

void setXY(vector inVector, const long xValue, const long yValue)
{
    inVector[0] = xValue;
    inVector[1] = yValue;
}

void free_vector(vector vec)
{
    free(vec);
}

void free_vectorArray(vector *matrix, const size_t size)
{
    int j;

    for(j = 0; j < size; j ++)
        free(matrix[j]);
    free(matrix);
}
