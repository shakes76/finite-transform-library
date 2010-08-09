/**
 * FRTW Mojette Library
 * \file mojette.c
 * \brief Mojette Transforms Source for the FRTW C Library.
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
#include <stdio.h>

#include <nttw/prime.h>
#include <nttw/image.h>

#include "mojette.h"
#include "noise.h"

mojetteProjection* mojetteArray_1D(vector *angles, const size_t mu, const size_t rows, const size_t cols)
{
    size_t j;
    mojetteProjection *projSet;

    projSet = (mojetteProjection *) malloc( mu*sizeof(mojetteProjection) );

    for(j = 0; j < mu; j ++)
    {
        projSet[j].binMax = 1 + (rows-1)*abs(getX(angles[j])) + (cols-1)*abs(getY(angles[j]));
        projSet[j].bins = array_1D(projSet[j].binMax);
        init_1D(projSet[j].bins, projSet[j].binMax, 0);
    }

    return projSet;
}

void free_mojetteArray(mojetteProjection *set, const size_t mu)
{
    size_t j;

    for(j = 0; j < mu; j ++)
        free_array(set[j].bins);

    free(set); ///\todo free_array instead?
}

int writeMojetteSet(mojetteProjection *set, const size_t mu, const char *filename)
{
    size_t j, k;
    FILE *outFile;
    int binary = FALSE;

    if( !openFile_Write(filename,&outFile,binary) )
        return FALSE;

    fprintf(outFile,"%lu\n",mu);
    for(j = 0; j < mu; j ++)
    {
        fprintf(outFile,"%lu\n",set[j].binMax);
        for(k = 0; k < set[j].binMax; k ++)
            fprintf(outFile,FORMAT_OUTPUT_STR,set[j].bins[k]);
        fprintf(outFile,"\n");
    }

    fclose(outFile);
    return TRUE;
}

int outputMojetteSet(mojetteProjection *set, const size_t mu)
{
    size_t j, k;

    fprintf(stdout,"%lu\n",mu);
    for(j = 0; j < mu; j ++)
    {
        fprintf(stdout,"%lu\n",set[j].binMax);
        for(k = 0; k < set[j].binMax; k ++)
            fprintf(stdout,FORMAT_OUTPUT_STR,set[j].bins[k]);
        fprintf(stdout,"\n");
    }

    return TRUE;
}

int readMojetteSet(mojetteProjection **set, size_t *mu, const char *filename)
{
    size_t j, k;
    FILE *inFile;
    int binary = FALSE, success = FALSE;

    if( !openFile_Read(filename, &inFile, binary) )
        return FALSE;

    success = fscanf(inFile, "%lu\n", mu);

    *set = (mojetteProjection *) malloc( (*mu)*sizeof(mojetteProjection) );

    for(j = 0; j < *mu; j ++)
    {
        success = fscanf(inFile, "%lu\n", &( (*set)[j].binMax ));
        (*set)[j].bins = array_1D((*set)[j].binMax);
        for(k = 0; k < (*set)[j].binMax; k ++)
            success = fscanf(inFile, FORMAT_INPUT_STR, &( (*set)[j].bins[k] ));
        success = fscanf(inFile, "\n");

        if(feof(inFile))
            break;
    }

    fclose(inFile);
    return TRUE;
}

int isKatzCriterion(const vector *projectionSet, const int projectionsNumber, const int rows, const int cols)
{
    int j, p, q;
    int sumOfP = 0, sumOfQ = 0;

    for(j = 0; j < projectionsNumber; j ++)
    {
        p = (int)getY(projectionSet[j]);
        q = (int)getX(projectionSet[j]);
        sumOfP += abs(p);
        sumOfQ += abs(q);
    }

    //printf("Sum of p's: %d, Sum of q's: %d\n",sumOfP,sumOfQ);

    if( (sumOfP + 1) > rows || (sumOfP + 1) > cols )
        return TRUE;
    else if( (sumOfQ + 1) > rows || (sumOfQ + 1) > cols )
        return TRUE;
    else
        return FALSE;
}

int totalFareyAngles(const int n)
{
    int N = (int) (0.304*n*n + 0.5);

    return N;
}

void nextFareyAngle(const int n, const vector angle1, const vector angle2, vector *result)
{
    long p1, q1, p2, q2, p3, q3;

    p1 = getY(angle1);
    q1 = getX(angle1);
    p2 = getY(angle2);
    q2 = getX(angle2);

    p3 = (long) floor( (q1+n)/(double)q2 )*p2 - p1;
    q3 = (long) floor( (q1+n)/(double)q2 )*q2 - q1;

    setY(*result,p3);
    setX(*result,q3);
}

void nextFareyAngle_Compact(const int n, const vector angle1, const vector angle2, vector *result)
{
    long p1, q1, p2, q2, p3, q3;

    p1 = getY(angle1);
    q1 = getX(angle1);
    p2 = getY(angle2);
    q2 = getX(angle2);

    p3 = (long) floor( (q1+p1+n)/(double)(q2+p2) )*p2 - p1;
    q3 = (long) floor( (q1+p1+n)/(double)(q2+p2) )*q2 - q1;

    setY(*result,p3);
    setX(*result,q3);
}

void allFareyAngles(const int n)
{
    vector angle1, angle2, nextAngle;

    angle1 = vector_2D();
    angle2 = vector_2D();
    nextAngle = vector_2D();

    setY(angle1,0);
    setX(angle1,1); // 0/1
    setY(angle2,1);
    setX(angle2,n); // 1/n

    fprintf(stderr,"Farey(%i): \n",n);
    while( !(getY(nextAngle) == 1 && getX(nextAngle) == 1) ) // 1/1
    {
        nextFareyAngle(n, angle1, angle2, &nextAngle);
        fprintf(stderr,"%lli/%lli,\t",getY(nextAngle),getX(nextAngle));

        setY(angle1, getY(angle2));
        setX(angle1, getX(angle2));
        setY(angle2, getY(nextAngle));
        setX(angle2, getX(nextAngle));
    }
    fprintf(stderr,"\n");
}

void allFareyAngles_Compact(const int n)
{
    vector angle1, angle2, nextAngle;

    angle1 = vector_2D();
    angle2 = vector_2D();
    nextAngle = vector_2D();

    setY(angle1,0);
    setX(angle1,1); // 0/1
    setY(angle2,1);
    setX(angle2,n); // 1/n

    fprintf(stderr,"Farey(%i) (Compact): \n",n);
    while( !(getY(nextAngle) == 1 && getX(nextAngle) == 1) ) // 1/1
    {
        nextFareyAngle_Compact(n, angle1, angle2, &nextAngle);
        fprintf(stderr,"%lli/%lli,\t",getY(nextAngle),getX(nextAngle));

        setY(angle1, getY(angle2));
        setX(angle1, getX(angle2));
        setY(angle2, getY(nextAngle));
        setX(angle2, getX(nextAngle));
    }
    fprintf(stderr,"\n");
}

nttw_big_integer finiteAngle(const vector fareyAngle, const nttw_big_integer modulus, nttw_integer *perp)
{
    nttw_big_integer pInv, qInv, p, q, sValue, mValue;

    if(getY(fareyAngle)%2 == 0)
    {
        pInv = minverse_big(getX(fareyAngle),modulus);
        q = getY(fareyAngle);
        sValue = (q*pInv)%modulus;
        sValue /= 2;

        *perp = TRUE;

        return sValue;
    }
    else
    {
        qInv = minverse_big(getY(fareyAngle),modulus);

        p = getX(fareyAngle);
        mValue = (p*qInv)%modulus;

        *perp = FALSE;

        return mValue;
    }
}

void markFiniteAngle(const nttw_big_integer value, const nttw_integer perp, nttw_integer *mLookUpTable, nttw_integer *sLookUpTable, nttw_integer *rejected)
{
    if(perp)
    {
        if(!sLookUpTable[value])
        {
            sLookUpTable[value] = TRUE;
            *rejected = TRUE;
        }
    }
    else
    {
        if(!mLookUpTable[value])
        {
            mLookUpTable[value] = TRUE;
            *rejected = TRUE;
        }
    }
}

void computeAngle(const size_t N, long p, long q, nttw_integer *mLookUp, nttw_integer *sLookUp, vector *farey, nttw_big_integer *finite, nttw_integer *perps, nttw_big_integer *angleCount)
{
    vector angle;
    nttw_integer found = FALSE, perp = FALSE;
    nttw_big_integer value;

    angle = vector_2D();

    ///Compute Angle by setting vector
    setY(angle,q);
    setX(angle,p); // q/p

    value = finiteAngle(angle, N, &perp); ///Determine the finite angle analytically
    markFiniteAngle(value,perp,mLookUp,sLookUp,&found); ///Check and Mark lookup tables

    if(found)
    {
        setXY(farey[*angleCount],getX(angle),getY(angle));
        perps[*angleCount] = perp;
        finite[*angleCount] = value;
        (*angleCount) ++;
    }
}

nttw_big_integer fmt_angleSet(const size_t N, const size_t P, const size_t Q, vector *farey, nttw_big_integer *finite, nttw_integer *perps)
{
    const size_t size = N+N/2;
    nttw_integer *mLookUp, *sLookUp;
    nttw_big_integer count = 0, total = 0;
    vector angle1, angle2, nextAngle;
    size_t n;

    if(P > Q)
        n = P;
    else
        n = Q;

    mLookUp = array_1D(N);
    sLookUp = array_1D(N/2);
    angle1 = vector_2D();
    angle2 = vector_2D();
    nextAngle = vector_2D();

    init_1D(mLookUp,N,FALSE);
    init_1D(sLookUp,N/2,FALSE);

    ///Compute Start and End angles, order important
    computeAngle(N, 0, 1, mLookUp, sLookUp, farey, finite, perps, &count);
    computeAngle(N, 1, (long) n, mLookUp, sLookUp, farey, finite, perps, &count); ///\todo Typecast safe?
    computeAngle(N, 1, 0, mLookUp, sLookUp, farey, finite, perps, &count);
    setXY(angle1,getX(farey[count-1]),getY(farey[count-1]));
    computeAngle(N, (long) n, 1, mLookUp, sLookUp, farey, finite, perps, &count); ///\todo Typecast safe?
    setXY(angle2,getX(farey[count-1]),getY(farey[count-1]));
    total = count;

    while( !(getY(nextAngle) == 1 && getX(nextAngle) == 1) && count < size) // 1/1
    {
        nextFareyAngle_Compact(n, angle1, angle2, &nextAngle);

        ///First Octant
        computeAngle(N, getX(nextAngle), getY(nextAngle), mLookUp, sLookUp, farey, finite, perps, &count);
        ///Second Octant
        computeAngle(N, getY(nextAngle), getX(nextAngle), mLookUp, sLookUp, farey, finite, perps, &count);
        ///Third Octant
        computeAngle(N, getX(nextAngle), -getY(nextAngle), mLookUp, sLookUp, farey, finite, perps, &count);
        ///Forth Octant
        computeAngle(N, getY(nextAngle), -getX(nextAngle), mLookUp, sLookUp, farey, finite, perps, &count);

        setXY(angle1,getX(angle2),getY(angle2));
        setXY(angle2,getX(nextAngle),getY(nextAngle));
        total ++;
    }

    ///Deallocate
    free_array(mLookUp);
    free_array(sLookUp);
    free_vector(angle1);
    free_vector(angle2);
    free_vector(nextAngle);

    return total;
}

nttw_big_integer fmt_angleSet_Simple(const size_t N, vector *farey, nttw_big_integer *finite, nttw_integer *perps)
{
    const size_t size = N+N/2;
    nttw_integer *mLookUp, *sLookUp;
    nttw_big_integer count = 0;
    size_t n = 0;

    mLookUp = array_1D(N);
    sLookUp = array_1D(N/2);

    init_1D(mLookUp,N,FALSE);
    init_1D(sLookUp,N/2,FALSE);

    ///Compute Start and End angles, order important
    while(count < size)
    {
        computeAngle(N, 1, (long) n, mLookUp, sLookUp, farey, finite, perps, &count); ///\todo Typecast safe?
        computeAngle(N, (long) n, 1, mLookUp, sLookUp, farey, finite, perps, &count); ///\todo Typecast safe?
        computeAngle(N, (long) n, -1, mLookUp, sLookUp, farey, finite, perps, &count); ///\todo Typecast safe?
        computeAngle(N, 1, - (long) n, mLookUp, sLookUp, farey, finite, perps, &count); ///\todo Typecast safe?

        n ++;
    }

    ///Deallocate
    free_array(mLookUp);
    free_array(sLookUp);

    return count;
}

nttw_big_integer fmt_angleSet_L1(const size_t N, const size_t P, const size_t Q, vector *farey, nttw_big_integer *finite, nttw_integer *perps)
{
    const size_t size = N+N/2;
    nttw_integer *mLookUp, *sLookUp, *L1;
    nttw_big_integer count = 0, total = 0, index = 0, min = 0;
    vector angle1, angle2, nextAngle, *anglesRaw, *anglesSorted;
    size_t j, k, n, angleSize = 0;

    if(P > Q)
        n = P;
    else
        n = Q;

    angleSize = 6*N*N/9; //Half of 12*N*N/pi^2, since angles in triangle

    mLookUp = array_1D(N);
    sLookUp = array_1D(N/2);
    L1 = array_1D(angleSize);
    anglesRaw = vectorArray_2D(angleSize);
    anglesSorted = vectorArray_2D(angleSize);
    angle1 = vector_2D();
    angle2 = vector_2D();
    nextAngle = vector_2D();

    init_1D(mLookUp,N,FALSE);
    init_1D(sLookUp,N/2,FALSE);

    ///Compute Start and End angles, order important
    setXY(anglesRaw[count],0,1);
    L1[count] = 1;
    count ++;

    setXY(anglesRaw[count],1,(long) n);
    L1[count] = n+1;
    count ++;

    setXY(anglesRaw[count],1,0);
    L1[count] = 1;
    setXY(angle1, 1, 0);
    count ++;

    setXY(anglesRaw[count],(long) n,1);
    L1[count] = n+1;
    setXY(angle2, (long) n, 1);
    count ++;

    total = count; //Compute angle changes count

    while( !(getY(nextAngle) == 1 && getX(nextAngle) == 1)) // 1/1
    {
        nextFareyAngle_Compact(N, angle1, angle2, &nextAngle);
        //fprintf(stderr,"%li/%li,\t",getY(nextAngle),getX(nextAngle));

        setY(angle1, getY(angle2));
        setX(angle1, getX(angle2));
        setY(angle2, getY(nextAngle));
        setX(angle2, getX(nextAngle));

        setXY(anglesRaw[count],getX(nextAngle),getY(nextAngle));
        L1[count] = abs(getX(nextAngle)) + abs(getY(nextAngle));
        count ++;
    }
    fprintf(stderr,"Done l1 set.\n");

    ///Sort Angles by L1, \todo Use qsort here
    index = 0;
    for(j = 0; j < count; j ++)
    {
        min = NTTW_BIGINT_MAX;
        for(k = 0; k < count; k ++)
        {
            if(L1[k] < min && L1[k] > 0)
            {
                min = L1[k];
                index = k;
            }
        }

        setXY(anglesSorted[j],getX(anglesRaw[index]),getY(anglesRaw[index]));
        L1[index] = -1;
    }
    free_vectorArray(anglesRaw,size);
    free_array(L1);

    index = 0;
    count = 0;
    ///Determine Finite angles for sorted L1 set
    while(index < angleSize && count < size) // 1/1
    {
        setXY(nextAngle,getX(anglesSorted[index]),getY(anglesSorted[index]));
        //fprintf(stderr,"%li/%li,\t",getY(nextAngle),getX(nextAngle));

        ///First Octant
        computeAngle(N, getX(nextAngle), getY(nextAngle), mLookUp, sLookUp, farey, finite, perps, &count);
        ///Second Octant
        computeAngle(N, getY(nextAngle), getX(nextAngle), mLookUp, sLookUp, farey, finite, perps, &count);
        ///Third Octant
        computeAngle(N, getX(nextAngle), -getY(nextAngle), mLookUp, sLookUp, farey, finite, perps, &count);
        ///Forth Octant
        computeAngle(N, getY(nextAngle), -getX(nextAngle), mLookUp, sLookUp, farey, finite, perps, &count);

        index ++;
    }

    ///Deallocate
    free_array(mLookUp);
    free_array(sLookUp);
    free_vectorArray(anglesSorted,angleSize);
    free_vector(angle1);
    free_vector(angle2);
    free_vector(nextAngle);

    return count;
}

void mt(const vector *angles, const nttw_integer *data, mojetteProjection *set, const size_t mu, const size_t rows, const size_t cols)
{
    size_t n, translateMojette = 0, offsetMojette = 0, x = 0, y = 0;

	///Compute Fast Mojette
	for(n = 0; n < mu; n ++)
	{
        if(getY(angles[n])*getX(angles[n]) >= 0) //If not negative slope
            offsetMojette = getX(angles[n])*(rows-1);
        else
            offsetMojette = 0;

        for(x = 0; x < rows; x ++)
            for(y = 0; y < cols; y ++)
            {
                if(getY(angles[n])*getX(angles[n]) >= 0)
                    translateMojette = getY(angles[n])*y - getX(angles[n])*x + offsetMojette; //GetY = q, GetX = p
                else
                    translateMojette = getX(angles[n])*x - getY(angles[n])*y; //GetY = q, GetX = p
                set[n].bins[translateMojette] += data[x*rows+y];
            }
	}
}

nttw_integer* fmt(const size_t P, const size_t Q, const nttw_integer *image, const vector *farey, const nttw_big_integer *finite, nttw_integer *perps, const size_t N)
{
    const size_t size = N+N/2, angleCount = N+N/2;
    size_t n, row, translateMojette = 0, translateFinite = 0, x = 0, y = 0;
    nttw_integer *frtSpace;
    nttw_big_integer inv;

    frtSpace = array_1D(size*N);
    init_1D(frtSpace,size*N,0);

	///Compute Fast Mojette
	for(n = 0; n < angleCount; n ++)
	{
	    if(perps[n])
	    {
            inv = minverse_big(abs(getX(farey[n])),N);
            row = N + finite[n];
	    }
	    else
	    {
            inv = minverse_big(abs(getY(farey[n])),N);
            row = finite[n];
	    }

        for(x = 0; x < Q; x ++)
            for(y = 0; y < P; y ++)
            {
                if(getY(farey[n])*getX(farey[n]) >= 0 && !perps[n])
                    translateMojette = getY(farey[n])*y - getX(farey[n])*x; //GetY = q, GetX = p
                else
                    translateMojette = getX(farey[n])*x - getY(farey[n])*y; //GetY = q, GetX = p
                translateFinite = (inv*translateMojette)%N;
                frtSpace[row*N+translateFinite] += image[x*Q+y];
                //frtSpace[row*N+translateFinite] += (nttw_integer) Normal(image[x*Q+y],(1-SNR)*image[x*Q+y]);
            }
	}

	return frtSpace;
}

nttw_integer* fmt_noise(const size_t P, const size_t Q, const nttw_integer *image, const vector *farey, const nttw_big_integer *finite, nttw_integer *perps, const size_t N, const float SNR)
{
    const size_t size = N+N/2, angleCount = N+N/2;
    size_t n, row, translateMojette = 0, translateFinite = 0, x = 0, y = 0;
    nttw_integer *frtSpace;
    nttw_big_integer inv;

    frtSpace = array_1D(size*N);
    init_1D(frtSpace,size*N,0);

	///Compute Fast Mojette
	for(n = 0; n < angleCount; n ++)
	{
	    if(perps[n])
	    {
            inv = minverse_big(abs(getX(farey[n])),N);
            row = N + finite[n];
	    }
	    else
	    {
            inv = minverse_big(abs(getY(farey[n])),N);
            row = finite[n];
	    }

        for(x = 0; x < Q; x ++)
            for(y = 0; y < P; y ++)
            {
                if(getY(farey[n])*getX(farey[n]) >= 0 && !perps[n])
                    translateMojette = getY(farey[n])*y - getX(farey[n])*x; //GetY = q, GetX = p
                else
                    translateMojette = getX(farey[n])*x - getY(farey[n])*y; //GetY = q, GetX = p
                translateFinite = (inv*translateMojette)%N;
                frtSpace[row*N+translateFinite] += (nttw_integer) Normal(image[x*Q+y],(1.0-SNR)*image[x*Q+y]);
            }
	}

	return frtSpace;
}

void mt2frt(const vector *farey, const nttw_big_integer *finite, const nttw_integer *perps, const mojetteProjection *set, const size_t mu, const size_t rows, nttw_integer *frtSpace, const size_t N)
{
    const size_t size = N+N/2;
    size_t n, j, proj, translateOffset = 0, translateMojette = 0, translateFinite = 0;
    nttw_big_integer inv;
    long angleSign;

    init_1D(frtSpace,size*N,0);

	///Compute Fast Mojette
	for(n = 0; n < mu; n ++)
	{
	    angleSign = getY(farey[n])*getX(farey[n]);

	    if(perps[n])
	    {
            inv = minverse_big(abs(getX(farey[n])),N);
            proj = N + finite[n];
	    }
	    else
	    {
            inv = minverse_big(abs(getY(farey[n])),N);
            proj = finite[n];
	    }

	    if(angleSign >= 0 && !perps[n]) //If not negative slope
            translateOffset = getX(farey[n])*(rows-1);
        else if(angleSign >= 0) //If not negative slope and perp
            translateOffset = getX(farey[n])*(rows-1);
        else
            translateOffset = 0;

        for(j = 0; j < set[n].binMax; j ++)
        {
            if(angleSign >= 0 && perps[n]) //Reverse for perp
                translateMojette = translateOffset - j;
            else
                translateMojette = j - translateOffset;
            translateFinite = (inv*translateMojette)%N;
            frtSpace[proj*N+translateFinite] += set[n].bins[j];
        }
	}
}

void mt2frt_noise(const vector *farey, const nttw_big_integer *finite, const nttw_integer *perps, const mojetteProjection *set, const size_t mu, const size_t rows, nttw_integer *frtSpace, const size_t N, const float SNR)
{
    const size_t size = N+N/2;
    size_t n, j, proj, translateOffset = 0, translateMojette = 0, translateFinite = 0;
    nttw_big_integer inv;
    long angleSign;

    init_1D(frtSpace,size*N,0);

	///Compute Fast Mojette
	for(n = 0; n < mu; n ++)
	{
	    angleSign = getY(farey[n])*getX(farey[n]);

	    if(perps[n])
	    {
            inv = minverse_big(abs(getX(farey[n])),N);
            proj = N + finite[n];
	    }
	    else
	    {
            inv = minverse_big(abs(getY(farey[n])),N);
            proj = finite[n];
	    }

	    if(angleSign >= 0 && !perps[n]) //If not negative slope and perp proj
            translateOffset = getX(farey[n])*(rows-1);
        else if(angleSign >= 0) //If not negative slope and perp
            translateOffset = getX(farey[n])*(rows-1);
        else
            translateOffset = 0;

        for(j = 0; j < set[n].binMax; j ++)
        {
            if(angleSign >= 0 && perps[n]) //Reverse for perp
                translateMojette = translateOffset - j;
            else
                translateMojette = j - translateOffset;
            translateFinite = (inv*translateMojette)%N;
            frtSpace[proj*N+translateFinite] += (nttw_integer) Normal(set[n].bins[j],(1.0-SNR)*set[n].bins[j]);
        }
	}
}
