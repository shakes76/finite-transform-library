/**
* Katz Program
* Computes the Farey angles up to Katz criterion.
* \author Shekhar S. Chandra, 2009
*/
#include <stdio.h>
#include <math.h>

//NTTW Library
#include <nttw/timing.h>
//FRTW Library
#include "mojette.h"

///Globals
const int debugMode = FALSE;

///Function Prototypes
int processArgs(int argcount, char *argvars[], size_t *p, size_t *n);

int main(int argc, char *argv[])
{
    int success = FALSE;
    size_t j, k, size, N, P, Q, angleCount = 0, index = 0, tick = 0;
    vector angle1, angle2, nextAngle, *angles, *anglesSorted;
    unsigned long long duration = 0;
    nttw_big_integer *bins, *binsSorted, min, sumX, sumY;
    FILE *outFile;

    printf(">| Katz Criterion Program.\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %lu bits\n",(long)BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    //--------------------------------------------
    ///Process arguments
	success = processArgs(argc,argv,&P,&Q);
	if(!success)
	{
	    fprintf(stderr,">| Check arguments and try again.\n");
        return EXIT_FAILURE;
	}
	printf(">| Calculating Katz for %lux%lu image and output will be to %s.\n",P,Q,argv[3]);

	//--------------------------------------------
	///Calculate Katz
	if(P > Q)
        N = P;
    else
        N = Q;
    size = totalFareyAngles(N);

    angle1 = vector_2D();
    angle2 = vector_2D();
    nextAngle = vector_2D();
    angles = vectorArray_2D(size);
    bins = array_1D_big(size);

    setY(angle1,0);
    setX(angle1,1); // 0/1
    setXY(angles[0],getX(angle1),getY(angle1));
    bins[0] = abs(getX(angle1)) + abs(getY(angle1));
    setY(angle2,1);
    setX(angle2,N); // 1/n
    setXY(angles[1],getX(angle2),getY(angle2));
    bins[1] = abs(getX(angle2)) + abs(getY(angle2));
    angleCount = 2;

    fprintf(stderr,"N(%i): %i, Farey(%i): \n",N,size,N);
    fprintf(stderr,"%li/%li,\t",getY(angles[0]),getX(angles[0]));
    fprintf(stderr,"%li/%li,\t",getY(angles[1]),getX(angles[1]));

    START_TIMER;
    //while( !(isKatzCriterion(angles,angleCount,P,Q)) && !(getY(nextAngle) == 1 && getX(nextAngle) == 1) ) // 1/1
    while( !(getY(nextAngle) == 1 && getX(nextAngle) == 1) ) // 1/1
    {
        nextFareyAngle_Compact(N, angle1, angle2, &nextAngle);
        fprintf(stderr,"%li/%li,\t",getY(nextAngle),getX(nextAngle));

        setY(angle1, getY(angle2));
        setX(angle1, getX(angle2));
        setY(angle2, getY(nextAngle));
        setX(angle2, getX(nextAngle));

        setXY(angles[angleCount],getX(nextAngle),getY(nextAngle));
        bins[angleCount] = abs(getX(nextAngle)) + abs(getY(nextAngle));
        angleCount ++;
    }
    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"\n>| Done. Took %llu usecs.\n",duration);
    printf(">| Total angles: %lu\n",angleCount);

    printf(">| l1: \n");
    for(j = 0; j < angleCount; j ++)
        printf(FORMAT_BIG_CSVOUTPUT_STR,bins[j]);
    printf("\n");

    printf(">| Sorting Angles...\n");
    anglesSorted = vectorArray_2D(angleCount);
    binsSorted = array_1D_big(angleCount);
    for(j = 0; j < angleCount; j ++)
    {
        min = NTTW_BIGINT_MAX;
        for(k = 0; k < angleCount; k ++)
        {
            if(bins[k] < min && bins[k] > 0)
            {
                min = bins[k];
                index = k;
            }
        }

        binsSorted[j] = bins[index];
        setXY(anglesSorted[j],getX(angles[index]),getY(angles[index]));
        bins[index] = -1;
    }
    free_vectorArray(angles,size);
    free_array(bins);

    printf(">| l1 Sorted: \n");
    for(j = 0; j < angleCount; j ++)
        printf(FORMAT_BIG_CSVOUTPUT_STR,binsSorted[j]);
    printf("\n");
    printf(">| Angles Sorted: \n");
    for(j = 0; j < angleCount; j ++)
        printf("%li/%li, ",getX(anglesSorted[j]),getY(anglesSorted[j]));
    printf("\n");

    //--------------------------------------------
    ///Do Katz
    index = 0;
    sumX = 0;
    sumY = 0;
    printf(">| Katz Angle Set:\n");
    while( (sumX < P && sumY < Q) && index < angleCount)
    {
        if(tick == 0)
        {
            sumX += getX(anglesSorted[index]);
            sumY += getY(anglesSorted[index]);
            min = (P-1)*abs(getX(anglesSorted[index])) + (Q-1)*abs(getY(anglesSorted[index])) + 1;
            printf("%li/%li, ",getX(anglesSorted[index]),getY(anglesSorted[index]));
            tick ++;
        }
        else if(tick == 1) //Do +45 deg
        {
            if(getY(anglesSorted[index]) != getX(anglesSorted[index]))
            {
                sumX += getY(anglesSorted[index]);
                sumY += getX(anglesSorted[index]);
                min = (P-1)*abs(getX(anglesSorted[index])) + (Q-1)*abs(getY(anglesSorted[index])) + 1;
                printf("%li/%li, ",getY(anglesSorted[index]),getX(anglesSorted[index]));
            }
            tick ++;
        }
        else if(tick == 2)
        {
            if(getY(anglesSorted[index]) != 0)
            {
                sumX += getX(anglesSorted[index]);
                sumY += getY(anglesSorted[index]);
                min = (P-1)*abs(getX(anglesSorted[index])) + (Q-1)*abs(getY(anglesSorted[index])) + 1;
                printf("%li/%li, ",getX(anglesSorted[index]),-getY(anglesSorted[index]));
            }
            tick ++;
        }
        else if(tick == 3) //Do +45 deg
        {
            if(getY(anglesSorted[index]) != getX(anglesSorted[index]) && getX(anglesSorted[index]) != 0)
            {
                sumX += getY(anglesSorted[index]);
                sumY += getX(anglesSorted[index]);
                min = (P-1)*abs(getX(anglesSorted[index])) + (Q-1)*abs(getY(anglesSorted[index])) + 1;
                printf("%li/%li, ",getY(anglesSorted[index]),-getX(anglesSorted[index]));
            }
            tick ++;
            index ++;
        }

        tick %= 4;
    }
    printf("\n>| Katz Angle Set Done (Sum X: %lu, Sum Y: %lu)\n",sumX,sumY);
    printf(">| Max #Bins: %lu\n",min);

    outFile = fopen(argv[3],"a");
    fprintf(outFile,FORMAT_BIG_OUTPUT_STR,N);
    fprintf(outFile,FORMAT_BIG_OUTPUT_STR,min);
    fprintf(outFile,"\n");
    fclose(outFile);

    //--------------------------------------------
    free_array(binsSorted);
    free_vectorArray(anglesSorted,angleCount);
    printf(">| Operation Complete\n");
	return EXIT_SUCCESS;
}

int processArgs(int argcount, char *argvars[], size_t *P, size_t *Q)
{
    if(argcount == 4)
    {
        *P = atoi(argvars[1]);
        *Q = atoi(argvars[2]);
        return TRUE;
    }
    else
    {
        fprintf(stderr,">| Usage: %s <P> <Q> <output>\n",argvars[0]);
        fprintf(stderr,">| Computes the Farey angles up to Katz criterion.\n");
        fprintf(stderr,">| Output will be a tuple (text) file.\n");
        return FALSE;
    }
}
