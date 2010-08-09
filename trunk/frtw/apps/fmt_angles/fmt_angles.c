/**
* Fast Mojette Angle Set Generate App
* Generates the FMT angle set.
* Outputs text result.
* \author Shekhar S. Chandra, 2008-10
*/
//NTTW
#include <nttw/array.h>
#include <nttw/timing.h>
#include <nttw/image.h>

//CGhosts 2 Lib
#include "mojette.h"

int main(int argc, char *argv[])
{
    size_t size, Q, n, N;
    int binaryFile = FALSE, forceL1 = FALSE;
    unsigned long long duration = 0;
    nttw_integer *perpFlags;
    nttw_big_integer *finiteAngles, totalAngles = 0;
    vector *angles;
    FILE *outFile;

    printf(">| Fast Mojette Angle Set Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 5)
    {
        printf(">| Usage: %s <n> <N> <Angle Type> <output>\n",argv[0]);
        printf(">| Angle set for an nxn image for NxN FMT space is generated.\n");
        printf(">| Angle Type is an integer 0-2 where 1 - l1, 2 - Simple.\n");
        printf(">| N must be dyadic and output to stdout is made when N < 64.\n");
        return EXIT_FAILURE;
    }
    Q = n = atoi(argv[1]);
    N = atoi(argv[2]);
    forceL1 = atoi(argv[3]);

    //--------------------------------------------
    ///Arrays to be used
    size = N+N/2;

    perpFlags = array_1D(size);
    finiteAngles = array_1D_big(size);
    angles = vectorArray_2D(size);

    init_1D(perpFlags,size,FALSE);

    fprintf(stderr,">| Generating... ");
    START_TIMER;
if(!forceL1)
{
    //--------------------------------------------
    ///Fast FRT Angles
    totalAngles = fmt_angleSet(N, n, n, angles, finiteAngles, perpFlags);
}
else if(forceL1 == 2)
{
    //--------------------------------------------
    ///Use simple Farey set for Finite Set
    ///Fast FRT Angles
    totalAngles = fmt_angleSet_Simple(N, angles, finiteAngles, perpFlags);
}
else
{
    //--------------------------------------------
    ///Force L1 while creating Finite Set
    ///Fast FRT Angles
    totalAngles = fmt_angleSet_L1(N, n, n, angles, finiteAngles, perpFlags);
}
    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);
    printf(">| Acquired %lu of %lu angles.\n", size, totalAngles);

    //--------------------------------------------
    ///Output Result to Std if small (N < 64)
    if(N < 64)
    {
        printf(">| Selected: \n");
        for(n = 0; n < size; n ++)
            printf("%li/%li, ",getY(angles[n]),getX(angles[n]));
        printf("\n");

        printf(">| Finite Angles: \n");
        for(n = 0; n < size; n ++)
            printf(FORMAT_BIG_CSVOUTPUT_STR,finiteAngles[n]);
        printf("\n");

        printf(">| Perp Flags: \n");
        for(n = 0; n < size; n ++)
            printf(FORMAT_CSVOUTPUT_STR,perpFlags[n]);
        printf("\n");
    }

    ///Output to file
    printf(">| Output to file: %s\n",argv[4]);
    printf(">| Format is \t m   p   q   s?\n");
    openFile_Write(argv[4],&outFile,binaryFile);
    fprintf(outFile,"%lu\n",size);
    for(n = 0; n < size; n ++)
        fprintf(outFile,"%li %li %li %i\n",finiteAngles[n],getX(angles[n]),getY(angles[n]),perpFlags[n]);
    fclose(outFile);

    ///Output Diophantine File
    printf(">| Diophantine Output is %s, which opens in DGV.\n","angles.dio");
    openFile_Write("angles.dio",&outFile,binaryFile);
    fprintf(outFile,"%lu %lu %i %i %i\n",N,N,0,0,0); //pri lattice
    fprintf(outFile,"%lu %lu %i %i %i\n",Q,Q,0,0,0); //sec lattice
    fprintf(outFile,"%i %i %i\n",255,255,255); //bgd
    fprintf(outFile,"%i %i %i\n",0,0,192); //change colour
    for(n = 0; n < size; n ++)
        fprintf(outFile,"%lu %lu %li %li\n",Q/2-1,Q/2-1,getY(angles[n]),getX(angles[n]));
    fclose(outFile);

    ///Output/Append times to global data file
    printf(">| Global result output to %s.\n","results.csv");
    openFile_Append("results.csv",&outFile,binaryFile);
    Q = (size_t) (N/(double)Q);
    fprintf(outFile,"%lu %lu %llu ", Q, forceL1, duration); //Write initial data and time
    ///Format will be k, angle type, angle time, mt time, mt2frt time, mse, rmse
    fclose(outFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(perpFlags);
    free_array(finiteAngles);
    free_vectorArray(angles,size);

    return EXIT_SUCCESS;
}
