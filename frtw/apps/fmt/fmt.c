/**
* Fast Mojette
* Computes the FMT of an image provided
* Outputs image result.
* \author Shekhar S. Chandra, 2008-10
*/
//NTTW
#include <nttw/array.h>
#include <nttw/timing.h>
#include <nttw/image.h>

//CGhosts 2 Lib
#include "mojette.h"
#include "radon.h"

const float SNR = 0.97;

int main(int argc, char *argv[])
{
    nttw_integer *image, *frtSpace, *perpFlags, type, addNoise = FALSE;
    nttw_big_integer *finiteAngles, totalAngles = 0;
    vector *angles;
    long *result;
    size_t n, N, size, Q, P, binaryFile = FALSE;
    unsigned long long durationAngle = 0, durationTransform = 0;

    printf(">| Fast Mojette Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 6)
    {
        printf(">| Usage: %s <filename> <N> <Angle Type> <Noise?> <output>\n",argv[0]);
        printf(">| filename is loaded and transformed to produce output.\n");
        printf(">| N is the size of the FFT space to be used.\n");
        printf(">| Angle Type is an integer 0-2 where 1 - l1, 2 - Simple.\n");
        printf(">| Noise is a Boolean and output is the FRT space.\n");
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }
    N = atoi(argv[2]);
    type = atoi(argv[3]);
    addNoise = atoi(argv[4]);
    size = N+N/2;

    //--------------------------------------------
    ///Load Image
    if(!readPGM(&image,&Q,&P,argv[1],binaryFile))
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }
    if(P > Q)
        n = P;
    else
        n = Q;

    perpFlags = array_1D(size);
    finiteAngles = array_1D_big(size);
    angles = vectorArray_2D(size);
    result = arraySigned_1D(N*N);

    init_1D(perpFlags,size,FALSE);

    //--------------------------------------------
    ///Fast FRT
    fprintf(stderr,">| Creating Angle Set... ");
    START_TIMER;
if(type == 0)
{
    ///Fast FRT Angles
    fprintf(stderr,"via a possible L1 set... ");
    totalAngles = fmt_angleSet(N, n, n, angles, finiteAngles, perpFlags);
}
else if(type == 2)
{
    ///Use simple Farey set for Finite Set
    ///Fast FRT Angles
    fprintf(stderr,"via the simple set... ");
    totalAngles = fmt_angleSet_Simple(N, angles, finiteAngles, perpFlags);
}
else
{
    ///Force L1 while creating Finite Set
    ///Fast FRT Angles
    fprintf(stderr,"via the L1 set... ");
    totalAngles = fmt_angleSet_L1(N, n, n, angles, finiteAngles, perpFlags);
}
    STOP_TIMER;
    durationAngle = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", durationAngle);
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

    //--------------------------------------------
    fprintf(stderr,">| Transforming... ");
    START_TIMER;

    if(addNoise)
        frtSpace = fmt_noise(P, Q, image, angles, finiteAngles, perpFlags, N, SNR);
    else
        frtSpace = fmt(P, Q, image, angles, finiteAngles, perpFlags, N);

    STOP_TIMER;
    durationTransform = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", durationTransform);
    fprintf(stderr,">| Total Time: %llu usecs.\n", durationAngle+durationTransform);

    ///TEST RESULT - iFRT
	fprintf(stderr,">| TEST: Computing Inversion... \n");
	START_TIMER;

    ifrt_dyadic_signed(frtSpace,result,N,TRUE);

    STOP_TIMER;
    durationTransform = MICROSECONDS_ELAPSED;
    fprintf(stderr,">| Done. Took %llu usecs.\n",durationTransform);

    ///Write result
    writePGM(frtSpace,size,N,255,argv[5],binaryFile);
    writeSignedPGM(result,N,N,255,"fmt_result.pgm",binaryFile);
    fprintf(stderr,">| Test Inversion output to fmt_result.pgm.\n");

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(image);
    free_array(perpFlags);
    free_array(finiteAngles);
    free_vectorArray(angles,size);
    free_array(frtSpace);
    free_array(result);

    return EXIT_SUCCESS;
}
