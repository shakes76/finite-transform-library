/**
* Inverse Fast Mojette
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
    nttw_integer *frtSpace, normResult = TRUE, deNoise = FALSE;
    nttw_big_integer Isum = 0, iterations = 1;
    long *result, *resultCropped, *samples, *frtSpaceSamples, *frtSpaceNoise;
    size_t j, k, n, N, size, binaryFile = FALSE;
    unsigned long long durationAngle = 0, durationTransform = 0;
    double scale;

    printf(">| Inverse Fast Mojette Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 6)
    {
        printf(">| Usage: %s <filename> <n> <DeNoise?> <Iterations> <output>\n",argv[0]);
        printf(">| filename is loaded and transformed to produce output.\n");
        printf(">| Iterations is an integer and DeNoise is a Boolean.\n");
        printf(">| Iterations is ignored if DeNoise is false.\n");
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }
    n = atoi(argv[2]);
    deNoise = atoi(argv[3]);
    iterations = atoi(argv[4]);

    //--------------------------------------------
    ///Load Image
    if(!readPGM(&frtSpace,&size,&N,argv[1],binaryFile))
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }

    result = arraySigned_1D(N*N);
    samples = arraySigned_1D(N*N);
    frtSpaceNoise = arraySigned_1D(size*N);
    frtSpaceSamples = arraySigned_1D(size*N);

    initSigned_1D(samples,N*N,1);

    //--------------------------------------------
//    ///Fast FRT
//    fprintf(stderr,">| Regularise Projections... ");
//    START_TIMER;
//
//    trueIsum = 0;
//    for(k = 0; k < N; k ++)
//        trueIsum += frtSpace[k];
//
//    ///Standardise Isum
//    for(j = 0; j < size; j ++)
//    {
//        Isum = 0;
//        for(k = 0; k < N; k ++)
//            Isum += frtSpace[j*N+k];
//        Isum -= trueIsum;
//        Isum /= N;
//        for(k = 0; k < N; k ++)
//            frtSpace[j*N+k] -= Isum;
//    }
//
//    STOP_TIMER;
//    durationAngle = MICROSECONDS_ELAPSED;
//    fprintf(stderr,"Done. Time: %llu usecs.\n", durationAngle);
//    printf(">| Acquired %lu of %lu angles.\n", size, totalAngles);

//    //--------------------------------------------
//    ///Output Result to Std if small (N < 64)
//    if(N < 64)
//    {
//        printf(">| Selected: \n");
//        for(n = 0; n < size; n ++)
//            printf("%li/%li, ",getY(angles[n]),getX(angles[n]));
//        printf("\n");
//
//        printf(">| Finite Angles: \n");
//        for(n = 0; n < size; n ++)
//            printf(FORMAT_BIG_CSVOUTPUT_STR,finiteAngles[n]);
//        printf("\n");
//
//        printf(">| Perp Flags: \n");
//        for(n = 0; n < size; n ++)
//            printf(FORMAT_CSVOUTPUT_STR,perpFlags[n]);
//        printf("\n");
//    }

    //--------------------------------------------
	fprintf(stderr,">| Computing Inversion... \n");
	START_TIMER;

    ifrt_dyadic_signed(frtSpace,result,N,normResult);

    STOP_TIMER;
    durationTransform = MICROSECONDS_ELAPSED;
    fprintf(stderr,">| Done. Took %llu usecs.\n",durationTransform);
    writeSignedPGM(result,N,N,255,"ifmt_original.pgm",binaryFile);

    //--------------------------------------------
    ///Denoise if requested
if(deNoise)
{
    while(iterations --)
    {
        for(j = 0; j < n; j ++)
        {
            for(k = 0; k < n; k ++)
                result[j*N+k] = samples[j*N+k] = 0;
        }

        frt_dyadic_signed(result,frtSpaceNoise,N);
        frt_dyadic_signed(samples,frtSpaceSamples,N);

        writeSignedPGM(result,N,N,255,"ifmt_redundant.pgm",binaryFile);
        writeSignedPGM(frtSpaceNoise,size,N,255,"ifmt_noise_est.pgm",binaryFile);
        writeSignedPGM(frtSpaceSamples,size,N,255,"ifmt_sample_est.pgm",binaryFile);
        writeSignedPGM(samples,N,N,255,"ifmt_samples.pgm",binaryFile);

        for(j = 0; j < size; j ++)
        {
            for(k = 0; k < N; k ++)
            {
//                if(frtSpaceSamples[j*N+k] > 0)
//                    scale = 1.0*frtSpaceSamples[j*N+k]/N;
//                else
//                    scale = 1.0;
//                scale = frtSpaceNoise[j*N+k]/scale;
//                scale = frtSpaceNoise[j*N+k];
//                frtSpaceSamples[j*N+k] = frtSpace[j*N+k];
//                scale = (double) frtSpaceSamples[j*N+k] - scale;
//                frtSpaceNoise[j*N+k] = scale + 0.5; //Round

                //frtSpace[j*N+k] = frtSpaceSamples[j*N+k];

                ///No scaling
                frtSpaceSamples[j*N+k] = frtSpace[j*N+k];
                frtSpaceSamples[j*N+k] -= frtSpaceNoise[j*N+k];
            }
        }

        writeSignedPGM(frtSpaceSamples,size,N,255,"ifmt_corrected.pgm",binaryFile);

        fprintf(stderr,">| Computing Inversion... \n");
        START_TIMER;

        ifrt_dyadic_signed2(frtSpaceSamples,result,N,normResult);

        STOP_TIMER;
        durationTransform = MICROSECONDS_ELAPSED;
        fprintf(stderr,">| Done. Took %llu usecs.\n",durationTransform);
    }
}

    ///Write result
    writeSignedPGM(result,N,N,255,argv[5],binaryFile);

    resultCropped = truncateImage(result,N,N,n,n); ///Crop result

    writeSignedPGM(resultCropped,n,n,255,"ifmt_cropped.pgm",binaryFile);
    fprintf(stderr,">| Cropped output name: ifmt_cropped.pgm\n");

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(frtSpace);
    free_array(frtSpaceNoise);
    free_array(frtSpaceSamples);
    free_array(result);
    free_array(resultCropped);
    free_array(samples);

    return EXIT_SUCCESS;
}
