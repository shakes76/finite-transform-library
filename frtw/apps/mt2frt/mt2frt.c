/**
* Mojette Transform to FRT
* Converts the MT of an image to the FRT of the same image
* Outputs data result.
* \author Shekhar S. Chandra, 2008-10
*/
//NTTW
#include <nttw/array.h>
#include <nttw/timing.h>
#include <nttw/image.h>

//CGhosts 2 Lib
#include "array_complex.h"
#include "radon.h"
#include "mojette.h"

///Prototypes
int readAngleFile(vector **angles, nttw_big_integer **finite, nttw_integer **perps, size_t *projNo, const char* filename);

int main(int argc, char *argv[])
{
    nttw_integer *frtSpace, *perpFlags;
    mojetteProjection *set;
    vector *angles;
    int binaryFile = FALSE, dyadicSize = TRUE;
    size_t n, N, size, mu;
    unsigned long long duration = 0;
    nttw_big_integer *finiteAngles;
    long *result, *resultCropped;
    float SNR = 1.0;
    FILE *outFile;

    printf(">| Mojette Transform to FRT Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");
    printf(">| Machine Integer Size of %u bits\n", BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 7)
    {
        printf(">| Usage: %s <filename> <angles> <n> <N> <SNR> <output>\n",argv[0]);
        printf(">| filename is loaded and converted FRT projections in output.\n");
        printf(">| n is the original image size and N is the overall image size.\n");
        printf(">| SNR is Signal-to-noise ratio. Set to less than 1.0 to add noise.\n");
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }
    n = atoi(argv[3]);
    N = atoi(argv[4]);
    SNR = atof(argv[5]);

    if(N % 2 == 1) ///Assume prime if odd
        dyadicSize = FALSE;

    //--------------------------------------------
    ///Load Image
    fprintf(stderr,">| Reading... ");
    if( !readMojetteSet(&set,&mu,argv[1]) )
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }

    if( !readAngleFile(&angles, &finiteAngles, &perpFlags, &mu, argv[2]) )
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[2]);
        return EXIT_FAILURE;
    }
    fprintf(stderr,"Done.\n");

    ///Output parameters to be used
    if(dyadicSize)
        size = N+N/2;
    else
        size = N+1;
    frtSpace = array_1D(N*size);
    result = arraySigned_1D(N*N);

    //--------------------------------------------
    ///Fast FRT
    fprintf(stderr,">| Transforming... ");
    START_TIMER;

    if(SNR < 1.0)
        mt2frt_noise(angles,finiteAngles,perpFlags,set,mu,n,frtSpace,N,SNR);
    else
        mt2frt(angles,finiteAngles,perpFlags,set,mu,n,frtSpace,N);

    ///Invert
    if(dyadicSize)
        ifrt_dyadic_signed(frtSpace, result, N, TRUE);
    else
        ifrt_signed(frtSpace, result, N, TRUE);

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);

    ///Write
    fprintf(stderr,">| Writing File... ");
    ///Save Result
    writePGM(frtSpace,size,N,255,argv[6],binaryFile);
    //fprintf(stderr,">| Writing File... ");
    //writeSignedPGM(result,N,N,255,argv[6],binaryFile);

    resultCropped = truncateImage(result,N,N,n,n); ///Crop result

    writeSignedPGM(resultCropped,n,n,255,"mt2frt_cropped.pgm",binaryFile);
    fprintf(stderr,">| Cropped output name: mt2frt_cropped.pgm\n");

    ///Output/Append times to global data file
    printf(">| Global result output to %s.\n","results.csv");
    openFile_Append("results.csv",&outFile,binaryFile);
    fprintf(outFile,"%llu ", duration); //Write initial data and time
    ///Format will be k, angle type, angle time, mt time, mt2frt time, mse, rmse
    fclose(outFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(frtSpace);
    free_array(perpFlags);
    free_array(finiteAngles);
    free_mojetteArray(set,mu);
    free_vectorArray(angles, mu);
    free_array(result);
    free_array(resultCropped);

    return EXIT_SUCCESS;
}

int readAngleFile(vector **angles, nttw_big_integer **finite, nttw_integer **perps, size_t *projNo, const char* filename)
{
    size_t j, success = FALSE;
    long p, q;
    nttw_big_integer m, flag;
    FILE *inFile;

    inFile = fopen(filename,"r");
    if(inFile == NULL)
        return FALSE;

    ///Read Angle File. Format
    ///m, p, q
    success = fscanf(inFile,"%lu",projNo);

    *angles = vectorArray_2D(*projNo);
    *finite = array_1D_big(*projNo);
    *perps = array_1D(*projNo);
    init_1D(*perps,*projNo,FALSE);

    for(j = 0; j < *projNo; j ++)
    {
        success = fscanf(inFile,"%lu %li %li %lu",&m, &p, &q, &flag);
        setXY( (*angles)[j], p, q );
        (*finite)[j] = m;
        if(flag == 1)
            (*perps)[j] = TRUE;
    }

    for(j = 0; j < *projNo; j ++)
        printf( "%li/%li - m:%lu s?:%lu, ", getY((*angles)[j]), getX((*angles)[j]), (*finite)[j], (*perps)[j] );
    printf("\n");

    fclose(inFile);
    return TRUE;
}
