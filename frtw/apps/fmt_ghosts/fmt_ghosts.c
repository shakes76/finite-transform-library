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

//FRTW Lib
#include "mojette.h"
#include "ghosts.h"

///Globals

///Prototypes
int readAngleFile(vector **angles, nttw_big_integer **finite, nttw_integer **perps, size_t *projNo, const char* filename);

int main(int argc, char *argv[])
{
    nttw_integer *ghostImage, *perpFlags, *eigenvalues;
    nttw_integer root, modulus, proot;
    vector *angles;
    int binaryFile = FALSE, dyadicSize = TRUE;
    size_t n, N, size, mu, noOfGhosts = 0, count = 0;
    unsigned long long duration = 0;
    nttw_big_integer *finiteAngles, *missingAngles, *missingAngleLookup;

    printf(">| Mojette Transform to FRT Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");
    printf(">| Machine Integer Size of %u bits\n", BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 4)
    {
        printf(">| Usage: %s <angles> <N> <output>\n",argv[0]);
        printf(">| filename is loaded and converted FRT projections in output.\n");
        printf(">| n is the original image size and N is the overall image size.\n");
        printf(">| SNR is Signal-to-noise ratio. Set to less than 1.0 to add noise.\n");
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }
    N = atoi(argv[2]);

    if(N % 2 == 1) ///Assume prime if odd
        dyadicSize = FALSE;

    //--------------------------------------------
    ///Load Image
    fprintf(stderr,">| Reading... ");
    if( !readAngleFile(&angles, &finiteAngles, &perpFlags, &mu, argv[1]) )
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }
    fprintf(stderr,"Done.\n");

    ///Output parameters to be used
    if(dyadicSize)
        size = N+N/2;
    else
        size = N+1;

    noOfGhosts = size-mu;

    ghostImage = array_1D(N*N);

    //--------------------------------------------
    ///Determine which finite angles are missing
    determineMissingAngles(N, finiteAngles, perpFlags, mu, &missingAngles, &missingAngleLookup);

    fprintf(stderr, "Missing Lookup: \n");
    for(n = 0; n < size; n ++)
        fprintf(stderr, "%lu, ", missingAngleLookup[n]);

    fprintf(stderr, "\nMissing Projections at: \n");
    for(n = 0; n < noOfGhosts; n ++)
        fprintf(stderr, "%lu, ", missingAngles[n]);

    //--------------------------------------------
    ///Ghost Eigenvalues
    fprintf(stderr,"\n>| Eigenvalues of Ghosts... ");
    START_TIMER;

    eigenvalues = ghostEigenvalues(N, &root, &modulus, &proot);

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);

    //--------------------------------------------
    ///Fast Ghost Construction
    fprintf(stderr,">| Forming Ghosts... ");
    START_TIMER;

    ///Ghost convolution
    ghostImage = convolveGhosts_1D(N, eigenvalues, root, modulus, proot, missingAngles, missingAngleLookup, noOfGhosts);

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);

    ///Write
    fprintf(stderr,">| Writing Eigenvalues File... ");
    ///Save Result
    writePGM(eigenvalues,N,N,255,"eigenvalues.pgm",binaryFile);

    fprintf(stderr,">| Writing Ghost Image File... ");
    writePGM(ghostImage,N,N,255,argv[3],binaryFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(perpFlags);
    free_array(finiteAngles);
    free_array(missingAngles);
    free_array(missingAngleLookup);
    free_vectorArray(angles, mu);
    free_array(eigenvalues);
    free_array(ghostImage);

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
