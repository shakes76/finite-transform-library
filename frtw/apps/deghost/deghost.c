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
    nttw_integer *deghostImage, *perpFlags, *frtSpace, *ghost;
    vector *angles;
    int binaryFile = FALSE, dyadicSize = TRUE, readN, readSize;
    size_t n, N, size, mu, noOfGhosts = 0;
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

    if(argc != 5)
    {
        printf(">| Usage: %s <angles> <FRT Space> <Ghost> <output>\n",argv[0]);
        printf(">| FRT Space and Ghost are PGM images.\n");
        printf(">| Output will contain the deghosted image.\n");
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }

    //--------------------------------------------
    ///Load Files
    fprintf(stderr,">| Reading... ");
    if( !readAngleFile(&angles, &finiteAngles, &perpFlags, &mu, argv[1]) )
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }
    fprintf(stderr,"Done.\n");

    if(!readPGM(&frtSpace,&readSize,&readN,argv[2],binaryFile))
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[2]);
        return EXIT_FAILURE;
    }
    size = (size_t) readSize;

    if(!readPGM(&ghost,&readN,&readN,argv[3],binaryFile))
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[3]);
        return EXIT_FAILURE;
    }
    N = (size_t) readN;

    noOfGhosts = size-mu;
    fprintf(stderr, ">| Number of Ghosts: %lu\n", noOfGhosts);

    if(N % 2 == 1) ///Assume prime if odd
        dyadicSize = FALSE;

    //--------------------------------------------
    ///Determine which finite angles are missing
    determineMissingAngles(N, finiteAngles, perpFlags, mu, &missingAngles, &missingAngleLookup);

//    printf("Missing Lookup: \n");
//    for(n = 0; n < size; n ++)
//        printf("%lu, ", missingAngleLookup[n]);

    fprintf(stderr, "\nMissing Projections at: \n");
    for(n = 0; n < noOfGhosts; n ++)
        fprintf(stderr, "%lu, ", missingAngles[n]);

    //--------------------------------------------
    ///Deghost
    fprintf(stderr,"\n>| De-Ghosting... ");
    START_TIMER;

    deghostImage = deghost(N, frtSpace, ghost, noOfGhosts);

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);

    fprintf(stderr,">| Writing Ghost Image File... ");
    writePGM(deghostImage,N,N,255,argv[4],binaryFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(perpFlags);
    free_array(finiteAngles);
    free_array(missingAngles);
    free_array(missingAngleLookup);
    free_vectorArray(angles, mu);
    free_array(frtSpace);
    free_array(deghostImage);
    free_array(ghost);

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
