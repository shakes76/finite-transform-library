/**
* Number RT
* Computes the Dyadic Number FRT of an image provided
* Outputs image result.
* \author Shekhar S. Chandra, 2008-9
*/
//NTTW
#include <nttw/array.h>
#include <nttw/timing.h>
#include <nttw/image.h>

//CGhosts 2 Lib
#include "radon.h"

int main(int argc, char *argv[])
{
    nttw_integer *image, *result, Isum = 0;
    int rows, projections, cols, binaryFile = FALSE, dyadicTransform = TRUE;
    unsigned long long duration = 0;

    printf(">| Number RT Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-9\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 3)
    {
        printf(">| Usage: %s <filename> <output>\n",argv[0]);
        printf(">| filename is loaded and transformed to produce output.\n");
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }

    //--------------------------------------------
    ///Load Image
    if(!readPGM(&image,&rows,&cols,argv[1],binaryFile))
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }
    if(rows % 2 == 1) ///Odd #rows then assume prime-length
    {
        dyadicTransform = FALSE;
        projections = rows+1;
    }
    else
        projections = rows+rows/2;

    ///Output parameters to be used
    result = array_1D(projections*cols);

    //--------------------------------------------
    ///Fast FRT
    fprintf(stderr,">| Transforming... ");
    START_TIMER;

    if(dyadicTransform)
        nrt_dyadic(image,result,rows);
    else
        nrt(image,result,rows);

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);

    ///Isum
    Isum = getIsum(result,rows);
    //Isum %= MODULUS;
    fprintf(stderr,">| Isum: %u.\n",Isum);

    ///Save Result
    writePGM(result,projections,cols,255,argv[2],binaryFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(image);
    free_array(result);

    return EXIT_SUCCESS;
}
