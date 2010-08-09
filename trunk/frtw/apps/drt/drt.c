/**
* Dyadic Discrete RT
* Computes the Dyadic FRT of an image provided
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
    int rows, cols, binaryFile = FALSE;
    unsigned long long duration = 0;

    printf(">| Dyadic Discrete RT Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");
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

    ///Output parameters to be used
    result = array_1D((rows+rows/2)*cols);

    //--------------------------------------------
    ///Fast FRT
    fprintf(stderr,">| Transforming... ");
    START_TIMER;

    drt_dyadic(image,result,rows);

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);

    ///Isum
    Isum = getIsum(result,rows);
    fprintf(stderr,">| Isum: %u.\n",Isum);

    ///Save Result
    writePGM(result,rows+rows/2,cols,255,argv[2],binaryFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(image);
    free_array(result);

    return EXIT_SUCCESS;
}
