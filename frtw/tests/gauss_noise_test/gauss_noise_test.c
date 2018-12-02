/**
* Gaussian Noise Test
* Computes the Dyadic FRT of an image provided
* Outputs image result.
* \author Shekhar S. Chandra, 2008-9
*/
#include <stdio.h>

//NTTW
#include <nttw/array.h>

//CGhosts 2 Lib
#include "noise.h"

const int Size = 128;

int main(int argc, char *argv[])
{
    nttw_integer *image;
    size_t j, k, mean, stdev;
    int binaryFile = FALSE;

    printf(">| Gaussian Noise Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 4)
    {
        printf(">| Usage: %s <mean> <stdev> <output>\n",argv[0]);
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }
    mean = atoi(argv[1]);
    stdev = atoi(argv[2]);

    //--------------------------------------------
    ///
    image = array_1D(Size*Size);

    for(j = 0; j < Size; j ++)
        for(k = 0; k < Size; k ++)
            image[j*Size+k] = (nttw_integer) Normal(mean,stdev);

    writePGM(image,Size,Size,255,argv[3],binaryFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(image);

    return EXIT_SUCCESS;
}
