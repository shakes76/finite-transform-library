/**
* Number Theoretic Transform
* Computes the NTT of an image provided
* Outputs image result.
* \author Shekhar S. Chandra, 2009
*/
//NTTW
#include "timing.h"
#include "prime.h"
#include "number32.h"
#include "image.h"

///Global
const int totalRows = 1000;
const int totalCols = 10;
const char *filename = "primes.csv";

int main(int argc, char *argv[])
{
    nttw_integer *image, *result, *primes, root = 0, modulus = 0, proot = 0;
    int rows, cols, binaryFile = FALSE, dyadicTransform = TRUE;
    unsigned long long duration = 0;

    printf(">| Number Theoretic Transform Program.\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %lu bits\n",BITS);
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
        ///Read primes list
        if( !readCSV(&primes,totalRows,totalCols,filename) )
        {
            printf(">| Could not open primes file.\n");
            exit(EXIT_FAILURE);
        }
        ///Compute roots
        root = findFirstPrimitiveRoot(rows);
        modulus = findAlternatePrime(primes,totalRows*totalCols,rows);
        proot = findFirstPrimitiveRoot(modulus);
        printf(">| Root: %u, Prime': %u, Primitive Root: %u.\n",root,modulus,proot);
    }

    ///Output parameters to be used
    result = array_1D(rows*cols);
    //init_1D(result,rows*cols,0);

    //--------------------------------------------
    ///Fast FRT
    fprintf(stderr,">| Transforming... ");
    START_TIMER;

    if(dyadicTransform)
        fntt_2D(image,result,rows,NTTW_FORWARD);
    else
        fntt_2D_prime(image,result,rows,root,modulus,proot,NTTW_FORWARD,TRUE);

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);

    ///Save Result
    writePGM(result,rows,cols,255,argv[2],binaryFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(image);
    free_array(result);
    if(!dyadicTransform)
        free_array(primes);

    return EXIT_SUCCESS;
}
