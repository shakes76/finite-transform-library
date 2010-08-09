/**
* Mojette Transform
* Computes the MT of an image provided using angles given
* Outputs data result.
* \author Shekhar S. Chandra, 2008-10
*/
//NTTW
#include <nttw/array.h>
#include <nttw/timing.h>
#include <nttw/image.h>

//CGhosts 2 Lib
#include "array_complex.h"
#include "mojette.h"

///Prototypes
int readAngleFile(vector **angles, size_t *projNo, const char* filename);

int main(int argc, char *argv[])
{
    nttw_integer *image;
    mojetteProjection *result;
    vector *angles;
    int rows, cols, binaryFile = FALSE;
    size_t mu;
    unsigned long long duration = 0;
    FILE *outFile;

    printf(">| Mojette Transform Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 4)
    {
        printf(">| Usage: %s <filename> <angles> <output>\n",argv[0]);
        printf(">| filename is loaded and transformed to produce output.\n");
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }

    //--------------------------------------------
    ///Load Image
    if( !readPGM(&image,&rows,&cols,argv[1],binaryFile) )
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }

    if( !readAngleFile(&angles, &mu, argv[2]) )
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[2]);
        return EXIT_FAILURE;
    }

    ///Output parameters to be used
    result = mojetteArray_1D(angles,mu,rows,cols);

    //--------------------------------------------
    ///Fast FRT
    fprintf(stderr,">| Transforming... ");
    START_TIMER;

    mt(angles,image,result,mu,rows,cols);

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);

    fprintf(stderr,">| Writing File... ");
    ///Save Result
    writeMojetteSet(result,mu,argv[3]);
    fprintf(stderr,"Done.\n");

    ///Output/Append times to global data file
    printf(">| Global result output to %s.\n","results.csv");
    openFile_Append("results.csv",&outFile,binaryFile);
    fprintf(outFile,"%llu ", duration); //Write initial data and time
    ///Format will be k, angle type, angle time, mt time, mt2frt time, mse, rmse
    fclose(outFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(image);
    free_mojetteArray(result,mu);
    free_vectorArray(angles, mu);

    return EXIT_SUCCESS;
}

int readAngleFile(vector **angles, size_t *projNo, const char* filename)
{
    int j, success = FALSE;
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

    for(j = 0; j < *projNo; j ++)
    {
        success = fscanf(inFile,"%lu %li %li %lu",&m, &p, &q, &flag);
        setXY( (*angles)[j], p, q );
    }

    for(j = 0; j < *projNo; j ++)
        printf( "%li/%li, ", getY((*angles)[j]), getX((*angles)[j]) );
    printf("\n");

    fclose(inFile);
    return TRUE;
}
