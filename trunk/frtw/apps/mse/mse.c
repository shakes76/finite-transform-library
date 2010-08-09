/**
* Dyadic Fast RT
* Computes the Dyadic FRT of an image provided
* Outputs image result.
* \author Shekhar S. Chandra, 2008-9
*/
//NTTW
#include <nttw/array.h>
#include <nttw/timing.h>
#include <nttw/image.h>

//CGhosts 2 Lib
#include "noise.h"

int main(int argc, char *argv[])
{
    long *image, *result, Isum = 0;
    int rows, cols, binaryFile = FALSE;
    unsigned long long duration = 0;
    double MSE = 0;
    FILE *outFile;

    printf(">| Dyadic Fast RT Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-9\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
#if defined (NTTW_64)
    printf(">| Using 64-bit mode.\n");
#else
    printf(">| Using 32-bit mode.\n");
#endif

    if(argc != 3)
    {
        printf(">| Usage: %s <computed> <actual>\n",argv[0]);
        printf(">| computed is loaded and compared with actual.\n");
        printf(">| Files should be PGM formats.\n");
        return EXIT_FAILURE;
    }

    //--------------------------------------------
    ///Load Image
    if(!readSignedPGM(&image,&rows,&cols,argv[1],binaryFile))
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }

    if(!readSignedPGM(&result,&rows,&cols,argv[2],binaryFile))
    {
        fprintf(stderr,">| Error Opening File: %s\n",argv[2]);
        return EXIT_FAILURE;
    }
    fprintf(stderr,">| Errors to be found from %lux%lu space.\n", rows, cols);

    //--------------------------------------------
    ///Fast FRT
    fprintf(stderr,">| Calculating Errors... ");
    START_TIMER;

    MSE = mse(image,result,rows,cols);

    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    fprintf(stderr,"Done. Time: %llu usecs.\n", duration);
    fprintf(stderr,">| MSE: %f.\n", MSE);
    fprintf(stderr,">| RMSE: %f.\n", RMSE(MSE));
    fprintf(stderr,">| PSNR: %f.\n", PSNR(MSE));

    ///Output/Append times to global data file
    printf(">| Global result output to %s.\n","results.csv");
    openFile_Append("results.csv",&outFile,binaryFile);
    fprintf(outFile,"%f %f %f\n", MSE, RMSE(MSE), PSNR(MSE)); //Write initial data and time
    ///Format will be k, angle type, angle time, mt time, mt2frt time, mse, rmse, psnr
    fclose(outFile);

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(image);
    free_array(result);

    return EXIT_SUCCESS;
}
