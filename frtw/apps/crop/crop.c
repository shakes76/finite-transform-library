/**
* Crop/Embed
* Crops/Embeds images
* Outputs image result.
* \author Shekhar S. Chandra, 2008-9
*/
//NTTW
#include <nttw/array.h>
#include <nttw/timing.h>
#include <nttw/image.h>

int main(int argc, char *argv[])
{
    long *image, *result;
    int binaryFile = FALSE, rows, cols, newRows, newCols, topLeft = FALSE;
    unsigned long long bits = 8*sizeof(nttw_integer); //1 << x is the x Power of 2

    printf(">| Crop Image Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");
    printf(">| Integer Size of %llu bits\n",(unsigned long long)bits);

    if(argc != 6)
    {
        printf(">| Usage: %s <filename> <rows> <cols> <top left?> <output>\n",argv[0]);
        printf(">| filename is loaded and cropped to produce output.\n");
        printf(">| If new size is larger, the original image is embedded.\n");
        printf(">| Top Left is a Boolean and is used whereever possible.\n");
        printf(">| Files should be PGM formats and rows and cols are the new sizes.\n");
        return EXIT_FAILURE;
    }
    newRows = atoi(argv[2]);
    newCols = atoi(argv[3]);
    topLeft = atoi(argv[4]);

    ///Load Image
    if(!readSignedPGM(&image,&rows,&cols,argv[1],binaryFile))
    {
        printf(">| Error Opening File: %s\n",argv[1]);
        return EXIT_FAILURE;
    }

    ///Crop Image if new rows/cols smaller
    if(newRows <= rows && newCols <= cols)
    {
        printf(">| Smaller size provided, image will be cropped.\n");
        if(!topLeft)
            result = cropImage(image,rows,cols,newRows,newCols);
        else
            result = truncateImage(image,rows,cols,newRows,newCols);
    }
    else ///Embed Image if rows/cols larger
    {
        printf(">| Larger size provided, image will be embedded.\n");
        result = embedImage(image,rows,cols,newRows,newCols);
    }

    ///Save Result
    writeSignedPGM(result,newRows,newCols,255,argv[5],binaryFile);

    ///Cleanup
    printf(">| Operation Complete\n");
    free_array(image);
    free_array(result);
    return EXIT_SUCCESS;
}
