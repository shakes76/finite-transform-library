/**
 * NTTW Image Module
 * \file image.h
 * \brief Image Source for the NTTW C Library.
 *
 * This file implements the functions for input and output of image files.
 *
 * This file is part of NTTW Library.
 *
 * NTTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NTTW is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NTTW. If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2008-9
*/
#include "image.h"

int openFile_Read(const char *filename, FILE **inFilePtr, int binary)
{
    if (binary)
    {
        printf("Opening File: %s as binary.\n", filename);
        *inFilePtr = fopen(filename, "rb"); //!< Open Binary File.
    }
    else
    {
        printf("Opening File: %s as ASCII.\n", filename);
        *inFilePtr = fopen(filename, "r"); //!< Open ASCII File.
    }

    if (*inFilePtr == NULL)
    {
        fprintf(stderr,"ERROR: %s could not be opened.\n",filename);
        return FALSE;
    }

    return TRUE;
}

int openFile_Write(const char *filename, FILE **inFilePtr, int binary)
{
    if (binary)
    {
        printf("Opening File: %s as binary.\n", filename);
        *inFilePtr = fopen(filename, "wb"); //!< Open Binary File.
    }
    else
    {
        printf("Opening File: %s as ASCII.\n", filename);
        *inFilePtr = fopen(filename, "w"); //!< Open ASCII File.
    }

    if (*inFilePtr == NULL)
    {
        fprintf(stderr,"ERROR: %s could not be opened.\n",filename);
        return FALSE;
    }

    return TRUE;
}

int readCSV(nttw_integer **data, const int rows, const int cols, const char *filename)
{
    int j, k, success, size, binary = FALSE;
    nttw_integer *tmptr;
    FILE *inFile = NULL;

    if( !openFile_Read(filename,&inFile,binary) )
        return FALSE;

    ///Allocate Memory for Image appropriately
    size = rows*cols;
    *data = array_1D(size);

    ///Read Data according to size and type provided in header
    tmptr = *data;
    for(j = 0; j < rows; j ++)
    {
        for(k = 0; k < cols; k ++)
            success = fscanf(inFile, FORMAT_CSVOUTPUT_STR, tmptr++);
        success = fscanf(inFile, "\n");
    }

    fclose(inFile); ///Close File
    return TRUE;
}

int readPGMHeader(FILE *inFile, char *type, int *rows, int *cols, int *maxGrey)
{
    char a, b, c[100];
    char *tmpStr;
    int width, height, greymax, success;

    ///Read PGM Header
    ///Read Filetype and File comment
    success = fscanf(inFile, "%c", &a);
    success = fscanf(inFile, "%c", &b);
    success = fscanf(inFile, "\n");
    tmpStr = fgets(c, 100, inFile);
    //printf("File Info: %s\n",a,b,c);

    if(c[0] != '#') ///See if Comment not present
    {
        rewind(inFile); ///If so, Start again and dont read comment
        ///Read Filetype
        success = fscanf(inFile, "%c", &a);
        success = fscanf(inFile, "%c", &b);
        success = fscanf(inFile, "\n");
        printf("File Type: %c%c\n",a,b);
    }
    else ///Else read comments until none left
    {
        printf("File Info: %s",c);
        //while(c[0] == '#')
            //tmpStr = fgets(c,100,inFile);
    }

    ///Read Image Characteristics
    success = fscanf(inFile, "%i", &width);
    success = fscanf(inFile, "%i", &height);
    success = fscanf(inFile, "%i", &greymax);
    success = fscanf(inFile, "\n");
    printf("File Size: %dx%d, Grey Scale Max: %d\n",width, height, greymax);

    ///Allocate Memory for Image appropriately
    *type = b;
    *rows = height;
    *cols = width;
    *maxGrey = greymax;

    return TRUE;
}

int readPGM(nttw_integer **data, int *rows, int *cols, const char *filename, int binary)
{
    char type;
    int width, height, greymax, size, success;
    nttw_integer *tmptr;
    FILE *inFile = NULL;

    if( !openFile_Read(filename,&inFile,binary) )
        return FALSE;

    readPGMHeader(inFile,&type,&height,&width,&greymax);

    ///Allocate Memory for Image appropriately
    size = width*height;
    *rows = height;
    *cols = width;
    *data = array_1D(size);

    ///Read Data according to size and type provided in header
    tmptr = *data;
    if ((type=='2') || (type=='5'))
    {
        while (!feof(inFile) && size--)
        {
            success = fscanf(inFile, FORMAT_INPUT_STR, tmptr++);
        }
    }
    else
        if (type=='3')
        {
            while (!feof(inFile) && size--)
            {
                success = fscanf(inFile, FORMAT_INPUT_STR, tmptr);
                success = fscanf(inFile, FORMAT_INPUT_STR, tmptr);
                success = fscanf(inFile, FORMAT_INPUT_STR, tmptr++);
            }
        }

    fclose(inFile); ///Close File
    return TRUE;
}

int readSignedPGM(long **data, int *rows, int *cols, const char *filename, int binary)
{
    char type;
    int width, height, greymax, size, success;
    long *tmptr;
    FILE *inFile = NULL;

    if( !openFile_Read(filename,&inFile,binary) )
        return FALSE;

    readPGMHeader(inFile,&type,&height,&width,&greymax);

    ///Allocate Memory for Image appropriately
    size = width*height;
    *rows = height;
    *cols = width;
    *data = arraySigned_1D(size);

    ///Read Data according to size and type provided in header
    tmptr = *data;
    if ((type=='2') || (type=='5'))
    {
        while (!feof(inFile) && size--)
        {
            success = fscanf(inFile, "%li", tmptr++);
        }
    }
    else
        if (type=='3')
        {
            while (!feof(inFile) && size--)
            {
                success = fscanf(inFile, "%li", tmptr);
                success = fscanf(inFile, "%li", tmptr);
                success = fscanf(inFile, "%li", tmptr++);
            }
        }

    fclose(inFile); ///Close File
    return TRUE;
}

int writePGM(nttw_integer *data, const int rows, const int cols, const int greyMax, const char *filename, int binary)
{
    int j, k, count = 0;
    FILE *outFile = NULL;

    if( !openFile_Write(filename,&outFile,binary) )
        return FALSE;

    fprintf(outFile,"P2\n# Generated PGM.\n");
    fprintf(outFile,"%d %d\n%d\n",cols,rows,greyMax);

    for (j = 0; j < rows; j ++)
    {
        for (k = 0; k < cols; k ++)
        {
            fprintf(outFile,FORMAT_OUTPUT_STR,data[count]);
            count ++;
        }
        fprintf(outFile,"\n");
    }

    fclose(outFile);
    return TRUE;
}

int writeSignedPGM(long *data, const int rows, const int cols, const int greyMax, const char *filename, int binary)
{
    int j, k, count = 0;
    FILE *outFile = NULL;

    if( !openFile_Write(filename,&outFile,binary) )
        return FALSE;

    fprintf(outFile,"P2\n# Generated PGM.\n");
    fprintf(outFile,"%d %d\n%d\n",cols,rows,greyMax);

    for (j = 0; j < rows; j ++)
    {
        for (k = 0; k < cols; k ++)
        {
            fprintf(outFile,"%li ",data[count]);
            count ++;
        }
        fprintf(outFile,"\n");
    }

    fclose(outFile);
    return TRUE;
}

//Operations
long* cropImage(long *initImage, const int initRows, const int initCols, const int newRows, const int newCols)
{
    int j, k, startx, endx, starty, endy;
    long *newImage = NULL;

    ///Set bounds of new image to extract
    startx = initRows/2 - newRows/2;
    starty = initCols/2 - newCols/2;
    if(initRows != newRows)
        endx = initRows/2 + newRows/2;
    else
        endx = initRows-1;
    if(initCols != newCols)
        endy = initCols/2 + newCols/2;
    else
        endy = initCols-1;

    if(startx < 0 || starty < 0 || endx > initRows || endy > initCols)
    {
        fprintf(stderr,"ERROR: Out of Bounds in Crop Image Function.\n");
        exit(EXIT_FAILURE);
    }

    ///Allocate new image memory
    newImage = arraySigned_1D((endx-startx+1)*(endy-starty+1));

    ///Extract relevant part of image to form cropped image.
    for(j = startx; j <= endx; j ++)
        for(k = starty; k <= endy; k ++)
            newImage[(j-startx)*newCols + (k-starty)] = initImage[j*initCols + k];

    return newImage;
}

long* truncateImage(long *initImage, const int initRows, const int initCols, const int newRows, const int newCols)
{
    int j, k;
    long *newImage = NULL;

    ///Allocate new image memory
    newImage = arraySigned_1D(newRows*newCols);

    ///Extract relevant part of image to form cropped image.
    for(j = 0; j < newRows; j ++)
        for(k = 0; k < newCols; k ++)
            newImage[j*newCols + k] = initImage[j*initCols + k];

    return newImage;
}

long* embedImage(long *initImage, const int initRows, const int initCols, const int newRows, const int newCols)
{
    int j, k, tmpRows, tmpCols, imgRows, imgCols;
    long *newImage = NULL;

    ///Ensure to only reference inside arrays
    if(newRows > initRows)
    {
        tmpRows = newRows;
        imgRows = initRows; //Safe to use image rows
    }
    else
    {
        tmpRows = initRows;
        imgRows = newRows;
    }

    if(newCols > initCols)
    {
        tmpCols = newCols;
        imgCols = initCols; //Safe to use image cols
    }
    else
    {
        tmpCols = initCols;
        imgCols = newCols;
    }

    newImage = arraySigned_1D(tmpRows*tmpCols);

    for(j = 0; j < imgRows; j ++)
        for(k = 0; k < imgCols; k ++)
            newImage[j*tmpCols+k] = initImage[j*imgCols+k];

    return newImage;
}

long* diffImage(long *image1, long *image2, const int rows, const int cols)
{
    int j, k;
    long *result = NULL;

    result = arraySigned_1D(rows*cols);

    for(j = 0; j < rows; j ++)
        for(k = 0; k < cols; k ++)
            result[j*cols+k] = image1[j*cols+k] - image2[j*cols+k];

    return result;
}
