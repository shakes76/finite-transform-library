/**
 * NTTW Image Module
 * \file image.h
 * \brief Image Header/Object for the NTTW C Library.
 *
 * This file defines the functions for input and output of image files.
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
#ifndef IMAGE_H_INCLUDED
#define IMAGE_H_INCLUDED

#include <stdio.h>

#include "array.h"

/**
 * \defgroup Imaging Image I/O and Helpers
 * \brief The input and output of PGM and CSV files.
 *
 * Read a PGM file by using readPGM() and write using writePGM(). CSV files can be read using readCSV().
 * All functions read and write 32-bit integers.
 */
//@{
/**
 * \fn openFile_Read(const char *filename, FILE **inFilePtr, int binary)
 * \brief Opens the file for read mode.
*/
NTTW_DLL_SYM int openFile_Read(const char *filename, FILE **inFilePtr, int binary);
/**
 * \fn openFile_Write(const char *filename, FILE **inFilePtr, int binary)
 * \brief Opens the file for write mode.
*/
NTTW_DLL_SYM int openFile_Write(const char *filename, FILE **inFilePtr, int binary);

/**
 * \fn readCSV(nttw_integer **data, const int rows, const int cols, const char *filename)
 * \brief Reads a Comma Separated Value (CSV) file. Memory is internally allocated to the correct size.
*/
NTTW_DLL_SYM int readCSV(nttw_integer **data, const int rows, const int cols, const char *filename);

/**
 * \fn readPGMHeader(FILE *inFile, char *type, int *rows, int *cols, int *maxGrey)
 * \brief Reads the parameters in the header of a PGM file.
*/
NTTW_DLL_SYM int readPGMHeader(FILE *inFile, char *type, int *rows, int *cols, int *maxGrey);

/**
 * \fn readPGM(nttw_integer **data, int *rows, int *cols, const char *filename, int binary)
 * \brief Reads in an image which has a PGM format of name given.
 *
 * Takes a reference to a 1D pointer, and reads the file as a 2D array into the 1D array.
 *
 * Memory is handled internally, no prior allocation is required.
 * \param data - This is an array passed by reference, or &(image)
 * where "image" is a 1D array declared with array_1D().
 * \param rows - number of rows that the file has. Passed by reference.
 * Prior knowledge of file size is not required. This will contain the
 * number of rows in file.
 * \param cols - number of columns that the file has. Passed by reference.
 * Prior knowledge of file size is not required. This will contain the
 * number of columns in file.
 * \param filename - The name of the PGM to be opened.
 * \param binary - Boolean for reading binary file.
 * \return TRUE if read is successful.
*/
NTTW_DLL_SYM int readPGM(nttw_integer **data, int *rows, int *cols, const char *filename, int binary);

/**
 * \fn readSignedPGM(long **data, int *rows, int *cols, const char *filename, int binary)
 * \brief Inputs data as a PGM file with signed values.
 * \return TRUE if successful.
*/
NTTW_DLL_SYM int readSignedPGM(long **data, int *rows, int *cols, const char *filename, int binary);

/**
 * \fn writePGM(nttw_integer *data, const int rows, const int cols, const int greyMax, const char *filename, int binary)
 * \brief Outputs data as a PGM file.
 * \return TRUE if successful.
*/
NTTW_DLL_SYM int writePGM(nttw_integer *data, const int rows, const int cols, const int greyMax, const char *filename, int binary);

/**
 * \fn writeSignedPGM(long *data, const int rows, const int cols, const int greyMax, const char *filename, int binary)
 * \brief Outputs data as a PGM file with signed values.
 * \return TRUE if successful.
*/
NTTW_DLL_SYM int writeSignedPGM(long *data, const int rows, const int cols, const int greyMax, const char *filename, int binary);
//@}

/**
* \defgroup Imaging_Ops Image Operations
*/
//@{
/**
 * \fn cropImage(long *initImage, const int initRows, const int initCols, const int newRows, const int newCols)
 * \brief Crops the image to the size provided from the centre. EXPERIMENTAL
*/
NTTW_DLL_SYM long* cropImage(long *initImage, const int initRows, const int initCols, const int newRows, const int newCols);
/**
 * \fn truncateImage(long *initImage, const int initRows, const int initCols, const int newRows, const int newCols)
 * \brief Crops the image to the size provided from the top left-hand side. EXPERIMENTAL
*/
NTTW_DLL_SYM long* truncateImage(long *initImage, const int initRows, const int initCols, const int newRows, const int newCols);
/**
 * \fn embedImage(long *initImage, const int initRows, const int initCols, const int newRows, const int newCols)
 * \brief Embeds the image to the size provided from the top left-hand corner. EXPERIMENTAL
*/
NTTW_DLL_SYM long* embedImage(long *initImage, const int initRows, const int initCols, const int newRows, const int newCols);
/**
 * \fn diffImage(long *image1, long *image2, const int rows, const int cols)
 * \brief Takes the differences of image1 and image2 and returns pointer to result. EXPERIMENTAL
*/
NTTW_DLL_SYM long* diffImage(long *image1, long *image2, const int rows, const int cols);
//@}
#endif // IMAGE_H_INCLUDED
