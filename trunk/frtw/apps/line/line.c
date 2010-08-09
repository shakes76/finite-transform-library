/**
* Finite Line Study Program
* Computes the finite line of given parameters
* Outputs Diophatine File of result.
* \author Shekhar S. Chandra, 2008-10
*/
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    size_t x, N = 0, y = 0, m = 0, t = 0;
    FILE *outFile;

    printf(">| Finite Line Study Program.\n");
    printf(">| Copyright Shekhar Chandra, 2008-10\n");

    if(argc != 5)
    {
        printf(">| Usage: %s <m> <t> <N> <output>\n",argv[0]);
        printf(">| A line of mx + t (mod N) is drawn.\n");
        printf(">| Files output is in Diophantine format for DGV.\n");
        return EXIT_FAILURE;
    }
    m = atoi(argv[1]);
    t = atoi(argv[2]);
    N = atoi(argv[3]);

    //--------------------------------------------
    ///Open outFile, writing dio file
    outFile = fopen(argv[4],"w");

    ///Write lattice parameters
    fprintf(outFile,"%lu %lu %lu %lu %lu #Primary Lattice\n",2*N,2*N,0,0,0);
    fprintf(outFile,"%lu %lu %lu %lu %lu #Secondary Lattice\n",N,N,255,255,255);
    fprintf(outFile,"%lu %lu %lu #Background Colour\n",26,52,104);

    ///Compute line
    for(x = 0; x < N; x ++)
    {
        y = (m*x + t)%N; ///Determine point
        fprintf(outFile,"%lu %lu %lu %lu %lu #Point\n",x,y,192,0,0);
    }

    //--------------------------------------------
    ///Cleanup
    printf(">| Operation Complete\n");
    fclose(outFile);

    return EXIT_SUCCESS;
}
