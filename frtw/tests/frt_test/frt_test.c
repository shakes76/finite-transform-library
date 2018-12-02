/**
* Test for the Fast Radon Transform.
* \author Shekhar S. Chandra, 2008-9
*/
#include <nttw/array.h>
#include <nttw/number32.h>
#include <nttw/image.h>

#include "radon.h"

const char *filename = "frt.pgm";

int main(int argc, char *argv[])
{
    nttw_integer *data, *result, *inversion;
    size_t j, k;
    size_t N = 0, size;
    int normTransform = TRUE, outputToStd = TRUE;

    printf(">| Fast RT Test\n");
    printf(">| Copyright Shekhar Chandra, 2009\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(data));
    if(argc != 2)
    {
        printf(">| %s <n>\n",argv[0]);
        return EXIT_FAILURE;
    }
    N = atoi(argv[1]);

    if(N > 32)
        outputToStd = FALSE;

    //-----------------------------------------
    printf(">| Allocating memory\n");
    size = N*N;
    data = array_1D(size);
    result = array_1D( N*(N+N/2) );
    inversion = array_1D(size);

    printf(">| Data: \n");
    for(j = 0; j < N; j ++)
        for(k = 0; k < N; k ++)
            data[j*N+k] = 1 + j*(nttw_integer)N + k;

if(outputToStd)
{
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
            printf(FORMAT_CSVOUTPUT_STR,data[j*N+k]);
        printf("\n");
    }
}
else
    printf(">| (0,0): %u,\t(0,1): %u\n",data[0],data[1]);

    //-----------------------------------------
    ///FNT
    printf(">| FRT TEST ----------- \n");
    printf(">| Forward FRT \n");
    frt_dyadic(data,result,N);

    printf(">| FRT: \n");
if(outputToStd)
{
    for(j = 0; j < N+N/2; j ++)
    {
        for(k = 0; k < N; k ++)
            printf(FORMAT_CSVOUTPUT_STR,result[j*N+k]);
        printf("\n");
    }
}
else
    printf(">| (0,0): %u,\t(0,1): %u\n",result[0],result[1]);

    ///iFNT
    printf(">| Inverse FRT \n");
    ifrt_dyadic(result, inversion, N, normTransform);

    printf(">| iFRT: \n");
if(outputToStd)
{
    for(j = 0; j < N; j ++)
    {
        for(k = 0; k < N; k ++)
            printf(FORMAT_CSVOUTPUT_STR,inversion[j*N+k]);
        printf("\n");
    }
}
else
{
    printf(">| (0,0): %u,\t(0,1): %u\n",inversion[0],inversion[1]);
    printf(">| Writing PGM: %s\n",filename);
    writePGM(inversion,N,N,255,filename,FALSE);
}

    //-----------------------------------------
    printf(">| Operation Complete\n");
    free_array(data);
    free_array(result);
    free_array(inversion);
    return EXIT_SUCCESS;
}
