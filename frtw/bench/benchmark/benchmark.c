/**
 * FRTW Benchmark Program
 * \file benchmark.c
 * \brief Transforms Benchmark Program for the NTTW/FRTW C Library.
 *
 * This file is part of FRTW Library.
 *
 * FRTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FRTW is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with FRTW.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2013
*/
#include <stdio.h>

#include <nttw/array.h>
#include <nttw/timing.h>
#include <nttw/image.h>
#include <nttw/number.h>

#include "array_complex.h"
#include "fourier.h"

///Defines to 'turn off' the transforms directly at compile time for most accurate timings
//#define FCT_BENCHMARK 1
//#define FNTT_BENCHMARK 1
#define FFT_BENCHMARK 1

int main(int argc, char *argv[])
{
    nttw_integer *data;
#ifdef FFT_BENCHMARK
    fftw_complex *dataComplex, *result, *inversion;
    nttw_big_integer *times3;
    int rank = 1, sizes[1];
#endif
    nttw_big_integer *times, *lengths, *times2;
    size_t j, k, l, N = 0, n = 1, mLower = 1, mUpper = 2, m = 1; //Use size_t instead of int for portability
    size_t count = 0, count2 = 0, count3 = 0;
    unsigned long long duration1 = 0, duration2 = 0;
    FILE *outFile = NULL;

    printf(">| Benchmark Application\n");
    printf(">| Copyright Shekhar Chandra, 2013\n");
    printf(">| Machine Integer Size of %u bits\n",BITS);
    printf(">| Memory Alignment at %u bytes\n",ALIGNOF(data));
    printf(">| Size of big integer: %u bytes\n",sizeof(nttw_big_integer));
    if(argc != 5)
    {
        printf(">| %s <n> <min dyadic> <max dyadic> <csv outputname>\n",argv[0]);
        printf(">| Computes the FCT, FNTT and FFT with n runs, from min up to max dyadic.\n");
        printf(">| Dyadic is power of 2, e.g. 14 = 2^14.\n");
        printf(">| Output CSV file with the timing results.\n");
        return EXIT_FAILURE;
    }
    n = atol(argv[1]);
    mLower = atol(argv[2]);
    mUpper = atol(argv[3]);
    m = (mUpper-mLower+1)*n;

    //-----------------------------------------
    printf(">| Allocating memory\n");
    lengths = array_1D_big(m);
    times = array_1D_big(m);
    times2 = array_1D_big(m);
#ifdef FFT_BENCHMARK
    times3 = array_1D_big(m);
#endif

    for(k = mLower, count = 0, count2 = 0, count3 = 0; k <= mUpper; k ++)
    {
        N = 1ULL << k;
        data = array_1D(N);
    #ifdef FFT_BENCHMARK
        dataComplex = fftw_array(N);
        result = fftw_array(N);
        inversion = fftw_array(N);
        sizes[0] = (int)N;
    #endif
//        printf(">| N: %llu.\n", N);

        ///FNTT
        for(j = 0; j < n; j ++)
        {
            //~ printf(">| Data: \n");
            for(l = 0; l < N; l ++)
            {
                data[l] = 1 + (nttw_integer)l;
                //~ printf(FORMAT_CSVOUTPUT_STR,data[l]);
            }
            //~ printf("\n");

            //-----------------------------------------
            //~ printf(">| FNTT TEST ----------- \n");
            START_TIMER;

            ///FNTT
            fntt(data, N, PROOT, NTTW_FORWARD);

            STOP_TIMER;
            duration1 = MICROSECONDS_ELAPSED;

            //~ printf(">| FNTT: \n");
            //~ for(j = 0; j < N; j ++)
                //~ printf(FORMAT_CSVOUTPUT_STR,data[j]);
            //~ printf("\n\n");

            RESTART_TIMER;

            ///iFNTT
            fntt(data, N, PROOT, NTTW_INVERSE);

            STOP_TIMER;
            duration2 = MICROSECONDS_ELAPSED;

            ///Norm
            ntt_norm(data,N);

//            printf(">| FNTT Took %llu usecs.\n", duration1);
//            printf(">| iFNTT Took %llu usecs.\n", duration2);
            times2[count2] = duration1 + duration2;
            count2 ++;
        }
    #ifdef FFT_BENCHMARK
        ///FFT
        for(j = 0; j < n; j ++)
        {
//            printf(">| Data: \n");
            for(l = 0; l < N; l ++)
            {
                dataComplex[l][0] = 1.0 + l;
                //~ printf(FORMAT_CSVOUTPUT_STR,data[l]);
            }
            //~ printf("\n");

            //-----------------------------------------
            //~ printf(">| FFT TEST ----------- \n");
            START_TIMER;

            ///FFT
            fft(rank,sizes,dataComplex,result);

            STOP_TIMER;
            duration1 = MICROSECONDS_ELAPSED;

            //~ printf(">| FFT: \n");
            //~ for(j = 0; j < N; j ++)
                //~ printf(FORMAT_CSVOUTPUT_STR,data[j]);
            //~ printf("\n\n");

            RESTART_TIMER;

            ///iFFT
            ifft(rank,sizes,result,inversion);

            STOP_TIMER;
            duration2 = MICROSECONDS_ELAPSED;

            ///Norm
            fft_norm(inversion,N);

//            printf(">| FFT Took %llu usecs.\n", duration1);
//            printf(">| iFFT Took %llu usecs.\n", duration2);
            times3[count3] = duration1 + duration2;
            count3 ++;
        }
    #endif

//        printf(">| FCT (%llu runs) for size %llu Took %llu usecs.\n", N, times[k]);

        free_array(data);
    #ifdef FFT_BENCHMARK
        fftw_free(dataComplex);
        fftw_free(result);
        fftw_free(inversion);
        fftw_cleanup();
    #endif
    }

    printf(">| Writing Benchmark Times\n");
    if( !openFile_Write(argv[4],&outFile,FALSE) )
        return FALSE;

    for (count = 0; count < m; count ++)
    #ifdef FFT_BENCHMARK
        fprintf(outFile,"%llu, %llu, %llu, %llu\n",lengths[count],times[count],times2[count],times3[count]);
    #else
        fprintf(outFile,"%llu, %llu, %llu\n",lengths[count],times[count],times2[count]);
    #endif
    fclose(outFile);

    printf(">| Operation Complete\n");
    free_array(times);
    free_array(lengths);
    free_array(times2);
#ifdef FFT_BENCHMARK
    free_array(times3);
#endif
    return EXIT_SUCCESS;
}
