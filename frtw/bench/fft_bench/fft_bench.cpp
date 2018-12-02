/**
* Template project. Does nothing but take arguments
* and print them out.
* \author Shekhar S. Chandra, 2007-9
*/
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

extern "C"
{
    //NTTW
    #include <nttw/timing.h>

    //C Ghosts 2 Library
    #include "array_complex.h"
    #include "fourier.h"
}

///Globals

///Function Prototypes
bool processArgs(int argcount, char *argvars[], size_t &n);
inline string numToString(double num);

int main(int argc, char *argv[])
{
    bool success = false;
    unsigned long long duration1 = 0, duration2 = 0;

    cout << ">| Benchmark FFT Algorithm Program -----------------------------" << endl;
	cout << ">| Copyright Shekhar S. Chandra, 2007-9" << endl;
	cout << ">| Machine Integer Size of " << BITS << " bits" << endl;

	///C CODE
	fftw_complex *data, *result, *inversion;
    ///END C

    ///Process arguments
    size_t N = 1;
	success = processArgs(argc,argv,N);
	if(!success)
	{
	    cout << ">| Check arguments and try again." << endl;
        return 1;
	}
	cout << ">| N: " << N << endl;
	cout << ">| Memory Alignment of arrays at " << ALIGNOF(data) << " bytes" << endl;

    //--------------------------------------------
    cerr << ">| Allocating and creating image... ";
    size_t size = N*N;
    int rank = 2, sizes[2] = {(int)N,(int)N};

	///C CODE
	data = fftw_array(size);
	result = fftw_array(size);
	inversion = fftw_array(size);
    ///END C

    for(size_t j = 0; j < N; j ++)
        for(size_t k = 0; k < N; k ++)
            data[j*N+k][0] = 1 + j*N + k;
    cerr << "Done." << endl;
    //cout << ">| Memory Alignment of values at " << ALIGNOF( data[0][0] ) << " bytes" << endl;
    cerr << ">| Data Element (0,0): " << data[0][0] << endl;

    //--------------------------------------------
    ///NFRT
    cerr << "\n>| FFT Benchmarking -------" << endl;
    cerr << ">| Forward Transform." << endl;
    ///C CODE
    START_TIMER;

    fft(rank,sizes,data,result); ///FFT

    STOP_TIMER;
    duration1 = MICROSECONDS_ELAPSED;
    ///END

    cerr << ">| Inverse Transform." << endl;
    ///C CODE
    START_TIMER;

    ifft(rank,sizes,result,inversion); ///iFFT
    for(int j = 0; j < sizes[0]; j ++) ///Norm
        for(int k = 0; k < sizes[1]; k ++)
        {
            inversion[j*N+k][0] /= N*N;
            inversion[j*N+k][1] /= N*N;
        }

    STOP_TIMER;
    duration2 = MICROSECONDS_ELAPSED;
    ///END
    cerr << ">| Transform Time Elapsed: " << duration1 << " usecs." << endl;
    cerr << ">| Inv. Transform Time Elapsed: " << duration2 << " usecs." << endl;
    cerr << ">| Inversion Element (0,0): " << inversion[0][0] << endl;

    cerr << ">| Writing Times" << endl;
    string filename = "times_fft_" + numToString(N) + ".csv";
    ofstream outFile(filename.c_str(),ios::app);

    outFile << duration1 << endl;
    outFile << duration2 << endl;

    outFile.close();

    //--------------------------------------------
    cerr << "\n>| De-Allocating ... ";
    ///C CODE
    fftw_free(data);
    fftw_free(result);
    fftw_free(inversion);
    ///END
    cerr << "Done." << endl;

    //--------------------------------------------
    cout << ">| Output Complete" << endl;
	return 0;
}

bool processArgs(int argcount, char *argvars[], size_t &n)
{
    if(argcount == 2)
    {
        stringstream arguments;

        //for(int j = 1; j < argcount-1; j ++) //-1 because dont need to cast last arg
            arguments << argvars[1] << " ";

        if (!(arguments >> n))
        {
            cerr << ">| Invalid argument N: " << argvars[1] << endl;
            return false;
        }

        return true;
    }
    else
    {
        cerr << ">| Usage: " << argvars[0] << " <N>" << endl;
        return false;
    }
}

/**
* \fn numToString(double num)
* Number to string converter
*/
inline string numToString(double num)
{
    ostringstream toString;
    if (!(toString << num))
        cerr << "ERROR: numToString(double)" << endl;
    return toString.str();
}
