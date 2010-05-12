/**
* Benchmark FNTT Algorithm
* \author Shekhar S. Chandra, 2007-9
*/
#include <iostream>
#include <sstream>

using namespace std;

extern "C"
{
    //C NTTW Library
    #include "timing.h"
    #include "number32.h"
    #include "array.h"
    #include "image.h"
}

///Globals
const string filename = "fntt_bench.pgm";

///Function Prototypes
bool processArgs(int argcount, char *argvars[], size_t &n);

int main(int argc, char *argv[])
{
    bool success = false;
    unsigned long long duration1 = 0, duration2 = 0;

    ///C CODE
	nttw_integer *data, *result, *inversion;
    ///END C

    cout << ">| Benchmark FNTT Algorithm Program -----------------------------" << endl;
	cout << ">| Copyright Shekhar S. Chandra, 2007-9" << endl;
	cout << ">| Machine Integer Size of " << BITS << " bits" << endl;
	cout << ">| Largest Library Integer Size of " << 8*sizeof(nttw_integer) << " bits" << endl;
	cout << ">| Array memory alignment at " << ALIGNOF(data) << endl;

    ///Process arguments
    size_t N = 1;
	success = processArgs(argc,argv,N);
	if(!success)
	{
	    cout << ">| Check arguments and try again." << endl;
        return 1;
	}
	cout << ">| N: " << N << endl;

    //--------------------------------------------
    cout << ">| Allocating and creating image... ";
    size_t size = N*N;
	///C CODE
	data = array_1D(size);
	result = array_1D(size);
	inversion = array_1D(size);
    ///END C

    for(size_t j = 0; j < N; j ++)
        for(size_t k = 0; k < N; k ++)
            data[j*N+k] = ( 1 + j*N + k )%MODULUS;
    cout << "Done." << endl;
    cout << ">| Data Element (0,0): " << data[0] << endl;

    //--------------------------------------------
    ///NFRT
    cout << "\n>| FNTT Benchmarking -------" << endl;
    START_TIMER;
    ///C CODE
    fntt_2D(data,result,N,NTTW_FORWARD);
    ///END
    STOP_TIMER;
    duration1 = MICROSECONDS_ELAPSED;

    START_TIMER;
    ///C CODE
    fntt_2D(result,inversion,N,NTTW_INVERSE);
    ///END
    STOP_TIMER;
    duration2 = MICROSECONDS_ELAPSED;
    cout << ">| Transform Time Elapsed: " << duration1 << " usecs." << endl;
    cout << ">| Inv. Transform Time Elapsed: " << duration2 << " usecs." << endl;
    cout << ">| Inversion Element (0,0): " << inversion[0] << endl;

    cout << ">| Writing PGM: " << filename << endl;
    writePGM(data,N,N,255,filename.c_str(),FALSE);

    //--------------------------------------------
    cout << "\n>| De-Allocating ... ";
    ///C CODE
    free_array(data);
    free_array(result);
    free_array(inversion);
    ///END
    cout << "Done." << endl;

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
            cout << ">| Invalid argument N: " << argvars[1] << endl;
            return false;
        }

        return true;
    }
    else
    {
        cout << ">| Usage: " << argvars[0] << " <N>" << endl;
        return false;
    }
}
