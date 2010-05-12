/**
* Benchmark FFRT Algorithm
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
    #include "radon.h"
}

///Globals

///Function Prototypes
bool processArgs(int argcount, char *argvars[], size_t &n);
inline string numToString(double num);

int main(int argc, char *argv[])
{
    bool success = false;
    unsigned long long duration1 = 0, duration2 = 0;

    cout << ">| Benchmark FFRT Algorithm Program -----------------------------" << endl;
	cout << ">| Copyright Shekhar S. Chandra, 2007-9" << endl;
	cout << ">| Machine Integer Size of " << BITS << " bits" << endl;

	///C CODE
	nttw_integer *data, *result, *inversion;
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
    cout << ">| Allocating and creating image... ";
    size_t size = N*N*sizeof(nttw_integer);
	///C CODE
	data = array_1D(size);
	result = array_1D(size);
	inversion = array_1D(size);
    ///END C

    for(size_t j = 0; j < N; j ++)
        for(size_t k = 0; k < N; k ++)
            data[j*N+k] = 1 + j*N + k;
    cout << "Done." << endl;
    cout << ">| Data Element (0,0): " << data[0] << endl;

    //--------------------------------------------
    ///NFRT
    cout << "\n>| FFRT Benchmarking -------" << endl;
    cout << ">| Forward Transform." << endl;
    ///C CODE
    START_TIMER;

    frt_dyadic(data,result,N);

    STOP_TIMER;
    duration1 = MICROSECONDS_ELAPSED;
    ///END

    cout << ">| Inverse Transform." << endl;
    ///C CODE
    START_TIMER;

    ifrt_dyadic(result,inversion,N,TRUE);

    STOP_TIMER;
    duration2 = MICROSECONDS_ELAPSED;
    ///END
    cout << ">| Transform Time Elapsed: " << duration1 << " usecs." << endl;
    cout << ">| Inv. Transform Time Elapsed: " << duration2 << " usecs." << endl;
    cout << ">| Inversion Element (0,0): " << inversion[0] << endl;

    cerr << ">| Writing Times" << endl;
    string filename = "times_ffrt_" + numToString(N) + ".csv";
    ofstream outFile(filename.c_str(),ios::app);

    outFile << duration1 << endl;
    outFile << duration2 << endl;

    outFile.close();

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
