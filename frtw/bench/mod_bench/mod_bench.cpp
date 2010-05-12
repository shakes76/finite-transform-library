/**
* Template project. Does nothing but take arguments
* and print them out.
* \author Shekhar S. Chandra, 2007-9
*/
#include <emmintrin.h>
//#include <xmmintrin.h>
#include <iostream>
#include <sstream>

using namespace std;

extern "C"
{
    //NTTW Library
    #include <nttw/timing.h>
    #include <nttw/number32.h>
}

///Globals

///Function Prototypes
bool processArgs(int argcount, char *argvars[], size_t &n);
static inline int round(double const x) ;
//inline nttw_integer Round(double a);
//inline nttw_integer RoundMod(double a);

int main(int argc, char *argv[])
{
    bool success = false;
    unsigned long long duration = 0, fpuTime;

    cout << ">| Benchmark FRT Algorithms Program -----------------------------" << endl;
	cout << ">| Copyright Shekhar S. Chandra, 2007-9" << endl;
	cout << ">| Machine Integer Size of " << BITS << " bits" << endl;

    ///Process arguments
    size_t N = 1;
	success = processArgs(argc,argv,N);
	if(!success)
	{
	    cout << ">| Check arguments and try again." << endl;
        return 1;
	}
	N = 1ULL << N;
	cout << ">| N: " << N << endl;

    //--------------------------------------------
    ///FPU MOD
    nttw_integer tmpVal = 3;
    nttw_real value = 1.7, recp;
    int intValue;
    //__m128d recp;
    //__m128i value64;
    //cerr << ">| Truncation test: 1.7 = " << TRUNCATE(value) << endl;

    START_TIMER;
    for(size_t j = 0; j < N; j ++)
    {
        value = tmpVal * j;
        intValue = round(value/MODULUS)*MODULUS;
        //recp = _mm_set_sd(value/MODULUS);
        //tmpVal = _mm_cvttsd_si32(recp);

        //tmpVal = FPUMOD(value,MODULUS,value/MODULUS);

        //value64 = _mm_cvtsi32_si128(tmpVal * j);
        //recp = _mm_set_sd(value/MODULUS);
        //tmpVal = _mm_cvttsd_si32(recp);
        //tmpVal *= MODULUS;

//        recp = _mm_set_sd(value/MODULUS);
//        tmpVal = _mm_cvttsd_si32(recp);
//        tmpVal *= MODULUS;
//        value -= tmpVal;

        //recp = _mm_set_sd(tmpVal);
        //recp = _mm_mul_sd(recp,modVal);
        //recp = _mm_sub_sd(recp,);

        //tmpVal *= MODULUS;
        //value -= tmpVal;
        //recp = _mm_set_sd(value);
        //tmpVal = _mm_cvttsd_si32(recp);
    }
    STOP_TIMER;
    fpuTime = MICROSECONDS_ELAPSED;
    cerr << ">| FPU Fixed Mod Time (us): " << fpuTime << endl;

    /*START_TIMER;
    for(size_t j = 0; j < N; j ++)
    {
        value = (nttw_real)tmpVal * j;
        recp = value/MODULUS;
        intValue = (nttw_integer)FPUMOD(value,MODULUS,recp);
    }
    STOP_TIMER;
    fputime = MICROSECONDS_ELAPSED;
    cerr << ">| FPU Standard Mod Time (us): " << fputime << endl;*/

    ///Mod Operator
    tmpVal = 3;
    START_TIMER;
    for(size_t j = 0; j < N; j ++)
        tmpVal = ((unsigned long long)tmpVal * j)%MODULUS;
    STOP_TIMER;
    duration = MICROSECONDS_ELAPSED;
    cerr << ">| Operator Mod Time (us): " << duration << endl;
    cerr << ">| Ratio: " << fpuTime/duration << endl;

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
        cerr << ">| Usage: " << argvars[0] << " <Power Of 2>" << endl;
        return false;
    }
}

static inline int round(double const x)
{
    //Truncate _mm_cvttsd_si32
    //Round _mm_cvtsd_si32
    return _mm_cvttsd_si32(_mm_load_sd(&x));
}

//inline nttw_integer Round(double a)
//{
//    /** x86 only
//	int retval;
//
//	__asm fld a
//	__asm fistp retval
//
//	return retval;
//	*/
//	a += typecastfix;
//	return ((nttw_integer*)&a)[IMAN] >> shiftamt;
//}
//
//inline nttw_integer RoundMod(double a)
//{
//    /** x86 only
//	int retval;
//
//	__asm fld a
//	__asm fistp retval
//
//	return retval;
//	*/
//	a += typecastfix;
//	return a - ( ((nttw_integer*)&a)[IMAN] >> shiftamt )*MODULUS;
//}
