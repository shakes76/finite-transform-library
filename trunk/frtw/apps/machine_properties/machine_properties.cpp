/**
* Machine Datatype Properties Program
* Prints the sizes and alignment of types for the machine.
* \author Shekhar S. Chandra, 2007-9
*/
#include <iostream>
#include <sstream>

using namespace std;

extern "C"
{
    //C Ghosts 2 Library
    #include <nttw/array.h>
}

///Globals

///Function Prototypes

int main(int argc, char *argv[])
{
    cout << ">| Machine Datatype Properties Program -----------------------------" << endl;
	cout << ">| Copyright Shekhar S. Chandra, 2007-9" << endl;
	cout << ">| Machine Integer Size of " << BITS << " bits" << endl;

    //--------------------------------------------
    CACHE_ALIGN nttw_integer value;
    CACHE_ALIGN nttw_real dblValue;
    CACHE_ALIGN nttw_integer *arrayPtr;

    arrayPtr = array_1D(2);

    cout << "\n>| Library Types ------" << endl;
    cout << ">| Memory Alignment of Library Type at " << ALIGNOF(value) << " bytes" << endl;
    cout << ">| Memory Alignment of Complex Library Type at " << ALIGNOF(dblValue) << " bytes" << endl;
    cout << ">| Memory Alignment of Array of Library Type at " << ALIGNOF(arrayPtr) << " bytes" << endl;
    cout << "\n>| Floating Point Types ------" << endl;
    cout << ">| Size of float:\t" << sizeof(float) << ",\t Alignment:\t" << ALIGNOF(float) << " bytes" << endl;
    cout << ">| Size of double:\t" << sizeof(double) << ",\t Alignment:\t" << ALIGNOF(double) << " bytes" << endl;
    cout << ">| Size of long double:\t" << sizeof(long double) << ",\t Alignment:\t" << ALIGNOF(long double) << " bytes" << endl;
    cout << "\n>| Signed Types ------" << endl;
    cout << ">| Size of bool:\t" << sizeof(bool) << ",\t Alignment:\t" << ALIGNOF(bool) << " bytes" << endl;
    cout << ">| Size of char:\t" << sizeof(char) << ",\t Alignment:\t" << ALIGNOF(char) << " bytes" << endl;
    cout << ">| Size of short:\t" << sizeof(short) << ",\t Alignment:\t" << ALIGNOF(short) << " bytes" << endl;
    cout << ">| Size of int:\t\t" << sizeof(int) << ",\t Alignment:\t" << ALIGNOF(int) << " bytes" << endl;
    cout << ">| Size of long:\t" << sizeof(long) << ",\t Alignment:\t" << ALIGNOF(long) << " bytes" << endl;
#if defined WIN32
        cout << ">| Size of __int64:\t" << sizeof(__int64) << ",\t Alignment:\t" << ALIGNOF(__int64) << " bytes" << endl;
#endif
    cout << ">| Size of long long:\t" << sizeof(long long) << ",\t Alignment:\t" << ALIGNOF(long long) << " bytes" << endl;
    cout << "\n>| Unsigned Types ------" << endl;
    cout << ">| Size of unsigned char:\t" << sizeof(unsigned char) << ",\t Alignment:\t" << ALIGNOF(unsigned char) << " bytes" << endl;
    cout << ">| Size of unsigned short:\t" << sizeof(unsigned short) << ",\t Alignment:\t" << ALIGNOF(unsigned short) << " bytes" << endl;
    cout << ">| Size of unsigned int:\t" << sizeof(unsigned int) << ",\t Alignment:\t" << ALIGNOF(unsigned int) << " bytes" << endl;
    cout << ">| Size of unsigned long:\t" << sizeof(unsigned long) << ",\t Alignment:\t" << ALIGNOF(unsigned long) << " bytes" << endl;
    cout << ">| Size of unsigned long long:\t" << sizeof(unsigned long long) << ",\t Alignment:\t" << ALIGNOF(unsigned long long) << " bytes" << endl;
    cout << ">| Size of size_t:\t\t" << sizeof(size_t) << ",\t Alignment:\t" << ALIGNOF(size_t) << " bytes" << endl;

    //--------------------------------------------
    free_array(arrayPtr);
    cout << ">| Output Complete" << endl;
	return 0;
}
