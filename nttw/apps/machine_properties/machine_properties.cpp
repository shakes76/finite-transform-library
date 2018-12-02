/**
 * Machine Datatype Properties Program
 * Prints the sizes and alignment of types for the machine.
 * NTTW Machine Test Program
 * \file machine_properties.cpp
 * \brief Machine Test Program for the NTTW C Library.
 *
 * This file is part of NTTW Library.
 *
 * NTTW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NTTW is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with NTTW.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2008-13
*/
#include <iostream>
#include <sstream>

using namespace std;

///C++ is used here since there are issues with printing size_t portably across platforms in C
extern "C"
{
    //NTTW C Library
    #include "nttw/array.h"
}

///Globals

///Function Prototypes

int main(int argc, char *argv[])
{
    cout << ">| Machine Datatype Properties Program -----------------------------" << endl;
    cout << ">| Copyright Shekhar S. Chandra, 2007-9" << endl;
    cout << ">| Machine Integer Size of " << BITS << " bits" << endl;

    //--------------------------------------------
    nttw_integer value;
    nttw_big_integer value2;
    nttw_real dblValue;
    nttw_integer *arrayPtr;

    arrayPtr = array_1D(2);

    cout << "\n>| Library Types ------" << endl;
    cout << ">| Library (nttw_integer) Type Size: " << sizeof(value) << ", Alignment: " << ALIGNOF(value) << " bytes" << endl;
    cout << ">| Library Big (nttw_big_integer) Type Size: " << sizeof(value2) << ", Alignment: " << ALIGNOF(value2) << " bytes" << endl;
    cout << ">| Complex Library (nttw_real)  Type Size: " << sizeof(dblValue) << ", Alignment: " << ALIGNOF(dblValue) << " bytes" << endl;
    cout << ">| Array Pointer of Library Type Size: " << sizeof(arrayPtr) << ", Alignment: " << ALIGNOF(arrayPtr) << " bytes" << endl;
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

    //--------------------------------------------
    free_array(arrayPtr);
    cout << ">| Output Complete" << endl;
    return 0;
}
