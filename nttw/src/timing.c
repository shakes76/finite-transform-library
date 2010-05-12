/**
 * NTTW High Resolution Timing Module
 * \file timing.c
 * \brief High Resolution Timing Source for the NTTW C Library.
 *
 * This file implements the functions that wrap all the OS dependent high resolution timing functions.
 * Only supports Windows and POSIX (including Linux) Operating Systems.
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
 * \author Shekhar S. Chandra, 2009
*/
#include "timing.h"

#if defined WIN32

	//Windows Differencing
    unsigned long long diff_times(const LARGE_INTEGER start, const LARGE_INTEGER end)
    {
        LARGE_INTEGER frequency;

        QueryPerformanceFrequency(&frequency);

        return (end.QuadPart - start.QuadPart)*1000000ULL / frequency.QuadPart;
    }

#else

	//Linux Differencing
    struct timespec diff_times(const struct timespec start, const struct timespec end)
    {
        struct timespec temp;
        if ((end.tv_nsec-start.tv_nsec)<0)
        {
            temp.tv_sec = end.tv_sec-start.tv_sec-1;
            temp.tv_nsec = 1000000000ULL+end.tv_nsec-start.tv_nsec;
        }
        else
        {
            temp.tv_sec = end.tv_sec-start.tv_sec;
            temp.tv_nsec = end.tv_nsec-start.tv_nsec;
        }
        return temp;
    }

#endif
