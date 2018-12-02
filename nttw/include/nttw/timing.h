/**
 * NTTW High Resolution Timing Module
 * \file timing.h
 * \brief High Resolution Timing Header/Object for the NTTW C Library.
 *
 * This file defines the functions that wrap all the OS dependent high resolution timing functions.
 * Only supports Windows and POSIX (including Linux) Operating Systems.
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
 * along with NTTW. If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Shekhar S. Chandra, 2009
*/
#ifndef TIMING_H_INCLUDED
#define TIMING_H_INCLUDED

/**
 * \defgroup Timing High Resolution (micro-second) Timing Module
 * \brief Timing functions for high resolution (micro-second) timing that can be used within any C code.
 * The performence is OS dependent and only supports Windows and POSIX/Linux.
 *
 * The main functions are START_TIMER, STOP_TIMER and MICROSECONDS_ELAPSED.
 * These commands require semi-colons at the end if used on its own.
*/
//@{
#if defined WIN32
    //Do Windows Equivalent
    #include <windows.h>
    #include "global.h"

	//Global Time variables
    static LARGE_INTEGER startTime; //!< Start Time Global Variable (Windows)
    static LARGE_INTEGER stopTime; //!< Stop Time Global Variable (Windows)

	/**
	 * \fn diff_times(const LARGE_INTEGER start, const LARGE_INTEGER end)
	 * Finds the difference in times. For Windows only.
	*/
    NTTW_DLL_SYM unsigned long long diff_times(const LARGE_INTEGER start, const LARGE_INTEGER end);

	//Timer Macros (to avoid function call overheads)
    #ifndef START_TIMER //!< Start Timer Macro (Windows)
        #define START_TIMER ( QueryPerformanceCounter(&startTime) )
    #endif

    #ifndef RESTART_TIMER //!< Restart Timer Macro (Windows)
        #define RESTART_TIMER ( QueryPerformanceCounter(&startTime) )
    #endif

    #ifndef STOP_TIMER //!< Stop Timer Macro (Windows)
        #define STOP_TIMER ( QueryPerformanceCounter(&stopTime) )
    #endif

    #ifndef TIME_ELAPSED //!< Time Elapsed Macro - EMPTY (Windows)
        #define TIME_ELAPSED
    #endif

    #ifndef SECONDS_ELAPSED //!< Time Elapsed Macro - EMPTY (Windows)
        #define SECONDS_ELAPSED
    #endif

    #ifndef NANOSECONDS_ELAPSED //!< Time Elapsed Macro - EMPTY (Windows)
        #define NANOSECONDS_ELAPSED
    #endif

    #ifndef MILLISECONDS_ELAPSED //!< Returns time in milli-seconds (Windows)
        #define MILLISECONDS_ELAPSED ( diff_times(startTime,stopTime)/1000 )
    #endif

    #ifndef MICROSECONDS_ELAPSED //!< Returns time in micro-seconds (Windows)
        #define MICROSECONDS_ELAPSED ( diff_times(startTime,stopTime) )
    #endif

#else
    //Do POSIX Timers
    #include <time.h>
    #include "global.h"

	//Global Time variables
    static struct timespec startTime; //!< Start Time Global Variable (Linux)
    static struct timespec stopTime; //!< Start Time Global Variable (Linux)

	/**
	 * \fn diff_times(const struct timespec start, const struct timespec end)
	 * Finds the difference in times. For POSIX/Linux only.
	*/
    NTTW_DLL_SYM struct timespec diff_times(const struct timespec start, const struct timespec end);

	//Timer Macros (to avoid function call overheads)
    #ifndef START_TIMER //!< Start Timer Macro (Linux)
        #define START_TIMER ( clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startTime) )
    #endif

    #ifndef RESTART_TIMER //!< Restart Timer Macro (Linux)
        #define RESTART_TIMER ( clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startTime) )
    #endif

    #ifndef STOP_TIMER //!< Stop Timer Macro (Linux)
        #define STOP_TIMER ( clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stopTime) )
    #endif

    #ifndef TIME_ELAPSED //!< Time elapsed so far, return timespec struct (Linux)
        #define TIME_ELAPSED ( diff_times(startTime, stopTime) )
    #endif

    #ifndef SECONDS_ELAPSED //!< Returns seconds elapsed (Linux)
        #define SECONDS_ELAPSED ( diff_times(startTime, stopTime).tv_sec )
    #endif

    #ifndef NANOSECONDS_ELAPSED //!< Return nano-seconds elapsed, use with SECONDS_ELAPSED to get total time (Linux)
        #define NANOSECONDS_ELAPSED ( diff_times(startTime, stopTime).tv_nsec )
    #endif

    #ifndef MILLISECONDS_ELAPSED //!< Returns time in milli-seconds (Linux)
        #define MILLISECONDS_ELAPSED ( diff_times(startTime, stopTime).tv_sec*1000ULL + diff_times(startTime, stopTime).tv_nsec/1000000ULL )
    #endif

    #ifndef MICROSECONDS_ELAPSED //!< Returns time in micro-seconds (Linux)
        #define MICROSECONDS_ELAPSED ( diff_times(startTime, stopTime).tv_sec*1000000ULL + diff_times(startTime, stopTime).tv_nsec/1000ULL )
    #endif

#endif //WIN32
//@}
#endif // TIMING_H_INCLUDED
