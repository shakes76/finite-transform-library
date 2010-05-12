/**
* Test Timing Module
* \author Shekhar S. Chandra, 2009
*/
#include <stdio.h>
#include "timing.h"

int main()
{
	int temp, i;

	START_TIMER;
	for (i = 0; i< 242000000; i++)
		temp+=temp;
    STOP_TIMER;

#ifndef WIN32
	printf("%li:%li\n",SECONDS_ELAPSED,NANOSECONDS_ELAPSED);
#endif
	printf("%llu usecs.\n",MICROSECONDS_ELAPSED);

	return 0;
}
