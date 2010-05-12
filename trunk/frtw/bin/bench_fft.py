#!/usr/bin/python

import os, time

limit = 13 #power of 2
runs = 100 #number of runs per length
for i in range(6,limit+1):
    for j in range(0,runs):
        print ">>>> Running FFT run " + str(j+1) + " for dyad " + str(1 << i)
        os.system( "./fft_bench " + str(1 << i) )
        time.sleep(2)
print ">>>> Done Benching FFT"