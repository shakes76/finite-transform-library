#!/usr/bin/python

import os

step = 2
limit = 10 #power of 2
for i in range(7,limit+1):
	print ">>>> Running NTT for dyadic " + str(pow(step,i))
	os.system( "./fntt_bench " + str(pow(step,i)) )
print ">>>> Done Benching NTT"