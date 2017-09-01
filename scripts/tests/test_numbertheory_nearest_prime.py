# -*- coding: utf-8 -*-
"""
Test if nearest prime function works in NT
Created on Tue Oct 28 10:48:41 2014

@author: uqscha22
"""
import _libpath #add custom libs
import finitetransform.numbertheory as nt

n = 98

print "Nearest prime to", n, "is", nt.nearestPrime(n)
