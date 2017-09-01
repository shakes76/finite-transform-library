# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 14:30:21 2015

@author: shakes
"""
import _libpath #add custom libs
import numpy as np
import finitetransform.farey as farey
                
n = 32

#Farey vectors
fareyVectors = farey.Farey()
#~ fareyVectors.compactOn()
fareyVectors.generate(n-1, 8)
#fareyVectors.generate2(n, 8)
vectors = fareyVectors.vectors
print vectors

countX = 0
countY = 0
vectorImage = np.zeros((2*n+1,2*n+1))
vectorImage[n, n] += 1
for vector in vectors:
    x = abs(vector.real)
    y = abs(vector.imag)
    vectorImage[int(vector.real+n), int(vector.imag+n)] += 1
    countX += x
    countY += y

print "Count x:", countX, "Count y:", countY

#plot
import matplotlib.pyplot as plt

plt.gray()

iax = plt.imshow(vectorImage, interpolation='nearest')
icbar = plt.colorbar(iax, cmap='gray')
plt.title("Farey Image for n="+str(n))

plt.show()