# -*- coding: utf-8 -*-
"""
Mask image test

Created on Tue Jul  5 21:55:15 2016

@author: shakes
"""
import _libpath #add custom libs
import finitetransform.imageio as imageio #local module

#parameters
N = 32
M = 64

#create test image
lena = imageio.lena(N,M)
mask = imageio.mask(N,M)
print mask

#mask and crop
print lena[mask].reshape((N,N))
cropLena = imageio.immask(lena, mask, N, N)
#print cropLena.shape

#plot
import matplotlib.pyplot as plt

plt.subplot(131)
plt.imshow(lena, interpolation='nearest', cmap='gray')
plt.title('Image')
plt.subplot(132)
plt.imshow(mask, interpolation='nearest', cmap='gray')
plt.title('Mask')
plt.subplot(133)
plt.imshow(cropLena, interpolation='nearest', cmap='gray')
plt.title('Cropped Image')

plt.show()
print "Complete"