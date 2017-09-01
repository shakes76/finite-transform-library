# -*- coding: utf-8 -*-
"""
Test Read PGM member
Created on Sun Jan 25 21:52:38 2015

@author: shakes
"""
import _libpath #add custom libs
import finitetransform.imageio as imageio

saveImage = False
image, depth = imageio.readPGM("../data/lena512.pgm")
print "Image dtype:", image.dtype

if saveImage:
    imageio.imsave("../data/lena.png", image)

#Plot
import matplotlib.pyplot as plt

plt.gray()
plt.imshow(image)
plt.axis('off')
plt.title('Image')

plt.show()
