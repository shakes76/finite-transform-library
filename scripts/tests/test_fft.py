'''
Test some properties of the DFT
'''
import _libpath #add custom libs
import numpy as np
import scipy.fftpack as fftpack
import finitetransform.imageio as imageio #local module

#create suitably sized Lena image
#zero pad to larger prime size to allow non-periodic rotations.
n = 256
image = imageio.lena(n)
imageFloat = image.astype('float32', copy=False)
print "Size:", image.shape
print "Data Type:", image.dtype

#FFT
fftLena = fftpack.fft2(image) #the '2' is important
fftLenaShifted = fftpack.fftshift(fftLena)

#power spectrum
powSpectLena = np.abs(fftLenaShifted)

#iFFT
ifftLena = fftpack.ifft2(fftLena)

#Plot
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(8, 5))

plt.gray()

ax[0].imshow(imageFloat)
ax[0].axis('off')
ax[0].set_title('Lena')
ax[1].imshow( np.log10(powSpectLena) ) # log plot
ax[1].axis('off')
ax[1].set_title('FFT (log10) of Lena')
ax[2].imshow( np.real(ifftLena) )
ax[2].axis('off')
ax[2].set_title('iFFT of Lena')

plt.show()