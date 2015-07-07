# Introduction #

FRTW is especially designed to do Discrete Radon Transforms (DRT). These include the Fast/Finite Radon Transform (FRT) of Bolker, Grigoryan, Gertner, Fill and others (1986-1989) and the Mojette Transform of Guedon et al. (1995). See my [pre-print paper "Fast Mojette Transform for Discrete Tomography"](http://arxiv.org/abs/1006.1965) on these transforms for more details.

The library depends on [FFTW](http://www.fftw.org/) and NTTW, the latter of which is also part of the FTL libraries.

# Details #

The library is structured as follows. The modules include are:
  * Radon module for DRTs (radon.h/.c).
  * Noise module using Random Number Generators for producing Gaussian/Poisson noise (noise.h/.c).
  * Mojette module for the (Fast) Mojette Transform operations (mojette.h/.c).
  * Vector (nD) module for simple vector operations (vector.h/.c).
  * Complex Array wrappers of the FFTW complex arrays (array\_complex.h/.c).
  * Fourier module for using FFTW (fourier.h/.c).
See the [FRTW API](http://finitetransform.sourceforge.net/frtw_api/index.html) documentation for more details.

# Building #
## QMake ##
Have qmake (qt4-qmake) installed and running on your system. This is a cross-platform build tool that comes with [Qt4](http://qt.nokia.com/). This is especially easy on Ubuntu, as you just install the package qt4-qmake. Then its just a matter of running
```
qmake
make
make install
```

for both NTTW and then FRTW. On Windows use nmake instead of make. You will have to repeat this commands to build the apps, tests or bench applications. If you do not have or cannot get qmake, then a sample Makefile is provided in the downloads section, which you must edit.

## CMake ##
Have cmake installed and running on your system. This is a cross-platform build tool that is an independent tool which you can grab from [CMake.org](http://www.cmake.org/) for every system. Cmake is advanteous because its always a small install. In the source directory
```
cmake .
make
make install
```
for both NTTW and then FRTW. On Windows use nmake instead of make. It is best to do an out of source build, so you can use the source tree multiple times for multiple architectures. Use the CMake GUI for this. On Linux you can use ccmake. See the cmake documentation [here](http://www.cmake.org/cmake/help/documentation.html).

# Apps #
## Fast Mojette Transform ##
There are two ways to compute the Fast Mojette Transform (FMT). Firstly, one may do the following, where the commands or options are in round brackets:
  1. Generate the Farey angle set (fmt\_angle) corresponding to either the simple (2) angle set or the L1-norm angle set (1).
  1. Mojette Transform (mt) the image using the generated angle set.
  1. Then convert the projections to the Fast/Finite Radon Transform (FRT) (mt2frt). This generates the FRT space which my be inverted using the iFRT app. as a test, the app automatically generates this for you at the time of writing this guide.

Secondly, one may to the following:
  1. Simply run the FMT app (fmt), where the angle set is generated, MT is taken and converted to FRT projections on the fly and a FRT space produced.
  1. Run the inverse FMT (ifmt) app to invert the result, which just happens to be similar to the inverse FRT app (ifrt).

## Fast/Finite Radon Transform ##
The apps included are drt, frt and ifrt for dyadic (power of two sizes). The DRT is a cubic complexity implementation of the FRT.

## Others ##
Apps mse, crop and line produce the Mean Squared Error, crop an image to specific size and produce a Diophantine line file for lines within the FRT. The latter is a file that can be viewed in [DGV](http://code.google.com/p/discrete-geometry-viewer/).

## Troubleshooting ##
### FRTW ###
  * If you get the error "nttw/global.h: No such file or directory", you havn't installed NTTW. NTTW is a prerequisite for FRTW.