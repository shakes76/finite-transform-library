![http://finitetransform.sourceforge.net/images/ftl_logo2.png](http://finitetransform.sourceforge.net/images/ftl_logo2.png)

Library for Finite Transforms such as the Number Theoretic Transform (NTT) and others. Current modules include NTTW which is the NTT library with high resolution (microsecond) timing, basic array and imaging. The transforms are optimised for performance. Finite Transform Library (FTL) has two modules that cover two transforms defined on the finite circle not commonly found in open source.

# NTTW #
The first is the NTTW library for Number Theoretic Transforms (NTTs). NTTs are pure digital transforms computed on integer fields. This allows convolutions to be done exactly and essentially without overflow. This module is inspired by the [APFloat library](http://www.apfloat.org/apfloat/) by Mikko Tommila et al where it has been utilised for an arbitrary precision library. See the NTT description by APFloat [here](http://www.apfloat.org/ntt.html). The NTTW library's goal is to focus only on NTTs and to achieve comparable performance to FFTs.

The NTTW module also features a simple Image Module and a High Resolution (micro-second) Timing module all of which are cross platform.

Version 1.31

---

  * Changed the license to LGPL
  * Added more tests and apps

# FRTW #
The second module is designed to compute the Discrete Radon Transform and their variants. A primary release has been made available that also covers the Mojette Transform. See my recent paper entitled "Fast Mojette Transform for Discrete Tomography" arXiv:1006.1965v1 [physics.med-ph].

Version 1.50

---

  * Added the Number Theoretic Radon Transform (NRT) to the Radon module (published in the DICTA 2009 paper and IET CV paper under review)
  * Added ghosts module and deghosting capabilty to library (published in the Ghosts paper of IEEE Trans. Image Proc. (TIP) October 2012)
  * Added applications building via CMake
  * Updated the library to the LPGL license
  * Updated NTTW detection in the CMake configs

# Dev Blog News #
<wiki:gadget url="http://google-code-feed-gadget.googlecode.com/svn/trunk/gadget.xml" up\_feeds="http://l3mmings.blogspot.com/feeds/posts/default" width="780" height="340" border="0" up\_showaddbutton="0"/>