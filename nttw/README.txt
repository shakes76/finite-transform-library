NTTW Library only needs C standard library, build projects are provided for CMake 2.6+
which is a cross-platform build tool. One can generate makefiles and Codeblocks, Visual Studio and Eclipse projects
from CMake.

BUILD
-----------------------------------------------
Simply create a build directory, then run ccmake in that directory or cmake-gui from anywhere. For ccmake, it looks like

 $ cd /path/to/nttw-source/
 $ mkdir build
 $ cd build
 $ ccmake ..
 
After configuring as desired, on Linux with GCC:
 $ make -j 4

Or on Windows, via the command line:
 $ jom /j 4 or nmake
If you dont have the open-source multi-threaded jom, use nmake.
     
For more help, consult any CMake guid or manual.

-----------------------------------------------
The Doxygen configuration is stored in file: nttw.doxy

MODULES
-----------------------------------------------
Imaging - image.h
The image library supports IO of 32-bit signed PGM images.
It also supports some basic arithmetic and operations such as flip and crop.

Array - array.h
Supports dynamically allocated arrays of different types.
Initialisation methods for the arrays are also provided.

Timing - timing.h
Cross-platform micro-second timing supprt

Prime Number Operations - prime.h
Finding closest prime from list, finding factors of a number etc.
Finding primitive roots, multiplicative inverses and GCDs etc.

Number Transforms - number.h
Number Theoretic Transforms for integer only signals