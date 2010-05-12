NTTW Library only needs C standard library, build projects are provided from qmake (Qt 4)
which is a cross-platform build tool. Projects are also provided for Code::Blocks a
cross platform IDE.

-----------------------------------------------
LINUX - gcc-4.2.4 with Qt 4 installed
Build Instructions:

1. Build the Library (nttw.pro)
	qmake nttw.pro
	make
2. Build the Benches
	cd bench
	qmake
	make
3. Build the Tests
	cd ../tests
	qmake
	make
4. Build the Applications
	cd ../apps
	qmake
	make
	
Note: To use older versions of gcc-4 comment out the compiler flags (-march and -mtune flags) in file cxxflags.pri.
-----------------------------------------------
WINDOWS - msvc-2005 or msvc-2008
Build Instructions:    

1. Build the Library (nttw.pro)
	qmake nttw.pro
	nmake
2. Build the Benches
	cd bench
	qmake
	nmake
3. Build the Tests
	cd ../tests
	qmake
	nmake
4. Build the Applications
	cd ../apps
	qmake
	nmake
    
-----------------------------------------------
The Library qmake configuration is stored in file: config.pri
The CXX compiler tweaks are stored in file: cxxflags.pri
The Doxygen configuration is stored in file: Doxyfile.doxy

