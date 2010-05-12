## C++ Compiler flags for performance enhancements for DGL, DGV and other branches

linux-g++ {
    #gcc-4.2.4
    FLAGS_CXX += -O3
    FLAGS_CXX += -march=native
    #FLAGS_CXX += -march=core2
    FLAGS_CXX += -mtune=native #Implied by the above flag
    #FLAGS_CXX += -mtune=core2 #Implied by the above flag
    #FLAGS_CXX += -mmmx #older instruction set
    FLAGS_CXX += -mfpmath=sse
    FLAGS_CXX += -msse2
    #FLAGS_CXX += -malign-double #Mostly for Pentiums and not for 64-bit
    #FLAGS_CXX += -m64 #64-Bit Only on x64 CPUs
    FLAGS_CXX += -funroll-all-loops #Un roll loops
    FLAGS_CXX += -momit-leaf-frame-pointer #Discard Debug help
    FLAGS_CXX += -fomit-frame-pointer #Discard Debug help
    FLAGS_CXX += -finline-functions
    FLAGS_CXX += -Wno-deprecated #VTK uses deprecated c libraries
    FLAGS_CXX += -fwrapv #Defines overflow of signed integers as wrapping
    FLAGS_CXX += -fopenmp #Multi Processor support
    #FLAGS_CXX += -funsafe-loop-optimizations #Ignore loop overflow checking
    #FLAGS_CXX += -foptimize-register-move
    #FLAGS_CXX += -fconserve-stack -floop-interchange 
    #FLAGS_CXX += -msahf #remainder 64-bit intructions
    #FLAGS_CXX += -mmovbe #byte reversal flag
    #FLAGS_CXX += -pg #Profiling

    QMAKE_CXXFLAGS_RELEASE -= -O2 #Ensure -O3 is used for Blitz++
    QMAKE_CXXFLAGS_RELEASE += $$FLAGS_CXX
}

win32-msvc2005 {
    #msvc 2008
    FLAGS_CXX += /Ox
    FLAGS_CXX += /Ot /Oi 
    FLAGS_CXX += /QIfist #Suppress float to int function call
    FLAGS_CXX += /fp:fast # Fast Floating points
    FLAGS_CXX += /GL
    FLAGS_CXX += /arch:SSE2
    FLAGS_CXX += /MP #Multi-processor switch for building
    FLAGS_CXX += /GS- #Turn off Security Buffer
    FLAGS_CXX += /favor:blend
    #FLAGS_CXX += /favor:INTEL64 #Intel x64
    #FLAGS_CXX += /favor:AMD64 #AMD x64
    FLAGS_LINKER += /LTCG /OPT:icf,ref

    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += $$FLAGS_CXX
    QMAKE_LFLAGS_WINDOWS += $$FLAGS_LINKER
    QMAKE_LFLAGS_WINDOWS_DLL += $$FLAGS_LINKER
}
