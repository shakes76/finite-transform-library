## C Compiler flags for performance enhancements for DGL, DGV and other branches

linux-g++ {
    #gcc-4.2.4
    FLAGS_C += -O3
    FLAGS_C += -march=native
    #FLAGS_C += -march=core2
    FLAGS_C += -mtune=native #Implied by the above flag
    #FLAGS_C += -mtune=core2 #Implied by the above flag
    #FLAGS_C += -mmmx #older instruction set
    FLAGS_C += -mfpmath=sse
    FLAGS_C += -msse2
    #FLAGS_C += -malign-double #Mostly for Pentiums and not for 64-bit
    #FLAGS_C += -m64 #64-Bit Only on x64 CPUs
    FLAGS_C += -ffast-math #Math Functions do not check arguments, care is needed here.
    FLAGS_C += -funroll-all-loops #Un roll loops
    FLAGS_C += -momit-leaf-frame-pointer #Discard Debug help
    FLAGS_C += -fomit-frame-pointer
    FLAGS_C += -finline-functions
    FLAGS_C += -fwrapv #Defines overflow of signed integers as wrapping
    FLAGS_C += -fopenmp #Multi Processor support
    #FLAGS_C += -funsafe-loop-optimizations #Ignore loop overflow checking
    #FLAGS_C += -foptimize-register-move
    #FLAGS_C += -fconserve-stack -floop-interchange 
    #FLAGS_C += -msahf #remainder 64-bit intructions
    #FLAGS_C += -mmovbe #byte reversal flag
    #FLAGS_C += -pg #Profiling

    QMAKE_CFLAGS_RELEASE -= -O2 #Ensure -O3 is used
    QMAKE_CFLAGS_RELEASE += $$FLAGS_C
}

win32-msvc2005 {
    #msvc 2008
    FLAGS_C += /Ox 
    FLAGS_C += /Ot /Oi
    FLAGS_C += /fp:fast # Fast Floating points
    FLAGS_C += /QIfist #Suppress float to int function call
    FLAGS_C += /GL
    FLAGS_C += /arch:SSE2
    FLAGS_C += /MP #Multi-processor switch for building
    FLAGS_C += /GS- #Turn off Security Buffer
    #FLAGS_C += /favor:blend
    #FLAGS_C += /favor:INTEL64 #Intel x64
    FLAGS_C += /favor:AMD64 #AMD x64
    FLAGS_LINKER += /LTCG /OPT:icf,ref

    QMAKE_CFLAGS_RELEASE -= -O2
    QMAKE_CFLAGS_RELEASE += $$FLAGS_C
    QMAKE_LFLAGS_WINDOWS += $$FLAGS_LINKER
    QMAKE_LFLAGS_WINDOWS_DLL += $$FLAGS_LINKER
}
