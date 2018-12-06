# Getting Started
This library is dependent on FFTW (http://www.fftw.org) and for Windows using Visual C++, you will need to setup lib files properly working [their guide](http://www.fftw.org/install/windows.html).
This project uses CMake.

# Windows
* I downloaded the fftw-3.3.5-dll64.zip file, extracted it to a directory (used later) and setup .lib files as explained by the [FFTW Windows guide](http://www.fftw.org/install/windows.html). You can see the log of what I did for Visual Studio 14 (VS 2015) [on the wiki](https://github.com/shakes76/finite-transform-library/wiki/FFTW-Setup-Guide-for-Windows).
* Using cmake-gui, run configure making sure to point to the FFTW libraries directory above by setting FFTW_LIBRARIES and FFTW_INCLUDES
* Be sure to point to the library you want, usually libfftw3-3.lib
* Generate from CMake and open solution file .sln
* Build and run!

# Linux
* Install FFTW3 from your distro's package manager
* Change to FTL's source directory and create a build directory
* Run ccmake
* Build using make (or whatever generator you decided on)

It will look something like:
* 'mkdir build' in the ftl directory
* 'cd build'
* 'ccmake ..' OR 'cmake ..'
* 'make'
* The binaries etc. will be in the relevant (frtw/nttw) directories

# ANDROID
This was done a while ago and might be out of date.

Add the following to your .bashrc
#Android Dev
export ANDROID_NDK=~/Dev/Install/android-toolchain
export ANDROID_SDK=~/Dev/android-sdk-linux
export ANDROID_NDK_TOOLCHAIN_ROOT=$ANDROID_NDK
export PATH=$ANDROID_SDK/tools:$ANDROID_SDK/platform-tools:$PATH
export PATH=$ANDROID_NDK/bin:$ANDROID_NDK/arm-linux-androideabi/bin:$PATH
export ANDROID_CMAKE=~/Dev/Install/android-cmake
export ANDTOOLCHAIN=$ANDROID_CMAKE/toolchain/android.toolchain.cmake
alias android-cmake='cmake -DCMAKE_TOOLCHAIN_FILE=$ANDTOOLCHAIN '

The above variables assume you have
1. Setup NDK standalone toolchain as
  export NDK=~/Dev/android-ndk-r9
  $NDK/build/tools/make-standalone-toolchain.sh --platform=android-4 --install-dir=$HOME/Dev/Install/android-toolchain
2. Downladed android-cmake from the Google Code hg repository into ~/Dev/Install/android-cmake
3. Downloaded the Android SDK to ~/Dev/android-sdk-linux and run the following to update and install
  android sdk
4. On a 64-bit Linux system, install the 32-bit libraries etc.
  sudo apt-get install ia32-libs
