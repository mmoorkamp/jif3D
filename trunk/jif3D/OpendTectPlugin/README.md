Jif3D OpendTect Plugin
----------------------

Jif3D comes with a plugin environment for OpendTect, an industry-standard geophysics modelling tool. This document gives an outline of how to compile and integrate the plugin.

The plugin has been tested on Linux (Ubuntu with g++ 5.2, cmake 3.2.2) and Windows 10 (with Visual Studio 2013 and 2015, cmake 3.2.2) using OpendTect version 6.0.2.

If you test successfully on other OSes, compilers or versions please do let us know!

## Building
Navigate to the OpendTectPlugin folder and run CMake with the generator of your choice.

You'll also need to define the variable -DOpendTect_DIR to the folder where your copy of OpendTect is installed or CMake is likely to fail.

Contrary to how CMake is usually used, you should run the CMake command in the same folder as the `CMakeLists.txt` file rather than in a separate file.

An example command would be:

```
cmake -DOpendTect_DIR=~/OpendTect/6.0.2 .
```

# Troubleshooting & Platform Specific Tips
## 32/64 Bit
Take care to use a 64 bit generator if you're running a 64-bit version of OpendTect - especially on Windows. This could be a source of hard to identify problems.

For example on Windows be sure to be specific with your generator on 64-bit; use:

```
cmake -DOpendTect_DIR=<path> -G "Visual Studio 14 2015 Win64 ."
```

for example. The `Win64` is important.

## Windows
Version 6.0.2 of OpendTect has, at the time of writing, no support for MSVC 2015 (2013 is supported fine).

The fix for this is very simple but requires editing your OpendTect install. Go to OpendTect_DIR and open CmakeModules/CreateLaunchers.cmake.

Patch line 80 `if(MSVC12)` to:

```
if(MSVC14)
	set(LAUNCHER_LINESEP "\n")
	set(USERFILE_VC_VERSION 14.00)
	set(USERFILE_EXTENSION user)
	set(VCPROJ_TYPE vcxproj)
elseif(MSVC12)
```

that is, replace the line with the lines above. A pull request has been submitted to the main OpendTect repository but may not have been integrated.

## Linux
### Compilers
Be sure to use gcc/g++ on Linux; clang/clang++ seems to produce issues. Use `-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++` if you've got a non-default compiler set up.

### ABI Versioning
As of versions 6.0.2 and below, OpendTect is compiled with the pre-C++11 ABI (using g++ from the 4.* range). If you have a more modern compiler (e.g. Ubuntu 15.10 or above) 
you **will** run into problems with ABI versioning if you try to change it. The CMakeLists.txt file provided already sets the ABI to the old version so it should be safe to 
use out of the box on any system.

If you're using a much newer future version of OpendTect, you may need to change this flag depending on your system. Likely for most current systems you won't need to change it.

# Other Details
## Manually Compiling OpendTect on Windows
If you want to compile OpendTect on Windows with Visual Studio 2015 for whatever reason, follow this guide for qt4: https://stackoverflow.com/questions/32848962/how-to-build-qt-4-8-6-with-visual-studio-2015-without-official-support