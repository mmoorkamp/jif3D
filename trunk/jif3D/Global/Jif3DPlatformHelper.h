//============================================================================
// Name        : Jif3DPlatformHelper.h
// Author      : Ashley Davis (ad299) https://github.com/SgtCoDFish
// Version     :
// Copyright   : 2016, AD
//============================================================================

/*
 * This file contains abstractions to help make Jif3D more portable across
 * different platforms and operating systems. Where a function might have
 * a difference across platforms (e.g. between GNU/Linux and Windows) it should
 * be defined here and implemented in Jif3DPlatformHelper.cpp, included in
 * "#define guards". For example, a Windows only implementation would check
 * "#ifdef _WIN32"
 *
 *
 */
 
#ifndef GLOBAL_JIF3DPLATFORMHELPER_H_
#define GLOBAL_JIF3DPLATFORMHELPER_H_

#include "Jif3DGlobal.h"

#if defined _WIN32
	#define JIF3D_PLATFORM "WINDOWS"
#elif defined __linux
	#define JIF3D_PLATFORM "LINUX"
#else
	#define JIF3D_PLATFORM "UNKNOWN"
	#warning Unknown platform detected in Global/Jif3DPlatformHelper.h - compilation unlikely to succeed.
#endif

namespace jif3D
{
namespace platform
{

J3DEXPORT int get_process_id();

J3DEXPORT double drand48();
J3DEXPORT void srand48(long seedval);

}
}

#endif
