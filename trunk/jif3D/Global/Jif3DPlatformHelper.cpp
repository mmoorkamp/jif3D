//============================================================================
// Name        : Jif3DPlatformHelper.cpp
// Author      : Ashley Davis (ad299) https://github.com/SgtCoDFish
// Version     :
// Copyright   : 2016, AD
//============================================================================

#include <cstdlib>
#include <ctime>

#include "Jif3DPlatformHelper.h"

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#else
/*
 * We assume that if we're not on Windows, we're likely to be on Linux.
 * Anything that works for Linux is more likely to work for other
 * OSes which are unsupported at this time (e.g. Mac OSX).
 */
#include <sys/types.h>
#include <unistd.h>

#endif


namespace jif3D
{
namespace platform
{

int get_process_id()
{
#ifdef _WIN32
	return (int) GetCurrentProcessId();
#else
	return (int) getpid();
#endif
}

double drand48()
{
#ifdef _WIN32
	return ((double)rand()/RAND_MAX);
#else
	return ::drand48();
#endif
}

void srand48(long seedval)
{
#ifdef _WIN32
	std::srand(seedval);
#else
	::srand48(seedval);
#endif
}

}
}
