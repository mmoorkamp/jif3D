//============================================================================
// Name        : Jif3DTesting.h
// Author      : Ashley Davis (ad299) https://github.com/SgtCoDFish
// Version     :
// Copyright   : 2016, AD
//============================================================================

/*
 * This file handles inclusion of testing libraries in a cross platform way.
 */

#ifndef GLOBAL_JIF3DTESTING_H_
#define GLOBAL_JIF3DTESTING_H_
#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#endif

#if BOOST_ALL_DYN_LINK
	#include <boost/test/unit_test.hpp>
#else
	#include <boost/test/included/unit_test.hpp>
#endif

#endif
