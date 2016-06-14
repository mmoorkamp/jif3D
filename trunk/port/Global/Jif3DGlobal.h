//============================================================================
// Name        : Jif3DGlobal.h
// Author      : Ashley Davis (ad299) https://github.com/SgtCoDFish
// Version     :
// Copyright   : 2016, AD
//============================================================================

/*
 * This file is mainly intended for Windows DLL support when exporting symbols
 * from a source file.
 * 
 * It should only have any effect on Windows.
 */

#ifndef GLOBAL_JIF3DGLOBAL_H_
#define GLOBAL_JIF3DGLOBAL_H_

#if defined(_WIN32) || defined(_WIN64) || defined(JIF3D_WINDOWS_BUILD)

#ifdef JIF3D_BUILDING
	#define J3DEXPORT __declspec(dllexport)
#else
	#define J3DEXPORT __declspec(dllimport)
#endif

#else
	#define J3DEXPORT 
#endif

#endif
