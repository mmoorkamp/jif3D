# Microsoft Developer Studio Project File - Name="meschach" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=meschach - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "meschach.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "meschach.mak" CFG="meschach - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "meschach - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "meschach - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "meschach - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x807 /d "NDEBUG"
# ADD RSC /l 0x807 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "meschach - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ  /c
# ADD BASE RSC /l 0x807 /d "_DEBUG"
# ADD RSC /l 0x807 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "meschach - Win32 Release"
# Name "meschach - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\bdfactor.c
# End Source File
# Begin Source File

SOURCE=..\bkpfacto.c
# End Source File
# Begin Source File

SOURCE=..\chfactor.c
# End Source File
# Begin Source File

SOURCE=..\copy.c
# End Source File
# Begin Source File

SOURCE=..\err.c
# End Source File
# Begin Source File

SOURCE=..\extras.c
# End Source File
# Begin Source File

SOURCE=..\fft.c
# End Source File
# Begin Source File

SOURCE=..\givens.c
# End Source File
# Begin Source File

SOURCE=..\hessen.c
# End Source File
# Begin Source File

SOURCE=..\hsehldr.c
# End Source File
# Begin Source File

SOURCE=..\init.c
# End Source File
# Begin Source File

SOURCE=..\iter0.c
# End Source File
# Begin Source File

SOURCE=..\iternsym.c
# End Source File
# Begin Source File

SOURCE=..\itersym.c
# End Source File
# Begin Source File

SOURCE=..\ivecop.c
# End Source File
# Begin Source File

SOURCE=..\lufactor.c
# End Source File
# Begin Source File

SOURCE=..\machine.c
# End Source File
# Begin Source File

SOURCE=..\matlab.c
# End Source File
# Begin Source File

SOURCE=..\matop.c
# End Source File
# Begin Source File

SOURCE=..\matrixio.c
# End Source File
# Begin Source File

SOURCE=..\meminfo.c
# End Source File
# Begin Source File

SOURCE=..\memory.c
# End Source File
# Begin Source File

SOURCE=..\memstat.c
# End Source File
# Begin Source File

SOURCE=..\mfunc.c
# End Source File
# Begin Source File

SOURCE=..\norm.c
# End Source File
# Begin Source File

SOURCE=..\otherio.c
# End Source File
# Begin Source File

SOURCE=..\pxop.c
# End Source File
# Begin Source File

SOURCE=..\qrfactor.c
# End Source File
# Begin Source File

SOURCE=..\schur.c
# End Source File
# Begin Source File

SOURCE=..\solve.c
# End Source File
# Begin Source File

SOURCE=..\sparse.c
# End Source File
# Begin Source File

SOURCE=..\sparseio.c
# End Source File
# Begin Source File

SOURCE=..\spbkp.c
# End Source File
# Begin Source File

SOURCE=..\spchfctr.c
# End Source File
# Begin Source File

SOURCE=..\splufctr.c
# End Source File
# Begin Source File

SOURCE=..\sprow.c
# End Source File
# Begin Source File

SOURCE=..\spswap.c
# End Source File
# Begin Source File

SOURCE=..\submat.c
# End Source File
# Begin Source File

SOURCE=..\svd.c
# End Source File
# Begin Source File

SOURCE=..\symmeig.c
# End Source File
# Begin Source File

SOURCE=..\update.c
# End Source File
# Begin Source File

SOURCE=..\vecop.c
# End Source File
# Begin Source File

SOURCE=..\version.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\dll_export.h
# End Source File
# Begin Source File

SOURCE=..\err.h
# End Source File
# Begin Source File

SOURCE=..\iter.h
# End Source File
# Begin Source File

SOURCE=..\machine.h
# End Source File
# Begin Source File

SOURCE=..\matlab.h
# End Source File
# Begin Source File

SOURCE=..\matrix.h
# End Source File
# Begin Source File

SOURCE=..\matrix2.h
# End Source File
# Begin Source File

SOURCE=..\meminfo.h
# End Source File
# Begin Source File

SOURCE=..\oldnames.h
# End Source File
# Begin Source File

SOURCE=..\sparse.h
# End Source File
# Begin Source File

SOURCE=..\sparse2.h
# End Source File
# Begin Source File

SOURCE=..\zmatrix.h
# End Source File
# Begin Source File

SOURCE=..\zmatrix2.h
# End Source File
# End Group
# End Target
# End Project
