
OUTDIR=.\
INTDIR=.\

ALL : "$(OUTDIR)\meschachd.dll" "$(OUTDIR)\meschachd.bsc" ".\Debug\regsvr32.trg"


CLEAN :
	-@erase "$(INTDIR)\bdfactor.obj"
	-@erase "$(INTDIR)\bdfactor.sbr"
	-@erase "$(INTDIR)\bkpfacto.obj"
	-@erase "$(INTDIR)\bkpfacto.sbr"
	-@erase "$(INTDIR)\chfactor.obj"
	-@erase "$(INTDIR)\chfactor.sbr"
	-@erase "$(INTDIR)\copy.obj"
	-@erase "$(INTDIR)\copy.sbr"
	-@erase "$(INTDIR)\err.obj"
	-@erase "$(INTDIR)\err.sbr"
	-@erase "$(INTDIR)\extras.obj"
	-@erase "$(INTDIR)\extras.sbr"
	-@erase "$(INTDIR)\fft.obj"
	-@erase "$(INTDIR)\fft.sbr"
	-@erase "$(INTDIR)\givens.obj"
	-@erase "$(INTDIR)\givens.sbr"
	-@erase "$(INTDIR)\hessen.obj"
	-@erase "$(INTDIR)\hessen.sbr"
	-@erase "$(INTDIR)\hsehldr.obj"
	-@erase "$(INTDIR)\hsehldr.sbr"
	-@erase "$(INTDIR)\init.obj"
	-@erase "$(INTDIR)\init.sbr"
	-@erase "$(INTDIR)\iter0.obj"
	-@erase "$(INTDIR)\iter0.sbr"
	-@erase "$(INTDIR)\iternsym.obj"
	-@erase "$(INTDIR)\iternsym.sbr"
	-@erase "$(INTDIR)\itersym.obj"
	-@erase "$(INTDIR)\itersym.sbr"
	-@erase "$(INTDIR)\ivecop.obj"
	-@erase "$(INTDIR)\ivecop.sbr"
	-@erase "$(INTDIR)\lufactor.obj"
	-@erase "$(INTDIR)\lufactor.sbr"
	-@erase "$(INTDIR)\machine.obj"
	-@erase "$(INTDIR)\machine.sbr"
	-@erase "$(INTDIR)\matlab.obj"
	-@erase "$(INTDIR)\matlab.sbr"
	-@erase "$(INTDIR)\matop.obj"
	-@erase "$(INTDIR)\matop.sbr"
	-@erase "$(INTDIR)\matrixio.obj"
	-@erase "$(INTDIR)\matrixio.sbr"
	-@erase "$(INTDIR)\meminfo.obj"
	-@erase "$(INTDIR)\meminfo.sbr"
	-@erase "$(INTDIR)\memory.obj"
	-@erase "$(INTDIR)\memory.sbr"
	-@erase "$(INTDIR)\memstat.obj"
	-@erase "$(INTDIR)\memstat.sbr"
	-@erase "$(INTDIR)\meschachd.obj"
	-@erase "$(INTDIR)\meschachd.res"
	-@erase "$(INTDIR)\meschachd.sbr"
	-@erase "$(INTDIR)\mfunc.obj"
	-@erase "$(INTDIR)\mfunc.sbr"
	-@erase "$(INTDIR)\norm.obj"
	-@erase "$(INTDIR)\norm.sbr"
	-@erase "$(INTDIR)\otherio.obj"
	-@erase "$(INTDIR)\otherio.sbr"
	-@erase "$(INTDIR)\pxop.obj"
	-@erase "$(INTDIR)\pxop.sbr"
	-@erase "$(INTDIR)\qrfactor.obj"
	-@erase "$(INTDIR)\qrfactor.sbr"
	-@erase "$(INTDIR)\schur.obj"
	-@erase "$(INTDIR)\schur.sbr"
	-@erase "$(INTDIR)\solve.obj"
	-@erase "$(INTDIR)\solve.sbr"
	-@erase "$(INTDIR)\sparse.obj"
	-@erase "$(INTDIR)\sparse.sbr"
	-@erase "$(INTDIR)\sparseio.obj"
	-@erase "$(INTDIR)\sparseio.sbr"
	-@erase "$(INTDIR)\spbkp.obj"
	-@erase "$(INTDIR)\spbkp.sbr"
	-@erase "$(INTDIR)\spchfctr.obj"
	-@erase "$(INTDIR)\spchfctr.sbr"
	-@erase "$(INTDIR)\splufctr.obj"
	-@erase "$(INTDIR)\splufctr.sbr"
	-@erase "$(INTDIR)\sprow.obj"
	-@erase "$(INTDIR)\sprow.sbr"
	-@erase "$(INTDIR)\spswap.obj"
	-@erase "$(INTDIR)\spswap.sbr"
	-@erase "$(INTDIR)\submat.obj"
	-@erase "$(INTDIR)\submat.sbr"
	-@erase "$(INTDIR)\svd.obj"
	-@erase "$(INTDIR)\svd.sbr"
	-@erase "$(INTDIR)\symmeig.obj"
	-@erase "$(INTDIR)\symmeig.sbr"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\update.sbr"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\vecop.obj"
	-@erase "$(INTDIR)\vecop.sbr"
	-@erase "$(INTDIR)\version.obj"
	-@erase "$(INTDIR)\version.sbr"
	-@erase "$(OUTDIR)\meschachd.bsc"
	-@erase "$(OUTDIR)\meschachd.dll"
	-@erase "$(OUTDIR)\meschachd.exp"
	-@erase "$(OUTDIR)\meschachd.ilk"
	-@erase "$(OUTDIR)\meschachd.lib"
	-@erase "$(OUTDIR)\meschachd.pdb"
	-@erase ".\meschachd.h"
	-@erase ".\meschachd.tlb"
	-@erase ".\meschachd_i.c"
	-@erase ".\Debug\regsvr32.trg"

CPP=cl.exe
CPP_PROJ=/nologo /MTd /W4 /Gm /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL HAVE_CONFIG_H ANSI_C HAVE_COMPLEX" /FR"$(INTDIR)\\" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

MTL=midl.exe
MTL_PROJ=
RSC=rc.exe
RSC_PROJ=/l 0x409 /fo"$(INTDIR)\meschachd.res" /d "_DEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\meschachd.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\bdfactor.sbr" \
	"$(INTDIR)\bkpfacto.sbr" \
	"$(INTDIR)\chfactor.sbr" \
	"$(INTDIR)\copy.sbr" \
	"$(INTDIR)\err.sbr" \
	"$(INTDIR)\extras.sbr" \
	"$(INTDIR)\fft.sbr" \
	"$(INTDIR)\givens.sbr" \
	"$(INTDIR)\hessen.sbr" \
	"$(INTDIR)\hsehldr.sbr" \
	"$(INTDIR)\init.sbr" \
	"$(INTDIR)\iter0.sbr" \
	"$(INTDIR)\iternsym.sbr" \
	"$(INTDIR)\itersym.sbr" \
	"$(INTDIR)\ivecop.sbr" \
	"$(INTDIR)\lufactor.sbr" \
	"$(INTDIR)\machine.sbr" \
	"$(INTDIR)\matlab.sbr" \
	"$(INTDIR)\matop.sbr" \
	"$(INTDIR)\matrixio.sbr" \
	"$(INTDIR)\meminfo.sbr" \
	"$(INTDIR)\memory.sbr" \
	"$(INTDIR)\memstat.sbr" \
	"$(INTDIR)\meschachd.sbr" \
	"$(INTDIR)\mfunc.sbr" \
	"$(INTDIR)\norm.sbr" \
	"$(INTDIR)\otherio.sbr" \
	"$(INTDIR)\pxop.sbr" \
	"$(INTDIR)\qrfactor.sbr" \
	"$(INTDIR)\schur.sbr" \
	"$(INTDIR)\solve.sbr" \
	"$(INTDIR)\sparse.sbr" \
	"$(INTDIR)\sparseio.sbr" \
	"$(INTDIR)\spbkp.sbr" \
	"$(INTDIR)\spchfctr.sbr" \
	"$(INTDIR)\splufctr.sbr" \
	"$(INTDIR)\sprow.sbr" \
	"$(INTDIR)\spswap.sbr" \
	"$(INTDIR)\submat.sbr" \
	"$(INTDIR)\svd.sbr" \
	"$(INTDIR)\symmeig.sbr" \
	"$(INTDIR)\update.sbr" \
	"$(INTDIR)\vecop.sbr" \
	"$(INTDIR)\version.sbr"

"$(OUTDIR)\meschachd.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /incremental:yes /pdb:"$(OUTDIR)\meschachd.pdb" /debug /machine:I386 /def:".\meschachd.def" /out:"$(OUTDIR)\meschachd.dll" /implib:"$(OUTDIR)\meschachd.lib" /pdbtype:sept 
DEF_FILE= \
	".\meschachd.def"
LINK32_OBJS= \
	"$(INTDIR)\bdfactor.obj" \
	"$(INTDIR)\bkpfacto.obj" \
	"$(INTDIR)\chfactor.obj" \
	"$(INTDIR)\copy.obj" \
	"$(INTDIR)\err.obj" \
	"$(INTDIR)\extras.obj" \
	"$(INTDIR)\fft.obj" \
	"$(INTDIR)\givens.obj" \
	"$(INTDIR)\hessen.obj" \
	"$(INTDIR)\hsehldr.obj" \
	"$(INTDIR)\init.obj" \
	"$(INTDIR)\iter0.obj" \
	"$(INTDIR)\iternsym.obj" \
	"$(INTDIR)\itersym.obj" \
	"$(INTDIR)\ivecop.obj" \
	"$(INTDIR)\lufactor.obj" \
	"$(INTDIR)\machine.obj" \
	"$(INTDIR)\matlab.obj" \
	"$(INTDIR)\matop.obj" \
	"$(INTDIR)\matrixio.obj" \
	"$(INTDIR)\meminfo.obj" \
	"$(INTDIR)\memory.obj" \
	"$(INTDIR)\memstat.obj" \
	"$(INTDIR)\meschachd.obj" \
	"$(INTDIR)\mfunc.obj" \
	"$(INTDIR)\norm.obj" \
	"$(INTDIR)\otherio.obj" \
	"$(INTDIR)\pxop.obj" \
	"$(INTDIR)\qrfactor.obj" \
	"$(INTDIR)\schur.obj" \
	"$(INTDIR)\solve.obj" \
	"$(INTDIR)\sparse.obj" \
	"$(INTDIR)\sparseio.obj" \
	"$(INTDIR)\spbkp.obj" \
	"$(INTDIR)\spchfctr.obj" \
	"$(INTDIR)\splufctr.obj" \
	"$(INTDIR)\sprow.obj" \
	"$(INTDIR)\spswap.obj" \
	"$(INTDIR)\submat.obj" \
	"$(INTDIR)\svd.obj" \
	"$(INTDIR)\symmeig.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\vecop.obj" \
	"$(INTDIR)\version.obj" \
	"$(INTDIR)\meschachd.res"

"$(OUTDIR)\meschachd.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

OutDir=.\Debug
TargetPath=.\Debug\meschachd.dll
InputPath=.\Debug\meschachd.dll
SOURCE="$(InputPath)"

"$(OUTDIR)\regsvr32.trg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	<<tempfile.bat 
	@echo off 
	regsvr32 /s /c "$(TargetPath)" 
	echo regsvr32 exec. time > "$(OutDir)\regsvr32.trg" 
<< 
	

!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"

OUTDIR=.\DebugU
INTDIR=.\DebugU
# Begin Custom Macros
OutDir=.\DebugU
# End Custom Macros

ALL : "$(OUTDIR)\meschachd.dll" ".\meschachd.tlb" ".\meschachd.h" ".\meschachd_i.c" ".\DebugU\regsvr32.trg"


CLEAN :
	-@erase "$(INTDIR)\bdfactor.obj"
	-@erase "$(INTDIR)\bkpfacto.obj"
	-@erase "$(INTDIR)\chfactor.obj"
	-@erase "$(INTDIR)\copy.obj"
	-@erase "$(INTDIR)\err.obj"
	-@erase "$(INTDIR)\extras.obj"
	-@erase "$(INTDIR)\fft.obj"
	-@erase "$(INTDIR)\givens.obj"
	-@erase "$(INTDIR)\hessen.obj"
	-@erase "$(INTDIR)\hsehldr.obj"
	-@erase "$(INTDIR)\init.obj"
	-@erase "$(INTDIR)\iter0.obj"
	-@erase "$(INTDIR)\iternsym.obj"
	-@erase "$(INTDIR)\itersym.obj"
	-@erase "$(INTDIR)\ivecop.obj"
	-@erase "$(INTDIR)\lufactor.obj"
	-@erase "$(INTDIR)\machine.obj"
	-@erase "$(INTDIR)\matlab.obj"
	-@erase "$(INTDIR)\matop.obj"
	-@erase "$(INTDIR)\matrixio.obj"
	-@erase "$(INTDIR)\meminfo.obj"
	-@erase "$(INTDIR)\memory.obj"
	-@erase "$(INTDIR)\memstat.obj"
	-@erase "$(INTDIR)\meschachd.obj"
	-@erase "$(INTDIR)\meschachd.res"
	-@erase "$(INTDIR)\mfunc.obj"
	-@erase "$(INTDIR)\norm.obj"
	-@erase "$(INTDIR)\otherio.obj"
	-@erase "$(INTDIR)\pxop.obj"
	-@erase "$(INTDIR)\qrfactor.obj"
	-@erase "$(INTDIR)\schur.obj"
	-@erase "$(INTDIR)\solve.obj"
	-@erase "$(INTDIR)\sparse.obj"
	-@erase "$(INTDIR)\sparseio.obj"
	-@erase "$(INTDIR)\spbkp.obj"
	-@erase "$(INTDIR)\spchfctr.obj"
	-@erase "$(INTDIR)\splufctr.obj"
	-@erase "$(INTDIR)\sprow.obj"
	-@erase "$(INTDIR)\spswap.obj"
	-@erase "$(INTDIR)\submat.obj"
	-@erase "$(INTDIR)\svd.obj"
	-@erase "$(INTDIR)\symmeig.obj"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\vecop.obj"
	-@erase "$(INTDIR)\version.obj"
	-@erase "$(OUTDIR)\meschachd.dll"
	-@erase "$(OUTDIR)\meschachd.exp"
	-@erase "$(OUTDIR)\meschachd.ilk"
	-@erase "$(OUTDIR)\meschachd.lib"
	-@erase "$(OUTDIR)\meschachd.pdb"
	-@erase ".\meschachd.h"
	-@erase ".\meschachd.tlb"
	-@erase ".\meschachd_i.c"
	-@erase ".\DebugU\regsvr32.trg"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MTd /W3 /Gm /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_USRDLL" /D "_UNICODE" /Fp"$(INTDIR)\meschachd.pch" /Yu"stdafx.h" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

MTL=midl.exe
MTL_PROJ=
RSC=rc.exe
RSC_PROJ=/l 0x409 /fo"$(INTDIR)\meschachd.res" /d "_DEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\meschachd.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /incremental:yes /pdb:"$(OUTDIR)\meschachd.pdb" /debug /machine:I386 /def:".\meschachd.def" /out:"$(OUTDIR)\meschachd.dll" /implib:"$(OUTDIR)\meschachd.lib" /pdbtype:sept 
DEF_FILE= \
	".\meschachd.def"
LINK32_OBJS= \
	"$(INTDIR)\bdfactor.obj" \
	"$(INTDIR)\bkpfacto.obj" \
	"$(INTDIR)\chfactor.obj" \
	"$(INTDIR)\copy.obj" \
	"$(INTDIR)\err.obj" \
	"$(INTDIR)\extras.obj" \
	"$(INTDIR)\fft.obj" \
	"$(INTDIR)\givens.obj" \
	"$(INTDIR)\hessen.obj" \
	"$(INTDIR)\hsehldr.obj" \
	"$(INTDIR)\init.obj" \
	"$(INTDIR)\iter0.obj" \
	"$(INTDIR)\iternsym.obj" \
	"$(INTDIR)\itersym.obj" \
	"$(INTDIR)\ivecop.obj" \
	"$(INTDIR)\lufactor.obj" \
	"$(INTDIR)\machine.obj" \
	"$(INTDIR)\matlab.obj" \
	"$(INTDIR)\matop.obj" \
	"$(INTDIR)\matrixio.obj" \
	"$(INTDIR)\meminfo.obj" \
	"$(INTDIR)\memory.obj" \
	"$(INTDIR)\memstat.obj" \
	"$(INTDIR)\meschachd.obj" \
	"$(INTDIR)\mfunc.obj" \
	"$(INTDIR)\norm.obj" \
	"$(INTDIR)\otherio.obj" \
	"$(INTDIR)\pxop.obj" \
	"$(INTDIR)\qrfactor.obj" \
	"$(INTDIR)\schur.obj" \
	"$(INTDIR)\solve.obj" \
	"$(INTDIR)\sparse.obj" \
	"$(INTDIR)\sparseio.obj" \
	"$(INTDIR)\spbkp.obj" \
	"$(INTDIR)\spchfctr.obj" \
	"$(INTDIR)\splufctr.obj" \
	"$(INTDIR)\sprow.obj" \
	"$(INTDIR)\spswap.obj" \
	"$(INTDIR)\submat.obj" \
	"$(INTDIR)\svd.obj" \
	"$(INTDIR)\symmeig.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\vecop.obj" \
	"$(INTDIR)\version.obj" \
	"$(INTDIR)\meschachd.res"

"$(OUTDIR)\meschachd.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

OutDir=.\DebugU
TargetPath=.\DebugU\meschachd.dll
InputPath=.\DebugU\meschachd.dll
SOURCE="$(InputPath)"

"$(OUTDIR)\regsvr32.trg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	<<tempfile.bat 
	@echo off 
	if "%OS%"=="" goto NOTNT 
	if not "%OS%"=="Windows_NT" goto NOTNT 
	regsvr32 /s /c "$(TargetPath)" 
	echo regsvr32 exec. time > "$(OutDir)\regsvr32.trg" 
	goto end 
	:NOTNT 
	echo Warning : Cannot register Unicode DLL on Windows 95 
	:end 
<< 
	

!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"

OUTDIR=.\ReleaseMinSize
INTDIR=.\ReleaseMinSize
# Begin Custom Macros
OutDir=.\ReleaseMinSize
# End Custom Macros

ALL : "$(OUTDIR)\meschachd.dll" ".\meschachd.tlb" ".\meschachd.h" ".\meschachd_i.c" ".\ReleaseMinSize\regsvr32.trg"


CLEAN :
	-@erase "$(INTDIR)\bdfactor.obj"
	-@erase "$(INTDIR)\bkpfacto.obj"
	-@erase "$(INTDIR)\chfactor.obj"
	-@erase "$(INTDIR)\copy.obj"
	-@erase "$(INTDIR)\err.obj"
	-@erase "$(INTDIR)\extras.obj"
	-@erase "$(INTDIR)\fft.obj"
	-@erase "$(INTDIR)\givens.obj"
	-@erase "$(INTDIR)\hessen.obj"
	-@erase "$(INTDIR)\hsehldr.obj"
	-@erase "$(INTDIR)\init.obj"
	-@erase "$(INTDIR)\iter0.obj"
	-@erase "$(INTDIR)\iternsym.obj"
	-@erase "$(INTDIR)\itersym.obj"
	-@erase "$(INTDIR)\ivecop.obj"
	-@erase "$(INTDIR)\lufactor.obj"
	-@erase "$(INTDIR)\machine.obj"
	-@erase "$(INTDIR)\matlab.obj"
	-@erase "$(INTDIR)\matop.obj"
	-@erase "$(INTDIR)\matrixio.obj"
	-@erase "$(INTDIR)\meminfo.obj"
	-@erase "$(INTDIR)\memory.obj"
	-@erase "$(INTDIR)\memstat.obj"
	-@erase "$(INTDIR)\meschachd.obj"
	-@erase "$(INTDIR)\meschachd.res"
	-@erase "$(INTDIR)\mfunc.obj"
	-@erase "$(INTDIR)\norm.obj"
	-@erase "$(INTDIR)\otherio.obj"
	-@erase "$(INTDIR)\pxop.obj"
	-@erase "$(INTDIR)\qrfactor.obj"
	-@erase "$(INTDIR)\schur.obj"
	-@erase "$(INTDIR)\solve.obj"
	-@erase "$(INTDIR)\sparse.obj"
	-@erase "$(INTDIR)\sparseio.obj"
	-@erase "$(INTDIR)\spbkp.obj"
	-@erase "$(INTDIR)\spchfctr.obj"
	-@erase "$(INTDIR)\splufctr.obj"
	-@erase "$(INTDIR)\sprow.obj"
	-@erase "$(INTDIR)\spswap.obj"
	-@erase "$(INTDIR)\submat.obj"
	-@erase "$(INTDIR)\svd.obj"
	-@erase "$(INTDIR)\symmeig.obj"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vecop.obj"
	-@erase "$(INTDIR)\version.obj"
	-@erase "$(OUTDIR)\meschachd.dll"
	-@erase "$(OUTDIR)\meschachd.exp"
	-@erase "$(OUTDIR)\meschachd.lib"
	-@erase ".\meschachd.h"
	-@erase ".\meschachd.tlb"
	-@erase ".\meschachd_i.c"
	-@erase ".\ReleaseMinSize\regsvr32.trg"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MT /W3 /O1 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "_ATL_DLL" /D "_ATL_MIN_CRT" /Fp"$(INTDIR)\meschachd.pch" /Yu"stdafx.h" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

MTL=midl.exe
MTL_PROJ=
RSC=rc.exe
RSC_PROJ=/l 0x409 /fo"$(INTDIR)\meschachd.res" /d "NDEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\meschachd.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /incremental:no /pdb:"$(OUTDIR)\meschachd.pdb" /machine:I386 /def:".\meschachd.def" /out:"$(OUTDIR)\meschachd.dll" /implib:"$(OUTDIR)\meschachd.lib" 
DEF_FILE= \
	".\meschachd.def"
LINK32_OBJS= \
	"$(INTDIR)\bdfactor.obj" \
	"$(INTDIR)\bkpfacto.obj" \
	"$(INTDIR)\chfactor.obj" \
	"$(INTDIR)\copy.obj" \
	"$(INTDIR)\err.obj" \
	"$(INTDIR)\extras.obj" \
	"$(INTDIR)\fft.obj" \
	"$(INTDIR)\givens.obj" \
	"$(INTDIR)\hessen.obj" \
	"$(INTDIR)\hsehldr.obj" \
	"$(INTDIR)\init.obj" \
	"$(INTDIR)\iter0.obj" \
	"$(INTDIR)\iternsym.obj" \
	"$(INTDIR)\itersym.obj" \
	"$(INTDIR)\ivecop.obj" \
	"$(INTDIR)\lufactor.obj" \
	"$(INTDIR)\machine.obj" \
	"$(INTDIR)\matlab.obj" \
	"$(INTDIR)\matop.obj" \
	"$(INTDIR)\matrixio.obj" \
	"$(INTDIR)\meminfo.obj" \
	"$(INTDIR)\memory.obj" \
	"$(INTDIR)\memstat.obj" \
	"$(INTDIR)\meschachd.obj" \
	"$(INTDIR)\mfunc.obj" \
	"$(INTDIR)\norm.obj" \
	"$(INTDIR)\otherio.obj" \
	"$(INTDIR)\pxop.obj" \
	"$(INTDIR)\qrfactor.obj" \
	"$(INTDIR)\schur.obj" \
	"$(INTDIR)\solve.obj" \
	"$(INTDIR)\sparse.obj" \
	"$(INTDIR)\sparseio.obj" \
	"$(INTDIR)\spbkp.obj" \
	"$(INTDIR)\spchfctr.obj" \
	"$(INTDIR)\splufctr.obj" \
	"$(INTDIR)\sprow.obj" \
	"$(INTDIR)\spswap.obj" \
	"$(INTDIR)\submat.obj" \
	"$(INTDIR)\svd.obj" \
	"$(INTDIR)\symmeig.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\vecop.obj" \
	"$(INTDIR)\version.obj" \
	"$(INTDIR)\meschachd.res"

"$(OUTDIR)\meschachd.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

OutDir=.\ReleaseMinSize
TargetPath=.\ReleaseMinSize\meschachd.dll
InputPath=.\ReleaseMinSize\meschachd.dll
SOURCE="$(InputPath)"

"$(OUTDIR)\regsvr32.trg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	<<tempfile.bat 
	@echo off 
	regsvr32 /s /c "$(TargetPath)" 
	echo regsvr32 exec. time > "$(OutDir)\regsvr32.trg" 
<< 
	

!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"

OUTDIR=.\ReleaseMinDependency
INTDIR=.\ReleaseMinDependency
# Begin Custom Macros
OutDir=.\ReleaseMinDependency
# End Custom Macros

ALL : "$(OUTDIR)\meschachd.dll" ".\meschachd.tlb" ".\meschachd.h" ".\meschachd_i.c" ".\ReleaseMinDependency\regsvr32.trg"


CLEAN :
	-@erase "$(INTDIR)\bdfactor.obj"
	-@erase "$(INTDIR)\bkpfacto.obj"
	-@erase "$(INTDIR)\chfactor.obj"
	-@erase "$(INTDIR)\copy.obj"
	-@erase "$(INTDIR)\err.obj"
	-@erase "$(INTDIR)\extras.obj"
	-@erase "$(INTDIR)\fft.obj"
	-@erase "$(INTDIR)\givens.obj"
	-@erase "$(INTDIR)\hessen.obj"
	-@erase "$(INTDIR)\hsehldr.obj"
	-@erase "$(INTDIR)\init.obj"
	-@erase "$(INTDIR)\iter0.obj"
	-@erase "$(INTDIR)\iternsym.obj"
	-@erase "$(INTDIR)\itersym.obj"
	-@erase "$(INTDIR)\ivecop.obj"
	-@erase "$(INTDIR)\lufactor.obj"
	-@erase "$(INTDIR)\machine.obj"
	-@erase "$(INTDIR)\matlab.obj"
	-@erase "$(INTDIR)\matop.obj"
	-@erase "$(INTDIR)\matrixio.obj"
	-@erase "$(INTDIR)\meminfo.obj"
	-@erase "$(INTDIR)\memory.obj"
	-@erase "$(INTDIR)\memstat.obj"
	-@erase "$(INTDIR)\meschachd.obj"
	-@erase "$(INTDIR)\meschachd.res"
	-@erase "$(INTDIR)\mfunc.obj"
	-@erase "$(INTDIR)\norm.obj"
	-@erase "$(INTDIR)\otherio.obj"
	-@erase "$(INTDIR)\pxop.obj"
	-@erase "$(INTDIR)\qrfactor.obj"
	-@erase "$(INTDIR)\schur.obj"
	-@erase "$(INTDIR)\solve.obj"
	-@erase "$(INTDIR)\sparse.obj"
	-@erase "$(INTDIR)\sparseio.obj"
	-@erase "$(INTDIR)\spbkp.obj"
	-@erase "$(INTDIR)\spchfctr.obj"
	-@erase "$(INTDIR)\splufctr.obj"
	-@erase "$(INTDIR)\sprow.obj"
	-@erase "$(INTDIR)\spswap.obj"
	-@erase "$(INTDIR)\submat.obj"
	-@erase "$(INTDIR)\svd.obj"
	-@erase "$(INTDIR)\symmeig.obj"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vecop.obj"
	-@erase "$(INTDIR)\version.obj"
	-@erase "$(OUTDIR)\meschachd.dll"
	-@erase "$(OUTDIR)\meschachd.exp"
	-@erase "$(OUTDIR)\meschachd.lib"
	-@erase ".\meschachd.h"
	-@erase ".\meschachd.tlb"
	-@erase ".\meschachd_i.c"
	-@erase ".\ReleaseMinDependency\regsvr32.trg"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MT /W3 /O1 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "_ATL_STATIC_REGISTRY" /D "_ATL_MIN_CRT" /Fp"$(INTDIR)\meschachd.pch" /Yu"stdafx.h" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

MTL=midl.exe
MTL_PROJ=
RSC=rc.exe
RSC_PROJ=/l 0x409 /fo"$(INTDIR)\meschachd.res" /d "NDEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\meschachd.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /incremental:no /pdb:"$(OUTDIR)\meschachd.pdb" /machine:I386 /def:".\meschachd.def" /out:"$(OUTDIR)\meschachd.dll" /implib:"$(OUTDIR)\meschachd.lib" 
DEF_FILE= \
	".\meschachd.def"
LINK32_OBJS= \
	"$(INTDIR)\bdfactor.obj" \
	"$(INTDIR)\bkpfacto.obj" \
	"$(INTDIR)\chfactor.obj" \
	"$(INTDIR)\copy.obj" \
	"$(INTDIR)\err.obj" \
	"$(INTDIR)\extras.obj" \
	"$(INTDIR)\fft.obj" \
	"$(INTDIR)\givens.obj" \
	"$(INTDIR)\hessen.obj" \
	"$(INTDIR)\hsehldr.obj" \
	"$(INTDIR)\init.obj" \
	"$(INTDIR)\iter0.obj" \
	"$(INTDIR)\iternsym.obj" \
	"$(INTDIR)\itersym.obj" \
	"$(INTDIR)\ivecop.obj" \
	"$(INTDIR)\lufactor.obj" \
	"$(INTDIR)\machine.obj" \
	"$(INTDIR)\matlab.obj" \
	"$(INTDIR)\matop.obj" \
	"$(INTDIR)\matrixio.obj" \
	"$(INTDIR)\meminfo.obj" \
	"$(INTDIR)\memory.obj" \
	"$(INTDIR)\memstat.obj" \
	"$(INTDIR)\meschachd.obj" \
	"$(INTDIR)\mfunc.obj" \
	"$(INTDIR)\norm.obj" \
	"$(INTDIR)\otherio.obj" \
	"$(INTDIR)\pxop.obj" \
	"$(INTDIR)\qrfactor.obj" \
	"$(INTDIR)\schur.obj" \
	"$(INTDIR)\solve.obj" \
	"$(INTDIR)\sparse.obj" \
	"$(INTDIR)\sparseio.obj" \
	"$(INTDIR)\spbkp.obj" \
	"$(INTDIR)\spchfctr.obj" \
	"$(INTDIR)\splufctr.obj" \
	"$(INTDIR)\sprow.obj" \
	"$(INTDIR)\spswap.obj" \
	"$(INTDIR)\submat.obj" \
	"$(INTDIR)\svd.obj" \
	"$(INTDIR)\symmeig.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\vecop.obj" \
	"$(INTDIR)\version.obj" \
	"$(INTDIR)\meschachd.res"

"$(OUTDIR)\meschachd.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

OutDir=.\ReleaseMinDependency
TargetPath=.\ReleaseMinDependency\meschachd.dll
InputPath=.\ReleaseMinDependency\meschachd.dll
SOURCE="$(InputPath)"

"$(OUTDIR)\regsvr32.trg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	<<tempfile.bat 
	@echo off 
	regsvr32 /s /c "$(TargetPath)" 
	echo regsvr32 exec. time > "$(OutDir)\regsvr32.trg" 
<< 
	

!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"

OUTDIR=.\ReleaseUMinSize
INTDIR=.\ReleaseUMinSize
# Begin Custom Macros
OutDir=.\ReleaseUMinSize
# End Custom Macros

ALL : "$(OUTDIR)\meschachd.dll" ".\meschachd.tlb" ".\meschachd.h" ".\meschachd_i.c" ".\ReleaseUMinSize\regsvr32.trg"


CLEAN :
	-@erase "$(INTDIR)\bdfactor.obj"
	-@erase "$(INTDIR)\bkpfacto.obj"
	-@erase "$(INTDIR)\chfactor.obj"
	-@erase "$(INTDIR)\copy.obj"
	-@erase "$(INTDIR)\err.obj"
	-@erase "$(INTDIR)\extras.obj"
	-@erase "$(INTDIR)\fft.obj"
	-@erase "$(INTDIR)\givens.obj"
	-@erase "$(INTDIR)\hessen.obj"
	-@erase "$(INTDIR)\hsehldr.obj"
	-@erase "$(INTDIR)\init.obj"
	-@erase "$(INTDIR)\iter0.obj"
	-@erase "$(INTDIR)\iternsym.obj"
	-@erase "$(INTDIR)\itersym.obj"
	-@erase "$(INTDIR)\ivecop.obj"
	-@erase "$(INTDIR)\lufactor.obj"
	-@erase "$(INTDIR)\machine.obj"
	-@erase "$(INTDIR)\matlab.obj"
	-@erase "$(INTDIR)\matop.obj"
	-@erase "$(INTDIR)\matrixio.obj"
	-@erase "$(INTDIR)\meminfo.obj"
	-@erase "$(INTDIR)\memory.obj"
	-@erase "$(INTDIR)\memstat.obj"
	-@erase "$(INTDIR)\meschachd.obj"
	-@erase "$(INTDIR)\meschachd.res"
	-@erase "$(INTDIR)\mfunc.obj"
	-@erase "$(INTDIR)\norm.obj"
	-@erase "$(INTDIR)\otherio.obj"
	-@erase "$(INTDIR)\pxop.obj"
	-@erase "$(INTDIR)\qrfactor.obj"
	-@erase "$(INTDIR)\schur.obj"
	-@erase "$(INTDIR)\solve.obj"
	-@erase "$(INTDIR)\sparse.obj"
	-@erase "$(INTDIR)\sparseio.obj"
	-@erase "$(INTDIR)\spbkp.obj"
	-@erase "$(INTDIR)\spchfctr.obj"
	-@erase "$(INTDIR)\splufctr.obj"
	-@erase "$(INTDIR)\sprow.obj"
	-@erase "$(INTDIR)\spswap.obj"
	-@erase "$(INTDIR)\submat.obj"
	-@erase "$(INTDIR)\svd.obj"
	-@erase "$(INTDIR)\symmeig.obj"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vecop.obj"
	-@erase "$(INTDIR)\version.obj"
	-@erase "$(OUTDIR)\meschachd.dll"
	-@erase "$(OUTDIR)\meschachd.exp"
	-@erase "$(OUTDIR)\meschachd.lib"
	-@erase ".\meschachd.h"
	-@erase ".\meschachd.tlb"
	-@erase ".\meschachd_i.c"
	-@erase ".\ReleaseUMinSize\regsvr32.trg"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MT /W3 /O1 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_USRDLL" /D "_UNICODE" /D "_ATL_DLL" /D "_ATL_MIN_CRT" /Fp"$(INTDIR)\meschachd.pch" /Yu"stdafx.h" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

MTL=midl.exe
MTL_PROJ=
RSC=rc.exe
RSC_PROJ=/l 0x409 /fo"$(INTDIR)\meschachd.res" /d "NDEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\meschachd.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /incremental:no /pdb:"$(OUTDIR)\meschachd.pdb" /machine:I386 /def:".\meschachd.def" /out:"$(OUTDIR)\meschachd.dll" /implib:"$(OUTDIR)\meschachd.lib" 
DEF_FILE= \
	".\meschachd.def"
LINK32_OBJS= \
	"$(INTDIR)\bdfactor.obj" \
	"$(INTDIR)\bkpfacto.obj" \
	"$(INTDIR)\chfactor.obj" \
	"$(INTDIR)\copy.obj" \
	"$(INTDIR)\err.obj" \
	"$(INTDIR)\extras.obj" \
	"$(INTDIR)\fft.obj" \
	"$(INTDIR)\givens.obj" \
	"$(INTDIR)\hessen.obj" \
	"$(INTDIR)\hsehldr.obj" \
	"$(INTDIR)\init.obj" \
	"$(INTDIR)\iter0.obj" \
	"$(INTDIR)\iternsym.obj" \
	"$(INTDIR)\itersym.obj" \
	"$(INTDIR)\ivecop.obj" \
	"$(INTDIR)\lufactor.obj" \
	"$(INTDIR)\machine.obj" \
	"$(INTDIR)\matlab.obj" \
	"$(INTDIR)\matop.obj" \
	"$(INTDIR)\matrixio.obj" \
	"$(INTDIR)\meminfo.obj" \
	"$(INTDIR)\memory.obj" \
	"$(INTDIR)\memstat.obj" \
	"$(INTDIR)\meschachd.obj" \
	"$(INTDIR)\mfunc.obj" \
	"$(INTDIR)\norm.obj" \
	"$(INTDIR)\otherio.obj" \
	"$(INTDIR)\pxop.obj" \
	"$(INTDIR)\qrfactor.obj" \
	"$(INTDIR)\schur.obj" \
	"$(INTDIR)\solve.obj" \
	"$(INTDIR)\sparse.obj" \
	"$(INTDIR)\sparseio.obj" \
	"$(INTDIR)\spbkp.obj" \
	"$(INTDIR)\spchfctr.obj" \
	"$(INTDIR)\splufctr.obj" \
	"$(INTDIR)\sprow.obj" \
	"$(INTDIR)\spswap.obj" \
	"$(INTDIR)\submat.obj" \
	"$(INTDIR)\svd.obj" \
	"$(INTDIR)\symmeig.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\vecop.obj" \
	"$(INTDIR)\version.obj" \
	"$(INTDIR)\meschachd.res"

"$(OUTDIR)\meschachd.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

OutDir=.\ReleaseUMinSize
TargetPath=.\ReleaseUMinSize\meschachd.dll
InputPath=.\ReleaseUMinSize\meschachd.dll
SOURCE="$(InputPath)"

"$(OUTDIR)\regsvr32.trg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	<<tempfile.bat 
	@echo off 
	if "%OS%"=="" goto NOTNT 
	if not "%OS%"=="Windows_NT" goto NOTNT 
	regsvr32 /s /c "$(TargetPath)" 
	echo regsvr32 exec. time > "$(OutDir)\regsvr32.trg" 
	goto end 
	:NOTNT 
	echo Warning : Cannot register Unicode DLL on Windows 95 
	:end 
<< 
	

!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"

OUTDIR=.\ReleaseUMinDependency
INTDIR=.\ReleaseUMinDependency
# Begin Custom Macros
OutDir=.\ReleaseUMinDependency
# End Custom Macros

ALL : "$(OUTDIR)\meschachd.dll" ".\meschachd.tlb" ".\meschachd.h" ".\meschachd_i.c" ".\ReleaseUMinDependency\regsvr32.trg"


CLEAN :
	-@erase "$(INTDIR)\bdfactor.obj"
	-@erase "$(INTDIR)\bkpfacto.obj"
	-@erase "$(INTDIR)\chfactor.obj"
	-@erase "$(INTDIR)\copy.obj"
	-@erase "$(INTDIR)\err.obj"
	-@erase "$(INTDIR)\extras.obj"
	-@erase "$(INTDIR)\fft.obj"
	-@erase "$(INTDIR)\givens.obj"
	-@erase "$(INTDIR)\hessen.obj"
	-@erase "$(INTDIR)\hsehldr.obj"
	-@erase "$(INTDIR)\init.obj"
	-@erase "$(INTDIR)\iter0.obj"
	-@erase "$(INTDIR)\iternsym.obj"
	-@erase "$(INTDIR)\itersym.obj"
	-@erase "$(INTDIR)\ivecop.obj"
	-@erase "$(INTDIR)\lufactor.obj"
	-@erase "$(INTDIR)\machine.obj"
	-@erase "$(INTDIR)\matlab.obj"
	-@erase "$(INTDIR)\matop.obj"
	-@erase "$(INTDIR)\matrixio.obj"
	-@erase "$(INTDIR)\meminfo.obj"
	-@erase "$(INTDIR)\memory.obj"
	-@erase "$(INTDIR)\memstat.obj"
	-@erase "$(INTDIR)\meschachd.obj"
	-@erase "$(INTDIR)\meschachd.res"
	-@erase "$(INTDIR)\mfunc.obj"
	-@erase "$(INTDIR)\norm.obj"
	-@erase "$(INTDIR)\otherio.obj"
	-@erase "$(INTDIR)\pxop.obj"
	-@erase "$(INTDIR)\qrfactor.obj"
	-@erase "$(INTDIR)\schur.obj"
	-@erase "$(INTDIR)\solve.obj"
	-@erase "$(INTDIR)\sparse.obj"
	-@erase "$(INTDIR)\sparseio.obj"
	-@erase "$(INTDIR)\spbkp.obj"
	-@erase "$(INTDIR)\spchfctr.obj"
	-@erase "$(INTDIR)\splufctr.obj"
	-@erase "$(INTDIR)\sprow.obj"
	-@erase "$(INTDIR)\spswap.obj"
	-@erase "$(INTDIR)\submat.obj"
	-@erase "$(INTDIR)\svd.obj"
	-@erase "$(INTDIR)\symmeig.obj"
	-@erase "$(INTDIR)\update.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vecop.obj"
	-@erase "$(INTDIR)\version.obj"
	-@erase "$(OUTDIR)\meschachd.dll"
	-@erase "$(OUTDIR)\meschachd.exp"
	-@erase "$(OUTDIR)\meschachd.lib"
	-@erase ".\meschachd.h"
	-@erase ".\meschachd.tlb"
	-@erase ".\meschachd_i.c"
	-@erase ".\ReleaseUMinDependency\regsvr32.trg"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

CPP=cl.exe
CPP_PROJ=/nologo /MT /W3 /O1 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_USRDLL" /D "_UNICODE" /D "_ATL_STATIC_REGISTRY" /D "_ATL_MIN_CRT" /Fp"$(INTDIR)\meschachd.pch" /Yu"stdafx.h" /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

MTL=midl.exe
MTL_PROJ=
RSC=rc.exe
RSC_PROJ=/l 0x409 /fo"$(INTDIR)\meschachd.res" /d "NDEBUG" 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\meschachd.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /dll /incremental:no /pdb:"$(OUTDIR)\meschachd.pdb" /machine:I386 /def:".\meschachd.def" /out:"$(OUTDIR)\meschachd.dll" /implib:"$(OUTDIR)\meschachd.lib" 
DEF_FILE= \
	".\meschachd.def"
LINK32_OBJS= \
	"$(INTDIR)\bdfactor.obj" \
	"$(INTDIR)\bkpfacto.obj" \
	"$(INTDIR)\chfactor.obj" \
	"$(INTDIR)\copy.obj" \
	"$(INTDIR)\err.obj" \
	"$(INTDIR)\extras.obj" \
	"$(INTDIR)\fft.obj" \
	"$(INTDIR)\givens.obj" \
	"$(INTDIR)\hessen.obj" \
	"$(INTDIR)\hsehldr.obj" \
	"$(INTDIR)\init.obj" \
	"$(INTDIR)\iter0.obj" \
	"$(INTDIR)\iternsym.obj" \
	"$(INTDIR)\itersym.obj" \
	"$(INTDIR)\ivecop.obj" \
	"$(INTDIR)\lufactor.obj" \
	"$(INTDIR)\machine.obj" \
	"$(INTDIR)\matlab.obj" \
	"$(INTDIR)\matop.obj" \
	"$(INTDIR)\matrixio.obj" \
	"$(INTDIR)\meminfo.obj" \
	"$(INTDIR)\memory.obj" \
	"$(INTDIR)\memstat.obj" \
	"$(INTDIR)\meschachd.obj" \
	"$(INTDIR)\mfunc.obj" \
	"$(INTDIR)\norm.obj" \
	"$(INTDIR)\otherio.obj" \
	"$(INTDIR)\pxop.obj" \
	"$(INTDIR)\qrfactor.obj" \
	"$(INTDIR)\schur.obj" \
	"$(INTDIR)\solve.obj" \
	"$(INTDIR)\sparse.obj" \
	"$(INTDIR)\sparseio.obj" \
	"$(INTDIR)\spbkp.obj" \
	"$(INTDIR)\spchfctr.obj" \
	"$(INTDIR)\splufctr.obj" \
	"$(INTDIR)\sprow.obj" \
	"$(INTDIR)\spswap.obj" \
	"$(INTDIR)\submat.obj" \
	"$(INTDIR)\svd.obj" \
	"$(INTDIR)\symmeig.obj" \
	"$(INTDIR)\update.obj" \
	"$(INTDIR)\vecop.obj" \
	"$(INTDIR)\version.obj" \
	"$(INTDIR)\meschachd.res"

"$(OUTDIR)\meschachd.dll" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

OutDir=.\ReleaseUMinDependency
TargetPath=.\ReleaseUMinDependency\meschachd.dll
InputPath=.\ReleaseUMinDependency\meschachd.dll
SOURCE="$(InputPath)"

"$(OUTDIR)\regsvr32.trg" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	<<tempfile.bat 
	@echo off 
	if "%OS%"=="" goto NOTNT 
	if not "%OS%"=="Windows_NT" goto NOTNT 
	regsvr32 /s /c "$(TargetPath)" 
	echo regsvr32 exec. time > "$(OutDir)\regsvr32.trg" 
	goto end 
	:NOTNT 
	echo Warning : Cannot register Unicode DLL on Windows 95 
	:end 
<< 
	

!ENDIF 


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("meschachd.dep")
!INCLUDE "meschachd.dep"
!ELSE 
!MESSAGE Warning: cannot find "meschachd.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "meschachd - Win32 Debug" || "$(CFG)" == "meschachd - Win32 Unicode Debug" || "$(CFG)" == "meschachd - Win32 Release MinSize" || "$(CFG)" == "meschachd - Win32 Release MinDependency" || "$(CFG)" == "meschachd - Win32 Unicode Release MinSize" || "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"
SOURCE=..\..\..\meschachd\visc\bdfactor.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\bdfactor.obj"	"$(INTDIR)\bdfactor.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\bdfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\bdfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\bdfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\bdfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\bdfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\bkpfacto.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\bkpfacto.obj"	"$(INTDIR)\bkpfacto.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\bkpfacto.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\bkpfacto.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\bkpfacto.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\bkpfacto.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\bkpfacto.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\chfactor.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\chfactor.obj"	"$(INTDIR)\chfactor.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\chfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\chfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\chfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\chfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\chfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\copy.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\copy.obj"	"$(INTDIR)\copy.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\copy.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\copy.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\copy.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\copy.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\copy.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\err.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\err.obj"	"$(INTDIR)\err.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\err.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\err.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\err.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\err.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\err.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\extras.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\extras.obj"	"$(INTDIR)\extras.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\extras.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\extras.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\extras.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\extras.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\extras.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\fft.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\fft.obj"	"$(INTDIR)\fft.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\fft.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\fft.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\fft.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\fft.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\fft.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\givens.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\givens.obj"	"$(INTDIR)\givens.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\givens.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\givens.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\givens.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\givens.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\givens.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\hessen.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\hessen.obj"	"$(INTDIR)\hessen.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\hessen.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\hessen.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\hessen.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\hessen.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\hessen.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\hsehldr.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\hsehldr.obj"	"$(INTDIR)\hsehldr.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\hsehldr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\hsehldr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\hsehldr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\hsehldr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\hsehldr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\init.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\init.obj"	"$(INTDIR)\init.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\init.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\init.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\init.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\init.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\init.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\iter0.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\iter0.obj"	"$(INTDIR)\iter0.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\iter0.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\iter0.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\iter0.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\iter0.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\iter0.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\iternsym.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\iternsym.obj"	"$(INTDIR)\iternsym.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\iternsym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\iternsym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\iternsym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\iternsym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\iternsym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\itersym.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\itersym.obj"	"$(INTDIR)\itersym.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\itersym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\itersym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\itersym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\itersym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\itersym.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\ivecop.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\ivecop.obj"	"$(INTDIR)\ivecop.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\ivecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\ivecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\ivecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\ivecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\ivecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\lufactor.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\lufactor.obj"	"$(INTDIR)\lufactor.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\lufactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\lufactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\lufactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\lufactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\lufactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\machine.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\machine.obj"	"$(INTDIR)\machine.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\machine.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\machine.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\machine.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\machine.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\machine.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\matlab.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\matlab.obj"	"$(INTDIR)\matlab.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\matlab.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\matlab.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\matlab.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\matlab.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\matlab.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\matop.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\matop.obj"	"$(INTDIR)\matop.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\matop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\matop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\matop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\matop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\matop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\matrixio.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\matrixio.obj"	"$(INTDIR)\matrixio.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\matrixio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\matrixio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\matrixio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\matrixio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\matrixio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\meminfo.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\meminfo.obj"	"$(INTDIR)\meminfo.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\meminfo.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\meminfo.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\meminfo.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\meminfo.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\meminfo.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\memory.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\memory.obj"	"$(INTDIR)\memory.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\memory.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\memory.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\memory.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\memory.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\memory.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\memstat.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\memstat.obj"	"$(INTDIR)\memstat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\memstat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\memstat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\memstat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\memstat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\memstat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=.\meschachd.cpp

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\meschachd.obj"	"$(INTDIR)\meschachd.sbr" : $(SOURCE) "$(INTDIR)" ".\meschachd_i.c" ".\meschachd.h"


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\meschachd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\meschachd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\meschachd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\meschachd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\meschachd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"


!ENDIF 

SOURCE=.\meschachd.idl

!IF  "$(CFG)" == "meschachd - Win32 Debug"

MTL_SWITCHES=/tlb ".\meschachd.tlb" /h "meschachd.h" /iid "meschachd_i.c" /Oicf 

".\meschachd.tlb"	".\meschachd.h"	".\meschachd_i.c" : $(SOURCE) "$(INTDIR)"
	$(MTL) @<<
  $(MTL_SWITCHES) $(SOURCE)
<<


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"

MTL_SWITCHES=/tlb ".\meschachd.tlb" /h "meschachd.h" /iid "meschachd_i.c" /Oicf 

".\meschachd.tlb"	".\meschachd.h"	".\meschachd_i.c" : $(SOURCE) "$(INTDIR)"
	$(MTL) @<<
  $(MTL_SWITCHES) $(SOURCE)
<<


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"

MTL_SWITCHES=/tlb ".\meschachd.tlb" /h "meschachd.h" /iid "meschachd_i.c" /Oicf 

".\meschachd.tlb"	".\meschachd.h"	".\meschachd_i.c" : $(SOURCE) "$(INTDIR)"
	$(MTL) @<<
  $(MTL_SWITCHES) $(SOURCE)
<<


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"

MTL_SWITCHES=/tlb ".\meschachd.tlb" /h "meschachd.h" /iid "meschachd_i.c" /Oicf 

".\meschachd.tlb"	".\meschachd.h"	".\meschachd_i.c" : $(SOURCE) "$(INTDIR)"
	$(MTL) @<<
  $(MTL_SWITCHES) $(SOURCE)
<<


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"

MTL_SWITCHES=/tlb ".\meschachd.tlb" /h "meschachd.h" /iid "meschachd_i.c" /Oicf 

".\meschachd.tlb"	".\meschachd.h"	".\meschachd_i.c" : $(SOURCE) "$(INTDIR)"
	$(MTL) @<<
  $(MTL_SWITCHES) $(SOURCE)
<<


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"

MTL_SWITCHES=/tlb ".\meschachd.tlb" /h "meschachd.h" /iid "meschachd_i.c" /Oicf 

".\meschachd.tlb"	".\meschachd.h"	".\meschachd_i.c" : $(SOURCE) "$(INTDIR)"
	$(MTL) @<<
  $(MTL_SWITCHES) $(SOURCE)
<<


!ENDIF 

SOURCE=.\meschachd.rc

"$(INTDIR)\meschachd.res" : $(SOURCE) "$(INTDIR)" ".\meschachd.tlb"
	$(RSC) $(RSC_PROJ) $(SOURCE)


SOURCE=..\..\..\meschachd\visc\mfunc.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\mfunc.obj"	"$(INTDIR)\mfunc.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\mfunc.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\mfunc.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\mfunc.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\mfunc.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\mfunc.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\norm.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\norm.obj"	"$(INTDIR)\norm.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\norm.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\norm.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\norm.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\norm.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\norm.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\otherio.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\otherio.obj"	"$(INTDIR)\otherio.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\otherio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\otherio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\otherio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\otherio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\otherio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\pxop.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\pxop.obj"	"$(INTDIR)\pxop.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\pxop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\pxop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\pxop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\pxop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\pxop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\qrfactor.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\qrfactor.obj"	"$(INTDIR)\qrfactor.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\qrfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\qrfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\qrfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\qrfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\qrfactor.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\schur.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\schur.obj"	"$(INTDIR)\schur.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\schur.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\schur.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\schur.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\schur.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\schur.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\solve.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\solve.obj"	"$(INTDIR)\solve.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\solve.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\solve.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\solve.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\solve.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\solve.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\sparse.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\sparse.obj"	"$(INTDIR)\sparse.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\sparse.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\sparse.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\sparse.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\sparse.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\sparse.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\sparseio.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\sparseio.obj"	"$(INTDIR)\sparseio.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\sparseio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\sparseio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\sparseio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\sparseio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\sparseio.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\spbkp.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\spbkp.obj"	"$(INTDIR)\spbkp.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\spbkp.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\spbkp.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\spbkp.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\spbkp.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\spbkp.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\spchfctr.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\spchfctr.obj"	"$(INTDIR)\spchfctr.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\spchfctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\spchfctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\spchfctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\spchfctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\spchfctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\splufctr.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\splufctr.obj"	"$(INTDIR)\splufctr.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\splufctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\splufctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\splufctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\splufctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\splufctr.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\sprow.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\sprow.obj"	"$(INTDIR)\sprow.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\sprow.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\sprow.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\sprow.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\sprow.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\sprow.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\spswap.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\spswap.obj"	"$(INTDIR)\spswap.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\spswap.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\spswap.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\spswap.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\spswap.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\spswap.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\submat.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\submat.obj"	"$(INTDIR)\submat.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\submat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\submat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\submat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\submat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\submat.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\svd.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\svd.obj"	"$(INTDIR)\svd.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\svd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\svd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\svd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\svd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\svd.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\symmeig.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\symmeig.obj"	"$(INTDIR)\symmeig.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\symmeig.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\symmeig.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\symmeig.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\symmeig.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\symmeig.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\update.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\update.obj"	"$(INTDIR)\update.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\update.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\update.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\update.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\update.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\update.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\vecop.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\vecop.obj"	"$(INTDIR)\vecop.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\vecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\vecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\vecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\vecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\vecop.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 

SOURCE=..\..\..\meschachd\visc\version.c

!IF  "$(CFG)" == "meschachd - Win32 Debug"


"$(INTDIR)\version.obj"	"$(INTDIR)\version.sbr" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Debug"


"$(INTDIR)\version.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinSize"


"$(INTDIR)\version.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Release MinDependency"


"$(INTDIR)\version.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinSize"


"$(INTDIR)\version.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "meschachd - Win32 Unicode Release MinDependency"


"$(INTDIR)\version.obj" : $(SOURCE) "$(INTDIR)" "$(INTDIR)\meschachd.pch"
	$(CPP) $(CPP_PROJ) $(SOURCE)


!ENDIF 


!ENDIF 

