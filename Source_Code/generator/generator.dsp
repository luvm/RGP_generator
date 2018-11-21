# Microsoft Developer Studio Project File - Name="generator" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=generator - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "generator.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "generator.mak" CFG="generator - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "generator - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "generator - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "generator - Win32 Release"

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
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x419 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "generator - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "generator - Win32 Release"
# Name "generator - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\fortran\COMMONS.INI
# End Source File
# Begin Source File

SOURCE=..\fortran\COMTRN.INI
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_ANC.FOR
DEP_F90_CORE_=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_BWR.FOR
DEP_F90_CORE_B=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_CPP.FOR
DEP_F90_CORE_C=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_DE.FOR
DEP_F90_CORE_D=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_ECS.FOR
DEP_F90_CORE_E=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_FEQ.FOR
DEP_F90_CORE_F=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_MLT.FOR
DEP_F90_CORE_M=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_PH0.FOR
DEP_F90_CORE_P=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_PR.FOR
DEP_F90_CORE_PR=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_QUI.FOR
DEP_F90_CORE_Q=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\CORE_STN.FOR
DEP_F90_CORE_S=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\FLASH2.FOR
DEP_F90_FLASH=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\FLSH_SUB.FOR
DEP_F90_FLSH_=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\IDEALGAS.FOR
DEP_F90_IDEAL=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\MIX_AGA8.FOR
DEP_F90_MIX_A=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\MIX_HMX.FOR
DEP_F90_MIX_H=\
	"..\fortran\COMMONS.INI"\
	"..\fortran\COMTRN.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\PASS_FTN.FOR
DEP_F90_PASS_=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\PROP_SUB.FOR
DEP_F90_PROP_=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\REALGAS.FOR
DEP_F90_REALG=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\..\..\RGP_gen_2\RGPgenerator.F90
# End Source File
# Begin Source File

SOURCE=..\fortran\SAT_SUB.FOR
DEP_F90_SAT_S=\
	"..\fortran\COMMONS.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\SETUP.FOR
DEP_F90_SETUP=\
	"..\fortran\COMMONS.INI"\
	"..\fortran\COMTRN.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\SETUP2.FOR
DEP_F90_SETUP2=\
	"..\fortran\COMMONS.INI"\
	"..\fortran\COMTRN.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\TRNS_ECS.FOR
DEP_F90_TRNS_=\
	"..\fortran\COMMONS.INI"\
	"..\fortran\COMTRN.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\TRNS_TCX.FOR
DEP_F90_TRNS_T=\
	"..\fortran\COMMONS.INI"\
	"..\fortran\COMTRN.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\TRNS_VIS.FOR
DEP_F90_TRNS_V=\
	"..\fortran\COMMONS.INI"\
	"..\fortran\COMTRN.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\TRNSP.FOR
DEP_F90_TRNSP=\
	"..\fortran\COMMONS.INI"\
	"..\fortran\COMTRN.INI"\
	
# End Source File
# Begin Source File

SOURCE=..\fortran\UTILITY.FOR
DEP_F90_UTILI=\
	"..\fortran\COMMONS.INI"\
	"..\fortran\COMTRN.INI"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
