# Microsoft Developer Studio Project File - Name="Betacool" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Betacool - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Betacool.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Betacool.mak" CFG="Betacool - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Betacool - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Betacool - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Betacool - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Yu"stdafx.h" /FD /c
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x419 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "Betacool - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /Yu"stdafx.h" /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /Yu"stdafx.h" /FD /GZ /c
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Betacool - Win32 Release"
# Name "Betacool - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\Betacool.cpp
# End Source File
# Begin Source File

SOURCE=.\bolideU.cpp
# End Source File
# Begin Source File

SOURCE=.\dataU.cpp
# End Source File
# Begin Source File

SOURCE=.\doubleU.cpp
# End Source File
# Begin Source File

SOURCE=.\matrixU.cpp
# End Source File
# Begin Source File

SOURCE=.\StdAfx.cpp
# ADD CPP /Yc"stdafx.h"
# End Source File
# Begin Source File

SOURCE=.\xBeam.cpp
# End Source File
# Begin Source File

SOURCE=.\xDistributor.cpp
# End Source File
# Begin Source File

SOURCE=.\xDraw.cpp
# End Source File
# Begin Source File

SOURCE=.\xDynamic.cpp
# End Source File
# Begin Source File

SOURCE=.\xEbeam.cpp
# End Source File
# Begin Source File

SOURCE=.\xEcool.cpp
# End Source File
# Begin Source File

SOURCE=.\xEffect.cpp
# End Source File
# Begin Source File

SOURCE=.\xForce.cpp
# End Source File
# Begin Source File

SOURCE=.\xhiroshi.cpp
# End Source File
# Begin Source File

SOURCE=.\xIBS.cpp
# End Source File
# Begin Source File

SOURCE=.\xLibrary.cpp
# End Source File
# Begin Source File

SOURCE=.\xOptics.cpp
# End Source File
# Begin Source File

SOURCE=.\xPowell.cpp
# End Source File
# Begin Source File

SOURCE=.\xrestgas.cpp
# End Source File
# Begin Source File

SOURCE=.\xRing.cpp
# End Source File
# Begin Source File

SOURCE=.\xRunge.cpp
# End Source File
# Begin Source File

SOURCE=.\xstoch.cpp
# End Source File
# Begin Source File

SOURCE=.\xTarget.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\bolideU.h
# End Source File
# Begin Source File

SOURCE=.\BTempinc.h
# End Source File
# Begin Source File

SOURCE=.\BTemplat.h
# End Source File
# Begin Source File

SOURCE=.\dataU.h
# End Source File
# Begin Source File

SOURCE=.\doubleU.h
# End Source File
# Begin Source File

SOURCE=.\matrixU.h
# End Source File
# Begin Source File

SOURCE=.\StdAfx.h
# End Source File
# Begin Source File

SOURCE=.\xBeam.h
# End Source File
# Begin Source File

SOURCE=.\xDistributor.h
# End Source File
# Begin Source File

SOURCE=.\xDraw.h
# End Source File
# Begin Source File

SOURCE=.\xDynamic.h
# End Source File
# Begin Source File

SOURCE=.\xEbeam.h
# End Source File
# Begin Source File

SOURCE=.\xEcool.h
# End Source File
# Begin Source File

SOURCE=.\xEffect.h
# End Source File
# Begin Source File

SOURCE=.\xForce.h
# End Source File
# Begin Source File

SOURCE=.\xhiroshi.h
# End Source File
# Begin Source File

SOURCE=.\xIBS.h
# End Source File
# Begin Source File

SOURCE=.\xLibrary.h
# End Source File
# Begin Source File

SOURCE=.\xOptics.h
# End Source File
# Begin Source File

SOURCE=.\xPowell.h
# End Source File
# Begin Source File

SOURCE=.\xrestgas.h
# End Source File
# Begin Source File

SOURCE=.\xRing.h
# End Source File
# Begin Source File

SOURCE=.\xRunge.h
# End Source File
# Begin Source File

SOURCE=.\XSTOCH.H
# End Source File
# Begin Source File

SOURCE=.\xTarget.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\ReadMe.txt
# End Source File
# End Target
# End Project
