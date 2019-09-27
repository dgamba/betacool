@echo off
set DICT=mad8.dic
if "%1" == "" goto noinput
del %2
mad8.exe<%1
if "%2" == "" goto endmad
if exist %2 goto endmad
ren print %2
goto endmad
:noinput
mad8.exe
:endmad
pause
