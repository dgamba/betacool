@echo off
:LOOP
xcopy o:\smirnov\betacool\nesr\*.cur /Y
xcopy o:\smirnov\betacool\nesr\betacool.war /Y
if "%1" == "" goto END
echo Press Ctrl-C to cancel or wait for %1 seconds...
delay.exe %1
goto LOOP
:END
