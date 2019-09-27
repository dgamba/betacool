@echo off
echo Do you want backup PLOT files? 
echo Press Ctrl-C to cancel 
pause
:LOOP
xcopy h:\smirnov\betacool\nica\*.cur /Y
rem xcopy h:\smirnov\betacool\nica\betacool.war /Y
echo Wait for 5 minutes delay...
echo Press Ctrl-C to cancel 
delay.exe 300
goto LOOP
