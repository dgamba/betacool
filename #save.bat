@echo off
if "%1" == "" goto EXIT
echo Use CTRL-C to stop of saving
mkdir #%1
copy %1 #%1 /-Y
copy betacool.war #%1 /Y
copy *.cur #%1 /Y
copy *.sur #%1 /Y
:EXIT