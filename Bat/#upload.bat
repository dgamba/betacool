@echo off
if not exist ..\..\..\betacool\backup\betacool\#upload.bat goto EXIT
echo Do you want upload BETACOOL files?
echo Press Ctrl-C to cancel 
pause
cd ..\..
echo ------------- SRC ---------------
xcopy backup\src\*.* /D /Y
echo ---------- GUI - DAT ------------
xcopy backup\betacool\*.* /D /Y
echo --------------------------------
echo BETACOOL files were uploaded
echo --------------------------------
pause
:EXIT