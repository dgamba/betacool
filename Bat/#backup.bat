@echo off
if exist ..\..\backup\betacool\#backup.bat goto EXIT
echo Do you want backup BETACOOL files? 
echo Press Ctrl-C to cancel 
pause
echo ------------- SRC ---------------
xcopy *.    backup\src\ /D /Y
xcopy *.cpp backup\src\ /D /Y
xcopy *.h   backup\src\ /D /Y
xcopy *.bpr backup\src\ /D /Y
xcopy *.dsp backup\src\ /D /Y
xcopy *.cbx backup\src\ /D /Y
xcopy *.vcproj backup\src\ /D /Y
xcopy *.cbproj backup\src\ /D /Y
xcopy *.sln backup\src\ /D /Y
echo ------------- GUI ---------------
xcopy *.bat backup\betacool\ /D /Y
xcopy *.exe backup\betacool\ /D /Y
xcopy *.dll backup\betacool\ /D /Y
xcopy *.dfm backup\betacool\ /D /Y
xcopy *.top backup\betacool\ /D /Y
xcopy *.grf backup\betacool\ /D /Y
xcopy *.srf backup\betacool\ /D /Y
xcopy *.cur backup\betacool\ /D /Y
xcopy *.sur backup\betacool\ /D /Y
xcopy *.dic backup\betacool\ /D /Y
xcopy *.tpt backup\betacool\ /D /Y
echo ------------- DAT ---------------              
xcopy *.bld backup\betacool%1\ /D /Y
xcopy *.mad backup\betacool%1\ /D /Y
xcopy *.lat backup\betacool%1\ /D /Y
xcopy *.tfs backup\betacool%1\ /D /Y
xcopy *.red backup\betacool%1\ /D /Y
xcopy *.err backup\betacool%1\ /D /Y
xcopy *.tvt backup\betacool%1\ /D /Y
xcopy *.lvt backup\betacool%1\ /D /Y
xcopy *.inj backup\betacool%1\ /D /Y
xcopy *.ela backup\betacool%1\ /D /Y
xcopy *.pat backup\betacool%1\ /D /Y
xcopy *.bar backup\betacool%1\ /D /Y
xcopy *.mov backup\betacool%1\ /D /Y
xcopy *.txt backup\betacool%1\ /D /Y
xcopy *.bump backup\betacool%1\ /D /Y
echo --------------------------------
echo BETACOOL files were saved
echo --------------------------------
pause
:EXIT