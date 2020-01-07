# Un-official version of BETACOOL code.

The present version is a copy of BETACOOL source code by A. Smirnov with corrections by N. Mityanina, A. Sidorin.

-- D. Gamba - Oct. 2019

In _betatrack_ branch I re-shuffled a bit the files around, and tried to build an adeguate makefile to build betacool (at least the core part)
The first steps to ordering the file structure was inspired from "#backup.bat" file by A. Smirnov. 
For reference, here the bash commands I used initially...:

```
> mkdir src
> mkdir GUI
> mkdir Dat
> mkdir Bat
> mv *.cpp *.h *.bpr *.dsp *.vcproj *.vcxproj *.cbx *.cbproj *.sln src/
> mv Betacool.vcxproj.user src/
> mv *.bat Bat/
> mv *.exe *.dll *.dfm *.top *.grf *.srf *.cur *.sur *.dic *.tpt GUI/
> mv *.bld *.mad *.lat *.tfs *.red *.err *.tvt *.lvt *.inj *.ela *.pat *.bar *.mov *.txt *.bump Dat/
> mv BRing-* Dat/
> mv survey Dat/
```

of course I had then to deviate, and start removing old log files trying to keep the essential.

Note that in the inital repository there were the follwing directoris:
- Betacool.tlog
- Debug
- Release
- ipch
- x64

## Modifications to Visual Studio project:
- since I moved all sources to the src/ folder, I had to add manually the "src\" relative path in Betacool.vcxproj file
- I used the "Build" directory also as "Intermediate Directory", which is then used to put logs of the compilation 
- Betacool.sln is (maybe) needed next to Betacool.vcxproj
- Betacool.v12.suo and Betacool.sdf are automatically generated, so probably not needed.

## For compilation:
In a new conda environment, I installed:
conda install -c anaconda gcc
conda install -c conda-forge openmp
conda install -c conda-forge stdlib-list
conda install libcxx

! To install stdlib:
cd /Library/Developer/CommandLineTools/Packages/
open macOS_SDK_headers_for_macOS_10.14.pkg 


