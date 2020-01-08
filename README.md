# Un-official version of BETACOOL code.

The present version is a copy of BETACOOL source code by A. Smirnov with corrections by N. Mityanina, A. Sidorin.

-- D. Gamba - Oct. 2019

---
# What is what

- README.md:         this file
- src:               all source files of BETACOOL (only!)
- Betacool.vcxproj:  Microsoft Visual Studio 2012 project file used to compile BETACOOL (only! no source for BOLIDE) for Windows (32 and 64 bit)
- Makefile:          UNIX Make file to compile BETACOOL (only! in Unix one has to deal with input files manually)
- BOLIDE:            minimal set of files to start BOLIDE interface, including the original "bolide.exe" executable for Win32
- Examples/BOLIDE:   a more extensive example of BOLIDE, with input/output files already generated for NICA (?) cooler
- Examples/BETACOOL: a minimalistic example for computing the cooling force with betacool from command line only 
- Docs:              folder with all documentation from original authors
- Dat:               folder with data files used to model previous coolers
- Bat:               folder with original *.bat files from authors
- Build:             normally and empty directory, it gets filled by temporary/log files when trying to compile BETACOOL for different architectures

---
# How to compile

## For compilation on Windows
Just load Betacool.vcxproj with Microsoft Visual Studio 2012 and "Build" it for your Windows architecture. 
This will create a betacool_<arc>.exe file 


## For compilation on UNIX
TODO, but will try to make it as simple as:
`
make
`
and this will create betacool_<arc> file executable on your architecture  


## For compilation on MacOSX:
In a new conda environment, I installed:
`
conda install -c anaconda gcc
conda install -c conda-forge openmp
conda install -c conda-forge stdlib-list
conda install libcxx
`

I installed the stdlib as:
`
cd /Library/Developer/CommandLineTools/Packages/
open macOS_SDK_headers_for_macOS_10.14.pkg 
`

Then something similar to UNIX compilation then....


---
# Logbook of things I did to obtain the present version

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

Note that in the inital repository there were the follwing directories:
- Betacool.tlog
- Debug
- Release
- ipch
- x64

## Modifications to Visual Studio project:
- since I moved all sources to the src/ folder, I had to add manually the "src\" relative path in Betacool.vcxproj file
- I used the "Build" directory also as "Intermediate Directory", which is then used to put logs of the compilation 
- several files could be cleaned up.
- didn't find the source code of the grafical interface. Will need to live with the binaries stored in BOLIDE folder for the time being.
- 


