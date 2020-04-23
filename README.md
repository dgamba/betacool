# Un-official version of BETACOOL code.

The present version is a copy of BETACOOL source code by A. Smirnov with corrections by N. Mityanina, A. Sidorin, with minor corrections to allow compiling it on UNIX.
The files have been re-shuffled to have more clean(?) structure.
Currently, it is able to compile in Windows with Microsoft Visual Studio 2013 starting from Betacool.vcxproj file.
For UNIX-like systems (Linux and MacOSX) it compiles using the provided Makefile.
For the time being no attempt to correct all warning during compilation has been made!

Note that here only the source code of BETACOOL is available. 
The source code of the Windows graphical interface "BOLIDE" is available only as a binary. 

Some details on how to use it are in [Docs/HowTo.md](./Docs/HowTo.md)

-- D. Gamba, A. Latina - Jan. 2020

---
# What is what

- README.md:         this file
- src:               all source files of BETACOOL (only!)
- Betacool.vcxproj:  Microsoft Visual Studio 2013 project file used to compile BETACOOL (only! no source for BOLIDE) for Windows (32 and 64 bit)
- Makefile:          UNIX Make file to compile BETACOOL (only! in Unix one has to deal with input files manually)
- BOLIDE:            minimal set of files to start BOLIDE interface, including the original "bolide.exe" executable for Win32
- Examples/BOLIDE:   a more extensive example of BOLIDE, with input/output files already generated for NICA (?) cooler
- Examples/BETACOOL: a minimalistic example for computing the cooling force with betacool from command line only 
- Examples/BETACOOL-ref: as Examples/BETACOOL, but including "reference" outputs generated with 2019 version of BETACOOL. 
- Docs:              folder with all documentation from original authors
- Dat:               folder with data files used to model previous coolers (not strictly needed)
- Bat:               folder with original *.bat files from authors (not strictly needed)
- Build:             normally and empty directory, it gets filled by temporary/log files when trying to compile BETACOOL for different architectures

---
# How to compile

## on Windows
Just load Betacool.vcxproj with Microsoft Visual Studio 2013 and "Build" it for your Windows architecture. 
This will create a betacool_<arc>.exe file 


## on UNIX
The provided Makefile should be good enough for a compilation on Linux and MacOSX, provided that the gcc and standard libraries are available.
To compile the code, simply execute:
```
make
```
and this will create Betacool_UNIX file executable on your architecture.
The compiled objects, will be placed into the "Build/UNIX" folder. They can be deleted (together with the Betacool_UNIX executable!) with
```
make clean
```

### For compilation on MacOSX Catalina, I used the following procedure:
1. Install MacPorts
2. Be sure everything is up do date
```
sudo port selfupdate
sudo port upgrade outdated
```
3. Add necessary packages
```
sudo port install gcc9
sudo port install libomp
sudo port install dos2unix
```
4. Select the right compiler
```
sudo port select --list gcc         # to list available options
sudo port select --set gcc mp-gcc9  # to select the one just installed above
```
:::Warning
Note that you might need to open a new shell to make the change effective
:::
5. Simply execute
```
make
```

---
# Logbook of things I did to obtain the present version

In _betatrack_ branch I re-shuffled a bit the files around, and tried to build an adeguate makefile to build betacool (at least the core part)
The first steps to ordering the file structure was inspired from "#backup.bat" file by A. Smirnov. 
For reference, here the bash commands I used initially:

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
I had then to deviate, and start removing old log files trying to keep the essential.

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

## Other modifications to the source code
- in a few files renamed the variable "gamma" to "gammaRel" to avoid re-definition of standard funciton in <math.h> which gives error compiling for MacOSX
- introduced a new "warning.h" and "warning.cpp" files containing "Warning()" function implementation which was duplicated in other source files
- adjusted the filename case of a few *.cpp files to be consistent with relative *.h files  

### Change case of files
```
for j in *.h; do b=$(echo $j | awk '{print tolower($0)}'); for i in *; do sed -i .bkp "s/"$b"/"$j"/g" $i; done; rm *.bkp; done;
```

### Make Git repository case sensitive
```
git config core.ignorecase false
```

## Check consistency of results from BETACOOL
For the time being I did a simple test like the following one using the BETACOOL example output run with Betacool_MacOSX
```
for i in *.cur; do echo $i; diff $i ../Examples/BETACOOL-ref/$i; done
```

### Convert DOS-like output files to UNIX-like standards
```
dos2unix ./Examples/BETACOOL-ref/*
```
