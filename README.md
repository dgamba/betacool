# Un-official version of BETACOOL code.

The present version is a copy of BETACOOL source code by A. Smirnov with corrections by N. Mityanina, A. Sidorin.

-- D. Gamba - Oct. 2019

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
- Dat:               folder with data files used to model previous coolers
- Bat:               folder with original *.bat files from authors
- Build:             normally and empty directory, it gets filled by temporary/log files when trying to compile BETACOOL for different architectures

---
# How to compile

## For compilation on Windows
Just load Betacool.vcxproj with Microsoft Visual Studio 2013 and "Build" it for your Windows architecture. 
This will create a betacool_<arc>.exe file 


## For compilation on UNIX
TODO, but will try to make it as simple as:
`
make
`
and this will create betacool_<arc> file executable on your architecture  


## For compilation on MacOSX Catalina:
1. Install MacPorts
2. Be sure everything is up do date
`
sudo port selfupdate
sudo port upgrade outdated
`
3. Add necessary packages
`
sudo port install gcc9
sudo port install libomp
sudo port install dos2unix
`
4. Select the right compiler
`
sudo port select --list gcc         # to list available options
sudo port select --set gcc mp-gcc9  # to select the one just installed above
`
:::Warning
Note that you might need to open a new shell to make the change effective
:::
5. Simply execute
`
make
`
and this will compile the source code under the Build folder, and create the Betacool_MacOSX executable in the main folder.


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

## Modifications to the source code
- in a few files renamed the variable "gamma" to "gammaRel" to avoid re-definition of standard funciton in <math.h> which gives error compiling for MacOSX
- introduced a new "warning.h" and "warning.cpp" files containing "Warning()" function implementation which was duplicated in other source files
- adjusted case of a few *.cpp files to be consistent with relative *.h files  

### Check consistency of results from BETACOOL
For the time being I did simple test like the following one using the BETACOOL example output run with Betacool_MacOSX
`
dos2unix ../Examples/BETACOOL-ref/*
for i in *.cur; do echo $i; diff $i ../Examples/BETACOOL-ref/$i; done
`
### Change case of files
`
for j in *.h; do b=$(echo $j | awk '{print tolower($0)}'); for i in *; do sed -i .bkp "s/"$b"/"$j"/g" $i; done; rm *.bkp; done;
`

### Make Git repository case sensitive
`
git config core.ignorecase false
`
