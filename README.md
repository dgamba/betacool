# Un-official version of BETACOOL code.

The present version is a copy of BETACOOL source code by A. Smirnov with corrections by N. Mityanina, A. Sidorin.

-- D. Gamba - Oct. 2019

In betatrack branch I re-shuffled a bit the files around, and tried to build an adeguate makefile to build betacool (at least the core part) ideally on 
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

