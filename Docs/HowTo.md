Simple documentation of BetaCOOL
===

# How to Execute

BetaCOOL uses a .bld configuration file to store all configuration/parameters to describe a situation.
To execute a task, simply run:
``` 
betacool.exe inputfile.bld /parameters
```
where `/parameters` is in reality can just be `/` followed by one letter only (the following ones are ignored!) - see Betacool.cpp:
- `d`: iDynamics.Dynamics();
- `m`: iDynamics.ModelBeam();
- `t`: iDynamics.Tracking();
- `l`: iDynamics.Lattice();
- `r`: iDynamics.Rates();
- `b`: iDynamics.BeamTest();
- `f`: It computes the friction force: iDynamics.FFTest();
- `g`: probably the same as `f`: iForce.FF(iTime, iEbeam, iRing);
- `c`: iDynamics.LuminosityCalculation();
- `s`: iDraw.SpaceChargeDraw(iBeam, iEbeam);
- `9` == `i`: several tests/computations ...
- ... see source for other options ...

# Source Code notes
Some useful notes
- BData: seems to be a generic class to host data (vectors)
    - Used to get and store data from/to files
- doubleU: is a double that carries a unit
    - U_xxxx: are units. They can be multiplied with other doubleU 
    - M_6, m_3, u_6, ....: power of 10: Mega, milli, micro, ...

# Physics notes

## Friction Force

### Models 

Several models are available in BETACOOL.
- All is handled by xForce class (see xForce.h) 
- At the end, most of the code is going to look into v[3]; and f[3] to look for velocity and force componets (x, y, z-longitudinal).
    - for transverse force, Vtr and Ftr variables are also used (for plotting, for example)

Models (with specifc parameters being used):
- 0: Budker - Budker
- 1: NonMag - Non-magnetized
    - obj.BLDContent[61][6,7,8]
- 2: DerSkr - Derbenev-Skrinsky-Meshkov
    - obj.BLDContent[61][0] = 2   # Smoothing coefficient (Derbenev model)
- 3: Parhom (default) - Parkhomchuk
    - obj.BLDContent[61][1] = 0   # Disable/Enable Effective of Angular Spread (if == 1, Teff computed from angular spread instead of from parameter 4 (which is then re-computed))
    - obj.BLDContent[61][2] = 0   # Disable/Enable Effective Temperature - NOT USED!!!
    - obj.BLDContent[61][3] = 0.0001   # Angular spread, rad - used only if BLDContent[61][1] == 1
    - obj.BLDContent[61][4] = 0.2 # Effective Temperature - used only if BLDContent[61][1] == 0
    - Additionally:
        - F.n_e
        - F.Z
        - F.V_eff_e  = (F.TempEff/U_me)^0.5;
            - F.TempEff   = Data.OnGet(61,5,1);
        - F.V_long_e = ((F.Ltemp/U_me)^0.5);
            - depends on e- model chosen Data.OnGet(55, 1, )
- 4: Toepffer - Toepffer
    - obj.BLDContent[61][9,10,11]
- 5: Table:
    - used to get data from data file, they should be two files, one with longitudinal and one with transverse forces
- 6: D3 - 3D model:
    - Integration over real coordinates of electrons 
- (7): D4 - 3D analytical model:
    - Integration with Maxwellian distribution ?!


