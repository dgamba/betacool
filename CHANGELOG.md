Version 7.2.beta
===

- Changed variable name `beta` into `betaRel` in `bpIBS.cpp` to avoid compilation errors
- Added some additional Warnings here and there for better debugging
- Changed variable name `InRangeData` into `StartTimestamp`. Before, this variable was computed at compilation time as unix timestamp (more or less) in seconds, and used to check that a simulation was not taking too long... This made an executable not-operational some time after compilation time. The variable now carries the execution time stamp in seconds (more or less). To be checked if this is more stable.
