## heat_spreading

### Structure
 - data - data sheets for process computation
 - include - common headers of implementation
 - legacy - deprecated realizations
 - out - folder for builds binaries output
 - src - `c` files that realized coputations methods
 - build.bat - script for build and run computation

### Next steps
 - Realize wave equation process (Implicit Euler)
 - Remove brute force calculation in Implicit Euler
 - Fix heat_map script

### Q:
1. Matrix methods
    - How we would implement matrix multiplication with 3 times?
### A:
1. Add prev time in separate state.

### Build commands

    .\build.bat euler input1.txt 4
    python ..\..\heat_map.py C:\work\Git\other\heat_spreading\src\kochin\data\input1.txt