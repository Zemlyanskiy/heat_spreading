## heat_spreading

### Structure
 - data - data sheets for process computation
 - include - common headers of implementation
 - legacy - deprecated realizations
 - out - folder for builds binaries output
 - src - `c` files that realized coputations methods
 - build.bat - script for build and run computation

### Next steps
 - Realize wave equation process (Runge-Kutta, Implicit Euler, Matrix Multiplication method)
 - Realize 4 levels borders sending in Runge-Kutta
 - Optimize Matrix methods - delete 1 redundant operation
 - Remove brute force calculation in Implicit Euler

### Questions:
1. Initial conditions
    - How we must initialize borders?
    - Should the `prev` and `currrent` be the same at the beginnig?
2. Matrix methods
    - How we would implement matrix multiplication with 3 times?

### Build commands

    .\build.bat euler input1.txt 4
    python ..\..\heat_map.py C:\work\Git\other\heat_spreading\src\kochin\data\input1.txt