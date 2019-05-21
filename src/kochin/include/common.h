/*
 * Common methods and includes for all methods
 */

#ifndef __COMMON
#define __COMMON

#include <math.h>
#include <omp.h>
#include <windows.h>

#include "io_parser.h"

// Array macro
#define Get(arr, x, y, z, Lx, Ly, Lz) (*((arr) + (z) + (y) * (Lz) + (x) * (Ly) * (Lz)))


// Common computation variables

double Xstep;
double Ystep;
double Zstep;

double Xdivider;
double Ydivider;
double Zdivider;

double   current_time;
unsigned out_freq;
unsigned count_threshold;

void initialize_common_variables( struct input data ){
    Xstep = (data.Xmax - data.Xmin) / data.Lx;
    Ystep = (data.Ymax - data.Ymin) / data.Ly;
    Zstep = (data.Zmax - data.Zmin) / data.Lz;

    Xdivider = 1 / (Xstep*Xstep);
    Ydivider = 1 / (Ystep*Ystep);
    Zdivider = 1 / (Zstep*Zstep);

    current_time = data.t;
    out_freq = (unsigned)((data.T - data.t) / data.deltaOut);
    count_threshold = (unsigned)((data.T - data.t) / data.deltaT);
    out_freq = count_threshold / out_freq;
}

#endif /* __COMMON */