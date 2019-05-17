/*
 * Common methods and includes for all methods
 */

#ifndef __COMMON
#define __COMMON

#include <math.h>
#include <omp.h>
#include <windows.h>

#include "io_parser.h"
#include "mpi_communications.h"

// Array macro
#define Get(arr, x, y, z, Lx, Ly, Lz) (*((arr) + (z) + (y) * (Lz) + (x) * (Ly) * (Lz)))

#endif /* __COMMON */