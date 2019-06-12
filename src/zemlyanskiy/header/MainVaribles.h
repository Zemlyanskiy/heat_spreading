#pragma once
#include "IOStream.h"
#include "MpiWorking.h"
//AD calculate data
double step;
double quad_step;
double sum;
double posX;
double posY;
double posZ;
int numbers_of_out;
int out_count;
int max_count;
int count;
double CurTime = 0;
double start, end;
double Restime;
double delta;
int num_threads;
char *file;
int prog_mode;
//Variables for CSR
int *first_line;
int *columns;
double *csr_values;
double *difference;
double *CsrDiag;
bool flag;
int result;
//Reading data
struct input_data data;
//Borders
int border_size = 2;
//Arrays for calculations
struct values3d points;      // array for caclulation results
struct values3d prev;        // array for previos state
struct values3d print_array; // array for previos state
//Constants
double Xstep;
double Ystep;
double Zstep;
double Xdivider;
double Ydivider;
double Zdivider;

static inline void InitializeMainVar(struct input_data data)
{
    //Calculate data
    numbers_of_out = (int)(data.T - data.t) / data.deltaOut;
    max_count = (int)(data.T - data.t) / data.deltaT;
    out_count = max_count / numbers_of_out;
    count = 0;
    InitValues3d(&points, Xlines, Ylines, Zlines);
    InitValues3d(&prev, Xlines, Ylines, Zlines);
    InitValues3d(&print_array, data.array.Lx, data.array.Ly, data.array.Lz);
    Xstep = (data.XMax - data.XMin) / data.array.Lx;
    Ystep = (data.YMax - data.YMin) / data.array.Ly;
    Zstep = (data.ZMax - data.ZMin) / data.array.Lz;
    Xdivider = 1 / (Xstep * Xstep);
    Ydivider = 1 / (Ystep * Ystep);
    Zdivider = 1 / (Zstep * Zstep);
}
