/*
 * Common methods and includes for all methods
 */

#ifndef __COMMON
#define __COMMON

#include <math.h>
#include <omp.h>
#include <mpi.h>
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

unsigned process_per_axis;
unsigned Xlines;
unsigned Ylines;
unsigned Xcube_pos;
unsigned Ycube_pos;

unsigned elements_per_process;

unsigned Xstart_copy_line;
unsigned Ystart_copy_line;

int *x_return_lines, *y_return_lines, *x_return_start, *y_return_start;
double** buffer;

MPI_Status status;

void initialize_common_variables( struct input data, int rank, int proc_num, int border_size ){
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

    process_per_axis = (unsigned)sqrt(proc_num);
    Xlines = data.Lx / process_per_axis;
    Ylines = data.Ly / process_per_axis;
    Xcube_pos = rank % process_per_axis;
    Ycube_pos = rank / process_per_axis;

    // Add remainder ---
    if (Xcube_pos == process_per_axis - 1)
        Xlines += data.Lx % process_per_axis;
    if (Ycube_pos == process_per_axis - 1)
        Ylines += data.Ly % process_per_axis;

    Xstart_copy_line = Xcube_pos * data.Lx / process_per_axis;
    Ystart_copy_line = Ycube_pos * data.Ly / process_per_axis;

    // Add borders ---
    if (Xcube_pos > 0) {
        Xlines+=border_size;
        Xstart_copy_line-=border_size;
    }
    if (Xcube_pos < process_per_axis - 1)
        Xlines+=border_size;
    if (Ycube_pos > 0) {
        Ylines+=border_size;
        Ystart_copy_line-=border_size;
    }
    if (Ycube_pos < process_per_axis - 1)
        Ylines+=border_size;

    // MPI communication variables ---
    elements_per_process = Xlines * Ylines * data.Lz;

    // Values for return collection (used only on root rank)
    x_return_lines = (int*)malloc(proc_num * sizeof(int));
    y_return_lines = (int*)malloc(proc_num * sizeof(int));
    x_return_start = (int*)malloc(proc_num * sizeof(int));
    y_return_start = (int*)malloc(proc_num * sizeof(int));
    buffer = (double**)malloc(proc_num * sizeof(double*));

    // Create on 1st process array of values for collect output data ---
    if (rank == 0) {
        x_return_lines[0] = Xlines;
        y_return_lines[0] = Ylines;
        x_return_start[0] = Xstart_copy_line;
        y_return_start[0] = Ystart_copy_line;
        buffer[0] = (double*)calloc(x_return_lines[0] * y_return_lines[0] * data.Lz, sizeof(double));

        for (int i = 1; i < proc_num; i++) {
            MPI_Recv(x_return_lines + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(y_return_lines + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(x_return_start + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(y_return_start + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            buffer[i] = (double*)calloc(x_return_lines[i] * y_return_lines[i] * data.Lz, sizeof(double));
        }

    }
    else
    {
        MPI_Send(&Xlines, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Ylines, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Xstart_copy_line, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Ystart_copy_line, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
}

// MPI Functions
inline void SendBorderValues(double* points, unsigned Lx, unsigned Ly, unsigned Lz, unsigned rank,
    unsigned process_per_axis, unsigned Xcube_pos, unsigned Ycube_pos, unsigned proc_num, unsigned border_size
) {
    if (proc_num > 1) {
        MPI_Status Status;
        int section_size = Ly * Lz;
        // Send X
        if (Xcube_pos > 0) {
            // send X border section to previos process
            for (unsigned n = 0; n < border_size; n++)
                MPI_Send(points + section_size * (border_size + n), section_size, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
        }
        if (Xcube_pos < process_per_axis - 1) {
            // send X border section to next process
            for (unsigned n = 0; n < border_size; n++)
                MPI_Send(points + (Lx - border_size*2 + n)*section_size, section_size, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
        }

        // Recieve X
        if (Xcube_pos > 0) {
            // get X border section from previos process
            for (unsigned n = 0; n < border_size; n++)
                MPI_Recv(points + section_size * n, section_size, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &Status);
        }
        if (Xcube_pos < process_per_axis - 1) {
            // get X border section from next process
            for (unsigned n = 0; n < border_size; n++)
                MPI_Recv(points + (Lx - border_size + n)*section_size,
                    section_size, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &Status);
        }

        // Send Y
        if (Ycube_pos > 0) {
            // send Y border section to previos process
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    MPI_Send(points + section_size * i + Lz * (border_size + n), Lz, MPI_DOUBLE, rank - process_per_axis, rank, MPI_COMM_WORLD);
        }
        if (Ycube_pos < process_per_axis - 1) {
            // send Y border section to next
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    MPI_Send(points + (i + 1)*section_size - Lz * (border_size*2 - n), Lz, MPI_DOUBLE, rank + process_per_axis, rank, MPI_COMM_WORLD);
        }

        // Recieve Y
        if (Ycube_pos > 0) {
            // get Y border section from previos process
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    MPI_Recv(points + section_size * i + Lz*n, Lz,
                        MPI_DOUBLE, rank - process_per_axis, rank - process_per_axis, MPI_COMM_WORLD, &Status);
        }

        if (Ycube_pos < process_per_axis - 1) {
            // get Y border section from next process
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    MPI_Recv(points + (i + 1)*section_size - Lz * (border_size - n),
                        Lz, MPI_DOUBLE, rank + process_per_axis, rank + process_per_axis, MPI_COMM_WORLD, &Status);
        }
    }
}

void receive_and_print_master( double* print_array, double* next, int Xlines, int Ylines,
                               int *x_return_lines, int* y_return_lines,
                               int* x_return_start, int* y_return_start,
                               double** buffer,
                               char* inputFile, struct input data, int proc_num) {

    //Can't send-receive from 0 process to himself
    for (unsigned i = 0; i < Xlines; i++)
        for (unsigned j = 0; j < Ylines; j++)
            for (unsigned k = 0; k < data.Lz; k++) {
                Get(print_array, i, j, k, data.Lx, data.Ly, data.Lz) = Get(next, i, j, k, Xlines, Ylines, data.Lz);
            }

    for (unsigned sender_number = 1; sender_number < proc_num; sender_number++) {
        MPI_Recv(buffer[sender_number], x_return_lines[sender_number] * y_return_lines[sender_number] * data.Lz,
            MPI_DOUBLE, sender_number, sender_number, MPI_COMM_WORLD, &status);

        for (int i = 0; i < x_return_lines[sender_number]; i++)
            for (int j = 0; j < y_return_lines[sender_number]; j++)
                for (int k = 0; k < data.Lz; k++) {
                    Get(print_array, x_return_start[sender_number] + i, y_return_start[sender_number] + j, k, data.Lx, data.Ly, data.Lz) =
                    Get(buffer[sender_number], i, j, k, x_return_lines[sender_number], y_return_lines[sender_number], data.Lz);
                }
    }

    OutputArr(inputFile, print_array, data.Lx*data.Ly*data.Lz);
    current_time += data.deltaOut;
    OutputCurrentTime(inputFile, current_time);
}

void send_to_print(double* next, int Xlines, int Ylines, int Zlines, int rank) {
    MPI_Send(next, Xlines*Ylines*Zlines, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
}

#endif /* __COMMON */