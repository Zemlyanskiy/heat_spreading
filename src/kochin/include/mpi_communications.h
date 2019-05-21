#ifndef __MPI_COMMUNICATIONS
#define __MPI_COMMUNICATIONS

#include <mpi.h>
#include "common.h"


MPI_Status status;

// MPI Functions
inline void SendBorderValues(double* points, unsigned Lx, unsigned Ly, unsigned Lz, unsigned rank,
    unsigned process_per_axis, unsigned Xcube_pos, unsigned Ycube_pos, unsigned proc_num
) {
    if (proc_num > 1) {
        MPI_Status Status;
        int section_size = Ly * Lz;
        // Send
        if (Xcube_pos > 0) {
            // send X border section to previos process
            MPI_Send(points + section_size, section_size, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
        }
        if (Xcube_pos < process_per_axis - 1) {
            // send X border section to next process
            MPI_Send(points + (Lx - 2)*section_size, section_size, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
        }

        if (Ycube_pos > 0) {
            // send Y border section to previos process
            for (unsigned i = 0; i < Lx; i++)
                MPI_Send(points + Lz + section_size * i, Lz, MPI_DOUBLE, rank - process_per_axis, rank, MPI_COMM_WORLD);
        }
        if (Ycube_pos < process_per_axis - 1) {
            // send Y border section to next process
            for (unsigned i = 0; i < Lx; i++)
                MPI_Send(points + (i + 1)*section_size - Lz * 2, Lz, MPI_DOUBLE, rank + process_per_axis, rank, MPI_COMM_WORLD);
        }

        // Recieve
        if (Xcube_pos > 0) {
            // get X border section from previos process
            MPI_Recv(points, section_size, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &Status);
        }
        if (Xcube_pos < process_per_axis - 1) {
            // get X border section from next process
            MPI_Recv(points + (Lx - 1)*section_size,
                section_size, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &Status);
        }

        if (Ycube_pos > 0) {
            // get Y border section from previos process
            for (unsigned i = 0; i < Lx; i++)
                MPI_Recv(points + section_size * i, Lz,
                    MPI_DOUBLE, rank - process_per_axis, rank - process_per_axis, MPI_COMM_WORLD, &Status);
        }
        if (Ycube_pos < process_per_axis - 1) {
            // get Y border section from next process
            for (unsigned i = 0; i < Lx; i++)
                MPI_Recv(points + (i + 1)*section_size - Lz,
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


int *x_return_lines, *y_return_lines, *x_return_start, *y_return_start;
double** buffer;

void initialize_return_variables( unsigned Xlines, unsigned Ylines, unsigned Zlines,
                                  unsigned Xstart_copy_line, unsigned Ystart_copy_line,
                                  int rank, int proc_num ){
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
        buffer[0] = (double*)calloc(x_return_lines[0] * y_return_lines[0] * Zlines, sizeof(double));

        for (int i = 1; i < proc_num; i++) {
            MPI_Recv(x_return_lines + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(y_return_lines + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(x_return_start + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(y_return_start + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            buffer[i] = (double*)calloc(x_return_lines[i] * y_return_lines[i] * Zlines, sizeof(double));
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

#endif /* __MPI_COMMUNICATIONS */