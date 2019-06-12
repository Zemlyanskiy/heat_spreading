#pragma once
#include <mpi.h>

#define Get(arr, x, y, z, Lx, Ly, Lz) (*((arr) + (z) + (y) * (Lz) + (x) * (Ly) * (Lz)))

unsigned process_per_axis;
unsigned Xlines;
unsigned Ylines;
unsigned Zlines;
unsigned Xcube_pos;
unsigned Ycube_pos;
unsigned Zcube_pos;

unsigned elements_per_process;

unsigned Xstart_copy_line;
unsigned Ystart_copy_line;
unsigned Zstart_copy_line;

int *x_return_lines, *y_return_lines, *z_return_lines, *x_return_start, *y_return_start, *z_return_start;
double **buffer;

MPI_Status status;

static inline void Initialize_MPI(struct input_data data, int rank, int proc_num, int border_size)
{
    process_per_axis = (unsigned)cbrt(proc_num);
    Xlines = data.array.Lx / process_per_axis;
    Ylines = data.array.Ly / process_per_axis;
    Zlines = data.array.Lz / process_per_axis;
    Xcube_pos = rank % process_per_axis;
    Ycube_pos = rank / process_per_axis;
    Zcube_pos = rank / process_per_axis * process_per_axis;
    // Add remainder ---
    if (Xcube_pos == process_per_axis - 1)
        Xlines += data.array.Lx % process_per_axis;
    if (Ycube_pos == process_per_axis - 1)
        Ylines += data.array.Ly % process_per_axis;
    if (Zcube_pos == process_per_axis - 1)
        Zlines += data.array.Lz % process_per_axis;

    Xstart_copy_line = Xcube_pos * data.array.Lx / process_per_axis;
    Ystart_copy_line = Ycube_pos * data.array.Ly / process_per_axis;
    Zstart_copy_line = Zcube_pos * data.array.Lz / process_per_axis;

    // Add borders ---
    if (Xcube_pos > 0)
    {
        Xlines += border_size;
        Xstart_copy_line -= border_size;
    }
    if (Xcube_pos < process_per_axis - 1)
        Xlines += border_size;
    if (Ycube_pos > 0)
    {
        Ylines += border_size;
        Ystart_copy_line -= border_size;
    }
    if (Ycube_pos < process_per_axis - 1)
        Ylines += border_size;
    if (Zcube_pos > 0)
    {
        Zlines += border_size;
        Zstart_copy_line -= border_size;
    }
    if (Zcube_pos < process_per_axis - 1)
        Zlines += border_size;

    // MPI communication variables ---
    elements_per_process = Xlines * Ylines * Zlines;

    // Values for return collection (used only on root rank)
    x_return_lines = (int *)malloc(proc_num * sizeof(int));
    y_return_lines = (int *)malloc(proc_num * sizeof(int));
    z_return_lines = (int *)malloc(proc_num * sizeof(int));

    x_return_start = (int *)malloc(proc_num * sizeof(int));
    y_return_start = (int *)malloc(proc_num * sizeof(int));
    z_return_start = (int *)malloc(proc_num * sizeof(int));

    buffer = (double **)malloc(proc_num * sizeof(double *));

    // Create on 1st process array of values for collect output data ---
    if (rank == 0)
    {
        x_return_lines[0] = Xlines;
        y_return_lines[0] = Ylines;
        z_return_lines[0] = Zlines;
        x_return_start[0] = Xstart_copy_line;
        y_return_start[0] = Ystart_copy_line;
        z_return_start[0] = Zstart_copy_line;

        buffer[0] = (double *)calloc(x_return_lines[0] * y_return_lines[0] * z_return_lines[0], sizeof(double));

        for (int i = 1; i < proc_num; i++)
        {
            MPI_Recv(x_return_lines + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(y_return_lines + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(z_return_lines + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(x_return_start + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(y_return_start + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(z_return_start + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            buffer[i] = (double *)calloc(x_return_lines[i] * y_return_lines[i] * z_return_lines[i], sizeof(double));
        }
    }
    else
    {
        MPI_Send(&Xlines, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Ylines, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Zlines, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Xstart_copy_line, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Ystart_copy_line, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Zstart_copy_line, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
}

// MPI Functions
static inline void UpdateBorders(double *points, unsigned Lx, unsigned Ly, unsigned Lz, unsigned rank,
                                 unsigned process_per_axis, unsigned Xcube_pos, unsigned Ycube_pos, unsigned Zcube_pos, unsigned proc_num, unsigned border_size)
{
    if (proc_num > 1)
    {
        //int section_size = Lx * Ly;
        //double *helpingbuf = calloc(Ly, sizeof(double));
        // Send X
        // if (Xcube_pos > 0)
        // {
        //     // send X border section to previos process
        //     for (unsigned n = 0; n < border_size; n++)
        //         MPI_Send(points + section_size * (border_size + n), section_size, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
        //     printf("Send X > 0\n");
        // }
//        if (Xcube_pos < process_per_axis - 1)
//        {
//            // send X border section to next process
//            for (unsigned n = 0; n < border_size; n++)
//                MPI_Send(points + section_size * (Lz - border_size + n), section_size, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
//            printf("Send X < process_per_axis - 1\n");
//        }
//
//        // Recieve X
//        if (Xcube_pos > 0)
//        {
//            // get X border section from previos process
//            for (unsigned n = 0; n < border_size; n++)
//                MPI_Recv(points + section_size * n, section_size, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &Status);
//            printf("Recieve X > 0\n");
//        }
        // if (Xcube_pos < process_per_axis - 1)
        // {
        //     // get X border section from next process
        //     for (unsigned n = 0; n < border_size; n++)
        //         MPI_Recv(points + section_size * (Lz + n), section_size, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &Status);
        //     printf("Recieve X < process_per_axis - 1\n");
        // }

        // Send Y
        /*if (Ycube_pos > 0)
        {
            printf("Send Y > 0\n");
            // send Y border section to previos process
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    MPI_Send(points + section_size * i + Lz * (border_size + n), Ly, MPI_DOUBLE, rank - 2, rank, MPI_COMM_WORLD);
        }
        if (Ycube_pos < process_per_axis - 1)
        {
            printf("Send Y < process_per_axis - 1\n");
            // send Y border section to next
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    MPI_Send(points + (i + 1) * section_size - Lz * (border_size * 2 - n), Ly, MPI_DOUBLE, rank + 2, rank, MPI_COMM_WORLD);
        }

        // Recieve Y
        if (Ycube_pos > 0)
        {
            printf("Recieve Y > 0\n");
            // get Y border section from previos process
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    MPI_Recv(points + section_size * i + Lz * n, Ly, MPI_DOUBLE, rank - 2, rank - 2, MPI_COMM_WORLD, &Status);
        }

        if (Ycube_pos < process_per_axis - 1)
        {
            printf("Recieve Y < process_per_axis - 1\n");
            // get Y border section from next process
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    MPI_Recv(points + (i + 1) * section_size - Lz * (border_size - n), Ly, MPI_DOUBLE, rank + 2, rank + 2, MPI_COMM_WORLD, &Status);
        }

        // Send Z
        if (Zcube_pos > 0)
        {
            printf("Send Z > 0\n");
            // send Z border section to previos process
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    for (unsigned j = 0; j < Ly; i++)
                        MPI_Send(points + section_size * i + Lz * j + (border_size + n), j, MPI_DOUBLE, rank - 3, rank, MPI_COMM_WORLD);
        }
        if (Zcube_pos < process_per_axis - 1)
        {
            printf("Send Z < process_per_axis - 1\n");
            // send Z border section to next
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    for (unsigned j = 0; j < Ly; i++)
                        MPI_Send(points + (i + 1) * section_size + Lz * j - (border_size * 2 - n), j, MPI_DOUBLE, rank + 3, rank, MPI_COMM_WORLD);
        }

        // Recieve Z
        if (Zcube_pos > 0)
        {
            printf("Recieve Z > 0\n");
            // get Z border section from previos process
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    for (unsigned j = 0; j < Ly; j++)
                        MPI_Recv(points + section_size * i + Lz * j + n, j, MPI_DOUBLE, rank - 3, rank - 3, MPI_COMM_WORLD, &Status);
        }

        if (Zcube_pos < process_per_axis - 1)
        {
            printf("Recieve Z < process_per_axis - 1\n");
            // get Z border section from next process
            for (unsigned n = 0; n < border_size; n++)
                for (unsigned i = 0; i < Lx; i++)
                    for (unsigned j = 0; j < Ly; j++)
                        MPI_Recv(points + (i + 1) * section_size + Lz * j - (border_size - n), j, MPI_DOUBLE, rank + 3, rank + 3, MPI_COMM_WORLD, &Status);
        }*/
    }
}

void ReceivePrint(double *print_array, double *points, int Xlines, int Ylines, int Zlines,
                              int *x_return_lines, int *y_return_lines, int *z_return_lines,
                              int *x_return_start, int *y_return_start, int *z_return_start,
                              double **buffer,
                              char *file, struct input_data data, int proc_num)
{

    //Can't send-receive from 0 process to himself
    for (unsigned i = 0; i < Xlines; i++)
        for (unsigned j = 0; j < Ylines; j++)
            for (unsigned k = 0; k < Zlines; k++)
            {
                Get(print_array, i, j, k, data.array.Lx, data.array.Ly, data.array.Lz) = Get(points, i, j, k, Xlines, Ylines, Zlines);
            }

    for (unsigned sender_number = 1; sender_number < proc_num; sender_number++)
    {
        MPI_Recv(buffer[sender_number], x_return_lines[sender_number] * y_return_lines[sender_number] * z_return_lines[sender_number],
                 MPI_DOUBLE, sender_number, sender_number, MPI_COMM_WORLD, &status);

        for (int i = 0; i < x_return_lines[sender_number]; i++)
            for (int j = 0; j < y_return_lines[sender_number]; j++)
                for (int k = 0; k < z_return_lines[sender_number]; k++)
                {
                    Get(print_array, x_return_start[sender_number] + i, y_return_start[sender_number] + j, k, data.array.Lx, data.array.Ly, data.array.Lz) =
                        Get(buffer[sender_number], i, j, k, x_return_lines[sender_number], y_return_lines[sender_number], z_return_lines[sender_number]);
                }
    }
}

void SendPrint(double *points, int Xlines, int Ylines, int Zlines, int rank)
{
    MPI_Send(points, Xlines * Ylines * Zlines, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
}
