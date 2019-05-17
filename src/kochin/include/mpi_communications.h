#ifndef __MPI_COMMUNICATIONS
#define __MPI_COMMUNICATIONS

#include <mpi.h>

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

#endif /* __MPI_COMMUNICATIONS */