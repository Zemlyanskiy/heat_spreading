#include "common.h"

int main(int argc, char* argv[])
{
    // MPI initialization ---
    unsigned rank, proc_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if ((sqrt(proc_num) - (unsigned)sqrt(proc_num)) != 0) {
        printf("The number of processes must be able to proccess sqrt() without remainder.\n");
        MPI_Finalize();
        return 1;
    }

    // Processing command line arguments ---
    char* inputFile = NULL;
    unsigned number_of_threads = 1;
    {
        if (argc > 1)
            inputFile = argv[1];
        if (argc > 2)
            number_of_threads = strtol(argv[2], NULL, 0);
        if (argc > 3) {
            printf("Incorrect number of arguments");
            return 4;
        }
    }

    // Read input data ---
    struct input data;
    data = ReadInput(inputFile);

    // Variables initialization ---
    if ( data.t >= data.T) {
        if (rank == 0) { printf("Start time >= end time, please clear io file\n"); }
        MPI_Finalize();
        return 0;
    }

    int border_size = 4;
    initialize_common_variables(data, rank, proc_num, border_size);

    // Standart computation arrays ---
    double *prev, *current, *next;
    prev = (double*)calloc(Xlines * Ylines * data.Lz, sizeof(double));
    current = (double*)calloc(Xlines * Ylines * data.Lz, sizeof(double));
    next = (double*)calloc(Xlines * Ylines * data.Lz, sizeof(double));

    double *k1, *k2, *k3, *k4, *medium, *derivative;
    k1 = (double*) calloc( Xlines* Ylines * data.Lz, sizeof(double));
    k2 = (double*) calloc( Xlines* Ylines * data.Lz, sizeof(double));
    k3 = (double*) calloc( Xlines* Ylines * data.Lz, sizeof(double));
    k4 = (double*) calloc( Xlines* Ylines * data.Lz, sizeof(double));
    medium = (double*) calloc(Xlines* Ylines* data.Lz, sizeof(double));
    derivative = (double*) calloc (Xlines * Ylines * data.Lz, sizeof(double));

    // Array for print result to file ---
    double* print_array = NULL;
    if (rank == 0) {
        print_array = (double*)calloc(data.Lx * data.Ly * data.Lz, sizeof(double));
    }

    // Env preparation for OpenMP ---
    omp_set_num_threads(number_of_threads);
    double time = 0;
    double out_time;

    // Read start state of array ---
    for (unsigned i = 0; i < Xlines; i++)
        for (unsigned j = 0; j < Ylines; j++)
            for (unsigned k = 0; k < data.Lz; k++) {
                Get(current, i, j, k, Xlines, Ylines, data.Lz) = Get(data.arr[0], Xstart_copy_line + i, Ystart_copy_line + j, k, data.Lx, data.Ly, data.Lz);
                Get(next, i, j, k, Xlines, Ylines, data.Lz) = Get(data.arr[1], Xstart_copy_line + i, Ystart_copy_line + j, k, data.Lx, data.Ly, data.Lz);
            }

    // Calculation process ---
    unsigned i, j, k, count;
    for (count = 0; count <= count_threshold; count++) {
        if (rank==0) out_time = omp_get_wtime();

#pragma omp parallel for
        for (i = 0; i < elements_per_process; i++) {
            prev[i] = current[i];
            current[i] = next[i];
        }
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(next, i, j, k, Xlines, Ylines, data.Lz) =
                        2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) - Get(prev, i, j, k, Xlines, Ylines, data.Lz) + pow(data.Sigma,2) * pow(data.deltaT,2)*(
                        (Get(current, i - 1, j, k, Xlines, Ylines, data.Lz) - 2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) + Get(current, i + 1, j, k, Xlines, Ylines, data.Lz)) * Xdivider +
                        (Get(current, i, j - 1, k, Xlines, Ylines, data.Lz) - 2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) + Get(current, i, j + 1, k, Xlines, Ylines, data.Lz)) * Ydivider +
                        (Get(current, i, j, k - 1, Xlines, Ylines, data.Lz) - 2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) + Get(current, i, j, k + 1, Xlines, Ylines, data.Lz)) * Zdivider);
                    // next[i][j][k] = 2*current[i][j][k] - prev[i,j,k] + data.Sigma^2*data.deltaT^2*((current[i-1][j][k] - 2*current[i][j][k] + current[i+1][j][k]) * Xdivider +
                    //                                                                                (current[i][j-1][k] - 2*current[i][j][k] + current[i][j+1][k]) * Ydivider +
                    //                                                                                (current[i][j][k-1] - 2*current[i][j][k] + current[i][j][k+1]) * Zdivider);
                }

#if RUNGE_KUTTA
#pragma omp parallel for
        for (i = 0; i < elements_per_process; i++)
            prev[i] = points[i];

        // k1 calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(k1, i, j, k, Xlines, Ylines, data.Lz) = pow(data.Sigma,2) * pow(data.deltaT,2)*(
                        (Get(current, i - 1, j, k, Xlines, Ylines, data.Lz) - 2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) + Get(current, i + 1, j, k, Xlines, Ylines, data.Lz)) * Xdivider +
                        (Get(current, i, j - 1, k, Xlines, Ylines, data.Lz) - 2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) + Get(current, i, j + 1, k, Xlines, Ylines, data.Lz)) * Ydivider +
                        (Get(current, i, j, k - 1, Xlines, Ylines, data.Lz) - 2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) + Get(current, i, j, k + 1, Xlines, Ylines, data.Lz)) * Zdivider);
                }

#pragma omp parallel for
        for (i = 0; i < Xlines; i++)
            for (j = 0; j < Ylines; j++)
                for (k = 0; k < data.Lz; k++) {
                    Get(medium, i, j, k, Xlines, Ylines, data.Lz) =  2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) - Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(k1, i, j, k, Xlines, Ylines, data.Lz) * 0.5;
                }

        // k2 calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(k2, i, j, k, Xlines, Ylines, data.Lz) = pow(data.Sigma,2) * pow(data.deltaT,2)*(
                        (Get(medium, i - 1, j, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i + 1, j, k, Xlines, Ylines, data.Lz)) * Xdivider +
                        (Get(medium, i, j - 1, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j + 1, k, Xlines, Ylines, data.Lz)) * Ydivider +
                        (Get(medium, i, j, k - 1, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j, k + 1, Xlines, Ylines, data.Lz)) * Zdivider);
                }

#pragma omp parallel for
        for (i = 0; i < Xlines; i++)
            for (j = 0; j < Ylines; j++)
                for (k = 0; k < data.Lz; k++) {
                    Get(medium, i, j, k, Xlines, Ylines, data.Lz) = 2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) - Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(k2, i, j, k, Xlines, Ylines, data.Lz) * 0.5;
                }

        // k3 calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(k3, i, j, k, Xlines, Ylines, data.Lz) = pow(data.Sigma,2) * pow(data.deltaT,2)*(
                        (Get(medium, i - 1, j, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i + 1, j, k, Xlines, Ylines, data.Lz)) * Xdivider +
                        (Get(medium, i, j - 1, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j + 1, k, Xlines, Ylines, data.Lz)) * Ydivider +
                        (Get(medium, i, j, k - 1, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j, k + 1, Xlines, Ylines, data.Lz)) * Zdivider);
                }

#pragma omp parallel for
        for (i = 0; i < Xlines; i++)
            for (j = 0; j < Ylines; j++)
                for (k = 0; k < data.Lz; k++) {
                    Get(medium, i, j, k, Xlines, Ylines, data.Lz) = 2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) - Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(k3, i, j, k, Xlines, Ylines, data.Lz);
                }

        // k4 calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(k4, i, j, k, Xlines, Ylines, data.Lz) = pow(data.Sigma,2) * pow(data.deltaT,2)*(
                        (Get(medium, i - 1, j, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i + 1, j, k, Xlines, Ylines, data.Lz)) * Xdivider +
                        (Get(medium, i, j - 1, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j + 1, k, Xlines, Ylines, data.Lz)) * Ydivider +
                        (Get(medium, i, j, k - 1, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j, k + 1, Xlines, Ylines, data.Lz)) * Zdivider);
                }

        // derivative calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(derivative, i, j, k, Xlines, Ylines, data.Lz) = (Get(k1, i, j, k, Xlines, Ylines, data.Lz) +
                                                                         Get(k2, i, j, k, Xlines, Ylines, data.Lz) * 2 +
                                                                         Get(k3, i, j, k, Xlines, Ylines, data.Lz) * 2 +
                                                                         Get(k4, i, j, k, Xlines, Ylines, data.Lz)) * 0.1666666666;
                }

        // final result ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(points, i, j, k, Xlines, Ylines, data.Lz) = 2 * Get(current, i, j, k, Xlines, Ylines, data.Lz) - Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(derivative, i, j, k, Xlines, Ylines, data.Lz);
                }
#endif /*RUNGE_KUTTA*/

        SendBorderValues(next, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num, border_size);

        // Delete output time from time calculating ---
        if (rank==0) time += omp_get_wtime() - out_time;

        // Print ---
        if (count % out_freq == 0 && count) {
            if (rank == 0) {
                receive_and_print_master( print_array, next, Xlines, Ylines,
                                          x_return_lines, y_return_lines, x_return_start, y_return_start,
                                          buffer, inputFile, data, proc_num);
            }
            else {
                send_to_print(next, Xlines, Ylines, data.Lz, rank);
            }
        }
    }

    // End ---
    if (!rank) { printf("Time is: %lf\n", time); }
    MPI_Finalize();
    return 0;
}