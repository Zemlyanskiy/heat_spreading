#include "common.h"
#include "csr_api.h"

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

    int border_size = 1;
    initialize_common_variables(data, rank, proc_num, border_size);

    // Standart computation arrays ---
    double *prev, *current, *next;
    prev = (double*)calloc(Xlines * Ylines * data.Lz, sizeof(double));
    current = (double*)calloc(Xlines * Ylines * data.Lz, sizeof(double));
    next = (double*)calloc(Xlines * Ylines * data.Lz, sizeof(double));

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

    // Yakobi method variables ---
    double * values;
    int* column_num;
    int* line_first;
    int flag, result;
    double delta, sum;
    double *MainDiag = NULL;
    double* discrepancy;
    discrepancy = (double*)calloc(elements_per_process,sizeof(double));
    MainDiag = (double*)calloc(elements_per_process, sizeof(double));
    delta = 0.00000000001;
    flag = 0;

    // Initialize matrix ---
    unsigned values_number = (Xlines - 2) * (Ylines - 2) * (data.Lz - 2) * 6 + Xlines*Ylines*data.Lz;
    values = (double*)calloc(values_number, sizeof(double));
    column_num = (unsigned*)calloc(values_number, sizeof(unsigned));
    line_first = (unsigned*)calloc(elements_per_process + 1, sizeof(unsigned));

    unsigned counter = 0;
    for (unsigned i = 0; i < Xlines; i++)
        for (unsigned j = 0; j < Ylines; j++)
            for (unsigned k = 0; k < data.Lz; k++) {
                line_first[i*Ylines*data.Lz + j * data.Lz + k] = counter;
                if (i == 0 || j == 0 || k == 0 || i == Xlines - 1 || j == Ylines - 1 || k == data.Lz - 1) {
                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k;
                    MainDiag[i*Ylines*data.Lz + j * data.Lz + k] = 1;
                    values[counter++] = 1;
                }
                else {
                    // x y z central*3 z y x
                    column_num[counter] = (i - 1) * Ylines*data.Lz + j * data.Lz + k;
                    values[counter++] = - pow(data.Sigma,2) * pow(data.deltaT,2)* Xdivider;
                    column_num[counter] = i * Ylines*data.Lz + (j - 1) * data.Lz + k;
                    values[counter++] = - pow(data.Sigma,2) * pow(data.deltaT,2)* Ydivider;
                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k - 1;
                    values[counter++] = - pow(data.Sigma,2) * pow(data.deltaT,2)* Zdivider;

                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k;
                    values[counter] = 1 + 2 * pow(data.Sigma,2) * pow(data.deltaT,2)* (Xdivider + Ydivider + Zdivider);
                    MainDiag[i*Ylines*data.Lz + j * data.Lz + k] = 1 / values[counter];
                    counter++;

                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k + 1;
                    values[counter++] = - pow(data.Sigma,2) * pow(data.deltaT,2)* Zdivider;
                    column_num[counter] = i * Ylines*data.Lz + (j + 1) * data.Lz + k;
                    values[counter++] = - pow(data.Sigma,2) * pow(data.deltaT,2)* Ydivider;
                    column_num[counter] = (i + 1) * Ylines*data.Lz + j * data.Lz + k;
                    values[counter++] = - pow(data.Sigma,2) * pow(data.deltaT,2)* Xdivider;
                }
            }
    line_first[elements_per_process] = counter;

    // Calculation process ---
    unsigned i, j, k, count;
    for (count = 0; count <= count_threshold; count++) {
        if (rank==0) out_time = omp_get_wtime();

        flag = 1;
#pragma omp parallel for
        for (i = 0; i < elements_per_process; i++) {
            prev[i] = current[i];
            current[i] = next[i];
        }

        while (flag) {
            flag = 0;
            CsrMult(values, column_num, line_first, next, elements_per_process, discrepancy);

#pragma omp parallel for
            for (i = 0; i < elements_per_process; i++) {
                discrepancy[i] = discrepancy[i] - 2 * current[i] + prev[i];
                discrepancy[i] = fabs(discrepancy[i]);
                if (discrepancy[i] > delta) {
                    flag = 1;
                    break;
                }
            }

            result = 0;
            MPI_Allreduce(&flag, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (result == 0) {
                break;
            }
            flag = 1;

            // We must set discrepancy to next to calculate sum from previous calculated points.
            for (i = 0; i < elements_per_process; i++)
                discrepancy[i] = next[i];

            SendBorderValues(discrepancy, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num, border_size);

            for (i = 0; i < elements_per_process; i++) {
                sum = 0;
                for (j = 0; j < elements_per_process; j++) {
                    // TODO: remove brute force sum calculating
                    if (i != j) sum += CsrAccess(values, column_num, line_first, i, j, elements_per_process)*discrepancy[j];
                }
                next[i] = MainDiag[i] * (2*current[i] - prev[i] - sum);
            }
        }

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