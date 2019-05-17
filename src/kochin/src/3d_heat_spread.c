#include "common.h"
#include "csr_api.h"

#define DEBUG 0

#define EILER 1
#define RUNGE_KUTTA 0
#define MATRIX_EILER 0
#define MATRIX_RUNGE_KUTTA 0
#define YAKOBI 0

int main(int argc, char* argv[])
{
    // MPI initialization ---
    unsigned rank, proc_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

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

    double Xstep = (data.Xmax - data.Xmin) / data.Lx;
    double Ystep = (data.Ymax - data.Ymin) / data.Ly;
    double Zstep = (data.Zmax - data.Zmin) / data.Lz;

    double Xdivider = 1 / (Xstep*Xstep);
    double Ydivider = 1 / (Ystep*Ystep);
    double Zdivider = 1 / (Zstep*Zstep);

    // Variables initialization ---
    if ( data.t >= data.T) {
        if (rank == 0) { printf("Start time >= end time, please clear io file\n"); }
        MPI_Finalize();
        return 0;
    }
    double current_time = data.t;
    unsigned out_freq = (unsigned)((data.T - data.t) / data.deltaOut);
    unsigned count_threshold = (unsigned)((data.T - data.t) / data.deltaT);
    out_freq = count_threshold / out_freq;

    unsigned process_per_axis = (unsigned)sqrt(proc_num);
    unsigned Xlines = data.Lx / process_per_axis;
    unsigned Ylines = data.Ly / process_per_axis;
    unsigned Xcube_pos = rank % process_per_axis;
    unsigned Ycube_pos = rank / process_per_axis;

    // Add remainder ---
    if (Xcube_pos == process_per_axis - 1)
        Xlines += data.Lx % process_per_axis;
    if (Ycube_pos == process_per_axis - 1)
        Ylines += data.Ly % process_per_axis;

    unsigned Xstart_copy_line = Xcube_pos * data.Lx / process_per_axis;
    unsigned Ystart_copy_line = Ycube_pos * data.Ly / process_per_axis;

    // Add borders ---
    if (Xcube_pos > 0) {
        Xlines++;
        Xstart_copy_line--;
    }
    if (Xcube_pos < process_per_axis - 1)
        Xlines++;
    if (Ycube_pos > 0) {
        Ylines++;
        Ystart_copy_line--;
    }
    if (Ycube_pos < process_per_axis - 1)
        Ylines++;

    // Standart calculation arrays ---
    double *points, *prev;
    points = (double*)calloc(Xlines * Ylines * data.Lz, sizeof(double));
    prev = (double*)calloc(Xlines * Ylines * data.Lz, sizeof(double));

    // MPI communication variables ---
    unsigned elements_per_process = Xlines * Ylines * data.Lz;

    // Values for return collection (used only on root rank)
    int *x_return_lines, *y_return_lines, *x_return_start, *y_return_start;
    x_return_lines = (int*)malloc(proc_num * sizeof(int));
    y_return_lines = (int*)malloc(proc_num * sizeof(int));
    x_return_start = (int*)malloc(proc_num * sizeof(int));
    y_return_start = (int*)malloc(proc_num * sizeof(int));
    double** buffer = (double**)malloc(proc_num * sizeof(double*));

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
                Get(points, i, j, k, Xlines, Ylines, data.Lz) = Get(data.arr, Xstart_copy_line + i, Ystart_copy_line + j, k, data.Lx, data.Ly, data.Lz);
            }

    // Runge Kutta variables ---
#if RUNGE_KUTTA || MATRIX_RUNGE_KUTTA
    double *k1, *k2, *k3, *k4, *medium, *derivative;

    k1 = (double*) calloc( Xlines* Ylines * data.Lz, sizeof(double));
    k2 = (double*) calloc( Xlines* Ylines * data.Lz, sizeof(double));
    k3 = (double*) calloc( Xlines* Ylines * data.Lz, sizeof(double));
    k4 = (double*) calloc( Xlines* Ylines * data.Lz, sizeof(double));
    medium = (double*) calloc(Xlines* Ylines* data.Lz, sizeof(double));
    derivative = (double*) calloc (Xlines * Ylines * data.Lz, sizeof(double));
#endif /*RUNGE_KUTTA*/

    // CSR variables ---
#if  MATRIX_EILER ||  MATRIX_RUNGE_KUTTA
    double * values;
    int* column_num;
    int* line_first;

    unsigned borders_number = Xlines * Ylines * 2 + Ylines * data.Lz * 2 + Xlines * data.Lz * 2 - 8 - Xlines * 4 - Ylines * 4 - data.Lz * 4;
    unsigned values_number = (elements_per_process - borders_number) * 7 + borders_number;
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
                    values[counter++] = 1;
                }
                else {
                    // x y z central*3 z y x
                    column_num[counter] = (i - 1) * Ylines*data.Lz + j * data.Lz + k;
                    values[counter++] = data.Sigma * data.deltaT * Xdivider;
                    column_num[counter] = i * Ylines*data.Lz + (j - 1) * data.Lz + k;
                    values[counter++] = data.Sigma * data.deltaT * Ydivider;
                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k - 1;
                    values[counter++] = data.Sigma * data.deltaT * Zdivider;

                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k;
#if MATRIX_EILER
                    values[counter++] = 1 - 2 * data.Sigma * data.deltaT * (Xdivider + Ydivider + Zdivider);
#else
                    values[counter++] = -2 * data.Sigma * data.deltaT * (Xdivider + Ydivider + Zdivider);
#endif
                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k + 1;
                    values[counter++] = data.Sigma * data.deltaT * Zdivider;
                    column_num[counter] = i * Ylines*data.Lz + (j + 1) * data.Lz + k;
                    values[counter++] = data.Sigma * data.deltaT * Ydivider;
                    column_num[counter] = (i + 1) * Ylines*data.Lz + j * data.Lz + k;
                    values[counter++] = data.Sigma * data.deltaT * Xdivider;
                }
            }
    line_first[elements_per_process] = counter;

    if (DEBUG) {
        Sleep(100 * rank);
        printf("\n");
        printf("\n");
        for (unsigned i = 0; i < values_number; i++)
            printf("%f ", values[i]);
        printf("\n");
        for (unsigned i = 0; i < values_number; i++)
            printf("%d ", column_num[i]);
        printf("\n");
        for (unsigned i = 0; i < elements_per_process + 1; i++)
            printf("%d ", line_first[i]);
        printf("\n");
    }
#endif /*MATRIX_EILER ||  MATRIX_RUNGE_KUTTA*/

    // Yakobi variables ---
#if YAKOBI
    double* values;
    int* column_num;
    int* line_first;
    int flag, result;
    double delta, sum;
    double *MainDiag = NULL;
    double* discrepancy;
    discrepancy = (double*)calloc(Xlines * Ylines * data.Lz,sizeof(double));

    unsigned borders_number = elements_per_process - (Xlines - 2)*(Ylines - 2)*(data.Lz - 2);
    unsigned values_number = (elements_per_process - borders_number) * 7 + borders_number;
    values = (double*)calloc(values_number, sizeof(double));
    column_num = (unsigned*)calloc(values_number, sizeof(unsigned));
    line_first = (unsigned*)calloc(elements_per_process + 1, sizeof(unsigned));
    MainDiag = (double*)calloc(elements_per_process, sizeof(double));

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
                    values[counter++] = -data.Sigma * data.deltaT * Xdivider;
                    column_num[counter] = i * Ylines*data.Lz + (j - 1) * data.Lz + k;
                    values[counter++] = -data.Sigma * data.deltaT * Ydivider;
                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k - 1;
                    values[counter++] = -data.Sigma * data.deltaT * Zdivider;

                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k;
                    values[counter] = 1 + 2 * data.Sigma * data.deltaT * (Xdivider + Ydivider + Zdivider);
                    MainDiag[i*Ylines*data.Lz + j * data.Lz + k] = 1 / values[counter];
                    counter++;

                    column_num[counter] = i * Ylines*data.Lz + j * data.Lz + k + 1;
                    values[counter++] = -data.Sigma * data.deltaT * Zdivider;
                    column_num[counter] = i * Ylines*data.Lz + (j + 1) * data.Lz + k;
                    values[counter++] = -data.Sigma * data.deltaT * Ydivider;
                    column_num[counter] = (i + 1) * Ylines*data.Lz + j * data.Lz + k;
                    values[counter++] = -data.Sigma * data.deltaT * Xdivider;
                }
            }
    line_first[elements_per_process] = counter;

    delta = 0.00001;
    flag = 0;

    if (DEBUG) {
        Sleep(100 * rank);
        printf("\n rank: %d", rank);
        printf("\n values: ");
        for (unsigned i = 0; i < values_number; i++)
            printf("%f ", values[i]);
        printf("\n column num: ");
        for (unsigned i = 0; i < values_number; i++)
            printf("%d ", column_num[i]);
        printf("\n line first: ");
        for (unsigned i = 0; i < elements_per_process + 1; i++)
            printf("%d ", line_first[i]);
        printf("\nMain Diag: ");

        for (unsigned i = 0; i < elements_per_process; i++)
            printf("%f ", MainDiag[i]);
        printf("\n");
        MPI_Finalize();
        return 0;
    }
#endif /*YAKOBI*/

    // Amount of computation per process
    if (DEBUG) {
        Sleep(1000 * rank);
        printf("rank: %d xpos: %d ypos: %d proc_per_ax: %d\n"
            "rank: %d xlines: %d ylines: %d zlines: %d\n"
            "rank: %d x_start_copy %d y_start_copy %d\n",
            rank, Xcube_pos, Ycube_pos, process_per_axis,
            rank, Xlines, Ylines, data.Lz,
            rank, Xstart_copy_line, Ystart_copy_line);

        if (rank == 0) {
            for (int i = 0; i < proc_num; i++) {
                printf("i: %d ret_x_l: %d, ret_y_l: %d ret_x_s: %d ret_y_s: %d\n",
                    i, x_return_lines[i], y_return_lines[i], x_return_start[i], y_return_start[i]);
            }
        }
    }

    // Calculation process ---
    unsigned i, j, k, count;
    for (count = 0; count <= count_threshold; count++) {
        if (rank==0) out_time = omp_get_wtime();

#if EILER
    //printf("count: %d rank: %d\n", count, rank);
    //fflush(stdout);
#pragma omp parallel for
        for (i = 0; i < elements_per_process; i++)
            prev[i] = points[i];
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(points, i, j, k, Xlines, Ylines, data.Lz) =
                        Get(prev, i, j, k, Xlines, Ylines, data.Lz) + data.Sigma*data.deltaT*(
                        (Get(prev, i - 1, j, k, Xlines, Ylines, data.Lz) - 2 * Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(prev, i + 1, j, k, Xlines, Ylines, data.Lz)) * Xdivider +
                        (Get(prev, i, j - 1, k, Xlines, Ylines, data.Lz) - 2 * Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(prev, i, j + 1, k, Xlines, Ylines, data.Lz)) * Ydivider +
                        (Get(prev, i, j, k - 1, Xlines, Ylines, data.Lz) - 2 * Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(prev, i, j, k + 1, Xlines, Ylines, data.Lz)) * Zdivider);
                    // points[i][j][k] = prev[i][j][k] + data.Sigma*data.deltaT*((prev[i-1][j][k] - 2*prev[i][j][k] + prev[i+1][j][k]) * Xdivider +
                    //                                                           (prev[i][j-1][k] - 2*prev[i][j][k] + prev[i][j+1][k]) * Ydivider +
                    //                                                           (prev[i][j][k-1] - 2*prev[i][j][k] + prev[i][j][k+1]) * Zdivider);
                }
#endif /*EILER*/
#if RUNGE_KUTTA
#pragma omp parallel for
        for (i = 0; i < elements_per_process; i++)
            prev[i] = points[i];

        // k1 calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(k1, i, j, k, Xlines, Ylines, data.Lz) = data.Sigma*data.deltaT*(
                        (Get(prev, i - 1, j, k, Xlines, Ylines, data.Lz) - 2 * Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(prev, i + 1, j, k, Xlines, Ylines, data.Lz)) * Xdivider +
                        (Get(prev, i, j - 1, k, Xlines, Ylines, data.Lz) - 2 * Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(prev, i, j + 1, k, Xlines, Ylines, data.Lz)) * Ydivider +
                        (Get(prev, i, j, k - 1, Xlines, Ylines, data.Lz) - 2 * Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(prev, i, j, k + 1, Xlines, Ylines, data.Lz)) * Zdivider);
                }
        SendBorderValues(k1, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num);

#pragma omp parallel for
        for (i = 0; i < Xlines; i++)
            for (j = 0; j < Ylines; j++)
                for (k = 0; k < data.Lz; k++) {
                    Get(medium, i, j, k, Xlines, Ylines, data.Lz) =  Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(k1, i, j, k, Xlines, Ylines, data.Lz) * 0.5;
                }

        // k2 calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(k2, i, j, k, Xlines, Ylines, data.Lz) = data.Sigma*data.deltaT*(
                        (Get(medium, i - 1, j, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i + 1, j, k, Xlines, Ylines, data.Lz)) * Xdivider +
                        (Get(medium, i, j - 1, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j + 1, k, Xlines, Ylines, data.Lz)) * Ydivider +
                        (Get(medium, i, j, k - 1, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j, k + 1, Xlines, Ylines, data.Lz)) * Zdivider);
                }
        SendBorderValues(k2, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num);

#pragma omp parallel for
        for (i = 0; i < Xlines; i++)
            for (j = 0; j < Ylines; j++)
                for (k = 0; k < data.Lz; k++) {
                    Get(medium, i, j, k, Xlines, Ylines, data.Lz) = Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(k2, i, j, k, Xlines, Ylines, data.Lz) * 0.5;
                }

        // k3 calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(k3, i, j, k, Xlines, Ylines, data.Lz) = data.Sigma*data.deltaT*(
                        (Get(medium, i - 1, j, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i + 1, j, k, Xlines, Ylines, data.Lz)) * Xdivider +
                        (Get(medium, i, j - 1, k, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j + 1, k, Xlines, Ylines, data.Lz)) * Ydivider +
                        (Get(medium, i, j, k - 1, Xlines, Ylines, data.Lz) - 2 * Get(medium, i, j, k, Xlines, Ylines, data.Lz) + Get(medium, i, j, k + 1, Xlines, Ylines, data.Lz)) * Zdivider);
                }
        SendBorderValues(k3, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num);

#pragma omp parallel for
        for (i = 0; i < Xlines; i++)
            for (j = 0; j < Ylines; j++)
                for (k = 0; k < data.Lz; k++) {
                    Get(medium, i, j, k, Xlines, Ylines, data.Lz) = Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(k3, i, j, k, Xlines, Ylines, data.Lz);
                }

        // k4 calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(k4, i, j, k, Xlines, Ylines, data.Lz) = data.Sigma*data.deltaT*(
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
                    Get(points, i, j, k, Xlines, Ylines, data.Lz) = Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(derivative, i, j, k, Xlines, Ylines, data.Lz);
                }
#endif /*RUNGE_KUTTA*/
#if MATRIX_EILER
#pragma omp parallel for
        for (i = 0; i < elements_per_process; i++)
            prev[i] = points[i];

        CsrMult(values, column_num, line_first, prev, elements_per_process, points);
#endif /*MATRIX_EILER*/
#if MATRIX_RUNGE_KUTTA
#pragma omp parallel for
        for (i = 0; i < elements_per_process; i++)
            prev[i] = points[i];

        // k1 calculation ---
        CsrMult(values, column_num, line_first, prev, elements_per_process, k1);
        SendBorderValues(k1, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num);
#pragma omp parallel for
        for (i = 0; i < Xlines; i++)
            for (j = 0; j < Ylines; j++)
                for (k = 0; k < data.Lz; k++) {
                    Get(medium, i, j, k, Xlines, Ylines, data.Lz) = Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(k1, i, j, k, Xlines, Ylines, data.Lz) * 0.5;
                }

        // k2 calculation ---
        CsrMult(values, column_num, line_first, medium, elements_per_process, k2);
        SendBorderValues(k2, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num);
#pragma omp parallel for
        for (i = 0; i < Xlines; i++)
            for (j = 0; j < Ylines; j++)
                for (k = 0; k < data.Lz; k++) {
                    Get(medium, i, j, k, Xlines, Ylines, data.Lz) = Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(k2, i, j, k, Xlines, Ylines, data.Lz) * 0.5;
                }

        // k3 calculation ---
        CsrMult(values, column_num, line_first, medium, elements_per_process, k3);
        SendBorderValues(k3, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num);
#pragma omp parallel for
        for (i = 0; i < Xlines; i++)
            for (j = 0; j < Ylines; j++)
                for (k = 0; k < data.Lz; k++) {
                    Get(medium, i, j, k, Xlines, Ylines, data.Lz) = Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(k3, i, j, k, Xlines, Ylines, data.Lz);
                }

        // k4 calculation ---
        CsrMult(values, column_num, line_first, medium, elements_per_process, k4);

        // derivative calculation ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(derivative, i, j, k, Xlines, Ylines, data.Lz) = (Get(k1, i, j, k, Xlines, Ylines, data.Lz) +
                                                                      2 * Get(k2, i, j, k, Xlines, Ylines, data.Lz) +
                                                                      2 * Get(k3, i, j, k, Xlines, Ylines, data.Lz) +
                                                                          Get(k4, i, j, k, Xlines, Ylines, data.Lz))* 0.1666666666;
                }

        // final result ---
#pragma omp parallel for
        for (i = 1; i < Xlines - 1; i++)
            for (j = 1; j < Ylines - 1; j++)
                for (k = 1; k < data.Lz - 1; k++) {
                    Get(points, i, j, k, Xlines, Ylines, data.Lz) = Get(prev, i, j, k, Xlines, Ylines, data.Lz) + Get(derivative, i, j, k, Xlines, Ylines, data.Lz);
                }
#endif /*MATRIX_RUNGE_KUTTA*/
#if YAKOBI
        flag = 1;
#pragma omp parallel for
        for (i = 0; i < elements_per_process; i++)
            prev[i] = points[i];

        while (flag) {
            flag = 0;
            CsrMult(values, column_num, line_first, points, elements_per_process, discrepancy,rank);

#pragma omp parallel for
            for (i = 0; i < elements_per_process; i++) {
                discrepancy[i] -= prev[i];
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

            // We must set discrepancy to points to calculate sum from previous calculated points.
            for (i = 0; i < elements_per_process; i++)
                discrepancy[i] = points[i];
            SendBorderValues(discrepancy, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num);

            for (i = 0; i < elements_per_process; i++) {
                sum = 0;
                for (j = 0; j < elements_per_process; j++) {
                    // TODO: remove brute force sum calculating
                    if (i != j) sum += CsrAccess(values, column_num, line_first, i, j, elements_per_process)*discrepancy[j];
                }
                points[i] = MainDiag[i] * (prev[i] - sum);
            }
        }
#endif

        SendBorderValues(points, Xlines, Ylines, data.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num);

        // Delete output time from time calculating ---
        if (rank==0) time += omp_get_wtime() - out_time;

        // Print ---
        if (count % out_freq == 0 && count) {
            if (rank == 0) {
                //Can't send-receive from 0 process to himself
                for (unsigned i = 0; i < Xlines; i++)
                    for (unsigned j = 0; j < Ylines; j++)
                        for (unsigned k = 0; k < data.Lz; k++) {
                            Get(print_array, i, j, k, data.Lx, data.Ly, data.Lz) = Get(points, i, j, k, Xlines, Ylines, data.Lz);
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
            else {
                MPI_Send(points, Xlines*Ylines*data.Lz, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
            }
        }
    }

    // End ---
    if (!rank) { printf("Time is: %lf\n", time); }
    MPI_Finalize();
    return 0;
}