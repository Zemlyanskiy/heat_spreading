#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <windows.h>

#define START_STRING 0
#define EILER 0
#define RUNGE_KUTTA 1
#define DEBUG 0

// Array functions

struct Array3D{
    // Number of points (x, y, z)
    unsigned Lx, Ly, Lz;
    double* values;
};

inline void initArray3D(struct Array3D* arr, int X, int Y, int Z) {
    arr->Lx = X;
    arr->Ly = Y;
    arr->Lz = Z;
    arr->values = (double*)calloc(X*Y*Z, sizeof(double));
}

inline double get(struct Array3D *arr, int x, int y, int z) {
    return arr->values[z + y*arr->Lz + x*arr->Ly*arr->Lz];
}

inline void set(struct Array3D *arr, int x, int y, int z, double value) {
    arr->values[z + y*arr->Lz + x * arr->Ly*arr->Lz] = value;
}

// Input functions

struct input {
    double t;// Start time (t)
    double T;// End time (T)
    double deltaT;// Gap between measurements (dt)
    double Xmin;
    double Xmax;
    double Ymin;
    double Ymax;
    double Zmin;
    double Zmax;
    double Sigma; // The coefficient of quation
    double deltaOut;// Delta of output
    struct Array3D arr;
};

inline struct input ReadInput(char* PATH) {
    FILE *file;
    struct input temp = { .t = 0,.T = 0,.deltaT = 0,
                          .Xmin = 0,.Xmax = 0,.arr.Lx = 0,
                          .Ymin = 0,.Ymax = 0,.arr.Ly = 0,
                          .Zmin = 0,.Zmax = 0,.arr.Lz = 0,
                          .Sigma = 0,.deltaOut = 0,.arr.values = 0 };
    char* token;
    int pos=0;
    double buf=0;
    char* string = (char*)calloc(10000000,sizeof(char));

    fopen_s(&file, PATH, "r");
    if (file == NULL) {
            printf("Can`t find file\n");
    }
    else {
        for (int i = 0; fgets(string, 10000000 * sizeof(char), file) != NULL; i++)
            switch (i) {
            case 0:
            {
                sscanf_s(string, "t=%lf", &temp.t);
                break;
            }
            case 1:
            {
                sscanf_s(string, "T=%lf", &temp.T);
                break;
            }
            case 2:
            {
                sscanf_s(string, "deltaT=%lf", &temp.deltaT);
                break;
            }
            case 3:
            {
                sscanf_s(string, "XMin=%lf", &temp.Xmin);
                break;
            }
            case 4:
            {
                sscanf_s(string, "XMax=%lf", &temp.Xmax);
                break;
            }
            case 5:
            {
                sscanf_s(string, "YMin=%lf", &temp.Ymin);
                break;
            }
            case 6:
            {
                sscanf_s(string, "YMax=%lf", &temp.Ymax);
                break;
            }
            case 7:
            {
                sscanf_s(string, "ZMin=%lf", &temp.Zmin);
                break;
            }
            case 8:
            {
                sscanf_s(string, "ZMax=%lf", &temp.Zmax);
                break;
            }
            case 9:
            {
                sscanf_s(string, "Lx=%d", &temp.arr.Lx);
                break;
            }
            case 10:
            {
                sscanf_s(string, "Ly=%d", &temp.arr.Ly);
                break;
            }
            case 11:
            {
                sscanf_s(string, "Lz=%d", &temp.arr.Lz);
                break;
            }
            case 12:
            {
                sscanf_s(string, "Sigma=%lf", &temp.Sigma);
                break;
            }
            case 13:
            {
                sscanf_s(string, "deltaOut=%lf", &temp.deltaOut);
                break;
            }
            default:
#if START_STRING
                printf("Mode: start string, Problem: file allready contain start string\n");
#endif
                initArray3D(&temp.arr, temp.arr.Lx, temp.arr.Ly, temp.arr.Lz);
                token = strtok(string, " ");
                while (token != NULL)
                {
                    buf = strtof(token, NULL);
                    temp.arr.values[pos] = buf;
                    token = strtok(NULL, " ");
                    pos++;
                }
            }
    }

    fclose(file);
    return temp;
}

// Output functions

inline void OutputArr(char* PATH, double* arr, int length) {
    FILE* file;
    fopen_s(&file, PATH, "a+");
    if ( file == NULL)
        printf("Cant open file\n");
    else {
        fprintf(file, "\n");
        for (int i = 0; i < length; i++)
            fprintf(file, "%lf ", arr[i]);
        fclose(file);
    }
}

#if START_STRING

double CalcFunc(double Xcoord, double Ycoord, double Zcoord) {
    if (Xcoord >= -0.5 && Xcoord <= 0.5 &&
        Ycoord >= -0.5 && Ycoord <= 0.5 &&
        Zcoord >= -0.5 && Zcoord <= 0.5) {
        return cos(Xcoord*3.141592) + cos(Ycoord*3.141592) + cos(Zcoord*3.141592);
    }
    return 0;
}

#endif /*START_STRING*/

inline int OutputCurrentTime(char* PATH, double currT) {
    FILE* file;
    int ERR_CODE, num_of_space = 15 - 8;
    int tmp = (int)currT;
    fopen_s(&file, PATH, "r+");
    if (file == NULL)
        printf("Cant open file\n");
    else {
        while (tmp > 10) {
            tmp /= 10;
            num_of_space--;
        }

        fseek(file, 0, SEEK_SET);
        fputs("t=", file);
        for (int i = 0; i<num_of_space; i++)
            fputs(" ", file);
        fprintf(file, "%f", currT);
        ERR_CODE = ftell(file);
        fclose(file);
    }
    return 0;
}

// MPI Functions

inline void SendBorderValues( struct Array3D* points , unsigned rank, 
    unsigned process_per_axis, unsigned Xcube_pos, unsigned Ycube_pos
) {
    MPI_Status Status;
    int section_size = points->Ly * points ->Lz;
    // Send
    if (Xcube_pos > 0) {
        // send X border section to previos process
        MPI_Send(points->values + section_size, section_size, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
    }
    if (Xcube_pos < process_per_axis - 1) {
        // send X border section to next process
        MPI_Send(points->values + (points->Lx - 2)*section_size, section_size, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
    }

    if (Ycube_pos > 0) {
        // send Y border section to previos process
        for (unsigned i = 0; i < points->Lx; i++)
            MPI_Send(points->values + points->Lz + section_size * i, points->Lz, MPI_DOUBLE, rank - process_per_axis, rank, MPI_COMM_WORLD);
    }
    if (Ycube_pos < process_per_axis - 1) {
        // send Y border section to next process
        for (unsigned i = 0; i < points->Lx; i++)
            MPI_Send(points->values + (i + 1)*section_size - points->Lz*2, points->Lz, MPI_DOUBLE, rank + process_per_axis, rank, MPI_COMM_WORLD);
    }

    // Recieve
    if (Xcube_pos > 0) {
        // get X border section from previos process
        MPI_Recv(points->values, section_size, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &Status);
    }
    if (Xcube_pos < process_per_axis - 1) {
        // get X border section from next process
        MPI_Recv(points->values + (points->Lx - 1)*section_size,
                 section_size, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &Status);
    }

    if (Ycube_pos > 0) {
        // get Y border section from previos process
        for (unsigned i = 0; i < points->Lx; i++)
            MPI_Recv(points->values + section_size * i, points->Lz,
                     MPI_DOUBLE, rank - process_per_axis, rank - process_per_axis, MPI_COMM_WORLD, &Status);
    }
    if (Ycube_pos < process_per_axis - 1) {
        // get Y border section from next process
        for (unsigned i = 0; i < points->Lx; i++)
            MPI_Recv(points->values + (i + 1)*section_size - points->Lz,
                     points->Lz, MPI_DOUBLE, rank + process_per_axis, rank + process_per_axis, MPI_COMM_WORLD, &Status);
    }
}


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

    double Xstep = (data.Xmax - data.Xmin) / data.arr.Lx;
    double Ystep = (data.Ymax - data.Ymin) / data.arr.Ly;
    double Zstep = (data.Zmax - data.Zmin) / data.arr.Lz;

    double Xdivider = 1 / (Xstep*Xstep);
    double Ydivider = 1 / (Ystep*Ystep);
    double Zdivider = 1 / (Zstep*Zstep);

#if START_STRING
    if (!rank) {
        double posX = data.Xmin;
        double posY = data.Ymin;
        double posZ = data.Zmin;

        struct Array3D points;
        initArray3D(&points, data.arr.Lx, data.arr.Ly, data.arr.Lz);

        for (unsigned i = 0; i < data.arr.Lx; i++) {
            posX += Xstep;
            posY = data.Ymin;
            for (unsigned j = 0; j < data.arr.Ly; j++) {
                posY += Ystep;
                posZ = data.Zmin;
                for (unsigned k = 0; k < data.arr.Lz; k++) {
                    posZ += Zstep;
                    set(&points, i, j, k, CalcFunc(posX, posY, posZ));
                }
            }
        }

        OutputArr(inputFile, points.values, points.Lx * points.Ly * points.Lz);
    }
    MPI_Finalize();
    return 0;
#endif /*START STRING*/

    // Variables initialization ---
    double current_time = data.t;
    unsigned out_freq = (unsigned)(data.T - data.t) / data.deltaOut;
    unsigned count_threshold = (unsigned)(data.T - data.t) / data.deltaT;
    out_freq = count_threshold / out_freq;

    unsigned process_per_axis = (unsigned)sqrt(proc_num);
    unsigned Xlines = data.arr.Lx / process_per_axis;
    unsigned Ylines = data.arr.Ly / process_per_axis;
    unsigned Xcube_pos = rank % process_per_axis;
    unsigned Ycube_pos = rank / process_per_axis;

    // Add borders ---
    if (Xcube_pos > 0)
        Xlines++;
    if (Xcube_pos < process_per_axis - 1)
        Xlines++;
    if (Ycube_pos > 0)
        Ylines++;
    if (Ycube_pos < process_per_axis - 1)
        Ylines++;

    // Add remainder ---
    if (Xcube_pos == process_per_axis - 1)
        Xlines += data.arr.Lx % process_per_axis;
    if (Ycube_pos == process_per_axis - 1)
        Ylines += data.arr.Ly % process_per_axis;
    
    // Standart calculation arrays ---
    struct Array3D points;
    struct Array3D prev;
    initArray3D(&points, Xlines, Ylines, data.arr.Lz);
    initArray3D(&prev, Xlines, Ylines, data.arr.Lz);

    // MPI communication variables ---
    unsigned elements_per_process = Xlines * Ylines * data.arr.Lz;

    unsigned Xstart_copy_line = Xcube_pos * data.arr.Lx / process_per_axis;
    unsigned Ystart_copy_line = Ycube_pos * data.arr.Ly / process_per_axis;

    if (Xcube_pos > 0)
        Xstart_copy_line -= 1;
    if (Ycube_pos > 0)
        Ystart_copy_line -= 1;

    // Values for return collection (used only on root rank)
    int* x_return_lines = (int*)malloc(proc_num * sizeof(int));
    int* y_return_lines = (int*)malloc(proc_num * sizeof(int));
    int* x_return_start = (int*)malloc(proc_num * sizeof(int));
    int* y_return_start = (int*)malloc(proc_num * sizeof(int));
    struct Array3D* buffer = (struct Array3D*)malloc(sizeof(struct Array3D)*proc_num);

    // Create on 1st process array of values for collect return data ---
    if (rank == 0) {
        x_return_lines[0] = Xlines;
        y_return_lines[0] = Ylines;
        x_return_start[0] = Xstart_copy_line;
        y_return_start[0] = Ystart_copy_line;
        initArray3D(&buffer[0], x_return_lines[0], y_return_lines[0], data.arr.Lz);

        for (int i = 1; i < proc_num; i++) {
            MPI_Recv(x_return_lines + i, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(y_return_lines + i, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(x_return_start + i, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            MPI_Recv(y_return_start + i, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            initArray3D(&buffer[i], x_return_lines[i], y_return_lines[i], data.arr.Lz);
        }

    }
    else
    {
        MPI_Send(&Xlines, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Ylines, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Xstart_copy_line, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
        MPI_Send(&Ystart_copy_line, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
    }

    // Array for print result to file ---
    struct Array3D print_array;
    if (rank == 0) {
        initArray3D(&print_array, data.arr.Lx, data.arr.Ly, data.arr.Lz);
    }

    // Env preparation for OpenMP ---
    omp_set_num_threads(number_of_threads);
    double time = 0;
    double out_time;

    // Read start state of array ---
    for (unsigned i = 0; i < Xlines; i++)
        for (unsigned j = 0; j < Ylines; j++)
            for (unsigned k = 0; k < data.arr.Lz; k++) {
                set(&points, i, j, k, get(&data.arr, Xstart_copy_line + i, Ystart_copy_line + j, k));
            }

    // Runge Kutta variables
    struct Array3D k1, k2, k3, k4, medium, derivative;

#if RUNGE_KUTTA
    initArray3D(&k1, Xlines, Ylines, data.arr.Lz);
    initArray3D(&k2, Xlines, Ylines, data.arr.Lz);
    initArray3D(&k3, Xlines, Ylines, data.arr.Lz);
    initArray3D(&k4, Xlines, Ylines, data.arr.Lz);
    initArray3D(&medium, Xlines, Ylines, data.arr.Lz);
    initArray3D(&derivative, Xlines, Ylines, data.arr.Lz);
#endif /*RUNGE_KUTTA*/

    if (DEBUG) {
        printf("rank: %d xpos: %d ypos: %d proc_per_ax: %d\n"
            "rank: %d xlines: %d ylines: %d zlines: %d\n"
            "rank: %d x_start_copy %d y_start_copy %d\n",
            rank, Xcube_pos, Ycube_pos, process_per_axis,
            rank, Xlines, Ylines, data.arr.Lz,
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
        if (!rank) out_time = omp_get_wtime();

#if EILER
#pragma omp parallel for
        for (i = 0; i < elements_per_process; i++)
             prev.values[i] = points.values[i];
 #pragma omp parallel for
         for (i = 1; i < Xlines - 1; i++)
             for (j = 1; j < Ylines - 1; j++)
                 for (k = 1; k < points.Lz - 1; k++) {
                     set(&points, i, j, k,
                         get(&prev, i, j, k) + data.Sigma*data.deltaT*(
                         (get(&prev, i - 1, j, k) - 2 * get(&prev, i, j, k) + get(&prev, i + 1, j, k)) * Xdivider +
                             (get(&prev, i, j - 1, k) - 2 * get(&prev, i, j, k) + get(&prev, i, j + 1, k)) * Ydivider +
                             (get(&prev, i, j, k - 1) - 2 * get(&prev, i, j, k) + get(&prev, i, j, k + 1)) * Zdivider));
                    // points[i][j][k] = prev[i][j][k] + data.Sigma*data.deltaT*((prev[i-1][j][k] - 2*prev[i][j][k] + prev[i+1][j][k]) * Xdivider +
                    //                                                           (prev[i][j-1][k] - 2*prev[i][j][k] + prev[i][j+1][k]) * Ydivider +
                    //                                                           (prev[i][j][k-1] - 2*prev[i][j][k] + prev[i][j][k+1]) * Zdivider);
                }
#endif /*EILER*/
#if RUNGE_KUTTA
#pragma omp parallel for
         for (i = 0; i < elements_per_process; i++)
             prev.values[i] = points.values[i];

// k1 calculation ---
#pragma omp parallel for
         for (i = 1; i < Xlines - 1; i++)
             for (j = 1; j < Ylines - 1; j++)
                 for (k = 1; k < points.Lz - 1; k++) {
                     set(&k1, i, j, k, data.Sigma*data.deltaT*(
                         (get(&prev, i - 1, j, k) - 2 * get(&prev, i, j, k) + get(&prev, i + 1, j, k)) * Xdivider +
                         (get(&prev, i, j - 1, k) - 2 * get(&prev, i, j, k) + get(&prev, i, j + 1, k)) * Ydivider +
                         (get(&prev, i, j, k - 1) - 2 * get(&prev, i, j, k) + get(&prev, i, j, k + 1)) * Zdivider));
                 }
#pragma omp parallel for
         for (i = 0; i < Xlines; i++)
             for (j = 0; j < Ylines; j++)
                 for (k = 0; k < points.Lz; k++) {
                     set(&medium, i, j, k, get(&prev, i, j, k) + get(&k1, i, j, k) * 0.5);
                 }

// k2 calculation ---
#pragma omp parallel for
         for (i = 1; i < Xlines - 1; i++)
             for (j = 1; j < Ylines - 1; j++)
                 for (k = 1; k < points.Lz - 1; k++) {
                     set(&k2, i, j, k, data.Sigma*data.deltaT*(
                         (get(&medium, i - 1, j, k) - 2 * get(&medium, i, j, k) + get(&medium, i + 1, j, k)) * Xdivider +
                         (get(&medium, i, j - 1, k) - 2 * get(&medium, i, j, k) + get(&medium, i, j + 1, k)) * Ydivider +
                         (get(&medium, i, j, k - 1) - 2 * get(&medium, i, j, k) + get(&medium, i, j, k + 1)) * Zdivider));
                 }
#pragma omp parallel for
         for (i = 0; i < Xlines; i++)
             for (j = 0; j < Ylines; j++)
                 for (k = 0; k < points.Lz; k++) {
                     set(&medium, i, j, k, get(&prev, i, j, k) + get(&k2, i, j, k) * 0.5);
                 }

// k3 calculation ---
#pragma omp parallel for
         for (i = 1; i < Xlines - 1; i++)
             for (j = 1; j < Ylines - 1; j++)
                 for (k = 1; k < points.Lz - 1; k++) {
                     set(&k3, i, j, k, data.Sigma*data.deltaT*(
                         (get(&medium, i - 1, j, k) - 2 * get(&medium, i, j, k) + get(&medium, i + 1, j, k)) * Xdivider +
                         (get(&medium, i, j - 1, k) - 2 * get(&medium, i, j, k) + get(&medium, i, j + 1, k)) * Ydivider +
                         (get(&medium, i, j, k - 1) - 2 * get(&medium, i, j, k) + get(&medium, i, j, k + 1)) * Zdivider));
                 }
#pragma omp parallel for
         for (i = 0; i < Xlines; i++)
             for (j = 0; j < Ylines; j++)
                 for (k = 0; k < points.Lz; k++) {
                     set(&medium, i, j, k, get(&prev, i, j, k) + get(&k2, i, j, k));
                 }

// k4 calculation ---
#pragma omp parallel for
         for (i = 1; i < Xlines - 1; i++)
             for (j = 1; j < Ylines - 1; j++)
                 for (k = 1; k < points.Lz - 1; k++) {
                     set(&k4, i, j, k, data.Sigma*data.deltaT*(
                         (get(&medium, i - 1, j, k) - 2 * get(&medium, i, j, k) + get(&medium, i + 1, j, k)) * Xdivider +
                         (get(&medium, i, j - 1, k) - 2 * get(&medium, i, j, k) + get(&medium, i, j + 1, k)) * Ydivider +
                         (get(&medium, i, j, k - 1) - 2 * get(&medium, i, j, k) + get(&medium, i, j, k + 1)) * Zdivider));
                 }

// derivative calculation ---
#pragma omp parallel for
         for (i = 1; i < Xlines - 1; i++)
             for (j = 1; j < Ylines - 1; j++)
                 for (k = 1; k < points.Lz - 1; k++) {
                     set(&derivative, i, j, k, (get(&k1, i, j, k)+ get(&k2, i, j, k)+ get(&k3, i, j, k)+ get(&k4, i, j, k))* 0.1666666666);
                 }

// final result ---
#pragma omp parallel for
         for (i = 1; i < Xlines - 1; i++)
             for (j = 1; j < Ylines - 1; j++)
                 for (k = 1; k < points.Lz - 1; k++) {
                     set(&points, i, j, k, get(&prev, i, j, k) + get(&derivative, i, j, k));
                 }
#endif /*RUNGE_KUTTA*/

        SendBorderValues(&points, rank, process_per_axis, Xcube_pos, Ycube_pos);

        // Delete output time from time calculating ---
        if (!rank) time += omp_get_wtime() - out_time;

        // Print ---
        if (count%out_freq == 0 && count) {
            if (rank == 0) {
                //Can't send-receive from 0 process to himself
                for (unsigned i = 0; i < points.Lx; i++)
                    for (unsigned j = 0; j < points.Ly; j++)
                        for (unsigned k = 0; k < points.Lz; k++) {
                            set(&print_array, i, j, k, get(&points, i, j, k));
                        }

                for (unsigned sender_number = 1; sender_number < proc_num; sender_number++) {
                    MPI_Recv((buffer + sender_number) -> values, x_return_lines[sender_number] * y_return_lines[sender_number] * data.arr.Lz,
                        MPI_DOUBLE, sender_number, sender_number, MPI_COMM_WORLD, &status);

                    for (int i = 0; i < x_return_lines[sender_number]; i++)
                        for (int j = 0; j < y_return_lines[sender_number]; j++)
                            for (int k = 0; k < data.arr.Lz; k++) {
                                set(&print_array, x_return_start[sender_number] + i, y_return_start[sender_number] + j, k, get(buffer + sender_number, i, j, k));
                            }
                }
            
                OutputArr(inputFile, print_array.values, print_array.Lx*print_array.Ly*print_array.Lz);
                current_time += data.deltaOut;
                OutputCurrentTime(inputFile, current_time);
            }
            else {
                MPI_Send(points.values, points.Lx*points.Ly*points.Lz, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
            }
        }
    }

    // End ---
    if (!rank) printf("Time is: %lf\n", time);
    MPI_Finalize();
    return 0;
}