#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include<mpi.h>

#define START_STRING 0
#define EILER 1
#define DEBUG 1

// Array functions

struct Array3D{
    // Number of points (x, y, z)
    int Lx;
    int Ly;
    int Lz;
    double* values;
};

inline void initArray3D(struct Array3D* arr, int X, int Y, int Z) {
    arr->Lx = X;
    arr->Ly = Y;
    arr->Lz = Z;
    arr->values = (double*)calloc(X*Y*Z, sizeof(double));
}

inline double get(struct Array3D *arr, int x, int y, int z) {
    return arr->values[z + y*arr->Ly + x*arr->Ly*arr->Lz];
}

inline void set(struct Array3D *arr, int x, int y, int z, double value) {
    arr->values[z + y * arr->Ly + x * arr->Ly*arr->Lz] = value;
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
                initArray3D( &temp.arr, temp.arr.Lx, temp.arr.Ly, temp.arr.Lz );
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

// CSR Functions
#if 0
inline void CsrMult(double* Values, int* ColumnNum, int* LineFirst, double * Array, const int size, double* result)
{
    int i;
#pragma omp parallel for
    for (i = 0; i < size; i++) {
        result[i] = 0;
        for (int j = LineFirst[i]; j < LineFirst[i + 1]; j++)
            result[i] += Values[j] * Array[ColumnNum[j]];
    }
}

inline double CrsAccess(double* Values, int* ColumnNum, int* LineFirst, int i, int j, int length) {
    if (i >= length)
        return 0;
    for (int k = LineFirst[i]; k<LineFirst[i + 1]; k++)
        if (ColumnNum[k] == j) {
            return Values[k];
        }
    return 0;
}
#endif
// MPI Functions

inline void SendBorderValues(struct Array3D* points, int sendsections_for_rank, int rank, int proc_num) {
    if (proc_num != 1) {
        MPI_Status Status;
        int message_size = points->Ly * points->Lz;\

        if (rank)
        {
            // send second section to previos process
            MPI_Send(points->values + message_size, message_size, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
        }

        if (rank != proc_num - 1) {
            // send pre last section to next process
            MPI_Send( points->values + (sendsections_for_rank - 2)*message_size,
                      message_size, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
        }

        if (rank) {
            // get first section from previos process
            MPI_Recv(points->values, message_size, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &Status);
        }

        if (rank != proc_num - 1) {
            // get last section from next process
            MPI_Recv( points->values + (sendsections_for_rank - 1)*message_size,
                      message_size, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &Status);
        }
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

    // Processing command line arguments ---
    char* inputFile = NULL;
    unsigned number_of_threads = 1;
    {

        if (argc > 1)
            inputFile = argv[1];
        if (argc > 2)
            number_of_threads = strtol(argv[3], NULL, 0);
        if (argc > 3) {
            printf("Incorrect number of arguments");
            return 4;
        }
    }

    // Read input data and intialize common variables ---
    struct input data;
    data = ReadInput(inputFile);

    unsigned sections_per_process = data.arr.Lx / proc_num;
    unsigned remainder = data.arr.Lx % proc_num;
    unsigned elements_per_process = sections_per_process * data.arr.Ly * data.arr.Lz;

    unsigned count = 0;

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

    // Time variables initializing ---
    double current_time = data.t;
    unsigned out_freq = (int)(data.T - data.t) / data.deltaOut;
    unsigned count_threshold = (int)(data.T - data.t) / data.deltaT;
    out_freq = count_threshold / out_freq;

    // Env preparation for MPI ---
    int * sendsections = (int*)malloc(sizeof(int)*proc_num);
    int * displssections = (int*)malloc(sizeof(int)*proc_num);

    int * sendcounts = (int*)malloc(sizeof(int)*proc_num);
    int * displscounts = (int*)malloc(sizeof(int)*proc_num);

    int * returncounts = (int*)malloc(sizeof(int)*proc_num);
    int * returndisplscounts = (int*)malloc(sizeof(int)*proc_num);

    // Scatter preparation ---
    for (unsigned proc_rank = 0; proc_rank < proc_num; proc_rank++) {
        sendsections[proc_rank] = sections_per_process;
        displssections[proc_rank] = sections_per_process * proc_rank;
    }

    if (proc_num != 1) {
        for (unsigned proc_rank = 0; proc_rank < proc_num; proc_rank++) {
            // Distribute remainder
            if (proc_rank < remainder) {
                sendsections[proc_rank] += 1;
                for (unsigned sub_count = proc_rank + 1; sub_count < proc_num; sub_count++)
                    displssections[sub_count] += 1;
            }

            // Init return values
            returncounts[proc_rank] = sendsections[proc_rank] * data.arr.Ly * data.arr.Lz;
            returndisplscounts[proc_rank] = displssections[proc_rank] * data.arr.Ly * data.arr.Lz;


            if (proc_rank != 0) {
                // And section to begin of every process (exclude first)
                displssections[proc_rank] -= 1;
                sendsections[proc_rank] += 1;
            }
            if (proc_rank != proc_num - 1) {
                // Add section to end of every process (exclude last)
                sendsections[proc_rank] += 1;
            }
        }
    } 
    else {
        returncounts[0] = elements_per_process;
        returndisplscounts[0] = 0;
    }

    for (int proc_rank = 0; proc_rank < proc_num; proc_rank++) {
        sendcounts[proc_rank] = sendsections[proc_rank] * data.arr.Ly * data.arr.Lz;
        displscounts[proc_rank] = displssections[proc_rank] * data.arr.Ly * data.arr.Lz;
    }

    // Env preparation for OpenMP ---
    omp_set_num_threads(number_of_threads);
    double time = 0;
    double out_time;

    // Arrays initialization ---
    struct Array3D points;// array for caclulation results
    struct Array3D prev;// array for previos state
    struct Array3D print_array;// array for previos state
    initArray3D(&points, sendsections[rank], data.arr.Ly, data.arr.Lz);
    initArray3D(&prev, sendsections[rank], data.arr.Ly, data.arr.Lz);
    initArray3D(&print_array, data.arr.Lx, data.arr.Ly, data.arr.Lz);

    if (DEBUG && rank == 0) {
        printf("Lx: %d, Ly: %d, Lz: %d, SectionsPerProcess %d, Remainder:%d\n\n",
            data.arr.Lx, data.arr.Ly, data.arr.Lz, sections_per_process, remainder);
        for (int proc_rank = 0; proc_rank < proc_num; proc_rank++) {
            printf("Process %d:\n", proc_rank);
            printf("sendcounts: %d, sendsections: %d, displscounts: %d, displssections: %d\n",
                sendcounts[proc_rank], sendsections[proc_rank], displscounts[proc_rank], displssections[proc_rank]);
            printf("returncounts: %d, returndisplscounts %d \n", returncounts[proc_rank], returndisplscounts[proc_rank]);
            printf("-----------------------------------------------------\n");
        }
    }

    // Distribute data between processes ---
    MPI_Scatterv(data.arr.values, sendcounts, displscounts, MPI_DOUBLE, points.values,
        sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculation process ---
    unsigned i, j, k;
    for (count = 0; count <= count_threshold; count++) {
        if (!rank) out_time = omp_get_wtime();

#if EILER
#pragma omp parallel for
        for (i = 0; i < sendcounts[rank]; i++)
            prev.values[i] = points.values[i];
#pragma omp parallel for
        for (i = 1; i < sendsections[rank] - 1; i++)
            for (j = 1; j < points.Ly - 1; j++)
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

        SendBorderValues(&points, sendsections[rank], rank, proc_num);

        // Delete output time from time calculating ---
        if (!rank) time += omp_get_wtime() - out_time;

        // Print ---
        if (count%out_freq == 0 && count) {
            if (!rank) {
                //Can't send-receive from 0 process to himself
                for (int count = 0; count < returncounts[rank]; count++)
                    print_array.values[count] = points.values[count];
                for (int sender_number = 1; sender_number < proc_num; sender_number++)
                    MPI_Recv(&print_array.values[returndisplscounts[sender_number]], returncounts[sender_number],
                        MPI_DOUBLE, sender_number, sender_number, MPI_COMM_WORLD, &status);

                OutputArr(inputFile, print_array.values, data.arr.Lx*data.arr.Ly*data.arr.Lz);
                current_time += data.deltaOut;
                OutputCurrentTime(inputFile, current_time);
            }
            else {
                MPI_Send(points.values + data.arr.Ly * data.arr.Lz, returncounts[rank], MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
            }
        }
    }

    // End ---
    if (!rank) printf("Time is: %lf\n", time);
    MPI_Finalize();
    return 0;
}