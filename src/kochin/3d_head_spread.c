#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include<mpi.h>

#define START_STRING 0
#define EILER 1

// Array functions
struct Array3D{
    // Number of points (x, y, z)
    int Lx;
    int Ly;
    int Lz;
    double* values; //LAST ARRAY WITH VALUES
};

inline void initArray3D(struct Array3D* arr, int X, int Y, int Z) {
    arr->Lx = X;
    arr->Ly = Y;
    arr->Lz = Z;
    arr->values = (double*)calloc(X*Y*Z, sizeof(double));
}

inline double get(struct Array3D *arr, int x, int y, int z) {
    return arr->values[z + y*arr->Ly + x*arr->Ly*arr->Lx];
}

inline void set(struct Array3D *arr, int x, int y, int z, double value) {
    arr->values[z + y * arr->Ly + x * arr->Ly*arr->Lx] = value;
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
    double Sigma; //THE COEFFICIENT OF QUATION
    double deltaOut;//DELTA OF OUTPUT
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
//    char string[1000];

    fopen_s(&file, PATH, "r");
    if (file == NULL) {
            printf("Can`t find file\n");
    }
    else {
        for (int i = 0; fgets(string, 10000000*sizeof(char), file) != NULL; i++)
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

// Output function
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

#else

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

inline void SendBorderValues(struct Array3D* points, int sections_per_process, int rank, int proc_num) {
    if (proc_num != 1) {
        MPI_Status Status;
        struct Array3D border_matrix;
        initArray3D(&border_matrix, 1, points->Ly, points->Lz);
        int message_size = points->Ly * points->Lz;

        if (rank != proc_num - 1) {
            for (int j = 0; j < points->Ly; j++)
                for (int k = 0; k < points->Lz; k++)
                    set(&border_matrix, 0, j, k, get(points, sections_per_process - 1, j, k));

            // send pre last section to next process
            MPI_Send(border_matrix.values, message_size, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
        }

        if (rank)
        {
            for (int j = 0; j < points->Ly; j++)
                for (int k = 0; k < points->Lz; k++)
                    set(&border_matrix, 0, j, k, get(points, 1, j, k));

            // send second section to previos process
            MPI_Send(border_matrix.values, message_size, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
        }

        if (rank != proc_num - 1) {
            // get last section from next process
            MPI_Recv(border_matrix.values, message_size, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &Status);

            for (int j = 0; j < points->Ly; j++)
                for (int k = 0; k < points->Lz; k++)
                    set(points, sections_per_process, j, k, get(&border_matrix, 0, j, k));
        }

        if (rank) {
            // get first section from previos process
            MPI_Recv(border_matrix.values, message_size, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &Status);

            for (int j = 0; j < points->Ly; j++)
                for (int k = 0; k < points->Lz; k++)
                    set(points, 0, j, k, get(&border_matrix, 0, j, k));
        }
    }
}

#endif

int main(int argc, char* argv[])
{
    // MPI initialization
    int rank, proc_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Processing command line arguments
    char* inputFile = NULL;
    int NUMDER_OF_THREADS = 1;
    {

        if (argc > 1)
            inputFile = argv[1];
        if (argc > 2)
            NUMDER_OF_THREADS = strtol(argv[3], NULL, 0);
        if (argc > 3) {
            printf("Incorrect number of arguments");
            return 4;
        }
    }

    //Read input data and intialize common variables
    struct input data;
    data = ReadInput(inputFile);

    unsigned sections_per_process = data.arr.Lx / proc_num;
    unsigned elements_per_process = sections_per_process * data.arr.Ly * data.arr.Lz;

    struct Array3D points;// array for caclulation results
    struct Array3D prev;// array for previos state
    struct Array3D print_array;// array for previos state
    initArray3D(&points, data.arr.Lx, data.arr.Ly, data.arr.Lz);
    initArray3D(&prev, data.arr.Lx, data.arr.Ly, data.arr.Lz);
    initArray3D(&print_array, data.arr.Lx, data.arr.Ly, data.arr.Lz);

    int count = 0;

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

        for (int i = 0; i < data.arr.Lx; i++) {
            posX += Xstep;
            posY = data.Ymin;
            for (int j = 0; j < data.arr.Ly; j++) {
                posY += Ystep;
                posZ = data.Zmin;
                for (int k = 0; k < data.arr.Lz; k++) {
                    posZ += Zstep;
                    set(&points, i, j, k, CalcFunc(posX, posY, posZ));
                }
            }
        }

        OutputArr(inputFile, points.values, points.Lx * points.Ly * points.Lz);
    }
    MPI_Finalize();
    return 0;
#endif

    double OutTime = data.t;
    int OutCount = (int)(data.T - data.t) / data.deltaOut;
    int maxCount = (int)(data.T - data.t) / data.deltaT;
    OutCount = maxCount / OutCount;

    //ENV PREPARATION FOR MPI
    int * sendsections = malloc(sizeof(int)*proc_num);
    int * sendcounts = malloc(sizeof(int)*proc_num);
    int * displs = malloc(sizeof(int)*proc_num);
    for (int i = 0; i < proc_num; i++) {
        sendcounts[i] = 0;
        sendsections[i] = sections_per_process;
        displs[i] = sections_per_process * i;
    }
    //remainder of the division
    sendsections[proc_num - 1] += data.arr.Lx % proc_num;
    //send one symbol from each side of sent array
    if (proc_num != 1) {
        //add to end of array in first process one section
        sendsections[0] += 1;
        //and to begin of array in last process 1 element
        displs[proc_num - 1] -= 1;
        sendsections[proc_num - 1] += 1;
        //add for 1 element to begin and end other processes
        for (int i = 1; i < proc_num - 1; i++) {
            displs[i] -= 1;
            sendsections[i] += 2;
        }
        for (int i = 0; i < proc_num; i++) {
            sendcounts[i] = sendsections[i] * data.arr.Lx * data.arr.Ly;
            displs[i] *= data.arr.Lx * data.arr.Ly;
        }
    }
    //ENV PREPARATION FOR OpenMP
    omp_set_num_threads(NUMDER_OF_THREADS);
    double time = 0;
    double out_time;
    if (!rank) {
        //for (int i = 0; i < data.arr.Lx * data.arr.Ly * data.arr.Lz; i++)
        //    printf("%f\n", data.arr.values[i]);
        printf("%d %d %d %d %d\n", sendsections[0], sendsections[1], sendsections[2], sendsections[3],
            sendsections[0] + sendsections[1] + sendsections[2] + sendsections[3]);
        printf("%d %d %d %d %d\n", sendcounts[0], sendcounts[1], sendcounts[2], sendcounts[3],
                                   sendcounts[0] + sendcounts[1] + sendcounts[2] + sendcounts[3]);
        printf("%d %d %d %d\n", displs[0], displs[1], displs[2], displs[3]);
        printf("%d\n", data.arr.Lx * data.arr.Ly * data.arr.Lz);
        printf("%d %d \n", elements_per_process, sections_per_process);
    }
    MPI_Status Status;
    MPI_Scatterv(data.arr.values, sendcounts, displs,
        MPI_DOUBLE, points.values, sendcounts[rank],
        MPI_DOUBLE,
        0, MPI_COMM_WORLD);
    //START CALCULATING

    for (count = 0; count <= maxCount; count++) {
        if (!rank) out_time = omp_get_wtime();
#if EILER
#pragma omp parallel for
        for (unsigned i = 0; i < sendcounts[rank]; i++)
            prev.values[i] = points.values[i];//make pointers swap
#pragma omp parallel for
        for (unsigned i = 1; i < sendsections[rank] - 1; i++)
            for (unsigned j = 1; j < points.Ly - 1; j++)
                for (unsigned k = 1; k < points.Lz - 1; k++)
                    set(&points, i, j, k,
                        get(&prev, i, j, k) + data.Sigma*data.deltaT*(
                        (get(&prev, i - 1, j, k) - 2 * get(&prev, i, j, k) + get(&prev, i + 1, j, k)) * Xdivider +
                        (get(&prev, i, j - 1, k) - 2 * get(&prev, i, j, k) + get(&prev, i, j + 1, k)) * Ydivider +
                        (get(&prev, i, j, k - 1) - 2 * get(&prev, i, j, k) + get(&prev, i, j, k + 1)) * Zdivider));
    
#endif /*EILER*/

        // Sending and recieve border points
        SendBorderValues(&points, sendsections[rank], rank, proc_num);

        // delete output time from time calculating
        if (!rank) time += omp_get_wtime() - out_time;

        if (count%OutCount == 0 && count) {
            if (!rank) {
                //Can't send-receive from 0 process to himself
                for (int i = 0; i < elements_per_process; i++)
                    print_array.values[i] = points.values[i];
                for (int sender_number = 1; sender_number < proc_num; sender_number++)
                    MPI_Recv(&print_array.values[elements_per_process*sender_number],
                        /* If we Recv from last process - getting remainder(ostatok) of the division */
                        sender_number == proc_num - 1 ?
                        elements_per_process + data.arr.Lx % proc_num * data.arr.Ly * data.arr.Lz :
                        elements_per_process
                        /* In other cases getting only elements_per_process count */
                        , MPI_DOUBLE, sender_number, sender_number, MPI_COMM_WORLD, &Status);
                OutputArr(inputFile, print_array.values, data.arr.Lx*data.arr.Ly*data.arr.Lz);
                OutTime += data.deltaOut;
                OutputCurrentTime(inputFile, OutTime);
            }
            // Send points from not root process to root for Outpt
            else if (rank == proc_num - 1)
                MPI_Send( points.values + data.arr.Ly * data.arr.Lz,
                          elements_per_process + data.arr.Lx % proc_num * data.arr.Ly * data.arr.Lz, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
            else
                MPI_Send( points.values + data.arr.Ly * data.arr.Lz, elements_per_process, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
    return 0;
}