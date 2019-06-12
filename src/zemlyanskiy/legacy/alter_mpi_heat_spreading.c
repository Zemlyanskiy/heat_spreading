#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <mpi.h>
#include "../header/3d_heat_spreading.h"
#include "../header/mpi_check.h"

#define EILER 0;
#define RUNGE_KUTTA 1;
#define UNSIGNED_EILER 2;
#define START_STRING 3;

int main(int argc, char *argv[])
{
    // MPI initialization ---
    unsigned rank, proc_num;

    char *file = argv[1];
    int prog_mode;
    if (argv[3][0] == 'e' && argv[2][0] == 'm')
        prog_mode = 0;
    if (argv[3][0] == 'r' && argv[2][0] == 'm')
        prog_mode = 1;
    if (argv[3][0] == 'u' && argv[2][0] == 'm')
        prog_mode = 2;
    if (argv[3][0] == 'e' && argv[2][0] == 'n')
        prog_mode = 3;
    if (argv[3][0] == 'r' && argv[2][0] == 'n')
        prog_mode = 4;
    if (argv[3][0] == 's')
        prog_mode = 5;
    //Reading data
    struct input_data data = Read(file, prog_mode);
    int num_threads = atoi(argv[4]);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    if ((sqrt(proc_num) - (unsigned)sqrt(proc_num)) != 0)
    {
        printf("The number of processes must be able to proccess sqrt() without remainder.\n");
        MPI_Finalize();
        return 1;
    }
    //AD calculate data
    double step;
    double quad_step;
    double sum;
    double posX;
    double posY;
    double posZ;
    int numbers_of_out;
    int out_count;
    int max_count;
    int count;
    double CurTime = 0;
    double start, end;
    double Restime;
    double delta;
    //variables for CSR
    int *first_line;
    int *columns;
    double *csr_values;
    double *difference;
    double *CsrDiag;
    bool flag;
    int result;
    int border_size = 1;
    initialize_common_variables(data, rank, proc_num, border_size);
    struct values3d points;      // array for caclulation results
    struct values3d prev;        // array for previos state
    struct values3d print_array; // array for previos state
    InitValues3d(&points, Xlines, Ylines, data.array.Lz);
    InitValues3d(&prev, Xlines, Ylines, data.array.Lz);
    InitValues3d(&print_array, data.array.Lx, data.array.Ly, data.array.Lz);

    double Xstep = (data.XMax - data.XMin) / data.array.Lx;
    double Ystep = (data.YMax - data.YMin) / data.array.Ly;
    double Zstep = (data.ZMax - data.ZMin) / data.array.Lz;

    double Xdivider = 1 / (Xstep * Xstep);
    double Ydivider = 1 / (Ystep * Ystep);
    double Zdivider = 1 / (Zstep * Zstep);

    double *k1, *k2, *k3, *k4, *middle;
    step = (data.XMax - data.XMin) / data.array.Lx;
    quad_step = 1 / pow(step, 2);
    //Calculate data
    //Start_string
    if (prog_mode == 5)
    {
        posX = data.XMin;
        posY = data.YMin;
        posZ = data.ZMin;
        for (int i = 0; i < data.array.Lx; i++)
        {
            posX += Xstep;
            posY = data.YMin;
            for (int j = 0; j < data.array.Ly; j++)
            {
                posY += Ystep;
                posZ = data.ZMin;
                for (int k = 0; k < data.array.Lz; k++)
                {
                    posZ += Zstep;
                    SetValue(&points, i, j, k, CalcFunc(posX, posY, posZ));
                }
            }
        }
        //First output result
        Output(file, "a", points.values, data.array.Lx * data.array.Ly * data.array.Lz);
        MPI_Finalize();
    }
    else
    {
        for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
        {
            points.values[i] = data.array.values[i];
        }
    }
    //Initialize CSR Matrix
    if (prog_mode == 0 || prog_mode == 1)
    {
        unsigned borders_number = data.array.Lx * data.array.Ly * 2 + data.array.Ly * data.array.Lz * 2 + data.array.Lx * data.array.Lz * 2 - 8 - data.array.Lx * 4 - data.array.Ly * 4 - data.array.Lz * 4;
        unsigned values_number = (data.array.Lx * data.array.Ly * data.array.Lz - borders_number) * 7 + borders_number;
        first_line = (unsigned *)calloc(data.array.Lx * data.array.Ly * data.array.Lz + 1, sizeof(unsigned));
        columns = (unsigned *)calloc(values_number, sizeof(unsigned));
        csr_values = (double *)calloc(values_number, sizeof(double));
        unsigned counter = 0;
        for (unsigned i = 0; i < data.array.Lx; i++)
            for (unsigned j = 0; j < data.array.Ly; j++)
                for (unsigned k = 0; k < data.array.Lz; k++)
                {
                    first_line[i * data.array.Ly * data.array.Lz + j * data.array.Lz + k] = counter;

                    if (i == 0 || j == 0 || k == 0 || i == data.array.Lx - 1 || j == data.array.Ly - 1 || k == data.array.Lz - 1)
                    {
                        columns[counter] = i * data.array.Ly * data.array.Lz + j * data.array.Lz + k;
                        csr_values[counter++] = 1;
                    }
                    else
                    {
                        columns[counter] = (i - 1) * data.array.Ly * data.array.Lz + j * data.array.Lz + k;
                        csr_values[counter++] = data.Sigma * data.deltaT * Xdivider;
                        columns[counter] = i * data.array.Ly * data.array.Lz + (j - 1) * data.array.Lz + k;
                        csr_values[counter++] = data.Sigma * data.deltaT * Ydivider;
                        columns[counter] = i * data.array.Ly * data.array.Lz + j * data.array.Lz + k - 1;
                        csr_values[counter++] = data.Sigma * data.deltaT * Zdivider;

                        columns[counter] = i * data.array.Ly * data.array.Lz + j * data.array.Lz + k;
                        if (prog_mode == 0)
                            csr_values[counter++] = 1 - 2 * data.Sigma * data.deltaT * (Xdivider + Ydivider + Zdivider);
                        else
                            csr_values[counter++] = -2 * data.Sigma * data.deltaT * (Xdivider + Ydivider + Zdivider);

                        columns[counter] = i * data.array.Ly * data.array.Lz + j * data.array.Lz + k + 1;
                        csr_values[counter++] = data.Sigma * data.deltaT * Zdivider;
                        columns[counter] = i * data.array.Ly * data.array.Lz + (j + 1) * data.array.Lz + k;
                        csr_values[counter++] = data.Sigma * data.deltaT * Ydivider;
                        columns[counter] = (i + 1) * data.array.Ly * data.array.Lz + j * data.array.Lz + k;
                        csr_values[counter++] = data.Sigma * data.deltaT * Xdivider;
                    }
                }
        first_line[data.array.Lx * data.array.Ly * data.array.Lz] = counter;
    }
    else if (prog_mode == 2)
    {
        difference = (double *)calloc(data.array.Lx * data.array.Ly * data.array.Lz, sizeof(double));
        unsigned borders_number = data.array.Lx * data.array.Ly * data.array.Lz - (data.array.Lx - 2) * (data.array.Ly - 2) * (data.array.Lz - 2);
        unsigned values_number = (data.array.Lx * data.array.Ly * data.array.Lz - borders_number) * 7 + borders_number;
        csr_values = (double *)calloc(values_number, sizeof(double));
        columns = (unsigned *)calloc(values_number, sizeof(unsigned));
        first_line = (unsigned *)calloc(data.array.Lx * data.array.Ly * data.array.Lz + 1, sizeof(unsigned));
        CsrDiag = (double *)calloc(data.array.Lx * data.array.Ly * data.array.Lz, sizeof(double));

        unsigned counter = 0;
        for (unsigned i = 0; i < data.array.Lx; i++)
            for (unsigned j = 0; j < data.array.Ly; j++)
                for (unsigned k = 0; k < data.array.Lz; k++)
                {

                    first_line[i * data.array.Ly * data.array.Lz + j * data.array.Lz + k] = counter;
                    if (i == 0 || j == 0 || k == 0 || i == data.array.Lx - 1 || j == data.array.Ly - 1 || k == data.array.Lz - 1)
                    {
                        columns[counter] = i * data.array.Ly * data.array.Lz + j * data.array.Lz + k;
                        CsrDiag[i * data.array.Ly * data.array.Lz + j * data.array.Lz + k] = 1;
                        csr_values[counter++] = 1;
                    }
                    else
                    {
                        // x y z central*3 z y x
                        columns[counter] = (i - 1) * data.array.Ly * data.array.Lz + j * data.array.Lz + k;
                        csr_values[counter++] = -data.Sigma * data.deltaT * Xdivider;
                        columns[counter] = i * data.array.Ly * data.array.Lz + (j - 1) * data.array.Lz + k;
                        csr_values[counter++] = -data.Sigma * data.deltaT * Ydivider;
                        columns[counter] = i * data.array.Ly * data.array.Lz + j * data.array.Lz + k - 1;
                        csr_values[counter++] = -data.Sigma * data.deltaT * Zdivider;

                        columns[counter] = i * data.array.Ly * data.array.Lz + j * data.array.Lz + k;
                        csr_values[counter] = 1 + 2 * data.Sigma * data.deltaT * (Xdivider + Ydivider + Zdivider);
                        CsrDiag[i * data.array.Ly * data.array.Lz + j * data.array.Lz + k] = 1 / csr_values[counter];
                        counter++;

                        columns[counter] = i * data.array.Ly * data.array.Lz + j * data.array.Lz + k + 1;
                        csr_values[counter++] = -data.Sigma * data.deltaT * Zdivider;
                        columns[counter] = i * data.array.Ly * data.array.Lz + (j + 1) * data.array.Lz + k;
                        csr_values[counter++] = -data.Sigma * data.deltaT * Ydivider;
                        columns[counter] = (i + 1) * data.array.Ly * data.array.Lz + j * data.array.Lz + k;
                        csr_values[counter++] = -data.Sigma * data.deltaT * Xdivider;
                    }
                }
        first_line[data.array.Lx * data.array.Ly * data.array.Lz] = counter;

        delta = 0.00001;
        flag = false;
    }

    //Calculate data
    numbers_of_out = (int)(data.T - data.t) / data.deltaOut;
    max_count = (int)(data.T - data.t) / data.deltaT;
    out_count = max_count / numbers_of_out;
    count = 0;
    //EulerCsr
    if (prog_mode == 0)
    {
        omp_set_num_threads(num_threads);
        start = omp_get_wtime();
        for (double time = data.t; time <= data.T; time += data.deltaT, count++)
        {
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
                prev.values[i] = points.values[i];

            points.values = MultiCSR(prev.values, csr_values, first_line, columns, data);

            if (count % out_count == 0)
            {
                Output(file, "a", points.values, data.array.Lx * data.array.Ly * data.array.Lz);
                CurrentTime(file, CurTime);
                CurTime += data.deltaOut;
            }
        }
        end = omp_get_wtime();
        Restime = end - start;
        printf("Lead time in sec: \n");
        printf("%lf\n", Restime);
    }
    //Runge_KuttaCsr
    if (prog_mode == 1)
    {
        struct values3d k1;
        struct values3d k2;
        struct values3d k3;
        struct values3d k4;
        struct values3d middle;

        InitValues3d(&k1, data.array.Lx, data.array.Ly, data.array.Lz);
        InitValues3d(&k2, data.array.Lx, data.array.Ly, data.array.Lz);
        InitValues3d(&k3, data.array.Lx, data.array.Ly, data.array.Lz);
        InitValues3d(&k4, data.array.Lx, data.array.Ly, data.array.Lz);
        InitValues3d(&middle, data.array.Lx, data.array.Ly, data.array.Lz);
        omp_set_num_threads(num_threads);
        start = omp_get_wtime();
        for (double time = data.t; time <= data.T; time += data.deltaT, count++)
        {
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
                prev.values[i] = points.values[i];

            k1.values = MultiCSR(prev.values, csr_values, first_line, columns, data);
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx; i++)
                for (int j = 0; j < data.array.Ly; j++)
                    for (int k = 0; k < data.array.Lz; k++)
                    {
                        SetValue(&middle, i, j, k, GetValue(&prev, i, j, k) + GetValue(&k1, i, j, k) * 0.5);
                    }
            k2.values = MultiCSR(middle.values, csr_values, first_line, columns, data);
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx; i++)
                for (int j = 0; j < data.array.Ly; j++)
                    for (int k = 0; k < data.array.Lz; k++)
                    {
                        SetValue(&middle, i, j, k, GetValue(&prev, i, j, k) + GetValue(&k2, i, j, k) * 0.5);
                    }
            k3.values = MultiCSR(middle.values, csr_values, first_line, columns, data);
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx; i++)
                for (int j = 0; j < data.array.Ly; j++)
                    for (int k = 0; k < data.array.Lz; k++)
                    {
                        SetValue(&middle, i, j, k, GetValue(&prev, i, j, k) + GetValue(&k3, i, j, k));
                    }
            k4.values = MultiCSR(middle.values, csr_values, first_line, columns, data);
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx; i++)
                for (int j = 0; j < data.array.Ly; j++)
                    for (int k = 0; k < data.array.Lz; k++)
                    {
                        SetValue(&middle, i, j, k, (GetValue(&k1, i, j, k) + 2 * GetValue(&k2, i, j, k) + 2 * GetValue(&k3, i, j, k) + GetValue(&k4, i, j, k)) * 0.166666667);
                    }
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx; i++)
                for (int j = 0; j < data.array.Ly; j++)
                    for (int k = 0; k < data.array.Lz; k++)
                    {
                        SetValue(&points, i, j, k, GetValue(&prev, i, j, k) + GetValue(&middle, i, j, k));
                    }
            if (count % out_count == 0)
            {
                Output(file, "a", points.values, data.array.Lx * data.array.Ly * data.array.Lz);
                CurrentTime(file, CurTime);
                CurTime += data.deltaOut;
            }
        }
        end = omp_get_wtime();
        Restime = end - start;
        printf("Lead time in sec: \n");
        printf("%lf\n", Restime);
    }
    //Unsigned Eiler
    if (prog_mode == 2)
    {
        omp_set_num_threads(num_threads);
        start = omp_get_wtime();
        for (double time = data.t; time <= data.T; time += data.deltaT, count++)
        {
            flag = true;
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
                prev.values[i] = points.values[i];

            while (flag)
            {
                flag = false;
                difference = MultiCSR(prev.values, csr_values, first_line, columns, data);

                for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
                {
                    difference[i] -= prev.values[i];
                    difference[i] = fabs(difference[i]);
                    if (difference[i] > delta)
                    {
                        flag = true;
                        break;
                    }
                }
                flag = false;

                for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
                    difference[i] = points.values[i];

                for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
                {
                    sum = 0;
                    for (int j = 0; j < data.array.Lx * data.array.Ly * data.array.Lz; j++)
                    {
                        if (i != j)
                            sum += GetCsrDiag(csr_values, columns, first_line, i, j, data.array.Lx * data.array.Ly * data.array.Lz) * difference[j];
                    }
                    points.values[i] = CsrDiag[i] * (prev.values[i] - sum);
                }
            }
            if (count % out_count == 0)
            {
                Output(file, "a", points.values, data.array.Lx * data.array.Ly * data.array.Lz);
                CurrentTime(file, CurTime);
                CurTime += data.deltaOut;
            }
        }
        end = omp_get_wtime();
        Restime = end - start;
        printf("Lead time in sec: \n");
        printf("%lf\n", Restime);
    }
    //Euler
    if (prog_mode == 3)
    {
        omp_set_num_threads(num_threads);
        if (rank == 0)
            start = omp_get_wtime();
        for (double time = data.t; time <= data.T; time += data.deltaT, count++)
        {
#pragma omp parallel for
            for (int i = 0; i < elements_per_process; i++)
                prev.values[i] = points.values[i];
#pragma omp parallel for
            for (int i = 1; i < Xlines - 1; i++)
                for (int j = 1; j < Ylines - 1; j++)
                    for (int k = 1; k < data.array.Lz - 1; k++)
                        SetValue(&points, i, j, k,
                                 GetValue(&prev, i, j, k) + data.Sigma * data.deltaT * ((GetValue(&prev, i - 1, j, k) - 2 * GetValue(&prev, i, j, k) + GetValue(&prev, i + 1, j, k)) * Xdivider + (GetValue(&prev, i, j - 1, k) - 2 * GetValue(&prev, i, j, k) + GetValue(&prev, i, j + 1, k)) * Ydivider + (GetValue(&prev, i, j, k - 1) - 2 * GetValue(&prev, i, j, k) + GetValue(&prev, i, j, k + 1)) * Zdivider));
            //SendingBorderValues(points.values, Xlines, Ylines, Zlines, rank, process_per_axis, Xcube_pos, Ycube_pos, Zcube_pos, proc_num, border_size);
            SendBorderValues(points.values, Xlines, Ylines, data.array.Lz, rank, process_per_axis, Xcube_pos, Ycube_pos, proc_num, border_size);
            if (count % out_count == 0)
            {
                if (rank == 0)
                {
                    receive_and_print_master(print_array.values, points.values, Xlines, Ylines,
                                             x_return_lines, y_return_lines, x_return_start, y_return_start,
                                             buffer, file, data, proc_num);
                    Output(file, "a", points.values, data.array.Lx * data.array.Ly * data.array.Lz);
                    CurrentTime(file, CurTime);
                    CurTime += data.deltaOut;
                }
                else
                {
                    send_to_print(points.values, Xlines, Ylines, data.array.Lz, rank);
                }
            }
        }
        if (!rank)
        {
            end = omp_get_wtime();
            Restime = end - start;
            printf("Lead time in sec: \n");
            printf("%lf\n", Restime);
        }
    }
    //Runge-Kutta
    if (prog_mode == 4)
    {
        struct values3d k1;
        struct values3d k2;
        struct values3d k3;
        struct values3d k4;
        struct values3d middle;

        InitValues3d(&k1, data.array.Lx, data.array.Ly, data.array.Lz);
        InitValues3d(&k2, data.array.Lx, data.array.Ly, data.array.Lz);
        InitValues3d(&k3, data.array.Lx, data.array.Ly, data.array.Lz);
        InitValues3d(&k4, data.array.Lx, data.array.Ly, data.array.Lz);
        InitValues3d(&middle, data.array.Lx, data.array.Ly, data.array.Lz);

        omp_set_num_threads(num_threads);
        start = omp_get_wtime();
        for (double time = data.t; time <= data.T; time += data.deltaT, count++)
        {
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
                prev.values[i] = points.values[i];
#pragma omp parallel for
            for (int i = 1; i < data.array.Lx - 1; i++)
            {
                for (int j = 1; j < data.array.Ly - 1; j++)
                {
                    for (int k = 1; k < data.array.Lz - 1; k++)
                    {
                        SetValue(&k1, i, j, k, data.Sigma * data.deltaT * ((GetValue(&prev, i - 1, j, k) - 2 * GetValue(&prev, i, j, k) + GetValue(&prev, i + 1, j, k)) * Xdivider + (GetValue(&prev, i, j - 1, k) - 2 * GetValue(&prev, i, j, k) + GetValue(&prev, i, j + 1, k)) * Ydivider + (GetValue(&prev, i, j, k - 1) - 2 * GetValue(&prev, i, j, k) + GetValue(&prev, i, j, k + 1)) * Zdivider));

                        SetValue(&k2, i, j, k, data.Sigma * data.deltaT * (((GetValue(&prev, i - 1, j, k) + GetValue(&k1, i - 1, j, k) * 0.5) - 2 * (GetValue(&prev, i, j, k) + GetValue(&k1, i, j, k) * 0.5) + (GetValue(&prev, i + 1, j, k) + GetValue(&k1, i + 1, j, k) * 0.5)) * Xdivider + ((GetValue(&prev, i, j - 1, k) + GetValue(&k1, i, j - 1, k) * 0.5) - 2 * (GetValue(&prev, i, j, k) + GetValue(&k1, i, j, k) * 0.5) + (GetValue(&prev, i, j + 1, k) + GetValue(&k1, i, j + 1, k) * 0.5)) * Ydivider + ((GetValue(&prev, i, j, k - 1) + GetValue(&k1, i, j, k - 1) * 0.5) - 2 * (GetValue(&prev, i, j, k) + GetValue(&k1, i, j, k) * 0.5) + (GetValue(&prev, i, j, k + 1) + GetValue(&k1, i, j, k + 1) * 0.5)) * Zdivider));

                        SetValue(&k3, i, j, k, data.Sigma * data.deltaT * (((GetValue(&prev, i - 1, j, k) + GetValue(&k2, i - 1, j, k) * 0.5) - 2 * (GetValue(&prev, i, j, k) + GetValue(&k2, i, j, k) * 0.5) + (GetValue(&prev, i + 1, j, k) + GetValue(&k2, i + 1, j, k) * 0.5)) * Xdivider + ((GetValue(&prev, i, j - 1, k) + GetValue(&k2, i, j - 1, k) * 0.5) - 2 * (GetValue(&prev, i, j, k) + GetValue(&k2, i, j, k) * 0.5) + (GetValue(&prev, i, j + 1, k) + GetValue(&k2, i, j + 1, k) * 0.5)) * Ydivider + ((GetValue(&prev, i, j, k - 1) + GetValue(&k2, i, j, k - 1) * 0.5) - 2 * (GetValue(&prev, i, j, k) + GetValue(&k2, i, j, k) * 0.5) + (GetValue(&prev, i, j, k + 1) + GetValue(&k2, i, j, k + 1) * 0.5)) * Zdivider));

                        SetValue(&k4, i, j, k, data.Sigma * data.deltaT * (((GetValue(&prev, i - 1, j, k) + GetValue(&k3, i - 1, j, k)) - 2 * (GetValue(&prev, i, j, k) + GetValue(&k3, i, j, k)) + (GetValue(&prev, i + 1, j, k) + GetValue(&k3, i + 1, j, k))) * Xdivider + ((GetValue(&prev, i, j - 1, k) + GetValue(&k3, i, j - 1, k)) - 2 * (GetValue(&prev, i, j, k) + GetValue(&k3, i, j, k)) + (GetValue(&prev, i, j + 1, k) + GetValue(&k3, i, j + 1, k))) * Ydivider + ((GetValue(&prev, i, j, k - 1) + GetValue(&k3, i, j, k - 1)) - 2 * (GetValue(&prev, i, j, k) + GetValue(&k3, i, j, k)) + (GetValue(&prev, i, j, k + 1) + GetValue(&k3, i, j, k + 1))) * Zdivider));

                        SetValue(&middle, i, j, k, (GetValue(&k1, i, j, k) + 2 * GetValue(&k2, i, j, k) + 2 * GetValue(&k3, i, j, k) + GetValue(&k4, i, j, k)) * 0.166666667);
                    }
                }
            }
#pragma omp parallel for
            for (int i = 0; i < data.array.Lx; i++)
                for (int j = 0; j < data.array.Ly; j++)
                    for (int k = 0; k < data.array.Lz; k++)
                    {
                        SetValue(&points, i, j, k, GetValue(&prev, i, j, k) + GetValue(&middle, i, j, k));
                    }
            if (count % out_count == 0)
            {
                Output(file, "a", points.values, data.array.Lx * data.array.Ly * data.array.Lz);
                CurrentTime(file, CurTime);
                CurTime += data.deltaOut;
            }
        }

        end = omp_get_wtime();
        Restime = end - start;
        printf("Lead time in sec: \n");
        printf("%lf\n", Restime);
    }
    printf("THE END\n");
    MPI_Finalize();
    return 0;
}
