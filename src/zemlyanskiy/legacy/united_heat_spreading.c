#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include "../header/heat_spreading.h"

#define EILER 0;
#define RUNGE_KUTTA 1;
#define UNSIGNED_EILER 2;
#define START_STRING 3;

int main(int argc, char *argv[])
{
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

    int num_threads = atoi(argv[4]);
    //AD calculate data
    double step;
    double quad_step;
    double sum;
    double *cur_points;
    double *prev_points;
    double posX;
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

    //Reading data
    struct input_data data = Read(file, prog_mode);

    double *k1, *k2, *k3, *k4, *middle;

    step = (data.XMax - data.XMin) / data.Lx;
    quad_step = 1 / pow(step, 2);
    cur_points = (double *)malloc(sizeof(double) * data.Lx);
    posX = data.XMin;

    //Calculate data
    //Start_string
    if (prog_mode == 5)
    {
        for (int i = 0; i < data.Lx; i++)
        {
            if ((posX >= -0.5) && (posX <= 0.5))
                cur_points[i] = cos(posX * 3.141592);
            else
                cur_points[i] = 0;
            posX += step;
        }
        //First output result
        Output(file, "a", cur_points, data.Lx);
    }
    else
    {
        for (int i = 0; i < data.Lx; i++)
        {
            cur_points[i] = data.values[i];
        }
    }
    //Initialize CSR Matrix
    if (prog_mode == 0 || prog_mode == 1)
    {
        first_line = (int *)malloc(sizeof(int) * data.Lx + 1);
        columns = (int *)malloc(sizeof(int) * (data.Lx - 2) * 3);
        csr_values = (double *)malloc(sizeof(double) * (data.Lx - 2) * 3);
        for (int i = 0; i < data.Lx - 2; i++)
        {
            csr_values[i * 3] = data.Sigma * data.deltaT * quad_step;
            csr_values[i * 3 + 1] = 1 - 2 * data.Sigma * data.deltaT * quad_step;
            csr_values[i * 3 + 2] = data.Sigma * data.deltaT * quad_step;
            first_line[i + 1] = i * 3;
            columns[i * 3] = i;
            columns[i * 3 + 1] = i + 1;
            columns[i * 3 + 2] = i + 2;
        }
        first_line[0] = first_line[1];
        first_line[data.Lx - 1] = first_line[data.Lx - 2];
        first_line[data.Lx] = (data.Lx - 2) * 3;
    }
    else if (prog_mode == 2)
    {
        first_line = (int *)malloc(sizeof(int) * data.Lx + 1);
        columns = (int *)malloc(sizeof(int) * (data.Lx - 2) * 3 + 2);
        csr_values = (double *)malloc(sizeof(double) * (data.Lx - 2) * 3 + 2);
        for (int i = 0; i < data.Lx - 2; i++)
        {
            csr_values[i * 3] = -data.Sigma * data.deltaT * quad_step;
            csr_values[i * 3 + 1] = 1 + 2 * data.Sigma * data.deltaT * quad_step;
            csr_values[i * 3 + 2] = -data.Sigma * data.deltaT * quad_step;
            first_line[i + 1] = i * 3 + 1;
            columns[i * 3] = i;
            columns[i * 3 + 1] = i + 1;
            columns[i * 3 + 2] = i + 2;
        }
        first_line[0] = first_line[1];
        first_line[data.Lx - 1] = first_line[data.Lx - 2];
        first_line[data.Lx] = (data.Lx - 2) * 3;

        delta = 0.0001;
        difference = (double *)malloc(sizeof(double) * data.Lx);
        flag = false;

        CsrDiag = (double *)malloc(sizeof(double) * data.Lx);
        for (int i = 0; i < data.Lx; i++)
        {
            CsrDiag[i] = 1 / GetCsrDiag(csr_values, columns, first_line, i, i, data.Lx);
        }
    }

    //Calculate data
    prev_points = (double *)malloc(sizeof(double) * data.Lx);
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
#pragma omp parallel
            {
#pragma omp for
                for (int i = 0; i < data.Lx; i++)
                    prev_points[i] = cur_points[i];

                cur_points = MultiCSR(prev_points, csr_values, first_line, columns, data);
            }
            if (count % out_count == 0)
            {
                Output(file, "a", cur_points, data.Lx);
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
        k1 = (double *)malloc(sizeof(double) * data.Lx);
        k2 = (double *)malloc(sizeof(double) * data.Lx);
        k3 = (double *)malloc(sizeof(double) * data.Lx);
        k4 = (double *)malloc(sizeof(double) * data.Lx);
        middle = (double *)malloc(sizeof(double) * data.Lx);
        for (int i = 0; i < data.Lx; i++)
        {
            k1[i] = 0;
            k2[i] = 0;
            k3[i] = 0;
            k4[i] = 0;
            middle[i] = 0;
        }
        omp_set_num_threads(num_threads);
        start = omp_get_wtime();
        for (double time = data.t; time <= data.T; time += data.deltaT, count++)
        {
#pragma omp parallel
            {
#pragma omp for
                for (int i = 0; i < data.Lx; i++)
                    prev_points[i] = cur_points[i];

                k1 = MultiCSR(prev_points, csr_values, first_line, columns, data);

#pragma omp for
                for (int i = 0; i < data.Lx; i++)
                    middle[i] = (prev_points[i] + k1[i]) * 0.5;

                k2 = MultiCSR(middle, csr_values, first_line, columns, data);
#pragma omp for
                for (int i = 0; i < data.Lx; i++)
                    middle[i] = (prev_points[i] + k2[i]) * 0.5;

                k3 = MultiCSR(middle, csr_values, first_line, columns, data);
#pragma omp for
                for (int i = 0; i < data.Lx; i++)
                    middle[i] = prev_points[i];

                k4 = MultiCSR(middle, csr_values, first_line, columns, data);

#pragma omp for
                for (int i = 1; i < data.Lx - 1; i++)
                    middle[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * 0.166666667;
// Calculate next steps
#pragma omp for
                for (int i = 1; i < data.Lx - 1; i++)
                    cur_points[i] = middle[i];
            }
            if (count % out_count == 0)
            {
                Output(file, "a", cur_points, data.Lx);
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
            for (int i = 0; i < data.Lx; i++)
            {
                prev_points[i] = cur_points[i];
            }
            while (flag)
            {
                flag = false;
                //checking difference by delta
                difference = MultiCSR(prev_points, csr_values, first_line, columns, data);
#pragma omp parallel for
                for (int i = 0; i < data.Lx; i++)
                {
                    difference[i] -= prev_points[i];
                    difference[i] = fabs(difference[i]);
                    if (difference[i] >= delta)
                    {
                        flag = true;
                    }
                }
                if (flag)
                    break;

#pragma omp parallel for
                for (int i = 0; i < data.Lx; i++)
                    difference[i] = cur_points[i];

                for (int i = 0; i < data.Lx; i++)
                {
                    sum = 0;
                    for (int j = 0; j < data.Lx; j++)
                    {
                        if (j != i)
                            sum += GetCsrDiag(csr_values, columns, first_line, j, j, data.Lx) * difference[j];
                    }
                    cur_points[i] = CsrDiag[i] * (prev_points[i] - sum);
                }
            }

            if (count % out_count == 0)
            {
                Output(file, "a", cur_points, data.Lx);
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
        start = omp_get_wtime();
        for (double time = data.t; time <= data.T; time += data.deltaT, count++)
        {
#pragma omp parallel
            {
#pragma omp for
                for (int i = 0; i < data.Lx; i++)
                    prev_points[i] = cur_points[i];
#pragma omp for
                for (int i = 1; i < data.Lx - 1; i++)
                    cur_points[i] = prev_points[i] + data.Sigma * data.deltaT * (prev_points[i - 1] - 2 * prev_points[i] + prev_points[i + 1]) * quad_step;
            }
            if (count % out_count == 0)
            {
                Output(file, "a", cur_points, data.Lx);
                CurrentTime(file, CurTime);
                CurTime += data.deltaOut;
            }
        }
        end = omp_get_wtime();
        Restime = end - start;
        printf("Lead time in sec: \n");
        printf("%lf\n", Restime);
    }
    if (prog_mode == 4)
    {
        k1 = (double *)malloc(sizeof(double) * data.Lx);
        k2 = (double *)malloc(sizeof(double) * data.Lx);
        k3 = (double *)malloc(sizeof(double) * data.Lx);
        k4 = (double *)malloc(sizeof(double) * data.Lx);
        middle = (double *)malloc(sizeof(double) * data.Lx);
        for (int i = 0; i < data.Lx; i++)
        {
            k1[i] = 0;
            k2[i] = 0;
            k3[i] = 0;
            k4[i] = 0;
            middle[i] = 0;
        }
        omp_set_num_threads(num_threads);
        start = omp_get_wtime();
        for (double time = data.t; time <= data.T; time += data.deltaT, count++)
        {
#pragma omp parallel
            {
#pragma omp for
                for (int i = 0; i < data.Lx; i++)
                    prev_points[i] = cur_points[i];
#pragma omp for
                for (int i = 1; i < data.Lx - 1; i++)
                {
                    k1[i] = data.Sigma * data.deltaT * (prev_points[i - 1] - 2 * prev_points[i] + prev_points[i + 1]) * quad_step;

                    k2[i] = data.Sigma * data.deltaT * ((prev_points[i - 1] + k1[i - 1] * 0.5) - 2 * (prev_points[i] + k1[i] * 0.5) + (prev_points[i + 1] + k1[i + 1] * 0.5)) * quad_step;

                    k3[i] = data.Sigma * data.deltaT * ((prev_points[i - 1] + k2[i - 1] * 0.5) - 2 * (prev_points[i] + k2[i] * 0.5) + (prev_points[i + 1] + k2[i + 1] * 0.5)) * quad_step;

                    k4[i] = data.Sigma * data.deltaT * ((prev_points[i - 1] + k3[i - 1]) - 2 * (prev_points[i] + k3[i]) + (prev_points[i + 1] + k3[i + 1])) * quad_step;

                    middle[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * 0.166666667;
                }
// Calculate next steps
#pragma omp for
                for (int i = 0; i < data.Lx; i++)
                    cur_points[i] = prev_points[i] + middle[i];
            }
            if (count % out_count == 0)
            {
                Output(file, "a", cur_points, data.Lx);
                CurrentTime(file, CurTime);
                CurTime += data.deltaOut;
            }
        }
        end = omp_get_wtime();
        Restime = end - start;
        printf("Lead time in sec: \n");
        printf("%lf\n", Restime);
    }
    free(cur_points);
    free(prev_points);
    printf("THE END\n");
    return 0;
}
