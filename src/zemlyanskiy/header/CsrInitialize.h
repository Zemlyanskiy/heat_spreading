#pragma once
#include "MpiWorking.h"
static inline void InitializeCsr(struct input_data data, int prog_mode)
{
    if (prog_mode == 2 || prog_mode == 3)
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
                        if (prog_mode == 2)
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
    else if (prog_mode == 4)
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
}