#pragma once
#include "MpiWorking.h"
static inline void Runge_Kutta(struct input_data data, unsigned rank, unsigned proc_num)
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
            Output(file, points.values, data.array.Lx * data.array.Ly * data.array.Lz);
            CurrentTime(file, CurTime);
            CurTime += data.deltaOut;
        }
    }

    end = omp_get_wtime();
    Restime = end - start;
    printf("Lead time in sec: \n");
    printf("%lf\n", Restime);
}