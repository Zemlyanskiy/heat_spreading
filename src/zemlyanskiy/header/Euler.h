#pragma once
#include "MpiWorking.h"
static inline void Euler(struct input_data data, unsigned rank, unsigned proc_num)
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
                for (int k = 1; k < Zlines - 1; k++)
                    SetValue(&points, i, j, k,
                             GetValue(&prev, i, j, k) + data.Sigma * data.deltaT * ((GetValue(&prev, i - 1, j, k) - 2 * GetValue(&prev, i, j, k) + GetValue(&prev, i + 1, j, k)) * Xdivider + (GetValue(&prev, i, j - 1, k) - 2 * GetValue(&prev, i, j, k) + GetValue(&prev, i, j + 1, k)) * Ydivider + (GetValue(&prev, i, j, k - 1) - 2 * GetValue(&prev, i, j, k) + GetValue(&prev, i, j, k + 1)) * Zdivider));
        UpdateBorders(points.values, Xlines, Ylines, Zlines, rank, process_per_axis, Xcube_pos, Ycube_pos, Zcube_pos, proc_num, border_size);
        if (count % out_count == 0)
        {
            if (rank == 0)
            {
                ReceivePrint(print_array.values, points.values, Xlines, Ylines, Zlines,
                                         x_return_lines, y_return_lines, z_return_lines, x_return_start, y_return_start, z_return_start,
                                         buffer, file, data, proc_num);
                Output(file, points.values, data.array.Lx * data.array.Ly * data.array.Lz);
                CurrentTime(file, CurTime);
                CurTime += data.deltaOut;
            }
            else
            {
                SendPrint(points.values, Xlines, Ylines, Zlines, rank);
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
