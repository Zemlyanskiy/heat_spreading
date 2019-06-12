#pragma once
#include "CsrInteraction.h"
#include "MpiWorking.h"
static inline void CsrEuler(struct input_data data, unsigned rank, unsigned proc_num)
{
    omp_set_num_threads(num_threads);
    start = omp_get_wtime();
    for (double time = data.t; time <= data.T; time += data.deltaT, count++)
    {
#pragma omp parallel for
        for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
            prev.values[i] = points.values[i];

        MultiCSR(prev.values, points.values, csr_values, first_line, columns, data.array.Lx * data.array.Ly * data.array.Lz);

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
