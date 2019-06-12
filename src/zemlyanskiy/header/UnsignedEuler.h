#pragma once
static inline void UnsignedEuler(struct input_data data, unsigned rank, unsigned proc_num)
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
            MultiCSR(prev.values, difference, csr_values, first_line, columns, data.array.Lx * data.array.Ly * data.array.Lz);

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