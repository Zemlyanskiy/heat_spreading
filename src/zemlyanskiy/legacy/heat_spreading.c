#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "../header/heat_spreading.h"

#define EILER 0;
#define RUNGE_KUTTA 1;
#define START_STRING 2;

int main(int argc, char *argv[])
{
  char *file = argv[1];
  int prog_mode;
  if (argv[2][0] == 'e')
    prog_mode = 0;
  if (argv[2][0] == 'r')
    prog_mode = 1;
  if (argv[2][0] == 's')
    prog_mode = 2;
  int num_threads = atoi(argv[3]);
  //AD calculate data
  double step;
  double quad_step;
  double *cur_points;
  double *prev_points;
  double posX;
  int numbers_of_out;
  int out_count;
  int count;
  double CurTime = 0;
  double start, end;
  double Restime;

  //Reading data
  struct input_data data = Read(file, prog_mode);

  double *k1, *k2, *k3, *k4, *middle;

  step = (data.XMax - data.XMin) / data.Lx;
  quad_step = 1 / pow(step, 2);
  cur_points = (double *)malloc(sizeof(double) * data.Lx);
  posX = data.XMin;

  //Calculate data
  if (prog_mode == 2)
  { //Start_string
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
  //Calculate data
  prev_points = (double *)malloc(sizeof(double) * data.Lx);
  numbers_of_out = (int)(data.T - data.t) / data.deltaOut;
  out_count = (int)(data.T - data.t) / data.deltaT;
  out_count = out_count / numbers_of_out;
  count = 0;
  if (prog_mode == 0)
  { //Euler
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
  //Runge_Kutta
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
