#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main() {
    FILE *fp = fopen("input1.txt", "r");

    //Input data
    double XMin, XMax;
    double Sigma;
    int Lx;
    double t,T;
    double deltaT;
    double deltaOut;

    //Calculate data
    double step;
    double *points;
    double *prev;
    double posX;
    int numbers_of_out;
    int out_count;
    int count;

    //Reading
    if (fp == NULL) {
        printf("File can`t open!\n");
        exit(-1);
    }

    fscanf(fp, "t=%lf\n", &t);
    fscanf(fp, "T=%lf\n", &T);
    fscanf(fp, "deltaT=%lf\n", &deltaT);
    fscanf(fp, "XMin=%lf\n", &XMin);
    fscanf(fp, "XMax=%lf\n", &XMax);
    fscanf(fp, "Lx=%d\n", &Lx);
    fscanf(fp, "Sigma=%lf\n", &Sigma);
    fscanf(fp, "deltaOut=%lf\n", &deltaOut);

    //Calculate
    step = (XMax - XMin) / Lx;
    points = (double*)malloc(sizeof(double)*Lx);
    prev = (double*)malloc(sizeof(double)*Lx);
    posX = XMin;
    for (int i = 0; i < Lx; i++) {
      if ((posX >= -0.5) && (posX <= 0.5))
    		points[i] = cos(posX*3.141592);
    	else
    		points[i] = 0;
  		posX += step;
  	}
    numbers_of_out =(int) (T - t) / deltaOut;
  	out_count = (int)(T - t) / deltaT;
  	out_count = out_count / numbers_of_out;
  	count = 0;
    //Writing
    //First step
    fp = fopen("input1.txt","a");
    fprintf(fp, "\n");
    for (int i = 0; i < Lx; i++) {
      fprintf(fp, "%f ", points[i]);
      //printf("%lf\n", points[i]);
    }
    fprintf(fp, "\n");
    //Next steps
    for (double time = t; time <= T; time += deltaT, count++) {
  		for (int i = 0; i < Lx; i++)
  			prev[i] = points[i];
  		for (int i = 1; i < Lx - 1; i++) {
  			points[i] = prev[i] + Sigma*deltaT*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) / pow(step, 2);
        if (count%out_count == 0) {
          fprintf(fp, "%f ", points[i]);
        }
      }
    }
    fclose(fp);
    free(points);
    free(prev);
    printf("THE END\n");
    return 0;
}
