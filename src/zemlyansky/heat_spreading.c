#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define EILER 0;
#define RUNGE_KUTTA 1;
#define START_STRING 2;

struct input {
  //Input data
  double t;//Start time
  double T;//End time
  double deltaT;//Gap between measurements
  double XMin;//Start OX
  double XMax;//End OX
  int Lx;//Number of measurements
  double Sigma;//Coefficient of thermal conductivity
  double deltaOut;//Delta of Output
  double *values;//Last array with values
};

struct input Read(char *DATA,int PROGRAMM_MODE) {
  FILE *file = fopen(DATA, "r");
  struct input temp;
  char str[5000];
  int pos;
  char *lexeme;
  double buffer;

  if(file == NULL) {
    printf("Can`t open file\n");
    exit(-1);
  } else {
    printf("Read input file:\n");
    for (int i = 0; fgets(str, sizeof(str), file) != NULL; i++) {
    fscanf(file, "t=%lf\n", &temp.t);
    fscanf(file, "T=%lf\n", &temp.T);
    fscanf(file, "deltaT=%lf\n", &temp.deltaT);
    fscanf(file, "XMin=%lf\n", &temp.XMin);
    fscanf(file, "XMax=%lf\n", &temp.XMax);
    fscanf(file, "Lx=%d\n", &temp.Lx);
    fscanf(file, "Sigma=%lf\n", &temp.Sigma);
    fscanf(file, "deltaOut=%lf\n", &temp.deltaOut);

    temp.values = malloc(sizeof(double)*temp.Lx);
    for(int i = 0; i < temp.Lx; i++)
      temp.values[i] = 0;

    if(PROGRAMM_MODE == 2)
      printf("Wrong line\n");
    else {
      pos = 0;
      lexeme = strtok(str, " ");
      while (lexeme != NULL) {
        buffer = strtof(lexeme, NULL);
        temp.values[pos] = buffer;
        lexeme = strtok(NULL, " ");
        pos++;
      }
    }
  }
}
printf("%f\n",temp.t);
printf("%f\n",temp.T);
printf("%f\n",temp.deltaT);
printf("%f\n",temp.XMin);
printf("%f\n",temp.XMax);
printf("%i\n",temp.Lx);
printf("%f\n",temp.Sigma);
printf("%f\n",temp.deltaOut);
  for (int i = 0; i < temp.Lx; i++)
		printf("%lf ", temp.values[i]);
	printf("\n");
  fclose(file);
  return temp;
}

int CurrentTime(char* DATA, double currentT) {
	FILE* file = fopen(DATA, "r+");
	int ERR_CODE,num_of_space=15-8;
	int tmp=(int)currentT;
	if (file == NULL)
		printf("Can't open file\n");
	else {
		while (tmp > 10) {
			tmp /= 10;
			num_of_space--;
		}
		fseek(file, 0, SEEK_SET);
		fputs("t=", file);
		for (int i=0;i<num_of_space;i++)
			fputs(" ", file);
		fprintf(file,"%f",currentT);
		ERR_CODE = ftell(file);
		fclose(file);
	}
  return 0;
}

int Output(char* DATA,char* mode,double* arr,int length) {
  FILE* file = fopen(DATA, mode);
  if (file == NULL)
		printf("Can`t open file\n");
	else {
		for (int i = 0; i < length; i++)
			fprintf(file, "%f ", arr[i]);
		fprintf(file, "\n");
		fclose(file);
	}
}

int main(int argc, char *argv[]) {
    //FILE *fp = fopen("input1.txt", "r");
    char* file = argv[1];
    int prog_mode;
    if(argv[2][0] == 'e')
      prog_mode = 0;
    if(argv[2][0] == 'r')
      prog_mode = 1;
    if(argv[2][0] == 's')
      prog_mode = 2;

    //AD calculate data
    double step;
    double *points;
    double *prev;
    double posX;
    int numbers_of_out;
    int out_count;
    int count;
    double OutTime = 0;

    //Reading data
    struct input data = Read(file, prog_mode);

    double *k1, *k2, *k3, *k4, *medium;
    double *derivative;

    step = (data.XMax - data.XMin) / data.Lx;
    points = (double*)malloc(sizeof(double)*data.Lx);
    posX = data.XMin;

    //Calculate data
    if(prog_mode == 2) {//Start_string
      for (int i = 0; i < data.Lx; i++) {
        if ((posX >= -0.5) && (posX <= 0.5))
    		  points[i] = cos(posX*3.141592);
    	  else
    		  points[i] = 0;
  		  posX += step;
  	    }
      //First output result
      Output(file, "a", points, data.Lx);
    } else {
        for(int i = 0; i < data.Lx; i++) {
          points[i] = data.values[i];
        }
        //Calculate data
        prev = (double*)malloc(sizeof(double)*data.Lx);
        numbers_of_out =(int) (data.T - data.t) / data.deltaOut;
      	out_count = (int)(data.T - data.t) / data.deltaT;
      	out_count = out_count / numbers_of_out;
      	count = 0;
        if(prog_mode == 0) {//Euler
          for (double time = data.t; time <= data.T; time += data.deltaT, count++) {
            for (int i = 0; i < data.Lx; i++)
              prev[i] = points[i];
            for (int i = 1; i < data.Lx - 1; i++)
              points[i] = prev[i] + data.Sigma*data.deltaT*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) / pow(step, 2);
            if (count%out_count == 0) {
              Output(file, "a", points, data.Lx);
              CurrentTime(file, OutTime);
              OutTime += data.deltaOut;
            }
          }
        }
        if(prog_mode == 1) {//Runge_Kutta
          k1 = (double*)malloc(sizeof(double)*data.Lx);
    			k2 = (double*)malloc(sizeof(double)*data.Lx);
    			k3 = (double*)malloc(sizeof(double)*data.Lx);
    			k4 = (double*)malloc(sizeof(double)*data.Lx);
    			medium = (double*)malloc(sizeof(double)*data.Lx);
    			derivative = (double*)malloc(sizeof(double)*data.Lx);
    			double* eiler_derivative = (double*)malloc(sizeof(double)*data.Lx);
          for (int i = 0; i < data.Lx; i++) {
    				k1[i] = 0;
    				k2[i] = 0;
    				k3[i] = 0;
    				k4[i] = 0;
    				medium[i] = 0;
    				derivative[i] = 0;
    				eiler_derivative[i] = 0;
    			}
        }
      }

    free(points);
    free(prev);
    printf("THE END\n");
    return 0;
}
