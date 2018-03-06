#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define EILER 0;
#define RUNGE_KUTTA 1;
#define START_STRING 2;

struct inputDATA {
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

struct inputDATA Read(char *DATA,int PROGRAMM_MODE) {
  FILE *file = fopen(DATA, "r");
  struct inputDATA temp;
  char str[5000];
  int pos;
  char *lexeme;
  double buffer;

  if(file == NULL) {
    printf("Can`t open file\n");
    exit(-1);
  } else {
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

    if(PROGRAMM_MODE == 2 && fgets(str, sizeof(str), file) == NULL)
      printf("Wrong line\n");
    else {
      pos = 0;
      lexeme = strtok(str, " ");
      while (lexeme != NULL) {
        buffer = strtod(lexeme, NULL);
        temp.values[pos] = buffer;
        lexeme = strtok(NULL, " ");
        pos++;
      }
    }
  }
  fclose(file);
  return temp;
}

int CurrentTime(char* DATA, double currentT) {
	FILE* file = fopen(DATA, "r");
	int ERR_CODE,num_of_space=15-8;
	int tmp=(int)currentT;
	if (file == NULL)
		printf("Cant open file\n");
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
    char* FILE = argv[1];
    int PROGRAMM_MODE;
    if(argv[2][0] == 'e')
      PROGRAMM_MODE = 0;
    if(argv[2][0] == 'r')
      PROGRAMM_MODE = 1;
    if(argv[2][0] == 's')
      PROGRAMM_MODE = 2;

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
    struct inputDATA DATA = Read(FILE, PROGRAMM_MODE);

    double *k1, *k2, *k3, *k4, *medium;
    double *derivative;

    step = (DATA.XMax - DATA.XMin) / DATA.Lx;
    points = (double*)malloc(sizeof(double)*DATA.Lx);
    posX = DATA.XMin;
    printf("%f\n",DATA.t);
    printf("%f\n",DATA.T);
    printf("%f\n",DATA.deltaT);
    printf("%f\n",DATA.XMin);
    printf("%f\n",DATA.XMax);
    printf("%i\n",DATA.Lx);
    printf("%f\n",DATA.Sigma);
    printf("%f\n",DATA.deltaOut);
    //Calculate data
    if(PROGRAMM_MODE == 2) {
      for (int i = 0; i < DATA.Lx; i++) {
        if ((posX >= -0.5) && (posX <= 0.5))
    		  points[i] = cos(posX*3.141592);
    	     else
    		     points[i] = 0;
  		         posX += step;
  	    }
    } else {
        for(int i = 0; i < DATA.Lx; i++) {
          points[i] = DATA.values[i];
        }
        //Calculate data
        prev = (double*)malloc(sizeof(double)*DATA.Lx);
        numbers_of_out =(int) (DATA.T - DATA.t) / DATA.deltaOut;
      	out_count = (int)(DATA.T - DATA.t) / DATA.deltaT;
      	out_count = out_count / numbers_of_out;
      	count = 0;
        if(PROGRAMM_MODE == 0) {//Euler
          for (double time = DATA.t; time <= DATA.T; time += DATA.deltaT, count++) {
            for (int i = 0; i < DATA.Lx; i++)
              prev[i] = points[i];
            for (int i = 1; i < DATA.Lx - 1; i++)
              points[i] = prev[i] + DATA.Sigma*DATA.deltaT*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) / pow(step, 2);
            if (count%out_count == 0) {
              Output(FILE, "a", points, DATA.Lx);
              CurrentTime(FILE, OutTime);
              OutTime += DATA.deltaOut;
            }
          }
        }
        if(PROGRAMM_MODE == 1) {//Runge_Kutta
          k1 = (double*)malloc(sizeof(double)*DATA.Lx);
    			k2 = (double*)malloc(sizeof(double)*DATA.Lx);
    			k3 = (double*)malloc(sizeof(double)*DATA.Lx);
    			k4 = (double*)malloc(sizeof(double)*DATA.Lx);
    			medium = (double*)malloc(sizeof(double)*DATA.Lx);
    			derivative = (double*)malloc(sizeof(double)*DATA.Lx);
    			double* eiler_derivative = (double*)malloc(sizeof(double)*DATA.Lx);
          for (int i = 0; i < DATA.Lx; i++) {
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

    //Writing
    //First step
    //fprintf(fp,"\n");
    //Output("input1.txt", "a", points, DATA.Lx);
    //Next steps

    //fclose(fp);
    free(points);
    free(prev);
    printf("THE END\n");
    return 0;
}
