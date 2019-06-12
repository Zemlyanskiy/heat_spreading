#pragma once
struct input_data
{
  //Input data
  double t;        //Start time
  double T;        //End time
  double deltaT;   //Gap between measurements
  double XMin;     //Start OX
  double XMax;     //End OX
  int Lx;          //Number of measurements
  double Sigma;    //Coefficient of thermal conductivity
  double deltaOut; //Delta of Output
  double *values;  //Last array with values
};

struct input_data Read(char *DATA, int PROGRAMM_MODE)
{
  FILE *file = fopen(DATA, "r");
  struct input_data temp;
  char str[5000];
  int pos;
  char *lexeme;
  double buffer;

  if (file == NULL)
  {
    printf("Can`t open file\n");
    exit(-1);
  }
  else
  {
    printf("Read input file:\n");
    for (int i = 0; fgets(str, sizeof(str), file) != NULL; i++)
    {
      fscanf(file, "t=%lf\n", &temp.t);
      fscanf(file, "T=%lf\n", &temp.T);
      fscanf(file, "deltaT=%lf\n", &temp.deltaT);
      fscanf(file, "XMin=%lf\n", &temp.XMin);
      fscanf(file, "XMax=%lf\n", &temp.XMax);
      fscanf(file, "Lx=%d\n", &temp.Lx);
      fscanf(file, "Sigma=%lf\n", &temp.Sigma);
      fscanf(file, "deltaOut=%lf\n", &temp.deltaOut);

      temp.values = malloc(sizeof(double) * temp.Lx);
      for (int i = 0; i < temp.Lx; i++)
        temp.values[i] = 0;

      if (PROGRAMM_MODE == 6)
        printf("Wrong!!!\n");
      else
      {
        pos = 0;
        lexeme = strtok(str, " ");
        while (lexeme != NULL)
        {
          buffer = strtof(lexeme, NULL);
          temp.values[pos] = buffer;
          lexeme = strtok(NULL, " ");
          pos++;
        }
      }
    }
  }
  printf("%f\n", temp.t);
  printf("%f\n", temp.T);
  printf("%f\n", temp.deltaT);
  printf("%f\n", temp.XMin);
  printf("%f\n", temp.XMax);
  printf("%i\n", temp.Lx);
  printf("%f\n", temp.Sigma);
  printf("%f\n", temp.deltaOut);
  for (int i = 0; i < temp.Lx; i++)
    printf("%lf ", temp.values[i]);
  printf("\n");
  fclose(file);
  return temp;
}

int CurrentTime(char *DATA, double currentT)
{
  FILE *file = fopen(DATA, "r+");
  int ERR_CODE, num_of_space = 15 - 8;
  int tmp = (int)currentT;
  if (file == NULL)
    printf("Can't open file\n");
  else
  {
    while (tmp > 10)
    {
      tmp /= 10;
      num_of_space--;
    }
    fseek(file, 0, SEEK_SET);
    fputs("t=", file);
    for (int i = 0; i < num_of_space; i++)
      fputs(" ", file);
    fprintf(file, "%f", currentT);
    ERR_CODE = ftell(file);
    fclose(file);
  }
  return 0;
}

int Output(char *DATA, char *mode, double *arr, int length)
{
  FILE *file = fopen(DATA, mode);
  if (file == NULL)
    printf("Can`t open file\n");
  else
  {
    for (int i = 0; i < length; i++)
      fprintf(file, "%f ", arr[i]);
    fprintf(file, "\n");
    fclose(file);
  }
}

double *MultiCSR(double *arr, double *csr_values, int *first_line, int *columns, struct input_data temp)
{
  double *composition = (double *)malloc(sizeof(double) * temp.Lx);
#pragma omp parallel for
  for (int i = 0; i < temp.Lx; i++)
  {
    composition[i] = 0;
    for (int j = first_line[i]; j < first_line[i + 1]; j++)
      composition[i] += csr_values[j] * arr[columns[j]];
  }
  return composition;
}

double GetCsrDiag(double *csr_values, int *columns, int *first_line, int iRaw, int jRaw, int Lx)
{
  if (iRaw >= Lx)
    return 0;
  int firstNum = first_line[iRaw];
  int secondNum = first_line[iRaw + 1];
  for (int j = firstNum; j < secondNum; j++)
  {
    if (columns[j] == jRaw)
    {
      return csr_values[j];
    }
  }

  return 0;
}