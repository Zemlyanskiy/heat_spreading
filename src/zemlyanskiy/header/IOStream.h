#pragma once
//Number of measurements (x, y, z)
struct values3d
{
  int Lx;
  int Ly;
  int Lz;
  double *values;
};

//Input data
struct input_data
{
  double t;        //Start time
  double T;        //End time
  double deltaT;   //Gap between measurements
  double XMin;     //Start OX
  double XMax;     //End OX
  double YMin;     //Start OY
  double YMax;     //End OY
  double ZMin;     //Start OZ
  double ZMax;     //End OZ
  double Sigma;    //Coefficient of thermal conductivity
  double deltaOut; //Delta of Output
  struct values3d array;
};

static inline void InitValues3d(struct values3d *array, int x, int y, int z)
{
  array->Lx = x;
  array->Ly = y;
  array->Lz = z;
  array->values = (double *)calloc(x * y * z, sizeof(double));
}

static inline double GetValue(struct values3d *array, int x, int y, int z)
{
  return array->values[z + y * array->Ly + x * array->Ly * array->Lx];
}

static inline void SetValue(struct values3d *array, int x, int y, int z, double value)
{
  array->values[z + y * array->Ly + x * array->Ly * array->Lx] = value;
}

static inline struct input_data Read(char *DATA, int PROGRAMM_MODE)
{
  FILE *file = fopen(DATA, "r");
  struct input_data temp = {.t = 0, .T = 0, .deltaT = 0, .XMin = 0, .XMax = 0, .array.Lx = 0, .YMin = 0, .YMax = 0, .array.Ly = 0, .ZMin = 0, .ZMax = 0, .array.Lz = 0, .Sigma = 0, .deltaOut = 0, .array.values = 0};
  char *str = (char *)calloc(10000000, sizeof(char));
  int pos = 0;
  char *lexeme;
  double buffer = 0;

  if (file == NULL)
  {
    printf("Can`t open file\n");
    exit(-1);
  }
  else
  {
    printf("Read input file:\n");
    for (int i = 0; fgets(str, 10000000*sizeof(char), file) != NULL; i++)
    {
      switch (i)
      {
      case 0:
      {
        sscanf(str, "t=%lf\n", &temp.t);
        break;
      }
      case 1:
      {
        sscanf(str, "T=%lf\n", &temp.T);
        break;
      }
      case 2:
      {
        sscanf(str, "deltaT=%lf\n", &temp.deltaT);
        break;
      }
      case 3:
      {
        sscanf(str, "XMin=%lf\n", &temp.XMin);
        break;
      }
      case 4:
      {
        sscanf(str, "XMax=%lf\n", &temp.XMax);
        break;
      }
      case 5:
      {
        sscanf(str, "YMin=%lf\n", &temp.YMin);
        break;
      }
      case 6:
      {
        sscanf(str, "YMax=%lf\n", &temp.YMax);
        break;
      }
      case 7:
      {
        sscanf(str, "ZMin=%lf\n", &temp.ZMin);
        break;
      }
      case 8:
      {
        sscanf(str, "ZMax=%lf\n", &temp.ZMax);
        break;
      }
      case 9:
      {
        sscanf(str, "Lx=%d\n", &temp.array.Lx);
        break;
      }
      case 10:
      {
        sscanf(str, "Ly=%d\n", &temp.array.Ly);
        break;
      }
      case 11:
      {
        sscanf(str, "Lz=%d\n", &temp.array.Lz);
        break;
      }
      case 12:
      {
        sscanf(str, "Sigma=%lf\n", &temp.Sigma);
        break;
      }
      case 13:
      {
        sscanf(str, "deltaOut=%lf\n", &temp.deltaOut);
        break;
      }
      default:
        if (PROGRAMM_MODE == 5)
          printf("Program mode: start string\n");
        else
        {
          lexeme = strtok(str, " ");
          if (temp.array.values == NULL)
            temp.array.values = (double *)calloc(temp.array.Lx * temp.array.Ly * temp.array.Lz, sizeof(double));
          while (pos < temp.array.Lx * temp.array.Ly * temp.array.Lz)
          {
            buffer = strtof(lexeme, NULL);
            temp.array.values[pos] = buffer;
            lexeme = strtok(NULL, " ");
            pos++;
          }
          pos = 0;
        }
      }
    }

    printf("%lf\n", temp.t);
    printf("%lf\n", temp.T);
    printf("%lf\n", temp.deltaT);
    printf("%lf\n", temp.XMin);
    printf("%lf\n", temp.XMax);
    printf("%lf\n", temp.YMin);
    printf("%lf\n", temp.YMax);
    printf("%lf\n", temp.ZMin);
    printf("%lf\n", temp.ZMax);
    printf("%d\n", temp.array.Lx);
    printf("%d\n", temp.array.Ly);
    printf("%d\n", temp.array.Lz);
    printf("%lf\n", temp.Sigma);
    printf("%lf\n", temp.deltaOut);
    printf("\n");
    fclose(file);
    return temp;
  }
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

static inline void Output(char *DATA, double *arr, int length)
{
  FILE *file = fopen(DATA, "a");
  if (file == NULL)
    printf("Can`t open file\n");
  else
  {
    fprintf(file, "\n");
    for (int i = 0; i < length; i++)
      fprintf(file, "%lf ", arr[i]);
    fclose(file);
  }
}