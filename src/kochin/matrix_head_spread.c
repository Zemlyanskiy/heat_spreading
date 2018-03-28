#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define EILER 0;
#define RUNGE_KUTTA 1;
#define START_STRING 2;

struct input {
	double t;//START TIME (t)
	double T;//END TIME (T)
	double deltaT;//GAP BETWEEN MEASUREMENTS (dt)
	double Xmin; //START X (XMIN)
	double Xmax; //END X (XMAX)
	int NumPoints; // NUMBER OF MEASUREMENTS (LX)
	double Sigma; //THE COEFFICIENT OF QUATION  (SIGMA)
	double deltaOut;//DELTA OF OUTPUT(deltaOUT)
	double* values; //LAST ARRAY WITH VALUES
};

struct input ReadInput(char* PATH, int mode) {
	struct input temp;
	FILE *file;

	char* token;
	int pos;
	double buf;

	fopen_s(&file, PATH, "r");
	char string[5000];
	if (file == NULL) {
		printf("CANT FIND FILE\n");
	}
	else {
		printf("READ INPUT FILE:\n");
		for (int i = 0; fgets(string, sizeof(string), file) != NULL; i++)
			switch (i) {
			case 0:
			{
				sscanf_s(string, "t=%lf", &temp.t);
				printf("t = %f\n", temp.t);
				break;
			}
			case 1:
			{
				sscanf_s(string, "T=%lf", &temp.T);
				printf("T = %f\n", temp.T);
				break;
			}
			case 2:
			{
				sscanf_s(string, "deltaT=%lf", &temp.deltaT);
				printf("deltaT = %f\n", temp.deltaT);
				break;
			}
			case 3:
			{
				sscanf_s(string, "XMin=%lf", &temp.Xmin);
				printf("XMin = %f\n", temp.Xmin);
				break;
			}
			case 4:
			{
				sscanf_s(string, "XMax=%lf", &temp.Xmax);
				printf("XMax = %f\n", temp.Xmax);
				break;
			}
			case 5:
			{
				sscanf_s(string, "Lx=%d", &temp.NumPoints);
				printf("NumPoints = %d\n", temp.NumPoints);
				temp.values = malloc(sizeof(double)*temp.NumPoints);
				for (int i = 0; i < temp.NumPoints; i++) {
					temp.values[i] = 0;
				}
				break;
			}
			case 6:
			{
				sscanf_s(string, "Sigma=%lf", &temp.Sigma);
				printf("Sigma = %f\n", temp.Sigma);
				break;
			}
			case 7:
			{
				sscanf_s(string, "deltaOut=%lf", &temp.deltaOut);
				printf("deltaOut = %f\n", temp.deltaOut);
				break;
			}
			default:
				if (mode == 2)
					printf("WRONG LINE INPUT NUMBER\n");
				else {
					pos = 0;
					token = strtok(string, " ");
					while (token != NULL)
					{
						buf = strtof(token, NULL);
						temp.values[pos] = buf;
						token = strtok(NULL, " ");
						pos++;
					}
				}
			}
	}
	for (int i = 0; i < temp.NumPoints; i++)
		printf("%lf ", temp.values[i]);
	printf("\n");
	fclose(file);
	return temp;
}

double CalcFunc(double point) {
	double result;
	if ((point >= -0.5) && (point <= 0.5))
		result = cos(point*3.141592);
	else
		result = 0;
	return result;
}

int OutputArr(char* PATH, double* arr, int length) {
	FILE* file;
	fopen_s(&file, PATH, "a+");
	if (file == NULL)
		printf("Cant open file\n");
	else {
		fprintf(file, "\n");
		for (int i = 0; i < length; i++)
			fprintf(file, "%lf ", arr[i]);
		fclose(file);
	}

}

char* FloatToString(double in) {
	int IntPart = (int)in;
	double DoublePart = in - IntPart;
	int size = 1;//char #','
	int buf;
	if (IntPart == 0)
		size++;
	if (DoublePart == 0)
		size++;
	while (IntPart > 0) {
		size++;
		IntPart /= 10;
	}
	while (DoublePart > 0) {
		DoublePart *= 10;
		buf = (int)DoublePart;
		DoublePart -= buf;
		size++;
	}

}

int OutputCurrentTime(char* PATH, double currT) {
	FILE* file;
	int ERR_CODE, num_of_space = 15 - 8;
	int tmp = (int)currT;
	fopen_s(&file, PATH, "r+");
	if (file == NULL)
		printf("Cant open file\n");
	else {
		while (tmp > 10) {
			tmp /= 10;
			num_of_space--;
		}

		fseek(file, 0, SEEK_SET);
		fputs("t=", file);
		for (int i = 0; i<num_of_space; i++)
			fputs(" ", file);
		fprintf(file, "%f", currT);
		ERR_CODE = ftell(file);
		fclose(file);
	}
	return 0;
}

void CsrMult(double* Values, int* ColumnNum, int* LineFirst, double * Array, const int size, double* result)
{
//	#pragma omp parallel threads_num(4)
//	{
//		printf("hello");
//	}
	#pragma omp parallel for
	for (int i = 0; i < size; i++) {
		result[i] = 0;
		for (int j = LineFirst[i]; j < LineFirst[i + 1]; j++)
			result[i] += Values[j] * Array[ColumnNum[j]];
	}
}

int main(int argc, char *argv[]) {
	char* inputFile;
	int PROGRAMM_MODE = 1;
	int NUMBER_OF_THREADS = 1;
	if (argc == 2) {
		inputFile = argv[1];
	}
	else if (argc == 3) {
		inputFile = argv[1];
		printf("%c", argv[2]);
		if (argv[2][0] == 'e')
			PROGRAMM_MODE = 1;
		if (argv[2][0] == 'r')
			PROGRAMM_MODE = 0;
		if (argv[2][0] == 's')
			PROGRAMM_MODE = 2;
	}
	else if (argc == 4) {
		inputFile = argv[1];
		printf("%c", argv[2]);
		if (argv[2][0] == 'e')
			PROGRAMM_MODE = 1;
		if (argv[2][0] == 'r')
			PROGRAMM_MODE = 0;
		if (argv[2][0] == 's')
			PROGRAMM_MODE = 2;
		NUMBER_OF_THREADS = strtol(argv[3], NULL, 0);
	}
	else {
		printf("Incorrect number of arguments");
		return 4;
	}
	//READ INPUT DATA AND CREATE START SITUATION
	struct input data;//struct for input file
	double* points;//array for caclulate results
	double* prev;//array to previos state
	double posX; //var for calculate first state
	double step;

	double *k1, *k2, *k3, *k4, *medium;
	double *derivative;


	data = ReadInput(inputFile, PROGRAMM_MODE);
	points = (double*)malloc(sizeof(double)*data.NumPoints);
	posX = data.Xmin;
	step = (data.Xmax - data.Xmin) / data.NumPoints;
	//CSR Compressed Sparse Row
	double * Values = (double*)malloc(sizeof(double)*(data.NumPoints - 2)*3);
	int* ColumnNum = (int*)malloc(sizeof(int)*(data.NumPoints - 2) * 3);
	int* LineFirst = (int*)malloc(sizeof(int)*(data.NumPoints + 1));
	
	
	double OutTime = data.t;
	int OutCount = (int)(data.T - data.t) / data.deltaOut;
	int maxCount = (int)(data.T - data.t) / data.deltaT;
	OutCount = maxCount / OutCount;
	int count;

	if (PROGRAMM_MODE == 2) {
		for (int i = 0; i < data.NumPoints; i++) {
			points[i] = CalcFunc(posX);
			posX += step;
		}
		//OUTPUT FIRST RESULTS
		OutputArr(inputFile, points, data.NumPoints);
	}
	else
	{	//START PROCESS(EILER OR KUTTA)
		omp_set_num_threads(NUMBER_OF_THREADS);
		
		int i;
		double divider = 1 / (step*step);
		for (int i = 0; i < data.NumPoints; i++) {
			points[i] = data.values[i];
		}
		//INITIALIZE CSR
		for (int i = 0; i < data.NumPoints-2; i++) {
				Values[i * 3] = data.Sigma * data.deltaT * divider;
				Values[i * 3 + 1] = 1 - 2 * data.Sigma * data.deltaT * divider;
				Values[i * 3 + 2] = data.Sigma * data.deltaT * divider;
				ColumnNum[i * 3] = i;
				ColumnNum[i * 3 + 1] = i + 1;
				ColumnNum[i * 3 + 2] = i + 2;
				LineFirst[i+1] = i * 3;
		}
		LineFirst[0] = LineFirst[1];
		LineFirst[data.NumPoints - 1] = LineFirst[data.NumPoints - 2];
		LineFirst[data.NumPoints] = (data.NumPoints-2) * 3;
		//END OF INITIALIZATION
		prev = (double*)malloc(sizeof(double)*data.NumPoints);

		double time;
		double out_time;
		if (PROGRAMM_MODE == 1)//EILER
		{
			time = omp_get_wtime();
			for (count = 0; count <= maxCount; count++) {
				#pragma omp parallel for
				for (i = 0; i < data.NumPoints; i++)
					prev[i] = points[i];//make pointers swap
				/*
				#pragma omp parallel for
				for (i = 1; i < data.NumPoints - 1; i++)
					points[i] = prev[i] + data.Sigma*data.deltaT*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) * divider;
				//delete output time from time calculating*/
				CsrMult(Values, ColumnNum, LineFirst, prev, data.NumPoints, points);
				out_time = omp_get_wtime();
				if (count%OutCount == 0) {
					OutputArr(inputFile, points, data.NumPoints);//if ERROR will be after this line, numbers of lines and t parameter will be wrong
					OutputCurrentTime(inputFile, OutTime);
					OutTime += data.deltaOut;
				}
				out_time = omp_get_wtime() - out_time;
				time = time + out_time;
			}
		}
		if (PROGRAMM_MODE == 0) {//RUNGE_KUTTA
			k1 = (double*)malloc(sizeof(double)*data.NumPoints);
			k2 = (double*)malloc(sizeof(double)*data.NumPoints);
			k3 = (double*)malloc(sizeof(double)*data.NumPoints);
			k4 = (double*)malloc(sizeof(double)*data.NumPoints);
			medium = (double*)malloc(sizeof(double)*data.NumPoints);
			derivative = (double*)malloc(sizeof(double)*data.NumPoints);
			for (int i = 0; i < data.NumPoints; i++) {
				k1[i] = 0;
				k2[i] = 0;
				k3[i] = 0;
				k4[i] = 0;
				medium[i] = 0;
				derivative[i] = 0;
			}
			time = omp_get_wtime();
			for (count = 0; count <= maxCount; count++) {
				#pragma omp parallel for
				for (i = 0; i < data.NumPoints; i++)
					prev[i] = points[i];//make pointers swap
				//DERIVATIVE CALCULATING
				//1
				CsrMult(Values, ColumnNum, LineFirst, prev, data.NumPoints, k1);
//				#pragma omp parallel for
//				for (i = 1; i < data.NumPoints - 1; i++)
//					//1 - write, 5 - reads
//					k1[i] = prev[i]+data.deltaT*data.Sigma*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) * divider;
				
				#pragma omp parallel for
				for (i = 0; i < data.NumPoints; i++)
					medium[i] = (prev[i]+k1[i]) * 0.5;
				//2
				CsrMult(Values, ColumnNum, LineFirst, medium, data.NumPoints, k2);
//				#pragma omp parallel for
//				for (i = 1; i < data.NumPoints - 1; i++)
//					k2[i] = prev[i] + data.deltaT*data.Sigma*(medium[i - 1] - 2 * medium[i] + medium[i + 1]) * divider;
				#pragma omp parallel for
				for (i = 0; i < data.NumPoints; i++)
					medium[i] = (prev[i] + k2[i]) * 0.5;
				//3
				CsrMult(Values, ColumnNum, LineFirst, medium, data.NumPoints, k3);
//				#pragma omp parallel for
//				for (i = 1; i < data.NumPoints - 1; i++)
//					k3[i] = prev[i] + data.deltaT*data.Sigma*(medium[i - 1] - 2 * medium[i] + medium[i + 1]) * divider;
				#pragma omp parallel for
				for (i = 0; i < data.NumPoints; i++)
					medium[i] = k2[i];
				//4
				CsrMult(Values, ColumnNum, LineFirst, medium, data.NumPoints, k4);
//				#pragma omp parallel for
//				for (i = 1; i < data.NumPoints - 1; i++)
//					k4[i] = prev[i]+data.deltaT*data.Sigma*(medium[i - 1] - 2 * medium[i] + medium[i + 1]) * divider;
				#pragma omp parallel for
				for (i = 1; i < data.NumPoints - 1; i++) {
					derivative[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * 0.1666666666;
				}
				//NEXT POINT CALCULATING
				#pragma omp parallel for
				for (i = 1; i < data.NumPoints - 1; i++)
					points[i] = derivative[i];
				out_time = omp_get_wtime();
				if (count%OutCount == 0) {
					OutputArr(inputFile, points, data.NumPoints);//if ERROR will be after this line, numbers of lines and t parameter will be wrong
					OutputCurrentTime(inputFile, OutTime);
					OutTime += data.deltaOut;
				}
				out_time = omp_get_wtime() - out_time;
				time = time + out_time;
			}
		}
		time = omp_get_wtime() - time;
		printf("\nTIME IS %f%\n", time);
		scanf_s("ENTER ANY KEY TO EXIT %f", &time);
	}
	printf("\nEND\n");

	return 0;
}