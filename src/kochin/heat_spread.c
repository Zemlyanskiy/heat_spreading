#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct input {
	double t;//START TIME (t)
	double T;//END TIME (T)
	double deltaT;//GAP BETWEEN MEASUREMENTS (dt)
	double Xmin; //START X (XMIN)
	double Xmax; //END X (XMAX)
	int NumPoints; // NUMBER OF MEASUREMENTS (LX)
	double Sigma; //THE COEFFICIENT OF QUATION  (SIGMA)
	double deltaOut;//DELTA OF OUTPUT(deltaOUT)
};

struct input ReadInput(char* PATH) {
	struct input temp;
	FILE *file;
	fopen_s(&file,PATH, "r");
	char string[100];
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
				printf("WRONG LINE INPUT NUMBER\n");
			}
	}
	fclose(file);
	return temp;
}

int Output(struct input temp, double** points, unsigned int len) {
	printf("t = %f\n", temp.t);
	printf("T = %f\n", temp.T);
	printf("deltaT = %f\n", temp.deltaT);
	printf("Xmin = %f\n", temp.Xmin);
	printf("Xmax = %f\n", temp.Xmax);
	printf("NumPoints = %d\n", temp.NumPoints);
	printf("Sigma = %f\n", temp.Sigma);
	printf("deltaOut=%f\n", temp.deltaOut);
	if (points != NULL) {
		for (int i = 0; i < len; i++)
			printf("%f ", points[0][i]);
		printf("\n");
		for (int i = 0; i < len; i++)
			printf("%f ", points[1][i]);
		printf("\n");
	}
}

double CalcFunc(double point) {
	double result;
	if ((point >= -0.5) && (point <= 0.5))
		result = cos(point*3.141592);
	else
		result = 0;
	return result;
}

int OutputArr(char* PATH, char* mode, double* arr, int length) {
	FILE* file;
	fopen_s(&file, PATH, mode);
	if (file == NULL)
		printf("Cant open file\n");
	else {
		for (int i = 0; i < length; i++)
			fprintf(file, "%f ", arr[i]);
		fprintf(file, "\n");
		fclose(file);
	}

}

int main() {

	//READ INPUT DATA AND CREATE START SITUATION
	struct input data;//struct for input file
	double* points;//array for caclulate results
	double* prev;//array to previos state
	double posX; //var for calculate first state
	double step;
	int numOutputs;
	int outCoun;
	int count;
	data = ReadInput("input1.txt");
	step = (data.Xmax - data.Xmin) / data.NumPoints;
	points = (double*)malloc(sizeof(double)*data.NumPoints);
	prev = (double*)malloc(sizeof(double)*data.NumPoints);
	posX = data.Xmin;
	for (int i = 0; i < data.NumPoints; i++) {
		points[i] = CalcFunc(posX);
		posX += step;
	}

	//OUTPUT FIRST RESULTS
	OutputArr("input1.txt", "a", points, data.NumPoints);
	//START PROCESS
	numOutputs =(int) (data.T - data.t) / data.deltaOut;
	outCoun = (int)(data.T - data.t) / data.deltaT;
	outCoun = outCoun / numOutputs;
	count = 0;
	for (double time = data.t; time <= data.T; time += data.deltaT, count++) {
		for (int i = 0; i < data.NumPoints; i++)
			prev[i] = points[i];
		for (int i = 1; i < data.NumPoints - 1; i++)
			points[i] = prev[i] + data.Sigma*data.deltaT*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) / pow(step, 2);
		if (count%outCoun == 0) {
			OutputArr("input1.txt", "a", points, data.NumPoints);
		}
	}
	printf("\n\nEND\n\n");
	return 0;
}