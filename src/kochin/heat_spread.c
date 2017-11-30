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
				sscanf_s(string, "%lf", &temp.t);
				printf("t = %f\n", temp.t);
				break;
			}
			case 1:
			{
				sscanf_s(string, "%lf", &temp.T);
				printf("T = %f\n", temp.T);
				break;
			}
			case 2:
			{
				sscanf_s(string, "%lf", &temp.deltaT);
				printf("deltaT = %f\n", temp.deltaT);
				break;
			}
			case 3:
			{
				sscanf_s(string, "%lf", &temp.Xmin);
				printf("Xmin = %f\n", temp.Xmin);
				break;
			}
			case 4:
			{
				sscanf_s(string, "%lf", &temp.Xmax);
				printf("Xmax = %f\n", temp.Xmax);
				break;
			}
			case 5:
			{
				temp.NumPoints = atoi(string);
				printf("NumPoints = %d\n", temp.NumPoints);
				break;
			}
			case 6:
			{
				sscanf_s(string, "%lf", &temp.Sigma);
				printf("Sigma = %f\n", temp.Sigma);
				break;
			}
			default:
				printf("WRONG LINE INPUT NUMBER\n");
			}
	}
	fclose(file);
	return temp;
}

int Output(struct input temp, double** points,unsigned int len) {
	printf("t = %f\n", temp.t);
	printf("T = %f\n", temp.T);
	printf("deltaT = %f\n", temp.deltaT);
	printf("Xmin = %f\n", temp.Xmin);
	printf("Xmax = %f\n", temp.Xmax);
	printf("NumPoints = %d\n", temp.NumPoints);
	printf("Sigma = %f\n", temp.Sigma);
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
		result = cos(point);
	else
		result = 0;
	return result;
}

int main() {

	//READ INPUT DATA AND CREATE START SITUATION
	struct input data;
	double* points[2];
	double* results;
	double* derivates;
	double* prev;
	double posX;
	double step;
	int time;
	int in; //buffer for user input

	data = ReadInput("input1.txt");
	step = (data.Xmax - data.Xmin) / data.NumPoints;
	results = (double*)malloc(sizeof(double)*data.NumPoints);
	derivates = (double*)malloc(sizeof(double)*data.NumPoints);
	prev = (double*)malloc(sizeof(double)*data.NumPoints);
	for (int i = 0; i < 2; i++)
		points[i] = (double*)malloc(sizeof(double)*data.NumPoints);
	posX = data.Xmin;

	for (int i = 0; i < data.NumPoints; i++) {
		points[0][i] = posX;
		points[1][i] = CalcFunc(posX);
		posX += step;
	}
	for (int i = 0; i < data.NumPoints; i++)
		printf("%.2f ", points[0][i]);
	printf("\n");
	for (int i = 0; i < data.NumPoints; i++)
		printf("%.2f ", points[1][i]);
	printf("\n");


/*	FILE* file;
	double tmp=1.0;
	fopen_s(&file, "input1.txt", "ab");
	if (file == NULL)
		printf("Cant open file\n");
	else {
		fprintf(file, "\n");
		for (int i = 0; i < data.NumPoints; i++)
			fprintf(file, "%f ", points[0][i]);
		fprintf(file, "\n");
		for (int i = 0; i < data.NumPoints; i++)
			fprintf(file, "%f ", points[1][i]);
		fclose(file);
	}
*/
	//START OF PROCESS
	for (int i = 0; i < data.NumPoints; i++)
		prev[i] = points[1][i];

	derivates[0] = 0;
	derivates[data.NumPoints] = 0;
	results[0] = 0;
	results[data.NumPoints] = 0;

	for (double posT = data.t; posT <= data.T; posT += data.deltaT) {
		for (int i = 1; i < data.NumPoints - 1; i++) {
			derivates[i] = data.Sigma*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) / pow(step, 2);
			results[i] = prev[i] + derivates[i];
		}
		for (int i = 0; i < data.NumPoints; i++)
			prev[i] = results[i];
		time = posT * 1000;
		if (time % 100 == 0) {
			for (int i = 0; i < data.NumPoints; i++)
				printf("%f ", prev[i]);
			printf("\n\n\n\n");
		}
		scanf_s("%d", &in);
	}
	return 0;
}