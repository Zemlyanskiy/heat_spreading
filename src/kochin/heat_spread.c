#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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

struct input ReadInput(char* PATH,int mode) {
	struct input temp;
	FILE *file;

	char* token;
	int pos;
	double buf;

	fopen_s(&file,PATH, "r");
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
						buf = strtof(token,NULL);
						temp.values[pos] = buf;
						token = strtok(NULL, " ");
						pos++;
					}
				}
			}
	}
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
			fprintf(file, "%f ", arr[i]);
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
	int ERR_CODE,num_of_space=15-8;
	int tmp=(int)currT;
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
		for (int i=0;i<num_of_space;i++)
			fputs(" ", file);
		fprintf(file,"%f",currT);
		ERR_CODE = ftell(file);
		fclose(file);
	}
	return 0;
}

int main(int argc, char *argv[]) {
	char* inputFile;
	int PROGRAMM_MODE = EILER;
	if (argc == 2) {
		inputFile = argv[1];
	}
	else if (argc == 3) {
		inputFile = argv[1];
		if (argv[2] == 'e')
			PROGRAMM_MODE = EILER;
		if (argv[2] == 'r')
			PROGRAMM_MODE = RUNGE_KUTTA;
		if (argv[2] == 's')
			PROGRAMM_MODE = START_STRING;
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
	double OutTime = 0;
	int numOutputs;
	int outCoun;
	int count;
	data = ReadInput(inputFile, PROGRAMM_MODE);
	points = (double*)malloc(sizeof(double)*data.NumPoints);
	posX = data.Xmin;
	step = (data.Xmax - data.Xmin) / data.NumPoints;
	
	
	if (PROGRAMM_MODE == 2) {
		for (int i = 0; i < data.NumPoints; i++) {
			points[i] = CalcFunc(posX);
			posX += step;
		}
		//OUTPUT FIRST RESULTS
		OutputArr(inputFile, points, data.NumPoints);
	}
	else
	{	//START PROCESS(EILER)
		for (int i = 0; i < data.NumPoints; i++) {
			points[i] = data.values[i];
		}
		prev = (double*)malloc(sizeof(double)*data.NumPoints);
		numOutputs = (int)(data.T - data.t) / data.deltaOut;
		outCoun = (int)(data.T - data.t) / data.deltaT;
		outCoun = outCoun / numOutputs;
		count = 0;
		for (double time = data.t; time <= data.T; time += data.deltaT, count++) {
			for (int i = 0; i < data.NumPoints; i++)
				prev[i] = points[i];//make pointers swap
			for (int i = 1; i < data.NumPoints - 1; i++)
				points[i] = prev[i] + data.Sigma*data.deltaT*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) / pow(step, 2);
			if (count%outCoun == 0) {
				OutputArr(inputFile, points, data.NumPoints);//if ERROR will be after this line, numbers of lines and t parameter will be wrong
				OutputCurrentTime(inputFile, OutTime);
				OutTime += data.deltaOut;
			}
		}
	}
	printf("\n\nEND\n\n");
	return 0;
}