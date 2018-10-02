#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include<mpi.h>

enum Method {
	EILER,
	RUNGE_KUTTA,
	START_STRING,
	MATRIX_EILER,
	MATRIX_RUNGE_KUTTA,
	YAKOBI
};

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
	FILE *file;
	struct input temp = { .t = 0, .T=0, .deltaT=0, .Xmin=0, .Xmax=0, .NumPoints=0, .Sigma=0, .deltaOut=0, .values=0 };

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
	int i;
#pragma omp parallel for
	for (i = 0; i < size; i++) {
		result[i] = 0;
		for (int j = LineFirst[i]; j < LineFirst[i + 1]; j++)
			result[i] += Values[j] * Array[ColumnNum[j]];
	}
}

double CrsAccess(double* Values, int* ColumnNum, int* LineFirst, int i, int j, int length) {
	int result;
	if (i >= length)
		return 0;
	int n1 = LineFirst[i];
	int n2 = LineFirst[i + 1];
	int k;
	for (k = n1; k<n2; k++)
		if (ColumnNum[k] == j) {
			return Values[k];
		}
	return 0;
}

int main(int argc, char *argv[]) {
	//PROCESSING CL ARGUMENTS
	char* inputFile = NULL;
	enum Method PROGRAMM_MODE = EILER;
	int NUMDER_OF_THREADS = 1;
	{

		if (argc > 1)
			inputFile = argv[1];
		if (argc > 2) {
			if (argv[2][0] == 'e' && argv[2][1] == '\0')
				PROGRAMM_MODE = EILER;
			else if (argv[2][0] == 'r' && argv[2][1] == '\0')
				PROGRAMM_MODE = RUNGE_KUTTA;
			else if (argv[2][0] == 's' && argv[2][1] == '\0')
				PROGRAMM_MODE = START_STRING;
			else if (argv[2][0] == 'm'&& argv[2][1] == 'e'&& argv[2][2] == '\0')
				PROGRAMM_MODE = MATRIX_EILER;
			else if (argv[2][0] == 'm'&& argv[2][1] == 'r'&& argv[2][2] == '\0')
				PROGRAMM_MODE = MATRIX_RUNGE_KUTTA;
			else if (argv[2][0] == 'y' && argv[2][1] == '\0')
				PROGRAMM_MODE = YAKOBI;
			else {
				printf("WRONG FORMAT\n");
				return 127;
			}
		}
		if (argc > 3)
			NUMDER_OF_THREADS = strtol(argv[3], NULL, 0);
		if (argc > 4) {
			printf("Incorrect number of arguments");
			return 4;
		}
	}

	//READ INPUT DATA AND INTIALIZE COMMON VARIABLES
	struct input data;//struct for input file
	double* points;//array for caclulate results
	double* prev;//array to previos state
	double posX; //var for calculate first state
	double step;

	int count = 0, i = 0;

	data = ReadInput(inputFile, PROGRAMM_MODE);
	points = (double*)malloc(sizeof(double)*data.NumPoints);
	prev = (double*)malloc(sizeof(double)*data.NumPoints);
	posX = data.Xmin;
	step = (data.Xmax - data.Xmin) / data.NumPoints;

	//IF PROGRAM MODE == START_STRING - GENERATE IT AND COMPLETE PROGRAMM
	if (PROGRAMM_MODE == START_STRING) {
		for (i = 0; i < data.NumPoints; i++) {
			points[i] = CalcFunc(posX);
			posX += step;
		}
		//OUTPUT FIRST RESULTS
		OutputArr(inputFile, points, data.NumPoints);
		return 0;
	}

	double *k1 = NULL, *k2 = NULL, *k3 = NULL, *k4 = NULL,
		*medium = NULL, *derivative = NULL,
		*discrepancy = NULL, *MainDiag = NULL;
	double *Values = NULL;
	int *ColumnNum = NULL, *LineFirst = NULL;

	double OutTime = data.t;
	int OutCount = (int)(data.T - data.t) / data.deltaOut;
	int maxCount = (int)(data.T - data.t) / data.deltaT;
	OutCount = maxCount / OutCount;

	double divider = 1 / (step*step);
	double delta, sum;
	int flag;

	//INITIALIZE CUSTOM VARIABLES
	//CSR Compressed Sparse Row
	if (PROGRAMM_MODE == MATRIX_EILER || PROGRAMM_MODE == MATRIX_RUNGE_KUTTA) {
		Values = (double*)malloc(sizeof(double)*(data.NumPoints - 2) * 3);
		ColumnNum = (int*)malloc(sizeof(int)*(data.NumPoints - 2) * 3);
		LineFirst = (int*)malloc(sizeof(int)*(data.NumPoints + 1));

		//INITIALIZE CSR
		for (i = 0; i < data.NumPoints - 2; i++) {
			Values[i * 3] = data.Sigma * data.deltaT * divider;
			Values[i * 3 + 1] = 1 - 2 * data.Sigma * data.deltaT * divider;
			Values[i * 3 + 2] = data.Sigma * data.deltaT * divider;
			ColumnNum[i * 3] = i;
			ColumnNum[i * 3 + 1] = i + 1;
			ColumnNum[i * 3 + 2] = i + 2;
			LineFirst[i + 1] = i * 3;
		}
		LineFirst[0] = LineFirst[1];
		LineFirst[data.NumPoints - 1] = LineFirst[data.NumPoints - 2];
		LineFirst[data.NumPoints] = (data.NumPoints - 2) * 3;
		//END OF INITIALIZATION
	}
	else if (PROGRAMM_MODE == YAKOBI) {
		Values = (double*)malloc(sizeof(double)*(data.NumPoints - 2) * 3 + 2);
		ColumnNum = (int*)malloc(sizeof(int)*(data.NumPoints - 2) * 3 + 2);
		LineFirst = (int*)malloc(sizeof(int)*(data.NumPoints + 1));


		//INITIALIZE CSR
		for (i = 0; i < data.NumPoints - 2; i++) {
			Values[i * 3 + 1] = -data.Sigma * data.deltaT * divider;
			Values[i * 3 + 2] = 1 + 2 * data.Sigma * data.deltaT * divider;
			Values[i * 3 + 3] = -data.Sigma * data.deltaT * divider;
			ColumnNum[i * 3 + 1] = i;
			ColumnNum[i * 3 + 2] = i + 1;
			ColumnNum[i * 3 + 3] = i + 2;
			LineFirst[i + 1] = i * 3 + 1;
		}

		Values[0] = 1;
		ColumnNum[0] = 0;
		LineFirst[0] = 0;

		Values[(data.NumPoints - 2) * 3 + 1] = 1;
		ColumnNum[(data.NumPoints - 2) * 3 + 1] = data.NumPoints - 1;
		LineFirst[data.NumPoints - 1] = LineFirst[data.NumPoints - 2] + 3;

		LineFirst[data.NumPoints] = (data.NumPoints - 2) * 3 + 2;

		delta = 0.00000000001;
		discrepancy = (double*)malloc(sizeof(double)*data.NumPoints);
		flag = 0;
		sum;

		MainDiag = (double*)malloc(sizeof(double)*data.NumPoints);
		for (i = 0; i < data.NumPoints; i++) {
			MainDiag[i] = 1 / CrsAccess(Values, ColumnNum, LineFirst, i, i, data.NumPoints);
		}

		//END OF INITIALIZATION

	}

	if (PROGRAMM_MODE == MATRIX_RUNGE_KUTTA || PROGRAMM_MODE == RUNGE_KUTTA) {
		k1 = (double*)malloc(sizeof(double)*data.NumPoints);
		k2 = (double*)malloc(sizeof(double)*data.NumPoints);
		k3 = (double*)malloc(sizeof(double)*data.NumPoints);
		k4 = (double*)malloc(sizeof(double)*data.NumPoints);
		medium = (double*)malloc(sizeof(double)*data.NumPoints);
		derivative = (double*)malloc(sizeof(double)*data.NumPoints);

		for (i = 0; i < data.NumPoints; i++) {
			k1[i] = 0;
			k2[i] = 0;
			k3[i] = 0;
			k4[i] = 0;
			medium[i] = 0;
			derivative[i] = 0;
		}
	}

	for (i = 0; i < data.NumPoints; i++)
		points[i] = data.values[i];

	//ENV PREPARATION FOR OpenMP
	omp_set_num_threads(NUMDER_OF_THREADS);
	double time = 0;
	double out_time;
	//TODO: ENV PREPARATION FOR MPI
	int rank, proc_num;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	unsigned step_length = data.NumPoints / proc_num;
	int* start_positions = (int*)malloc(sizeof(int)*proc_num);
	int* steps_length = (int*)malloc(sizeof(int)*proc_num);
	for (i = 0; i < proc_num; i++) {
		start_positions[i] = i * step_length;
		steps_length[i] = steps_length + 1;
	}
	MPI_Scatterv(data.values, steps_length, start_positions,
		MPI_INT, points, step_length + 1,
		MPI_INT, 0, MPI_COMM_WORLD);

	/*int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
					MPI_Datatype sendtype, void *recvbuf, int recvcount,
					MPI_Datatype recvtype,
					int root, MPI_Comm comm)
	Input Parameters
			sendbuf - address of send buffer(choice, significant only at root)
			sendcounts - integer array(of length group size) specifying the number of elements to send to each processor
			displs - integer array(of length group size).Entry i specifies the displacement(relative to sendbuf from which to take the outgoing data to process i
			sendtype - data type of send buffer elements(handle)
			recvcount - number of elements in receive buffer(integer)
			recvtype - data type of receive buffer elements(handle)
			root - rank of sending process(integer)
			comm - communicator(handle)
	*/

	//START CALCULATING
	if (PROGRAMM_MODE == EILER)
	{
		printf("I am process number %d and my array is:\n", rank);
		for (int i = 0; i < step_length + 1; i++)
			printf("%d ", points[i]);
		printf("\n");
		rank = getchar();
	}
	else
	{
		for (count = 0; count <= maxCount; count++) {
			out_time = omp_get_wtime();
			#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				prev[i] = points[i];//make pointers swap
			#pragma omp parallel for
			for (i = 1; i < data.NumPoints - 1; i++)
				points[i] = prev[i] + data.Sigma*data.deltaT*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) * divider;
			//delete output time from time calculating
			time += omp_get_wtime() - out_time;
			if (count%OutCount == 0) {
				OutputArr(inputFile, points, data.NumPoints);//if ERROR will be after this line, numbers of lines and t parameter will be wrong
				OutputCurrentTime(inputFile, OutTime);
				OutTime += data.deltaOut;
			}
		}
	}
	if (PROGRAMM_MODE == RUNGE_KUTTA) {
		for (count = 0; count <= maxCount; count++) {
			out_time = omp_get_wtime();
#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				prev[i] = points[i];//make pointers swap
									//DERIVATIVE CALCULATING
									//1
#pragma omp parallel for
			for (i = 1; i < data.NumPoints - 1; i++)
				//1 - write, 4 - reads, store 1/pow(step, 2) in separate variable
				k1[i] = data.deltaT*data.Sigma*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) * divider;
#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				medium[i] = prev[i] + k1[i] * 0.5;
			//2
#pragma omp parallel for
			for (i = 1; i < data.NumPoints - 1; i++)
				k2[i] = data.deltaT*data.Sigma*(medium[i - 1] - 2 * medium[i] + medium[i + 1]) * divider;
#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				medium[i] = prev[i] + k2[i] * 0.5;
			//3
#pragma omp parallel for
			for (i = 1; i < data.NumPoints - 1; i++)
				k3[i] = data.deltaT*data.Sigma*(medium[i - 1] - 2 * medium[i] + medium[i + 1]) * divider;
#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				medium[i] = prev[i] + k2[i];
			//4
#pragma omp parallel for
			for (i = 1; i < data.NumPoints - 1; i++)
				k4[i] = data.deltaT*data.Sigma*(medium[i - 1] - 2 * medium[i] + medium[i + 1]) * divider;
#pragma omp parallel for
			for (i = 1; i < data.NumPoints - 1; i++) {
				derivative[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * 0.1666666666;
			}
			//NEXT POINT CALCULATING
#pragma omp parallel for
			for (i = 1; i < data.NumPoints - 1; i++)
				points[i] = prev[i] + derivative[i];
			time += omp_get_wtime() - out_time;
			if (count%OutCount == 0) {
				OutputArr(inputFile, points, data.NumPoints);//if ERROR will be after this line, numbers of lines and t parameter will be wrong
				OutputCurrentTime(inputFile, OutTime);
				OutTime += data.deltaOut;
			}
		}
	}
	if (PROGRAMM_MODE == MATRIX_EILER) {
		for (count = 0; count <= maxCount; count++) {
			out_time = omp_get_wtime();
#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				prev[i] = points[i];//make pointers swap

			CsrMult(Values, ColumnNum, LineFirst, prev, data.NumPoints, points);
			time += omp_get_wtime() - out_time;
			if (count%OutCount == 0) {
				OutputArr(inputFile, points, data.NumPoints);//if ERROR will be after this line, numbers of lines and t parameter will be wrong
				OutputCurrentTime(inputFile, OutTime);
				OutTime += data.deltaOut;
			}
		}
	}
	if (PROGRAMM_MODE == MATRIX_RUNGE_KUTTA) {
		for (count = 0; count <= maxCount; count++) {
			out_time = omp_get_wtime();
#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				prev[i] = points[i];//make pointers swap
									//DERIVATIVE CALCULATING
									//1
			CsrMult(Values, ColumnNum, LineFirst, prev, data.NumPoints, k1);

#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				medium[i] = (prev[i] + k1[i]) * 0.5;
			//2
			CsrMult(Values, ColumnNum, LineFirst, medium, data.NumPoints, k2);
#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				medium[i] = (prev[i] + k2[i]) * 0.5;
			//3
			CsrMult(Values, ColumnNum, LineFirst, medium, data.NumPoints, k3);
#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				medium[i] = k2[i];
			//4
			CsrMult(Values, ColumnNum, LineFirst, medium, data.NumPoints, k4);
#pragma omp parallel for
			for (i = 1; i < data.NumPoints - 1; i++) {
				derivative[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * 0.1666666666;
			}
			//NEXT POINT CALCULATING
#pragma omp parallel for
			for (i = 1; i < data.NumPoints - 1; i++)
				points[i] = derivative[i];
			time += omp_get_wtime() - out_time;
			if (count%OutCount == 0) {
				OutputArr(inputFile, points, data.NumPoints);//if ERROR will be after this line, numbers of lines and t parameter will be wrong
				OutputCurrentTime(inputFile, OutTime);
				OutTime += data.deltaOut;
			}
		}
	}
	if (PROGRAMM_MODE == YAKOBI) {
		for (count = 0; count <= maxCount; count++) {
			out_time = omp_get_wtime();
			flag = 1;
#pragma omp parallel for
			for (i = 0; i < data.NumPoints; i++)
				prev[i] = points[i];//make pointers swap
			while (flag) {
				flag = 0;
				//TRY TO CHECK DISCREPANCY DELTA
				CsrMult(Values, ColumnNum, LineFirst, points, data.NumPoints, discrepancy);
#pragma omp parallel for
				for (i = 0; i < data.NumPoints; i++) {
					discrepancy[i] -= prev[i];
					discrepancy[i] = fabs(discrepancy[i]);
					if (discrepancy[i] >= delta) {
						flag = 1;
					}
				}
				if (!flag)
					break;

#pragma omp parallel for
				for (i = 0; i < data.NumPoints; i++)
					discrepancy[i] = points[i];

				for (i = 0; i < data.NumPoints; i++) {
					sum = 0;
					for (int j = 0; j < data.NumPoints; j++)
						if (j != i)
							sum += CrsAccess(Values, ColumnNum, LineFirst, i, j, data.NumPoints)*discrepancy[j];
					points[i] = MainDiag[i] * (prev[i] - sum);
				}
			}
			time += omp_get_wtime() - out_time;
			if (count%OutCount == 0) {
				OutputArr(inputFile, points, data.NumPoints);//if ERROR will be after this line, numbers of lines and t parameter will be wrong
				OutputCurrentTime(inputFile, OutTime);
				OutTime += data.deltaOut;
			}
		}
	}


	printf("\nTIME IS %f%\n", time);
	scanf_s("ENTER ANY KEY TO EXIT %f", &time);
	printf("\nEND");

	return 0;
}