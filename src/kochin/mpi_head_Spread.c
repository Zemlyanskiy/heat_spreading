#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include<mpi.h>

#define DEBUG 0
#define START_STRING 0
#define EILER 0
#define RUNGE_KUTTA 0
#define MATRIX_EILER 0
#define MATRIX_RUNGE_KUTTA 0
#define YAKOBI 1

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

inline struct input ReadInput(char* PATH, int proc_rank) {
	FILE *file;
	struct input temp = { .t = 0, .T=0, .deltaT=0, .Xmin=0, .Xmax=0, .NumPoints=0, .Sigma=0, .deltaOut=0, .values=0 };

	char* token;
	int pos;
	double buf;

	fopen_s(&file, PATH, "r");
	char string[5000];
	if (file == NULL) {
		if (!proc_rank)
		printf("CANT FIND FILE\n");
	}
	else {
		if (!proc_rank)
		printf("READ INPUT FILE:\n");
		for (int i = 0; fgets(string, sizeof(string), file) != NULL; i++)
			switch (i) {
			case 0:
			{
				sscanf_s(string, "t=%lf", &temp.t);
				if (!proc_rank)
				printf("t = %f\n", temp.t);
				break;
			}
			case 1:
			{
				sscanf_s(string, "T=%lf", &temp.T);
				if (!proc_rank)
				printf("T = %f\n", temp.T);
				break;
			}
			case 2:
			{
				sscanf_s(string, "deltaT=%lf", &temp.deltaT);
				if (!proc_rank)
				printf("deltaT = %f\n", temp.deltaT);
				break;
			}
			case 3:
			{
				sscanf_s(string, "XMin=%lf", &temp.Xmin);
				if (!proc_rank)
				printf("XMin = %f\n", temp.Xmin);
				break;
			}
			case 4:
			{
				sscanf_s(string, "XMax=%lf", &temp.Xmax);
				if (!proc_rank)
				printf("XMax = %f\n", temp.Xmax);
				break;
			}
			case 5:
			{
				sscanf_s(string, "Lx=%d", &temp.NumPoints);
				if (!proc_rank)
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
				if (!proc_rank)
				printf("Sigma = %f\n", temp.Sigma);
				break;
			}
			case 7:
			{
				sscanf_s(string, "deltaOut=%lf", &temp.deltaOut);
				if (!proc_rank)
				printf("deltaOut = %f\n", temp.deltaOut);
				break;
			}
			default:
#if START_STRING
				if (!proc_rank)
				    printf("Mode: start string, Problem: file allready contain start string\n");
#endif
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
	if (!proc_rank) {
		for (int i = 0; i < temp.NumPoints; i++)
			printf("%lf ", temp.values[i]);
		printf("\n");
	}
	fclose(file);
	return temp;
}

inline int OutputArr(char* PATH, double* arr, int length) {
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

#if START_STRING

double CalcFunc(double point) {
	double result;
	if ((point >= -0.5) && (point <= 0.5))
		result = cos(point*3.141592);
	else
		result = 0;
	return result;
}

#else

inline char* FloatToString(double in) {
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

inline int OutputCurrentTime(char* PATH, double currT) {
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

inline void CsrMult(double* Values, int* ColumnNum, int* LineFirst, double * Array, const int size, double* result)
{
	int i;
#pragma omp parallel for
	for (i = 0; i < size; i++) {
		result[i] = 0;
		for (int j = LineFirst[i]; j < LineFirst[i + 1]; j++)
			result[i] += Values[j] * Array[ColumnNum[j]];
	}
}

inline void SendBorderValues(double * points, int* sendcounts, int rank, int proc_num) {
    if (proc_num != 1) {
        MPI_Status Status;
        //Send Recv via 0 and 1 process
        if (!rank) {
            // send pre last element to 1 process
            MPI_Send(points + sendcounts[rank] - 2, 1, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
            // get last element from 1 process
            MPI_Recv(points + sendcounts[rank] - 1, 1, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &Status);
        }
        //Send Recv via last and pre last process
        else if (rank == proc_num - 1)
        {
            // send second element to pre last process
            MPI_Send(points + 1, 1, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
            // get first element from pre last process
            MPI_Recv(points, 1, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &Status);
        }
        //Other Threads
        else
        {
            // send second and pre last elements to rank - 1 and rank + 1 process
            MPI_Send(points + 1, 1, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
            MPI_Send(points + sendcounts[rank] - 2, 1, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
            // get first and last elements from rank - 1 and rank + 1 process
            MPI_Recv(points, 1, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, &Status);
            MPI_Recv(points + sendcounts[rank] - 1, 1, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, &Status);
        }
    }
}

inline double CrsAccess(double* Values, int* ColumnNum, int* LineFirst, int i, int j, int length) {
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

#endif

int main(int argc, char *argv[]) {
    int rank, proc_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //PROCESSING CL ARGUMENTS
    char* inputFile = NULL;
    int NUMDER_OF_THREADS = 1;
    {

        if (argc > 1)
            inputFile = argv[1];
        if (argc > 2)
            NUMDER_OF_THREADS = strtol(argv[3], NULL, 0);
        if (argc > 3) {
            printf("Incorrect number of arguments");
            return 4;
        }
    }

    //READ INPUT DATA AND INTIALIZE COMMON VARIABLES
    struct input data;//struct for input file
    double* points;//array for caclulate results
    double* prev;//array to previos state
    double* print_array;
    double posX; //var for calculate first state
    double step;

    int count = 0, i = 0;

    data = ReadInput(inputFile, 1);
    unsigned elements_per_process = data.NumPoints / proc_num;
    points = (double*)malloc(sizeof(double)* (elements_per_process + proc_num));
    prev = (double*)malloc(sizeof(double)* (elements_per_process + proc_num));
    print_array = (double*)malloc(sizeof(double)* data.NumPoints);
    posX = data.Xmin;
    step = (data.Xmax - data.Xmin) / data.NumPoints;


#if START_STRING
    if (!rank) {
        for (i = 0; i < data.NumPoints; i++) {
            points[i] = CalcFunc(posX);
            posX += step;
        }
        //OUTPUT FIRST RESULTS
        OutputArr(inputFile, points, data.NumPoints);
    }
    MPI_Finalize();
    return 0;
}
#else
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

    //ENV PREPARATION FOR MPI
    int * sendcounts = malloc(sizeof(int)*proc_num);
    int * displs = malloc(sizeof(int)*proc_num);
    for (int i = 0; i < proc_num; i++) {
        sendcounts[i] = elements_per_process;
        displs[i] = elements_per_process * i;
    }
    //remainder(ostatok) of the division
    sendcounts[proc_num - 1] += data.NumPoints % proc_num;
    //send one symbol from each side of sent array
    if (proc_num != 1) {
        //add to end of array in first process 1 element
        sendcounts[0] += 1;
        //and to begin of array in last process 1 element
        displs[proc_num - 1] -= 1;
        sendcounts[proc_num - 1] += 1;
        //add for 1 element to begin and end other processes
        for (int i = 1; i < proc_num - 1; i++) {
            displs[i] -= 1;
            sendcounts[i] += 2;
        }
    }

	//INITIALIZE CUSTOM VARIABLES
	//CSR Compressed Sparse Row
#if  MATRIX_EILER ||  MATRIX_RUNGE_KUTTA
        Values = (double*)malloc(sizeof(double)*((sendcounts[rank] - 2) * 3 + 2));
        ColumnNum = (int*)malloc(sizeof(int)*((sendcounts[rank] - 2) * 3 + 2));
        LineFirst = (int*)malloc(sizeof(int)*(sendcounts[rank] + 1));
        for (i = 0; i < (sendcounts[rank] - 2) * 3+2; i++){
            Values[i] = 0;
            ColumnNum[i] = 0;
        }
        for (i = 0; i < sendcounts[rank] +1; i++) {
            LineFirst[i] = 0;
        }
        //INITIALIZE CSR
        for (i = 0; i < sendcounts[rank] - 2; i++) {
            Values[i * 3 + 1] = data.Sigma * data.deltaT * divider;
            Values[i * 3 + 2] = 1 - 2 * data.Sigma * data.deltaT * divider;
            Values[i * 3 + 3] = data.Sigma * data.deltaT * divider;
            ColumnNum[i * 3 + 1] = i;
            ColumnNum[i * 3 + 2] = i +1;
            ColumnNum[i * 3 + 3] = i +2;
            LineFirst[i + 1] = i * 3 + 1;
        }
        Values[0] = 1;
        LineFirst[0] = 0;
        ColumnNum[0] = 0;
        
        Values[(sendcounts[rank] - 2) * 3 + 1] = 1;
        ColumnNum[(sendcounts[rank] - 2) * 3 + 1] = sendcounts[rank]-1;
        LineFirst[sendcounts[rank] - 1] = LineFirst[sendcounts[rank] - 2] + 3;
        LineFirst[sendcounts[rank]] = (sendcounts[rank] - 2) * 3 + 2;
#endif /*MATRIX_EILER ||  MATRIX_RUNGE_KUTTA*/

#if YAKOBI
		Values = (double*)malloc(sizeof(double)*((sendcounts[rank] - 2) * 3 + 2));
		ColumnNum = (int*)malloc(sizeof(int)*((sendcounts[rank] - 2) * 3 + 2));
		LineFirst = (int*)malloc(sizeof(int)*(sendcounts[rank] + 1));


		//INITIALIZE CSR
		for (i = 0; i < sendcounts[rank] - 2; i++) {
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

		Values[(sendcounts[rank] - 2) * 3 + 1] = 1;
		ColumnNum[(sendcounts[rank] - 2) * 3 + 1] = sendcounts[rank] - 1;
		LineFirst[sendcounts[rank] - 1] = LineFirst[sendcounts[rank] - 2] + 3;

		LineFirst[sendcounts[rank]] = (sendcounts[rank] - 2) * 3 + 2;

		delta = 0.00000000001;
		discrepancy = (double*)malloc(sizeof(double)*sendcounts[rank]);
		flag = 0;
		sum;

		MainDiag = (double*)malloc(sizeof(double)*sendcounts[rank]);
		for (i = 0; i < data.NumPoints; i++) {
			MainDiag[i] = 1 / CrsAccess(Values, ColumnNum, LineFirst, i, i, sendcounts[rank]);
		}

		//END OF INITIALIZATION

#endif /*YAKOBI*/
#if RUNGE_KUTTA || MATRIX_RUNGE_KUTTA
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
#endif /*RUNGE_KUTTA || MATRIX_RUNGE_KUTTA*/
#if DEBUG
    // set values to 1,2,3,4...
	if (!rank)
		for (i = 0; i < data.NumPoints; i++)
			data.values[i] = i;
#endif
	//ENV PREPARATION FOR OpenMP
	omp_set_num_threads(NUMDER_OF_THREADS);
	double time = 0;
	double out_time;

    MPI_Status Status;
	MPI_Scatterv(data.values, sendcounts, displs,
		MPI_DOUBLE, points, elements_per_process+proc_num,
		MPI_DOUBLE,
		0, MPI_COMM_WORLD);

	//START CALCULATING

        for (count = 0; count <=maxCount; count++) {
            if (!rank) out_time = omp_get_wtime();
#if EILER
            #pragma omp parallel for
            for (i = 0; i < sendcounts[rank]; i++)
                prev[i] = points[i];//make pointers swap
            #pragma omp parallel for
            for (i = 1; i < sendcounts[rank] - 1; i++)
                points[i] = prev[i] + data.Sigma*data.deltaT*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) * divider;
#endif /*EILER*/
#if RUNGE_KUTTA
#pragma omp parallel for
            for (i = 0; i < sendcounts[rank]; i++)
                prev[i] = points[i];// make pointers swap
                                    // DERIVATIVE CALCULATING
                                    // 1
#pragma omp parallel for
            for (i = 1; i < sendcounts[rank] - 1; i++)
                // 1 - write, 4 - reads, store 1/pow(step, 2) in separate variable
                k1[i] = data.deltaT*data.Sigma*(prev[i - 1] - 2 * prev[i] + prev[i + 1]) * divider;
#pragma omp parallel for
            for (i = 0; i < sendcounts[rank]; i++)
                medium[i] = prev[i] + k1[i] * 0.5;
            // 2
#pragma omp parallel for
            for (i = 1; i < sendcounts[rank] - 1; i++)
                k2[i] = data.deltaT*data.Sigma*(medium[i - 1] - 2 * medium[i] + medium[i + 1]) * divider;
#pragma omp parallel for
            for (i = 0; i < sendcounts[rank]; i++)
                medium[i] = prev[i] + k2[i] * 0.5;
            // 3
#pragma omp parallel for
            for (i = 1; i < sendcounts[rank] - 1; i++)
                k3[i] = data.deltaT*data.Sigma*(medium[i - 1] - 2 * medium[i] + medium[i + 1]) * divider;
#pragma omp parallel for
            for (i = 0; i < sendcounts[rank]; i++)
                medium[i] = prev[i] + k2[i];
            // 4
#pragma omp parallel for
            for (i = 1; i < sendcounts[rank] - 1; i++)
                k4[i] = data.deltaT*data.Sigma*(medium[i - 1] - 2 * medium[i] + medium[i + 1]) * divider;
#pragma omp parallel for
            for (i = 1; i < sendcounts[rank] - 1; i++) {
                derivative[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * 0.1666666666;
            }
            // NEXT POINT CALCULATING
#pragma omp parallel for
            for (i = 1; i < sendcounts[rank] - 1; i++)
                points[i] = prev[i] + derivative[i];

#endif /*RUNGE_KUTTA*/
#if MATRIX_EILER
            #pragma omp parallel for
            for (i = 0; i < sendcounts[rank]; i++)
                prev[i] = points[i];//make pointers swap
            //if(rank==1)
            CsrMult(Values, ColumnNum, LineFirst, prev, sendcounts[rank], points);
#endif
#if MATRIX_RUNGE_KUTTA
#pragma omp parallel for
                for (i = 0; i < sendcounts[rank]; i++)
                    prev[i] = points[i];//make pointers swap
                                        //DERIVATIVE CALCULATING
                                        //1
                CsrMult(Values, ColumnNum, LineFirst, prev, sendcounts[rank], k1);
                SendBorderValues(k1, sendcounts, rank, proc_num);
#pragma omp parallel for
                for (i = 0; i < sendcounts[rank]; i++)
                    medium[i] = (prev[i] + k1[i]) * 0.5;
                //2
                CsrMult(Values, ColumnNum, LineFirst, medium, sendcounts[rank], k2);
                SendBorderValues(k2, sendcounts, rank, proc_num);
#pragma omp parallel for
                for (i = 0; i < sendcounts[rank]; i++)
                    medium[i] = (prev[i] + k2[i]) * 0.5;
                //3
                CsrMult(Values, ColumnNum, LineFirst, medium, sendcounts[rank], k3);
                SendBorderValues(k3, sendcounts, rank, proc_num);
#pragma omp parallel for
                for (i = 0; i < sendcounts[rank]; i++)
                    medium[i] = k2[i];
                //4
                CsrMult(Values, ColumnNum, LineFirst, medium, sendcounts[rank], k4);
                SendBorderValues(k4, sendcounts, rank, proc_num);
#pragma omp parallel for
                for (i = 1; i < sendcounts[rank] - 1; i++) {
                    derivative[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * 0.1666666666;
                }
                //NEXT POINT CALCULATING
#pragma omp parallel for
                for (i = 1; i < sendcounts[rank] - 1; i++)
                    points[i] = derivative[i];
#endif
#if YAKOBI
                flag = 1;
#pragma omp parallel for
                for (i = 0; i < sendcounts[rank]; i++)
                    prev[i] = points[i];//make pointers swap
                while (flag) {
                    flag = 0;
                    //TRY TO CHECK DISCREPANCY DELTA
                    CsrMult(Values, ColumnNum, LineFirst, points, sendcounts[rank], discrepancy);
                    SendBorderValues(discrepancy, sendcounts, rank, proc_num);
#pragma omp parallel for
                    for (i = 0; i < sendcounts[rank]; i++) {
                        discrepancy[i] -= prev[i];
                        discrepancy[i] = fabs(discrepancy[i]);
                        if (discrepancy[i] >= delta) {
                            flag = 1;
                        }
                    }
                    if (!flag)
                        break;

                    #pragma omp parallel for
                    for (i = 0; i < sendcounts[rank]; i++)
                        discrepancy[i] = points[i];

                    for (i = 0; i < sendcounts[rank]; i++) {
                        sum = 0;
                        for (int j = 0; j < sendcounts[rank]; j++)
                            if (j != i)
                                sum += CrsAccess(Values, ColumnNum, LineFirst, i, j, sendcounts[rank])*discrepancy[j];
                        points[i] = MainDiag[i] * (prev[i] - sum);
                    }
                }
#endif
            // Sending and recieve border points
            SendBorderValues(points, sendcounts, rank, proc_num);

            // delete output time from time calculating
            if(!rank) time += omp_get_wtime() - out_time;
            if (count%OutCount == 0 && count) {
                if (!rank) {

                    //Can't send-receive from 0 process to himself
                    for (int i = 0; i < elements_per_process; i++)
                        print_array[i] = points[i];
                    for (int sender_number = 1; sender_number<proc_num; sender_number++)
                        MPI_Recv(&print_array[elements_per_process*sender_number],
                            /* If we Recv from last process - getting remainder(ostatok) of the division */
                            sender_number == proc_num - 1 ?
                            elements_per_process + data.NumPoints % proc_num :
                            elements_per_process
                            /* Other way getting only elements_per_process count */
                            , MPI_DOUBLE, sender_number, sender_number, MPI_COMM_WORLD, &Status);
                    OutputArr(inputFile, print_array, data.NumPoints);//if ERROR will be after this line, numbers of lines and t parameter will be wrong
                    OutTime += data.deltaOut;
                    OutputCurrentTime(inputFile, OutTime);
                }
                // Send points from not root process to root for Outpt
                else if (rank == proc_num - 1)
                    MPI_Send(points + 1, elements_per_process + data.NumPoints % proc_num, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
                else
                    MPI_Send(points + 1, elements_per_process, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
            }
        }

#if 0
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
#endif /*0*/
    if(!rank)printf("\ntime: %f%\n", time);
#endif /*!START_STRING*/
	MPI_Finalize();

	return 0;

}
