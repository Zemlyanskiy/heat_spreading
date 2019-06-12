#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <mpi.h>
#include "../header/IOStream.h"
#include "../header/MainVaribles.h"
#include "../header/CsrInitialize.h"
#include "../header/CsrInteraction.h"
#include "../header/MpiWorking.h"
#include "../header/StartString.h"
#include "../header/Euler.h"
#include "../header/RungeKutta.h"
#include "../header/CsrEuler.h"
#include "../header/CsrRungeKutta.h"
#include "../header/UnsignedEuler.h"

int main(int argc, char *argv[])
{
    // MPI initialization
    unsigned rank, proc_num;

    file = argv[1];
    if (argv[3][0] == 'e' && argv[2][0] == 'n')
        prog_mode = 0;
    if (argv[3][0] == 'r' && argv[2][0] == 'n')
        prog_mode = 1;
    if (argv[3][0] == 'e' && argv[2][0] == 'm')
        prog_mode = 2;
    if (argv[3][0] == 'r' && argv[2][0] == 'm')
        prog_mode = 3;
    if (argv[3][0] == 'u' && argv[2][0] == 'm')
        prog_mode = 4;
    if (argv[3][0] == 's')
        prog_mode = 5;
    num_threads = atoi(argv[4]);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //Reading data
    data = Read(file, prog_mode);
    Initialize_MPI(data, rank, proc_num, border_size);
    InitializeMainVar(data);
    //Start-String
    if (prog_mode == 5)
    {
        StartString(data);
    }
    else
    {
        for (int i = 0; i < data.array.Lx * data.array.Ly * data.array.Lz; i++)
        {
            points.values[i] = data.array.values[i];
        }
    }
    //Euler
    if (prog_mode == 0)
    {
        Euler(data, rank, proc_num);
    }
    //Runge-Kutta
    if (prog_mode == 1)
    {
        Runge_Kutta(data, rank, proc_num);
    }
    //CsrEuler
    if (prog_mode == 2)
    {
        //Initialize CSR Matrix
        InitializeCsr(data, prog_mode);
        CsrEuler(data, rank, proc_num);
    }
    //CsrRungeKutta
    if (prog_mode == 3)
    {
        //Initialize CSR Matrix
        InitializeCsr(data, prog_mode);
        CsrRungeKutta(data, rank, proc_num);
    }
    if (prog_mode == 4)
    {
        //Initialize CSR Matrix
        InitializeCsr(data, prog_mode);
        UnsignedEuler(data, rank, proc_num);
    }
    if (!rank)
        printf("THE END\n");
    MPI_Finalize();
    return 0;
}