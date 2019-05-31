#include "common.h"


double CalcFunc(double Xcoord, double Ycoord, double Zcoord) {
    if (Xcoord >= -0.5 && Xcoord <= 0.5 &&
        Ycoord >= -0.5 && Ycoord <= 0.5 &&
        Zcoord >= -0.5 && Zcoord <= 0.5) {
        return cos(Xcoord*3.141592) + cos(Ycoord*3.141592) + cos(Zcoord*3.141592);
    }
    return 0;
}

double CalcFuncDebug(double x, double y, double z) {
    static unsigned counter = 0;
    return counter++;
}

int main(int argc, char* argv[])
{
    // MPI initialization ---
    unsigned rank, proc_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        // Processing command line arguments ---
        char* inputFile = NULL;
        unsigned number_of_threads = 1;
        {
            if (argc > 1)
                inputFile = argv[1];
            if (argc > 2)
                number_of_threads = strtol(argv[2], NULL, 0);
            if (argc > 3) {
                printf("Incorrect number of arguments");
                return 4;
            }
        }

        // Read input data ---
        struct input data;
        data = ReadInput(inputFile);

        double Xstep = (data.Xmax - data.Xmin) / data.Lx;
        double Ystep = (data.Ymax - data.Ymin) / data.Ly;
        double Zstep = (data.Zmax - data.Zmin) / data.Lz;

        double posX = data.Xmin;
        double posY = data.Ymin;
        double posZ = data.Zmin;

        double* points;
        points = (double*)calloc(data.Lx * data.Ly * data.Lz, sizeof(double));

        for (unsigned i = 0; i < data.Lx; i++) {
            posX += Xstep;
            posY = data.Ymin;
            for (unsigned j = 0; j < data.Ly; j++) {
                posY += Ystep;
                posZ = data.Zmin;
                for (unsigned k = 0; k < data.Lz; k++) {
                    posZ += Zstep;
                    Get(points, i, j, k, data.Lx, data.Ly, data.Lz) = CalcFunc(posX, posY, posZ);
                }
            }
        }

        OutputArr(inputFile, points, data.Lx * data.Ly * data.Lz);
        OutputArr(inputFile, points, data.Lx * data.Ly * data.Lz);
    }
    MPI_Finalize();
    return 0;
}