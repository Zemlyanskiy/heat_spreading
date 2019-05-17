/*
 * Input/ouput file handling
 */

#ifndef __IO_PARSER
#define __IO_PARSER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Input functions
struct input {
    double t;// Start time (t)
    double T;// End time (T)
    double deltaT;// Gap between measurements (dt)
    double Xmin;
    double Xmax;
    double Ymin;
    double Ymax;
    double Zmin;
    double Zmax;
    double Sigma; // The coefficient of quation
    double deltaOut;// Delta of output
    unsigned Lx, Ly, Lz;
    double* arr;
};

inline struct input ReadInput(char* PATH) {
    FILE *file;
    struct input temp = { .t = 0,.T = 0,.deltaT = 0,
                          .Xmin = 0,.Xmax = 0, .Lx = 0,
                          .Ymin = 0,.Ymax = 0, .Ly = 0,
                          .Zmin = 0,.Zmax = 0, .Lz = 0,
                          .Sigma = 0,.deltaOut = 0, .arr = NULL };
    char* token;
    int pos = 0;
    double buf = 0;
    char* string = (char*)calloc(10000000, sizeof(char));

    fopen_s(&file, PATH, "r");
    if (file == NULL) {
        printf("Can`t find file\n");
    }
    else {
        for (int i = 0; fgets(string, 10000000 * sizeof(char), file) != NULL; i++)
            switch (i) {
            case 0:
            {
                sscanf_s(string, "t=%lf", &temp.t);
                break;
            }
            case 1:
            {
                sscanf_s(string, "T=%lf", &temp.T);
                break;
            }
            case 2:
            {
                sscanf_s(string, "deltaT=%lf", &temp.deltaT);
                break;
            }
            case 3:
            {
                sscanf_s(string, "XMin=%lf", &temp.Xmin);
                break;
            }
            case 4:
            {
                sscanf_s(string, "XMax=%lf", &temp.Xmax);
                break;
            }
            case 5:
            {
                sscanf_s(string, "YMin=%lf", &temp.Ymin);
                break;
            }
            case 6:
            {
                sscanf_s(string, "YMax=%lf", &temp.Ymax);
                break;
            }
            case 7:
            {
                sscanf_s(string, "ZMin=%lf", &temp.Zmin);
                break;
            }
            case 8:
            {
                sscanf_s(string, "ZMax=%lf", &temp.Zmax);
                break;
            }
            case 9:
            {
                sscanf_s(string, "Lx=%d", &temp.Lx);
                break;
            }
            case 10:
            {
                sscanf_s(string, "Ly=%d", &temp.Ly);
                break;
            }
            case 11:
            {
                sscanf_s(string, "Lz=%d", &temp.Lz);
                break;
            }
            case 12:
            {
                sscanf_s(string, "Sigma=%lf", &temp.Sigma);
                break;
            }
            case 13:
            {
                sscanf_s(string, "deltaOut=%lf", &temp.deltaOut);
                break;
            }
            default:
#if START_STRING
                printf("Mode: start string, Problem: file allready contain start string\n");
#endif
                token = strtok(string, " ");
                if (temp.arr == NULL)
                    temp.arr = (double*)calloc(temp.Lx * temp.Ly * temp.Lz, sizeof(double));
                while (token != NULL)
                {
                    buf = strtof(token, NULL);
                    temp.arr[pos] = buf;
                    token = strtok(NULL, " ");
                    pos++;
                }
                pos = 0;
            }
    }

    fclose(file);
    return temp;
}

// Output functions
inline void OutputArr(char* PATH, double* arr, int length) {
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

#endif /* __IO_PARSER */