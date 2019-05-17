/*
 * Compressed sparce row API
 */

#ifndef __CSR_API
#define __CSR_API

inline double CsrAccess(double* Values, int* ColumnNum, int* LineFirst, int i, int j, int length) {
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

inline void CsrMult(double* Values, int* ColumnNum, int* LineFirst, double * Array, const int size, double* result, unsigned rank)
{
    int i;
#pragma omp parallel for
    for (i = 0; i < size; i++) {
        result[i] = 0;
        for (int j = LineFirst[i]; j < LineFirst[i + 1]; j++) {
            result[i] += Values[j] * Array[ColumnNum[j]];
        }
    }
}

#endif /* __CSR_API */