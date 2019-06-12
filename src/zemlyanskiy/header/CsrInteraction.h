#pragma once
static inline void MultiCSR(double *arr, double *composition, double *csr_values, int *first_line, int *columns, const int size)
{
#pragma omp parallel for
  for (int i = 0; i < size; i++)
  {
    composition[i] = 0;
    for (int j = first_line[i]; j < first_line[i + 1]; j++)
      composition[i] += csr_values[j] * arr[columns[j]];
  }
}

static inline double GetCsrDiag(double *csr_values, int *columns, int *first_line, int iRaw, int jRaw, int Length)
{
  if (iRaw >= Length)
    return 0;
  int firstNum = first_line[iRaw];
  int secondNum = first_line[iRaw + 1];
  for (int j = firstNum; j < secondNum; j++)
  {
    if (columns[j] == jRaw)
    {
      return csr_values[j];
    }
  }

  return 0;
}
