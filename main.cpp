#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;
typedef struct sparseMatrix
{
    int i;
    int j;
    double val;
} sparseMatrix;

sparseMatrix *X;
int N;
int m;
void convertCRStoCCS(int *rptrs, int *columns, double *rvals, int NumRows, int NumElements, int **cptrs, int **rows,
                     double **cvals)
{
    // Allocate memory for CCS arrays
    *cptrs = (int *)malloc((NumRows + 1) * sizeof(int));
    *rows = (int *)malloc(NumElements * sizeof(int));
    *cvals = (double *)malloc(NumElements * sizeof(double));

    // Count non-zero elements in each column
    int *colCount = (int *)calloc(NumRows, sizeof(int));
    for (int i = 0; i < NumElements; i++)
    {
        colCount[columns[i]]++;
    }

    // Compute cptrs
    (*cptrs)[0] = 0;
    for (int i = 0; i < NumRows; i++)
    {
        (*cptrs)[i + 1] = (*cptrs)[i] + colCount[i];
        colCount[i] = 0; // Reset for reuse
    }

    // Fill rows and cvals
    for (int i = 0; i < NumRows; i++)
    {
        for (int j = rptrs[i]; j < rptrs[i + 1]; j++)
        {
            int col = columns[j];
            int index = (*cptrs)[col] + colCount[col];
            (*rows)[index] = i;
            (*cvals)[index] = rvals[j];
            colCount[col]++;
        }
    }

    free(colCount);
}
int loadDataSparse(char *fileName, sparseMatrix **data, int *size1)
{
    FILE *myfile;
    if ((myfile = fopen(fileName, "r")) == NULL)
    {
        printf("Error: file cannot be found\n");
        exit(1);
    }

    int s1, numel, result;
    if ((result = fscanf(myfile, "%d %d", &s1, &numel)) <= 0)
    {
        printf("error while reading file %d\n", result);
        exit(1);
    }
    *data = (sparseMatrix *)malloc(numel * sizeof(sparseMatrix));
    for (int i = 0; i < numel; i++)
    {
        int tempi, tempj;
        double tempval;
        if ((result = fscanf(myfile, "%d %d %lf", &tempi, &tempj, &tempval)) <= 0)
        {
            printf("error while reading file - %d\n", result);
            exit(1);
        }

        (*data)[i].i = tempi;
        (*data)[i].j = tempj;
        (*data)[i].val = tempval;
    }

    fclose(myfile);

    *size1 = s1;
    return numel;
}
double seqSpaRyserNzero(int *rptrs, int *columns, double *rvals, int *cptrs, int *rows, double *cvals)
{
    int nzeros = 0;
    double *X = (double *)malloc((N) * sizeof(double));
    for (size_t i = 0; i < N; i++)
    {
        double sum = 0;
        for (size_t ptr = rptrs[i]; ptr < rptrs[i + 1]; ptr++)
        {
            sum += rvals[ptr];
        }
        if (columns[rptrs[i + 1] - 1] == N - 1)
        {
            X[i] = rvals[rptrs[i + 1] - 1] - (sum / 2);
        }
        else
        {
            X[i] = 0 - (sum / 2);
        }

        if (X[i] == 0)
        {
            nzeros += 1;
        }
    }
    double p;
    if (nzeros <= 0)
    {
        p = 1.0;
        for (size_t i = 0; i < N; i++)
        {
            p *= X[i];
        }
    }
    else
    {
        p = 0.0;
    }

    for (int g = 1; g < (1 << (N - 1)); g++)
    {

        int y = (g >> 1) ^ g;                   // gray-code order
        int yy = ((g - 1) >> 1) ^ (g - 1);      // i-1's gray-code order
        int j = __builtin_ctz(y ^ yy);          // get the changing bit
        int s = ((y >> (j)) & 1) == 1 ? 1 : -1; // find changing bit
        int prodsign = (g & 1) == 0 ? 1 : -1;   // get the prodsign

        for (int ptr = cptrs[j]; ptr < cptrs[j + 1]; ptr++)
        {
            int row = rows[ptr];
            double val = cvals[ptr];
            if (X[row] == 0)
            {
                nzeros -= 1;
            }
            X[row] = X[row] + (s * val);
            if (X[row] == 0)
            {
                nzeros += 1;
            }
        }
        if (nzeros == 0)
        {
            double prod = 1.0;
            for (size_t i = 0; i < N; i++)
            {
                prod *= X[i];
            }
            p += (double)(prodsign * prod);
        }
    }
    free(X);
    return p * (4 * (N % 2) - 2);
}
double seqSpaRyser(int *rptrs, int *columns, double *rvals, int *cptrs, int *rows, double *cvals)
{
    double *X = (double *)malloc((N) * sizeof(double));
    for (size_t i = 0; i < N; i++)
    {
        double sum = 0;
        for (size_t ptr = rptrs[i]; ptr < rptrs[i + 1]; ptr++)
        {
            sum += rvals[ptr];
        }
        if (columns[rptrs[i + 1] - 1] == N - 1)
        {
            X[i] = rvals[rptrs[i + 1] - 1] - (sum / 2);
        }
        else
        {
            X[i] = 0 - (sum / 2);
        }
    }
    double p = 1.0;
    for (size_t i = 0; i < N; i++)
    {
        p *= X[i];
    }
    for (int g = 1; g < (1 << (N - 1)); g++)
    {

        int y = (g >> 1) ^ g;                   // gray-code order
        int yy = ((g - 1) >> 1) ^ (g - 1);      // i-1's gray-code order
        int j = __builtin_ctz(y ^ yy);          // get the changing bit
        int s = ((y >> (j)) & 1) == 1 ? 1 : -1; // find changing bit
        int prodsign = (g & 1) == 0 ? 1 : -1;   // get the prodsign

        for (int ptr = cptrs[j]; ptr < cptrs[j + 1]; ptr++)
        {
            int row = rows[ptr];
            double val = cvals[ptr];
            X[row] = X[row] + (s * val);
        }

        double prod = 1.0;
        for (size_t i = 0; i < N; i++)
        {
            prod *= X[i];
        }
        p += (double)(prodsign * prod);
    }
    free(X);
    return p * (4 * (N % 2) - 2);
}
int main()
{
    char fileName[80] = "out3.txt";
    int n; // Size of square matrix
    int m; // Number of non-zero items
    int *rowPtr = nullptr;
    int *colIdx = nullptr;
    double *values = nullptr;

    m = loadDataSparse(fileName, &X, &N);
    int *rptrs = (int *)malloc((N + 1) * sizeof(int));
    int *columns = (int *)malloc(m * sizeof(int));
    double *rvals = (double *)malloc(m * sizeof(double));
    memset(rptrs, 0, (N + 1) * sizeof(int));
    for (int i = 0; i < m; i++)
    {
        int rowid = X[i].i;
        if (rowid < 0 || rowid >= N)
        {
            printf("problem in X, quitting - %d\n", rowid);
            exit(1);
        }
        rptrs[rowid + 1]++;
    }
    for (int i = 1; i <= N; i++)
    {
        rptrs[i] += rptrs[i - 1];
    }
    for (int i = 0; i < m; i++)
    {
        int rowid = X[i].i;
        int index = rptrs[rowid];

        columns[index] = X[i].j;
        rvals[index] = X[i].val;

        rptrs[rowid] = rptrs[rowid] + 1;
    }

    for (int i = N; i > 0; i--)
    {
        rptrs[i] = rptrs[i - 1];
    }
    rptrs[0] = 0;
    if (rptrs[N - 1] == m)
    {
        printf("The result is : %lf\n", 0);
        exit(0);
    }

    for (size_t i = 1; i < N; i++)
    {
        if (rptrs[i - 1] >= rptrs[i])
        {
            printf("The result is : %lf\n", 0);
            exit(0);
        }
    }
    int *cptrs;
    int *rows;
    double *cvals;
    convertCRStoCCS(rptrs, columns, rvals, N, m, &cptrs, &rows, &cvals);
    double trueResult = seqSpaRyserNzero(rptrs, columns, rvals, cptrs, rows, cvals);
    printf("The result is : %lf \n", trueResult);
    return 0;
}
