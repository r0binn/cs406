#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
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
// Function to compare two pairs based on the second value (non-zero count)
bool comparePairs(const pair<int, int> &a, const pair<int, int> &b)
{
    return a.second < b.second;
}

// Function to count non-zero items per column
void countNonZeroPerColumn(int *rptrs, int *columns, int N, int m, int *colCounts)
{
    for (int row = 0; row < N; ++row)
    {
        for (int idx = rptrs[row]; idx < rptrs[row + 1]; ++idx)
        {
            int col = columns[idx];
            colCounts[col]++;
        }
    }
}

// Function to create a mapping from old column indices to new sorted column indices
void createColumnMapping(int *colCounts, int N, int *colMap, int *sortedCols)
{
    pair<int, int> *colWithCount = (pair<int, int> *)malloc(N * sizeof(pair<int, int>));
    for (int col = 0; col < N; ++col)
    {
        colWithCount[col] = {col, colCounts[col]};
    }

    sort(colWithCount, colWithCount + N, comparePairs);

    for (int i = 0; i < N; ++i)
    {
        sortedCols[i] = colWithCount[i].first;
        colMap[colWithCount[i].first] = i;
    }

    free(colWithCount);
}

// Function to rearrange CSR arrays according to the sorted columns
void rearrangeCSR(int *rptrs, int *columns, double *rvals, int *sortedCols, int *colMap, int *newRptrs, int *newColumns,
                  double *newRvals, int N, int m)
{
    int *tempCounts = (int *)calloc(N, sizeof(int));

    for (int row = 0; row < N; ++row)
    {
        for (int idx = rptrs[row]; idx < rptrs[row + 1]; ++idx)
        {
            int col = columns[idx];
            int newCol = colMap[col];
            newColumns[newRptrs[row] + tempCounts[row]] = newCol;
            newRvals[newRptrs[row] + tempCounts[row]] = rvals[idx];
            tempCounts[row]++;
        }
    }

    free(tempCounts);
}

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
    if (nzeros == 0)
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
    unsigned long int maxG = (1ULL << (N - 1));
    for (int g = 1; g < maxG; g++)
    {

        int grayNow = (g >> 1) ^ g;
        int grayPrev = ((g - 1) >> 1) ^ (g - 1);
        int j = __builtin_ctz(grayNow ^ grayPrev);
        int s = ((grayNow >> (j)) & 1) == 1 ? 1 : -1;
        int sign = (g & 1) == 0 ? 1 : -1;

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
            p += (double)(sign * prod);
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
    unsigned long int maxG = (1ULL << (N - 1));

    for (int g = 1; g < maxG; g++)
    {

        int grayNow = (g >> 1) ^ g;
        int grayPrev = ((g - 1) >> 1) ^ (g - 1);
        int j = __builtin_ctz(grayNow ^ grayPrev);
        int s = ((grayNow >> (j)) & 1) == 1 ? 1 : -1;
        int sign = (g & 1) == 0 ? 1 : -1;

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
        p += (double)(sign * prod);
    }
    free(X);
    return p * (4 * (N % 2) - 2);
}

int main(int argc, char *argv[])
{
    char fileName[80] = "out3.txt";
    int n; // Size of square matrix
    int m; // Number of non-zero items
    int *rowPtr = nullptr;
    int *colIdx = nullptr;
    double *values = nullptr;
    int numThreads = atoi(argv[1]);
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
        printf("The result is : %d\n", 0);
        exit(0);
    }
    for (size_t i = 1; i < N; i++)
    {
        if (rptrs[i - 1] >= rptrs[i])
        {
            printf("The result is : %d\n", 0);
            exit(0);
        }
    }
    // Step 1: Count non-zero items per column
    int *colCounts = (int *)calloc(N, sizeof(int));
    countNonZeroPerColumn(rptrs, columns, N, m, colCounts);

    // Step 2: Sort columns by the number of non-zero items
    int *sortedCols = (int *)malloc(N * sizeof(int));
    int *colMap = (int *)malloc(N * sizeof(int));
    createColumnMapping(colCounts, N, colMap, sortedCols);

    // Step 3: Rearrange CSR arrays according to the sorted columns
    int *newRptrs = (int *)malloc((N + 1) * sizeof(int));
    int *newColumns = (int *)malloc(m * sizeof(int));
    double *newRvals = (double *)malloc(m * sizeof(double));
    memset(newRptrs, 0, (N + 1) * sizeof(int));

    for (int i = 0; i < N; ++i)
    {
        newRptrs[i + 1] = rptrs[i + 1] - rptrs[i];
    }
    for (int i = 1; i <= N; ++i)
    {
        newRptrs[i] += newRptrs[i - 1];
    }

    rearrangeCSR(rptrs, columns, rvals, sortedCols, colMap, newRptrs, newColumns, newRvals, N, m);

    free(rvals);
    free(columns);
    free(rptrs);
    rvals = newRvals;
    columns = newColumns;
    rptrs = newRptrs;
    unsigned long int maxG = (1ULL << (N - 1));
    int *cptrs;
    int *rows;
    double *cvals;
    convertCRStoCCS(rptrs, columns, rvals, N, m, &cptrs, &rows, &cvals);
    double startS = omp_get_wtime();
    double trueResult = seqSpaRyserNzero(rptrs, columns, rvals, cptrs, rows, cvals);
    double endS = omp_get_wtime();
    printf("Sequential result is : %e seconds %lf\n", trueResult, endS - startS);
    omp_set_num_threads(numThreads);
    unsigned long long int chunksize = (1ULL << (N - 1));
    printf("Chunksize %lu\n", chunksize);

    chunksize = floor((chunksize - 1) / numThreads);
    chunksize = (chunksize > 1) ? chunksize : 1;
    printf("Chunksize %llu\n", chunksize);
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
    double tempResult;
    printf("Parallel has started\n");
    double startT = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, chunksize) reduction(+ : p) proc_bind(spread)
    for (unsigned long long int g = 1; g < maxG; g++)
    {
        int thread_id = omp_get_thread_num();
        int nzeros = 0;
        double *myX = (double *)malloc((N) * sizeof(double));
        for (size_t i = 0; i < N; i++)
        {
            myX[i] = X[i];
        }
        int myStart = g;
        int prevGray = (myStart - 1) ^ ((myStart - 1) >> 1);

        double myP = 0;
        int colIdx = 0;
        if (g != 1)
        {
            while (prevGray != 0)
            {
                if (prevGray & 1 == 1)
                {
                    for (int ptr = cptrs[colIdx]; ptr < cptrs[colIdx + 1]; ptr++)
                    {
                        int row = rows[ptr];
                        double val = cvals[ptr];
                        myX[row] = myX[row] + val;
                    }
                }
                prevGray = prevGray >> 1;
                colIdx += 1;
            }
        }
        for (size_t i = 0; i < N; i++)
        {
            if (myX[i] == 0)
            {
                nzeros += 1;
            }
        }
        unsigned long long int myEnd = g + chunksize;
        myEnd = (myEnd < maxG) ? myEnd : maxG;
        for (; g < myEnd; g++)
        {

            int grayNow = (g >> 1) ^ g;
            int grayPrev = ((g - 1) >> 1) ^ (g - 1);
            int j = __builtin_ctz(grayNow ^ grayPrev);
            int s = ((grayNow >> (j)) & 1) == 1 ? 1 : -1;
            int prodsign = (g & 1) == 0 ? 1 : -1;

            for (int ptr = cptrs[j]; ptr < cptrs[j + 1]; ptr++)
            {
                int row = rows[ptr];
                double val = cvals[ptr];
                if (myX[row] == 0)
                {
                    nzeros -= 1;
                }
                myX[row] = myX[row] + (s * val);
                if (myX[row] == 0)
                {
                    nzeros += 1;
                }
            }
            if (nzeros == 0)
            {
                double prod = 1.0;
                for (size_t i = 0; i < N; i++)
                {
                    prod *= myX[i];
                }
                myP += (double)(prodsign * prod);
            }
        }
        free(myX);
        p += myP;
    }

    tempResult = p * (4 * (N % 2) - 2);
    double endT = omp_get_wtime();
    printf("parallel ends in %f seconds with %d threads ", endT - startT, numThreads);
    printf("The given result is : %e\n", tempResult);
    printf("Error :%lf\n", abs(trueResult - tempResult) / trueResult);
    printf("Speedup :%lf\n", (endS - startS) / (endT - startT));
    return 0;
}
