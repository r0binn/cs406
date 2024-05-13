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
int NumCols, NumRows;
int NumElements;
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

    printf("The number of rows, columns, and nonzeros are %d, and %d respectively\n", s1, numel);

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

int main()
{
    char fileName[80] = "out1.txt";
    int n; // Size of square matrix
    int m; // Number of non-zero items
    int *rowPtr = nullptr;
    int *colIdx = nullptr;
    double *values = nullptr;

    NumElements = loadDataSparse(fileName, &X, &NumRows);
    int *crs_ptrs = (int *)malloc((NumRows + 1) * sizeof(int));
    int *crs_colids = (int *)malloc(NumElements * sizeof(int));
    double *crs_values = (double *)malloc(NumElements * sizeof(double));
    memset(crs_ptrs, 0, (NumRows + 1) * sizeof(int));
    for (int i = 0; i < NumElements; i++)
    {
        int rowid = X[i].i;
        if (rowid < 0 || rowid >= NumRows)
        {
            printf("problem in X, quitting - %d\n", rowid);
            exit(1);
        }
        crs_ptrs[rowid + 1]++;
    }
    for (int i = 1; i <= NumRows; i++)
    {
        crs_ptrs[i] += crs_ptrs[i - 1];
    }
    for (int i = 0; i < NumElements; i++)
    {
        int rowid = X[i].i;
        int index = crs_ptrs[rowid];

        crs_colids[index] = X[i].j;
        crs_values[index] = X[i].val;

        crs_ptrs[rowid] = crs_ptrs[rowid] + 1;
    }

    for (int i = NumRows; i > 0; i--)
    {
        crs_ptrs[i] = crs_ptrs[i - 1];
    }
    crs_ptrs[0] = 0;
    printf("CRS Pointers: ");
    for (int i = 0; i < NumRows; i++)
    {
        printf("%d ", crs_ptrs[i]);
    }
    printf("\n");
    printf("CRS Colids: ");
    for (int i = 0; i < NumElements; i++)
    {
        printf("%d ", crs_colids[i]);
    }
    printf("\n");
    printf("CRS Values: ");
    for (int i = 0; i < NumElements; i++)
    {
        printf("%lf ", crs_values[i]);
    }
    printf("\n");
    int nzeros = 0;
    double *X = (double *)malloc((NumRows) * sizeof(double));
    for (size_t i = 0; i < NumRows; i++)
    {
        int sum = 0;
        for (size_t j = crs_ptrs[i]; j < crs_ptrs[i + 1]; j++)
        {
            sum += crs_values[j];
        }
        X[i] = crs_values[crs_ptrs[i + 1] - 1] - (sum / 2);
        if (X[i] == 0)
            nzeros += 1;
    }

    printf("Nzeros: %d\n", nzeros);
    double p;
    if (nzeros > 0)
    {
        p = 1.0;
        for (size_t i = 0; i < NumRows; i++)
        {
            p *= X[i];
        }
    }
    else
    {
        p = 0.0;
    }
    for (long long int i = 1; i < (1 << (NumRows - 1)); i++)
    {
        int y = (i >> 1) ^ i;                 // gray-code order
        int yPrev = ((i - 1) >> 1) ^ (i - 1); // i-1's gray-code order
        int s = (int)__builtin_ctz(y ^ yPrev) + 1;

        // int prodsign = (i & 1) == 0 ? 1 : -1;   // get the prodsign
        // double dd = 1.0;

        // #pragma omp simd reduction(* : dd)
        //         for (int j = 0; j < N; j++)
        //         {
        //             x_specul[j] += (double)(s * MT[z][j]);
        //             dd *= x_specul[j];
        //         }
        //         p += (double)(prodsign * dd);
        //     }
        return 0;
    }
}