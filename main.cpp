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
    printf("This number should be equal to the number of nonzeros %d\n", crs_ptrs[NumRows]);
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
    printf("CRS Pointers:\n");
    for (int i = 0; i <= NumElements; i++)
    {
        printf("%d ", crs_ptrs[i]);
    }

    printf("CRS Colids:\n");
    for (int i = 0; i <= NumElements; i++)
    {
        printf("%d ", crs_colids[i]);
    }
    printf("CRS Values:\n");
    for (int i = 0; i <= NumElements; i++)
    {
        printf("%lf ", crs_values[i]);
    }
    return 0;
}
