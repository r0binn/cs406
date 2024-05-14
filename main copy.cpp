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
    char fileName[80] = "out3.txt";
    int n; // Size of square matrix
    int m; // Number of non-zero items
    int *rowPtr = nullptr;
    int *colIdx = nullptr;
    double *values = nullptr;

    NumElements = loadDataSparse(fileName, &X, &NumRows);
    int *rptrs = (int *)malloc((NumRows + 1) * sizeof(int));
    int *columns = (int *)malloc(NumElements * sizeof(int));
    double *rvals = (double *)malloc(NumElements * sizeof(double));
    memset(rptrs, 0, (NumRows + 1) * sizeof(int));
    for (int i = 0; i < NumElements; i++)
    {
        int rowid = X[i].i;
        if (rowid < 0 || rowid >= NumRows)
        {
            printf("problem in X, quitting - %d\n", rowid);
            exit(1);
        }
        rptrs[rowid + 1]++;
    }
    for (int i = 1; i <= NumRows; i++)
    {
        rptrs[i] += rptrs[i - 1];
    }
    for (int i = 0; i < NumElements; i++)
    {
        int rowid = X[i].i;
        int index = rptrs[rowid];

        columns[index] = X[i].j;
        rvals[index] = X[i].val;

        rptrs[rowid] = rptrs[rowid] + 1;
    }

    for (int i = NumRows; i > 0; i--)
    {
        rptrs[i] = rptrs[i - 1];
    }
    rptrs[0] = 0;

    int *cptrs;
    int *rows;
    double *cvals;
    convertCRStoCCS(rptrs, columns, rvals, NumRows, NumElements, &cptrs, &rows, &cvals);
    printf("rptrs: ");
    for (int i = 0; i < NumRows; i++)
    {
        printf("%d ", rptrs[i]);
    }
    printf("\n");
    printf("columns: ");
    for (int i = 0; i < NumElements; i++)
    {
        printf("%d ", columns[i]);
    }
    printf("\n");
    printf("rvals: ");
    for (int i = 0; i < NumElements; i++)
    {
        printf("%lf ", rvals[i]);
    }
    printf("\n");

    printf("cptrs: ");
    for (int i = 0; i < NumRows; i++)
    {
        printf("%d ", cptrs[i]);
    }
    printf("\n");
    printf("rows: ");
    for (int i = 0; i < NumElements; i++)
    {
        printf("%d ", rows[i]);
    }
    printf("\n");
    printf("cvals: ");
    for (int i = 0; i < NumElements; i++)
    {
        printf("%lf ", cvals[i]);
    }
    printf("\n");

    int nzeros = 0;
    double *X = (double *)malloc((NumRows) * sizeof(double));
    for (size_t i = 0; i < NumRows; i++)
    {
        double sum = 0;
        printf("%d %d %d\n", i, rptrs[i], rptrs[i + 1]);
        for (size_t ptr = rptrs[i]; ptr < rptrs[i + 1] - 1; ptr++)
        {
            sum += rvals[ptr];
        }
        X[i] = rvals[rptrs[i + 1] - 1] - (sum / 2);
        if (X[i] == 0)
            nzeros += 1;
    }
    printf("X: ");
    for (int i = 0; i < NumRows; i++)
    {
        printf("%lf ", X[i]);
    }
    printf("\n");
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
        X[0] = 0.5;
        X[1] = -0.5;
    }
    int tmp = 1 << (NumRows - 1);
    for (int g = 1; g < (1 << (NumRows - 1)); g++)
    {

        int y, yPrime, j, c, s;
        double pr_Sign;

        y = ((g) >> 1) ^ (g);
        yPrime = ((g - 1) >> 1) ^ (g - 1);
        j = __builtin_ctz(y ^ yPrime);
        s = ((y >> (j)) & 1) == 1 ? 1 : -1;
        pr_Sign = (g & 1) == 0 ? 1 : -1;
        printf("%d %d %d %d\n", g, y, j, s);
        printf("%d \n", rptrs[j + 1] - 1);
        for (int ptr = rptrs[j]; ptr < rptrs[j + 1] - 1; ptr++)
        {
            int row = rows[ptr];
            double val = cvals[ptr];
            if (X[row] == 0)
            {
                nzeros -= 1;
            }
            X[row] = X[row] + s * val;
            if (X[row] == 0)
            {
                nzeros += 1;
            }
        }
        if (nzeros == 0)
        {
            double prod = 1;
            for (size_t i = 0; i < NumRows; i++)
            {
                prod *= X[i];
            }
            p = p + (prod * pr_Sign);
        }
    }
    double result = p * (4 * (NumRows % 2) - 2);
    printf("The result is : %lf\n", result);
    return 0;
}
