#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 係数行列AをLU分解(成分はA自身に保存)し、連立方程式を解く //
void solver_LU_decomposition(double *A, double *x, double *b, int DoF_free)
{
    double *y; // 補助変数

    if ((y = (double *)calloc(DoF_free, sizeof(double))) == NULL)
    {
        printf("Error:Memory is not enough\n");
        exit(-1);
    }

    // AをLU分解(Lの対角成分を1とする) //
    for (int j = 0; j < DoF_free; j++)
    {
        for (int i = 0; i < j + 1; i++)
        {
            for (int k = 0; k < i; k++)
                A[DoF_free * i + j] -= A[DoF_free * i + k] * A[DoF_free * k + j];
        }
        for (int i = j + 1; i < DoF_free; i++)
        {
            for (int k = 0; k < j; k++)
                A[DoF_free * i + j] -= A[DoF_free * i + k] * A[DoF_free * k + j];
            A[DoF_free * i + j] /= A[DoF_free * j + j];
        }
    }

    // 方程式 Ly=bを解く //
    for (int i = 0; i < DoF_free; i++)
    {
        y[i] = b[i];
        for (int j = 0; j < i; j++)
            y[i] -= A[DoF_free * i + j] * y[j];
    }
    // 方程式 Ux=yを解く //
    for (int i = DoF_free - 1; -1 < i; i--)
    {
        x[i] = y[i];
        for (int j = i + 1; j < DoF_free; j++)
            x[i] -= A[DoF_free * i + j] * x[j];
        x[i] /= A[DoF_free * i + i];
    }

    free(y);
}
