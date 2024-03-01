#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"tensor.h"
#include"type.h"
#include"scalar.h"
#include"vector.h"
#include"matrix.h"

extern Option option;
static const double identity_tensor[3][3] = {{1.0, 0.0, 0.0},
                                             {0.0, 1.0, 0.0},
                                             {0.0, 0.0, 1.0}};

void calculateTensorExponent(double (*tensor_out)[3], double (*tensor_in)[3])
#if 0
{
    calculateIsotropicTensorFunction(tensor_out, tensor_in,
                                     exp);
}
#else
{
    const double taylor_series_tolerance = 1.0e-4;
    const int taylor_series_max_iteration_count = 1000;

    double power_tensor[3][3], previous_power_tensor[3][3];
    double factorial;
    double nn;
    int n;
    int i, j, k;

    identify3x3Matrix(power_tensor);
    identify3x3Matrix(tensor_out);

    for (n = 1, nn = 1.0, factorial = 1.0;
         n < taylor_series_max_iteration_count;
         n++, nn += 1.0, factorial *= nn)
    {
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                previous_power_tensor[i][j] = power_tensor[i][j];

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            {
                double temp = 0.0;

                for (k = 0; k < 3; k++)
                    temp += previous_power_tensor[i][k] * tensor_in[k][j];

                power_tensor[i][j] = temp;
            }

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                tensor_out[i][j] += power_tensor[i][j] / factorial;

        if (calc_3x3matrix_norm(power_tensor)
            <= taylor_series_tolerance * factorial)
            break;
    }

    if (n == taylor_series_max_iteration_count)
    {
#if 0
        printError("Warning: tensor exponent calculation iteration not converged\n");
#endif

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                tensor_out[i][j] = nan("");
    }
}
#endif

/*
 * Calculate tensor logarithm
 */
void calculateTensorLogarithm(double (*tensor_out)[3], double (*tensor_in)[3])
#if 0
{
    calculateIsotropicTensorFunction(tensor_out, tensor_in,
                                     log);
}
#else
{
    const double taylor_series_tolerance = 1.0E-4;
    const int taylor_series_max_iteration_count = 1000;

    double tensor[3][3], square_tensor[3][3];
    double power_tensor[3][3], previous_power_tensor[3][3];
    double nn;
    int n;
    int i, j, k;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            square_tensor[i][j] = tensor_in[i][j];
    for (i = 0; i < 3; i++)
        square_tensor[i][i] += 1.0;
    inverse_mat3x3(option.dim, square_tensor,  tensor);
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            tensor[i][j] *= -2.0;
    for (i = 0; i < 3; i++)
        tensor[i][i] += 1.0;

    calculate3x3MatrixSquare(square_tensor, tensor);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            power_tensor[i][j] = tensor[i][j];

    zeroize3x3Matrix(tensor_out);

    for (n = 0, nn = 0.0;
         n < taylor_series_max_iteration_count;
         n++, nn += 1.0)
    {
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                tensor_out[i][j] += power_tensor[i][j] / (2.0 * nn + 1.0);

        if (calc_3x3matrix_norm(power_tensor)
            <= taylor_series_tolerance * (2.0 * nn + 1.0))
            break;

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                previous_power_tensor[i][j] = power_tensor[i][j];

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            {
                double temp = 0.0;

                for (k = 0; k < 3; k++)
                    temp += previous_power_tensor[i][k] * square_tensor[k][j];

                power_tensor[i][j] = temp;
            }
    }

    if (n == taylor_series_max_iteration_count)
    {
#if 0
        printError("Warning: tensor logarithm calculation iteration not converged\n");
#endif

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                tensor_out[i][j] = nan("");
    }

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            tensor_out[i][j] *= 2.0;
}
#endif

void convertSymmetric4thOrderMatrixToTensor(double tensor[3][3][3][3],
                                            double matrix[6][6])
{
    tensor[0][0][0][0] = matrix[0][0];
    tensor[0][0][0][1] = 0.5 * (matrix[0][3] + matrix[3][0]);
    tensor[0][0][0][2] = 0.5 * (matrix[0][5] + matrix[5][0]);
    tensor[0][0][1][0] = 0.5 * (matrix[0][3] + matrix[3][0]);
    tensor[0][0][1][1] = 0.5 * (matrix[0][1] + matrix[1][0]);
    tensor[0][0][1][2] = 0.5 * (matrix[0][4] + matrix[4][0]);
    tensor[0][0][2][0] = 0.5 * (matrix[0][5] + matrix[5][0]);
    tensor[0][0][2][1] = 0.5 * (matrix[0][4] + matrix[4][0]);
    tensor[0][0][2][2] = 0.5 * (matrix[0][2] + matrix[2][0]);
    tensor[0][1][0][0] = 0.5 * (matrix[3][0] + matrix[0][3]);
    tensor[0][1][0][1] = matrix[3][3];
    tensor[0][1][0][2] = 0.5 * (matrix[3][5] + matrix[5][3]);
    tensor[0][1][1][0] = 0.5 * (matrix[3][3] + matrix[3][3]);
    tensor[0][1][1][1] = 0.5 * (matrix[3][1] + matrix[1][3]);
    tensor[0][1][1][2] = 0.5 * (matrix[3][4] + matrix[4][3]);
    tensor[0][1][2][0] = 0.5 * (matrix[3][5] + matrix[5][3]);
    tensor[0][1][2][1] = 0.5 * (matrix[3][4] + matrix[4][3]);
    tensor[0][1][2][2] = 0.5 * (matrix[3][2] + matrix[2][3]);
    tensor[0][2][0][0] = 0.5 * (matrix[5][0] + matrix[0][5]);
    tensor[0][2][0][1] = 0.5 * (matrix[5][3] + matrix[3][5]);
    tensor[0][2][0][2] = matrix[5][5];
    tensor[0][2][1][0] = 0.5 * (matrix[5][3] + matrix[3][5]);
    tensor[0][2][1][1] = 0.5 * (matrix[5][1] + matrix[1][5]);
    tensor[0][2][1][2] = 0.5 * (matrix[5][4] + matrix[4][5]);
    tensor[0][2][2][0] = 0.5 * (matrix[5][5] + matrix[5][5]);
    tensor[0][2][2][1] = 0.5 * (matrix[5][4] + matrix[4][5]);
    tensor[0][2][2][2] = 0.5 * (matrix[5][2] + matrix[2][5]);
    tensor[1][0][0][0] = 0.5 * (matrix[3][0] + matrix[0][3]);
    tensor[1][0][0][1] = 0.5 * (matrix[3][3] + matrix[3][3]);
    tensor[1][0][0][2] = 0.5 * (matrix[3][5] + matrix[5][3]);
    tensor[1][0][1][0] = matrix[3][3];
    tensor[1][0][1][1] = 0.5 * (matrix[3][1] + matrix[1][3]);
    tensor[1][0][1][2] = 0.5 * (matrix[3][4] + matrix[4][3]);
    tensor[1][0][2][0] = 0.5 * (matrix[3][5] + matrix[5][3]);
    tensor[1][0][2][1] = 0.5 * (matrix[3][4] + matrix[4][3]);
    tensor[1][0][2][2] = 0.5 * (matrix[3][2] + matrix[2][3]);
    tensor[1][1][0][0] = 0.5 * (matrix[0][1] + matrix[1][0]);
    tensor[1][1][0][1] = 0.5 * (matrix[3][1] + matrix[1][3]);
    tensor[1][1][0][2] = 0.5 * (matrix[5][1] + matrix[1][5]);
    tensor[1][1][1][0] = 0.5 * (matrix[3][1] + matrix[1][3]);
    tensor[1][1][1][1] = matrix[1][1];
    tensor[1][1][1][2] = 0.5 * (matrix[1][4] + matrix[4][1]);
    tensor[1][1][2][0] = 0.5 * (matrix[5][1] + matrix[1][5]);
    tensor[1][1][2][1] = 0.5 * (matrix[1][4] + matrix[4][1]);
    tensor[1][1][2][2] = 0.5 * (matrix[1][2] + matrix[2][1]);
    tensor[1][2][0][0] = 0.5 * (matrix[0][4] + matrix[4][0]);
    tensor[1][2][0][1] = 0.5 * (matrix[3][4] + matrix[4][3]);
    tensor[1][2][0][2] = 0.5 * (matrix[5][4] + matrix[4][5]);
    tensor[1][2][1][0] = 0.5 * (matrix[3][4] + matrix[4][3]);
    tensor[1][2][1][1] = 0.5 * (matrix[4][1] + matrix[1][4]);
    tensor[1][2][1][2] = matrix[4][4];
    tensor[1][2][2][0] = 0.5 * (matrix[5][4] + matrix[4][5]);
    tensor[1][2][2][1] = 0.5 * (matrix[4][4] + matrix[4][4]);
    tensor[1][2][2][2] = 0.5 * (matrix[4][2] + matrix[2][4]);
    tensor[2][0][0][0] = 0.5 * (matrix[5][0] + matrix[0][5]);
    tensor[2][0][0][1] = 0.5 * (matrix[5][3] + matrix[3][5]);
    tensor[2][0][0][2] = 0.5 * (matrix[5][5] + matrix[5][5]);
    tensor[2][0][1][0] = 0.5 * (matrix[5][3] + matrix[3][5]);
    tensor[2][0][1][1] = 0.5 * (matrix[5][1] + matrix[1][5]);
    tensor[2][0][1][2] = 0.5 * (matrix[5][4] + matrix[4][5]);
    tensor[2][0][2][0] = matrix[5][5];
    tensor[2][0][2][1] = 0.5 * (matrix[5][4] + matrix[4][5]);
    tensor[2][0][2][2] = 0.5 * (matrix[5][2] + matrix[2][5]);
    tensor[2][1][0][0] = 0.5 * (matrix[0][4] + matrix[4][0]);
    tensor[2][1][0][1] = 0.5 * (matrix[3][4] + matrix[4][3]);
    tensor[2][1][0][2] = 0.5 * (matrix[5][4] + matrix[4][5]);
    tensor[2][1][1][0] = 0.5 * (matrix[3][4] + matrix[4][3]);
    tensor[2][1][1][1] = 0.5 * (matrix[4][1] + matrix[1][4]);
    tensor[2][1][1][2] = 0.5 * (matrix[4][4] + matrix[4][4]);
    tensor[2][1][2][0] = 0.5 * (matrix[5][4] + matrix[4][5]);
    tensor[2][1][2][1] = matrix[4][4];
    tensor[2][1][2][2] = 0.5 * (matrix[4][2] + matrix[2][4]);
    tensor[2][2][0][0] = 0.5 * (matrix[0][2] + matrix[2][0]);
    tensor[2][2][0][1] = 0.5 * (matrix[3][2] + matrix[2][3]);
    tensor[2][2][0][2] = 0.5 * (matrix[5][2] + matrix[2][5]);
    tensor[2][2][1][0] = 0.5 * (matrix[3][2] + matrix[2][3]);
    tensor[2][2][1][1] = 0.5 * (matrix[1][2] + matrix[2][1]);
    tensor[2][2][1][2] = 0.5 * (matrix[4][2] + matrix[2][4]);
    tensor[2][2][2][0] = 0.5 * (matrix[5][2] + matrix[2][5]);
    tensor[2][2][2][1] = 0.5 * (matrix[4][2] + matrix[2][4]);
    tensor[2][2][2][2] = matrix[2][2];
}

void convertSymmetric4thOrderTensorToMatrix(double matrix[6][6],
                                            double tensor[3][3][3][3])
{
    matrix[0][0] = tensor[0][0][0][0];
    matrix[0][1] = 0.5 * (tensor[0][0][1][1] + tensor[1][1][0][0]);
    matrix[0][2] = 0.5 * (tensor[0][0][2][2] + tensor[2][2][0][0]);
    matrix[0][3] = 0.25 * (tensor[0][0][0][1] + tensor[0][0][1][0] + tensor[0][1][0][0] + tensor[1][0][0][0]);
    matrix[0][4] = 0.25 * (tensor[0][0][1][2] + tensor[0][0][2][1] + tensor[1][2][0][0] + tensor[2][1][0][0]);
    matrix[0][5] = 0.25 * (tensor[0][0][2][0] + tensor[0][0][0][2] + tensor[2][0][0][0] + tensor[0][2][0][0]);
    matrix[1][0] = 0.5 * (tensor[1][1][0][0] + tensor[0][0][1][1]);
    matrix[1][1] = tensor[1][1][1][1];
    matrix[1][2] = 0.5 * (tensor[1][1][2][2] + tensor[2][2][1][1]);
    matrix[1][3] = 0.25 * (tensor[1][1][0][1] + tensor[1][1][1][0] + tensor[0][1][1][1] + tensor[1][0][1][1]);
    matrix[1][4] = 0.25 * (tensor[1][1][1][2] + tensor[1][1][2][1] + tensor[1][2][1][1] + tensor[2][1][1][1]);
    matrix[1][5] = 0.25 * (tensor[1][1][2][0] + tensor[1][1][0][2] + tensor[2][0][1][1] + tensor[0][2][1][1]);
    matrix[2][0] = 0.5 * (tensor[2][2][0][0] + tensor[0][0][2][2]);
    matrix[2][1] = 0.5 * (tensor[2][2][1][1] + tensor[1][1][2][2]);
    matrix[2][2] = tensor[2][2][2][2];
    matrix[2][3] = 0.25 * (tensor[2][2][0][1] + tensor[2][2][1][0] + tensor[0][1][2][2] + tensor[1][0][2][2]);
    matrix[2][4] = 0.25 * (tensor[2][2][1][2] + tensor[2][2][2][1] + tensor[1][2][2][2] + tensor[2][1][2][2]);
    matrix[2][5] = 0.25 * (tensor[2][2][2][0] + tensor[2][2][0][2] + tensor[2][0][2][2] + tensor[0][2][2][2]);
    matrix[3][0] = 0.25 * (tensor[0][1][0][0] + tensor[1][0][0][0] + tensor[0][0][0][1] + tensor[0][0][1][0]);
    matrix[3][1] = 0.25 * (tensor[0][1][1][1] + tensor[1][0][1][1] + tensor[1][1][0][1] + tensor[1][1][1][0]);
    matrix[3][2] = 0.25 * (tensor[0][1][2][2] + tensor[1][0][2][2] + tensor[2][2][0][1] + tensor[2][2][1][0]);
    matrix[3][3] = 0.25 * (tensor[0][1][0][1] + tensor[0][1][1][0] + tensor[1][0][0][1] + tensor[1][0][1][0]);
    matrix[3][4] = 0.125 * (tensor[0][1][1][2] + tensor[0][1][2][1] + tensor[1][0][1][2] + tensor[1][0][2][1] + tensor[1][2][0][1] + tensor[1][2][1][0] + tensor[2][1][0][1] + tensor[2][1][1][0]);
    matrix[3][5] = 0.125 * (tensor[0][1][2][0] + tensor[0][1][0][2] + tensor[1][0][2][0] + tensor[1][0][0][2] + tensor[2][0][0][1] + tensor[2][0][1][0] + tensor[0][2][0][1] + tensor[0][2][1][0]);
    matrix[4][0] = 0.25 * (tensor[1][2][0][0] + tensor[2][1][0][0] + tensor[0][0][1][2] + tensor[0][0][2][1]);
    matrix[4][1] = 0.25 * (tensor[1][2][1][1] + tensor[2][1][1][1] + tensor[1][1][1][2] + tensor[1][1][2][1]);
    matrix[4][2] = 0.25 * (tensor[1][2][2][2] + tensor[2][1][2][2] + tensor[2][2][1][2] + tensor[2][2][2][1]);
    matrix[4][3] = 0.125 * (tensor[1][2][0][1] + tensor[1][2][1][0] + tensor[2][1][0][1] + tensor[2][1][1][0] + tensor[0][1][1][2] + tensor[0][1][2][1] + tensor[1][0][1][2] + tensor[1][0][2][1]);
    matrix[4][4] = 0.25 * (tensor[1][2][1][2] + tensor[1][2][2][1] + tensor[2][1][1][2] + tensor[2][1][2][1]);
    matrix[4][5] = 0.125 * (tensor[1][2][2][0] + tensor[1][2][0][2] + tensor[2][1][2][0] + tensor[2][1][0][2] + tensor[2][0][1][2] + tensor[2][0][2][1] + tensor[0][2][1][2] + tensor[0][2][2][1]);
    matrix[5][0] = 0.25 * (tensor[2][0][0][0] + tensor[0][2][0][0] + tensor[0][0][2][0] + tensor[0][0][0][2]);
    matrix[5][1] = 0.25 * (tensor[2][0][1][1] + tensor[0][2][1][1] + tensor[1][1][2][0] + tensor[1][1][0][2]);
    matrix[5][2] = 0.25 * (tensor[2][0][2][2] + tensor[0][2][2][2] + tensor[2][2][2][0] + tensor[2][2][0][2]);
    matrix[5][3] = 0.125 * (tensor[2][0][0][1] + tensor[2][0][1][0] + tensor[0][2][0][1] + tensor[0][2][1][0] + tensor[0][1][2][0] + tensor[0][1][0][2] + tensor[1][0][2][0] + tensor[1][0][0][2]);
    matrix[5][4] = 0.125 * (tensor[2][0][1][2] + tensor[2][0][2][1] + tensor[0][2][1][2] + tensor[0][2][2][1] + tensor[1][2][2][0] + tensor[1][2][0][2] + tensor[2][1][2][0] + tensor[2][1][0][2]);
    matrix[5][5] = 0.25 * (tensor[2][0][2][0] + tensor[2][0][0][2] + tensor[0][2][2][0] + tensor[0][2][0][2]);
}

void conver4thOrderTensorToMatrix(double matrix[9][9], double tensor[3][3][3][3]){
    matrix[0][0] = tensor[0][0][0][0];
    matrix[0][1] = tensor[0][0][0][1];
    matrix[0][2] = tensor[0][0][0][2];
    matrix[0][3] = tensor[0][0][1][0];
    matrix[0][4] = tensor[0][0][1][1];
    matrix[0][5] = tensor[0][0][1][2];
    matrix[0][6] = tensor[0][0][2][0];
    matrix[0][7] = tensor[0][0][2][1];
    matrix[0][8] = tensor[0][0][2][2];
    matrix[1][0] = tensor[0][1][0][0];
    matrix[1][1] = tensor[0][1][0][1];
    matrix[1][2] = tensor[0][1][0][2];
    matrix[1][3] = tensor[0][1][1][0];
    matrix[1][4] = tensor[0][1][1][1];
    matrix[1][5] = tensor[0][1][1][2];
    matrix[1][6] = tensor[0][1][2][0];
    matrix[1][7] = tensor[0][1][2][1];
    matrix[1][8] = tensor[0][1][2][2];
    matrix[2][0] = tensor[0][2][0][0];
    matrix[2][1] = tensor[0][2][0][1];
    matrix[2][2] = tensor[0][2][0][2];
    matrix[2][3] = tensor[0][2][1][0];
    matrix[2][4] = tensor[0][2][1][1];
    matrix[2][5] = tensor[0][2][1][2];
    matrix[2][6] = tensor[0][2][2][0];
    matrix[2][7] = tensor[0][2][2][1];
    matrix[2][8] = tensor[0][2][2][2];
    matrix[3][0] = tensor[1][0][0][0];
    matrix[3][1] = tensor[1][0][0][1];
    matrix[3][2] = tensor[1][0][0][2];
    matrix[3][3] = tensor[1][0][1][0];
    matrix[3][4] = tensor[1][0][1][1];
    matrix[3][5] = tensor[1][0][1][2];
    matrix[3][6] = tensor[1][0][2][0];
    matrix[3][7] = tensor[1][0][2][1];
    matrix[3][8] = tensor[1][0][2][2];
    matrix[4][0] = tensor[1][1][0][0];
    matrix[4][1] = tensor[1][1][0][1];
    matrix[4][2] = tensor[1][1][0][2];
    matrix[4][3] = tensor[1][1][1][0];
    matrix[4][4] = tensor[1][1][1][1];
    matrix[4][5] = tensor[1][1][1][2];
    matrix[4][6] = tensor[1][1][2][0];
    matrix[4][7] = tensor[1][1][2][1];
    matrix[4][8] = tensor[1][1][2][2];
    matrix[5][0] = tensor[1][2][0][0];
    matrix[5][1] = tensor[1][2][0][1];
    matrix[5][2] = tensor[1][2][0][2];
    matrix[5][3] = tensor[1][2][1][0];
    matrix[5][4] = tensor[1][2][1][1];
    matrix[5][5] = tensor[1][2][1][2];
    matrix[5][6] = tensor[1][2][2][0];
    matrix[5][7] = tensor[1][2][2][1];
    matrix[5][8] = tensor[1][2][2][2];
    matrix[6][0] = tensor[2][0][0][0];
    matrix[6][1] = tensor[2][0][0][1];
    matrix[6][2] = tensor[2][0][0][2];
    matrix[6][3] = tensor[2][0][1][0];
    matrix[6][4] = tensor[2][0][1][1];
    matrix[6][5] = tensor[2][0][1][2];
    matrix[6][6] = tensor[2][0][2][0];
    matrix[6][7] = tensor[2][0][2][1];
    matrix[6][8] = tensor[2][0][2][2];
    matrix[7][0] = tensor[2][1][0][0];
    matrix[7][1] = tensor[2][1][0][1];
    matrix[7][2] = tensor[2][1][0][2];
    matrix[7][3] = tensor[2][1][1][0];
    matrix[7][4] = tensor[2][1][1][1];
    matrix[7][5] = tensor[2][1][1][2];
    matrix[7][6] = tensor[2][1][2][0];
    matrix[7][7] = tensor[2][1][2][1];
    matrix[7][8] = tensor[2][1][2][2];
    matrix[8][0] = tensor[2][2][0][0];
    matrix[8][1] = tensor[2][2][0][1];
    matrix[8][2] = tensor[2][2][0][2];
    matrix[8][3] = tensor[2][2][1][0];
    matrix[8][4] = tensor[2][2][1][1];
    matrix[8][5] = tensor[2][2][1][2];
    matrix[8][6] = tensor[2][2][2][0];
    matrix[8][7] = tensor[2][2][2][1];
    matrix[8][8] = tensor[2][2][2][2];
}

void calculateTensorLogarithmDerivative(double tensor_out[3][3][3][3],
                                        double tensor_in[3][3])
{
    calculateIsotropicTensorFunctionDerivative(tensor_out, tensor_in,
                                               log, calculateInverse);
}

void calculateIsotropicTensorFunctionDerivative(double tensor_out[3][3][3][3],
                                                double tensor_in[3][3],
                                                double (*function)(const double variable),
                                                double (*function_derivative)(const double variable))
{
    double eigenvalue_tolerance = 1.0e-3;

    double new_eigenvalues[3];
    double new_eigenvalue_derivatives[3];
    double eigenvalues[3];
    double eigenprojections[3][3][3];
    double tolerance;
    int i, j, k, l, m;

    /* Calculate eigenvalues and eigenprojections */
    calculateEigenvalues(eigenvalues, tensor_in);
    calculateEigenprojections(eigenprojections, eigenvalues, tensor_in);

    /* Calculate new eigenvalues */
    for (i = 0; i < 3; i++)
        new_eigenvalues[i]
            = function(eigenvalues[i]);
    for (i = 0; i < 3; i++)
        new_eigenvalue_derivatives[i]
            = function_derivative(eigenvalues[i]);

    /* Calculate eigenvalue tolerance */
    tolerance = eigenvalue_tolerance * (fabs(eigenvalues[0])
                                        + fabs(eigenvalues[1])
                                        + fabs(eigenvalues[2])) / 3.0;

    /* Calculate new tensor when x_1 == x_2 == x_3 */
    if (fabs(eigenvalues[0] - eigenvalues[1]) <= tolerance
        && fabs(eigenvalues[1] - eigenvalues[2]) <= tolerance
        && fabs(eigenvalues[2] - eigenvalues[0]) <= tolerance)
    {
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                    for (l = 0; l < 3; l++)
                        tensor_out[i][j][k][l]
                            = new_eigenvalue_derivatives[0]
                            * 0.5 * (identity_tensor[i][k] * identity_tensor[j][l] + identity_tensor[i][l] * identity_tensor[j][k]);

        return;
    }

    /* Calculate new tensor when x_a != x_b == x_c */
    for (m = 0; m < 3; m++)
        if (fabs(eigenvalues[(m + 1) % 3] - eigenvalues[(m + 2) % 3]) <= tolerance)
        {
            const double s1
                = (new_eigenvalues[m] - new_eigenvalues[(m + 2) % 3])
                / calculateSquare(eigenvalues[m] - eigenvalues[(m + 2) % 3])
                - new_eigenvalue_derivatives[(m + 2) % 3]
                / (eigenvalues[m] - eigenvalues[(m + 2) % 3]);
            const double s2
                = 2.0 * eigenvalues[(m + 2) % 3]
                * (new_eigenvalues[m] - new_eigenvalues[(m + 2) % 3])
                / calculateSquare(eigenvalues[m] - eigenvalues[(m + 2) % 3])
                - (eigenvalues[m] + eigenvalues[(m + 2) % 3])
                / (eigenvalues[m] - eigenvalues[(m + 2) % 3])
                * new_eigenvalue_derivatives[(m + 2) % 3];
            const double s3
                = 2.0
                * (new_eigenvalues[m] - new_eigenvalues[(m + 2) % 3])
                / calculateCube(eigenvalues[m] - eigenvalues[(m + 2) % 3])
                - (new_eigenvalue_derivatives[m] + new_eigenvalue_derivatives[(m + 2) % 3])
                / calculateSquare(eigenvalues[m] - eigenvalues[(m + 2) % 3]);
            const double s4
                = eigenvalues[(m + 2) % 3] * s3;
            const double s5
                = s4;
            const double s6
                = calculateSquare(eigenvalues[(m + 2) % 3]) * s3;

            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    for (k = 0; k < 3; k++)
                        for (l = 0; l < 3; l++)
                            tensor_out[i][j][k][l]
                                = s1 * 0.5 * (identity_tensor[i][k] * tensor_in[l][j]
                                              + identity_tensor[i][l] * tensor_in[k][j]
                                              + identity_tensor[j][l] * tensor_in[i][k]
                                              + identity_tensor[k][j] * tensor_in[i][l])
                                - s2 * 0.5 * (identity_tensor[i][k] * identity_tensor[j][l]
                                              + identity_tensor[i][l] * identity_tensor[j][k])
                                - s3 * tensor_in[i][j] * tensor_in[k][l]
                                + s4 * tensor_in[i][j] * identity_tensor[k][l]
                                + s5 * identity_tensor[i][j] * tensor_in[k][l]
                                - s6 * identity_tensor[i][j] * identity_tensor[k][l];

            return;
        }

    /* Calculate new tensor when x_1 != x_2 != x_3 */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                for (l = 0; l < 3; l++)
                {
                    double temp = 0.0;

                    for (m = 0; m < 3; m++)
                    {
                        const double s1
                            = new_eigenvalues[m]
                            / ((eigenvalues[m] - eigenvalues[(m + 1) % 3])
                               * (eigenvalues[m] - eigenvalues[(m + 2) % 3]));
                        const double s2
                            = s1
                            * (eigenvalues[(m + 1) % 3] + eigenvalues[(m + 2) % 3]);
                        const double s3
                            = s1
                            * ((eigenvalues[m] - eigenvalues[(m + 1) % 3])
                               + (eigenvalues[m] - eigenvalues[(m + 2) % 3]));
                        const double s4
                            = s1
                            * (eigenvalues[(m + 1) % 3] - eigenvalues[(m + 2) % 3]);

                        temp
                            += s1 * 0.5 * (identity_tensor[i][k] * tensor_in[l][j]
                                           + identity_tensor[i][l] * tensor_in[k][j]
                                           + identity_tensor[j][l] * tensor_in[i][k]
                                           + identity_tensor[k][j] * tensor_in[i][l])
                            -  s2 * 0.5 * (identity_tensor[i][k] * identity_tensor[j][l]
                                           + identity_tensor[i][l] * identity_tensor[j][k])
                            -  s3 * eigenprojections[m][i][j] * eigenprojections[m][k][l]
                            -  s4 * (eigenprojections[(m + 1) % 3][i][j] * eigenprojections[(m + 1) % 3][k][l]
                                     - eigenprojections[(m + 2) % 3][i][j] * eigenprojections[(m + 2) % 3][k][l])
                            +  new_eigenvalue_derivatives[m] * eigenprojections[m][i][j] * eigenprojections[m][k][l];
                    }

                    tensor_out[i][j][k][l]
                        = temp;
                }
}

void calculateEigenvalues(double eigenvalues[3],
                          double tensor[3][3])
#if 1
{
    const double discriminant_tolerance = 1.0e-3;

    double square_tensor[3][3];
    double i1, i2, i3;
    double q, sqrt_q3, r;

    /* Calculate square of tensor */
    calculate3x3MatrixSquare(square_tensor, tensor);

    /* Calculate invariants */
    i1 = calculate3x3MatrixTrace(tensor);
    i2 = 0.5 * (i1 * i1 - calculate3x3MatrixTrace(square_tensor));
    i3 = calc_3x3matrix_determinant(tensor);

    /* Calculate Q and R */
    q  = (i1 * i1 - 3.0 * i2) / 9.0;
    sqrt_q3 = sqrt(q * q * q);
    r  = (-2.0 * i1 * i1 * i1 + 9.0 * i1 * i2 - 27.0 * i3) / 54.0;

    /* Calculate eigenvalues */
    if (q < 0.0
        || fabs(q) <= discriminant_tolerance * (i1 * i1 / 9.0))
    {
        const double b = i1 / 3.0;

        eigenvalues[0] = b;
        eigenvalues[1] = b;
        eigenvalues[2] = b;
    }
    else if (r - sqrt_q3 > 0
             || fabs(r - sqrt_q3) <= discriminant_tolerance * fabs(sqrt_q3))
    {
        const double a = sqrt(q);
        const double b = i1 / 3.0;

        eigenvalues[0] = -2.0 * a + b;
        eigenvalues[1] = a + b;
        eigenvalues[2] = a + b;
    }
    else if (r + sqrt_q3 < 0
             || fabs(r + sqrt_q3) <= discriminant_tolerance * fabs(sqrt_q3))
    {
        const double a = sqrt(q);
        const double b = i1 / 3.0;

        eigenvalues[0] = -a + b;
        eigenvalues[1] = 2.0 * a + b;
        eigenvalues[2] = -a + b;
    }
    else
    {
        const double theta = acos(r / sqrt_q3) / 3.0;
        const double cos_theta = cos(theta);
        const double sin_theta = sin(theta);
        const double a = sqrt(q);
        const double b = i1 / 3.0;

        eigenvalues[0] = -2.0 * a * cos_theta + b;
        eigenvalues[1] = a * (cos_theta + sqrt(3.0) * sin_theta) + b;
        eigenvalues[2] = a * (cos_theta - sqrt(3.0) * sin_theta) + b;
    }

    /* Sort and reverse eigenvalues */
    sortVector(eigenvalues, 3);
    reverseVector(eigenvalues, 3);
}
#endif

void calculateEigenprojections(double eigenprojections[3][3][3],
                               const double eigenvalues[3],
                               double tensor[3][3])
{
    double square_tensor[3][3];
    int i, j, k;

    /* Calculate square of tensor */
    calculate3x3MatrixSquare(square_tensor, tensor);

    /* Calculate eigenprojections */
    for (i = 0; i < 3; i++)
    {
        const double i1
            = calculate3x3MatrixTrace(tensor);
        const double i3
            = calc_3x3matrix_determinant(tensor);
        const double coefficient
            = eigenvalues[i]
            / (2.0 * calculateCube(eigenvalues[i])
               - i1 * calculateSquare(eigenvalues[i])
               + i3);

        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                eigenprojections[i][j][k]
                    = coefficient
                    * (square_tensor[j][k]
                       - (i1 - eigenvalues[i]) * tensor[j][k]
                       + i3 / eigenvalues[i] * identity_tensor[j][k]);
    }
}