#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"matrix.h"

extern Option option;


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

    calc_3x3_matrix_square(square_tensor, tensor);

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