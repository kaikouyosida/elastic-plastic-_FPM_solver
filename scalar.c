#include "scalar.h"
double calculateInverse(const double value)
{
    return 1.0 / value;
}

/*
 * Calculate square
 */
double calculateSquare(const double value)
{
    return value * value;
}

/*
 * Calculate cube
 */
double calculateCube(const double value)
{
    return value * value * value;
}

/*
 * Swap reals
 */
void swapReals(double *value1, double *value2)
{
    double temp = *value1;
    *value1  = *value2;
    *value2  = temp;
}

/*
 * Swap integers
 */
void swapIntegers(int *value1, int *value2)
{
    int temp = *value1;
    *value1  = *value2;
    *value2  = temp;
}
