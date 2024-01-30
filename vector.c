#include<stdio.h>
#include<math.h>

#include"scalar.h"
#include"vector.h"

void sortVector(double *vector, const int num)
{
    int i, j;

    for (i = 0; i < num - 1; i++)
        for (j = num - 1; j > i; j--)
            if (vector[j - 1] > vector[j])
                swapReals(&vector[j - 1], &vector[j]);
}

void reverseVector(double *vector, const int num)
{
    int i;

    for (i = 0; i < num / 2; i++)
        swapReals(&vector[i], &vector[num - 1 - i]);
}