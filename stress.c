#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

#include<math.h>
#include"stress.h"
#include"scalar.h"

double calc_equivalent_stress(const double *stresses)
{
    /* Calculate vonMises' equivalent stress */
    return sqrt(0.5 * (calculateSquare(stresses[0] - stresses[1])
                       + calculateSquare(stresses[1] - stresses[2])
                       + calculateSquare(stresses[2] - stresses[0])
                       + 6.0 * (calculateSquare(stresses[3])
                                + calculateSquare(stresses[4])
                                + calculateSquare(stresses[5]))));
}