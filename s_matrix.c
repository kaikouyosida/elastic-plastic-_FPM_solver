#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

#include<stdio.h>
#include"s_matrix.h"
#include"type.h"

extern Option option;

void generateSMatrix(double s_matrix[9][9],
                     const double current_stresses[6])
{
    //Sマトリクスのゼロ処理
    
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            s_matrix[i][j] = 0.0;
    
    //Sマトリクスの計算
    for (int i = 0; i < option.dim; i++)
    {
        s_matrix[option.dim  *  i][option.dim  *  i]     = current_stresses[0];
        s_matrix[option.dim  *  i][option.dim  *  i + 1] = current_stresses[3];
        s_matrix[option.dim  *  i][option.dim  *  i + 2] = current_stresses[5];

        s_matrix[option.dim  *  i + 1][option.dim  *  i]     = current_stresses[3];
        s_matrix[option.dim  *  i + 1][option.dim  *  i + 1]= current_stresses[1];
        s_matrix[option.dim  *  i + 1][option.dim  *  i + 2] = current_stresses[4];

        s_matrix[option.dim  *  i + 2][option.dim  *  i]     = current_stresses[5];
        s_matrix[option.dim  *  i + 2][option.dim  *  i + 1]= current_stresses[4];
        s_matrix[option.dim  *  i + 2][option.dim  *  i + 2] = current_stresses[2];
    }
    
}