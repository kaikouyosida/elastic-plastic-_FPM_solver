#include<stdio.h>
#include"s_matrix.h"
#include"type.h"

extern Option option;

void generateSMatrix(double s_matrix[9][9],
                     double current_stresses[6])
{
    //Sマトリクスのゼロ処理
    #if 0
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            s_matrix[i][j] = 0.0;
    
    for(int i = 0; i < 6; i++){
        printf("%+15.14e    ", current_stresses[i]);
    }
    printf("\n");
    #endif
    //Sマトリクスの計算
    for (int i = 0; i < option.dim; i++)
    {
        s_matrix[i][i]     = current_stresses[0];
        s_matrix[i][i + 3] = current_stresses[3];
        s_matrix[i][i + 6] = current_stresses[5];

        s_matrix[i + 3][i]     = current_stresses[3];
        s_matrix[i + 3][i + 3] = current_stresses[1];
        s_matrix[i + 3][i + 6] = current_stresses[4];

        s_matrix[i + 6][i]     = current_stresses[5];
        s_matrix[i + 6][i + 3] = current_stresses[4];
        s_matrix[i + 6][i + 6] = current_stresses[2];
    }
}