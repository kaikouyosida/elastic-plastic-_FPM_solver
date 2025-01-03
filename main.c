#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

#include<stdio.h>
#include<stdlib.h>
#include"model.h"
#include"type.h"
#include"fpm.h"
#include"Output_data.h"
#include"matrix.h"
#include"stress.h"

Global global;
Option option;
SS_CURVE ss_curve;

#define NON_LINEAR_SOLVER 1
#define INFINITESIMAL 2

int main(){
    
    Init_model();
    option.solver_type = 1;
    
    switch(option.solver_type){
        case NON_LINEAR_SOLVER:
            analize_by_NewtonRaphson();
        break;
        case INFINITESIMAL:
            infinitesimal_analization();
        break;
    }


    
    //nit_model()で確保したメモリの開放
    break_model_memory();

    printf("complete!\n");
    return 0;
}