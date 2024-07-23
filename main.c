#include<stdio.h>
#include<stdlib.h>
#include"model.h"
#include"type.h"
#include"fpm.h"
#include"Output.h"

Global global;
Option option;
SS_CURVE ss_curve;

#define NON_LINEAR_SOLVER 1
#define LINEAR_SOLVER 2

int main(){
    int solver_type = 1;

    Init_model();

    switch(solver_type){
        case NON_LINEAR_SOLVER:
            analize_by_NewtonRaphson();
        break;
        case LINEAR_SOLVER:
            Linear_analization();
        break;
    }
    
    //nit_model()で確保したメモリの開放
    break_model_memory();

    printf("complete!\n");
    return 0;
}