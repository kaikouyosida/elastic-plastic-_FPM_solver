#include<stdio.h>
#include<stdlib.h>
#include"model.h"
#include"type.h"
#include"fpm.h"
Global global;
Option option;
SS_CURVE ss_curve;

int main(){
    Init_model();
    analize_by_NewtonRapdon();
    
    printf("complete!\n");
    break_model_memory();
    return 0;
}