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
    int solver_type = 2;
    double *du;

    Init_model();

    switch(solver_type){
        case NON_LINEAR_SOLVER:
            analize_by_NewtonRapdon();
        break;
        case LINEAR_SOLVER:
            if((du = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
                printf("Error: du's Memory is not enough\n");
                exit(-1);
            }
            Linear_analization(du);
            Output_Linear_data(du);
            free(du);
        break;
    }
    
    printf("complete!\n");
    break_model_memory();
    return 0;
}