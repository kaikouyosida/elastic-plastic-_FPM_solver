#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"field.h"
#include"matrix.h"
#include"ss_curve.h"
extern Global global;
extern SS_CURVE ss_curve;
extern Option option;

void init_field(){
    //変位ベクトル、変位増分ベクトルをゼロ処理
    global.subdomain.displacement = matrix(global.subdomain.N_point, option.dim);
    global.subdomain.displacement_increment = matrix(global.subdomain.N_point, option.dim);

    //全体剛性マトリクスをゼロ処理
    if((global.subdomain.Global_K = (double *)calloc(option.dim * global.subdomain.N_point * option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:Global_K's memory is not enough\n");
        exit(-1);
    }

    //残差ベクトル{r} = {Fint} - λ{Fext}の計算
    if((global.subdomain.global_residual_force = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:Global_residual_force's memory is not enough\n");
        exit(-1);
    }

    //外力ベクトルをゼロ処理
    global.subdomain.external_force = matrix(global.subdomain.N_point, option.dim);
    global.subdomain.global_external_force = matrix(global.subdomain.N_point, option.dim);
    global.subdomain.previous_global_external_force = matrix(global.subdomain.N_point, option.dim);
    
    //内力ベクトルをゼロ処理
    global.subdomain.internal_force = matrix(global.subdomain.N_point, option.dim);
    global.subdomain.global_internal_force = matrix(global.subdomain.N_point, option.dim);

    //変形勾配テンソルをゼロ処理
    global.subdomain.deformation_gradients = threetimes_tensor(3, 3, global.subdomain.N_point);
    global.subdomain.current_deformation_gradients = threetimes_tensor(3, 3, global.subdomain.N_point);
    //変形勾配テンソルを単位テンソルにする
    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                if(i == j){
                    global.subdomain.deformation_gradients[i][j][point] = 1.0;
                    global.subdomain.current_deformation_gradients[i][j][point] = 1.0;
                }
            }
        }
    }

    //弾性ひずみ、現配置の弾性ひずみ、試行弾性ひずみをゼロ処理
    global.subdomain.elastic_strains = matrix(global.subdomain.N_point, 6);
    global.subdomain.current_elastic_strains = matrix(global.subdomain.N_point, 6);
    global.subdomain.trial_elastic_strains = matrix(global.subdomain.N_point, 6);

    //応力、現配置の応力をゼロ処理
    global.subdomain.stresses = matrix(global.subdomain.N_point, 6);
    global.subdomain.current_stresses = matrix(global.subdomain.N_point, 6);

    //相当応力、相当塑性ひずみ、相当塑性ひずみ増分をゼロ処理
    if((global.subdomain.equivalent_stresses = (double *)calloc(global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error: subdomain_eq_stresses's memory is not enough\n");
        exit(-1);
    }
    if((global.subdomain.equivalent_plastic_strains = (double *)calloc(global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error: subdomain_eq_strains's memory is not enough\n");
        exit(-1);
    }
    if((global.subdomain.equivalent_plastic_strain_increments = (double *)calloc(global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error: subdomain_eq_strains_increment's memory is not enough\n");
        exit(-1);
    }

    //降伏応力をゼロ処理
    if((global.subdomain.yield_stresses = (double *)calloc(global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error: yield_stresses's memory is not enough\n");
        exit(-1);
    }
    //初期の降伏応力値をyield_stressesに格納
    for(int point = 0; point < global.subdomain.N_point; point++){
        global.subdomain.yield_stresses[point] = get_hardening_stress(0.0);
    }
    
    //現配置の降伏応力をゼロ処理
    if((global.subdomain.current_yield_stresses = (double *)calloc(global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error: current_yield_stresses's memory is not enough\n");
        exit(-1);
    }
    for(int point = 0; point < global.subdomain.N_point; point++){
        global.subdomain.current_yield_stresses[point] = global.subdomain.yield_stresses[point];
    }
    
    //背応力をゼロ処理
    global.subdomain.back_stresses = matrix(global.subdomain.N_point, 6);
    global.subdomain.current_back_stresses = matrix(global.subdomain.N_point, 6);

    //節点値をゼロ処理
    global.subdomain.nodal_displacements = matrix(global.subdomain.N_node, option.dim);
    global.subdomain.nodal_displacement_increments = matrix(global.subdomain.N_node, option.dim);
    global.subdomain.nodal_stresses = matrix(global.subdomain.N_node, option.dim);

    if((global.subdomain.nodal_equivalent_plastic_strains = (double *)calloc(global.subdomain.N_node, sizeof(double))) == NULL){
        printf("Error: nodal_eq_plastic_strains's memory is not enough\n");
        exit(-1);
    }
    if((global.subdomain.nodal_yield_stresses = (double *)calloc(global.subdomain.N_node, sizeof(double))) == NULL){
        printf("Error: nodal_yield_stresses's memory is not enough\n");
        exit(-1);
    }
    global.subdomain.nodal_back_stresses = matrix(global.subdomain.N_node, 6);

}
void break_field(){
    printf("check0\n");
    free_matrix(global.subdomain.nodal_back_stresses);
    free(global.subdomain.nodal_yield_stresses);
    free(global.subdomain.nodal_equivalent_plastic_strains);
    free_matrix(global.subdomain.nodal_stresses);
    free_matrix(global.subdomain.nodal_displacement_increments);
    free_matrix(global.subdomain.nodal_displacements);
    free_matrix(global.subdomain.current_back_stresses);
    free_matrix(global.subdomain.back_stresses);
    free(global.subdomain.current_yield_stresses);
    free(global.subdomain.yield_stresses);
    free(global.subdomain.equivalent_plastic_strain_increments);
    free(global.subdomain.equivalent_plastic_strains);
    free(global.subdomain.equivalent_stresses);
    free_matrix(global.subdomain.current_stresses);
    free_matrix(global.subdomain.stresses);
    free_matrix(global.subdomain.trial_elastic_strains);
    free_matrix(global.subdomain.current_elastic_strains);
    free_matrix(global.subdomain.elastic_strains);
    free_tensor(global.subdomain.current_deformation_gradients);
    free_tensor(global.subdomain.deformation_gradients);
    free_matrix(global.subdomain.global_internal_force);
    free_matrix(global.subdomain.internal_force);
    free_matrix(global.subdomain.previous_global_external_force);
    free_matrix(global.subdomain.global_external_force);
    free_matrix(global.subdomain.external_force);
    free(global.subdomain.global_residual_force);
    free(global.subdomain.Global_K);
    free_matrix(global.subdomain.displacement_increment);
    free_matrix(global.subdomain.displacement);
}