#include<stdio.h>
#include"type.h"
extern Global global;
extern Option option;


void update_field_and_internal_forces(){

    double subdomain_internal_force[20][3];
    double d_matrix[6][6];          //Dマトリクス
    double inverse_relative_deformation_gradient[3][3];         //相対変形勾配テンソルの逆テンソル
    double relative_deformation_gradient[3][3];                 //相対変形勾配テンソル
    double elastic_strain_tensor[3][3];                         //弾性ひずみテンソル
    double elastic_left_cauchy_green_deformations[3][3];        //弾性左コーシーグリーンテンソル
    double trial_elastic_left_cauchy_green_deformations[3][3];  //試行弾性左コーシーグリーンテンソル
    double trial_relative_stresses[6];                          //試行相対応力
    double X[27][3];
    double w[27];                                               //ガウス求積に使う正規化座標と重み関数

    //internal_forceをゼロ処理
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            global.subdomain.internal_force[i][j] = 0.;
        }
    }

    for(int point = 0; point < global.subdomain.N_point; point++){
        double kinematic_hardening_fractions = global.material.kinematic_hardening_fractions;
        double N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];
        //サブドメインの内力ベクトルをゼロ処理
        for(int i = 0; i < N_support; i++){
            for(int j = 0; j < option.dim; j++){
                subdomain_internal_force[i][j] = 0.;
            }
        }

    }




}

void zero_fill_displacement_increments(){
    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i  = 0; i < option.dim; i++){
            global.subdomain.displacement[point][i] = 0.;
        }
    }
}