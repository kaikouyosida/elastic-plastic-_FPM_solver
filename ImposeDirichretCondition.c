#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "type.h"
#include "ImposeDirichretCondition.h"


extern Option option;
extern Global global;

//強制変位を設定
double fixed_deformation(double time, double time_end, double x1, double x2, double x3, int type){
  double fixed_u;

  if(type == 0){
    // 変位を固定 //
		fixed_u = 0.0;
	}
  
  else if(type == 1){
    // 変位を固定 //
		fixed_u = 0.01 * time / time_end;
	}
  
  else if(type == 2){
    // x軸方向のTimoshenko梁の変位固定 //
		double E_mod = global.material.E_mod;
    double nu_mod = global.material.nu_mod;
    double I = 1.0 / 12.0;
    double P = 1.0;
    double L = 10.0;
    fixed_u = P / (6.0 * E_mod * I) * (3.0 * x1 * (2.0 * L - x1) + (2.0 + nu_mod) * (x2 * x2 - 1.0 / 4.0)) * x2;
    
     if(time == 0)
      fixed_u = 0;
	}
  else if(type == 2){
    // y軸方向のTimoshenko梁の変位固定 //
		double E_mod = global.material.E_mod;
    double nu_mod = global.material.nu_mod;
    double I = 1.0 / 12.0;
    double P = 1.0;
    double L = 10.0;
    fixed_u = -P / (6.0 * E_mod * I) * (x1 * x1 * (3.0 * L - x1) + 3.0 * nu_mod * (L - x1) * x2 * x2 + (4.0 + 5.0 * nu_mod)/ 4.0 * x1);

    if(time == 0)
      fixed_u = 0;
	}
  else if(type == 3){
    // 時間に比例した変形 //
		fixed_u = (0.1*x2)*(time/time_end);
	}
  else if(type == 4){
    // 時間に比例した変形 //
		fixed_u = (1.0*x1)*(time/time_end);
	}
  return fixed_u;
}

//残差ベクトルにディリクレ条件を反映
void ImposeDirichretResidual(const int iteration_step)
{
    int count = 0;                                         //カウンタ
    int DoF_free = option.dim * global.subdomain.N_point;  // 拘束を考慮しない自由度
    double fixed_u_inc = 0.; // 規定された変位の増分量
    int type_num[3];         // ディリクレ条件の番号
    double fixed_xyz[3];     // 変位を固定するポイントの位置
    
    if (iteration_step == 0)
    {
        for (int i = 0; i < global.subdomain.N_point; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                type_num[j] = 0;
                fixed_xyz[j] = 0.;
            }
            for (int j = 0; j < 3; j++)
                fixed_xyz[j] = global.subdomain.point_XYZ[3 * i + j];

            // x方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 2 == 0)
            {
                type_num[0] = global.bc.Dirichlet_type[i] % 100;
                fixed_u_inc = fixed_deformation(option.time, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[0]);
                fixed_u_inc -= fixed_deformation(option.time_old, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[0]);
                global.subdomain.global_residual_force[i][0] = fixed_u_inc;

                for(int j = 0; j < global.subdomain.N_point; j++){
                  for(int k = 0; k < option.dim; k++){
                    for(int l = 0; l < global.bc.N_D_DoF; l++)
                      if(option.dim * j + k == global.bc.fixed_dof[l]) count++;
                    if(count == 0)
                       global.subdomain.global_residual_force[j][k] -= global.subdomain.Global_K[DoF_free * (option.dim * j + k) + option.dim * i] * fixed_u_inc;
                    count = 0;
                  }
                }
            }
            // y方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 3 == 0)
            {
                type_num[1] = global.bc.Dirichlet_type[i] % 10000;
                type_num[1] = (type_num[1] - type_num[0]) / 100;
                fixed_u_inc = fixed_deformation(option.time, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[1]);
                fixed_u_inc -= fixed_deformation(option.time_old, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[1]);
                global.subdomain.global_residual_force[i][1] = fixed_u_inc;
                #if 1
                for(int j = 0; j < global.subdomain.N_point; j++){
                  for(int k = 0; k < option.dim; k++){
                    for(int l = 0; l < global.bc.N_D_DoF; l++)
                      if(option.dim * j + k == global.bc.fixed_dof[l]) count++;
                    if(count == 0)
                       global.subdomain.global_residual_force[j][k] -= global.subdomain.Global_K[DoF_free * (option.dim * j + k) + option.dim * i + 1] * fixed_u_inc;
                    count = 0;
                  }
                }
                #endif
            }
            // z方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 5 == 0)
            {
                type_num[2] = (global.bc.Dirichlet_type[i] - type_num[1] * 100 - type_num[0]) / 10000;
                fixed_u_inc = fixed_deformation(option.time, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[2]);
                fixed_u_inc -= fixed_deformation(option.time_old, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[2]);
                global.subdomain.global_residual_force[i][2] = fixed_u_inc;
                #if 1
                for(int j = 0; j < global.subdomain.N_point; j++){
                  for(int k = 0; k < option.dim; k++){
                    for(int l = 0; l < global.bc.N_D_DoF; l++)
                      if(option.dim * j + k == global.bc.fixed_dof[l]) count++;
                    if(count == 0)
                       global.subdomain.global_residual_force[j][k] -= global.subdomain.Global_K[DoF_free * (option.dim * j + k) + option.dim * i + 2] * fixed_u_inc;
                    count = 0;
                  }
                }
                #endif
            }
        }
    }
    else
    {
        for (int i = 0; i < global.subdomain.N_point; i++)
        {
            // x方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 2 == 0)
                global.subdomain.global_residual_force[i][0] = 0.0;
            // y方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 3 == 0)
                global.subdomain.global_residual_force[i][1] = 0.0;
            // z方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 5 == 0)
                global.subdomain.global_residual_force[i][2] = 0.0;
        }
    }
}

// ディリクレ条件に応じて接線剛性マトリクスと残差を書き換える (full-matrix形式) //
void ImposeDirichletTangentialMatrix(){
	int DoF_free = option.dim * global.subdomain.N_point;  // 拘束を考慮しない自由度

	for(int i=0;i<global.subdomain.N_point;i++){
    // x方向の変位が固定されているとき //
    if(global.bc.fixed_dir[i] % 2 == 0){
      // 第dim*i行の第dim*i列目を1,それ以外を0とする //
      for(int j=0;j<DoF_free;j++){
          global.subdomain.Global_K[DoF_free * (option.dim * i) + j] = 0.;
          global.subdomain.Global_K[DoF_free * j + (option.dim * i)] = 0.;
      }
      global.subdomain.Global_K[DoF_free * (option.dim * i) + option.dim * i] = 1.0;
    }
    // y方向の変位が固定されているとき //
    if(global.bc.fixed_dir[i] % option.dim == 0){
      // 第dim*i+1行の第dim*i+1列目を1,それ以外を0とする //
      for(int j=0;j<DoF_free;j++){
        global.subdomain.Global_K[DoF_free * (option.dim * i + 1) + j] = 0.;
        global.subdomain.Global_K[DoF_free * j + (option.dim * i + 1)] = 0.;
      }
      global.subdomain.Global_K[DoF_free*(option.dim * i + 1) + option.dim * i + 1] = 1.0;
    }
    // z方向の変位が固定されているとき //
    if(global.bc.fixed_dir[i] % 5 == 0){
      // 第dim*i+2行の第dim*i+2列目を1,それ以外を0とする //
       for(int j=0;j<DoF_free;j++){
        global.subdomain.Global_K[DoF_free * (option.dim * i + 2) + j] = 0.;
        global.subdomain.Global_K[DoF_free * j + (option.dim * i + 2)] = 0.;
      }
      global.subdomain.Global_K[DoF_free * (option.dim * i + 2) + option.dim * i + 2] = 1.0;
    }
	}
}



void assemble_matrix_and_vector_for_Dirichlet(double *K_u, double *residual){
  int *checker_vector;
  int matrix_count = 0;
  int residual_count = 0;
  int DoF_free = option.dim * global.subdomain.N_point;

  if((checker_vector = (int *)calloc(DoF_free, sizeof(int)))== NULL){
    printf("Checker_vector's Memory is not enough\n");
    exit(-1);
  }
  
  //Dirichlet境界条件に相当する自由度に-1をチェック
  for(int i = 0; i < global.bc.N_D_DoF; i++){
    checker_vector[global.bc.fixed_dof[i]] = -1;
  }
  
  //全体の残差ベクトルを求解用に縮退
  for(int i = 0; i < global.subdomain.N_point; i++){
    for(int j = 0; j < option.dim; j++){
      if(checker_vector[option.dim * i + j] != -1){
        residual[residual_count] += global.subdomain.global_residual_force[i][j];
        residual_count++;
      }
    }
  }

  //変位増分に残差ベクトルの値を代入
  for(int i = 0; i < global.bc.N_D_DoF; i++){
      global.subdomain.displacement_increment[global.bc.fixed_dof[i] / 3][global.bc.fixed_dof[i] % 3] += global.subdomain.global_residual_force[global.bc.fixed_dof[i] / 3][global.bc.fixed_dof[i] % 3];
  }

  //全体剛性マトリクスを求解用に縮退
  for(int i = 0; i < DoF_free; i++){
      for(int j = 0; j < DoF_free; j++){
        if(checker_vector[i] != -1 && checker_vector[j] != -1){
          K_u[matrix_count] = global.subdomain.Global_K[DoF_free * i + j];
          matrix_count++;
        }
      }
  }


  free(checker_vector);

}

//残差ベクトルと係数マトリクスを求解用にアセンブリ(線形弾性用)
void assemble_matrix_and_vector_for_Dirichlet_Linear(double *K_u, double *deformation, double *residual){
  int DoF_free = option.dim * global.subdomain.N_point;
  int flag = 0;
  int count = 0;

  //残差ベクトルからディリクレ境界条件が反映される自由度を削除
  for(int i = 0; i < global.subdomain.N_point; i++){
    for(int j = 0; j < option.dim; j++){
      for(int k = 0; k < global.bc.N_D_DoF; k++)
        if(option.dim * i + j == global.bc.fixed_dof[k]){
          flag = 1;
          break;
        }

        if(flag != 1){
          residual[count] = global.subdomain.global_residual_force[i][j];
          count++;
        }else{
          deformation[option.dim * i + j] += global.subdomain.global_residual_force[i][j];
        }

        flag = 0;
    }
  }
  count = 0;

  //係数マトリクスからディリクレ境界条件が反映される自由度を削除
  for(int i = 0; i < DoF_free; i++){
    for(int j = 0; j < DoF_free; j++){
      for(int k = 0; k < global.bc.N_D_DoF; k++){
        if(i == global.bc.fixed_dof[k] || j == global.bc.fixed_dof[k]){ 
          flag = 1;
          break;
        }
      }
      if(flag != 1){
        K_u[count] = global.subdomain.Global_K[DoF_free * i + j];
        count++;
      }
      flag = 0;
    }
  }
}