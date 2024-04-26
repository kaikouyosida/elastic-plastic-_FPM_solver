#include <stdio.h>
#include <stdlib.h>
#include "type.h"

extern Option option;
extern Global global;

double fixed_deformation(double time, double time_end, double x1, double x2, double x3, int type){
  double fixed_u = 0.;

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
void ImposeDirichretResidual(int NR_step)
{
    double fixed_u_inc = 0.; // 規定された変位の増分量
    int type_num[3];         // ディリクレ条件の番号
    double fixed_xyz[3];     // 変位を固定するポイントの位置

    if (NR_step == 1)
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
                global.subdomain.global_residual_force[3 * i] = fixed_u_inc;
            }
            // y方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 3 == 0)
            {
                type_num[1] = global.bc.Dirichlet_type[i] % 10000;
                type_num[1] = (type_num[1] - type_num[0]) / 100;
                fixed_u_inc = fixed_deformation(option.time, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[1]);
                fixed_u_inc -= fixed_deformation(option.time_old, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[1]);
                global.subdomain.global_residual_force[3 * i + 1] = fixed_u_inc;
            }
            // z方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 5 == 0)
            {
                type_num[2] = (global.bc.Dirichlet_type[i] - type_num[1] * 100 - type_num[0]) / 10000;
                fixed_u_inc = fixed_deformation(option.time, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[2]);
                fixed_u_inc -= fixed_deformation(option.time_old, option.time_end, fixed_xyz[0], fixed_xyz[1], fixed_xyz[2], type_num[2]);
                global.subdomain.global_residual_force[3 * i + 2] = fixed_u_inc;
            }
        }
    }
    else
    {
        for (int i = 0; i < global.subdomain.N_point; i++)
        {
            // x方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 2 == 0)
                global.subdomain.global_residual_force[3 * i] = 0.0;
            // y方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 3 == 0)
                global.subdomain.global_residual_force[3 * i + 1] = 0.0;
            // z方向の変位が固定されているとき //
            if (global.bc.fixed_dir[i] % 5 == 0)
                global.subdomain.global_residual_force[3 * i + 2] = 0.0;
        }
    }
}
// ディリクレ条件に応じて接線剛性マトリクスと残差を書き換える (full-matrix形式) //
void ImposeDirichletTangentialMatrix(){
	int DoF_free = 3*global.subdomain.N_point;  // 拘束を考慮しない自由度

	for(int i=0;i<global.subdomain.N_point;i++){
    // x方向の変位が固定されているとき //
    if(global.bc.fixed_dir[i] % 2 == 0){
      // 第dim*i行の第dim*i列目を1,それ以外を0とする //
      for(int j=0;j<DoF_free;j++) global.subdomain.Global_K[DoF_free*3*i + j] = 0.;
      global.subdomain.Global_K[DoF_free*3*i + 3*i] = 1.0;
    }
    // y方向の変位が固定されているとき //
    if(global.bc.fixed_dir[i] % 3 == 0){
      // 第dim*i+1行の第dim*i+1列目を1,それ以外を0とする //
      for(int j=0;j<DoF_free;j++) global.subdomain.Global_K[DoF_free*(3*i+1) + j] = 0.;
      global.subdomain.Global_K[DoF_free*(3*i+1) + 3*i+1] = 1.0;
    }
    // z方向の変位が固定されているとき //
    if(global.bc.fixed_dir[i] % 5 == 0){
      // 第dim*i+2行の第dim*i+2列目を1,それ以外を0とする //
      for(int j=0;j<DoF_free;j++) global.subdomain.Global_K[DoF_free*(3*i+2) + j] = 0.;
      global.subdomain.Global_K[DoF_free*(3*i+2) + 3*i+2] = 1.0;
    }
	}
}