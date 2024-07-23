#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"external_force.h"
#include"scalar.h"
#include"vector.h"
#include"GetGaussPoints.h"
#include"b_matrix.h"

extern Global global;
extern Option option;

#define NUMBER_OF_NODE_IN_SUBDOMAIN 8
#define NUMBER_OF_NODE_IN_FACE 4

//外力ベクトルの更新
void update_external_force(int time){
    FILE *fp_debug;
    char FILE_name[128];

    int N_qu = 1;                           //積分点数
    double xyz[3];                          //積分点の座標
    double t_force[3];                      //トラクションのベクトル    
    double NT[60][3];                       //形状関数
    double face_node_XYZ[4][3];             //面を構成する頂点座標
    double X[27], w[27];                    //ガウス点と重み関数
    double mapping_parameter;               //座標変換のための因子
    double *current_point_XYZ;              //現配置のポイントの座標
    double subdomain_external_force[60];    //サブドメイン単位での外力ベクトル
    double lambda;                          //(ステップ数) / (最終ステップ数)     

    lambda = ((double)time + 1.0) / (double)option.N_timestep;

    //外力ベクトルをゼロ処理
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                global.subdomain.global_external_force[i][j] = 0.;
            }
        }

    //形状関数を計算するためのpoint現在座標を計算
    if((current_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:current_point_XYZ's memory is not enough\n");
        exit(-1);
    }
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                        + global.subdomain.displacement[i][j]
                                        + global.subdomain.displacement_increment[i][j];
        }
    }
   
   //ガウス積分点と重みを計算
    Gauss_points_and_weighting_factors(N_qu, X, w);

    //要素ごとの外力ベクトルを計算 ([N]^T{t})
    for(int face = 0; face < global.bc.N_t_face; face++){
        int N_support = global.subdomain.support_offset[global.bc.traction_point[face] + 1] - global.subdomain.support_offset[global.bc.traction_point[face]];

        generate_current_node_of_face(face_node_XYZ, global.bc.traction_face[face], global.bc.traction_point[face]);

        //要素外力ベクトルの計算
        for(int s = 0; s < N_qu; s++){
            for(int t = 0; t < N_qu; t++){
                //物理空間座標→正規化座標に変換するためのスカラー値を計算
                mapping_parameter = calc_mapping_parameter(global.bc.traction_face[face], global.bc.traction_point[face], s, t, X);

                //物理座標におけるガウス点の座標を計算
                generate_gauss_point_coordinate(s, t, face_node_XYZ, X, xyz);
            
                //トラクションを計算
                traction(xyz[0], xyz[1], xyz[2], t_force, global.bc.traction_type[face]);

                //形状関数を計算
                calc_shape(xyz, option.dim, global.bc.traction_point[face], current_point_XYZ, global.subdomain.support_offset, NT);

                //要素外力ベクトルの計算
                for(int i = 0; i < option.dim * (N_support + 1); i++){
                    double subdomain_external_force_i = 0.;
                    for(int j = 0; j < option.dim; j++){
                        subdomain_external_force_i += NT[i][j] * t_force[j];
                    }
                    subdomain_external_force[i] = subdomain_external_force_i * mapping_parameter * w[s] * w[t] * lambda;
                }
                
            }
        }

        //要素ごとの外力ベクトルを全体の外力ベクトルにアセンブル
        assemble_vector(global.bc.traction_point[face], global.subdomain.global_external_force, subdomain_external_force);
        
    }

    #if 0
            fp_debug = fopen("global_external_vector.dat","w");
            
            for(int i = 0; i < global.subdomain.N_point; i++){
                fprintf(fp_debug,"%5d    ", i);
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+15.14e  ", global.subdomain.global_external_force[i][j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
    #endif

    free(current_point_XYZ);
}

//トラクションの計算
void traction(double x1, double x2, double x3, double* t, int type){
	if(type == 0){
    // x=1の面に課す荷重 //
		t[0] = 0.0;  t[1] = 0.0;  t[2] = 1.0;
	}
  if(type == 1){
    // x=0の面に課す荷重 //
		t[0] = 0.0000067/1.499;  t[1] = -0.00020167/1.499;  t[2] = 0.0;
	}
  else if(type == 2){
    // y=1の面に課す荷重 //
		t[0] = 0.000201/1.499;  t[1] = -0.0000067/1.499;  t[2] = 0.0;
	}
  else if(type == 3){
    // y=0の面に課す荷重 //
		t[0] = -0.000201/1.499;  t[1] = 0.0000067/1.499;  t[2] = 0.0;
	}
}

