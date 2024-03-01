#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"external_force.h"
#include"scalar.h"
#include"GetGaussPoints.h"
#include"b_matrix.h"

extern Global global;
extern Option option;

void update_external_force(int time){

    FILE *fp_debug;
    char FILE_name[128];
    double xyz[3];
    double t_force[3];
    double NT[60][3];
    int face_node[4];
    double face_node_XYZ[4][3]; 
    double X[27], w[27];
    double jacobian;    
    int N_qu = 1;
    double *point_XYZ;
    double *node_XYZ;
    double subdomain_external_force[60];
    double lambda;

    lambda = ((double)time + 1.0) / (double)option.N_timestep;

    //形状関数を計算するためのpoint, node現在座標を計算
    if((point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:point_XYZ's memory is not enough\n");
        exit(-1);
    }
    if((node_XYZ = (double *)calloc(option.dim * global.subdomain.N_node, sizeof(double))) == NULL){
        printf("Error:node_XYZ's memory is not enough\n");
        exit(-1);
    }
    
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                        + global.subdomain.displacement[i][j]
                                        + global.subdomain.displacement_increment[i][j];
        }
    }
    for(int i = 0; i < global.subdomain.N_node; i++){
        for(int j = 0; j < option.dim; j++){
            node_XYZ[option.dim * i + j] = global.subdomain.node_XYZ[option.dim * i + j]
                                        + global.subdomain.nodal_displacements[i][j]
                                        + global.subdomain.nodal_displacement_increments[i][j];
        }
    }
    
    //外力ベクトルをゼロ処理
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            global.subdomain.global_external_force[i][j] = 0.;
        }
    }

    Gauss_points_and_weighting_factors(N_qu, X, w);

    for(int face = 0; face < global.bc.N_t_face; face++){
        int N_support = global.subdomain.support_offset[global.bc.traction_point[face] + 1] - global.subdomain.support_offset[global.bc.traction_point[face]];
        
        jacobian = calc_surface_area(global.bc.traction_face[face]) / 4.0;

        for(int i = 0; i < 4; i++)
            face_node[i] = global.subdomain.node[global.subdomain.vertex_offset[global.bc.traction_face[face]] + i];
        
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < option.dim; j++)
                face_node_XYZ[i][j] = node_XYZ[option.dim * face_node[i] + j];
        
        for(int s = 0; s < N_qu; s++){
            for(int t = 0; t < N_qu; t++){
            
                for(int i = 0; i < option.dim; i++)
                    xyz[i] = 0.25 * (1.0 - X[s]) * (1.0 - X[t]) * face_node_XYZ[0][i]
                            + 0.25 * (1.0 - X[s]) * (1.0 + X[t]) * face_node_XYZ[1][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 + X[t]) * face_node_XYZ[2][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 - X[t]) * face_node_XYZ[3][i];
            
                traction(xyz[0], xyz[1], xyz[2], t_force, global.bc.traction_type[face]);
                calc_shape(xyz, option.dim, global.bc.traction_point[face], point_XYZ, global.subdomain.support_offset, NT);

                for(int i = 0; i < option.dim * (N_support + 1); i++){
                    double subdomain_external_force_i = 0.;
                    for(int j = 0; j < option.dim; j++){
                        subdomain_external_force_i += NT[i][j] * t_force[j];
                    }
                    subdomain_external_force[i] = subdomain_external_force_i * jacobian * w[s] * w[t] * lambda;
                }
                
            }
        }
        //要素ごとの外力ベクトルを全体の外力ベクトルにアセンブル
                for(int i = 0; i < N_support; i++){
                    for(int j = 0; j < option.dim; j++){
                        global.subdomain.global_external_force[global.subdomain.support[global.subdomain.support_offset[global.bc.traction_point[face]] + i]][j]
                            += subdomain_external_force[option.dim * (i + 1) + j];
                    }
                }
            
                for(int i = 0; i < option.dim; i++)
                    global.subdomain.global_external_force[global.bc.traction_point[face]][i] 
                        += subdomain_external_force[i];
    }
    #if 0
            global.count++;
            snprintf(FILE_name, 128,"Data_Files_Output/debag%d.dat", global.count);
            fp_debug = fopen(FILE_name,"w");
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+4.3e  ", global.subdomain.global_external_force[i][j]);
                }
                fprintf(fp_debug, "\n");
            }
            fprintf(fp_debug, "\n");
    #endif
    
    
    free(node_XYZ);
    free(point_XYZ);

}
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

