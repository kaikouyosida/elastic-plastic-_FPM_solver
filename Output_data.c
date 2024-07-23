#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"Output_data.h"
#include"scalar.h"
#include"vector.h"

#define NUMBER_OF_NODE_IN_FACE 4
#define NUMBER_OF_NODE_IN_SUBDOMAIN 8

extern Global global;
extern Option option;

//応力節点データ等の出力
void Output_data(int time_step){
    FILE *fp_deformation;
    FILE *fp_stress;
    FILE *fp_strain;
    FILE *fp_F_grad;
    FILE *fp_volume;
    char FILE_name[128];

    snprintf(FILE_name, 128, "Data_Files_Output/Output_deformation_time%d.dat", time_step);
    fp_deformation = fopen(FILE_name,"w");
    if(fp_deformation == NULL){
        printf("file not open\n");
        exit(-1);
    }
    fprintf(fp_deformation,"point number / deformation  X Y Z\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_deformation, "%7d    ", i);
        for(int j = 0; j < option.dim; j++)
            fprintf(fp_deformation, "%+15.14e    ", global.subdomain.displacement[i][j] + global.subdomain.displacement_increment[i][j]);
        fprintf(fp_deformation, "\n");
    }
    fclose(fp_deformation);

    snprintf(FILE_name, 128, "Data_Files_Output/Output_Cauchy_stress_time%d.dat", time_step);
    fp_stress = fopen(FILE_name,"w");
    if(fp_stress == NULL){
        printf("file not open\n");
        exit(-1);
    }
    fprintf(fp_stress,"point number / stress xx yy zz zy yz zx\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_stress, "%7d    ", i);
        for(int j = 0; j < 6; j++)
            fprintf(fp_stress, "%+15.14e    ", global.subdomain.current_stresses[i][j]);
        fprintf(fp_stress, "\n");
    }
    fclose(fp_stress);

    snprintf(FILE_name, 128, "Data_Files_Output/Output_elastic_strain_time%d.dat", time_step);
    fp_strain = fopen(FILE_name,"w");
    if(fp_strain == NULL){
        printf("file not open\n");
        exit(-1);
    }
    fprintf(fp_strain,"point number / elastic strain xx yy zz zy yz zx\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_strain, "%7d    ", i);
        for(int j = 0; j < 6; j++)
            fprintf(fp_strain, "%+15.14e    ", global.subdomain.current_elastic_strains[i][j]);
        fprintf(fp_strain, "\n");
    }
    fclose(fp_strain);

    snprintf(FILE_name, 128, "Data_Files_Output/Output_deformation_gradient_time%d.dat", time_step);
    fp_F_grad = fopen(FILE_name, "w");
    if(fp_F_grad == NULL){
        printf("File is not open\n");
        exit(-1);
    }
    fprintf(fp_F_grad, "point number / deformation gradient xx xy xz yx yy yz zx zy zz\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_F_grad,  "%7d    ", i);
        for(int j = 0; j < option.dim ; j++){
            for(int k = 0; k < option.dim; k++){
                fprintf(fp_F_grad, "%+15.14e    ", global.subdomain.current_deformation_gradients[j][k][i]);
            }
        }
        fprintf(fp_F_grad, "\n");
    }
    fclose(fp_F_grad);
    snprintf(FILE_name, 128, "Data_Files_Output/Output_volumes_time%d.dat", time_step);
    fp_volume = fopen(FILE_name, "w");
    if(fp_volume == NULL){
        printf("File is not open\n");
        exit(-1);
    }
    fprintf(fp_F_grad, "point number / subdomain's volume\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        double volume = calc_subdomain_volume(i);
        fprintf(fp_volume, "%7d     %+15.14e\n", i, volume);
    }
    fclose(fp_volume);
}

//出力用の節点変位の計算
void update_nodal_coordinate(){
    double face_node_XYZ[4][3];
    double nodal_displacement[3];
    
    for(int node = 0; node < global.subdomain.N_node; node++){
        for(int i = 0; i < option.dim; i++)
            nodal_displacement[i] = 0.;

        for(int point = 0; point < global.subdomain.N_point; point++){
            int N_ar_point = global.subdomain.ar_node_offset[node + 1] - global.subdomain.ar_node_offset[node];
            for(int i = 0; i < NUMBER_OF_NODE_IN_SUBDOMAIN; i++){
                if(node == global.subdomain.subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN * point + i]){
                    for(int j = 0; j < option.dim; j++){
                        nodal_displacement[j] += global.subdomain.nodal_displacement_sd[point][i][j] / (double)N_ar_point;
                    }
                }
            }
        }

        for(int i = 0; i < option.dim; i++)
            global.subdomain.nodal_displacements[node][i] = nodal_displacement[i];
    }

}

//paraviewデータの出力
void paraview_node_data(int time_step){
    FILE *fp_paraview;
    double stresses[9];
    char FILE_name[128];
    
    fp_paraview = fopen("paraview/Outline_for_paraview.vtk",  "w");
    if(fp_paraview == NULL){
        printf("File is not open\n");
        exit(-1);
    }
    fprintf(fp_paraview, "# vtk DataFile Version 4.1\n");
    fprintf(fp_paraview,"FPM / 3D / hexa elements\n");
    fprintf(fp_paraview, "ASCII\n");
    fprintf(fp_paraview, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp_paraview, "POINTS %d double\n", global.subdomain.N_node);

    for(int i = 0; i < global.subdomain.N_node; i++){
        for(int j = 0; j < option.dim ; j++){
            fprintf(fp_paraview, "%+15.14e  ", global.subdomain.node_XYZ[option.dim * i + j]);
        }  
        fprintf(fp_paraview, "\n");
    }
    fprintf(fp_paraview, "CELLS %d %d\n", global.subdomain.N_point, global.subdomain.N_point * (NUMBER_OF_NODE_IN_SUBDOMAIN + 1));
    
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_paraview,  "%5d ", NUMBER_OF_NODE_IN_SUBDOMAIN);
        for(int j = 0; j < NUMBER_OF_NODE_IN_SUBDOMAIN; j++){
            fprintf(fp_paraview,  "%5d ", global.subdomain.subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN * i + j]);
        }
        fprintf(fp_paraview,  "\n");
    }
    fclose(fp_paraview);

    snprintf(FILE_name, 128,"paraview/paraview_time_deformation%d.vtk", time_step);
    fp_paraview = fopen(FILE_name,"w");
    if(fp_paraview == NULL){
        printf("File is not open\n");
        exit(-1);
    }
    fprintf(fp_paraview, "# vtk DataFile Version 4.1\n");
    fprintf(fp_paraview,"FPM / 3D / hexa elements\n");
    fprintf(fp_paraview, "ASCII\n");
    fprintf(fp_paraview, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp_paraview, "POINTS %d double\n", global.subdomain.N_node);

    for(int i = 0; i < global.subdomain.N_node; i++){
        for(int j = 0; j < option.dim ; j++){
            fprintf(fp_paraview, "%+15.14e  ", global.subdomain.node_XYZ[option.dim * i + j] + global.subdomain.nodal_displacements[i][j]);
        }  
        fprintf(fp_paraview, "\n");
    }
    fprintf(fp_paraview, "CELLS %d %d\n", global.subdomain.N_point, global.subdomain.N_point * (NUMBER_OF_NODE_IN_SUBDOMAIN + 1));
    
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_paraview,  "%5d ", NUMBER_OF_NODE_IN_SUBDOMAIN);
        for(int j = 0; j < NUMBER_OF_NODE_IN_SUBDOMAIN; j++){
            fprintf(fp_paraview,  "%5d ", global.subdomain.subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN * i + j]);
        }
        fprintf(fp_paraview,  "\n");
    }
    
    fprintf(fp_paraview,"CELL_TYPES %d\n", global.subdomain.N_point);
    for(int i = 0 ; i < global.subdomain.N_point; i++)
        fprintf(fp_paraview,  "12\n");
    fprintf(fp_paraview,"POINT_DATA %d\n", global.subdomain.N_node);
    fprintf(fp_paraview, "VECTORS node_deformation double\n");
    for(int i = 0; i < global.subdomain.N_node; i++){
        for(int j = 0; j < option.dim ; j++){
            fprintf(fp_paraview, "%+15.14e  ",  global.subdomain.nodal_displacements[i][j]);
        }  
        fprintf(fp_paraview, "\n");
    }
    fclose(fp_paraview);


    snprintf(FILE_name, 128,"paraview/paraview_time_stress%d.vtk", time_step);
    fp_paraview = fopen(FILE_name,"w");
    if(fp_paraview == NULL){
        printf("File is not open\n");
        exit(-1);
    }
    fprintf(fp_paraview, "# vtk DataFile Version 4.1\n");
    fprintf(fp_paraview,"FPM / 3D / hexa elements\n");
    fprintf(fp_paraview, "ASCII\n");
    fprintf(fp_paraview, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp_paraview, "POINTS %d double\n", global.subdomain.N_node);

    for(int i = 0; i < global.subdomain.N_node; i++){
        for(int j = 0; j < option.dim ; j++){
            fprintf(fp_paraview, "%+15.14e  ", global.subdomain.node_XYZ[option.dim * i + j] + global.subdomain.nodal_displacements[i][j]);
        }  
        fprintf(fp_paraview, "\n");
    }
    fprintf(fp_paraview, "CELLS %d %d\n", global.subdomain.N_point, global.subdomain.N_point * (NUMBER_OF_NODE_IN_SUBDOMAIN + 1));
    
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_paraview,  "%5d ", NUMBER_OF_NODE_IN_SUBDOMAIN);
        for(int j = 0; j < NUMBER_OF_NODE_IN_SUBDOMAIN; j++){
            fprintf(fp_paraview,  "%5d ", global.subdomain.subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN * i + j]);
        }
        fprintf(fp_paraview,  "\n");
    }
    
    fprintf(fp_paraview,"CELL_TYPES %d\n", global.subdomain.N_point);
    for(int i = 0 ; i < global.subdomain.N_point; i++)
        fprintf(fp_paraview,  "12\n");
    fprintf(fp_paraview,"POINT_DATA %d\n", global.subdomain.N_node);
    fprintf(fp_paraview, "TENSORS node_stress double\n");
    for(int node = 0; node < global.subdomain.N_node; node++){
        for(int i = 0; i < 9; i++)
            stresses[i] = 0.;
        int N_ar_node = global.subdomain.ar_node_offset[node + 1] - global.subdomain.ar_node_offset[node];
        for(int i = 0; i < N_ar_node; i++){
            int num_ar_node = global.subdomain.ar_node[global.subdomain.ar_node_offset[node] + i];
            stresses[0] += global.subdomain.stresses[num_ar_node][0];
            stresses[1] += global.subdomain.stresses[num_ar_node][3];
            stresses[2] += global.subdomain.stresses[num_ar_node][5];
            stresses[3] += global.subdomain.stresses[num_ar_node][3];
            stresses[4] += global.subdomain.stresses[num_ar_node][1];
            stresses[5] += global.subdomain.stresses[num_ar_node][4];
            stresses[6] += global.subdomain.stresses[num_ar_node][5];
            stresses[7] += global.subdomain.stresses[num_ar_node][4];
            stresses[8] += global.subdomain.stresses[num_ar_node][2];
        }
        for(int i = 0; i < 9; i++)
            stresses[i] /= (double)N_ar_node;
        fprintf(fp_paraview,"%+15.14e %+15.14e %+15.14e\n", stresses[0], stresses[1], stresses[2]);
        fprintf(fp_paraview,"%+15.14e %+15.14e %+15.14e\n", stresses[3], stresses[4], stresses[5]);
        fprintf(fp_paraview,"%+15.14e %+15.14e %+15.14e\n", stresses[6], stresses[7], stresses[8]);
    }
    fclose(fp_paraview);
}