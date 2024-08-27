#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"Output_data.h"
#include"scalar.h"
#include"vector.h"
#include"stress.h"
#include"d_matrix.h"

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

    snprintf(FILE_name, 128, "Data_Files_Output/Output_Deformation_time%d.dat", time_step);
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

    snprintf(FILE_name, 128, "Data_Files_Output/Output_cauchy_stress_time%d.dat", time_step);
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

    snprintf(FILE_name, 128, "Data_Files_Output/Output_Elastic_strain_time%d.dat", time_step);
    fp_strain = fopen(FILE_name,"w");
    if(fp_strain == NULL){
        printf("file not open\n");
        exit(-1);
    }
    fprintf(fp_strain,"point number / Elastic strain xx yy zz zy yz zx\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_strain, "%7d    ", i);
        for(int j = 0; j < 6; j++)
            fprintf(fp_strain, "%+15.14e    ", global.subdomain.current_elastic_strains[i][j]);
        fprintf(fp_strain, "\n");
    }
    fclose(fp_strain);

    #if 1
    snprintf(FILE_name, 128, "Data_Files_Output/Output_Plastic_strain_time%d.dat", time_step);
    fp_strain = fopen(FILE_name,"w");
    if(fp_strain == NULL){
        printf("file not open\n");
        exit(-1);
    }
    fprintf(fp_strain,"point number / plastic strain xx yy zz zy yz zx\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_strain, "%7d    ", i);
        for(int j = 0; j < 6; j++)
            fprintf(fp_strain, "%+15.14e    ", global.subdomain.current_plastic_strains[i][j]);
        fprintf(fp_strain, "\n");
    }
    fclose(fp_strain);
    #endif

    snprintf(FILE_name, 128, "Data_Files_Output/Output_total_strain_time%d.dat", time_step);
    fp_strain = fopen(FILE_name,"w");
    if(fp_strain == NULL){
        printf("file not open\n");
        exit(-1);
    }
    fprintf(fp_strain,"point number / total strain xx yy zz zy yz zx\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_strain, "%7d    ", i);
        for(int j = 0; j < 6; j++)
            fprintf(fp_strain, "%+15.14e    ", global.subdomain.current_elastic_strains[i][j] + global.subdomain.current_plastic_strains[i][j]);
        fprintf(fp_strain, "\n");
    }
    fclose(fp_strain);

    if(option.solver_type == 1){
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
}

//応力を抽出して出力
void calc_extract_component(){
    FILE *fp_extract;
    char ss[128];
    double stress[1000][6];
    double strain[1000][6];
    double u[1000][3];
    double d_matrix[6][6];
    double equivarent_stress[10000];
    double nodal_equivarent_stress[10000];

    fp_extract = fopen("Data_Files_Output/Output_cauchy_stress_time49.dat", "r");
    if(fp_extract == NULL){
        printf("File is not open\n");
        exit(-1);
    }
    fscanf(fp_extract, "%*[^\n]\n", ss);

    for(int i = 0; i < global.subdomain.N_point; i++){
        fscanf(fp_extract, "%*d");
        for(int j = 0; j < 6; j++){
            fscanf(fp_extract, "%lf", &stress[i][j]);
        }
        fscanf(fp_extract, "\n");
    }
    fclose(fp_extract);

    fp_extract = fopen("Data_Files_Output/Output_Plastic_strain_time49.dat", "r");
    if(fp_extract == NULL){
        printf("File is not enough\n");
        exit(-1);
    }
    fscanf(fp_extract, "%*[^\n]\n", ss);

    for(int i = 0; i < global.subdomain.N_point; i++){
        fscanf(fp_extract, "%*d");
        for(int j = 0; j < 6; j++){
            fscanf(fp_extract, "%lf", &strain[i][j]);
        }
        fscanf(fp_extract, "\n");
    }
    fclose(fp_extract);

    fp_extract = fopen("Data_Files_Output/Output_Deformation_time49.dat", "r");
    if(fp_extract == NULL){
        printf("Error:Memory is not enough\n");
        exit(-1);
    }
    fscanf(fp_extract, "%*[^\n]\n", ss);
    for(int i = 0; i < global.subdomain.N_point; i++){
        fscanf(fp_extract, "%*d");
        for(int j = 0; j < option.dim; j++){
            fscanf(fp_extract, "%lf", &u[i][j]);
        }
        fscanf(fp_extract, "\n");
    }
    fclose(fp_extract);
    #if 0
    for(int i = 0; i < global.subdomain.N_point; i++){
        double hydro_static_pressure = 1.0 / 3.0 * (stress[i][0] + stress[i][1] + stress[i][2]);
        if(fabs(global.subdomain.point_XYZ[3*i+1]) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*i+2] - 1.0) < 1.0e-5){
            printf("%5d %+15.14e %+15.14e\n", i, global.subdomain.point_XYZ[3*i], hydro_static_pressure);
        }
    }
    #endif
    #if 1
    for(int i = 0; i < global.subdomain.N_point; i++){
        equivarent_stress[i] = calc_equivalent_stress(strain[i]);
    }
    for(int i = 0; i < global.subdomain.N_point; i++){
        if(fabs(global.subdomain.point_XYZ[3*i+1]) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*i+2] - 1.0) < 1.0e-5){
            printf("%5d %+15.14e %+15.14e\n", i, global.subdomain.point_XYZ[3*i], equivarent_stress[i]);
        }
    }
    #endif

    #if 0
    double error = 0.;
    generateElasticDMatrix(d_matrix);

    double stress_ext[6] = {0, 0, 1.0, 0, 0, 0};
    double strain_ext[6] = {-0.00015, -0.00015, 0.0005, 0, 0, 0};
    for(int i = 0; i < global.subdomain.N_point; i++){
        double volume = calc_initial_subdomain_volume(i);
        double integrated_part = 0.;
        for(int j = 0; j < 6; j++){
            double strain_D_j = 0.;
            for(int k = 0; k < 6; k++){
                strain_D_j += (strain_ext[k]) * d_matrix[k][j];
            }
            integrated_part += strain_D_j * (strain_ext[j]);
        }
        error += integrated_part * volume;
    }
     
    error = sqrt(error);
    printf("%+15.14e\n", error);
    #endif
    #if 0
    for(int i = 0; i < global.subdomain.N_point; i++){
        equivarent_stress[i] = calc_equivalent_stress(stress[i]);
    }
    for(int i = 0; i < global.subdomain.N_node; i++){
        int N_ar_node = global.subdomain.ar_node_offset[i + 1] - global.subdomain.ar_node_offset[i];
        double temp = 0.;
        for(int j = 0; j < N_ar_node; j++){
            temp += equivarent_stress[global.subdomain.ar_node[global.subdomain.ar_node_offset[i] + j]] / N_ar_node;
        }
        printf("%+15.14e\n", temp);
    }
    #endif

}

//出力用の節点変位の計算
void update_nodal_coordinate(){
    int ar_point[27];
    double nodal_displacement[3];
    
    if(option.solver_type == 1){
        for(int node = 0; node < global.subdomain.N_node; node++){
            int N_ar_point = global.subdomain.ar_node_offset[node + 1] - global.subdomain.ar_node_offset[node];
            for(int i = 0; i < option.dim; i++)
                nodal_displacement[i] = 0.;

            for(int point = 0; point < global.subdomain.N_point; point++){
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
    }else{
        for(int node = 0; node < global.subdomain.N_node; node++){
            int N_ar_point = global.subdomain.ar_node_offset[node + 1] - global.subdomain.ar_node_offset[node];
            for(int i = 0; i < option.dim; i++)
                nodal_displacement[i] = 0.;   
            for(int i = 0; i < N_ar_point; i++)
                ar_point[i] = global.subdomain.ar_node[global.subdomain.ar_node_offset[node] + i];
            for(int i = 0; i < N_ar_point; i++){
                for(int j = 0; j < option.dim; j++){
                    nodal_displacement[j] += global.subdomain.displacement[ar_point[i]][j];
                }
            }
            for(int i = 0; i < option.dim; i++)
                global.subdomain.nodal_displacements[node][i] = nodal_displacement[i] / (double)N_ar_point;
        }
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

    snprintf(FILE_name, 128,"paraview/Paraview_time_deformation%d.vtk", time_step);
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


    snprintf(FILE_name, 128,"paraview/Paraview_time_stress%d.vtk", time_step);
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

    snprintf(FILE_name, 128, "paraview/Paraview_time_point_deformation%d.vtk", time_step);
    fp_paraview = fopen(FILE_name, "w");
    if(fp_paraview == NULL){
        printf("File is not open\n");
        exit(-1);
    }
    fprintf(fp_paraview, "# vtk DataFile Version 4.1\n");
    fprintf(fp_paraview,"FPM / 3D / point clowd\n");
    fprintf(fp_paraview, "ASCII\n");
    fprintf(fp_paraview, "DATASET POLYDATA\n");
    fprintf(fp_paraview, "POINTS %d double\n", global.subdomain.N_point);
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim ; j++){
            fprintf(fp_paraview, "%+15.14e  ", global.subdomain.point_XYZ[option.dim * i + j] + global.subdomain.displacement[i][j]);
        }  
        fprintf(fp_paraview, "\n");
    }

    fprintf(fp_paraview,"POINT_DATA %d\n", global.subdomain.N_point);
    fprintf(fp_paraview, "VECTORS point_deformation double\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim ; j++){
            fprintf(fp_paraview, "%+15.14e  ",  global.subdomain.displacement[i][j]);
        }  
        fprintf(fp_paraview, "\n");
    }
    fclose(fp_paraview);

}