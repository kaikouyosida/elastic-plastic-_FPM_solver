#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

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
    double node_displacement[100][3];
    int num;
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
            fprintf(fp_deformation, "%+15.14e    ", global.subdomain.displacement[i][j]);
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
            fprintf(fp_stress, "%+15.14e    ", global.subdomain.stresses[i][j]);
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
            fprintf(fp_strain, "%+15.14e    ", global.subdomain.elastic_strains[i][j]);
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
    snprintf(FILE_name, 128, "Data_Files_Output/Output_Equivalent_plastic_strain_time%d.dat", time_step);
    fp_strain = fopen(FILE_name,"w");
    if(fp_strain == NULL){
        printf("file not open\n");
        exit(-1);
    }
    fprintf(fp_strain,"point number / strain\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_strain, "%5d     %+15.14e\n", i, global.subdomain.equivalent_plastic_strains[i]);
    }
    fclose(fp_strain);

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

    snprintf(FILE_name, 128, "Data_Files_Output/Output_yield_stress%d.dat", time_step);
    fp_stress = fopen(FILE_name, "w");
    if(fp_stress == NULL){
        printf("File is not open\n");
        exit(-1);
    }
    fprintf(fp_stress, "point num / yield stress\n");
    for(int i = 0; i < global.subdomain.N_point; i++){
        fprintf(fp_stress, "%7d    %+15.14e\n", i, global.subdomain.yield_stresses[i]);
    }
    fclose(fp_stress);

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

    snprintf(FILE_name, 128, "Data_Files_Output/Output_nodal_deformation%d.dat", time_step);
    fp_deformation = fopen(FILE_name, "w");
    if(fp_deformation == NULL){
        printf("File is not open\n");
        exit(-1);
    }
    fprintf(fp_deformation, "node number / node_xyz\n");
    for(int node = 0; node < global.subdomain.N_node; node++){
        fprintf(fp_deformation, "%5d    ", node);
        num = extract_node_number(node, node_displacement);
        for(int i = 0; i < num; i++){
            for(int j = 0; j < option.dim; j++){
                fprintf(fp_deformation, "%+15.14e   ", node_displacement[i][j]);
            }
            fprintf(fp_deformation, "||");
        }
        fprintf(fp_deformation, "\n");
    }
    fclose(fp_deformation);

}

//応力を抽出して出力
void calc_extract_component(){
    FILE *fp_extract;
    FILE *fp_gnuData;
    char ss[128];
    double stress[10000][6];
    double strain[10000][6];
    double u[10000][3];
    double d_matrix[6][6];
    double equivarent_stress[10000];
    double yield_stress[10000];
    double nodal_equivarent_stress[10000];
   
   #if 1
    fp_extract = fopen("Data_Files_Output/circle_0905/Output_cauchy_stress_time39.dat", "r");
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
    fp_extract = fopen("Data_Files_Output/circle_0905/Output_Plastic_strain_time39.dat", "r");
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

    fp_extract = fopen("Data_Files_Output/circle_0905/Output_Deformation_time39.dat", "r");
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

    fp_extract = fopen("Data_Files_Output/circle_0905/Output_yield_stress39.dat", "r");
    if(fp_extract == NULL){
        printf("Error:File is not open\n");
        exit(-1);
    }
    fscanf(fp_extract, "%*[^\n]\n", ss);
    for(int i = 0; i < global.subdomain.N_point; i++){
        fscanf(fp_extract, "%*d %lf\n", &yield_stress[i]);
    }
    fclose(fp_extract);
    #endif

    #if 1
    for(int i = 0; i < global.subdomain.N_point; i++){
        double hydro_static_pressure = 1.0 / 3.0 * (stress[i][0] + stress[i][1] + stress[i][2]);
        if(fabs(global.subdomain.point_XYZ[3*i+1]) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*i+2] - 1.0) < 1.0e-5){
            printf("%5d %+15.14e %+15.14e %+15.14e %+15.14e\n", i, global.subdomain.point_XYZ[3*i], stress[i][0], stress[i][1], stress[i][2]);
        }
    }
    #endif
    #if 1
    for(int i = 0; i < global.subdomain.N_point; i++){
        equivarent_stress[i] = calc_equivalent_stress(stress[i]);
    }
    for(int i = 0; i < global.subdomain.N_point; i++){
        if(fabs(global.subdomain.point_XYZ[3*i+1]) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*i+2] - 1.0) < 1.0e-5){
            printf("%5d %+15.14e %+15.14e\n", i, global.subdomain.point_XYZ[3*i],  equivarent_stress[i]);
        }
    }
    #endif
    #if 1
    for(int i = 0; i < global.subdomain.N_point; i++){
        equivarent_stress[i] = calc_equivalent_stress(strain[i]);
    }
    for(int i = 0; i < global.subdomain.N_point; i++){
        if(fabs(global.subdomain.point_XYZ[3*i+1]) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*i+2] - 1.0) < 1.0e-5){
            printf("%5d %+15.14e %+15.14e\n", i, global.subdomain.point_XYZ[3*i], 3.0 / 2.0 * equivarent_stress[i]);
        }
    }
    #endif

    fp_gnuData = fopen("gnuplot/gnuData.txt", "w");
    for(int i = 0; i < global.subdomain.N_point; i++){
        if(fabs(global.subdomain.point_XYZ[3*i+1]) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*i+2] - 1.0) < 1.0e-5){
            fprintf(fp_extract, "%+15.14e        %+15.14e\n", global.subdomain.point_XYZ[3*i], global.subdomain.point_XYZ[3*i], stress[i][1]);
        }
    }
    fclose(fp_gnuData);
    
 


}

//Paraview出力用の節点変位の計算（各サブドメインの節点変位の平均を計算）
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

int extract_node_number(int node , double (*nodal_displacement)[3]){
    int count = 0;
    int flag = -1;
    int temp;
    int num;
    
    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i = 0; i < NUMBER_OF_NODE_IN_SUBDOMAIN; i++){
            if(node == global.subdomain.subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN * point + i]){
               flag = 1; temp = i;
               break; 
            }
        }
        
        if(flag == 1){
            for(int i = 0; i < option.dim; i++){
                nodal_displacement[count][i] = global.subdomain.nodal_displacement_sd[point][temp][i];
            }
            count++; 
        }
        flag = -1;
    }
    num = count + 1;
    return num;
}

//paraviewデータの出力
void paraview_node_data(int time_step){
    FILE *fp_paraview;
    double stresses[9];
    double strain[9];
    double F_grad[9];
    char FILE_name[128];
    double Equivalent_strain[30000];

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

    snprintf(FILE_name, 128,"paraview/Paraview_time_deformation_gradient%d.vtk", time_step);
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
    fprintf(fp_paraview, "TENSORS node_deformation_gradient double\n");
    for(int node = 0; node < global.subdomain.N_node; node++){
        for(int i = 0; i < 9; i++)
            F_grad[i] = 0.;
        int N_ar_node = global.subdomain.ar_node_offset[node + 1] - global.subdomain.ar_node_offset[node];
        for(int i = 0; i < N_ar_node; i++){
            for(int j = 0; j < 3; j++){
                for(int k = 0; k < 3; k++){
                    F_grad[3 * j + k] += global.subdomain.deformation_gradients[j][k][global.subdomain.ar_node[global.subdomain.ar_node_offset[node] + i]] / (double)N_ar_node;
                }
            }
        }
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                fprintf(fp_paraview, "%+15.14e  ", F_grad[3 * i + j]);
            }
            fprintf(fp_paraview, "\n");
        }
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

    snprintf(FILE_name, 128,"paraview/Paraview_time_strain%d.vtk", time_step);
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
    fprintf(fp_paraview, "TENSORS node_strain double\n");
    for(int node = 0; node < global.subdomain.N_node; node++){
        for(int i = 0 ; i < 9; i++)
            strain[i] = 0.;
        int N_ar_node = global.subdomain.ar_node_offset[node + 1] - global.subdomain.ar_node_offset[node];
        for(int i = 0; i < N_ar_node; i++){
            
            int num_ar_node = global.subdomain.ar_node[global.subdomain.ar_node_offset[node] + i];
            strain[0] += global.subdomain.elastic_strains[num_ar_node][0] + global.subdomain.current_plastic_strains[num_ar_node][0];
            strain[1] += global.subdomain.elastic_strains[num_ar_node][3] + global.subdomain.current_plastic_strains[num_ar_node][3];   
            strain[2] += global.subdomain.elastic_strains[num_ar_node][5] + global.subdomain.current_plastic_strains[num_ar_node][5];
            strain[3] += global.subdomain.elastic_strains[num_ar_node][3] + global.subdomain.current_plastic_strains[num_ar_node][3];
            strain[4] += global.subdomain.elastic_strains[num_ar_node][1] + global.subdomain.current_plastic_strains[num_ar_node][1];
            strain[5] += global.subdomain.elastic_strains[num_ar_node][4] + global.subdomain.current_plastic_strains[num_ar_node][4];
            strain[6] += global.subdomain.elastic_strains[num_ar_node][5] + global.subdomain.current_plastic_strains[num_ar_node][5];
            strain[7] += global.subdomain.elastic_strains[num_ar_node][4] + global.subdomain.current_plastic_strains[num_ar_node][4];
            strain[8] += global.subdomain.elastic_strains[num_ar_node][2] + global.subdomain.current_plastic_strains[num_ar_node][2];
        }
        for(int i = 0; i < 9; i++)
            strain[i] /= (double)N_ar_node;
        
        fprintf(fp_paraview,"%+15.14e %+15.14e %+15.14e\n", strain[0], strain[1], strain[2]);
        fprintf(fp_paraview,"%+15.14e %+15.14e %+15.14e\n", strain[3], strain[4], strain[5]);
        fprintf(fp_paraview,"%+15.14e %+15.14e %+15.14e\n", strain[6], strain[7], strain[8]);
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

    snprintf(FILE_name, 128,"paraview/Paraview_time_Equivalent_stress%d.vtk", time_step);
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
    fprintf(fp_paraview, "SCALARS Equivalent_stress double\n");
    fprintf(fp_paraview, "LOOKUP_TABLE default\n");
    for(int i = 0; i < global.subdomain.N_node; i++){
        double equivalent_stress = 0;
        int N_ar_node = global.subdomain.ar_node_offset[i+1] - global.subdomain.ar_node_offset[i];
        for(int j = 0; j < N_ar_node; j++){
            equivalent_stress += calc_equivalent_stress(global.subdomain.stresses[global.subdomain.ar_node[global.subdomain.ar_node_offset[i]+j]]);
        }
        equivalent_stress /= N_ar_node;
        fprintf(fp_paraview, "%+15.14e\n", equivalent_stress);
    }
    fclose(fp_paraview);

    snprintf(FILE_name, 128, "paraview/Paraview_time_Equivalent_strain%d.vtk", time_step);
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
    fprintf(fp_paraview, "SCALARS Equivalent_strain double\n");
    fprintf(fp_paraview, "LOOKUP_TABLE default\n");

    for(int i = 0; i < global.subdomain.N_point; i++){
        Equivalent_strain[i] = 2.0 / 3.0 * calc_equivalent_stress(global.subdomain.current_plastic_strains[i]);
    }
    for(int i = 0 ;i < global.subdomain.N_node; i++){
        int N_ar_point = global.subdomain.ar_node_offset[i+1] - global.subdomain.ar_node_offset[i];
        double Equivalent_strain_i = 0.;
        for(int j = 0; j < N_ar_point; j++){
            Equivalent_strain_i += Equivalent_strain[global.subdomain.ar_node[global.subdomain.ar_node_offset[i] + j]] / N_ar_point;
        }
        fprintf(fp_paraview, "%+15.14e\n", Equivalent_strain_i);
    }

    fclose(fp_paraview);

    snprintf(FILE_name, 128, "paraview/Paraview_point_deformation_gradient%d.vtk", time_step);
    fp_paraview = fopen(FILE_name, "w");
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

    fprintf(fp_paraview,"CELL_DATA %d\n", global.subdomain.N_point);
    fprintf(fp_paraview, "TENSORS point_deformation_gradient double\n");
    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i = 0; i < 9; i++)
            F_grad[i] = 0.;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                F_grad[3 * i + j] = global.subdomain.deformation_gradients[i][j][point];
            }
        }
        
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                fprintf(fp_paraview, "%+15.14e  ", F_grad[3 * i + j]);
            }
            fprintf(fp_paraview, "\n");
        }
    }

    fclose(fp_paraview);
}
