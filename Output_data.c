#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"Output_data.h"

extern Global global;
extern Option option;

void Output_data(int time_step){
    FILE *fp_deformation;
    FILE *fp_stress;
    FILE *fp_strain;
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

    snprintf(FILE_name, 128, "Data_Files_Output/Output_strain_time%d.dat", time_step);
    fp_stress = fopen(FILE_name,"w");
    if(fp_stress == NULL){
        printf("file not open\n");
        exit(-1);
    }
    fprintf(fp_stress,"point number / strain xx yy zz zy yz zx\n");
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
            fprintf(fp_strain, "%+15.14e    ", global.subdomain.trial_elastic_strains[i][j]);
        fprintf(fp_strain, "\n");
    }
    fclose(fp_strain);
}