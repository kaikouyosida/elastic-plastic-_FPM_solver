#include <stdio.h>
#include <stdlib.h>
#include "type.h"
#include "model.h"

extern Global global;
extern Option option;
extern SS_CURVE ss_curve;

//*Inputファイルのデータを読み込むファイル*/

void Init_model()
{
    FILE *fp_analysis_set;
    fp_analysis_set = fopen("Analysis_condition_setting.dat", "r");
    if (fp_analysis_set == NULL)
    {
        printf("Error1: setting data is not open\n");
        exit(-1);
    }
    fscanf(fp_analysis_set, "%*[^\n]\n");
    fscanf(fp_analysis_set, "%d\n", &option.dim);
    fgetc(fp_analysis_set);
    fscanf(fp_analysis_set, "%*[^\n]\n");
    fscanf(fp_analysis_set, "%lf\n", &option.magni);
    fgetc(fp_analysis_set);
    fscanf(fp_analysis_set, "%*[^\n]\n");
    fscanf(fp_analysis_set, "%lf\n", &option.Delta_time);
    fgetc(fp_analysis_set);
    fscanf(fp_analysis_set, "%*[^\n]\n");
    fscanf(fp_analysis_set, "%d\n", &option.N_timestep);
    fgetc(fp_analysis_set);
    fscanf(fp_analysis_set, "%*[^\n]\n");
    fscanf(fp_analysis_set, "%d\n", &option.time_output);
    fgetc(fp_analysis_set);
    fscanf(fp_analysis_set, "%*[^\n]\n");
    fscanf(fp_analysis_set, "%lf\n", &option.NR_tol);
    fgetc(fp_analysis_set);
    fscanf(fp_analysis_set, "%*[^\n]\n");
    fscanf(fp_analysis_set, "%d\n", &option.stop_count);
    fclose(fp_analysis_set);
    // 終了時刻を計算 //
    option.time_end = option.Delta_time * (double)option.N_timestep;

    FILE *fp_material;
    fp_material = fopen("material_constant_setting.dat", "r");
    if (fp_material == NULL)
    {
        printf("Error2: material data is not open\n");
        exit(-1);
    }
    fscanf(fp_material, "%*[^\n]\n");
    fscanf(fp_material, "%*[^\n]\n");
    fscanf(fp_material, "%d\n", &global.material.material_type);
    fgetc(fp_material);
    fscanf(fp_material, "%*[^\n]\n");
    fscanf(fp_material, "%lf\n", &global.material.rho);
    fgetc(fp_material);
    fscanf(fp_material, "%*[^\n]\n");
    fscanf(fp_material, "%lf\n", &global.material.E_mod);
    fgetc(fp_material);
    fscanf(fp_material, "%*[^\n]\n");
    fscanf(fp_material, "%lf\n", &global.material.nu_mod);
    fgetc(fp_material);
    fscanf(fp_material, "%*[^\n]\n");
    fscanf(fp_material, "%d\n", &ss_curve.num);
    if ((ss_curve.plastic_strain = (double *)calloc(ss_curve.num, sizeof(double))) == NULL)
    {
        printf("Error:Plastic strain's memory is not enough\n");
        exit(-1);
    }
    if ((ss_curve.stress = (double *)calloc(ss_curve.num, sizeof(double))) == NULL)
    {
        printf("Error:Stress's memory is not enough\n");
        exit(-1);
    }
    for (int i = 0; i < 3; i++)
    fscanf(fp_material, "%lf %lf\n", &ss_curve.plastic_strain[i], &ss_curve.stress[i]);
    fscanf(fp_material, "%lf\n", &global.material.kinematic_hardening_fractions);
    fclose(fp_material);

    global.material.penalty = option.magni * global.material.E_mod;

    FILE *fp_subdomain;
    fp_subdomain = fopen("Data_Files_Input/Input_subdomain.dat", "r");
    if (fp_subdomain == NULL)
    {
        printf("Error:Subdomain file is not open\n");
        exit(-1);
    }

    fscanf(fp_subdomain, "%*s %d\n", &global.subdomain.N_point);
    fscanf(fp_subdomain, "%*s %d\n", &global.subdomain.sum_N_node);
    fscanf(fp_subdomain, "%*s %d\n", &global.subdomain.N_node);
    fclose(fp_subdomain);

    if ((global.subdomain.node_XYZ = (double *)calloc(option.dim * global.subdomain.N_node, sizeof(double))) == NULL)
    {
        printf("Error:Node_XYZ's memory is not enough\n");
        exit(-1);
    }
    if ((global.subdomain.subdomain_node = (int *)calloc(option.dim * global.subdomain.sum_N_node, sizeof(int))) == NULL)
    {
        printf("Error:Subdomain node's memory is not enough\n");
        exit(-1);
    }

    fp_subdomain = fopen("Data_Files_Input/Input_subdomain.dat", "r");
    if (fp_subdomain == NULL)
    {
        printf("Error: subdomain file not open\n");
        exit(-1);
    }
    for (int i = 0; i < 4; i++)
        fscanf(fp_subdomain, "%*[^\n]\n");
    for (int i = 0; i < global.subdomain.N_point; i++)
    {
        fscanf(fp_subdomain, "%*d %d\n", &global.buf);
        for (int j = 0; j < global.buf; j++)
            fscanf(fp_subdomain, "%d", &global.subdomain.subdomain_node[global.buf * i + j]);
        fgetc(fp_subdomain);
    }
    fgetc(fp_subdomain);
    fscanf(fp_subdomain, "%*[^\n]\n");

    for (int i = 0; i < global.subdomain.N_node; i++)
    {
        fscanf(fp_subdomain, "%*d");
        for (int j = 0; j < option.dim; j++)
        {
            fscanf(fp_subdomain, "%lf", &global.subdomain.node_XYZ[option.dim * i + j]);
        }
        fgetc(fp_subdomain);
    }
    fclose(fp_subdomain);

    // ポイントの初期配置における座標を読み込む //
    FILE *fp_point;
    if ((global.subdomain.point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL)
    {
        printf("Error:Point_XYZ's memory is not enough\n");
        exit(-1);
    }
    fp_point = fopen("Data_Files_Input/Input_point_coordinate.dat", "r");
    if (fp_point == NULL)
    {
        printf("Error: point file not open\n");
        exit(-1);
    }
    for (int i = 0; i < 2; i++)
        fscanf(fp_point, "%*[^\n]\n");

    for (int i = 0; i < global.subdomain.N_point; i++)
    {
        fscanf(fp_point, "%*d");
        for (int j = 0; j < option.dim; j++)
        {
            fscanf(fp_point, "%lf", &global.subdomain.point_XYZ[option.dim * i + j]);
        }
        fgetc(fp_point);
    }
    fclose(fp_point);

    // サブドメインの中心の座標を読み込む //
    FILE *fp_center;
    if ((global.subdomain.center_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL)
    {
        printf("Error:center_XYZ's memory is not enough\n");
        exit(-1);
    }
    fp_center = fopen("Data_Files_Input/Input_center_of_subdomain.dat", "r");
    if (fp_center == NULL)
    {
        printf("Error: center file not open\n");
        exit(-1);
    }
    for (int i = 0; i < 2; i++)
        fscanf(fp_center, "%*[^\n]\n");

    for (int i = 0; i < global.subdomain.N_point; i++)
    {
        fscanf(fp_center, "%*d");
        for (int j = 0; j < option.dim; j++)
        {
            fscanf(fp_center, "%lf", &global.subdomain.center_XYZ[option.dim * i + j]);
        }
        fgetc(fp_center);
    }
    fclose(fp_center);

    // 面のデータを開き、サブドメインと面の関係を読み込む //
    FILE *fp_face;
    fp_face = fopen("Data_Files_Input/Input_face_of_subdomain.dat", "r");
    if (fp_face == NULL)
    {
        printf("Error: face of subdomain file not open\n");
        exit(-1);
    }

    fscanf(fp_face, "%*[^\n]\n");
    fscanf(fp_face, "%*s %d\n", &global.subdomain.sum_N_face);
    fscanf(fp_face, "%*[^\n]\n");
    for (int i = 0; i < global.subdomain.N_point; i++)
        fscanf(fp_face, "%*[^\n]\n");
    fgetc(fp_face);
    fscanf(fp_face, "%*s %d\n", &global.subdomain.N_face);
    fscanf(fp_face, "%*s %d\n", &global.subdomain.sum_N_vertex);
    fclose(fp_face);

    if ((global.subdomain.face_offset = (int *)calloc(global.subdomain.N_point + 1, sizeof(int))) == NULL)
    {
        printf("Error:Face_offset's memory is not enough\n");
        exit(-1);
    }
    if ((global.subdomain.face = (int *)calloc(global.subdomain.sum_N_face, sizeof(int))) == NULL)
    {
        printf("Error:Face's memory is not enough\n");
        exit(-1);
    }
    if ((global.subdomain.vertex_offset = (int *)calloc(global.subdomain.N_face + 1, sizeof(int))) == NULL)
    {
        printf("Error:Vertex_offset's memory is not enough\n");
        exit(-1);
    }
    if ((global.subdomain.node = (int *)calloc(global.subdomain.sum_N_vertex, sizeof(int))) == NULL)
    {
        printf("Error:Node's memory is not enough\n");
        exit(-1);
    }

    fp_face = fopen("Data_Files_Input/Input_face_of_subdomain.dat", "r");
    if (fp_face == NULL)
    {
        printf("Error: Face_of_subdomain file not open\n");
        exit(-1);
    }
    for (int i = 0; i < 3; i++)
        fscanf(fp_face, "%*[^\n]\n");

    global.subdomain.face_offset[0] = 0;
    for (int i = 0; i < global.subdomain.N_point; i++)
    {
        fscanf(fp_face, "%*d %d", &global.buf);
        for (int j = 0; j < global.buf; j++)
        {
            fscanf(fp_face, "%d", &global.subdomain.face[global.count]);
            global.count += 1;
        }
        fgetc(fp_face);

        global.subdomain.face_offset[i + 1] = global.subdomain.face_offset[i] + global.buf;
    }
    global.count = 0;
    fgetc(fp_face);
    for (int i = 0; i < 3; i++)
        fscanf(fp_face, "%*[^\n]\n");

    global.subdomain.vertex_offset[0] = 0;
    for (int i = 0; i < global.subdomain.N_face; i++)
    {
        fscanf(fp_face, "%*d %d", &global.buf);
        for (int j = 0; j < global.buf; j++)
        {
            fscanf(fp_face, "%d", &global.subdomain.node[global.count]);
            global.count += 1;
        }
        fgetc(fp_face);

        global.subdomain.vertex_offset[i + 1] = global.subdomain.vertex_offset[i] + global.buf;
    }
    fclose(fp_face);
    global.count = 0;

    // 内部境界の情報を読み込む //
    FILE *fp_int_boundary;
    fp_int_boundary = fopen("Data_Files_Input/Input_internal_boundary.dat", "r");
    if (fp_int_boundary == NULL)
    {
        printf("Error: int_boundary file not open\n");
        exit(-1);
    }
    fscanf(fp_int_boundary, "%*s %d\n", &global.subdomain.N_int_boundary);
    fclose(fp_int_boundary);

    if ((global.subdomain.pair_point_ib = (int *)calloc(2 * global.subdomain.N_int_boundary, sizeof(int))) == NULL)
    {
        printf("Error:Pair_point_ib's memory is not enough\n");
        exit(-1);
    }
    if ((global.subdomain.shared_face = (int *)calloc(global.subdomain.N_int_boundary, sizeof(int))) == NULL)
    {
        printf("Error:Shared_face's memory is not enough\n");
        exit(-1);
    }

    fp_int_boundary = fopen("Data_Files_Input/Input_internal_boundary.dat", "r");
    if (fp_int_boundary == NULL)
    {
        printf("Error :internal boundary file not open\n");
        exit(-1);
    }
    for (int i = 0; i < 2; i++)
        fscanf(fp_int_boundary, "%*[^\n]\n");

    for (int i = 0; i < global.subdomain.N_int_boundary; i++)
        fscanf(fp_int_boundary, "%d %d %d\n", &global.subdomain.pair_point_ib[2 * i], &global.subdomain.pair_point_ib[2 * i + 1], &global.subdomain.shared_face[i]);
    fclose(fp_int_boundary);

    // サポートドメインを読み込む //
    FILE *fp_support;
    fp_support = fopen("Data_Files_Input/Input_support_domain.dat", "r");
    if (fp_support == NULL)
    {
        printf("Error: support_domain file is not open\n");
        exit(-1);
    }

    fscanf(fp_support, "%*[^\n]\n");
    fscanf(fp_support, "%*s %d\n", &global.subdomain.sum_N_support);
    fclose(fp_support);

    if ((global.subdomain.support_offset = (int *)calloc(global.subdomain.N_point + 1, sizeof(int))) == NULL)
    {
        printf("Error:Support_offset's memory is not enough\n");
        exit(-1);
    }
    if ((global.subdomain.support = (int *)calloc(global.subdomain.sum_N_support, sizeof(int))) == NULL)
    {
        printf("Error:Support's memory is not enough\n");
        exit(-1);
    }

    fp_support = fopen("Data_Files_Input/Input_support_domain.dat", "r");
    if (fp_support == NULL)
    {
        printf("Error:support_domain file is not open\n");
        exit(-1);
    }
    for (int i = 0; i < 3; i++)
        fscanf(fp_support, "%*[^\n]\n");

    global.subdomain.support_offset[0] = 0;
    for (int i = 0; i < global.subdomain.N_point; i++)
    {
        fscanf(fp_support, "%*d %d", &global.buf);

        global.subdomain.support_offset[i + 1] = global.subdomain.support_offset[i] + global.buf;
        for (int j = 0; j < global.buf; j++)
        {
            fscanf(fp_support, "%d", &global.subdomain.support[global.count]);
            global.count += 1;
        }
        fgetc(fp_support);
    }

    fclose(fp_support);
    global.count = 0;

    // 境界条件に関する情報を読み込む//
    // ディリクレ条件を読み込む //
    FILE *fp_Dirichlet;
    if ((global.bc.fixed_dir = (int *)calloc(global.subdomain.N_point, sizeof(int))) == NULL)
    {
        printf("Memory error:fixed_dir memory is not enough\n");
        exit(-1);
    }
    for (int i = 0; i < global.subdomain.N_point; i++)
        global.bc.fixed_dir[i] = 1;
    if ((global.bc.Dirichlet_type = (int *)calloc(global.subdomain.N_point, sizeof(int))) == NULL)
    {
        printf("Error:Dirichlet_type's memory is not enough\n");
        exit(-1);
    }

    fp_Dirichlet = fopen("Data_Files_Input/Input_Dirichlet_condition.dat", "r");
    if (fp_Dirichlet == NULL)
    {
        printf("Error: Dhirichlet file is not open\n");
        exit(-1);
    }
    fscanf(fp_Dirichlet, "%*s %d\n", &global.bc.N_D_DoF);
    if((global.bc.fixed_dof = (int *)calloc(global.bc.N_D_DoF, sizeof(int))) == NULL){
        printf("Error: Dhirichet dof's memory is not enough\n");
        exit(-1);
    }
    fscanf(fp_Dirichlet, "%*[^\n]\n");

    for (int i = 0; i < global.bc.N_D_DoF; i++)
    {
        fscanf(fp_Dirichlet, "%d %d %d %*d\n", &global.buf, &global.bc.buf_dir, &global.bc.buf_type);
        global.bc.fixed_dof[i] = option.dim * global.buf + global.bc.buf_dir;
        if (global.bc.buf_dir == 0)
        {
            global.bc.fixed_dir[global.buf] *= 2;
            global.bc.Dirichlet_type[global.buf] += global.bc.buf_type;
        }
        if (global.bc.buf_dir == 1)
        {
            global.bc.fixed_dir[global.buf] *= 3;
            global.bc.Dirichlet_type[global.buf] += 100 * global.bc.buf_type;
        }
        if (global.bc.buf_dir == 2)
        {
            global.bc.fixed_dir[global.buf] *= 5;
            global.bc.Dirichlet_type[global.buf] += 10000 * global.bc.buf_type;
        }
    }
    fclose(fp_Dirichlet);
    //for(int i = 0; i < global.bc.N_D_DoF; i++)
        //printf("%5d %5d\n", i, global.bc.fixed_dof[i]);
    //printf("ok!\n");
    //exit(-1);

    // トラクション条件を読み込む //
    FILE *fp_traction;
    fp_traction = fopen("Data_Files_Input/Input_traction_condition.dat", "r");
    if (fp_traction == NULL)
    {
        printf("Error:traction file not open\n");
        exit(-1);
    }
    fscanf(fp_traction, "%*s %d\n", &global.bc.N_t_face);
    fclose(fp_traction);

    if ((global.bc.traction_point = (int *)calloc(global.bc.N_t_face, sizeof(int))) == NULL)
    {
        printf("Error:Traction_point's memory is not enough\n");
        exit(-1);
    }
    if ((global.bc.traction_type = (int *)calloc(global.bc.N_t_face, sizeof(int))) == NULL)
    {
        printf("Error:Traction_type's memory is not enough\n");
        exit(-1);
    }
    if ((global.bc.traction_face = (int *)calloc(global.bc.N_t_face, sizeof(int))) == NULL)
    {
        printf("Error:Traction_face's memory is not enough\n");
        exit(-1);
    }

    fp_traction = fopen("Data_Files_Input/Input_traction_condition.dat", "r");
    if (fp_traction == NULL)
    {
        printf("Error:traction file not open\n");
        exit(-1);
    }
    for (int i = 0; i < 2; i++)
        fscanf(fp_traction, "%*[^\n]\n");

    for (int i = 0; i < global.bc.N_t_face; i++)
        fscanf(fp_traction, "%d %d %d\n", &global.bc.traction_point[i], &global.bc.traction_type[i], &global.bc.traction_face[i]);

    fclose(fp_traction);

    // ポイント周りの頂点のデータを読み込む
    //  各ノードが属するサブドメインの番号を読み込む //
    if ((global.subdomain.ar_node_offset = (int *)calloc(global.subdomain.N_node + 1, sizeof(int))) == NULL)
    {
        printf("Error:Around_node_offset's memory is not enough\n");
        exit(-1);
    }
    if ((global.subdomain.ar_node = (int *)calloc(global.subdomain.sum_N_vertex, sizeof(int))) == NULL)
    {
        printf("Error:Around_node's Memory is not enough\n");
        exit(-1);
    }
    
    FILE *fp_point_ar_node;
    fp_point_ar_node = fopen("Data_Files_Input/Input_point_around_node.dat", "r");
    if (fp_point_ar_node == NULL)
    {
        printf("Error:Around_node file is not open\n");
        exit(-1);
    }
    for (int i = 0; i < 3; i++)
        fscanf(fp_point_ar_node, "%*[^\n]\n");

    global.subdomain.ar_node_offset[0] = 0;
    for (int i = 0; i < global.subdomain.N_node; i++)
    {
        fscanf(fp_point_ar_node, "%*d %d", &global.buf);
        for (int j = 0; j < global.buf; j++)
        {
            fscanf(fp_point_ar_node, "%d", &global.subdomain.ar_node[global.count]);
            global.count += 1;
        }
        fgetc(fp_point_ar_node);

        global.subdomain.ar_node_offset[i + 1] = global.subdomain.ar_node_offset[i] + global.buf;
    }

    fclose(fp_point_ar_node);
    global.count = 0;
}

void break_model_memory(){
    free(global.subdomain.ar_node);
    free(global.subdomain.ar_node_offset);
    free(global.bc.traction_face);
    free(global.bc.traction_type);
    free(global.bc.traction_point);
    free(global.bc.fixed_dof);
    free(global.bc.Dirichlet_type);
    free(global.bc.fixed_dir);
    free(global.subdomain.support);
    free(global.subdomain.support_offset);
    free(global.subdomain.shared_face);
    free(global.subdomain.pair_point_ib);
    free(global.subdomain.node);
    free(global.subdomain.vertex_offset);
    free(global.subdomain.face);
    free(global.subdomain.support);
    free(global.subdomain.face_offset);
    free(global.subdomain.center_XYZ);
    free(global.subdomain.point_XYZ);
    free(global.subdomain.subdomain_node);
    free(global.subdomain.node_XYZ);
    free(ss_curve.stress);
    free(ss_curve.plastic_strain);   
}