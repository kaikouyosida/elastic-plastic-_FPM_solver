#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"internal_force.h"
#include"d_matrix.h"
#include"b_matrix.h"
#include"matrix.h"
#include"tensor.h"
#include"stress.h"
#include"ss_curve.h"
#include"scalar.h"
#include"GetGaussPoints.h"
#include"ImposeDirichretCondition.h"

extern Global global;
extern Option option;

void update_field_and_internal_forces(){

    FILE *fp_debug;
    char FILE_name[128];
    //fp_debug = fopen("debug_stress.dat", "w");
    //fprintf(fp_debug, "stress      /        σxx     /       σyy     /       σzz     /       σxy     /       σyz     /       σzx     \n");
    
    double subdomain_internal_force[20][3];
    double d_matrix[6][6];          //Dマトリクス
    double inverse_relative_deformation_gradient[3][3];         //相対変形勾配テンソルの逆テンソル
    double relative_deformation_gradient[3][3];                 //相対変形勾配テンソル
    double elastic_strain_tensor[3][3];                         //弾性ひずみテンソル
    double deformation_gradients[3][3];                         //変形勾配テンソル
    double current_deformation_gradients[3][3];                 //現配置の変形勾配テンソル
    double elastic_left_cauchy_green_deformations[3][3];        //弾性左コーシーグリーンテンソル
    double trial_elastic_left_cauchy_green_deformations[3][3];  //試行弾性左コーシーグリーンテンソル
    double trial_relative_stresses[6];                          //試行相対応力
    double b_t_matrix[60][6];                                   //Bマトリクスの転置
    double *latest_point_XYZ;                                   //現配置のポイント配置
    double **G;                                                 //(u = {G}uE）につかうGマトリクス
    int support[60];                                             //サポートドメイン内のポイント数
    double displacement_increment[3];                           //サポートの変位増分
    double elastic_strains[6];                                  //弾性ひずみ
    double current_elastic_strains[6];                          //現配置の弾性ひずみ
    double trial_elastic_strains[6];                            //試行弾性ひずみ
    double current_stresses[6];                                 //現配置の応力
    double back_stresses[6];                                    //背応力
    double current_back_stresses[6];                            //現配置の背応力
    double **all_stress;                                        //各サブドメインにおける応力
    double trial_relative_equivalent_stress;
    double factor;

    all_stress = matrix(global.subdomain.N_point, 6);
    
    //internal_forceをゼロ処理
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            global.subdomain.global_internal_force[i][j] = 0.;
        }
    }
    #if 0
    FILE *fp_debug_du;
    if((fp_debug_du = fopen("Data_Files_Output/Output_displacement.dat",  "r")) == NULL){
        printf("File is not enough\n");
        exit(-1);
    }
    fscanf(fp_debug_du, "%*[^\n]\n");
    for(int i = 0; i < global.subdomain.N_point; i++)
        fscanf(fp_debug_du,  "%*d %lf %lf %lf\n", &global.subdomain.displacement_increment[i][0], &global.subdomain.displacement_increment[i][1], &global.subdomain.displacement_increment[i][2]);
    //for(int i = 0; i < global.subdomain.N_point; i++)
        //printf("%+15.14e %+15.14e %+15.14e\n", global.subdomain.displacement_increment[i][0], global.subdomain.displacement_increment[i][1], global.subdomain.displacement_increment[i][2]);
        
    #endif

    for(int point = 0; point < global.subdomain.N_point; point++){
       
        double kinematic_hardening_fractions = global.material.kinematic_hardening_fractions;
        int N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];

        for(int i = 0 ; i < N_support; i++)
            support[i] = global.subdomain.support[global.subdomain.support_offset[point] + i];

        //サブドメインの内力ベクトルをゼロ処理
        for(int i = 0; i < N_support + 1; i++){
            for(int j = 0; j < option.dim; j++){
                subdomain_internal_force[i][j] = 0.;
            }
        }

        //配列の時数を一つ落とす
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                deformation_gradients[i][j] = global.subdomain.deformation_gradients[i][j][point];
                current_deformation_gradients[i][j] = global.subdomain.current_deformation_gradients[i][j][point];           
            }
        }
        for(int i = 0; i < 6; i++){
            elastic_strains[i] = global.subdomain.elastic_strains[point][i];
            current_elastic_strains[i] = global.subdomain.current_elastic_strains[point][i];
            trial_elastic_strains[i] = global.subdomain.trial_elastic_strains[point][i];
            current_stresses[i] = global.subdomain.current_stresses[point][i];
            back_stresses[i] = global.subdomain.back_stresses[point][i];
            current_back_stresses[i] = global.subdomain.current_back_stresses[point][i];
        }
        double *equivalent_plastic_strain = &global.subdomain.equivalent_plastic_strains[point];
        double *equivalent_plastic_strain_increment = &global.subdomain.equivalent_plastic_strain_increments[point];
        double *yield_stress = &global.subdomain.yield_stresses[point];
        double *current_yield_stress = &global.subdomain.current_yield_stresses[point];

        generate_linear_b_matrix(b_t_matrix, point);
        #if 0
        double varepsiron[6];
        for(int i = 0; i < 6; i++){
            double varepsiron_i = 0.;
            for(int j = 0; j < N_support; j++){
                for(int k = 0 ; k < 3; k++){
                    varepsiron_i += b_t_matrix[3*(j + 1)+k][i] * global.subdomain.displacement_increment[support[j]][k];
                }
            }
            for(int j = 0; j < 3; j++)
                varepsiron_i += b_t_matrix[j][i] * global.subdomain.displacement_increment[point][j];
            varepsiron[i] = varepsiron_i;
            printf("%+8.7e  ", varepsiron[i]);
        }
        printf("\n");
        #endif
        #if 0
        for(int i = 0; i < 6; i++){
            for(int j = 0 ;j < 3*(N_support+1); j++){
                printf("%+4.3e ", b_t_matrix[j][i]);
            }
            printf("\n");
        }
        printf("\n\n");

        for(int j = 0 ;j < 3; j++){
                printf("%+4.3e ", global.subdomain.displacement_increment[point][j]);
        }
        for(int i = 0; i < N_support; i++){
            for(int j = 0 ;j < 3; j++){
                printf("%+4.3e ", global.subdomain.displacement_increment[support[i]][j]);
            }
        }
        printf("\n\n");
        //exit(-1);
        #endif
        
        generateElasticDMatrix(d_matrix);
        
        //相対変形勾配テンソルを計算。dF=(I-d(du)/d(x+u+du))^-1//
        identify3x3Matrix(inverse_relative_deformation_gradient);
        
        if((latest_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
            printf("Error: Latest_point_XYZ's memory is not enough\n");
            exit(-1);
        }
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                latest_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                        + global.subdomain.displacement[i][j]
                        + global.subdomain.displacement_increment[i][j];
            }
        }

        G = matrix(9, option.dim * (N_support + 1));
        calc_G(option.dim, point, latest_point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);
        
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                for(int k = 0; k < N_support; k++){
                    inverse_relative_deformation_gradient[i][j] -= G[3 * i + j][3 * (k + 1) + i] * global.subdomain.displacement_increment[support[k]][i];
                }
            }
            for(int j = 0; j < 3; j++){
                inverse_relative_deformation_gradient[i][j] -= G[3 * i + j][i] * global.subdomain.displacement_increment[point][i];
            }
        }
        free_matrix(G);
        free(latest_point_XYZ);
        invert3x3Matrix(relative_deformation_gradient, inverse_relative_deformation_gradient);
        #if 0
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                if(i == j) 
                    relative_deformation_gradient[i][j] -= 1.0;
                printf("%+15.14e ", relative_deformation_gradient[i][j]);
            }
            printf("\n");
        }
        printf("\n");
        #endif

        //変形勾配テンソル（初期配置に対してのテンソル）
        for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        double deformation_gradient_i_j = 0.0;

                        for (int k = 0; k < 3; k++)
                            deformation_gradient_i_j
                                += relative_deformation_gradient[i][k]
                                *  deformation_gradients[k][j];

                        current_deformation_gradients[i][j]
                            = deformation_gradient_i_j;
                    }
        #if 0
        for(int i = 0; i < 3; i++){
            for(int j = 0;  j < 3; j++){
                printf("%+15.14e    ", current_deformation_gradients[i][j]);
            }
            printf("\n");
        }
        printf("\n");
        #endif
              
        //弾性左コーシーグリーンテンソルの計算([B]^e = exp(2 * {epsilon}^e))
        elastic_strain_tensor[0][0] = 2.0 * elastic_strains[0];
        elastic_strain_tensor[0][1] = 2.0 * 0.5 * elastic_strains[3];
        elastic_strain_tensor[0][2] = 2.0 * 0.5 * elastic_strains[5];
        elastic_strain_tensor[1][0] = 2.0 * 0.5 * elastic_strains[3];
        elastic_strain_tensor[1][1] = 2.0 * elastic_strains[1];
        elastic_strain_tensor[1][2] = 2.0 * 0.5 * elastic_strains[4];
        elastic_strain_tensor[2][0] = 2.0 * 0.5 * elastic_strains[5];
        elastic_strain_tensor[2][1] = 2.0 * 0.5 * elastic_strains[4];
        elastic_strain_tensor[2][2] = 2.0 * elastic_strains[2];
        
        calculateTensorExponent(elastic_left_cauchy_green_deformations, elastic_strain_tensor);
        #if 0
        double A[3][3];
        double B[3][3];
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                A[i][j] = 0.;
            }
        }
        A[0][0] = 1.0; A[1][1] = 1.0; A[2][2] = 9.0 / 4.0;
        calculateTensorExponent(B, A);
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                printf("%+7.6e  ", B[i][j]);
            }
            printf("\n");
        }
        exit(-1);
        #endif

        //試行弾性左コーシーグリーンテンソルの計算 [B]^trial = [dF] * [B]^e * [dF]^T 
        for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        double left_cauchy_green_deformation_i_j = 0.0;

                        for (int k = 0; k < 3; k++)
                        {
                            double left_cauchy_green_deformation_k_l_times_relative_deformation_gradient_j_l = 0.0;

                            for (int l = 0; l < 3; l++)
                                left_cauchy_green_deformation_k_l_times_relative_deformation_gradient_j_l
                                    += elastic_left_cauchy_green_deformations[k][l]
                                    *  relative_deformation_gradient[j][l];

                            left_cauchy_green_deformation_i_j
                                += relative_deformation_gradient[i][k]
                                *  left_cauchy_green_deformation_k_l_times_relative_deformation_gradient_j_l;
                        }

                        trial_elastic_left_cauchy_green_deformations[i][j]
                            = left_cauchy_green_deformation_i_j;
                    }
        
        //試行弾性ひずみを計算（ {epsilon}^trial = 1/2 * ln([B]^trial)　
        calculateTensorLogarithm(elastic_strain_tensor, trial_elastic_left_cauchy_green_deformations);
        #if 0
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                printf("%+6.5e ", elastic_strain_tensor[i][j]);
            }
            printf("\n");
        }
        printf("\n");
        #endif

        #if 0
        double A[3][3];
        double B[3][3];
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                A[i][j] = 0.;
            }
        }
        A[0][0] = 1.0; A[1][1] = 1.0; A[2][2] = 9.0 / 4.0;
        calculateTensorLogarithm(B, A);
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                printf("%+7.6e  ", B[i][j]);
            }
            printf("\n");
        }
        exit(-1);
        #endif
        current_elastic_strains[0] = 0.5 * elastic_strain_tensor[0][0];
        current_elastic_strains[1] = 0.5 * elastic_strain_tensor[1][1];
        current_elastic_strains[2] = 0.5 * elastic_strain_tensor[2][2];
        current_elastic_strains[3] = 0.5 * (elastic_strain_tensor[0][1]
                                            + elastic_strain_tensor[1][0]);
        current_elastic_strains[4] = 0.5 * (elastic_strain_tensor[1][2]
                                            + elastic_strain_tensor[2][1]);
        current_elastic_strains[5] = 0.5 * (elastic_strain_tensor[2][0]
                                            + elastic_strain_tensor[0][2]);
        
        for (int i = 0; i < 6; i++)
                trial_elastic_strains[i]
                    = current_elastic_strains[i];
        
        //試行応力の計算({sigma}^trial = [D] * {epsilon}^trial)
        for (int i = 0; i < 6; i++)
            {
                double stress_i = 0.0;

                for (int j = 0; j < 6; j++)
                    stress_i += d_matrix[i][j] * current_elastic_strains[j];

                current_stresses[i] = stress_i;
            }
        //printf("%5d  %+15.14e %+15.14e %+15.14e %+15.14e %+15.14e %+15.14e\n", point, current_elastic_strains[0], current_elastic_strains[1], current_elastic_strains[2], current_elastic_strains[3], current_elastic_strains[4], current_elastic_strains[5]);
        //printf("%5d  %+15.14e %+15.14e %+15.14e %+15.14e %+15.14e %+15.14e\n", point, current_stresses[0], current_stresses[1], current_stresses[2], current_stresses[3], current_stresses[4], current_stresses[5]);


        //試行相対応力の計算({sigma}^trial = {sigma}^trial - {beta})
        for (int i = 0; i < 6; i++)
                trial_relative_stresses[i]
                    = current_stresses[i] - back_stresses[i];
        //試行相対相当応力の計算(sigma^trial_e)
        trial_relative_equivalent_stress
                = calc_equivalent_stress(trial_relative_stresses);
        
        if (trial_relative_equivalent_stress <= (*yield_stress))
            {   
                
                *equivalent_plastic_strain_increment = 0.0;
                *current_yield_stress = *yield_stress;
                for (int i = 0; i < 6; i++)
                    current_back_stresses[i] = back_stresses[i];
            }else{
                printf("stress is too large!!\n");exit(-1);
                double hardening_stress_increment;
                double current_relative_hydrostatic_stress;

                //相当塑性ひずみ増分の計算
                *equivalent_plastic_strain
                    = calc_equivalent_plastic_strain_increment(trial_relative_equivalent_stress,
                                                                *equivalent_plastic_strain,
                                                                *yield_stress);
                //硬化応力増分の計算
                hardening_stress_increment
                    = get_hardening_stress((*equivalent_plastic_strain) + (*equivalent_plastic_strain_increment))
                    - get_hardening_stress((*equivalent_plastic_strain));
                
                //試行相対偏差応力の計算
                current_relative_hydrostatic_stress
                    = (1.0 / 3.0)
                    * (trial_relative_stresses[0]
                       + trial_relative_stresses[1]
                       + trial_relative_stresses[2]);
                trial_relative_stresses[0] -= current_relative_hydrostatic_stress;
                trial_relative_stresses[1] -= current_relative_hydrostatic_stress;
                trial_relative_stresses[2] -= current_relative_hydrostatic_stress;


                //最終的な弾性応力の計算
                factor
                    = (*equivalent_plastic_strain_increment)
                    * 1.5
                    / trial_relative_equivalent_stress;
                for (int i = 0; i < 3; i++)
                    current_elastic_strains[i]
                        -= factor * trial_relative_stresses[i];
                for (int i = 3; i < 6; i++)
                    current_elastic_strains[i]
                        -= 2.0 * factor * trial_relative_stresses[i];
                
                //最終的な降伏応力
                *current_yield_stress
                    = (*yield_stress)
                    + (1.0 - kinematic_hardening_fractions)
                    * hardening_stress_increment;
                //最終的な背応力の計算
                factor
                    = kinematic_hardening_fractions
                    * hardening_stress_increment
                    / trial_relative_equivalent_stress;
                for (int i = 0; i < 6; i++)
                    current_back_stresses[i]
                        = back_stresses[i]
                        + factor * trial_relative_stresses[i];   
                //最終的な応力の計算
                factor
                        = (*current_yield_stress)
                        / trial_relative_equivalent_stress;
                    for (int i = 0; i < 6; i++)
                        current_stresses[i]
                            = current_back_stresses[i]
                            + factor * trial_relative_stresses[i];
                current_stresses[0] += current_relative_hydrostatic_stress;
                current_stresses[1] += current_relative_hydrostatic_stress;
                current_stresses[2] += current_relative_hydrostatic_stress; 
            }

            //Kirchhoff応力からCauchy応力の計算
            double inverse_volume_change
                = 1.0 / calc_3x3matrix_determinant(current_deformation_gradients);
            
            for (int i = 0; i < 6; i++)
                    current_stresses[i] *= inverse_volume_change;

            //printf("%5d  %+15.14e %+15.14e %+15.14e %+15.14e %+15.14e %+15.14e\n", point, current_stresses[0], current_stresses[1], current_stresses[2], current_stresses[3], current_stresses[4], current_stresses[5]);
            //printf("inverse = %+15.14e\n", inverse_volume_change);

            #if 0
            fp_debug = fopen("Data_Files_Output/Output_stress.dat", "r");
            if(fp_debug == NULL){
                printf("File is not open\n");
                exit(-1);
            }
            fscanf(fp_debug, "%*[^\n]\n");
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < option.dim; j++){
                    if(i == point){
                        fscanf(fp_debug,"%*d %lf %lf %lf %lf %lf %lf\n", &current_stresses[0], &current_stresses[1], &current_stresses[2], &current_stresses[3], &current_stresses[4], &current_stresses[5]);
                    }else{
                        fscanf(fp_debug, "\n");
                    }
                }
            }
            fclose(fp_debug);
            #endif
            //printf("%lf %lf %lf %lf %lf %lf\n", current_stresses[0], current_stresses[1], current_stresses[2], current_stresses[3], current_stresses[4], current_stresses[5]);
            
            double volume = calc_subdomain_volume(point);
            //printf("%+15.14e\n",  volume);
            #if 1
            //内力ベクトル一項目を計算（[B]^T * {sigma} * dV）
            for(int i = 0; i < N_support + 1; i++){
                for(int j = 0; j < option.dim; j++){
                    double force_j = 0.;
                    for(int k = 0; k < 6; k++)
                        force_j += b_t_matrix[option.dim * i + j][k] * current_stresses[k];
                    
                    subdomain_internal_force[i][j] = force_j * volume;
                }
            }

            //各サブドメインで内力ベクトルをアセンブル
            for(int i = 0; i < N_support; i++){
                for(int j = 0; j < option.dim; j++){
                    global.subdomain.global_internal_force[support[i]][j] += subdomain_internal_force[i + 1][j];
                }
            }
            for(int i = 0; i < option.dim; i++)
                global.subdomain.global_internal_force[point][i] += subdomain_internal_force[0][i];
            #endif

            //各サブドメインにおける応力とひずみを記録（ペナルティ項の計算に用いる）
            for(int i = 0; i < 6; i++){
                all_stress[point][i] = current_stresses[i];
                global.subdomain.current_stresses[point][i] = current_stresses[i];
                global.subdomain.trial_elastic_strains[point][i] = trial_elastic_strains[i];
                global.subdomain.current_elastic_strains[point][i] = current_elastic_strains[i];
            }
            for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++){
                    global.subdomain.current_deformation_gradients[i][j][point] = current_deformation_gradients[i][j];
                }
            }
    }
    #if 0
    for(int i = 0; i < global.subdomain.N_point; i++){
        printf("%5d    ", i);
        for(int j = 0; j < 3; j++){
            printf("%+15.14e  ", global.subdomain.global_internal_force[i][j]);
        }
        printf("\n");
    }
    #endif 
    //fclose(fp_debug);
    
    #if 0
            global.count++;
            snprintf(FILE_name, 128,"debug_for_internal/1st_internal_vector/global_internal_vector%d.dat", global.count);
            fp_debug = fopen(FILE_name,"w");
            
            for(int i = 0; i < global.subdomain.N_point; i++){
                fprintf(fp_debug,"%5d    ", i);
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+15.14e  ", global.subdomain.global_internal_force[i][j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
    #endif 
    
    //内力ベクトルのペナルティ項を計算
    calc_internal_force_penalty(all_stress, 1);
    #if  0
    for(int i = 0; i < global.subdomain.N_point; i++){
        printf("%5d    ", i);
        for(int j = 0; j < 3; j++){
            printf("%+15.14e  ", global.subdomain.global_internal_force[i][j]);
        }
        printf("\n");
    }
    #endif 
    #if 0
            global.count++;
            snprintf(FILE_name, 128,"debug_for_internal/2nd_internal_vector/global_internal_vector%d.dat", global.count);
            fp_debug = fopen(FILE_name,"w");
            
            for(int i = 0; i < global.subdomain.N_point; i++){
                fprintf(fp_debug,"%5d    ", i);
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+15.14e  ", global.subdomain.global_internal_force[i][j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
    #endif
    #if 0
    FILE *fp_debug_du;
    if((fp_debug_du = fopen("Data_Files_Output/Output_displacement.dat",  "r")) == NULL){
        printf("File is not enough\n");
        exit(-1);
    }
    fscanf(fp_debug_du, "%*[^\n]\n");
    for(int i = 0; i < global.subdomain.N_point; i++)
        fscanf(fp_debug_du,  "%*d %lf %lf %lf\n", &global.subdomain.displacement_increment[i][0], &global.subdomain.displacement_increment[i][1], &global.subdomain.displacement_increment[i][2]);
    //for(int i = 0; i < global.subdomain.N_point; i++)
        //printf("%+15.14e %+15.14e %+15.14e\n", global.subdomain.displacement_increment[i][0], global.subdomain.displacement_increment[i][1], global.subdomain.displacement_increment[i][2]);
    #endif     
    calc_internal_force_penalty_stabilization(2);
    #if 0
            global.count++;
            snprintf(FILE_name, 128,"debug_for_internal/3rd_internal_vector/global_internal_vector%d.dat", global.count);
            fp_debug = fopen(FILE_name,"w");
            
            for(int i = 0; i < global.subdomain.N_point; i++){
                fprintf(fp_debug,"%5d    ", i);
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+15.14e  ", global.subdomain.global_internal_force[i][j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
    #endif
    #if  0
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < 3; j++){
            printf("%5d %+15.14e\n", 3*i+j, global.subdomain.global_internal_force[i][j]);
        }
    }
    #endif
    #if 0
            global.count++;
            snprintf(FILE_name, 128,"debug_for_internal/global_internal_vector%d.dat", global.count);
            fp_debug = fopen(FILE_name,"w");
            
            for(int i = 0; i < global.subdomain.N_point; i++){
                fprintf(fp_debug,"%5d    ", i);
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+15.14e  ", global.subdomain.global_internal_force[i][j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
    #endif
    //fclose(fp_debug_du);
    //exit(-1);

    free_matrix(all_stress);
}

double calc_equivalent_plastic_strain_increment(double trial_relative_equivalent_stress,
                                                double equivalent_plastic_strain,
                                                double yield_stress){
        double tolerance = 1.0E-8 * yield_stress; /* Piecewise linear SS curve should give the exact solution */
        int max_iteration_count = 100;

        double three_times_shear_modulus
        = 3.0 * 0.5 * global.material.E_mod / (1.0 + global.material.nu_mod);

        double hardening_stress = get_hardening_stress(equivalent_plastic_strain);

        double current_hardening_stress, current_hardening_modulus;
        double residual, residual_gradient;

        double dep = 0.0;

        //Solve 3 * G * dep - (se - sy - (sh (ep + dep) - sh (ep))) = 0の求解
    for (int i = 0; i < max_iteration_count; i++){
        current_hardening_stress
            = get_hardening_stress(equivalent_plastic_strain + dep);

        residual
            = three_times_shear_modulus * dep
            - (trial_relative_equivalent_stress - yield_stress
            - (current_hardening_stress - hardening_stress));

        if (fabs(residual) <= tolerance)
            return dep;

        current_hardening_modulus
            = get_hardening_stress(equivalent_plastic_strain + dep);

        residual_gradient = three_times_shear_modulus + current_hardening_modulus;

        dep -= residual / residual_gradient;
    }

#if 0
    printError("Warning: equivalent plastic strain increment calculation iteration not converged\n");
#endif

    return nan("");
}

void zero_fill_displacement_increments(){
    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i  = 0; i < option.dim; i++){
            global.subdomain.displacement_increment[point][i] = 0.;
        }
    }
}


void calc_internal_force_penalty(double **all_stress,int N_qu){
    double X[27];
    double w[27];                                               //ガウス求積に使う正規化座標と重み関数
    double xyz[3];                                              //求積点の座標
    double Ne[3][6];                                            //内部境界の法線ベクトル
    double N1T[60][3], N2T[60][3];                              //形状関数の転置
    double N1Tne[60][6], N2Tne[60][6];                          //形状関数と法線ベクトルの積
    double subdomain_internal_force[60];                        //サブドメインごとの内力ベクトル
    int face_node[4];                                           //面を構成する節点
    double face_node_XYZ[4][3];                                 //面を構成する節点の座標
    double jacobian;                                            //ヤコビアン
    double *point_XYZ;                                          //ポイントの現在の座標
    double *node_XYZ;                                           //節点の現在の座標
    
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
    Gauss_points_and_weighting_factors(N_qu, X, w);

    for(int face = 0; face < global.subdomain.N_int_boundary; face++){
        //jacobian = calc_surface_area(global.subdomain.shared_face[face]) / 4.0;

        int N1_support = global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face] + 1]
                        - global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face]];
        int N2_support = global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + 1] + 1] 
                        - global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + 1]]; 

        for(int i = 0; i < 4; i++)
            face_node[i] = global.subdomain.node[global.subdomain.vertex_offset[global.subdomain.shared_face[face]] + i];
        
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < option.dim; j++)
                face_node_XYZ[i][j] = node_XYZ[option.dim * face_node[i] + j];
                
        //内部境界上での積分
        for(int s = 0; s < N_qu; s++){
            for(int t = 0; t < N_qu; t++){
                jacobian = calc_area_change(global.subdomain.shared_face[face], s, t, X);

                //ノードの座標からガウス積分点の位置を近似、決定
                for(int i = 0; i < option.dim; i++)
                    xyz[i] = 0.25 * (1.0 - X[s]) * (1.0 - X[t]) * face_node_XYZ[0][i]
                            + 0.25 * (1.0 - X[s]) * (1.0 + X[t]) * face_node_XYZ[1][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 + X[t]) * face_node_XYZ[2][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 - X[t]) * face_node_XYZ[3][i];
                
                //形状関数の転置を計算
                calc_shape(xyz, option.dim, global.subdomain.pair_point_ib[2 * face], point_XYZ, global.subdomain.support_offset, N1T);
                calc_shape(xyz, option.dim, global.subdomain.pair_point_ib[2 * face + 1], point_XYZ, global.subdomain.support_offset, N2T);
                
                //法線ベクトルを計算
                calc_Ne(option.dim, global.subdomain.pair_point_ib[2 * face], global.subdomain.pair_point_ib[2 * face + 1]
                        , global.subdomain.shared_face[face], global.subdomain.vertex_offset, global.subdomain.node, node_XYZ, point_XYZ, Ne);
            
                //形状関数と法線ベクトルの積を計算
                for(int i = 0; i < option.dim * (N1_support + 1); i++){
                    for(int j = 0; j < 6; j++){
                        double N1Tne_ij = 0.;
                        for(int k = 0; k < option.dim; k++){
                            N1Tne_ij += N1T[i][k] * Ne[k][j];
                        }
                        N1Tne[i][j] = N1Tne_ij;
                    }
                }

                for(int i = 0; i < option.dim * (N2_support + 1); i++){
                    for(int j = 0; j < 6; j++){
                        double N2Tne_ij = 0.;
                        for(int k = 0; k < option.dim; k++) N2Tne_ij += N2T[i][k] * Ne[k][j];
                        N2Tne[i][j] = N2Tne_ij;
                    }
                }

                //[NT]{neσ}の計算→全体の内力ベクトルにアセンブル
                for(int i = 0; i < option.dim * (N1_support + 1); i++){
                    double subdomain_internal_force_i = 0.;
                    for(int j = 0; j < 6; j++){
                        subdomain_internal_force_i += -0.5 * (N1Tne[i][j] * all_stress[global.subdomain.pair_point_ib[2 * face]][j]
                                                        + N1Tne[i][j] * all_stress[global.subdomain.pair_point_ib[2 * face + 1]][j]);
                    }
                    subdomain_internal_force[i] = subdomain_internal_force_i * jacobian * w[s] * w[t];
                }
                for(int i = 0; i < N1_support; i++){
                    for(int j = 0; j < option.dim; j++){
                        global.subdomain.global_internal_force[global.subdomain.support[global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face]] + i]][j]
                        += subdomain_internal_force[option.dim * (i + 1) + j];   
                    }
                }
                for(int i = 0 ; i < option.dim; i++)
                    global.subdomain.global_internal_force[global.subdomain.pair_point_ib[2 * face]][i]
                        += subdomain_internal_force[i];

                for(int i = 0; i < option.dim * (N2_support + 1); i++){
                    double subdomain_internal_force_i = 0.;
                    for(int j = 0; j < 6; j++){
                        subdomain_internal_force_i += 0.5 * (N2Tne[i][j] * all_stress[global.subdomain.pair_point_ib[2 * face]][j]
                                                        + N2Tne[i][j] * all_stress[global.subdomain.pair_point_ib[2 * face + 1]][j]);
                    }
                    subdomain_internal_force[i] = subdomain_internal_force_i * jacobian * w[s] * w[t];
                }
    
                for(int i = 0; i < N2_support; i++){
                    for(int j = 0; j < option.dim; j++){
                        global.subdomain.global_internal_force[global.subdomain.support[global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + 1]] + i]][j]
                        += subdomain_internal_force[option.dim * (i + 1) + j];
                    }
                }

                for(int i = 0 ; i < option.dim; i++)
                    global.subdomain.global_internal_force[global.subdomain.pair_point_ib[2 * face + 1]][i]
                        += subdomain_internal_force[i];
 
            }
        }
    }
    
    free(node_XYZ);
    free(point_XYZ);
}

void calc_internal_force_penalty_stabilization(int N_qu){
    double X[27];
    double w[27];                                               //ガウス求積に使う正規化座標と重み関数
    double xyz[3];                                              //求積点の座標
    double u1[3], u2[3];                                        //試行関数
    double N1T[60][3], N2T[60][3];                              //形状関数の転置
    double N1Tne[60][6], N2Tne[60][6];                          //形状関数と法線ベクトルの積
    double subdomain_internal_force[60];                        //サブドメインごとの内力ベクトル
    int face_node[4];                                           //面を構成する節点
    double face_node_XYZ[4][3];                                 //面を構成する節点の座標
    double jacobian;                                            //ヤコビアン
    double *point_XYZ;                                          //ポイントの現在の座標
    double *node_XYZ;                                           //節点の現在の座標
    double he;                                                  //ポイント間の距離
    double eta = global.material.penalty;                       //ペナルティパラメータ

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
            point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j];
                                        + global.subdomain.displacement[i][j]
                                        + global.subdomain.displacement_increment[i][j];
        }
    }
    for(int i = 0; i < global.subdomain.N_node; i++){
        for(int j = 0; j < option.dim; j++){
            node_XYZ[option.dim * i + j] = global.subdomain.node_XYZ[option.dim * i + j];
                                        + global.subdomain.nodal_displacements[i][j]
                                        + global.subdomain.nodal_displacement_increments[i][j];
        }
    }
    
    Gauss_points_and_weighting_factors(N_qu, X, w);

    for(int face = 0; face < global.subdomain.N_int_boundary; face++){
        //jacobian = calc_surface_area(global.subdomain.shared_face[face]) / 4.0;
    
        he = distance(option.dim, global.subdomain.pair_point_ib[2 * face], global.subdomain.pair_point_ib[2 * face + 1], point_XYZ);

        int N1_support = global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face] + 1]
                        - global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face]];
        int N2_support = global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + 1] + 1] 
                        - global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + 1]]; 

        for(int i = 0; i < 4; i++)
            face_node[i] = global.subdomain.node[global.subdomain.vertex_offset[global.subdomain.shared_face[face]] + i];
        
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < option.dim; j++)
                face_node_XYZ[i][j] = node_XYZ[option.dim * face_node[i] + j];   

        for(int s = 0; s < N_qu; s++){
            for(int t = 0; t < N_qu; t++){
                jacobian = calc_area_change(global.subdomain.shared_face[face], s,t,X);
                for(int i = 0; i < option.dim; i++)
                    xyz[i] = 0.25 * (1.0 - X[s]) * (1.0 - X[t]) * face_node_XYZ[0][i]
                            + 0.25 * (1.0 - X[s]) * (1.0 + X[t]) * face_node_XYZ[1][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 + X[t]) * face_node_XYZ[2][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 - X[t]) * face_node_XYZ[3][i];
                
                //形状関数の転置を計算
                calc_shape(xyz, option.dim, global.subdomain.pair_point_ib[2 * face], point_XYZ, global.subdomain.support_offset, N1T);
                calc_shape(xyz, option.dim, global.subdomain.pair_point_ib[2 * face + 1], point_XYZ, global.subdomain.support_offset, N2T);

                //試行関数を計算
                trial_u(xyz, global.subdomain.pair_point_ib[2 * face], point_XYZ, u1);
                trial_u(xyz, global.subdomain.pair_point_ib[2 * face + 1], point_XYZ, u2);

                for(int i = 0; i < option.dim * (N1_support + 1); i++){
                    double subdomain_internal_force_i = 0.;
                    for(int j = 0; j < option.dim; j++){
                        subdomain_internal_force_i += N1T[i][j] * (u1[j] - u2[j]);
                    }
                    subdomain_internal_force[i] = subdomain_internal_force_i;
                }
                for(int i = 0; i < N1_support; i++){
                    for(int j = 0; j < option.dim; j++){
                        global.subdomain.global_internal_force[global.subdomain.support[global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face]] + i]][j]
                        += eta / he * subdomain_internal_force[option.dim * (i + 1) + j] * jacobian * w[s] * w[t];   
                    }
                }
                for(int i = 0 ; i < option.dim; i++)
                    global.subdomain.global_internal_force[global.subdomain.pair_point_ib[2 * face]][i]
                        +=  eta / he * subdomain_internal_force[i] * jacobian * w[s] * w[t];

                for(int i = 0; i < option.dim * (N2_support + 1); i++){
                    double subdomain_internal_force_i = 0.;
                    for(int j = 0; j < option.dim; j++){
                        subdomain_internal_force_i += -N2T[i][j] * (u1[j] - u2[j]);
                    }
                    subdomain_internal_force[i] = subdomain_internal_force_i;
                }
                for(int i = 0; i < N2_support; i++){
                    for(int j = 0; j < option.dim; j++){
                        global.subdomain.global_internal_force[global.subdomain.support[global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + 1]] + i]][j]
                        += eta / he * subdomain_internal_force[option.dim * (i + 1) + j] * jacobian * w[s] * w[t];   
                    }
                }
                for(int i = 0 ; i < option.dim; i++)
                    global.subdomain.global_internal_force[global.subdomain.pair_point_ib[2 * face + 1]][i]
                        +=  eta / he * subdomain_internal_force[i] * jacobian * w[s] * w[t];
            }
        }
    }

    free(node_XYZ);
    free(point_XYZ);
}

double calc_global_force_residual_norm(int iteration_step){
    double global_f_norm = 0., global_r_norm = 0.;
    char FILE_name[128];
    FILE *fp_debug;

    for(int i = 0; i < global.subdomain.N_point; i++)
        for(int j = 0; j < option.dim; j++)
            global.subdomain.global_residual_force[option.dim * i + j]
                 =  global.subdomain.global_external_force[i][j] - global.subdomain.global_internal_force[i][j];
    
    
            ImposeDirichretResidual(iteration_step);
    
    #if 0
            snprintf(FILE_name, 128,"debug_for_residual/residual_vector%d.dat", iteration_step);
            fp_debug = fopen(FILE_name,"w");
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+15.14e  ", global.subdomain.global_residual_force[option.dim * i + j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
    #endif
    
    //ノルムの計算
    for(int i = 0; i < global.subdomain.N_point; i++)
        for(int j = 0; j < option.dim; j++){
            global_f_norm += global.subdomain.global_external_force[i][j] * global.subdomain.global_external_force[i][j];
            global_r_norm += global.subdomain.global_residual_force[option.dim * i + j] * global.subdomain.global_residual_force[option.dim * i + j];
        }
    global_f_norm = sqrt(global_f_norm);
    global_r_norm = sqrt(global_r_norm);
 
    if(iteration_step == 0){
        global.temp = global_r_norm;
        //printf("%+15.14e\n", global.temp);
    }

    if(global_f_norm == 0){
        //printf("norm: %+15.14e\n", global.temp-0.25);
        return global_r_norm / global.temp;
        
    }else{
        //printf("status\n");
        return global_r_norm / global_f_norm;
    }
}
double calc_du_norm(){
    
}
void update_nodal_displacement_increment(double *latest_point_xyz){
    double *latest_node_xyz;
    double u_h[3];
    #if 0
    for(int node = 0; node < global.subdomain.N_node; node++){
        int N_ar_point = global.subdomain.ar_node_offset[node + 1] - global.subdomain.ar_node_offset[node];
        for(int i = 0; i < option.dim; i++){
            for(int j = 0; j < N_ar_point; j++){
                global.subdomain.nodal_displacement_increments[node][i]
                    += global.subdomain.displacement_increment[global.subdomain.ar_node[global.subdomain.ar_node_offset[node] + j]][i] / (double)N_ar_point;
            }
        }
    }
    #endif
    #if 1
    if((latest_node_xyz = (double *)calloc(option.dim * global.subdomain.N_node, sizeof(double))) == NULL){
        printf("latest_point_xyz's is not neough\n");
        exit(-1);
    }
    for(int i = 0; i < global.subdomain.N_node; i++){
        for(int j = 0; j < option.dim; j++){
            latest_node_xyz[option.dim * i + j] = global.subdomain.node_XYZ[option.dim * i + j]
                                                + global.subdomain.nodal_displacements[i][j]
                                                + global.subdomain.nodal_displacement_increments[i][j];
        }
    }
    for(int node = 0; node < global.subdomain.N_node; node++){
        int N_ar_point = global.subdomain.ar_node_offset[node + 1] - global.subdomain.ar_node_offset[node];
        for(int i = 0; i < N_ar_point; i++){
            trial_u(latest_node_xyz, global.subdomain.ar_node[global.subdomain.ar_node_offset[node] + i], latest_point_xyz, u_h);
            for(int j = 0; j < option.dim; j++){
                global.subdomain.nodal_displacement_increments[node][j] += u_h[j];
            }
        }
        for(int i = 0; i < option.dim; i++)
            global.subdomain.nodal_displacement_increments[node][i] /= (double)N_ar_point;
    }
    free(latest_node_xyz);
    #endif

}

void increment_field(){
    double stresses[6];
    //変位を更新し、変位増分をゼロ処理
    for(int point = 0;  point < global.subdomain.N_point; point++){
        for(int i = 0; i < option.dim; i++){
            global.subdomain.displacement[point][i] += global.subdomain.displacement_increment[point][i];
            global.subdomain.displacement_increment[point][i] = 0.;
        }
    }

    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i = 0; i < option.dim; i++){
            for(int j = 0; j < option.dim; j++){
                global.subdomain.deformation_gradients[i][j][point] = global.subdomain.current_deformation_gradients[i][j][point];
            }
        }

        for(int i = 0; i < 6; i++){
            global.subdomain.elastic_strains[point][i] = global.subdomain.current_elastic_strains[point][i];
            global.subdomain.stresses[point][i] = global.subdomain.current_stresses[point][i];
            global.subdomain.back_stresses[point][i] = global.subdomain.current_back_stresses[point][i];
            stresses[i] = global.subdomain.stresses[point][i];  //相当応力の計算のためにポイントごとの応力ベクトルを用意
        }

        global.subdomain.equivalent_stresses[point] = calc_equivalent_stress(stresses);

        global.subdomain.equivalent_plastic_strains[point] += global.subdomain.equivalent_plastic_strain_increments[point];

        global.subdomain.equivalent_plastic_strain_increments[point] = 0.;

        global.subdomain.yield_stresses[point] = global.subdomain.current_yield_stresses[point];
    }

    for(int node = 0; node < global.subdomain.N_node; node++){
        for(int i = 0; i < option.dim; i++){
            global.subdomain.nodal_displacements[node][i] += global.subdomain.nodal_displacement_increments[node][i];
            global.subdomain.nodal_displacement_increments[node][i] = 0.;
        }
    }
}