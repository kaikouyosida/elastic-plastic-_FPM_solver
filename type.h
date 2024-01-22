typedef struct 
{
    int dim;                //解析次元
    double magni;           //ペナルティパラメータとヤング率の比
    double Delta_time;      //時間刻み幅（Newton-Raphson法）
    int N_timestep;         //タイムステップ数
    int time_output;        //データを出力する時間刻み幅
    double NR_tol;          //Newton-Raphson法の収束閾値
    int stop_count;         //Newton_raphson法の収束が著しく悪い時に反復をストップする回数
    double time;            //Newton-Raphson法の仮想的な時間
    double time_old;        //前ステップの仮想時間
    double time_end;        //反復を終了する仮想的な時間
    double norm_comp;       //サブドメインごとの残差ノルム
    double R_max;           //norm_compの最大値
    double R_first;         //R_maxの最大値
    double Error_NR;        //ニュートンラプソン法のエラー
    double norm_Fext;       //外力ベクトルのノルム

} Option;

typedef struct
{
    int material_type;                       //材料モデルも種類を選ぶフラグ
    double rho;                              //材料の密度
    double E_mod;                            //ヤング率
    double nu_mod;                           //ポアソン比
    double penalty;                          //ペナルティパラメータ
    double kinematic_hardening_fractions;    //移動硬化の影響を示すパラメータ
} Material;

typedef struct
{
    int N_point;                        //ポイント数
    int N_node;                         //サブドメインの頂点の数
    double *node_XYZ;                   //頂点の座標
    double *point_XYZ;                  //ポイントの座標
    double *center_XYZ;                 //サブドメインの中心までの距離（積分計算に使う）
    int sum_N_face;                     //面のアドレス配列の要素合計
    int N_face;                         //面の数
    int sum_N_vertex;                   //頂点のアドレス配列の要素合計
    int *face_offset;                   //面の先頭番号の配列
    int *face;                          //面の番号
    int *vertex_offset;                 //頂点の先頭番号の配列
    int *node;                          //頂点の番号
    int *pair_point_ib;                 //隣接するポイント番号
    int *shared_face;                   //共有する内部境界の面番号
    int N_int_boundary;                 //内部境界の数
    int *support;                       //サポートドメインに含まれるポイント番号
    int *support_offset;                //サポート番号の先頭番号の配列
    int sum_N_support;                  //サポート番号のアドレス配列の要素合計
    int sum_N_node;                     //サブドメインの頂点番号の配列の要素尾久系
    int *subdomain_node;                //サブドメインの頂点番号の配列
    int *ar_node_offset;                //サブドメイン内に配置された頂点の先頭番号
    int *ar_node;                       //サブドメイン内に配置された頂点番号

    double **displacement;              //ポイントごとの変位
    double **displacement_increment;    //ポイントごとの変位の増分

    double **external_force;                    //外力ベクトル
    double **global_external_force;             //全体の外力ベクトル
    double **previous_global_external_force;    //前ステップの外力ベクトル

    double **internal_force;                    //内力ベクトル
    double **global_internal_force;             //全体の内力ベクトル
    
    double ***deformation_gradients;            //変形勾配テンソル
    double ***current_deformation_gradients;    //現配置に対する変形勾配テンソル

    double **elastic_strains;                   //弾性ひずみ
    double **current_elastic_strains;           //現配置に対する弾性ひずみ
    double **trial_elastic_strains;             //試行弾性ひずみ

    double **stresses;                          //応力テンソル（フォークト表記）
    double **current_stresses;                  //現配置に対する応力

    double *equivalent_stresses;                    //相当応力
    double *equivalent_plastic_strains;             //相当塑性ひずみ
    double *equivalent_plastic_strain_increments;   //相当塑性ひずみ増分（塑性乗数）

    double *yield_stresses;                         //降伏応力
    double *current_yield_stresses;                 //現配置における降伏応力

    double **back_stresses;                         //背応力
    double **current_back_stresses;                 //現配置における背応力

    double **nodal_displacements;                   //節点の変位
    double **nodal_displacement_increments;         //節点の変位増分
    double **nodal_stresses;                        //応力の節点値
    double *nodal_equivalent_plastic_strains;       //相当塑性ひずみの節点値
    double *nodal_yield_stresses;                   //降伏応力の節点値
    double **nodal_back_stresses;                    //背応力の節点値

} Subdomain;

typedef struct
{
    int *fixed_dir;                     //ディリクレ境界条件をかける方向
    int *Dirichlet_type;                //ディリクレ境界条件のタイプ
    int N_D_DoF;                        //ディリクレ境界条件の自由度数
    int buf_dir;                        //バッファ
    int buf_type;                        //バッファ
    int N_t_face;                       //トラクションをかける面の数
    int *traction_point;                //トラクションをかけるポイント番号
    int *traction_type;                 //トラクションをかけるタイプ
    int *traction_face;                 //トラクションをかける面番号
} BC;

typedef struct{
    double *plastic_strain;
    double *stress;
    int num;
} SS_CURVE;

typedef struct{
    Material material;
    Subdomain subdomain;
    BC bc;

    int buf;                                        //バッファ
    int count;                                      //カウンタ
    
} Global;
