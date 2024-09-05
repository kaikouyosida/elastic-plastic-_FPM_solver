//ベクトルを降順にする
void sortVector(double *vector, const int num);

//ベクトルの成分をreverseする
void reverseVector(double *vector, const int num);

//ベクトルのアセンブル
void assemble_vector(int point_n, double **global_vecter, const double *element_vector);

//頂点2から頂点1に伸びるベクトルの生成
void generate_current_node_vector(int node_1, int node_2, double *vector);
void generate_current_edge_vector(double *vector ,int point_n, int node_id_1, int node_id_2, int *subdomain_node);

//頂点から任意の座標に伸びるベクトルの生成
void generate_current_node_to_point_vector(int node,double *point,double *vector);
void generate_current_points_vector(double *vector, double *point_xyz, int point_n, int node_id,  int *subdomain_node);

//point側のサブドメインにおける形状関数から得た節点の現在座標
void generate_current_face_node(double face_node_XYZ[4][3], int *node_id ,int subdomain_node[8], int point);
void generate_current_node_of_face(double face_node_XYZ[4][3], int face_n, int point);

//サブドメインがもつ節点の番号を格納する
void generate_subdomain_node(int point_n, int subdomain_node[8]);

//内部境界面のノード番号のアドレスを格納
void generate_node_id(int face_n, int point_n, int subdomain_node[8], int node_id[4]);

//物理座標におけるガウス点の座標を計算
void generate_gauss_point_coordinate(const int s, const int t, const double face_node_XYZ[4][3], const double *X, double xyz[3]);
