//ベクトルを降順にする
void sortVector(double *vector, const int num);

//ベクトルの成分をreverseする
void reverseVector(double *vector, const int num);

//ベクトルのアセンブル
void assemble_vector(int point_n, double **global_vecter, double *element_vector);

//頂点2から頂点1に伸びるベクトルの生成
void generate_current_node_vector(int node_1, int node_2, double *vector);

//頂点から任意の座標に伸びるベクトルの生成
void generate_current_node_to_point_vector(int node,double *point,double *vector);