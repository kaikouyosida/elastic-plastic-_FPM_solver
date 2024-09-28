//応力節点データ等の出力
void Output_data(int time_step);

//出力用の節点変位の計算
void update_nodal_coordinate();

//抽出データの出力
void calc_extract_component();

//paraviewデータの出力
void paraview_node_data(int time_step);

//頂点を共有するサブドメイン数
int extract_node_number(int node , double (*nodal_displacement)[3]);