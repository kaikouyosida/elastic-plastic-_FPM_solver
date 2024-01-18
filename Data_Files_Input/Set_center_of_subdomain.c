#include<stdio.h>
#include<stdlib.h>

int main(){
  int dim = 0;  // モデルの次元

  int N_point = 0;     // ポイントの数
  double *center_XYZ;  // サブドメインの中心の座標

  int N_node = 0;        // ノードの数
  int sum_N_vertex = 0;  // 各サブドメインがもつ頂点の数の総和
  // 各要素の頂点が記録されたnode配列の先頭番号 //
  // 最後の要素にはsum_N_vertexを格納 //
  int *vertex_offset;
  int *node;             // 各サブドメインの頂点の番号
	double *node_XYZ;      // ノードの初期配置の座標



  // モデルの設定ファイルを開く (処理1) //
  FILE *fp_model;  // モデルの設定ファイル

  // 各サブドメインの頂点番号と、頂点の座標を読み込む (処理2) //
  int buf = 0;         // 整数値を格納する変数
  int count = 0;       // カウンタ
  FILE *fp_subdomain;  // サブドメインのデータ

  // 各サブドメインの中心座標(頂点座標の算術平均)を計算 (処理3) //

  // 各サブドメインの中心座標を出力 (処理4) //
  FILE *fp_center;  // サブドメイン中心のデータ



  // モデルの設定ファイルを開く (処理1) //
  fp_model = fopen("Model_settings.dat","r");
  if(fp_model == NULL){
    printf("file not open\n");
    return -1;
  }
  fscanf(fp_model,"%*[^\n]\n");
  fscanf(fp_model,"%d\n", &dim);
  fclose(fp_model);



  // 各サブドメインの頂点番号と、頂点の座標を読み込む (処理2) //
  fp_subdomain = fopen("Input_subdomain.dat","r");
  if(fp_subdomain == NULL){
    printf("file not open\n");
    return -1;
  }

  fscanf(fp_subdomain,"%*s %d\n", &N_point);
  fscanf(fp_subdomain,"%*s %d\n", &sum_N_vertex);
  fscanf(fp_subdomain,"%*s %d\n", &N_node);
  fclose(fp_subdomain);

  if((center_XYZ = (double *)calloc(dim*N_point,sizeof(double))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((node_XYZ = (double *)calloc(dim*N_node,sizeof(double))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((vertex_offset = (int *)calloc(N_point+1,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((node = (int *)calloc(sum_N_vertex,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}

  fp_subdomain = fopen("Input_subdomain.dat","r");
  if(fp_subdomain == NULL){
    printf("file not open\n");
    return -1;
  }
  for(int i=0;i<4;i++) fscanf(fp_subdomain,"%*[^\n]\n");

  vertex_offset[0] = 0;
  for(int i=0;i<N_point;i++){
    fscanf(fp_subdomain,"%*d %d", &buf);
    vertex_offset[i+1] = vertex_offset[i] + buf;

    for(int j=0;j<buf;j++){
      fscanf(fp_subdomain,"%d", &node[count]);
      count += 1;
    }
    fgetc(fp_subdomain);
  }

  fgetc(fp_subdomain);
  fscanf(fp_subdomain,"%*[^\n]\n");

  for(int i=0;i<N_node;i++){
    fscanf(fp_subdomain,"%*d");
    for(int j=0;j<dim;j++){
      fscanf(fp_subdomain,"%lf", &node_XYZ[dim*i+j]);
    }
    fgetc(fp_subdomain);
  }
  fclose(fp_subdomain);



  // 各サブドメインの中心座標(頂点座標の算術平均)を計算 (処理3) //
  for(int i=0;i<N_point;i++){
    for(int j=0;j<dim;j++){
      for(int k=vertex_offset[i];k<vertex_offset[i+1];k++) center_XYZ[dim*i+j] += node_XYZ[dim*node[k]+j];
      center_XYZ[dim*i+j] /= (double)(vertex_offset[i+1] - vertex_offset[i]);
    }
  }



  // 各サブドメインの中心座標を出力 (処理4) //
  fp_center = fopen("Input_center_of_subdomain.dat","w");
  if(fp_center == NULL){
    printf("file not open\n");
    return -1;
  }

  fprintf(fp_center,"N_point  %d\n", N_point);
  if(dim == 2) fprintf(fp_center,"subdomain number / coordinates of center  X Y\n");
  else if(dim == 3) fprintf(fp_center,"subdomain number / coordinates of center  X Y Z\n");
  for(int i=0;i<N_point;i++){
    fprintf(fp_center,"%7d", i);
    for(int j=0;j<dim;j++) fprintf(fp_center,"  %+.15e", center_XYZ[dim*i+j]);
    fprintf(fp_center,"\n");
  }

  fclose(fp_center);
  free(node);
  free(vertex_offset);
  free(node_XYZ);
  free(center_XYZ);

  return 0;
}
