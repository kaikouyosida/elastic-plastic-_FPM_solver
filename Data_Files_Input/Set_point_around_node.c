#include<stdio.h>
#include<stdlib.h>

int main(){
  int N_point = 0;  // ポイント(サブドメイン)の数

  int N_node = 0;        // ノードの数
  int sum_N_vertex = 0;  // 各サブドメインがもつ頂点の数の総和
  // 各要素の頂点が記録されたnode配列の先頭番号 //
  // 最後の要素にはsum_N_vertexを格納 //
  int *vertex_offset;
	int *node;             // 各サブドメインの頂点の番号

  int *subdomain;         // 各ノードが属するサブドメインの番号
  int *N_subdomain;       // 各ノードが属するサブドメインの数
  int *subdomain_offset;  // subdomain配列のオフセット



  // 各サブドメインの頂点番号を読み込む (処理1) //
  int buf = 0;         // 整数値を格納する変数
  int count = 0;       // カウンタ
  FILE *fp_subdomain;  // サブドメインのデータ

  // 各頂点が属するサブドメインの番号を記録して出力 (処理2) //
  int node_num = 0;        // 頂点の番号
  FILE *fp_point_ar_node;  // 頂点周りのポイントの番号のファイル



  // 各サブドメインの頂点番号を読み込む (処理1) //
  fp_subdomain = fopen("Input_subdomain.dat","r");
  if(fp_subdomain == NULL){
    printf("file not open\n");
    return -1;
  }

  fscanf(fp_subdomain,"%*s %d\n", &N_point);
  fscanf(fp_subdomain,"%*s %d\n", &sum_N_vertex);
  fscanf(fp_subdomain,"%*s %d\n", &N_node);
  fclose(fp_subdomain);

  if((vertex_offset = (int *)calloc(N_point+1,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((node = (int *)calloc(sum_N_vertex,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((subdomain = (int *)calloc(sum_N_vertex,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  for(int i=0;i<sum_N_vertex;i++) subdomain[i] = -1;
  if((subdomain_offset = (int *)calloc(N_node,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((N_subdomain = (int *)calloc(N_node,sizeof(int))) == NULL){
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
      N_subdomain[node[count]] += 1;

			count += 1;
		}
		fgetc(fp_subdomain);
	}

  fclose(fp_subdomain);
  count = 0;

  // subdomain_offset配列を作る //
  subdomain_offset[0] = 0;
  for(int i=1;i<N_node;i++){
    subdomain_offset[i] = subdomain_offset[i-1] + N_subdomain[i-1];
  }



  // 各頂点が属するサブドメインの番号を調べる (処理2) //
  for(int i=0;i<N_point;i++){
    for(int j=vertex_offset[i];j<vertex_offset[i+1];j++){
      node_num = node[j];

      for(int k=0;k<N_subdomain[node_num];k++){
        if(subdomain[subdomain_offset[node_num]+k] == -1){
          subdomain[subdomain_offset[node_num]+k] = i;
          break;
        }
      }
    }
  }

  fp_point_ar_node = fopen("Input_point_around_node.dat","w");
  if(fp_point_ar_node == NULL){
    printf("file not open\n");
    return -1;
  }

  fprintf(fp_point_ar_node,"N_node  %d\n", N_node);
  fprintf(fp_point_ar_node,"sum_N_point  %d\n", sum_N_vertex);
  fprintf(fp_point_ar_node,"node number / number of points / point number\n");
  for(int i=0;i<N_node;i++){
    fprintf(fp_point_ar_node,"%7d  %2d ", i, N_subdomain[i]);
    for(int j=0;j<N_subdomain[i];j++) fprintf(fp_point_ar_node," %7d", subdomain[subdomain_offset[i]+j]);
    fprintf(fp_point_ar_node,"\n");
  }

  fclose(fp_point_ar_node);
  free(N_subdomain);
  free(subdomain_offset);
  free(subdomain);
  free(node);
  free(vertex_offset);

  return 0;
}
