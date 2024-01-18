#include<stdio.h>
#include<stdlib.h>

void sort_NodeNum(int vertex_f, int* ex_num, int* buf, int* node_sort);
void set_number_tetra(int j, int* ex_num);
void set_number_hexa(int j, int* ex_num);

int main(){
  // サブドメインの形は凸多角形または凸面体と約束する //
  // 2次元の場合、頂点番号が時計/反時計回りに記述されてる必要あり //
  // 四面体または六面体の場合、頂点番号が決まった順で記述されている必要あり //

  int dim = 0;  // モデルの次元

  int N_subdomain = 0;     // サブドメインの数
  int type_subdomain = 0;  // サブドメインの種類
  int vertex_f = 0;        // 1つの面の頂点の数
  int face = 0;            // 1つのサブドメインの面の数
  int sum_N_vertex = 0;    // 頂点数の総和



  // モデルの設定ファイルを開く (処理1) //
  FILE *fp_model;  // モデルの設定ファイル

  // 各サブドメインの頂点番号を読み取り、面の情報を出力 (処理2) //
  int size_node = 0;   // node配列のサイズ
  int N_vertex = 0;    // サブドメインの頂点数
  int N_face = 0;      // 面の数
  int flag = 0;        // 新しい面か否かを判定するフラグ
  int *buf;            // 頂点番号を一時的に格納する配列
  int *ex_num;         // buf配列から取り出す要素の番号
  int *node_sort;      // 頂点番号をソートした配列
  int *node;           // 各面がもつ頂点の番号
  FILE *fp_subdomain;  // サブドメインのデータ
  FILE *fp_face;       // 面の情報を出力するファイル



  // モデルの設定ファイルを開く (処理1) //
  fp_model = fopen("Model_settings.dat","r");
  if(fp_model == NULL){
    printf("file not open\n");
    return -1;
  }
  fscanf(fp_model,"%*[^\n]\n");
  fscanf(fp_model,"%d\n", &dim);
  fgetc(fp_model);
  fscanf(fp_model,"%*[^\n]\n");
  fscanf(fp_model,"%*[^\n]\n");
  fscanf(fp_model,"%d\n", &type_subdomain);

  fclose(fp_model);

  if(dim == 2) vertex_f = 2;
  else if(type_subdomain == 2){
    vertex_f = 3;
    face = 4;
  }
  else if(type_subdomain == 3){
    vertex_f = 4;
    face = 6;
  }



  // 各サブドメインの頂点番号を読み取り、面の情報を出力 (処理2) //
  fp_subdomain = fopen("Input_subdomain.dat","r");
  if(fp_subdomain == NULL){
    printf("file not open\n");
    return -1;
  }

  fscanf(fp_subdomain,"%*s %d\n", &N_subdomain);
  fscanf(fp_subdomain,"%*s %d\n", &sum_N_vertex);
  fclose(fp_subdomain);

  if(dim == 2) size_node = 2*sum_N_vertex;
  else if(dim == 3) size_node = N_subdomain*vertex_f*face;

  if((node = (int *)calloc(size_node,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((node_sort = (int *)calloc(vertex_f,sizeof(int))) == NULL){
    printf("Error:Memory is not enough\n");
    return -1;
  }
  if((ex_num = (int *)calloc(vertex_f,sizeof(int))) == NULL){
    printf("Error:Memory is not enough\n");
    return -1;
  }

  fp_subdomain = fopen("Input_subdomain.dat","r");
  if(fp_subdomain == NULL){
    printf("file not open\n");
    return -1;
  }
  fp_face = fopen("Input_face_of_subdomain.dat","w");
  if(fp_face == NULL){
    printf("file not open\n");
    return -1;
  }

  // 各サブドメインがもつ面の番号を書き出す //
  fprintf(fp_face,"N_subdomain  %d\n", N_subdomain);
  if(dim == 2) fprintf(fp_face,"sum_N_face  %d\n", sum_N_vertex);
  else if(dim == 3) fprintf(fp_face,"sum_N_face  %d\n", N_subdomain*face);

  for(int i=0;i<4;i++) fscanf(fp_subdomain,"%*[^\n]\n");
  fprintf(fp_face,"subdomain number / face / face number\n");

  if(dim == 2){
    for(int i=0;i<N_subdomain;i++){
      fscanf(fp_subdomain,"%*d %d",&N_vertex);
      face = N_vertex;
      fprintf(fp_face,"%7d  %2d ", i, face);

      if((buf = (int *)calloc(N_vertex,sizeof(int))) == NULL){
    		printf("Error:Memory is not enough\n");
    		return -1;
      }
      for(int j=0;j<N_vertex;j++){
        fscanf(fp_subdomain,"%d", &buf[j]);
      }
      fgetc(fp_subdomain);

      // buf配列から取り出す要素の番号 //
      ex_num[0] = 0; ex_num[1] = N_vertex-1;

      sort_NodeNum(2, ex_num, buf, node_sort);

      // 既に記録している辺と比較 //
      for(int k=0;k<N_face;k++){
        if(node[2*k] == node_sort[0]){
          if(node[2*k+1] == node_sort[1]){
            fprintf(fp_face," %7d", k);
            flag = 1;
            break;
          }
        }
      }
      // 新しい辺ならば頂点番号を記録 //
      if(flag == 0){
        fprintf(fp_face," %7d", N_face);
        for(int k=0;k<2;k++) node[2*N_face+k] = node_sort[k];
        N_face += 1;
      }

      flag = 0;

      for(int j=0;j<N_vertex-1;j++){
        // buf配列から取り出す要素の番号 //
        ex_num[0] = j; ex_num[1] = j+1;

        sort_NodeNum(2, ex_num, buf, node_sort);

        // 既に記録している辺と比較 //
        for(int k=0;k<N_face;k++){
          if(node[2*k] == node_sort[0]){
            if(node[2*k+1] == node_sort[1]){
              fprintf(fp_face," %7d", k);
              flag = 1;
              break;
            }
          }
        }
        // 新しい辺ならば頂点番号を記録 //
        if(flag == 0){
          fprintf(fp_face," %7d", N_face);
          for(int k=0;k<2;k++) node[2*N_face+k] = node_sort[k];
          N_face += 1;
        }

        flag = 0;
      }

      fprintf(fp_face,"\n");
      free(buf);
    }
  }
  else if(type_subdomain == 2){
    for(int i=0;i<N_subdomain;i++){
      fscanf(fp_subdomain,"%*d %d", &N_vertex);
      fprintf(fp_face,"%7d  %2d ", i, face);
      if((buf = (int *)calloc(N_vertex,sizeof(int))) == NULL){
    		printf("Error:Memory is not enough\n");
    		return -1;
      }
      for(int j=0;j<N_vertex;j++){
        fscanf(fp_subdomain,"%d", &buf[j]);
      }
      fgetc(fp_subdomain);

      for(int j=0;j<4;j++){
        set_number_tetra(j, ex_num);

        sort_NodeNum(3, ex_num, buf, node_sort);

        // 既に記録している辺と比較 //
        for(int k=0;k<N_face;k++){
          if(node[3*k] == node_sort[0]){
            if(node[3*k+1] == node_sort[1]){
              if(node[3*k+2] == node_sort[2]){
                fprintf(fp_face," %7d", k);
                flag = 1;
                break;
              }
            }
            else if(node[3*k+2] == node_sort[1]){
              if(node[3*k+1] == node_sort[2]){
                fprintf(fp_face," %7d", k);
                flag = 1;
                break;
              }
            }
          }
        }
        // 新しい辺ならば頂点番号を記録 //
        if(flag == 0){
          fprintf(fp_face," %7d", N_face);
          for(int k=0;k<3;k++) node[3*N_face+k] = node_sort[k];
          N_face += 1;
        }

        flag = 0;
      }

      fprintf(fp_face,"\n");
      free(buf);
    }
  }
  else if(type_subdomain == 3){
    for(int i=0;i<N_subdomain;i++){
      fscanf(fp_subdomain,"%*d %d", &N_vertex);
      fprintf(fp_face,"%7d  %2d ", i, face);
      if((buf = (int *)calloc(N_vertex,sizeof(int))) == NULL){
    		printf("Error:Memory is not enough\n");
    		return -1;
      }
      for(int j=0;j<N_vertex;j++){
        fscanf(fp_subdomain,"%d", &buf[j]);
      }
      fgetc(fp_subdomain);

      for(int j=0;j<6;j++){
        set_number_hexa(j, ex_num);

        sort_NodeNum(4, ex_num, buf, node_sort);

        // 既に記録している辺と比較 //
        for(int k=0;k<N_face;k++){
          if(node[4*k] == node_sort[0]){
            if(node[4*k+1] == node_sort[1]){
              if(node[4*k+2] == node_sort[2]){
                fprintf(fp_face," %7d", k);
                flag = 1;
                break;
              }
            }
            else if(node[4*k+3] == node_sort[1]){
              if(node[4*k+2] == node_sort[2]){
                fprintf(fp_face," %7d", k);
                flag = 1;
                break;
              }
            }
          }
        }
        // 新しい辺ならば頂点番号を記録 //
        if(flag == 0){
          fprintf(fp_face," %7d", N_face);
          for(int k=0;k<4;k++) node[4*N_face+k] = node_sort[k];
          N_face += 1;
        }

        flag = 0;
      }

      fprintf(fp_face,"\n");
      free(buf);
    }
  }

  fprintf(fp_face,"\n");
  fprintf(fp_face,"N_face  %d\n", N_face);
  fprintf(fp_face,"sum_N_vertex  %d\n", N_face*vertex_f);
  fprintf(fp_face,"face number / vertex / node number\n");
  for(int i=0;i<N_face;i++){
    fprintf(fp_face,"%7d  %2d ", i, vertex_f);
    for(int j=0;j<vertex_f;j++) fprintf(fp_face," %7d", node[i*vertex_f+j]);
    fprintf(fp_face,"\n");
  }

  fclose(fp_face);
  fclose(fp_subdomain);
  free(ex_num);
  free(node_sort);
  free(node);

  return 0;
}

void sort_NodeNum(int vertex_f, int* ex_num, int* buf, int* node_sort){
  int min_i = 0;  // 頂点番号が最小の要素
  int count = 0;  // カウンタ
  int *temp;      // ソート時に使う配列

  if((temp = (int *)calloc(vertex_f,sizeof(int))) == NULL){
    printf("Error:Memory is not enough\n");
    exit(-1);
  }

  // 1つの面を構成する頂点番号を取り出す //
  for(int i=0;i<vertex_f;i++) temp[i] = buf[ex_num[i]];
  // 格納した頂点番号のうち最も小さいものを先頭にする //
  for(int i=1;i<vertex_f;i++){
    if(temp[i] < temp[min_i]) min_i = i;
  }

  for(int i=min_i;i<vertex_f;i++){
    node_sort[count] = temp[i];
    count += 1;
  }
  for(int i=0;i<min_i;i++){
    node_sort[count] = temp[i];
    count += 1;
  }

  free(temp);
}
void set_number_tetra(int j, int* ex_num){
  if(j==0){
    ex_num[0] = 0; ex_num[1] = 1; ex_num[2] = 2;
  }
  if(j==1){
    ex_num[0] = 2; ex_num[1] = 1; ex_num[2] = 3;
  }
  if(j==2){
    ex_num[0] = 2; ex_num[1] = 3; ex_num[2] = 0;
  }
  if(j==3){
    ex_num[0] = 0; ex_num[1] = 3; ex_num[2] = 1;
  }
}
void set_number_hexa(int j, int* ex_num){
  if(j==0){
    ex_num[0] = 0; ex_num[1] = 1; ex_num[2] = 2; ex_num[3] = 3;
  }
  if(j==1){
    ex_num[0] = 4; ex_num[1] = 5; ex_num[2] = 6; ex_num[3] = 7;
  }
  if(j==2){
    ex_num[0] = 0; ex_num[1] = 3; ex_num[2] = 7; ex_num[3] = 4;
  }
  if(j==3){
    ex_num[0] = 1; ex_num[1] = 2; ex_num[2] = 6; ex_num[3] = 5;
  }
  if(j==4){
    ex_num[0] = 1; ex_num[1] = 0; ex_num[2] = 4; ex_num[3] = 5;
  }
  if(j==5){
    ex_num[0] = 2; ex_num[1] = 3; ex_num[2] = 7; ex_num[3] = 6;
  }
}
