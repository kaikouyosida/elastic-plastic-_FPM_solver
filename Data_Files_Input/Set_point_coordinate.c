#include<stdio.h>
#include<stdlib.h>

int main(){
  // 2次元ならば、サブドメインの形は三角形or四角形と約束 //
  // 3次元ならば、サブドメインの形は四面体or六面体 //

  int dim = 0;             // モデルの次元
  int type_subdomain = 0;  // サブドメインの種類

  int N_point = 0;     // ポイントの数
  double *point_XYZ;   // ポイントの初期配置の座標

  int N_face = 0;      // 面の数
  int sum_N_face = 0;  // サブドメインがもつ面の数の総和
  int face_s = 0;      // 1つのサブドメインの面の数
  int *face;           // 各サブドメインがもつ面の番号

  int N_node = 0;        // ノードの数
  int sum_N_vertex = 0;  // 各面がもつ頂点の数の総和
  int vertex_f = 0;      // 1つの面の頂点の数
  int *node;             // 各面の頂点の番号
	double *node_XYZ;      // ノードの初期配置の座標

  int sum_N_neighbor = 0;  // 隣接点の数の総和
  int *neighbor;           // 各ポイントの隣接点の番号
  // 各ポイントの隣接点が記録されたneighbor配列の先頭番号 //
  // 最後の要素にはsum_N_neighborを格納 //
  int *neighbor_offset;

  int N_int_boundary = 0;  // 内部境界の数

  // 複数の処理で使い回す変数 //
  int count = 0;  // カウンタ



  // モデルの設定ファイルを開く (処理1) //
  FILE *fp_model;  // モデルの設定ファイル

  // 頂点の座標を読み込む (処理2) //
  FILE *fp_subdomain;  // サブドメインのデータ

  // サブドメインの中心をポイントの座標として読み込む (処理3) //
  FILE *fp_center;  // サブドメイン中心のデータ

  // 面のデータを開き、サブドメインと面の関係を読み込む (処理4) //
  FILE *fp_face;   // 面のデータ

  // 各ポイントの隣接点の情報を読み込む (処理5) //
  int buf = 0;        // 整数値を格納する変数
  FILE *fp_neighbor;  // 隣接点のデータ

  // 内部境界の情報を読み込み、face配列から内部境界を消去 (処理6) //
  int num1 = 0;                // ポイントの番号
  int num2 = 0;                // 上に同じ
  int int_boundary = 0;        // 内部境界の番号
  FILE *fp_internal_boundary;  // 内部境界のデータ

  // 外部境界を有するサブドメインは、ポイントの位置を変える (処置7) //
  int *face_ex;  // 外部境界の面の番号
  int *node_ex;  // 外部境界のノードの番号

  // ポイントの座標を出力 (処置8) //
  FILE *fp_point;  // ポイントの座標を出力するファイル



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

  if(type_subdomain == 0){
    vertex_f = 2;
    face_s = 3;
  }
  else if(type_subdomain == 1){
    vertex_f = 2;
    face_s = 4;
  }
  else if(type_subdomain == 2){
    vertex_f = 3;
    face_s = 4;
  }
  else if(type_subdomain == 3){
    vertex_f = 4;
    face_s = 6;
  }



  // 頂点の座標を読み込む (処理2) //
  fp_subdomain = fopen("Input_subdomain.dat","r");
  if(fp_subdomain == NULL){
    printf("file not open\n");
    return -1;
  }

  fscanf(fp_subdomain,"%*s %d\n", &N_point);
  fscanf(fp_subdomain,"%*[^\n]\n");
  fscanf(fp_subdomain,"%*s %d\n", &N_node);
  fclose(fp_subdomain);

  if((node_XYZ = (double *)calloc(dim*N_node,sizeof(double))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}

  fp_subdomain = fopen("Input_subdomain.dat","r");
  if(fp_subdomain == NULL){
    printf("file not open\n");
    return -1;
  }
  for(int i=0;i<4;i++) fscanf(fp_subdomain,"%*[^\n]\n");
  for(int i=0;i<N_point;i++) fscanf(fp_subdomain,"%*[^\n]\n");
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



  // サブドメインの中心をポイントの座標として読み込む (処理3) //
  if((point_XYZ = (double *)calloc(dim*N_point,sizeof(double))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  fp_center = fopen("Input_center_of_subdomain.dat","r");
  if(fp_center == NULL){
    printf("file not open\n");
    return -1;
  }
  for(int i=0;i<2;i++) fscanf(fp_center,"%*[^\n]\n");

  for(int i=0;i<N_point;i++){
    fscanf(fp_center,"%*d");
    for(int j=0;j<dim;j++){
      fscanf(fp_center,"%lf", &point_XYZ[dim*i+j]);
    }
    fgetc(fp_center);
  }

  fclose(fp_center);



  // 面のデータを開き、サブドメインと面の関係を読み込む (処理4) //
  fp_face = fopen("Input_face_of_subdomain.dat","r");
  if(fp_face == NULL){
    printf("file not open\n");
    return -1;
  }

  fscanf(fp_face,"%*[^\n]\n");
  fscanf(fp_face,"%*s %d\n", &sum_N_face);
  fscanf(fp_face,"%*[^\n]\n");
  for(int i=0;i<N_point;i++) fscanf(fp_face,"%*[^\n]\n");
  fgetc(fp_face);
  fscanf(fp_face,"%*s %d\n", &N_face);
  fscanf(fp_face,"%*s %d\n", &sum_N_vertex);
  fclose(fp_face);

	if((face = (int *)calloc(sum_N_face,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((node = (int *)calloc(sum_N_vertex,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}

  fp_face = fopen("Input_face_of_subdomain.dat","r");
  if(fp_face == NULL){
    printf("file not open\n");
    return -1;
  }
  for(int i=0;i<3;i++) fscanf(fp_face,"%*[^\n]\n");

  for(int i=0;i<N_point;i++){
    fscanf(fp_face,"%*d %*d");
    for(int j=0;j<face_s;j++){
      fscanf(fp_face,"%d", &face[count]);
      count += 1;
    }
    fgetc(fp_face);
  }
  count = 0;

  fgetc(fp_face);
  for(int i=0;i<3;i++) fscanf(fp_face,"%*[^\n]\n");

  for(int i=0;i<N_face;i++){
    fscanf(fp_face,"%*d %*d");
    for(int j=0;j<vertex_f;j++){
      fscanf(fp_face,"%d", &node[count]);
      count += 1;
    }
    fgetc(fp_face);
  }

  fclose(fp_face);
  count = 0;



  // 各ポイントの隣接点の情報を読み込む (処理5) //
  fp_neighbor = fopen("Input_neighboring_point.dat","r");
  if(fp_neighbor == NULL){
    printf("file not open\n");
    return -1;
  }

  fscanf(fp_neighbor,"%*[^\n]\n");
  fscanf(fp_neighbor,"%*s %d\n", &sum_N_neighbor);
  fclose(fp_neighbor);

	if((neighbor_offset = (int *)calloc(N_point+1,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
	if((neighbor = (int *)calloc(sum_N_neighbor,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}

  fp_neighbor = fopen("Input_neighboring_point.dat","r");
  if(fp_neighbor == NULL){
    printf("file not open\n");
    return -1;
  }
  for(int i=0;i<3;i++) fscanf(fp_neighbor,"%*[^\n]\n");

  neighbor_offset[0] = 0;
  for(int i=0;i<N_point;i++){
    fscanf(fp_neighbor,"%*d %d", &buf);
    neighbor_offset[i+1] = neighbor_offset[i] + buf;

    for(int j=0;j<buf;j++){
      fscanf(fp_neighbor,"%d", &neighbor[count]);
      count += 1;
    }
    fgetc(fp_neighbor);
  }

  fclose(fp_neighbor);
  count = 0;



  // 内部境界の情報を読み込み、face配列から内部境界を消去 (処理6) //
  fp_internal_boundary = fopen("Input_internal_boundary.dat", "r");
  if(fp_internal_boundary == NULL){
    printf("file is not open\n");
    return -1;
  }
  fscanf(fp_internal_boundary,"%*s %d\n", &N_int_boundary);
  fscanf(fp_internal_boundary,"%*[^\n]\n");
  for(int i=0;i<N_int_boundary;i++){
    fscanf(fp_internal_boundary,"%d %d %d\n", &num1, &num2, &int_boundary);

    for(int j=face_s*num1;j<face_s*(num1+1);j++){
      // 内部境界に対応するface配列の要素を-1にする //
      if(face[j] == int_boundary){
        face[j] = -1;
        break;
      }
    }
    for(int j=face_s*num2;j<face_s*(num2+1);j++){
      // 内部境界に対応するface配列の要素を-1にする //
      if(face[j] == int_boundary){
        face[j] = -1;
        break;
      }
    }
  }

  fclose(fp_internal_boundary);



  // 外部境界を有するサブドメインは、ポイントの位置を変える (処置7) //
  if(dim == 2){
    for(int i=0;i<N_point;i++){
      if(face_s - (neighbor_offset[i+1] - neighbor_offset[i]) == 1){
        for(int j=face_s*i;j<face_s*(i+1);j++){
          if(face[j] != -1){
            for(int k=0;k<dim;k++) point_XYZ[dim*i+k] = 0.;

            for(int k=vertex_f*face[j];k<vertex_f*(face[j]+1);k++){
              for(int l=0;l<dim;l++) point_XYZ[dim*i+l] += node_XYZ[dim*node[k]+l];
            }

            for(int k=0;k<dim;k++) point_XYZ[dim*i+k] /= (double)vertex_f;
          }
        }
      }
      // サブドメインが角にあるとき //
      else if(face_s - (neighbor_offset[i+1] - neighbor_offset[i]) == 2){
        if((face_ex = (int *)calloc(2,sizeof(int))) == NULL){
      		printf("Error:Memory is not enough\n");
      		return -1;
      	}
        if((node_ex = (int *)calloc(1,sizeof(int))) == NULL){
      		printf("Error:Memory is not enough\n");
      		return -1;
      	}

        for(int j=face_s*i;j<face_s*(i+1);j++){
          if(face[j] != -1){
            face_ex[count] = face[j];
            count += 1;
          }
        }
        count = 0;

        for(int j=vertex_f*face_ex[0];j<vertex_f*(face_ex[0]+1);j++){
          for(int k=vertex_f*face_ex[1];k<vertex_f*(face_ex[1]+1);k++){
            if(node[j] == node[k]) node_ex[0] = node[j];
          }
        }

        for(int j=0;j<dim;j++) point_XYZ[dim*i+j] = node_XYZ[dim*node_ex[0]+j];

        free(node_ex);
        free(face_ex);
      }
    }
  }
  else if(dim == 3){
    for(int i=0;i<N_point;i++){
      if(face_s - (neighbor_offset[i+1] - neighbor_offset[i]) == 1){
        for(int j=face_s*i;j<face_s*(i+1);j++){
          if(face[j] != -1){
            for(int k=0;k<dim;k++) point_XYZ[dim*i+k] = 0.;

            for(int k=vertex_f*face[j];k<vertex_f*(face[j]+1);k++){
              for(int l=0;l<dim;l++) point_XYZ[dim*i+l] += node_XYZ[dim*node[k]+l];
            }

            for(int k=0;k<dim;k++) point_XYZ[dim*i+k] /= (double)vertex_f;
          }
        }
      }
      // サブドメインが辺上にあるとき //
      else if(face_s - (neighbor_offset[i+1] - neighbor_offset[i]) == 2){
        if((face_ex = (int *)calloc(2,sizeof(int))) == NULL){
      		printf("Error:Memory is not enough\n");
      		return -1;
      	}
        if((node_ex = (int *)calloc(2,sizeof(int))) == NULL){
      		printf("Error:Memory is not enough\n");
      		return -1;
      	}

        for(int j=face_s*i;j<face_s*(i+1);j++){
          if(face[j] != -1){
            face_ex[count] = face[j];
            count += 1;
          }
        }
        count = 0;

        for(int j=vertex_f*face_ex[0];j<vertex_f*(face_ex[0]+1);j++){
          for(int k=vertex_f*face_ex[1];k<vertex_f*(face_ex[1]+1);k++){
            if(node[j] == node[k]){
              node_ex[count] = node[j];
              count += 1;
            }
          }
        }
        count = 0;

        for(int j=0;j<dim;j++) point_XYZ[dim*i+j] = 0.5*(node_XYZ[dim*node_ex[0]+j] + node_XYZ[dim*node_ex[1]+j]);

        free(node_ex);
        free(face_ex);
      }
      // サブドメインが角にあるとき //
      else if(face_s - (neighbor_offset[i+1] - neighbor_offset[i]) == 3){
        if((face_ex = (int *)calloc(3,sizeof(int))) == NULL){
      		printf("Error:Memory is not enough\n");
      		return -1;
      	}
        if((node_ex = (int *)calloc(2,sizeof(int))) == NULL){
      		printf("Error:Memory is not enough\n");
      		return -1;
      	}

        for(int j=face_s*i;j<face_s*(i+1);j++){
          if(face[j] != -1){
            face_ex[count] = face[j];
            count += 1;
          }
        }
        count = 0;

        for(int j=vertex_f*face_ex[0];j<vertex_f*(face_ex[0]+1);j++){
          for(int k=vertex_f*face_ex[1];k<vertex_f*(face_ex[1]+1);k++){
            if(node[j] == node[k]){
              node_ex[count] = node[j];
              count += 1;
            }
          }
        }
        count = 0;

        for(int j=vertex_f*face_ex[2];j<vertex_f*(face_ex[2]+1);j++){
          for(int k=0;k<2;k++){
            if(node[j] == node_ex[k]){
              for(int l=0;l<dim;l++) point_XYZ[dim*i+l] = node_XYZ[dim*node[j]+l];
            }
          }
        }

        free(node_ex);
        free(face_ex);
      }
    }
  }



  // ポイントの座標を出力 (処置8) //
  fp_point = fopen("Input_point_coordinate.dat","w");
  if(fp_point == NULL){
    printf("file not open\n");
    return -1;
  }

  fprintf(fp_point,"N_point  %d\n", N_point);
  if(dim == 2) fprintf(fp_point,"point number / coordinate  X Y\n");
  else if(dim == 3) fprintf(fp_point,"point number / coordinate  X Y Z\n");
  for(int i=0;i<N_point;i++){
    fprintf(fp_point,"%7d", i);
    for(int j=0;j<dim;j++) fprintf(fp_point,"  %+.15e", point_XYZ[dim*i+j]);
    fprintf(fp_point,"\n");
  }

  fclose(fp_point);
  free(neighbor);
  free(neighbor_offset);
  free(node);
  free(face);
  free(point_XYZ);
  free(node_XYZ);
  return 0;
}
