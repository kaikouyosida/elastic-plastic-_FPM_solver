#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){
  int dim = 0;  // モデルの次元

  int N_point = 0;  // ポイントの数

  int N_face = 0;      // 面の数
  int sum_N_face = 0;  // サブドメインがもつ面の数の総和
  int *face;           // 各サブドメインがもつ面の番号
  // 各サブドメインの面が記録されたface配列の先頭番号 //
  // 最後の要素にはsum_N_faceを格納 //
  int *face_offset;

  int N_node = 0;        // ノードの数
  int sum_N_vertex = 0;  // 各面がもつ頂点の数の総和
  int *node;             // 各面の頂点の番号
  // 各面の頂点が記録されたnode配列の先頭番号 //
  // 最後の要素にはsum_N_vertexを格納 //
  int *vertex_offset;
	double *node_XYZ;      // ノードの初期配置の座標

  int sum_N_neighbor = 0;  // 隣接点の数の総和
  int *neighbor;           // 各ポイントの隣接点の番号
  // 各ポイントの隣接点が記録されたneighbor配列の先頭番号 //
  // 最後の要素にはsum_N_neighborを格納 //
  int *neighbor_offset;

  int N_int_boundary = 0;  // 内部境界の数

  // 複数の処理で使い回す変数 //
  int buf = 0;    // 整数値を格納する変数
  int count = 0;  // カウンタ



  // モデルの設定ファイルを開く (処理1) //
  FILE *fp_model;  // モデルの設定ファイル

  // 各サブドメインの頂点座標を読み込む (処理2) //
  FILE *fp_subdomain;  // サブドメインのデータ

  // 面のデータを開き、サブドメインと面の関係を読み込む (処理3) //
  FILE *fp_face;  // 面のデータ

  // 各ポイントの隣接点の情報を読み込む (処理4) //
  FILE *fp_neighbor;  // 隣接点のデータ

  // 内部境界の情報を読み込み、face配列から内部境界を消去 (処理5) //
  int num1 = 0;                // ポイントの番号
  int num2 = 0;                // 上に同じ
  int int_boundary = 0;        // 内部境界の番号
  FILE *fp_internal_boundary;  // 内部境界のデータ

  // 境界条件を課す面のデータを読み込む (処理6) //
  int N_condition_D = 0;  // ディリクレ条件の数
  int N_condition_t = 0;  // トラクション条件の数
  int temp_int = 0;       // ソート時に使用する変数
  double temp_lf = 0.;
  int *face_dir_D;        // ディリクレ条件を課す面の方向
  double *face_XYZ_D;     // ディリクレ条件を課す面の座標
  int *fixed_dir;         // 変位を指定する成分
  int *type_D;            // ディリクレ条件の種類(詳しくは解析プログラムで決める)
  int *face_dir_t;        // トラクション条件を課す面の方向
  double *face_XYZ_t;     // トラクション条件を課す面の座標
  int *type_t;            // トラクション条件の種類(詳しくは解析プログラムで決める)
  FILE *fp_face_D;        // ディリクレ条件を課す面のデータ
  FILE *fp_face_t;        // トラクション条件を課す面のデータ

  // 境界条件を課す面を調べ、データファイルを作る (処理7) //
  int flag_D = 0;           // ディリクレ条件を課すフラグ
  int flag_dir[3];          // 変位指定された方向に関するフラグ
  int flag_t = 0;           // トラクション条件を課すフラグ
  int N_Dirichlet_DoF = 0;  // ディリクレ条件を課す自由度の数
  int N_traction_face = 0;  // トラクション条件を課す面の数
  FILE *fp_Dirichlet;       // ディリクレ条件を出力するファイル
  FILE *fp_traction;        // トラクション条件を出力するファイル



  // モデルの設定ファイルを開く (処理1) //
  fp_model = fopen("Model_settings.dat","r");
  if(fp_model == NULL){
    printf("file not open\n");
    return -1;
  }
  fscanf(fp_model,"%*[^\n]\n");
  fscanf(fp_model,"%d\n", &dim);
  fclose(fp_model);



  // 各サブドメインの頂点座標を読み込む (処理2) //
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



  // 面のデータを開き、サブドメインと面の関係を読み込む (処理3) //
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

	if((face_offset = (int *)calloc(N_point+1,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
	if((face = (int *)calloc(sum_N_face,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((vertex_offset = (int *)calloc(N_face+1,sizeof(int))) == NULL){
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

  face_offset[0] = 0;
  for(int i=0;i<N_point;i++){
    fscanf(fp_face,"%*d %d", &buf);
    for(int j=0;j<buf;j++){
      fscanf(fp_face,"%d", &face[count]);
      count += 1;
    }
    fgetc(fp_face);

    face_offset[i+1] = face_offset[i] + buf;
  }
  count = 0;
  fgetc(fp_face);
  for(int i=0;i<3;i++) fscanf(fp_face,"%*[^\n]\n");

  vertex_offset[0] = 0;
  for(int i=0;i<N_face;i++){
    fscanf(fp_face,"%*d %d", &buf);
    for(int j=0;j<buf;j++){
      fscanf(fp_face,"%d", &node[count]);
      count += 1;
    }
    fgetc(fp_face);

    vertex_offset[i+1] = vertex_offset[i] + buf;
  }

  fclose(fp_face);
  count = 0;



  // 各ポイントの隣接点の情報を読み込む (処理4) //
  fp_neighbor = fopen("Input_neighboring_point.dat","r");
  if(fp_neighbor == NULL){
    printf("file not open\n");
    return -1;
  }

  fscanf(fp_neighbor,"%*s %d\n", &N_point);
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



  // 内部境界の情報を読み込み、face配列から内部境界を消去 (処理5) //
  fp_internal_boundary = fopen("Input_internal_boundary.dat", "r");
  if(fp_internal_boundary == NULL){
    printf("file is not open\n");
    return -1;
  }
  fscanf(fp_internal_boundary,"%*s %d\n", &N_int_boundary);
  fscanf(fp_internal_boundary,"%*[^\n]\n");
  for(int i=0;i<N_int_boundary;i++){
    fscanf(fp_internal_boundary,"%d %d %d\n", &num1, &num2, &int_boundary);

    for(int j=face_offset[num1];j<face_offset[num1+1];j++){
      // 内部境界に対応するface配列の要素を-1にする //
      if(face[j] == int_boundary){
        face[j] = -1;
        break;
      }
    }
    for(int j=face_offset[num2];j<face_offset[num2+1];j++){
      // 内部境界に対応するface配列の要素を-1にする //
      if(face[j] == int_boundary){
        face[j] = -1;
        break;
      }
    }
  }

  fclose(fp_internal_boundary);



  // 境界条件を課す面のデータを読み込む (処理6) //
  fp_face_D = fopen("Boundary_condition_Dirichlet.dat","r");
  if(fp_face_D == NULL){
    printf("file not open\n");
    return -1;
  }
  fscanf(fp_face_D,"%*s %d\n", &N_condition_D);
  fclose(fp_face_D);

  if((face_dir_D = (int *)calloc(N_condition_D,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((face_XYZ_D = (double *)calloc(N_condition_D,sizeof(double))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((fixed_dir = (int *)calloc(N_condition_D,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((type_D = (int *)calloc(N_condition_D,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}

  fp_face_D = fopen("Boundary_condition_Dirichlet.dat","r");
  if(fp_face_D == NULL){
    printf("file not open\n");
    return -1;
  }
  for(int i=0;i<2;i++) fscanf(fp_face_D,"%*[^\n]\n");
  for(int i=0;i<N_condition_D;i++) fscanf(fp_face_D,"%d %lf %d %d\n", &face_dir_D[i], &face_XYZ_D[i], &fixed_dir[i], &type_D[i]);

  fclose(fp_face_D);

  // データを自由度番号の順に出力するため、fixed_dirが昇順になるようソートする //
  for(int i=0;i<N_condition_D-1;i++){
    for(int j=i+1;j<N_condition_D;j++){
      if(fixed_dir[j]<fixed_dir[i]){
        temp_int = fixed_dir[i];  fixed_dir[i] = fixed_dir[j];  fixed_dir[j] = temp_int;
        temp_int = face_dir_D[i];  face_dir_D[i] = face_dir_D[j];  face_dir_D[j] = temp_int;
        temp_lf = face_XYZ_D[i];  face_XYZ_D[i] = face_XYZ_D[j];  face_XYZ_D[j] = temp_lf;
        temp_int = type_D[i];  type_D[i] = type_D[j];  type_D[j] = temp_int;
      }
    }
  }

  fp_face_t = fopen("Boundary_condition_traction.dat","r");
  if(fp_face_t == NULL){
    printf("file not open\n");
    return -1;
  }
  fscanf(fp_face_t,"%*s %d\n", &N_condition_t);
  fclose(fp_face_t);

  if((face_dir_t = (int *)calloc(N_condition_t,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((face_XYZ_t = (double *)calloc(N_condition_t,sizeof(double))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((type_t = (int *)calloc(N_condition_t,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}

  fp_face_t = fopen("Boundary_condition_traction.dat","r");
  if(fp_face_t == NULL){
    printf("file not open\n");
    return -1;
  }
  for(int i=0;i<2;i++) fscanf(fp_face_t,"%*[^\n]\n");
  for(int i=0;i<N_condition_t;i++) fscanf(fp_face_t,"%d %lf %d\n", &face_dir_t[i], &face_XYZ_t[i], &type_t[i]);

  fclose(fp_face_t);



  // 境界条件を課す面を調べ、データファイルを作る (処理7) //
  // ディリクレ条件が課せられる自由度の数とトラクション条件が課せられる面の数を数える //
  for(int i=0;i<N_point;i++){
    for(int i=0;i<dim;i++) flag_dir[i] = 0;
    for(int j=face_offset[i];j<face_offset[i+1];j++){
      if(face[j] == -1) continue;
      for(int k=0;k<N_condition_D;k++){
        if(flag_dir[fixed_dir[k]] == 1) continue;
        flag_D = 0;
        for(int l=vertex_offset[face[j]];l<vertex_offset[face[j]+1];l++){
          if(fabs(node_XYZ[dim*node[l]+face_dir_D[k]] - face_XYZ_D[k]) > 1.0e-6){
            flag_D = 1;
            break;
          }
        }
        if(flag_D == 0){
          N_Dirichlet_DoF += 1;
          flag_dir[fixed_dir[k]] = 1;
        }
      }
    }
  }
  for(int i=0;i<N_point;i++){
    for(int j=face_offset[i];j<face_offset[i+1];j++){
      if(face[j] == -1) continue;
      for(int k=0;k<N_condition_t;k++){
        flag_t = 0;
        for(int l=vertex_offset[face[j]];l<vertex_offset[face[j]+1];l++){
          if(fabs(node_XYZ[dim*node[l]+face_dir_t[k]] - face_XYZ_t[k]) > 1.0e-6){
            flag_t = 1;
            break;
          }
        }
        if(flag_t == 0){
          N_traction_face += 1;
          break;
        }
      }
    }
  }

  fp_Dirichlet = fopen("Input_Dirichlet_condition.dat","w");
  if(fp_Dirichlet == NULL){
    printf("file not open\n");
    return -1;
  }
  fprintf(fp_Dirichlet,"N_Dirichlet_DoF  %d\n", N_Dirichlet_DoF);
  fprintf(fp_Dirichlet,"point number / fixed direction / type of displacement / surface\n");

  fp_traction = fopen("Input_traction_condition.dat","w");
  if(fp_traction == NULL){
    printf("file not open\n");
    return -1;
  }
  fprintf(fp_traction,"N_traction_face  %d\n", N_traction_face);
  fprintf(fp_traction,"point number / type of traction / surface\n");

  for(int i=0;i<N_point;i++){
    for(int i=0;i<dim;i++) flag_dir[i] = 0;
    for(int j=face_offset[i];j<face_offset[i+1];j++){
      if(face[j] == -1) continue;
      for(int k=0;k<N_condition_D;k++){
        if(flag_dir[fixed_dir[k]] == 1) continue;
        flag_D = 0;
        for(int l=vertex_offset[face[j]];l<vertex_offset[face[j]+1];l++){
          if(fabs(node_XYZ[dim*node[l]+face_dir_D[k]] - face_XYZ_D[k]) > 1.0e-6){
            flag_D = 1;
            break;
          }
        }
        if(flag_D == 0){
          fprintf(fp_Dirichlet,"%7d  %d  %2d  %7d\n", i, fixed_dir[k], type_D[k], face[j]);
          flag_dir[fixed_dir[k]] = 1;
        }
      }
    }
  }
  for(int i=0;i<N_point;i++){
    for(int j=face_offset[i];j<face_offset[i+1];j++){
      if(face[j] == -1) continue;
      for(int k=0;k<N_condition_t;k++){
        flag_t = 0;
        for(int l=vertex_offset[face[j]];l<vertex_offset[face[j]+1];l++){
          if(fabs(node_XYZ[dim*node[l]+face_dir_t[k]] - face_XYZ_t[k]) > 1.0e-6){
            flag_t = 1;
            break;
          }
        }
        if(flag_t == 0){
          fprintf(fp_traction,"%7d  %2d  %7d\n", i, type_t[k], face[j]);
          break;
        }
      }
    }
  }

  fclose(fp_traction);
  fclose(fp_Dirichlet);

  free(type_t);
  free(face_XYZ_t);
  free(face_dir_t);
  free(type_D);
  free(fixed_dir);
  free(face_XYZ_D);
  free(face_dir_D);
  free(neighbor);
  free(neighbor_offset);
  free(node);
  free(vertex_offset);
  free(face);
  free(face_offset);
  free(node_XYZ);

  return 0;
}
