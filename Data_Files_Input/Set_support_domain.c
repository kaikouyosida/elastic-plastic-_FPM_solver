#include<stdio.h>
#include<stdlib.h>

int main(){
  // サポートドメインの範囲 //
  // 1ならば、隣あった要素をサポートドメインとする //
  // 2ならば、内部境界を2つ跨いで到達する要素までサポートドメインに含める //
  int support_range = 1;

  int N_point = 0;  // ポイント(サブドメイン)の数

  int sum_N_neighbor = 0;  // 隣接点の数の総和
  int *neighbor;           // 各ポイントの隣接点の番号
  // 各ポイントの隣接点が記録されたneighbor配列の先頭番号 //
  // 最後の要素にはsum_N_neighborを格納 //
  int *neighbor_offset;

  // 複数の処理で使い回す変数 //
  int temp = 0;   // ソート時などに使用する変数
  int count = 0;  // カウンタ



  // モデルの設定ファイルを開く (処理1) //
  FILE *fp_model;  // モデルの設定ファイル

  // 各ポイントの隣接点の情報を読み込む (処理2) //
  int buf = 0;        // 整数値を格納する変数
  FILE *fp_neighbor;  // 隣接点のデータ

  // 各ポイントのサポートドメインを記録 (処理3) //
  int size_support = 0;   // support配列のサイズ
  int num_j = 0;          // 隣接点の番号
  int sum_N_support = 0;  // サポートドメインの数の総和
  int num_k = 0;          // 隣接点に隣接するポイントの番号
  int *support;           // サポートドメインを記録する配列
  int *N_support;         // サポートドメインの数
  // 各ポイントのサポートドメインが記録されたsupport配列の先頭番号 //
  // 最後の要素にはsize_supportを格納 //
  int *support_offset;

  // サポートドメインの情報を出力 (処理4) //
  FILE *fp_support;  // サポートドメインの数と番号を出力するファイル



  // モデルの設定ファイルを開く (処理1) //
  fp_model = fopen("Model_settings.dat","r");
  if(fp_model == NULL){
    printf("file not open\n");
    return -1;
  }
  for(int i=0;i<2;i++) fscanf(fp_model,"%*[^\n]\n");
  fgetc(fp_model);
  for(int i=0;i<3;i++) fscanf(fp_model,"%*[^\n]\n");
  fgetc(fp_model);
  fscanf(fp_model,"%*[^\n]\n");
  fscanf(fp_model,"%d\n", &support_range);

  fclose(fp_model);



  // 各ポイントの隣接点の情報を読み込む (処理2) //
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



  // 各ポイントのサポートドメインを記録 (処理3) //
  if(support_range == 1) size_support = sum_N_neighbor;
  else if(support_range == 2){
    // support配列のサイズを大きめにとる //
    for(int i=0;i<N_point;i++){
      size_support += neighbor_offset[i+1] - neighbor_offset[i];

      for(int j=neighbor_offset[i];j<neighbor_offset[i+1];j++){
        num_j = neighbor[j];
        size_support += neighbor_offset[num_j+1] - neighbor_offset[num_j] - 1;
      }
    }
  }

  if((support = (int *)calloc(size_support,sizeof(int))) == NULL){
    printf("Error:Memory is not enough\n");
    return -1;
  }
  for(int i=0;i<size_support;i++) support[i] = -1;
  if((N_support = (int *)calloc(N_point,sizeof(int))) == NULL){
    printf("Error:Memory is not enough\n");
    return -1;
  }
  if((support_offset = (int *)calloc(N_point+1,sizeof(int))) == NULL){
    printf("Error:Memory is not enough\n");
    return -1;
  }
  support_offset[0] = 0;

  if(support_range == 1){
    for(int i=0;i<N_point;i++){
      // 隣接点をサポートドメインとして記録 //
      for(int j=neighbor_offset[i];j<neighbor_offset[i+1];j++) support[j] = neighbor[j];

      N_support[i] = neighbor_offset[i+1] - neighbor_offset[i];
      support_offset[i+1] = neighbor_offset[i+1];
      sum_N_support += N_support[i];
    }
  }
  else if(support_range == 2){
    for(int i=0;i<N_point;i++){
      count = 0;

      // 隣接点をサポートドメインとして記録 //
      for(int j=0;j<neighbor_offset[i+1]-neighbor_offset[i];j++) support[support_offset[i]+j] = neighbor[neighbor_offset[i]+j];
      count += neighbor_offset[i+1] - neighbor_offset[i];

      // 隣接点の隣接点をサポートドメインとして記録 //
      for(int j=neighbor_offset[i];j<neighbor_offset[i+1];j++){
        num_j = neighbor[j];

        for(int k=neighbor_offset[num_j];k<neighbor_offset[num_j+1];k++){
          num_k = neighbor[k];
          if(num_k == i) continue;

          support[support_offset[i]+count] = num_k;
          count += 1;
        }
      }

      // support配列をソート //
      for(int j=0;j<count-1;j++){
        for(int k=j+1;k<count;k++){
          if(support[support_offset[i]+k] < support[support_offset[i]+j]){
            temp = support[support_offset[i]+j];
            support[support_offset[i]+j] = support[support_offset[i]+k];
            support[support_offset[i]+k] = temp;
          }
        }
      }

      // support配列の重複要素に注意して、サポートドメインの数を調べる //
      temp = -1;
      for(int j=0;j<count;j++){
        if(support[support_offset[i]+j] > temp){
          temp = support[support_offset[i]+j];
          N_support[i] += 1;
        }
      }

      support_offset[i+1] = support_offset[i] + count;
      sum_N_support += N_support[i];
    }
  }
  count = 0;



  // サポートドメインの情報を出力 (処理4) //
  fp_support = fopen("Input_support_domain.dat","w");
  if(fp_support == NULL){
    printf("file is not open\n");
    return -1;
  }
  fprintf(fp_support,"N_point  %d\n", N_point);
  fprintf(fp_support,"sum_N_support  %d\n", sum_N_support);
  fprintf(fp_support,"point number / number of support / support domain number\n");
  for(int i=0;i<N_point;i++){
    fprintf(fp_support,"%7d  %2d ", i, N_support[i]);

    temp = -1;
    for(int j=support_offset[i];j<support_offset[i+1];j++){
      if(support[j] > temp){
        temp = support[j];
        fprintf(fp_support," %7d", support[j]);
      }
    }
    fprintf(fp_support,"\n");
  }

  fclose(fp_support);
  free(support_offset);
  free(N_support);
  free(support);
  free(neighbor);
  free(neighbor_offset);

  return 0;
}
