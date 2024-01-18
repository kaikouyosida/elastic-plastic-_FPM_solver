#include<stdio.h>
#include<stdlib.h>

int main(){
  int N_point = 0;  // ポイントの数

  int sum_N_support = 0;  // サポートドメインの数の総和
  int *support;           // 各ポイントのサポートドメインの番号
  // 各ポイントのサポートドメインが記録されたsupport配列の先頭番号 //
  // 最後の要素にはsum_N_supportを格納 //
  int *support_offset;

  int sum_N_neighbor = 0;  // 隣接点の数の総和
  int *neighbor;           // 各ポイントの隣接点の番号
  // 各ポイントの隣接点が記録されたneighbor配列の先頭番号 //
  // 最後の要素にはsum_N_neighborを格納 //
  int *neighbor_offset;

  // 複数の処理で使い回す変数 //
  int buf = 0;    // 整数値を格納する変数
  int count = 0;  // カウンタ
  int temp = 0;   // ソート時などに使用する変数



  // サポートドメインの情報を読み込む (処理1) //
  FILE *fp_support;   // サポートドメインの番号のデータ

  // 各ポイントの隣接点の情報を読み込む (処理2) //
  FILE *fp_neighbor;   // 隣接点のデータ

  // 各ポイントの間接サポートドメインを記録 (処理3) //
  int size_ind_support = 0;   // ind_support配列のサイズ
  int num_j = 0;              // サポートドメインの番号
  int num_k = 0;              // サポートドメインの隣接点の番号
  int swit = 0;               // 番号を記録する判定
  int sum_N_ind_support = 0;  // 間接サポートドメインの数の総和
  int *ind_support;           // 間接サポートドメインを記録する配列
  int *N_ind_support;         // 間接サポートドメインの数
  // 各ポイントのサポートドメインが記録されたind_support配列の先頭番号 //
  // 最後の要素にはsize_ind_supportを格納 //
  int *ind_support_offset;

  // 間接サポートドメインの情報を出力 (処理4) //
  FILE *fp_ind_support;  // 間接サポートドメインの数と番号を出力するファイル



  // サポートドメインの情報を読み込む (処理1) //
  fp_support = fopen("Input_support_domain.dat","r");
  if(fp_support == NULL){
    printf("file not open\n");
    return -1;
  }

  fscanf(fp_support,"%*s %d\n", &N_point);
  fscanf(fp_support,"%*s %d\n", &sum_N_support);
  fclose(fp_support);

  if((support_offset = (int *)calloc(N_point+1,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
	if((support = (int *)calloc(sum_N_support,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}

  fp_support = fopen("Input_support_domain.dat","r");
  if(fp_support == NULL){
    printf("file not open\n");
    return -1;
  }
  for(int i=0;i<3;i++) fscanf(fp_support,"%*[^\n]\n");

  support_offset[0] = 0;
  for(int i=0;i<N_point;i++){
    fscanf(fp_support,"%*d %d", &buf);

    support_offset[i+1] = support_offset[i] + buf;
    for(int j=0;j<buf;j++){
      fscanf(fp_support,"%d", &support[count]);
      count += 1;
    }
    fgetc(fp_support);
  }

  fclose(fp_support);
  count = 0;



  // 各ポイントの隣接点の情報を読み込む (処理2) //
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



  // 各ポイントの間接サポートドメインを記録 (処理3) //
  // ind_support配列のサイズを大きめにとる //
  for(int i=0;i<N_point;i++){
    for(int j=support_offset[i];j<support_offset[i+1];j++){
      num_j = support[j];

      size_ind_support += support_offset[num_j+1] - support_offset[num_j] - 1;

      for(int k=neighbor_offset[num_j];k<neighbor_offset[num_j+1];k++){
        num_k = neighbor[k];

        size_ind_support += support_offset[num_k+1] - support_offset[num_k] - 1;
      }
    }
  }

  if((ind_support = (int *)calloc(size_ind_support,sizeof(int))) == NULL){
    printf("Error:Memory is not enough\n");
    return -1;
  }
  for(int i=0;i<size_ind_support;i++) ind_support[i] = -1;
  if((N_ind_support = (int *)calloc(N_point,sizeof(int))) == NULL){
    printf("Error:Memory is not enough\n");
    return -1;
  }
  if((ind_support_offset = (int *)calloc(N_point+1,sizeof(int))) == NULL){
    printf("Error:Memory is not enough\n");
    return -1;
  }
  ind_support_offset[0] = 0;

  for(int i=0;i<N_point;i++){
    count = 0;

    // サポートドメインのサポートドメインを記録 //
    for(int j=support_offset[i];j<support_offset[i+1];j++){
      num_j = support[j];
      for(int k=support_offset[num_j];k<support_offset[num_j+1];k++){
        if(support[k] == i) continue;
        for(int l=support_offset[i];l<support_offset[i+1];l++){
          if(support[k] == support[l]){
            swit = 1;
            break;
          }
        }
        if(swit == 1){
          swit = 0;
          continue;
        }

        ind_support[ind_support_offset[i]+count] = support[k];
        count += 1;
      }
    }

    // サポートドメインの隣接点のサポートドメインを記録 //
    for(int j=support_offset[i];j<support_offset[i+1];j++){
      num_j = support[j];
      for(int k=neighbor_offset[num_j];k<neighbor_offset[num_j+1];k++){
        num_k = neighbor[k];
        for(int l=support_offset[num_k];l<support_offset[num_k+1];l++){
          if(support[l] == i) continue;
          for(int m=support_offset[i];m<support_offset[i+1];m++){
            if(support[l] == support[m]){
              swit = 1;
              break;
            }
          }
          if(swit == 1){
            swit = 0;
            continue;
          }

          ind_support[ind_support_offset[i]+count] = support[l];
          count += 1;
        }
      }
    }

    // ind_support配列をソート //
    for(int j=0;j<count-1;j++){
      for(int k=j+1;k<count;k++){
        if(ind_support[ind_support_offset[i]+k] < ind_support[ind_support_offset[i]+j]){
          temp = ind_support[ind_support_offset[i]+j];
          ind_support[ind_support_offset[i]+j] = ind_support[ind_support_offset[i]+k];
          ind_support[ind_support_offset[i]+k] = temp;
        }
      }
    }

    // ind_support配列の重複要素に注意して、間接サポートドメインの数を調べる //
    temp = -1;
    for(int j=0;j<count;j++){
      if(ind_support[ind_support_offset[i]+j] > temp){
        temp = ind_support[ind_support_offset[i]+j];
        N_ind_support[i] += 1;
      }
    }

    ind_support_offset[i+1] = ind_support_offset[i] + count;
    sum_N_ind_support += N_ind_support[i];
  }
  count = 0;



  // 間接サポートドメインの情報を出力 (処理4) //
  fp_ind_support = fopen("Input_indirect_support_domain.dat","w");
  if(fp_ind_support == NULL){
    printf("file is not open\n");
    return -1;
  }
  fprintf(fp_ind_support,"N_point  %d\n", N_point);
  fprintf(fp_ind_support,"sum_N_indirect_support  %d\n", sum_N_ind_support);
  fprintf(fp_ind_support,"point number / number of indirect support / indirect support domain number\n");
  for(int i=0;i<N_point;i++){
    fprintf(fp_ind_support,"%7d  %2d ", i, N_ind_support[i]);

    temp = -1;
    for(int j=ind_support_offset[i];j<ind_support_offset[i+1];j++){
      if(ind_support[j] > temp){
        temp = ind_support[j];
        fprintf(fp_ind_support," %7d", ind_support[j]);
      }
    }
    fprintf(fp_ind_support,"\n");
  }

  fclose(fp_ind_support);
  free(ind_support_offset);
  free(N_ind_support);
  free(ind_support);
  free(neighbor);
  free(neighbor_offset);
  free(support);
  free(support_offset);

  return 0;
}
