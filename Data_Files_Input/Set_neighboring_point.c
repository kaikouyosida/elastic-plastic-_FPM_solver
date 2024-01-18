#include<stdio.h>
#include<stdlib.h>

int main(){
  int N_subdomain = 0;  // サブドメイン(ポイント)の数

  int sum_N_face = 0;  // サブドメインがもつ面の数の総和
  // 各サブドメインがもつ面が記録されたface配列の先頭番号 //
  // 最後の要素にはsum_N_faceを格納 //
  int *face_offset;
  int *face;           // 各サブドメインがもつ面の番号

  int sum_N_neighbor = 0;  // 隣接点の数の総和
  int *N_neighbor;         // 各ポイントに隣接しているポイントの数
  int *neighbor;           // 隣接点の番号

  int N_int_boundary = 0;  // 内部境界の数
  int *Data_int_boundary;  // ポイント番号、隣接点番号、内部境界(面)の番号を順に入れた配列

  // 複数の処理で使い回す変数 //
  int count = 0;   // カウンタ



  // 面のデータを開き、サブドメインと面の関係を読み込む (処理1) //
  int buf = 0;    // 整数値を格納する変数
  FILE *fp_face;  // 面のデータ

  // 各サブドメインに隣接するサブドメインと、共有している面を調べる (処理2) //
  int count_lim = 0;  // カウンタの上限

  // 内部境界の情報を出力 (処理3) //
  FILE *fp_internal_boundary;  // 内部境界の情報を出力するファイル

  // 隣接点の情報を出力 (処理4) //
  FILE *fp_neighbor;  // 隣接点の数と番号を出力するファイル



  // 面のデータを開き、サブドメインと面の関係を読み込む (処理1) //
  fp_face = fopen("Input_face_of_subdomain.dat","r");
  if(fp_face == NULL){
    printf("file not open\n");
    return -1;
  }

  fscanf(fp_face,"%*s %d\n", &N_subdomain);
  fscanf(fp_face,"%*s %d\n", &sum_N_face);
  fclose(fp_face);

	if((face_offset = (int *)calloc(N_subdomain+1,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
	if((face = (int *)calloc(sum_N_face,sizeof(int))) == NULL){
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
  for(int i=0;i<N_subdomain;i++){
    fscanf(fp_face,"%*d %d", &buf);
    for(int j=0;j<buf;j++){
      fscanf(fp_face,"%d", &face[count]);
      count += 1;
    }
    fgetc(fp_face);

    face_offset[i+1] = face_offset[i] + buf;
  }

  fclose(fp_face);
  count = 0;



  // 各サブドメインに隣接するサブドメインと、共有している面を調べる (処理2) //
	if((N_neighbor = (int *)calloc(N_subdomain,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
	if((neighbor = (int *)calloc(sum_N_face,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}
  if((Data_int_boundary = (int *)calloc(sum_N_face*3/2,sizeof(int))) == NULL){
		printf("Error:Memory is not enough\n");
		return -1;
	}

  for(int i=0;i<N_subdomain-1;i++){
    count = 0;
    count_lim = face_offset[i+1] - face_offset[i];

    // サブドメインiより番号が大きいサブドメインの中から、面を共有しているものを探す //
    for(int j=i+1;j<N_subdomain && count<count_lim;j++){
      for(int k=face_offset[i];k<face_offset[i+1];k++){
        for(int l=face_offset[j];l<face_offset[j+1];l++){
          if(face[k] == face[l]){
            Data_int_boundary[3*N_int_boundary] = i;
            Data_int_boundary[3*N_int_boundary+1] = j;
            Data_int_boundary[3*N_int_boundary+2] = face[k];
            N_int_boundary += 1;

            neighbor[face_offset[i]+N_neighbor[i]] = j;
            N_neighbor[i] += 1;
            neighbor[face_offset[j]+N_neighbor[j]] = i;
            N_neighbor[j] += 1;

            sum_N_neighbor += 2;
            count += 1;
          }
        }
      }
    }
  }



  // 内部境界の情報を出力 (処理3) //
  fp_internal_boundary = fopen("Input_internal_boundary.dat", "w");
  if(fp_internal_boundary == NULL){
    printf("file is not open\n");
    return -1;
  }
  fprintf(fp_internal_boundary,"N_internal_boundary  %d\n", N_int_boundary);
  fprintf(fp_internal_boundary,"main point / neighbor point / shared face\n");
  for(int i=0;i<N_int_boundary;i++) fprintf(fp_internal_boundary,"%7d  %7d  %7d\n", Data_int_boundary[3*i], Data_int_boundary[3*i+1], Data_int_boundary[3*i+2]);

  fclose(fp_internal_boundary);
  free(Data_int_boundary);
  count = 0;



  // 隣接点の情報を出力 (処理4) //
  fp_neighbor = fopen("Input_neighboring_point.dat","w");
  if(fp_neighbor == NULL){
    printf("file is not open\n");
    return -1;
  }
  fprintf(fp_neighbor,"N_point  %d\n", N_subdomain);
  fprintf(fp_neighbor,"sum_N_neighbor  %d\n", sum_N_neighbor);
  fprintf(fp_neighbor,"point number / number of neighbor / neighboring point number\n");

  for(int i=0;i<N_subdomain;i++){
    fprintf(fp_neighbor,"%7d  %2d ", i, N_neighbor[i]);
    for(int j=0;j<N_neighbor[i];j++) fprintf(fp_neighbor," %7d", neighbor[face_offset[i]+j]);
    fprintf(fp_neighbor,"\n");
  }

  fclose(fp_neighbor);
  free(neighbor);
  free(N_neighbor);
  free(face);
  free(face_offset);

  return 0;
}
