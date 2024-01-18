#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(){
  int dim = 0;  // モデルの次元

  int N_element = 0;     // 要素(サブドメイン)の数
  int type_element = 0;  // 要素の種類
  // 2次元モデルの場合、要素の頂点数は3(三角形)か4(四角形) //
  // 3次元モデルの場合、要素の頂点数は4(四面体)か8(六面体) //
  int vertex = 0;  // 1つの要素の頂点の数

  int N_node = 0;  // 節点(頂点)の数



  // モデルの設定ファイルを開く (処理1) //
  FILE *fp_model;  // モデルの設定ファイル

  // Marcのモデルデータを開く (処理2) //
  int buf = 0;         // 整数を格納する変数
  char ss1[512];       // 文字列を格納する変数
  char ss2[] = "connectivity";
  double X = 0.; int Xe = 0;  // x座標の仮数部と指数部
  double Y = 0.; int Ye = 0;  // y座標の仮数部と指数部
  double Z = 0.; int Ze = 0;  // z座標の仮数部と指数部
  FILE *fp_Marc;       // Marcのデータファイル
  FILE *fp_subdomain;  // サブドメインの情報を出力するファイル



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
  fscanf(fp_model,"%d\n", &type_element);

  fclose(fp_model);

  if(type_element == 0) vertex = 3;
  else if(type_element == 1) vertex = 4;
  else if(type_element == 2) vertex = 4;
  else vertex = 8;



  // Marcのモデルデータを開く (処理2) //
  fp_Marc = fopen("Marc_data.dat","r");
  if(fp_Marc == NULL){
    printf("file not open\n");
    return -1;
  }
  fp_subdomain = fopen("Input_subdomain.dat","w");
  if(fp_subdomain == NULL){
    printf("file not open\n");
    return -1;
  }

  // 要素数と節点数が書かれている場所まで空読み //
  for(int i=0;i<6;i++) fscanf(fp_Marc,"%*[^\n]\n");
  // 要素数と節点数を格納 //
  fscanf(fp_Marc,"%*s %*d %d %d %*d\n", &N_element, &N_node);

  // Connectivityが書かれている場所まで空読み //
  while(strncmp(ss1,ss2,12) != 0){
    fscanf(fp_Marc,"%[^\n]\n",ss1);
  }

  fprintf(fp_subdomain,"N_subdomain  %d\n", N_element);
  fprintf(fp_subdomain,"sum_N_vertex  %d\n", vertex*N_element);
  fprintf(fp_subdomain,"N_node  %d\n", N_node);
  fprintf(fp_subdomain,"subdomain number / vertex / node number\n");
  fscanf(fp_Marc,"%*[^\n]\n");
  // 節点番号を格納 //
  for(int i=0;i<N_element;i++){
    fscanf(fp_Marc,"%*d %*d");
    fprintf(fp_subdomain,"%7d  %2d ", i, vertex);
    for(int j=0;j<vertex-1;j++){
      fscanf(fp_Marc,"%d", &buf);
      fprintf(fp_subdomain," %7d", buf-1);
    }
    // 改行文字単独の読み込みができなかったのでこの処理になってる //
    fscanf(fp_Marc,"%d\n", &buf);
    fprintf(fp_subdomain," %7d\n", buf-1);
  }

  fscanf(fp_Marc,"%*[^\n]\n");
  fscanf(fp_Marc,"%*[^\n]\n");

  fprintf(fp_subdomain,"\n");
  // 座標を読み込んで出力 //
  if(dim == 2){
    fprintf(fp_subdomain,"node number / coordinate  X Y\n");
    for(int i=0;i<N_node;i++){
      fscanf(fp_Marc,"%*d %lf%d %lf%d %lf%d\n", &X, &Xe, &Y, &Ye, &Z, &Ze);
      fprintf(fp_subdomain,"%7d  %+.15lfe%+d  %+.15lfe%+d\n", i, X, Xe, Y, Ye);
    }
  }
  else if(dim == 3){
    fprintf(fp_subdomain,"node number / coordinate  X Y Z\n");
    for(int i=0;i<N_node;i++){
      fscanf(fp_Marc,"%*d %lf%d %lf%d %lf%d\n", &X, &Xe, &Y, &Ye, &Z, &Ze);
      fprintf(fp_subdomain,"%7d  %+.15lfe%+d  %+.15lfe%+d  %+.15lfe%+d\n", i, X, Xe, Y, Ye, Z, Ze);
    }
  }

  fclose(fp_subdomain);
  fclose(fp_Marc);

  return 0;
}
