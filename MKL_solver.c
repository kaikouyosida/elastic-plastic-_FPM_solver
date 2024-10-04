#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>

void LU_dens(double *global_k, double *du, double *Q, int DoF_free){
    // 行列のサイズ
    int N = DoF_free;   //マトリクスの自由度数
    int LDA = N;        //global_kの先頭次元
    int LDB = N;        //Qの先頭次元
    int NHRS = 1;       //右辺ベクトルの数（複数パターンの方程式を一回で解く場合は1から変更）
    int INFO;           //実行結果

    int *IPIV = (int *)calloc(N, sizeof(int));      //ピボット結果を格納する配列
    
    if (IPIV == NULL) {
        printf("Error: Memory is not enough\n");
        exit(-1);
    }

    // global_k, Q, duのポインタがNULLでないか確認
    if (global_k == NULL || Q == NULL || du == NULL) {
        printf("Error: One of the input arrays is NULL\n");
        free(IPIV);
        exit(-1);
    }
    
    // LAPACK の dgesv 関数を使用して連立一次方程式を解く
    dgesv_(&N, &NHRS, global_k, &LDA, IPIV, Q, &LDB, &INFO);

    // INFOの値を確認
    if (INFO == 0) {
        for (int i = 0; i < N; i++) {
            du[i] = Q[i];
        }
    } else {
        printf("An error occurred: %d (INFO code)\n", INFO);
    }

    free(IPIV);
}

void Paradiso(int n, int nnz, double *a, int *ia, int *ja, double *b, double *x){
    void *pt[64] = {0};
    int iparm[64] = {0};
    int maxfct = 1, mnum = 1, phase = 13, error = 0, msglvl = 0;
    int nrhs = 1;  // 右辺ベクトルの数
    int mtype = 11;  // 実数の非対称行列

    // PARDISOパラメータの初期化
    for (int i = 0; i < 64; i++) {
        iparm[i] = 0;
    }
    // PARDISOパラメータの初期化
    for (int i = 0; i < 64; i++) {
        iparm[i] = 0;
    }

    // スレッド数の設定
    mkl_set_num_threads(4);  // 4スレッドを使用
    
    iparm[0] = 1;  // デフォルト以外の値を使用
    iparm[1] = 2;  // フィル減少順序付け
    iparm[3] = 0;  // 反復改良の回数
    iparm[4] = 0;  // ユーザー並べ替えを使用しない
    iparm[5] = 0;  // 解ベクトルを書き出す
    iparm[7] = 2;  // 反復改良の最大回数
    iparm[9] = 13; // ピボット摂動の大きさ
    iparm[10] = 1; // スケーリングベクトル
    iparm[12] = 1; // 改良されたウェイト付きマッチング
    iparm[17] = -1; // 報告される非ゼロ要素数
    iparm[18] = -1; // 報告される反復回数
    iparm[24] = 1;  // 並列数値因数分解
    iparm[25] = 1;  // 並列前進/後退代入
    iparm[26] = 1;  //チェッカーを有効にする
    iparm[34] = 1;  // 0ベースインデックスを使用

    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, b, x, &error);

    if (error != 0) {
        printf("PARDISO solver error: %d\n", error);
        exit(-1);
    }

    // PARDISOの終了処理
    phase = -1;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, NULL, &nrhs, iparm, &msglvl, b, x, &error);
}