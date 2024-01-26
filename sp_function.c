#include<stdio.h>

void make_sp_der(double s,double t,double (*shape_func_der)[4]){
    shape_func_der[0][0]=(1.0+t)/4.0;shape_func_der[0][1]=(1.0-t)/4.0;shape_func_der[0][2]=(-1.0+t)/4.0;shape_func_der[0][3]=(-1.0-t)/4.0;
    shape_func_der[1][0]=(1.0+s)/4.0;shape_func_der[1][1]=(-1.0-s)/4.0;shape_func_der[1][2]=(-1.0+s)/4.0;shape_func_der[1][3]=(1.0-s)/4.0;
}

void make_sp(double s, double t, double (*shape)[8]){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 8; j++) shape[i][j] =0.;
    }
    double f[4];
     f[0] = 0.25 * (1.0 + s) * (1.0 + t);
     f[1] = 0.25 * (1.0 + s) * (1.0 - t);
     f[2] = 0.25 * (1.0 - s) * (1.0 - t);
     f[3] = 0.25 * (1.0 - s) * (1.0 + t);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 2; j++){
            shape[j][2 * i + j] = f[i];
        }
    }
}