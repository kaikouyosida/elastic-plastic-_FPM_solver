//残差ベクトルにDirichlet条件を付与
void ImposeDirichretResidual(int NR_step);

//係数マトリクスにDirichlet条件を付与
void ImposeDirichletTangentialMatrix();

//固定する変位の設定
double fixed_deformation(double time, double time_end, double x1, double x2, double x3, int type);