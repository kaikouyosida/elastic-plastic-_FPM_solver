//テンソルの指数テンソルを計算する
void calculateTensorExponent(double (*tensor_out)[3], double (*tensor_in)[3]);

//テンソルの対数テンソルを計算する
void calculateTensorLogarithm(double (*tensor_out)[3], double (*tensor_in)[3]);

//4階のテンソルをマトリクスからテンソルへ変換する
void convertSymmetric4thOrderMatrixToTensor(double tensor[3][3][3][3],
                                            double matrix[6][6]);

//4階の対象テンソルからマトリクスへ変換する                                           
void convertSymmetric4thOrderTensorToMatrix(double matrix[6][6],
                                            double tensor[3][3][3][3]);

//4階のテンソルからマトリクスへ変換する                                           
void conver4thOrderTensorToMatrix(double matrix[9][9], double tensor[3][3][3][3]);

//対数テンソルの微分を計算する
void calculateTensorLogarithmDerivative(double tensor_out[3][3][3][3],
                                        double tensor_in[3][3]);

//等方テンソルの微分を計算
void calculateIsotropicTensorFunctionDerivative(double tensor_out[3][3][3][3],
                                                double tensor_in[3][3],
                                                double (*function)(double variable),
                                                double (*function_derivative)(double variable));

//固有値を計算する                                              
void calculateEigenvalues(double eigenvalues[3],
                          double tensor[3][3]);

//固有写像を計算する
void calculateEigenprojections(double eigenprojections[3][3][3],
                               double eigenvalues[3],
                               double tensor[3][3]);