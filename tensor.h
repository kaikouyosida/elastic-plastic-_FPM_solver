void calculateTensorExponent(double (*tensor_out)[3], double (*tensor_in)[3]);
void calculateTensorLogarithm(double (*tensor_out)[3], double (*tensor_in)[3]);
void convertSymmetric4thOrderMatrixToTensor(double tensor[3][3][3][3],
                                            double matrix[6][6]);
void convertSymmetric4thOrderTensorToMatrix(double matrix[6][6],
                                            double tensor[3][3][3][3]);
void conver4thOrderTensorToMatrix(double matrix[9][9], double tensor[3][3][3][3]);
void calculateTensorLogarithmDerivative(double tensor_out[3][3][3][3],
                                        double tensor_in[3][3]);
                                        
void calculateIsotropicTensorFunctionDerivative(double tensor_out[3][3][3][3],
                                                double tensor_in[3][3],
                                                double (*function)(const double variable),
                                                double (*function_derivative)(const double variable));
                                                
void calculateEigenvalues(double eigenvalues[3],
                          double tensor[3][3]);
void calculateEigenprojections(double eigenprojections[3][3][3],
                               const double eigenvalues[3],
                               double tensor[3][3]);

