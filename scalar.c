#include<stdio.h>
#include<math.h>
#include"scalar.h"
#include"type.h"

extern Global global;
extern Option option;

double calculateInverse(const double value)
{
    return 1.0 / value;
}

/*
 * Calculate square
 */
double calculateSquare(const double value)
{
    return value * value;
}

/*
 * Calculate cube
 */
double calculateCube(const double value)
{
    return value * value * value;
}

/*
 * Swap reals
 */
void swapReals(double *value1, double *value2)
{
    double temp = *value1;
    *value1  = *value2;
    *value2  = temp;
}

/*
 * Swap integers
 */
void swapIntegers(int *value1, int *value2)
{
    int temp = *value1;
    *value1  = *value2;
    *value2  = temp;
}
/*
* Calculate subdomain volume
*/
double calc_subdomain_volume(int point_n){
    double volume = 0.;
    double center[3];
    int N_face = global.subdomain.face_offset[point_n + 1] - global.subdomain.face_offset[point_n];
    int subdomain_node[60];

    for(int i = 0; i < 8; i++)
        subdomain_node[i] = global.subdomain.subdomain_node[8 * point_n + i];
    
    for(int i = 0; i < 3; i++){
        double center_i = 0.;
        for(int j = 0; j < 8; j++){
            center_i += (global.subdomain.node_XYZ[option.dim * subdomain_node[8 * point_n + j] + i]
                     + global.subdomain.nodal_displacements[subdomain_node[8 * point_n + j]][i]
                     + global.subdomain.nodal_displacement_increments[subdomain_node[8 * point_n + j]][i]);
        }
        center[i] = center_i / 8.0;
    }
    for(int i = 0; i < N_face; i++){
        int ref_num = global.subdomain.vertex_offset[global.subdomain.face[global.subdomain.face_offset[point_n] + i]];
        
    }
}
double dot_product(int N, double *vec1, double *vec2)
{
	double X = 0.;
	for (int i = 0; i < N; i++)
		X += vec1[i] * vec2[i];
	return X;
}
void cross_product(int dim, double *vecA, double *vecB, double *AcrossB)
{
	if (dim == 2)
		AcrossB[0] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
	else if (dim == 3)
	{
		AcrossB[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
		AcrossB[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
		AcrossB[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
	}
}
