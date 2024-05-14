#include<stdio.h>
#include<math.h>
#include"scalar.h"
#include"tensor.h"
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

double calculate3x3MatrixTrace(double matrix[3][3])
{
    return matrix[0][0] + matrix[1][1] + matrix[2][2];
}
/*
* Calculate subdomain volume(only the case for hexahedral...)
*/
double calc_subdomain_volume(int point_n){
    double volume = 0.;
    double center[3];
    int N_face = global.subdomain.face_offset[point_n + 1] - global.subdomain.face_offset[point_n];
    int subdomain_node[60];
    double edge1[3],edge2[3],edge3[3];   //辺のベクトル
    double edge1_cross_edge2[3];         //辺のベクトルの外積

    for(int i = 0; i < 8; i++)
        subdomain_node[i] = global.subdomain.subdomain_node[8 * point_n + i];
   
    for(int i = 0; i < 3; i++){
        double center_i = 0.;
        for(int j = 0; j < 8; j++){
            center_i += (global.subdomain.node_XYZ[option.dim * subdomain_node[j] + i]
                     + global.subdomain.nodal_displacements[subdomain_node[j]][i]
                     + global.subdomain.nodal_displacement_increments[subdomain_node[j]][i]);
        }
        center[i] = center_i / 8.0;
    }
    for(int i = 0; i < N_face; i++){
        int ref_num = global.subdomain.vertex_offset[global.subdomain.face[global.subdomain.face_offset[point_n] + i]];
        int node1 = global.subdomain.node[ref_num];
        int node2 = global.subdomain.node[ref_num + 1];
        int node3 = global.subdomain.node[ref_num + 2];
        int node4 = global.subdomain.node[ref_num + 3];

        for(int j = 0; j < option.dim; j++){
            edge1[j] = (global.subdomain.node_XYZ[option.dim * node4 + j] - global.subdomain.node_XYZ[option.dim * node1 + j]
                     + global.subdomain.nodal_displacements[node4][j] - global.subdomain.nodal_displacements[node1][j]
                     + global.subdomain.nodal_displacement_increments[node4][j] - global.subdomain.nodal_displacement_increments[node1][j]);
            edge2[j] = (global.subdomain.node_XYZ[option.dim * node2 + j] - global.subdomain.node_XYZ[option.dim * node1 + j]
                     + global.subdomain.nodal_displacements[node2][j] - global.subdomain.nodal_displacements[node1][j]
                     + global.subdomain.nodal_displacement_increments[node2][j] - global.subdomain.nodal_displacement_increments[node1][j]);
            edge3[j] = center[j] - global.subdomain.node_XYZ[option.dim * node1 + j]
                     - global.subdomain.nodal_displacements[node1][j]
                     - global.subdomain.nodal_displacement_increments[node1][j];    
        }
        cross_product(option.dim, edge1, edge2, edge1_cross_edge2);
        volume += fabs(dot_product(option.dim,edge1_cross_edge2, edge3));
        for(int j = 0; j < option.dim; j++){
            edge1[j] = (global.subdomain.node_XYZ[option.dim * node4 + j] - global.subdomain.node_XYZ[option.dim * node3 + j]
                     + global.subdomain.nodal_displacements[node4][j] - global.subdomain.nodal_displacements[node3][j]
                     + global.subdomain.nodal_displacement_increments[node4][j] - global.subdomain.nodal_displacement_increments[node3][j]);
            edge2[j] = (global.subdomain.node_XYZ[option.dim * node2 + j] - global.subdomain.node_XYZ[option.dim * node3 + j]
                     + global.subdomain.nodal_displacements[node2][j] - global.subdomain.nodal_displacements[node3][j]
                     + global.subdomain.nodal_displacement_increments[node2][j] - global.subdomain.nodal_displacement_increments[node3][j]);
            edge3[j] = center[j] - global.subdomain.node_XYZ[option.dim * node3 + j]
                     - global.subdomain.nodal_displacements[node3][j]
                     - global.subdomain.nodal_displacement_increments[node3][j]; 
        }
        cross_product(option.dim, edge1, edge2, edge1_cross_edge2);
        volume += fabs(dot_product(option.dim,edge1_cross_edge2, edge3));
    }

        volume /= 6.0;
        return volume;
}
#if 0
//内部境界を二分割して和をとることで面積を求める
double calc_surface_area(int face_n){
    int ref_num = global.subdomain.vertex_offset[face_n];
    int subdomain_node[60];
    double edge1[3],edge2[3];           //辺のベクトル
    double edge1_cross_edge2[3];        //ベクトル同士の外積
    double area1, area2;                //2つの分割したそれぞれの面積

    int node1 = global.subdomain.node[ref_num];
    int node2 = global.subdomain.node[ref_num + 1];
    int node3 = global.subdomain.node[ref_num + 2];
    int node4 = global.subdomain.node[ref_num + 3];

    for(int i = 0; i < option.dim; i++){
        edge1[i] = (global.subdomain.node_XYZ[option.dim * node4 + i] - global.subdomain.node_XYZ[option.dim * node1 + i]
                    + global.subdomain.nodal_displacements[node4][i] - global.subdomain.nodal_displacements[node1][i]
                    + global.subdomain.nodal_displacement_increments[node4][i] - global.subdomain.nodal_displacement_increments[node1][i]);
        edge2[i] = (global.subdomain.node_XYZ[option.dim * node2 + i] - global.subdomain.node_XYZ[option.dim * node1 + i]
                    + global.subdomain.nodal_displacements[node2][i] - global.subdomain.nodal_displacements[node1][i]
                    + global.subdomain.nodal_displacement_increments[node2][i] - global.subdomain.nodal_displacement_increments[node1][i]);
    }
    cross_product(option.dim, edge1, edge2, edge1_cross_edge2);
    area1 = 0.5 * norm(edge1_cross_edge2, option.dim);

    for(int i = 0; i < option.dim; i++){
        edge1[i] = (global.subdomain.node_XYZ[option.dim * node4 + i] - global.subdomain.node_XYZ[option.dim * node3 + i]
                    + global.subdomain.nodal_displacements[node4][i] - global.subdomain.nodal_displacements[node3][i]
                    + global.subdomain.nodal_displacement_increments[node4][i] - global.subdomain.nodal_displacement_increments[node3][i]);
        edge2[i] = (global.subdomain.node_XYZ[option.dim * node2 + i] - global.subdomain.node_XYZ[option.dim * node3 + i]
                    + global.subdomain.nodal_displacements[node2][i] - global.subdomain.nodal_displacements[node3][i]
                    + global.subdomain.nodal_displacement_increments[node2][i] - global.subdomain.nodal_displacement_increments[node3][i]);
    }
    cross_product(option.dim, edge1, edge2, edge1_cross_edge2);
    area2 = 0.5 * norm(edge1_cross_edge2, option.dim);

    return area1 + area2;
}
#endif
//物理空間座標→正規化座標に変換するためのスカラー値を計算
double calc_area_change(int face_n, int s, int t, double *X){
    int ref_num = global.subdomain.vertex_offset[face_n];
    double scalar = 0.;
    double area_vector[3];
    double e[3][3][3];
    double dx_ds[3];
    double dx_dt[3];
    double x1[3];
    double x2[3];
    double x3[3];
    double x4[3];
    int node1 = global.subdomain.node[ref_num];
    int node2 = global.subdomain.node[ref_num + 1];
    int node3 = global.subdomain.node[ref_num + 2];
    int node4 = global.subdomain.node[ref_num + 3];

    alternating_matrix(e);
    for(int i = 0; i < option.dim; i++){
        x1[i] = global.subdomain.node_XYZ[option.dim * node1 + i] + global.subdomain.nodal_displacements[node1][i] + global.subdomain.nodal_displacement_increments[node1][i];
        x2[i] = global.subdomain.node_XYZ[option.dim * node2 + i] + global.subdomain.nodal_displacements[node2][i] + global.subdomain.nodal_displacement_increments[node2][i];
        x3[i] = global.subdomain.node_XYZ[option.dim * node3 + i] + global.subdomain.nodal_displacements[node3][i] + global.subdomain.nodal_displacement_increments[node3][i];
        x4[i] = global.subdomain.node_XYZ[option.dim * node4 + i] + global.subdomain.nodal_displacements[node4][i] + global.subdomain.nodal_displacement_increments[node4][i];
    }
    for(int i = 0; i < option.dim ; i++){
        dx_ds[i] = -0.25 * (1.0 - X[t]) * x1[i] - 0.25 * (1.0 + X[t]) * x2[i] + 0.25 * (1.0 + X[t]) * x3[i] + 0.25 * (1.0 - X[t]) * x4[i];
        dx_dt[i] = -0.25 * (1.0 - X[s]) * x1[i] + 0.25 * (1.0 - X[s]) * x2[i] + 0.25 * (1.0 + X[s]) * x3[i] - 0.25 * (1.0 + X[s]) * x4[i];
    }
    for(int i = 0; i < option.dim; i++){
        double area_vector_i = 0.;
        for(int j = 0; j < option.dim; j++){
            for(int k = 0; k < option.dim; k++){
                area_vector_i += dx_ds[j] * dx_dt[k] * e[j][k][i];
            }
        }
        area_vector[i] = area_vector_i;
    }
    for(int i = 0; i < option.dim; i++)
        scalar += area_vector[i] * area_vector[i];
    scalar = sqrt(scalar);
    
    return scalar;
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
double norm(double *vec, int n)
{
	double s = 0.;
	s = dot_product(n, vec, vec);
	return sqrt(s);
}
double distance(int dim, int i, int j, double *point_xyz)
{
	double d = 0.;
	for (int k = 0; k < dim; k++)
		d += point_xyz[dim * i + k] * point_xyz[dim * i + k] + point_xyz[dim * j + k] * point_xyz[dim * j + k];
	for (int k = 0; k < dim; k++)
		d -= 2.0 * point_xyz[dim * i + k] * point_xyz[dim * j + k];
	return sqrt(d);
}