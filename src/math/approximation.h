#ifndef APPROXIMATION_H
#define APPROXIMATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "Peridigm_Constants.hpp"
using namespace Eigen;
namespace APPROXIMATION {

void caller();
void fit_points(const MatrixXd& points, int p, int n_control_points, VectorXd& U, MatrixXd& control_points);
void basis_functions(
    int p, 
    double u, 
    const VectorXd& U, 
    VectorXd& N
    );
void create_approximation(
const int node,
const int nneighbors,
const int* neighborhoodlist,
const double* coordinates,
const int num_control_points,
const int degree,
const bool twoD,
double* AMatrix
);
double basis_func(
    const int i, 
    const int p, 
    const double* U, 
    const double u
    );
int get_field_size(
    const int nnodes,
    const int* neighborhoodList,
    const int num_control_points,
    const bool twoD
    );
double deriv_basis_func(
    const int i, 
    const int p, 
    const double* U, 
    const double u
    );
void knots(
    const int num_control_points,
    const int degree,
    const bool inter,
    double* knots
    );
double get_sample_weighted(
    const double coor,
    const double minVal,
    const double maxVal
    );

void get_approximation(
    const int nnodes,
    const int* neighborhoodList,
    const double* coordinates,
    const int num_control_points,
    const int degree,
    const bool twoD,
    double* AMatrix
    );
void get_jacobians(
    const int nnodes,
    const double* contP,
    const int num_control_points,
    const int degree,
    const bool twoD,
    double* jacobians    
    );
void get_jacobian(
    const int p,
    const int num_control_points,
    const double* contP,
    const double* U,
    const double u,
    const double* V,
    const double v,
    const double* W,
    const double w,
    const bool twoD,
    double* jacobian
    );
void get_control_point(
    const int numNode,
    const int nneighbors,
    const int* neighborhoodlist,
    const int num_control_points,
    const double* coor,
    const double* AMatrix,
    const bool twoD,
    double* contP
    );
void get_gradient(
    const int p,
    const int num_control_points,
    const double* contP,
    const double* U,
    const double u,
    const double* V,
    const double v,
    const double* W,
    const double w,
    const bool twoD,
    MatrixXd& gradientMxM
    );
void get_local_gradient(
    const int p,
    const int num_control_points,
    const double* contP,
    const bool twoD,
    const double* jacobian,
    double* gradient
    );
  
void get_control_points(
    const int numNode,
    const int* neighborhoodlist,
    const int num_control_points,
    const double* coor,
    const double* AMatrix,
    const bool twoD,
    double* contP
    );
}
#endif // APPROXIMATION_H