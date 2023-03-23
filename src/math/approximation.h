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
    const int num_control_points
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
    const double maxVal,
    const double minVal
    );
}
#endif // APPROXIMATION_H