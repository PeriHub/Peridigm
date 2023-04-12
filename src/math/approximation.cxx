#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include "approximation.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace std;
using namespace Eigen;
namespace APPROXIMATION {


double deriv_basis_func(
    const int i, 
    const int p, 
    const double* U, 
    const double u
    )
    {
       /*
        Computes the ith NURBS basis function of degree p on the knot vector U at the point u.
        
        Args:
        i (int): Index of the basis function to compute.
        p (int): Degree of the basis functions.
        U (ndarray): Knot vector.
        u (float): Point to evaluate the basis function at.
        
        Returns:
        float: The value of the ith NURBS basis function at u.
        */
                    
       // dev comment
       // function.cxx to include all functions
       //
        double denom1 = U[i+p] - U[i];
        double denom2 = U[i+p+1] - U[i+1];
        double term1, term2;

        if (denom1 == 0.0) term1 = 0.0;
        else term1 = p/denom1 * basis_func(i, p-1, U, u);
        if (denom2 == 0.0) term2 = 0.0;
        else term2 = p/denom2 * basis_func(i+1, p-1, U, u);
        
        return term1 - term2;    

    }

double  basis_func(
    const int i, 
    const int p, 
    const double* U, 
    const double u
    )
    {
        /*
        Computes the ith NURBS basis function of degree p on the knot vector U at the point u.
        
        Args:
        i (int): Index of the basis function to compute.
        p (int): Degree of the basis functions.
        U (ndarray): Knot vector.
        u (float): Point to evaluate the basis function at.
        
        Returns:
        double: The value of the ith NURBS basis function at u.
        */
        if (p == 0){
            
            // work around -> to be checked if it works for inner cruve interpolation
            // self.maxP assumes right now equal p degrees
            if (U[i] <= u && u < U[i+1] && i<=p)return 1.0;
            if (U[i] < u && u <= U[i+1] && i>=p)return 1.0;
            else return 0.0;
        }
            
        double denom1 = U[i+p] - U[i];
        double denom2 = U[i+p+1] - U[i+1];
        double numer1 = u - U[i];
        double numer2 = U[i+p+1] - u;
        
        double term1;
        double term2;
        if (denom1 == 0.0){term1 = 0.0;}
        else{term1 = numer1/denom1 * basis_func(i, p-1, U, u);}
        if (denom2 == 0.0){term2 = 0.0;} else {term2  = numer2/denom2 * basis_func(i+1, p-1, U, u);}
        
        return term1 + term2;
    }


void knots(
    const int num_control_points,
    const int degree,
    const bool inter,
    double* knots
    )
    {
        if (num_control_points<degree) std::cout<<"you need more control points; at least p+1"<<std::endl;
        int num_knots = int(num_control_points + degree + 1);
        
        // Generate uniform knots
        if (inter){
            int n = 0;
            for(int i=0 ; i<degree ; ++i,++n)knots[n] = 0; 
            for(int i=0 ; i<num_knots-2*degree-1 ; ++i,++n)knots[n] = float(i) / (num_knots - 2*degree - 1); 
            for(int i=0 ; i<degree+1 ; ++i,++n)knots[n] = 1; 
        }
        else{
            for(int i=0 ; i<num_knots ; ++i){
                knots[i] = float(i) /  (num_knots - 1); 
            }
        }
    }


void get_approximation(
    const int nnodes,
    const int* neighborhoodList,
    const double* coordinates,
    const int num_control_points,
    const int degree,
    const bool twoD,
    double* AMatrix
)
{
    const int *neighborListPtr = neighborhoodList;
    int nsquare = num_control_points * num_control_points;
    double* AMatrixLocal;

    for(int iID=0 ; iID<nnodes ; ++iID){
        int numNeighbors = *neighborListPtr; neighborListPtr++;
        AMatrix += numNeighbors * nsquare;
        create_approximation(iID, numNeighbors, neighborListPtr,coordinates,num_control_points,degree,twoD,AMatrix);
      
    }
}


void create_approximation(
    const int node,
    const int nneighbors,
    const int* neighborhoodlist,
    const double* coordinates,
    const int num_control_points,
    const int degree,
    const bool twoD,
    double* AMatrix
)
{
    // vectoren müssen klar sein. Da das alles in der Init passiert ist es unkritisch
    const int dof = PeridigmNS::dof();
    const int dof2D = 2;
    int neighborID;
    int nsquare = num_control_points * num_control_points;
    std::vector<double> UVector(num_control_points + degree + 1);
    double* U = &UVector[0];
    std::vector<double> VVector(num_control_points + degree + 1);
    double* V = &VVector[0];
    Eigen::MatrixXd A(nsquare,nneighbors+1);
    Eigen::MatrixXd AAT(nsquare,nsquare);
    Eigen::MatrixXd approx(nsquare,nneighbors+1);
    std::vector<double> PVector(dof2D*(nneighbors+1));
    double* P = &PVector[0];
    double u, v;
    
    APPROXIMATION::knots(num_control_points,degree,true,U);
    APPROXIMATION::knots(num_control_points,degree,true,V);
    
    double minVal[dof];
    double maxVal[dof];
    // create a list of all points within the horizon of point A including point A
    for(int i=0 ; i<dof2D ; ++i){
        P[i] = coordinates[node * dof + i];
        minVal[i] = coordinates[node * dof + i];
        maxVal[i] = coordinates[node * dof + i];
    }
    for(int i=1 ; i<nneighbors+1 ; ++i){   
        neighborID = neighborhoodlist[i - 1];    
        for(int j=0 ; j<dof2D ; ++j){
            P[i * dof2D + j] = coordinates[neighborID * dof + j];
            if (maxVal[j]<P[i * dof2D + j])maxVal[j]=P[i * dof2D + j];
            if (minVal[j]>P[i * dof2D + j])minVal[j]=P[i * dof2D + j];
        }
    } 

    for(int pos=0 ; pos<nneighbors+1 ; ++pos){  
        for(int i=0 ; i<num_control_points ; ++i){       
            for(int j=0 ; j<num_control_points ; ++j){            
                u = APPROXIMATION::get_sample_weighted(P[pos*dof2D],minVal[0],maxVal[0]);
                v = APPROXIMATION::get_sample_weighted(P[pos*dof2D+1],minVal[1],maxVal[1]);
                A(i*num_control_points+j,pos) = basis_func(i,degree,U,u)*basis_func(j,degree,V,v);    
          
            }
        }
    }
    AAT = A*A.transpose();
  
    approx = AAT.inverse()*A;

    for(int i=0 ; i<nsquare ; ++i){    
        for(int j=0 ; j<nneighbors+1; ++j){ 
            AMatrix[i*(nneighbors+1) + j] = approx(i,j);
        }
    }

}

double get_sample_weighted(
    const double coor,
    const double minVal,
    const double maxVal
    )
    {
       return (coor-minVal)/(maxVal - minVal);
    }


int get_field_size(
    const int nnodes,
    const int* neighborhoodList,
    const int num_control_points
    )
    {
       
       const int *neighborListPtr = neighborhoodList;
       int val = 0;
       for(int iID=0 ; iID<nnodes ; ++iID){
        int numNeighbors = *neighborListPtr; neighborListPtr+=numNeighbors+1;
        val += numNeighbors + 1 + num_control_points; // +1 because of the iID node which is not in the neighborhoodlist
    }
    return val;
    }
}