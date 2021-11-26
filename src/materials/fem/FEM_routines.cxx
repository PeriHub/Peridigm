//! \file FEM_routines.cxx

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//
// funded by dfg project reference
// Licence agreement
//
// Christian Willberg    christian.willberg@dlr.de
//@HEADER

#include "matrices.h"
#include <Sacado.hpp>
#include <math.h>
#include <cmath> 

//#include <Teuchos_Assert.hpp>
//#include <Epetra_SerialComm.h>


using namespace std;
namespace FEM {

void weightsAndIntegrationPoints
(
const int order, 
double* elCoor,
double* weights
)
{
//https://de.wikipedia.org/wiki/Gau%C3%9F-Quadratur
if (order == 0)
{
    elCoor[0] = 0;
    weights[0] = 2;
}
if (order == 1)
{
    elCoor[0] = -sqrt(1.0/3.0);
    weights[0] = 1;
    elCoor[1] =  sqrt(1.0/3.0);
    weights[1] = 1;
   }
else if (order==2)
{
    elCoor[0] = -sqrt(3.0/5.0);
    weights[0] = 5.0/9.0;
    elCoor[1] =  0.0;
    weights[1] = sqrt(1.0/3.0);
    elCoor[2] =  sqrt(3.0/5.0);
    weights[2] = 5.0/9.0;
}
//else if (order == 3)
//{}
else
{
    //hier muss eine Fehlermeldung rein
}

}



void getLagrangeElementData
(
const int order[3], 
const double* elCoorx,
const double* elCoory,
const double* elCoorz,
double* Nxi,
double* Neta,
double* Npsi,
double* Bxi,
double* Beta,
double* Bpsi
)
{

// --> EQ (2.25) Willberg Diss; Differential Operator
// B is stored rowise Bxx, Byy, Bzz, Bxz, Bzx, Byz, Bzy, Bxy, Byx
// coor is added later --> issue
// testcase
//Eq 9.12 Zienkiewicz
    std::vector<double> xiVector(order[0]+1), etaVector(order[1]+1), psiVector(order[2]+1);
    double* xi = &xiVector[0], eta = &xiVector[0], psi = &xiVector[0];
    FEM::defineLagrangianGridSpace(order[0], xi);
    FEM::defineLagrangianGridSpace(order[1], eta);
    FEM::defineLagrangianGridSpace(order[2], psi);
    FEM::shapeFunctionsLagrangeRecursive(Nxi,  order[0],xi, elCoorx);
    FEM::shapeFunctionsLagrangeRecursive(Neta, order[1],eta,elCoory);
    FEM::shapeFunctionsLagrangeRecursive(Npsi, order[2],psi,elCoorz);
    FEM::derivativeShapeFunctionsLagrangeRecursive(Bxi,  Nxi,  order[0], xi, elCoorx);
    FEM::derivativeShapeFunctionsLagrangeRecursive(Beta, Neta, order[1], eta,elCoory);
    FEM::derivativeShapeFunctionsLagrangeRecursive(Bpsi, Npsi, order[2], psi,elCoorz);

}

void createBMatrix
(
const int order[3], 
const double elCoor[3]
double* Nmatrix,
double* Bmatrix
)
{
    std::vector<double> xiVector(order[0]+1), etaVector(order[1]+1), psiVector(order[2]+1);
    double* xi = &xiVector[0], eta = &xiVector[0], psi = &xiVector[0];
    std::vector<double> NxiVector(order[0]+1), NetaVector(order[1]+1), NpsiVector(order[2]+1);
    double* Nxi = &xiVector[0], Neta = &xiVector[0], Npsi = &xiVector[0];
    std::vector<double> BxiVector(order[0]+1), BetaVector(order[1]+1), BpsiVector(order[2]+1);
    double* Bxi = &xiVector[0], Beta = &xiVector[0], Bpsi = &xiVector[0];
    FEM::getLagrangeElementData()
}

void shapeFunctionsLagrangeRecursive
    (
        double* N, 
        const int order, 
        const double* xi,
        const double elCoor
    )

}

void defineLagrangianGridSpace
(
const int order,
double* xi
)
{
    double len = 2.0 / order;
    for(int i=0;i<order+1;i++){
        xi[i] = -1 + i*len;
    }

}


void derivativeShapeFunctionsLagrangeRecursive
(
    double* B, 
    const double* N,
    const int order,
    const double* xi,
    const double elCoor
)
{
// https://math.stackexchange.com/questions/809927/first-derivative-of-lagrange-polynomial

    for(int k=0;k<order+1;k++){
        for(int i=0;i<order+1;i++){
            if (i/=k){
                B[k] = 1.0/(elCoor-xi[i])*N[k];
            }
        }
}

void Jacobian
(
    const double* Nxi,
    const double* Neta,
    const double* Npsi,
    const double* Bxi,
    const double* Beta,
    const double* Bpsi,
    const int nELnodes, 
    const double coor,
    double J[3][3],
    double Jinv[3][3],
    double detJ
)
{// EQ. 9.11 - 9.12 The finite element method. The Basis (2000) Zienkievicz, Taylor
    for(int i=0;k<nELnodes+1;i++){
        J[0][0] += Bxi[i]*Neta[i]*Npsi[i]*coor[3*i];
        J[0][1] += Bxi[i]*Neta[i]*Npsi[i]*coor[3*i];
        J[0][2] += Bxi[i]*Neta[i]*Npsi[i]*coor[3*i];
        J[1][0] += Nxi[i]*Beta[i]*Npsi[i]*coor[3*i+1];
        J[1][1] += Nxi[i]*Beta[i]*Npsi[i]*coor[3*i+1];
        J[1][2] += Nxi[i]*Beta[i]*Npsi[i]*coor[3*i+1];
        J[2][0] += Nxi[i]*Neta[i]*Bpsi[i]*coor[3*i+2];
        J[2][1] += Nxi[i]*Neta[i]*Bpsi[i]*coor[3*i+2];
        J[2][2] += Nxi[i]*Neta[i]*Bpsi[i]*coor[3*i+2];
    }
        
    MATRICES::Invert3by3Matrix(J, detJ, Jinv)
}

void shapeFunctionsLagrangeRecursive
    (
        double* N, 
        const int order, 
        const double* xi,
        const double elCoor
    )
{

    for(int k=0;k<order+1;k++){
        N[k] = 1;
        for(int i=0;i<order+1;i++){
            if (i/=k){
                N[k] = N[k] * (elCoor-xi[i]) / (xi[k]-xi[i]);
            }
        }
    }
}



void shapeFunctionsLagrange
(
double* Nmatrix, 
const int order[3], 
const double elCoor[3]
)
{
    if (order[0] == 1){
        Nmatrix[0] = (1-elCoor[0])*(1-elCoor[1])*(1-elCoor[2])/8.0;
        Nmatrix[1] = (1+elCoor[0])*(1-elCoor[1])*(1-elCoor[2])/8.0;
        Nmatrix[2] = (1+elCoor[0])*(1+elCoor[1])*(1-elCoor[2])/8.0;
        Nmatrix[3] = (1-elCoor[0])*(1+elCoor[1])*(1+elCoor[2])/8.0;
        Nmatrix[4] = (1-elCoor[0])*(1-elCoor[1])*(1+elCoor[2])/8.0;
        Nmatrix[5] = (1+elCoor[0])*(1-elCoor[1])*(1+elCoor[2])/8.0;
        Nmatrix[6] = (1+elCoor[0])*(1+elCoor[1])*(1+elCoor[2])/8.0;
        Nmatrix[7] = (1-elCoor[0])*(1+elCoor[1])*(1+elCoor[2])/8.0;
    }

}
void computeStrain
(
const double B[6][],
const double* u, 
const int dof,
double strain[6]
)
{
    for (int iID = 0; iID < 6; ++iID){
        strain[iID] = 0.0;
        for(int jID=0 ; jID < dof ; ++jID){
            strain[iID] += B[iID][jID] * u[jID];
        }
    }
}


}
