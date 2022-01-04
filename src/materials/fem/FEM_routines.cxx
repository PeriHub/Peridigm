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
#include "FEM_routines.h"
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


void getElementTopo
(
const bool twoD,
const int order[3], 
int topo[][3]
)
{
    int count = 0;
    if (twoD){
        for (int kID=0 ; kID<order[2]+1 ; ++kID){
            for (int jID=0 ; jID<order[1]+1 ; ++jID){
                for (int iID=0 ; iID<order[0]+1 ; ++iID){
                
                topo[count][0] = iID;
                topo[count][1] = jID;
                topo[count][2] = kID;
                count += 1;
        
                }
            }
        }
    }
    else{
        
        for (int jID=0 ; jID<order[1]+1 ; ++jID){
            for (int iID=0 ; iID<order[0]+1 ; ++iID){
                
            topo[count][0] = iID;
            topo[count][1] = jID;
            count += 1;
        
                
            }
        }


    }

}
void getLagrangeElementData
(
const int order, 
const double elCoor,
double* Nxi,
double* Bxi
)
{

// --> EQ (2.25) Willberg Diss; Differential Operator
// B is stored rowise Bxx, Byy, Bzz, Bxz, Bzx, Byz, Bzy, Bxy, Byx
// coor is added later --> issue
// testcase
//Eq 9.12 Zienkiewicz
    std::vector<double> xiVector(order+1);
    double* xi = &xiVector[0];
    //std::vector<double> etaVector(order[1]+1);
    //double* eta = &etaVector[0];
    //std::vector<double> psiVector(order[2]+1);
    //double* psi = &psiVector[0];
    FEM::defineLagrangianGridSpace(order, xi);
    //FEM::defineLagrangianGridSpace(order[1], eta);
    //FEM::defineLagrangianGridSpace(order[2], psi);
    FEM::shapeFunctionsLagrangeRecursive(Nxi,  order, xi, elCoor);
    //FEM::shapeFunctionsLagrangeRecursive(Neta, order[1], eta,elCoory);
    //FEM::shapeFunctionsLagrangeRecursive(Npsi, order[2], psi,elCoorz);
    FEM::derivativeShapeFunctionsLagrangeRecursive(Bxi,  Nxi,  order, xi, elCoor);
    //FEM::derivativeShapeFunctionsLagrangeRecursive(Beta, Neta, order[1], eta,elCoory);
    //FEM::derivativeShapeFunctionsLagrangeRecursive(Bpsi, Npsi, order[2], psi,elCoorz);

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
            if (i!=k){
                B[k] = 1.0/(elCoor-xi[i])*N[k];
            }
        }
    }
}

int getNumberOfIntegrationPoints
(
const bool twoD, 
const int order[3]
)
{
    int numInt;
    if  (twoD){
      numInt = (order[0]+1) * (order[1]+1) ;
    }
    else{
      numInt = (order[0]+1) * (order[1]+1) * (order[2]+1);
    }
    return numInt;
}

void getNodelForce
(
const double* Nxi,
const double* Neta,
const double* Npsi,
const double* Bxi,
const double* Beta,
const double* Bpsi,
const int topo[][3],
const double* sigmaInt, 
const int dof,
const double detJ,
const double* Jinv, 
const bool twoD,
double* elNodalForces
)
{   


    double BxiTemp, BetaTemp, BpsiTemp;
    for(int iID=0 ; iID < dof/3 ; ++iID){
        if (twoD){
            BxiTemp  = Bxi[topo[iID][0]]*Neta[topo[iID][1]] * detJ;
            BetaTemp = Nxi[topo[iID][0]]*Beta[topo[iID][1]] * detJ;
            BpsiTemp = Nxi[topo[iID][0]]*Neta[topo[iID][1]] * detJ;
        }
        else{
            BxiTemp  = Bxi[topo[iID][0]]*Neta[topo[iID][1]]*Npsi[topo[iID][2]] * detJ;
            BetaTemp = Nxi[topo[iID][0]]*Beta[topo[iID][1]]*Npsi[topo[iID][2]] * detJ;
            BpsiTemp = Nxi[topo[iID][0]]*Neta[topo[iID][1]]*Bpsi[topo[iID][2]] * detJ;
        }
        // determined by python sympy
        // Bxi_i*J00*s11 + s12*(Beta_i*J00 + Bxi_i*J01) + s23*(Bpsi_i*J00 + Bxi_i*J02)
        // Beta_i*J11*s22 + s12*(Beta_i*J10 + Bxi_i*J11) + s13*(Beta_i*J12 + Bpsi_i*J11)
        // Bpsi_i*J22*s33 + s13*(Beta_i*J22 + Bpsi_i*J21) + s23*(Bpsi_i*J20 + Bxi_i*J22)
        elNodalForces[3*iID]   = *(Jinv) * *(sigmaInt)*BxiTemp + *(sigmaInt+1)*(*(Jinv)*BetaTemp + *(Jinv+1)*BxiTemp) + *(sigmaInt+5)*(*(Jinv)*BpsiTemp + *(Jinv+2)*BxiTemp);
        elNodalForces[3*iID+1] = *(Jinv+4) * *(sigmaInt+4)*BetaTemp + *(sigmaInt+1)*(*(Jinv+3)*BetaTemp + *(Jinv+4)*BxiTemp) + *(sigmaInt+2)*(*(Jinv+4)*BpsiTemp + *(Jinv+5)*BetaTemp);
        elNodalForces[3*iID+2] = *(Jinv+8) * *(sigmaInt+8)*BpsiTemp + *(sigmaInt+2)*(*(Jinv+7)*BpsiTemp + *(Jinv+8)*BetaTemp) + *(sigmaInt+5)*(*(Jinv+6)*BpsiTemp + *(Jinv+8)*BxiTemp);
    }

}




void getJacobian
(
    const double* Nxi,
    const double* Neta,
    const double* Npsi,
    const double* Bxi,
    const double* Beta,
    const double* Bpsi,
    const int dof, 
    const int topo[][3],
    const double* coor,
    const bool twoD,
    double* J,
    double detJ,
    double* Jinv
)
{
    int returnCode;
    // EQ. 9.11 - 9.12 The finite element method. The Basis (2000) Zienkievicz, Taylor
    if (twoD){
        for(int i=0;i<dof/3;i++){
            *(J)   += Bxi[topo[i][0]]*Neta[topo[i][1]]*coor[3*i];
            *(J+1) += Bxi[topo[i][0]]*Neta[topo[i][1]]*coor[3*i];
            *(J+2) += 0.0;
            *(J+3) += Nxi[topo[i][0]]*Beta[topo[i][1]]*coor[3*i+1];
            *(J+4) += Nxi[topo[i][0]]*Beta[topo[i][1]]*coor[3*i+1];
            *(J+5) += 0.0;
            *(J+6) += 0.0;
            *(J+7) += 0.0;
            *(J+8) += 0.0;
        }
        returnCode = MATRICES::Invert2by2Matrix(J, detJ, Jinv);
    }
    else{
        for(int i=0;i<dof/3;i++){
            *(J)   += Bxi[topo[i][0]]*Neta[topo[i][1]]*Npsi[topo[i][2]]*coor[3*i];
            *(J+1) += Bxi[topo[i][0]]*Neta[topo[i][1]]*Npsi[topo[i][2]]*coor[3*i];
            *(J+2) += Bxi[topo[i][0]]*Neta[topo[i][1]]*Npsi[topo[i][2]]*coor[3*i];
            *(J+3) += Nxi[topo[i][0]]*Beta[topo[i][1]]*Npsi[topo[i][2]]*coor[3*i+1];
            *(J+4) += Nxi[topo[i][0]]*Beta[topo[i][1]]*Npsi[topo[i][2]]*coor[3*i+1];
            *(J+5) += Nxi[topo[i][0]]*Beta[topo[i][1]]*Npsi[topo[i][2]]*coor[3*i+1];
            *(J+6) += Nxi[topo[i][0]]*Neta[topo[i][1]]*Bpsi[topo[i][2]]*coor[3*i+2];
            *(J+7) += Nxi[topo[i][0]]*Neta[topo[i][1]]*Bpsi[topo[i][2]]*coor[3*i+2];
            *(J+8) += Nxi[topo[i][0]]*Neta[topo[i][1]]*Bpsi[topo[i][2]]*coor[3*i+2];
        }
       
        returnCode = MATRICES::Invert3by3Matrix(J, detJ, Jinv);
    }
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
            if (i!=k){
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
const double* Nxi,
const double* Neta,
const double* Npsi,
const double* Bxi,
const double* Beta,
const double* Bpsi,
const int topo[][3],
const double* u, 
const int dof,
const double* Jinv,
const bool twoD,
double strain[3][3]
)
// zienkiewicz Basics EQ 6.11 with different Voigt notation
// determined with python sympy
//J00*(Bxi_i*u1_i)
//J11*(Beta_i*u2_i)
//J22*(Bpsi_i*u3_i)
//u1_i*(Bpsi_i*J00 + Bxi_i*J02) + u3_i*(Bpsi_i*J20 + Bxi_i*J22)
//u2_i*(Beta_i*J12 + Bpsi_i*J11) + u3_i*(Beta_i*J22 + Bpsi_i*J21)
//u1_i*(Beta_i*J00 + Bxi_i*J01) + u2_i*(Beta_i*J10 + Bxi_i*J11)

{
    double BxiTemp, BetaTemp, BpsiTemp;
    for (int iID = 0; iID < 6; ++iID){
        strain[iID][0] = 0.0;
        strain[iID][1] = 0.0;
        strain[iID][2] = 0.0;
        strain[iID][3] = 0.0;
        strain[iID][4] = 0.0;
        strain[iID][5] = 0.0;
    }
    for(int iID=0 ; iID < dof/3 ; ++iID){
        if (twoD){
            BxiTemp  = Bxi[topo[iID][0]]*Neta[topo[iID][1]];
            BetaTemp = Nxi[topo[iID][0]]*Beta[topo[iID][1]];
            BpsiTemp = Nxi[topo[iID][0]]*Neta[topo[iID][1]];
        }
        else{
            BxiTemp  = Bxi[topo[iID][0]]*Neta[topo[iID][1]]*Npsi[topo[iID][2]];
            BetaTemp = Nxi[topo[iID][0]]*Beta[topo[iID][1]]*Npsi[topo[iID][2]];
            BpsiTemp = Nxi[topo[iID][0]]*Neta[topo[iID][1]]*Bpsi[topo[iID][2]];
        }
        // determined with python sympy
        //J00*(Bxi_i*u1_i)
        //J11*(Beta_i*u2_i)
        //J22*(Bpsi_i*u3_i)
        
        
        
        strain[0][0] += *(Jinv)   * BxiTemp  * u[3*iID];
        strain[1][1] += *(Jinv+4) * BetaTemp * u[3*iID+1];
        strain[2][2] += *(Jinv+8) * BpsiTemp * u[3*iID+2];
        //u1_i*(Bpsi_i*J00 + Bxi_i*J02) + u3_i*(Bpsi_i*J20 + Bxi_i*J22)
        
        strain[1][2] += (*(Jinv) * BpsiTemp + *(Jinv+2) * BxiTemp)  * u[3*iID] + (*(Jinv+6) * BxiTemp  + *(Jinv+8) * BxiTemp) * u[3*iID+2];
        //u2_i*(Beta_i*J12 + Bpsi_i*J11) + u3_i*(Beta_i*J22 + Bpsi_i*J21)
        strain[0][2] += (*(Jinv+5) * BetaTemp + *(Jinv+4) * BpsiTemp) * u[3*iID+1] + (*(Jinv+7) * BpsiTemp + *(Jinv+8) * BetaTemp)* u[3*iID+2];
        //u1_i*(Beta_i*J00 + Bxi_i*J01) + u2_i*(Beta_i*J10 + Bxi_i*J11)
        strain[0][1] += (*(Jinv) * BetaTemp + *(Jinv+1) * BxiTemp)  * u[3*iID] + (*(Jinv+3) * BetaTemp + *(Jinv+4) * BxiTemp) * u[3*iID+1];
      
        
    }
    strain[1][0] = strain[0][1];
    strain[2][0] = strain[0][2];
    strain[2][1] = strain[1][2];
}

void getDisplacements
(
    int numOwnedPoints,
    const double* modelCoordinates,
    const double* coordinatesNP1,
    double* displacements
)
{
    const double* modelCoord = modelCoordinates;
    const double* coorNP1 = coordinatesNP1;
    double* disp = displacements;
    for(int i=0 ; i<numOwnedPoints ; ++i){
        *(disp+i) = *(coorNP1+i)-*(modelCoord+i);
    }
}

}
