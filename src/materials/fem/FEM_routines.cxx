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
#include "elastic_correspondence.h"
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
int *topo
)
{
    int count = 0;
    if (twoD){
        for (int jID=0 ; jID<order[1]+1 ; ++jID){
            for (int iID=0 ; iID<order[0]+1 ; ++iID){
                
                topo[3*count] = iID;
                topo[3*count+1] = jID;
                count += 1;
                       
            }
        }

    }
    
    else{
              for (int kID=0 ; kID<order[2]+1 ; ++kID){
            for (int jID=0 ; jID<order[1]+1 ; ++jID){
                for (int iID=0 ; iID<order[0]+1 ; ++iID){
                
                    topo[3*count] = iID;
                    topo[3*count+1] = jID;
                    topo[3*count+2] = kID;
                    count += 1;
        
                }
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

//Eq 9.12 Zienkiewicz
    std::vector<double> xiVector(order+1);
    double* xi = &xiVector[0];
    FEM::defineLagrangianGridSpace(order, xi);
    FEM::shapeFunctionsLagrangeRecursive(Nxi,  order, xi, elCoor);
    FEM::derivativeShapeFunctionsLagrangeRecursive(Bxi,  Nxi,  order, xi, elCoor);

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
const int numIntDir[3]
)
{
    int numInt;
    if  (twoD){
      numInt = numIntDir[0] * numIntDir[1] ;
    }
    else{
      numInt = numIntDir[0] * numIntDir[1] * numIntDir[2] ;
    }
    return numInt;
}

void getNodalForce
(
const double* Bx,
const double* By,
const double* Bz,
const int intPointPtr,
const int nnode,
const double detJ,
const double* Jinv, 
const bool twoD,
const double* sigmaInt, 
double* elNodalForces
)
{   
    double BxiTemp = 0.0, BetaTemp = 0.0, BpsiTemp = 0.0;
 
    for(int nID=0 ; nID<nnode ; ++nID){ 
        BxiTemp  = Bx[intPointPtr + nID];
        BetaTemp = By[intPointPtr + nID];
        if (twoD==false)  {
            BpsiTemp = Bz[intPointPtr + 3*nID];
        }
        // determined by python sympy
        // Bxi_i*J00*s11  + s12*(Beta_i*J00 + Bxi_i*J01) + s23*(Bpsi_i*J00 + Bxi_i*J02)
        // Beta_i*J11*s22 + s12*(Beta_i*J10 + Bxi_i*J11) + s13*(Beta_i*J12 + Bpsi_i*J11)
        // Bpsi_i*J22*s33 + s13*(Beta_i*J22 + Bpsi_i*J21) + s23*(Bpsi_i*J20 + Bxi_i*J22)
        for(int nID=0;nID<nnode;nID++){ 
            elNodalForces[3*nID]   += *(Jinv)   * *(sigmaInt)*  BxiTemp  + *(sigmaInt+1)*(*(Jinv)*BetaTemp   + *(Jinv+1)*BxiTemp)  + *(sigmaInt+5)*(*(Jinv)*BpsiTemp   + *(Jinv+2)*BxiTemp);
            elNodalForces[3*nID+1] += *(Jinv+4) * *(sigmaInt+4)*BetaTemp + *(sigmaInt+1)*(*(Jinv+3)*BetaTemp + *(Jinv+4)*BxiTemp)  + *(sigmaInt+2)*(*(Jinv+4)*BpsiTemp + *(Jinv+5)*BetaTemp);
            elNodalForces[3*nID+2] += *(Jinv+8) * *(sigmaInt+8)*BpsiTemp + *(sigmaInt+2)*(*(Jinv+7)*BpsiTemp + *(Jinv+8)*BetaTemp) + *(sigmaInt+5)*(*(Jinv+6)*BpsiTemp + *(Jinv+8)*BxiTemp);
         }

    }
}    

double getJacobian
(
    const double* Bx,
    const double* By,
    const double* Bz,
    const int nnode, 
    const int intPointPtr,
    const double* coor,
    const bool twoD,
    const double weight,
    double* J, 
    double* Jinv
)
{   
    double detJ;
    int returnCode;
    // EQ. 9.11 - 9.12 The finite element method. The Basis (2000) Zienkievicz, Taylor
    for(int i=0;i<9;i++){
        *(J+i) = 0.0;
    }
    
    if (twoD){
        for(int i=0;i<nnode;i++){
            *(J)   += Bx[intPointPtr + 3*i]*coor[3*i];
            *(J+1) += Bx[intPointPtr + 3*i+1]*coor[3*i+1];
            *(J+2) += 0.0;
            *(J+3) += By[intPointPtr + 3*i]*coor[3*i];
            *(J+4) += By[intPointPtr + 3*i+1]*coor[3*i+1];
            *(J+5) += 0.0;
            *(J+6) += 0.0;
            *(J+7) += 0.0;
            *(J+8) += 0.0;
        }
        returnCode = MATRICES::Invert2by2Matrix(J, detJ, Jinv);
        
        //detJ *= weightsx*weightsy;
    }
    else{
        for(int i=0;i<nnode;i++){
            *(J)   += Bx[intPointPtr + 3*i]*coor[3*i];
            *(J+1) += Bx[intPointPtr + 3*i+1]*coor[3*i+1];
            *(J+2) += Bx[intPointPtr + 3*i+2]*coor[3*i+2];
            *(J+3) += By[intPointPtr + 3*i]*coor[3*i];
            *(J+4) += By[intPointPtr + 3*i+1]*coor[3*i+1];
            *(J+5) += By[intPointPtr + 3*i+2]*coor[3*i+2];
            *(J+6) += Bz[intPointPtr + 3*i]*coor[3*i];            
            *(J+7) += Bz[intPointPtr + 3*i+1]*coor[3*i+1]; 
            *(J+8) += Bz[intPointPtr + 3*i+2]*coor[3*i+2]; 
            
        }
       
        returnCode = MATRICES::Invert3by3Matrix(J, detJ, Jinv);
        //detJ *= weightsx*weightsy*weightsz;
    }

    return detJ*weight;
}
void setElementMatrices
(
    const bool twoD,
    const int offset,
    const int order[3],
    const double* Nxi,
    const double* Neta,
    const double* Npsi,
    const double* Bxi,
    const double* Beta,
    const double* Bpsi,
    double* Bx,
    double* By,
    double* Bz
)
{
  int count = 0;
  if (twoD){

    for(int j=0;j<order[1];j++){
        for(int i=0;i<order[0];i++){
            Bx[offset + count] = Bxi[i]*Neta[j];
            By[offset + count] = Nxi[i]*Beta[j];

            count++;
        }
    }
  }
  else
  {
    for(int k=0;k<order[2];k++){
        for(int j=0;j<order[1];j++){
            for(int i=0;i<order[0];i++){
                Bx[offset + count] = Bxi[i]*Neta[j]*Npsi[k];
                By[offset + count] = Nxi[i]*Beta[j]*Npsi[k];
                Bz[offset + count] = Nxi[i]*Neta[j]*Bpsi[k];
                count++;

            }
        }
    }
  }
}

void setWeights
(
    const int numIntDir[3],
    const bool twoD, 
    const double* weightsx,
    const double* weightsy,
    const double* weightsz,
    double* weights
)
{
    int count = 0;
    if (twoD){
        for(int j=0;j<numIntDir[1];j++){
            for(int i=0;i<numIntDir[0];i++){
               weights[count] = weightsx[i]*weightsy[j];
               count++;
            }
        }

    }   
    else{
        for(int k=0;k<numIntDir[2];k++){
            for(int j=0;j<numIntDir[1];j++){
                for(int i=0;i<numIntDir[0];i++){
                    weights[count] = weightsx[i]*weightsy[j]*weightsz[k];
                    count++;
                }
            }
        }
    }

}
void setGlobalForces
(
    const int nnode,
    const int topoPtr,
    const int* topology,
    const double* elNodalForces,
    const double detJ,
    double* force
)
{
    int localId;
    
    for(int nID=0 ; nID<nnode ; ++nID){
      localId = topology[topoPtr + nID];
      force[3*localId]   += elNodalForces[3*nID]*detJ; 
      force[3*localId+1] += elNodalForces[3*nID+1]*detJ;
      force[3*localId+2] += elNodalForces[3*nID+2]*detJ;
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
const double* Bx,
const double* By,
const double* Bz,
const double* u, 
const int intPointPtr,
const int nnode,
const double* Jinv,
const bool twoD,
double* strain
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
    double BxiTemp = 0.0, BetaTemp = 0.0, BpsiTemp = 0.0;

    for(int i=0 ; i<9 ; ++i)*(strain+i) = 0.0;   
    
    for(int nID=0 ; nID<nnode ; ++nID){ 
        BxiTemp  = Bx[intPointPtr + nID];
        BetaTemp = By[intPointPtr + nID];

        if (twoD==false)  {
            BpsiTemp = Bz[intPointPtr + nID];
 
        }
        *(strain)   += *(Jinv)   * BxiTemp  * u[3*nID];
        *(strain+4) += *(Jinv+4) * BetaTemp * u[3*nID+1];
        *(strain+8) += *(Jinv+8) * BpsiTemp * u[3*nID+2];
        //u1_i*(Beta_i*J00 + Bxi_i*J01) + u2_i*(Beta_i*J10 + Bxi_i*J11)
        *(strain+1) += (*(Jinv) * BetaTemp + *(Jinv+1) * BxiTemp)  * u[3*nID] + (*(Jinv+3) * BetaTemp + *(Jinv+4) * BxiTemp) * u[3*nID+1];
        //u2_i*(Beta_i*J12 + Bpsi_i*J11) + u3_i*(Beta_i*J22 + Bpsi_i*J21)
        *(strain+2) += (*(Jinv+5) * BetaTemp + *(Jinv+4) * BpsiTemp) * u[3*nID+1] + (*(Jinv+7) * BpsiTemp + *(Jinv+8) * BetaTemp)* u[3*nID+2];
        //u1_i*(Bpsi_i*J00 + Bxi_i*J02) + u3_i*(Bpsi_i*J20 + Bxi_i*J22)
        *(strain+5) += (*(Jinv) * BpsiTemp + *(Jinv+2) * BxiTemp)  * u[3*nID] + (*(Jinv+6) * BxiTemp  + *(Jinv+8) * BxiTemp) * u[3*nID+2];
    }
    // 0 1 2
    // 3 4 5
    // 6 7 8
    *(strain+3) = *(strain+1);
    *(strain+6) = *(strain+2);
    *(strain+7) = *(strain+5);
}
int getNnode
(
    const int order[3], 
    const bool twoD
)
{
    int nnode;
    if (twoD)nnode = (order[0]+1)*(order[1]+1);
    else nnode = (order[0]+1)*(order[1]+1)*(order[2]+1);
    return nnode;

}
void getDisplacements
(
    const int numOwnedPoints,
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
