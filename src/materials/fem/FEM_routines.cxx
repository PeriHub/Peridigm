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
#include <Teuchos_Assert.hpp>
//#include <Epetra_SerialComm.h>


using namespace std;
namespace FEM {

bool weightsAndIntegrationPoints
(
const int order, 
double* elCoor,
double* weights
)
{
//https://de.wikipedia.org/wiki/Gau%C3%9F-Quadratur
bool success = true;
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
    weights[1] = 8.0/9.0;
    elCoor[2] =  sqrt(3.0/5.0);
    weights[2] = 5.0/9.0;
}
else
{
    success = false;
}
    return success;
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
    FEM::derivativeShapeFunctionsLagrangeRecursive(Bxi, order, xi, elCoor);

}
std::vector<int> getTopology
(
const int numOwnedPoints,
const int* ownedIDs,
const int* neighborhoodList,
const double* nodeType,
double* detachedNodes
)
{
  std::vector<int> topology;
  topology.push_back(0);
  int topoPtr = 0;
  int nodeId;
  int numberOfFeNodes;
  int numberOfPdNodes;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    nodeId = ownedIDs[iID];
    int numNodes = neighborhoodList[topoPtr++];

    detachedNodes[nodeId]=1;

    if (nodeType[nodeId]==2){
        numberOfFeNodes = 4;
        topology[0] += 1; // Number of elements
        topology.push_back(nodeId); // elementID
        topology.push_back(numberOfFeNodes);
        numberOfPdNodes = numNodes - numberOfFeNodes;
        for(int j=0 ; j<numberOfFeNodes ; ++j){
            topology.push_back(neighborhoodList[topoPtr++]);
        }
        topology.push_back(numberOfPdNodes);
        for(int j=0 ; j<numberOfPdNodes ; ++j){
            topology.push_back(neighborhoodList[topoPtr++]);
        }
    }
    else{topoPtr+=numNodes;}
  }
  return topology;
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
    const int order,
    const double* xi,
    const double elCoor
)
{
// https://en.wikipedia.org/wiki/Lagrange_polynomial#Derivation[6]
// sympy calculated

    for(int k=0;k<order+1;k++){
        B[k] = 0.0;
        for(int i=0;i<order+1;i++){
            if (i!=k){
                double temp = 1 / (xi[i]-xi[k]);
                for(int m=0;m<order+1;m++){
                    if (m!=i && m!=k){
                        temp *= (elCoor-xi[m]) / (xi[k] - xi[m]);
                    }
                }
                B[k] += temp;
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

void nodalForce
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
double* nforces
)
{   
    double BxiTemp = 0.0, BetaTemp = 0.0, BpsiTemp = 0.0;
    // int elId;
    for(int nID=0 ; nID<nnode ; ++nID){ 
        BxiTemp  = Bx[intPointPtr + nID];
        BetaTemp = By[intPointPtr + nID];
        if (twoD==false)  {
            BpsiTemp = Bz[intPointPtr + nID];
        }   
        nforces[3*nID]  -=( *(sigmaInt)  *(*(Jinv)*BxiTemp + *(Jinv+1)*BetaTemp + *(Jinv+2)*BpsiTemp) + *(sigmaInt+1)*(*(Jinv+3)*BxiTemp + *(Jinv+4)*BetaTemp + *(Jinv+5)*BpsiTemp) + *(sigmaInt+2)*(*(Jinv+6)*BxiTemp + *(Jinv+7)*BetaTemp + *(Jinv+8)*BpsiTemp) )*detJ;
        nforces[3*nID+1]-=( *(sigmaInt+1)*(*(Jinv)*BxiTemp + *(Jinv+1)*BetaTemp + *(Jinv+2)*BpsiTemp) + *(sigmaInt+4)*(*(Jinv+3)*BxiTemp + *(Jinv+4)*BetaTemp + *(Jinv+5)*BpsiTemp) + *(sigmaInt+5)*(*(Jinv+6)*BxiTemp + *(Jinv+7)*BetaTemp + *(Jinv+8)*BpsiTemp) )*detJ;
        nforces[3*nID+2]-=( *(sigmaInt+2)*(*(Jinv)*BxiTemp + *(Jinv+1)*BetaTemp + *(Jinv+2)*BpsiTemp) + *(sigmaInt+5)*(*(Jinv+3)*BxiTemp + *(Jinv+4)*BetaTemp + *(Jinv+5)*BpsiTemp) + *(sigmaInt+8)*(*(Jinv+6)*BxiTemp + *(Jinv+7)*BetaTemp + *(Jinv+8)*BpsiTemp) )*detJ;
        
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
    const double weight,
    const bool twoD,
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
            *(J)   += Bx[intPointPtr + i]*coor[3*i];
            *(J+1) += Bx[intPointPtr + i]*coor[3*i+1];
            *(J+2) += 0.0;
            *(J+3) += By[intPointPtr + i]*coor[3*i];
            *(J+4) += By[intPointPtr + i]*coor[3*i+1];
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
            *(J)   += Bx[intPointPtr + i]*coor[3*i];
            *(J+1) += Bx[intPointPtr + i]*coor[3*i+1];
            *(J+2) += Bx[intPointPtr + i]*coor[3*i+2];
            *(J+3) += By[intPointPtr + i]*coor[3*i];
            *(J+4) += By[intPointPtr + i]*coor[3*i+1];
            *(J+5) += By[intPointPtr + i]*coor[3*i+2];
            *(J+6) += Bz[intPointPtr + i]*coor[3*i];            
            *(J+7) += Bz[intPointPtr + i]*coor[3*i+1]; 
            *(J+8) += Bz[intPointPtr + i]*coor[3*i+2]; 
            
        }
       
        returnCode = MATRICES::Invert3by3Matrix(J, detJ, Jinv);
        //detJ *= weightsx*weightsy*weightsz;
    }
    
    TEUCHOS_TEST_FOR_EXCEPT_MSG(returnCode>0, "**** Jacobi Matrix is not invertable.\n");
  
    
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

        for(int j=0;j<order[1]+1;j++){
            for(int i=0;i<order[0]+1;i++){
                Bx[offset + count]  = Bxi[i]*Neta[j];
                By[offset + count] = Nxi[i]*Beta[j];

                count++;
            }
        }
    }
else
    {
        for(int k=0;k<order[2]+1;k++){
            for(int j=0;j<order[1]+1;j++){
                for(int i=0;i<order[0]+1;i++){
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
void setCouplingStiffnessMatrix
(
    const double* undeformedCoor,
    int topoPtr,
    const int* topology,
    double* couplingStiffnessMatrix
)
{
    double kappa = 10e4;
    int globalId;
    int elementID;
    int numElemNodes, numPdNodes;

    for(int iID=0 ; iID<topology[0] ; ++iID){
        // for averaging the element number to which the node is connected has to be known
        elementID = topology[topoPtr++];
        numElemNodes = topology[topoPtr++];

        double x[4];
        double y[4];
        double z[4];

        for(int nID=0 ; nID<numElemNodes ; nID++, topoPtr++){
            globalId = topology[topoPtr];

            // Extract coordinates of ith element nodes
            x[nID] = undeformedCoor[3*globalId+0];
            y[nID] = undeformedCoor[3*globalId+1];
            z[nID] = undeformedCoor[3*globalId+2];

        }

        numPdNodes = topology[topoPtr++];

        double* coord = new double[numPdNodes * 3];

        for(int nID=0 ; nID<numPdNodes ; nID++, topoPtr++){
            globalId = topology[topoPtr];
            for(int j = 0; j < 3; j++){
                coord[nID * 3 + j] =undeformedCoor[3*globalId+j];
            }
        }

        // Compute distances from point 1 to the other three points
        double dist13 = sqrt(pow((x[2]-x[0]),2) + pow((y[2]-y[0]),2));
        double dist12 = sqrt(pow((x[1]-x[0]),2) + pow((y[1]-y[0]),2));
        double dist14 = sqrt(pow((x[3]-x[0]),2) + pow((y[3]-y[0]),2));
        
        // Find the farthest point
        double max_dist = std::max({dist13, dist12, dist14});
        double min_dist = std::min({dist13, dist12, dist14});
        
        // std::cout << coord[0][0] <<std::endl;
        // std::cout << coord[0][1] <<std::endl;
        // std::cout << coord[0][2] <<std::endl;
        std::cout << x[0] <<std::endl;
        std::cout << x[1] <<std::endl;
        std::cout << x[2] <<std::endl;
        std::cout << x[3] <<std::endl;
        std::cout << y[0] <<std::endl;
        std::cout << y[1] <<std::endl;
        std::cout << y[2] <<std::endl;
        std::cout << y[3] <<std::endl;
        std::cout << max_dist <<std::endl;
        std::cout << min_dist <<std::endl;

        double ksi = 0, eta = 0;
        
        int ipd = 0;

        // Compute ksi and eta accordingly
        if (max_dist == dist13) {
            ksi = 2 * (coord[ipd * 3 + 0] - (x[2]+x[0])/2) / min_dist;
            eta = 2 * (coord[ipd * 3 + 1] - (y[2]+y[0])/2) / min_dist;
        }
        else if (max_dist == dist12) {
            ksi = 2 * (coord[ipd * 3 + 0] - (x[1]+x[0])/2) / min_dist;
            eta = 2 * (coord[ipd * 3 + 1] - (y[1]+y[0])/2) / min_dist;
        }
        else if (max_dist == dist14) {
            ksi = 2 * (coord[ipd * 3 + 0] - (x[3]+x[0])/2) / min_dist;
            eta = 2 * (coord[ipd * 3 + 1] - (y[3]+y[0])/2) / min_dist;
        }

        std::cout << ksi <<std::endl;
        std::cout << eta <<std::endl;

        // Define shape functions for PD point in local coordinates
        double N1p = 0.25 * (1-ksi) * (1-eta);
        double N2p = 0.25 * (1+ksi) * (1-eta);
        double N3p = 0.25 * (1+ksi) * (1+eta);
        double N4p = 0.25 * (1-ksi) * (1+eta);
        
        // Create matrix of shape function values at PD point
        double Np[1][4] = {{N1p, N2p, N3p, N4p}};
        
        // Define coupling stiffness matrix
        int I = 1;

        double Np0 = Np[0][0];

        couplingStiffnessMatrix[0] = kappa * I;
        couplingStiffnessMatrix[1] = kappa * -Np0;
        couplingStiffnessMatrix[2] = kappa * -Np[0][1];
        couplingStiffnessMatrix[3] = kappa * -Np[0][2];
        couplingStiffnessMatrix[4] = kappa * -Np[0][3];
        couplingStiffnessMatrix[5] = kappa * -Np0;
        couplingStiffnessMatrix[6] = kappa * pow(Np0, 2);
        couplingStiffnessMatrix[7] = kappa * Np0 * Np[0][1];
        couplingStiffnessMatrix[8] = kappa * Np0 * Np[0][2];
        couplingStiffnessMatrix[9] = kappa * Np0 * Np[0][3];
        couplingStiffnessMatrix[10] = kappa * -Np[0][1];
        couplingStiffnessMatrix[11] = kappa * Np0 * Np[0][1];
        couplingStiffnessMatrix[12] = kappa * pow(Np[0][1], 2);
        couplingStiffnessMatrix[13] = kappa * Np[0][1] * Np[0][2];
        couplingStiffnessMatrix[14] = kappa * Np[0][1] * Np[0][3];
        couplingStiffnessMatrix[15] = kappa * -Np[0][2];
        couplingStiffnessMatrix[16] = kappa * Np0 * Np[0][2];
        couplingStiffnessMatrix[17] = kappa * Np[0][1] * Np[0][2];
        couplingStiffnessMatrix[18] = kappa * Np[0][2] * Np[0][2];
        couplingStiffnessMatrix[19] = kappa * Np[0][2] * Np[0][3];
        couplingStiffnessMatrix[20] = kappa * -Np[0][3];
        couplingStiffnessMatrix[21] = kappa * Np0 * Np[0][3];
        couplingStiffnessMatrix[22] = kappa * Np[0][1] * Np[0][3];
        couplingStiffnessMatrix[23] = kappa * Np[0][2] * Np[0][3];
        couplingStiffnessMatrix[24] = kappa * Np[0][3] * Np[0][3];

        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                std::cout << *(couplingStiffnessMatrix + (i*5+j))<<std::endl;
            }
        }

    }
}
void setGlobalStresses
(
    const int elementID,
    const double* sigmaInt,
    double* sigmaNP1
)

{
    for (int i=0 ; i<9 ; ++i){
       sigmaNP1[9 * elementID + i] += *(sigmaInt+i) ;      
    }
}
void setToZero(
    double* A,
    int len
)
{
    MATRICES::setToZero(A, len);
}

void setNodalStresses
(
    const int nnode,
    const int elementID,
    int topoPtr,
    const int* topology,
    double* sigmaNP1
)

{
    int globalId;
    for(int nID=0 ; nID<nnode ; ++nID){
        globalId = topology[topoPtr + nID];
        for (int i=0 ; i<9 ; ++i){
            // set the element stresses equal to the element nodel stress
            // to have stresses at the nodes for visualisation
            sigmaNP1[9 * globalId + i] = sigmaNP1[9 * elementID + i];
            
        }
        
    }
}
void setGlobalForces
(
    const int nnode,
    const int elementID,
    int topoPtr,
    const int* topology,
    const double* elNodalForces,
    const double* volume, //needed for PD solver
    double* force
)
{
    int globalId;

    setToZero(&force[3 * elementID], 3);

    for(int nID=0 ; nID<nnode ; ++nID){
      globalId = topology[topoPtr + nID];
      for (int i=0 ; i<3 ; ++i){  
          force[3*globalId+i]      += elNodalForces[3*nID+i]  * volume[globalId];
      }
     }

}
void setGlobalCouplingForces
(
    const int nnode,
    const int elementID,
    int topoPtr,
    const int* topology,
    const double* elNodalForces,
    const double* volume, //needed for PD solver
    const double* deformedCoor,
    const double* undeformedCoor,
    const double* couplingStiffnessMatrix,
    double* force
)
{
    double ux[5], uy[5], uz[5];

    int globalId;

    for(int nID=0 ; nID<nnode ; ++nID){
        globalId = topology[topoPtr + nID];
        ux[nID + 1] = deformedCoor[3 * globalId + 0] - undeformedCoor[3 * globalId + 0];
        uy[nID + 1] = deformedCoor[3 * globalId + 1] - undeformedCoor[3 * globalId + 1];
        uz[nID + 1] = deformedCoor[3 * globalId + 2] - undeformedCoor[3 * globalId + 2];
    }
    
    //pd nodes
    globalId = topology[topoPtr + nnode];
    ux[0] = deformedCoor[3 * globalId + 0] - undeformedCoor[3 * globalId + 0];
    uy[0] = deformedCoor[3 * globalId + 1] - undeformedCoor[3 * globalId + 1];
    uz[0] = deformedCoor[3 * globalId + 2] - undeformedCoor[3 * globalId + 2];

    double* couplingForces = new double[5 * 3];

    for (int i=0;i<5*3;i++){
         couplingForces[i]=0;
    }
    for (int i=0;i<5;i++){
        for (int j=0;j<5;j++){
            couplingForces[i]+=( couplingStiffnessMatrix[i * 3 + j]*ux[j]);
            couplingForces[i+1]+=( couplingStiffnessMatrix[i * 3 + j]*uy[j]);
            couplingForces[i+2]+=( couplingStiffnessMatrix[i * 3 + j]*uz[j]);
            std::cout<<"ux["<<j<<"]="<<ux[j]<<std::endl;
            std::cout<<"uy["<<j<<"]="<<uy[j]<<std::endl;
            std::cout<<"uz["<<j<<"]="<<uz[j]<<std::endl;
        }
        std::cout<<couplingForces[i]<<" "<<couplingForces[i+1]<<" "<<couplingForces[i+2]<<std::endl;
    }

    setToZero(&force[3 * elementID], 3);

    for(int nID=0 ; nID<nnode ; ++nID){
        globalId = topology[topoPtr + nID];

        for (int i=0 ; i<3 ; i++){  
            force[3*globalId+i] += elNodalForces[3*nID+i]  * volume[globalId];
            // force[3*globalId+i] -= couplingForces[3*nID+i];
        }
    }

    globalId = topology[topoPtr + nnode];

    for (int i=0 ; i<3 ; i++){  
        force[3*globalId+i] -= couplingForces[3 * nnode + i];
        std::cout<< force[3*globalId+i] <<std::endl;
    }
}


void setElementCoordinates
(
    const int nnode,
    const int elementID,
    int topoPtr,
    const int* topology,
    const double* nodalCoor,
    double* elCoor
)
{
    int globalId;

    //setToZero(&elCoor[3 * elementID], 3);

    for(int nID=0 ; nID<nnode ; ++nID){
      globalId = topology[topoPtr + nID];
      
      for (int i=0 ; i<3 ; ++i){    
          elCoor[3 * elementID + i] += nodalCoor[3 * globalId + i] / nnode;
      }
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

{
    double BxiTemp = 0.0, BetaTemp = 0.0, BpsiTemp = 0.0;

    for(int i=0 ; i<9 ; ++i)*(strain+i) = 0.0;   
    
    for(int nID=0 ; nID<nnode ; ++nID){ 
        BxiTemp  = Bx[intPointPtr + nID];
        BetaTemp = By[intPointPtr + nID];

        if (twoD==false)  {
            BpsiTemp = Bz[intPointPtr + nID];
 
        }
    *(strain)   +=  u[3*nID]*(*(Jinv)*BxiTemp + *(Jinv+1)*BetaTemp + *(Jinv+2)*BpsiTemp) ;
    *(strain+4) +=  u[3*nID+1]*(*(Jinv+3)*BxiTemp + *(Jinv+4)*BetaTemp + *(Jinv+5)*BpsiTemp) ;
    *(strain+8) +=  u[3*nID+2]*(*(Jinv+6)*BxiTemp + *(Jinv+7)*BetaTemp + *(Jinv+8)*BpsiTemp) ;
    *(strain+5) +=  0.5*(u[3*nID+1]*(*(Jinv+6)*BxiTemp + *(Jinv+7)*BetaTemp + *(Jinv+8)*BpsiTemp) + u[3*nID+2]*(*(Jinv+3)*BxiTemp + *(Jinv+4)*BetaTemp + *(Jinv+5)*BpsiTemp)) ;
    *(strain+2) +=  0.5*(u[3*nID+2]*(*(Jinv)*BxiTemp + *(Jinv+1)*BetaTemp + *(Jinv+2)*BpsiTemp) + u[3*nID]*(*(Jinv+6)*BxiTemp + *(Jinv+7)*BetaTemp + *(Jinv+8)*BpsiTemp)) ;
    *(strain+1) +=  0.5*(u[3*nID+1]*(*(Jinv)*BxiTemp + *(Jinv+1)*BetaTemp + *(Jinv+2)*BpsiTemp) + u[3*nID]*(*(Jinv+3)*BxiTemp + *(Jinv+4)*BetaTemp + *(Jinv+5)*BpsiTemp)) ;       // 0 1 2
    // 3 4 5
    // 6 7 8
    *(strain+3) = *(strain+1);
    *(strain+6) = *(strain+2);
    *(strain+7) = *(strain+5);


}

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
    for(int i=0 ; i<3*numOwnedPoints ; ++i){
        *(disp+i) = *(coorNP1+i)-*(modelCoord+i);    
    }
}
template<typename ScalarT>
void tensorRotation
(
    const double* angles,
    const ScalarT* tensorIn,
    const bool globToLoc,
    ScalarT* tensorOut
){
    MATRICES::tensorRotation(angles, tensorIn, globToLoc, tensorOut);
}


template void tensorRotation<double>
(
    const double* angles,
    const double* tensorIn,
    const bool globToLoc,
    double* tensorOut
);

template void tensorRotation<Sacado::Fad::DFad<double> >
(
    const double* angles,
    const Sacado::Fad::DFad<double>* tensorIn,
    const bool globToLoc,
    Sacado::Fad::DFad<double>* tensorOut
);

double vectorNorm
(
    const double* vector,
    const int len
){
    return MATRICES::vectorNorm(vector, len);
}

}
