//! \file temperature_diffusion.cxx

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
//@HEADER

#include "temperature_diffusion.h"
#include "matrices.h"
#include <vector> 
#include <cmath> 
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
using namespace std;
namespace DIFFUSION {
  void computeFlux(
    const double* modelCoord,
    const double* temperature,
    const int* neighborhoodList,
    const double* quadratureWeights,
    const int numOwnedPoints,
    const bool useImprovedQuadrature,
    const double* horizon,
    const double coefficient,
    const double* volume,
    double* fluxDivergence
){
    int neighborhoodListIndex(0);
    int numNeighbors, neighborID, iID, iNID, bondListIndex;
    double nodeInitialPosition[3], initialDistance, quadWeight;
    double kernel, nodeTemperature, temperatureDifference, nodeFluxDivergence;//, neighborFluxDivergence;

    const double pi = PeridigmNS::value_of_pi();

    bondListIndex = 0;
    for(iID=0 ; iID<numOwnedPoints ; ++iID){
        nodeTemperature = temperature[iID];
        nodeInitialPosition[0] = modelCoord[iID*3];
        nodeInitialPosition[1] = modelCoord[iID*3+1];
        nodeInitialPosition[2] = modelCoord[iID*3+2];
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        for(iNID=0 ; iNID<numNeighbors ; ++iNID){
        neighborID = neighborhoodList[neighborhoodListIndex++];
        quadWeight = volume[neighborID];
        if (useImprovedQuadrature) {
            quadWeight = quadratureWeights[bondListIndex++];
        }
        initialDistance = MATRICES::distance(nodeInitialPosition[0], nodeInitialPosition[1], nodeInitialPosition[2], modelCoord[neighborID*3], modelCoord[neighborID*3+1], modelCoord[neighborID*3+2]);
        kernel = 6.0/(pi*horizon[iID]*horizon[iID]*horizon[iID]*horizon[iID]*initialDistance);
        temperatureDifference = temperature[neighborID] - nodeTemperature;
        nodeFluxDivergence = coefficient*kernel*temperatureDifference*quadWeight; 
        //TEUCHOS_TEST_FOR_TERMINATION(!std::isfinite(nodeFluxDivergence), "**** NaN detected in DiffusionMaterial::computeFluxDivergence().\n");
        fluxDivergence[iID] += nodeFluxDivergence;

        }
    }
    
  }
  void computeHeatFlowState_correspondence(    
    const double* modelCoord,
    const int numOwnedPoints,
    const int* neighborhoodList,
    const double* shapeTensorInverse,
    const double* temperature,
    const double* horizon,
    const double* kappa,
    const double* volume,
    const double* detachedNodes,
    const double* bondDamage,
    const bool twoD,
    double* heatFlowState
    )
  {
    // based on "A Review of Peridynamics (PD) Theory of Diffusion Based Problems"
    // DOI: 10.1155/2021/7782326
    int neighborhoodListIndex(0), secondNeighborhoodListIndex(0);
    int numNeighbors, neighborID(0), iID, iNID;


    const double *KInv = shapeTensorInverse;
    std::vector<double> Hvector(3), Qvector(3), nablaVector(3), Xvector(3), XpVector(3);
    double* H = &Hvector[0];
    double* X = &Xvector[0];
    double* q = &Qvector[0];
    double* Xp = &XpVector[0];
    double* nablaT = &nablaVector[0];
    double tempState = 0.0; //Eq. (2)
    double temp[3];
    for(iID=0 ; iID<numOwnedPoints ; ++iID, KInv+=9){
      numNeighbors = neighborhoodList[neighborhoodListIndex++];
      if (detachedNodes[iID]!=0){
        neighborhoodListIndex += numNeighbors;
        bondDamage += numNeighbors;
        continue; 
      }
      
      secondNeighborhoodListIndex = neighborhoodListIndex;
      
      for (int i=0 ; i<3 ; ++i) {
        H[i] = 0.0;
        Xp[i] = modelCoord[3*iID+i];
      }

      for(iNID=0 ; iNID<numNeighbors ; ++iNID, bondDamage++){
        neighborID = neighborhoodList[neighborhoodListIndex++];
        for (int i=0 ; i<3 ; ++i) X[i] = modelCoord[3*neighborID+i] - Xp[i];
        tempState = temperature[neighborID] - temperature[iID];
        // sum_j (Tj-Ti)*rij*Vj -> EQ. (8)
        for (int i=0 ; i<3 ; ++i) H[i] += tempState * X[i] * volume[neighborID] * (1 - *bondDamage);
        // std::cout<<*bondDamage<<std::endl;
      }
      // Ki * H -> EQ. (8)
      for (int i=0 ; i<3 ; ++i) {
        nablaT[i] = 0.0;
        for (int j=0 ; j<3 ; ++j) {
          nablaT[i] += KInv[3*i + j] * H[j];
          //std::cout<< nablaT[i]<<std::endl;
        }
      }
      for (int i=0 ; i<3 ; ++i) q[i] = kappa[i] * nablaT[i]; // heat state

      for(iNID=0 ; iNID<numNeighbors ; ++iNID){
        neighborID = neighborhoodList[secondNeighborhoodListIndex++];
        for (int i=0 ; i<3 ; ++i)X[i] = modelCoord[3*neighborID+i] - Xp[i];
        for (int i=0 ; i<3 ; ++i){
          temp[i] = 0.0;
          for (int j=0 ; j<3 ; ++j) {
            temp[i] += KInv[3*i + j] * X[j]; // K * rij -> Eq. (7)
          }
          //if (nablaT[i]!=0)std::cout<<nablaT[i]<<" " << temp<<std::endl;
          
        }

        for (int i=0 ; i<3 ; ++i)  heatFlowState[iID] -= temp[i] * q[i] * volume[neighborID]; // qj * temp -> Eq (7)
        for (int i=0 ; i<3 ; ++i)  heatFlowState[neighborID] += temp[i] * q[i] * volume[iID];
      }
      
    }
  }
  void computeHeatTransfer_correspondence(    
    const int numOwnedPoints,
    const int* neighborhoodList,
    const double* volume,
    const double* temperature,
    const double* horizon,
    const double* detachedNodes,
    const double* bondDamage,
    const bool twoD,
    const double alpha,
    const double Tenv,
    const double factor,
    const double surfaceCorrection,
    double* specificVolume,
    double* heatFlowState
    )
    {
 








    }
  // Peridigm.cpp -> synchro von Detached_nodes
  
  //
  // Peridigm_Correspondence -> Funktionsaufruf + Datenaqkuise
  // (V_neigbor/ V_Horizon) -> als output einpflegen
  //---------------------
  // hier muss Konvektion rein  x x x 0 0 0 
  // V_neighbor => V[i] + summe V[neighbor]*(1-detachedNodes[neighbor])
  // (V_neigbor/ V_Horizon) if < 0.8 -> die 0.8 sind default und können definiert werden
  // V_Horizon = 4/3*pi*delta^3 oder pi*delta^2*h
  //  then   (T - T_umgebung)  * alpha * A * (1 - (V_neigbor/ V_Horizon)^100) 
  // surface coorection factor
  // sqrt^3(V[i])^2 -> alpha*A*surfaceCorrect
  // surfCorrect = 1 als default
  // else q = 0
}