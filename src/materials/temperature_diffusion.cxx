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


using namespace std;
namespace DIFFUSION {
  void computeFlux(
    const double* modelCoord,
    const double* temperatureN,
    double* temperatureNP1,
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

    const double pi = 3.1415; //::value_of_pi();

    bondListIndex = 0;
    for(iID=0 ; iID<numOwnedPoints ; ++iID){
        nodeTemperature = temperatureN[iID];
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
        temperatureDifference = temperatureN[neighborID] - nodeTemperature;
        nodeFluxDivergence = coefficient*kernel*temperatureDifference*quadWeight; 
        //TEUCHOS_TEST_FOR_TERMINATION(!std::isfinite(nodeFluxDivergence), "**** NaN detected in DiffusionMaterial::computeFluxDivergence().\n");
        fluxDivergence[iID] += nodeFluxDivergence;

        }
    }
    neighborhoodListIndex = 0;
    for(iID=0 ; iID<numOwnedPoints ; ++iID){
      nodeTemperature = temperatureN[iID];
      numNeighbors = neighborhoodList[neighborhoodListIndex++];
      for(iNID=0 ; iNID<numNeighbors ; ++iNID){
        neighborID = neighborhoodList[neighborhoodListIndex++];
                quadWeight = volume[neighborID];
        if (useImprovedQuadrature) {
            quadWeight = quadratureWeights[bondListIndex++];
        }
        initialDistance = MATRICES::distance(nodeInitialPosition[0], nodeInitialPosition[1], nodeInitialPosition[2], modelCoord[neighborID*3], modelCoord[neighborID*3+1], modelCoord[neighborID*3+2]);
        kernel = 6.0/(pi*horizon[iID]*horizon[iID]*horizon[iID]*horizon[iID]*initialDistance);
        temperatureNP1[neighborID]  += fluxDivergence[iID] / (coefficient*kernel*quadWeight) + nodeTemperature ;
        }
    }
  }
  void computeHeatFlux_correspondence(    
    const double* modelCoord,
    const int numOwnedPoints,
    const int* neighborhoodList,
    const double* shapeTensorInverse,
    const double* temperature,
    const double* horizon,
    const double* kappa,
    const double* volume,
    const double* bondDamage,
    double* fluxDivergence
    )
  {

    int neighborhoodListIndex(0);
    int numNeighbors, neighborID, iID, iNID, bondListIndex;
    double undeformedBondX[3], initialDistance, quadWeight;
    double kernel[3], nodeTemperature, temperatureDifference, nodeFluxDivergence;//, neighborFluxDivergence;

    const double pi = 3.1415; //::value_of_pi();
    double H[3];
    bondListIndex = 0;
    for(iID=0 ; iID<numOwnedPoints ; ++iID){
        nodeTemperature = temperature[iID];
        undeformedBondX[0] = modelCoord[iID*3];
        undeformedBondX[1] = modelCoord[iID*3+1];
        undeformedBondX[2] = modelCoord[iID*3+2];
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        H[0] = 0.0;H[1] = 0.0;H[2] = 0.0;
        for(iNID=0 ; iNID<numNeighbors ; ++iNID, bondDamage++){
          neighborID = neighborhoodList[neighborhoodListIndex++];
          

          initialDistance = MATRICES::distance(undeformedBondX[0], undeformedBondX[1], undeformedBondX[2], modelCoord[neighborID*3], modelCoord[neighborID*3+1], modelCoord[neighborID*3+2]);
          for(int i=0 ; i<3 ; ++i){ 
            kernel[i] = 6.0*kappa[i]/(pi*horizon[iID]*horizon[iID]*horizon[iID]*horizon[iID]*initialDistance);
            }
          
          
          H[0] += (1-*bondDamage)*(temperature[neighborID] - nodeTemperature)*(modelCoord[neighborID*3] -undeformedBondX[0])*volume[neighborID];
          H[1] += (1-*bondDamage)*(temperature[neighborID] - nodeTemperature)*(modelCoord[neighborID*3+1]-undeformedBondX[1])*volume[neighborID+2];
          H[2] += (1-*bondDamage)*(temperature[neighborID] - nodeTemperature)*(modelCoord[neighborID*3+2]-undeformedBondX[2])*volume[neighborID];
          
          //nodeFluxDivergence = coefficient*kernel*temperatureDifference*quadWeight; 
          //TEUCHOS_TEST_FOR_TERMINATION(!std::isfinite(nodeFluxDivergence), "**** NaN detected in DiffusionMaterial::computeFluxDivergence().\n");
          fluxDivergence[iID] += nodeFluxDivergence;

        }
        //MATRICES::MatrixMultiply3x3toVector(H,K,nablaT);
        //MATRICES::MatrixMultiply3x3toVector(nablaT,kappa,flux);
        //heatconductionstate = flux*Kinv*dist; // rausgehen
    }
  }
}