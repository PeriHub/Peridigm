/*! \file Peridigm_CriticalTimeStep.cpp */

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

#include "Peridigm_CriticalTimeStep.hpp"
#include "Peridigm_HorizonManager.hpp"
#include "Peridigm_Field.hpp"
#include "Peridigm_Constants.hpp"
#include "Peridigm_DegreesOfFreedomManager.hpp"
#include <cmath>
// #include <math.h>

double PeridigmNS::ComputeCriticalTimeStep(const Epetra_Comm& comm, PeridigmNS::Block& block){

  PeridigmNS::DegreesOfFreedomManager &dofManager = PeridigmNS::DegreesOfFreedomManager::self();
  bool solveForDisplacement = dofManager.displacementTreatedAsUnknown();
  bool solveForTemperature = dofManager.temperatureTreatedAsUnknown();

  Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = block.getNeighborhoodData();
  const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
  const int* ownedIDs = neighborhoodData->OwnedIDs();
  const int* neighborhoodList = neighborhoodData->NeighborhoodList();
  Teuchos::RCP<const PeridigmNS::Material> materialModel = block.getMaterialModel();

  double density = materialModel()->Density();
  double bulkModulus = materialModel()->BulkModulus();
  double heatCapacity(0.0), lambda11(0.0), lambda22(0.0), lambda33(0.0), maxLambda(0.0);
  if(!solveForDisplacement && solveForTemperature){
    heatCapacity = materialModel()->lookupMaterialProperty("Specific Heat Capacity");
    lambda11 = materialModel()->lookupMaterialProperty("Thermal Conductivity 11");
    lambda22 = materialModel()->lookupMaterialProperty("Thermal Conductivity 22");
    lambda33 = materialModel()->lookupMaterialProperty("Thermal Conductivity 33");
    maxLambda = std::max(lambda11,std::max(lambda22,lambda33));
  }

  double horizon(0.0);
  std::string blockName = block.getName();
  PeridigmNS::HorizonManager& horizonManager = PeridigmNS::HorizonManager::self();
  bool blockHasConstantHorizon = horizonManager.blockHasConstantHorizon(blockName);
  if(blockHasConstantHorizon)
    horizon = horizonManager.getBlockConstantHorizonValue(blockName);

  double *cellVolume, *x;
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  block.getData(fieldManager.getFieldId("Volume"), PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
  block.getData(fieldManager.getFieldId("Model_Coordinates"), PeridigmField::STEP_NONE)->ExtractView(&x);

  const double pi = value_of_pi();
  double springConstant(0.0);
  if(blockHasConstantHorizon)
    springConstant = 18.0*bulkModulus/(pi*horizon*horizon*horizon*horizon);

  double minCriticalTimeStep = 1.0e50;

  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    double timestepDenominator = 0.0;

    int nodeID = ownedIDs[iID];
    double X[3] = { x[nodeID*3], x[nodeID*3+1], x[nodeID*3+2] };
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];

    if(!blockHasConstantHorizon){
      double delta = horizonManager.evaluateHorizon(blockName, X[0], X[1], X[2]);
      springConstant = 18.0*bulkModulus/(pi*delta*delta*delta*delta);
    }

    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      double neighborVolume = cellVolume[neighborID];
      double initialDistance = std::sqrt( (X[0] - x[neighborID*3  ])*(X[0] - x[neighborID*3  ]) +
                                          (X[1] - x[neighborID*3+1])*(X[1] - x[neighborID*3+1]) +
                                          (X[2] - x[neighborID*3+2])*(X[2] - x[neighborID*3+2]) );

      // Issue a warning if the bond length is very very small (as in zero)
      static bool warningGiven = false;
      if(!warningGiven && initialDistance < 1.0e-50){
        std::cout << "\nWarning:  Possible zero length bond detected (length = " << initialDistance << ")." << std::endl;
        std::cout << "            Bonds of length zero are not valid, the input mesh may contain coincident nodes.\n" << std::endl;
        warningGiven = true;
      }

      if(!solveForDisplacement && solveForTemperature){
        timestepDenominator += ( maxLambda / initialDistance) * neighborVolume;
      }
      else{
        timestepDenominator += neighborVolume*springConstant/initialDistance;
      }
    }

    double criticalTimeStep = 1.0e50;
    if(numNeighbors > 0){
      if(!solveForDisplacement && solveForTemperature){
        // Selda Oterkus, Erdogan Madenci, and Abigail G. Agwai.  Fully coupled peridynamic thermomechanics
        criticalTimeStep = (density * heatCapacity) / timestepDenominator;
      }
      else{
        criticalTimeStep = sqrt(2.0*density/timestepDenominator);
      }
    }
    if(criticalTimeStep < minCriticalTimeStep)
      minCriticalTimeStep = criticalTimeStep;
  }

  // Find the minimum time step for this block across all processors
  double globalMinCriticalTimeStep;
  comm.MinAll(&minCriticalTimeStep, &globalMinCriticalTimeStep, 1);

  return globalMinCriticalTimeStep;
}
