/*! \file Peridigm_ModelEvaluator.hpp */

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
// Adapted by C. Willberg 
// Christian Willberg    christian.willberg@dlr.de
// ************************************************************************
//@HEADER

#include "Peridigm_ModelEvaluator.hpp"
#include "Peridigm_Timer.hpp"

PeridigmNS::ModelEvaluator::ModelEvaluator(){
}

PeridigmNS::ModelEvaluator::~ModelEvaluator(){
}
void 
PeridigmNS::ModelEvaluator::evalDamageModel(Teuchos::RCP<Workset> workset) const
{
  const double dt = workset->timeStep;
  const double currentTime = workset->currentTime;
  std::vector<PeridigmNS::Block>::iterator blockIt;

  // ---- Evaluate Damage ---

  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<const PeridigmNS::DamageModel> damageModel = blockIt->getDamageModel();
    if(!damageModel.is_null() && blockIt->getDamageEnabled()){
      Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
      const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
      const int* ownedIDs = neighborhoodData->OwnedIDs();
      const int* neighborhoodList = neighborhoodData->NeighborhoodList();
      int blockInterfaceId = blockIt->getBlockInterfaceID();
      Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
      
      PeridigmNS::Timer::self().startTimer("Evaluate Damage Model:Compute Damage");
      damageModel->computeDamage(dt, 
                                 numOwnedPoints,
                                 ownedIDs,
                                 neighborhoodList,
                                 *dataManager,
                                 blockInterfaceId,
                                 currentTime);
      PeridigmNS::Timer::self().stopTimer("Evaluate Damage Model:Compute Damage");
    }
  }
}
void
PeridigmNS::ModelEvaluator::evalModel(Teuchos::RCP<Workset> workset, bool damageExist) const
{
  const double dt = workset->timeStep;
  const double currentTime = workset->currentTime;
  std::vector<PeridigmNS::Block>::iterator blockIt;

  // ---- Evaluate Precompute ----
  
  PeridigmNS::Timer::self().startTimer("Internal Force:Evaluate Precompute");
  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* ownedIDs = neighborhoodData->OwnedIDs();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::Material> materialModel = blockIt->getMaterialModel();

    materialModel->precompute(dt,
                              numOwnedPoints,
                              ownedIDs,
                              neighborhoodList,
                              *dataManager);
  }
  PeridigmNS::Timer::self().stopTimer("Internal Force:Evaluate Precompute");

  // ---- Synchronize data computed in precompute ----
 PeridigmNS::DataManagerSynchronizer::self().synchronizeDataAfterPrecompute(workset->blocks);

 // PeridigmNS::DataManagerSynchronizer::self().synchronizeDataAfterPrecompute(workset->blocks);

  // ---- Evaluate Internal Force ----

  PeridigmNS::Timer::self().startTimer("Internal Force:Evaluate Internal Force");
  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* ownedIDs = neighborhoodData->OwnedIDs();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::Material> materialModel = blockIt->getMaterialModel();

    bool runEval = true;
    
    if (materialModel->Name().find("Correspondence")!=std::string::npos) runEval = damageExist;
    
    
    if(runEval){
      PeridigmNS::Timer::self().startTimer("Internal Force:Evaluate Internal Force:Compute Force");
      materialModel->computeForce(dt,
                                  numOwnedPoints,
                                  ownedIDs,
                                  neighborhoodList,
                                  *dataManager,
                                  currentTime);
      PeridigmNS::Timer::self().stopTimer("Internal Force:Evaluate Internal Force:Compute Force");
    }
    
    PeridigmNS::Timer::self().startTimer("Internal Force:Evaluate Internal Force:Flux Divergence");
    materialModel->computeFluxDivergence(dt,
                                         numOwnedPoints,
                                         ownedIDs,
                                         neighborhoodList,
                                         *dataManager);
    PeridigmNS::Timer::self().stopTimer("Internal Force:Evaluate Internal Force:Flux Divergence");
  }
  PeridigmNS::Timer::self().stopTimer("Internal Force:Evaluate Internal Force");

  // ---- Evaluate Contact ----

  PeridigmNS::Timer::self().startTimer("Internal Force:Evaluate Contact");
  if(!workset->contactManager.is_null())
    workset->contactManager->evaluateContactForce(dt);
  PeridigmNS::Timer::self().stopTimer("Internal Force:Evaluate Contact");
}

void
PeridigmNS::ModelEvaluator::evalJacobian(Teuchos::RCP<Workset> workset) const
{
  const double dt = workset->timeStep;
  std::vector<PeridigmNS::Block>::iterator blockIt;
  PeridigmNS::Material::JacobianType jacobianType = *(workset->jacobianType);
  PeridigmNS::SerialMatrix& jacobian = *(workset->jacobian);

  // ---- Compute the Tangent Stiffness Matrix ----

  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* ownedIDs = neighborhoodData->OwnedIDs();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::Material> materialModel = blockIt->getMaterialModel();

    materialModel->computeJacobian(dt,
                                   numOwnedPoints,
                                   ownedIDs,
                                   neighborhoodList,
                                   *dataManager,
                                   jacobian,
                                   jacobianType);
  }
}

void 
PeridigmNS::ModelEvaluator::computeVelocityGradient(Teuchos::RCP<Workset> workset) const
{
  const double dt = workset->timeStep;
  std::vector<PeridigmNS::Block>::iterator blockIt;

  // ---- Compute the Velocity Gradient for Hypoelastic Model ----

  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* ownedIDs = neighborhoodData->OwnedIDs();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::Material> materialModel = blockIt->getMaterialModel();

    materialModel->computeNodeLevelVelocityGradient(dt, 
                                                    numOwnedPoints,
                                                    ownedIDs,
                                                    neighborhoodList,
                                                    *dataManager);
  }
}

void 
PeridigmNS::ModelEvaluator::computeBondVelocityGradient(Teuchos::RCP<Workset> workset) const
{
  //const double dt = workset->timeStep;
  std::vector<PeridigmNS::Block>::iterator blockIt;

  // ---- Compute the Velocity Gradient for Hypoelastic Model ----

  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    //const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    //const int* ownedIDs = neighborhoodData->OwnedIDs();
    //const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::Material> materialModel = blockIt->getMaterialModel();

 //   materialModel->computeBondVelocityGradient(dt, 
 //                                              numOwnedPoints,
 //                                              ownedIDs,
 //                                              neighborhoodList,
 //                                              *dataManager);
  }
}
void 
PeridigmNS::ModelEvaluator::updateDilatation(Teuchos::RCP<Workset> workset) const
{

  std::vector<PeridigmNS::Block>::iterator blockIt;

  // All blocks, because nodes outside the block definition are needed for the interface bond calculations  
  // what if no damage is included in a block
  // test 
  const double dt = workset->timeStep;
  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* ownedIDs = neighborhoodData->OwnedIDs();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::Material> materialModel = blockIt->getMaterialModel();
    
    if (materialModel->Name() == "Elastic"){
        
        materialModel->evalDilatation(dt, 
                                numOwnedPoints,
                                ownedIDs,
                                neighborhoodList,
                                *dataManager);
    
    }
  }
}

void 
PeridigmNS::ModelEvaluator::updateCauchyStress(Teuchos::RCP<Workset> workset) const
{
  const double dt = workset->timeStep;
  std::vector<PeridigmNS::Block>::iterator blockIt;

  // All blocks, because nodes outside the block definition are needed for the interface bond calculations  
  // what if no damage is included in a block
  // test 

  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* ownedIDs = neighborhoodData->OwnedIDs();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::Material> materialModel = blockIt->getMaterialModel();
    Teuchos::RCP<const PeridigmNS::DamageModel> damageModel = blockIt->getDamageModel();

            if (materialModel->Name().find("Correspondence")!=std::string::npos){
                PeridigmNS::Timer::self().startTimer("Update Cauchy Stress:Compute Force");
                materialModel->computeForce(dt, 
                                numOwnedPoints,
                                ownedIDs,
                                neighborhoodList,
                                *dataManager);
                PeridigmNS::Timer::self().stopTimer("Update Cauchy Stress:Compute Force");
                                
            }

  }
}