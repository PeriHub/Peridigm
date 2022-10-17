/*! \file Peridigm_SimpleAdditiveModel.cpp */

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

#include "Peridigm_SimpleAdditiveModel.hpp"
#include "Peridigm_Field.hpp"
#include "additive_utilities.h"
using namespace std;

PeridigmNS::SimpleAdditiveModel::SimpleAdditiveModel(const Teuchos::ParameterList& params)
  : AdditiveModel(params), m_applyThermalStrains(false), m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_bondDamageFieldId(-1), m_deltaTemperatureFieldId(-1)
{
  
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  
  m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::NODE, PeridigmNS::PeridigmField::VECTOR, PeridigmNS::PeridigmField::CONSTANT,"Model_Coordinates");
  m_coordinatesFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::NODE, PeridigmField::VECTOR, PeridigmNS::PeridigmField::TWO_STEP, "Coordinates");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bond_Damage");
  m_volumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_detachedNodesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Detached_Nodes");
  if(m_applyThermalStrains)
    m_deltaTemperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");


  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_detachedNodesFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_volumeFieldId);
  if(m_applyThermalStrains)
    m_fieldIds.push_back(m_deltaTemperatureFieldId);

}

PeridigmNS::SimpleAdditiveModel::~SimpleAdditiveModel()
{
}

void
PeridigmNS::SimpleAdditiveModel::initialize(const double dt,
                                                   const int numOwnedPoints,
                                                   const int* ownedIDs,
                                                   const int* neighborhoodList,
                                                   PeridigmNS::DataManager& dataManager) const
{
  double *bondDamage, *detachedNodes;
  
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamage);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
  // Initialize damage to zero
  ADDITIVE_UTILITIES::deleteAllBonds(numOwnedPoints, ownedIDs, neighborhoodList, bondDamage, detachedNodes);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)=dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N);
}

void
PeridigmNS::SimpleAdditiveModel::computeAdditive(const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,
                                                      PeridigmNS::DataManager& dataManager) const
{
  
  }


    
 



