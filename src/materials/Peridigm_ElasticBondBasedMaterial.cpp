/*! \file Peridigm_ElasticBondBasedMaterial.cpp */

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

#include "Peridigm_ElasticBondBasedMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "Peridigm_Logging.hpp"
#include "elastic_bond_based.h"
#include <Teuchos_Assert.hpp>

PeridigmNS::ElasticBondBasedMaterial::ElasticBondBasedMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_bulkModulus(0.0), m_density(0.0), m_horizon(0.0), m_volumeFieldId(-1), m_damageFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1),m_partialStressFieldId(-1)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");
  if(params.isParameter("Young's Modulus") || params.isParameter("Poisson's Ratio") || params.isParameter("Shear Modulus")){
    LOG(LogLevel::WARNING,"The Elastic bond based material model supports only one elastic constant, the bulk modulus. It is calculated from the other elasticity constants if not defined.");
  }
  m_applyThermalStrains = getThermalExpansionCoefficient(params,alpha);
  if (params.isParameter("Apply Thermal Flow")){
    m_lambda = getThermalFlowAndConductivityCoefficients(params, m_C, m_lambdaBed, m_Tbed, m_applyThermalFlow, m_applyThermalPrintBedFlow);
    m_applyHeatTransfer = getHeatTransferCoefficients(params, m_kappa, m_Tenv, m_factor, m_surfaceCorrection, m_limit);
    if (m_applyHeatTransfer){     
      materialProperties["Specific Heat Capacity"] = m_C;
      materialProperties["Thermal Conductivity 11"] = m_lambda[0];
      materialProperties["Thermal Conductivity 22"] = m_lambda[4];
      materialProperties["Thermal Conductivity 33"] = m_lambda[8];
    }
  }
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Volume");
  m_horizonFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_damageFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Damage");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Force_Density");
  m_bondDamageFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Bond_Damage");
  m_partialStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Partial_Stress");
  
  if (m_applyThermalFlow || m_applyThermalStrains){
    m_temperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature");
    m_fieldIds.push_back(m_temperatureFieldId);
    m_detachedNodesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Detached_Nodes");
    m_fieldIds.push_back(m_detachedNodesFieldId);
    if (m_applyThermalFlow){
      m_thermalFlowStateFieldId          = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Flux_Divergence"); // reuse of the field
      m_fieldIds.push_back(m_thermalFlowStateFieldId); 
    }
  }

  
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_partialStressFieldId);
  m_fieldIds.push_back(m_horizonFieldId);
}

PeridigmNS::ElasticBondBasedMaterial::~ElasticBondBasedMaterial()
{
}

void
PeridigmNS::ElasticBondBasedMaterial::initialize(const double dt,
                                                 const int numOwnedPoints,
                                                 const int* ownedIDs,
                                                 const int* neighborhoodList,
                                                 PeridigmNS::DataManager& dataManager)
{
}

void
PeridigmNS::ElasticBondBasedMaterial::computeForce(const double dt,
                                          const int numOwnedPoints,
                                          const int* ownedIDs,
                                          const int* neighborhoodList,
                                          PeridigmNS::DataManager& dataManager,
                                          const double currentTime) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *cellVolume, *bondDamage, *force;

  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);
  double *partialStress;
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&partialStress);
  double *temperature;
  double *detachedNodes;
  double *horizon;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
  if (m_applyThermalFlow || m_applyThermalStrains){
    dataManager.getData(m_temperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&temperature);
    
    dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
  }
  if (m_applyThermalFlow){
    double *thermalFlow; // quadratureWeights not included yet
    
    dataManager.getData(m_thermalFlowStateFieldId, PeridigmField::STEP_NP1)->ExtractView(&thermalFlow);
    dataManager.getData(m_thermalFlowStateFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    MATERIAL_EVALUATION::computeHeatFlowState(
                                  x,
                                  numOwnedPoints,
                                  neighborhoodList,
                                  temperature,
                                  horizon,
                                  m_lambda[0],
                                  cellVolume,
                                  detachedNodes,
                                  bondDamage,
                                  false,
                                  m_applyThermalPrintBedFlow,
                                  m_lambdaBed,
                                  m_Tbed,
                                  thermalFlow);
    if (m_applyHeatTransfer){
      double *specificVolume;
      dataManager.getData(m_specificVolumeFieldId, PeridigmField::STEP_NP1)->ExtractView(&specificVolume);
  
      MATERIAL_EVALUATION::computeHeatTransfer(
                                  x,
                                  numOwnedPoints,
                                  neighborhoodList,
                                  cellVolume,
                                  temperature,
                                  horizon,
                                  detachedNodes,
                                  bondDamage,
                                  false,
                                  m_kappa,
                                  m_Tenv,
                                  m_factor,
                                  m_surfaceCorrection,
                                  m_limit,
                                  m_applyThermalPrintBedFlow,
                                  specificVolume,
                                  thermalFlow);
    }
  }
  MATERIAL_EVALUATION::computeInternalForceElasticBondBased(x,y,cellVolume,bondDamage,neighborhoodList,numOwnedPoints,m_bulkModulus,m_horizon,m_applyThermalStrains,temperature,force,partialStress);
}
