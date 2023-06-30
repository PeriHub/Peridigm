/*! \file Peridigm_ElasticLinearCorrespondenceMaterial.cpp */

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

#include "Peridigm_ElasticLinearCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ElasticLinearCorrespondenceMaterial::ElasticLinearCorrespondenceMaterial(const Teuchos::ParameterList &params)
    : CorrespondenceMaterial(params),
      m_modelAnglesId(-1),
      m_deformationGradientFieldId(-1),
      m_cauchyStressFieldId(-1),
      m_strain(-1),
      m_temperatureFieldId(-1)
{
  bool m_planeStrain = false, m_planeStress = false;
  m_type = 0;
  m_density = params.get<double>("Density");
  if (params.isParameter("Plane Strain"))
  {
    m_planeStrain = params.get<bool>("Plane Strain");
  }
  if (params.isParameter("Plane Stress"))
  {
    m_planeStress = params.get<bool>("Plane Stress");
  }
  if (m_planeStrain == true)
    m_type = 1;
  if (m_planeStress == true)
    m_type = 2;

  m_hencky = false;
  if (params.isParameter("Hencky Strain"))
  {
    m_hencky = params.get<bool>("Hencky Strain");
  }
  getStiffnessmatrix(params, C, m_planeStrain, m_planeStress);
  m_applyThermalStrains = getThermalExpansionCoefficient(params,alpha,m_Tref);
  
  PeridigmNS::FieldManager &fieldManager = PeridigmNS::FieldManager::self();

  m_cauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_deformationGradientFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Deformation_Gradient");
  m_strain = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Strain");
  m_modelAnglesId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Local_Angles");
  if (m_applyThermalStrains)
  {
    m_temperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature");
    m_fieldIds.push_back(m_temperatureFieldId);
  }
  m_fieldIds.push_back(m_modelAnglesId);
  m_fieldIds.push_back(m_deformationGradientFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);
  m_fieldIds.push_back(m_strain);
  
  m_vonMisesStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Von_Mises_Stress");
  m_fieldIds.push_back(m_vonMisesStressFieldId);
}

PeridigmNS::ElasticLinearCorrespondenceMaterial::~ElasticLinearCorrespondenceMaterial()
{
}

void PeridigmNS::ElasticLinearCorrespondenceMaterial::initialize(const double dt,
                                                                 const int numOwnedPoints,
                                                                 const int *ownedIDs,
                                                                 const int *neighborhoodList,
                                                                 PeridigmNS::DataManager &dataManager)
{
  PeridigmNS::CorrespondenceMaterial::initialize(dt,
                                                 numOwnedPoints,
                                                 ownedIDs,
                                                 neighborhoodList,
                                                 dataManager);
                                                 
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);


}

void PeridigmNS::ElasticLinearCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                                          const int numOwnedPoints,
                                                                          PeridigmNS::DataManager &dataManager,
                                                                          const double time) const
{

  double *CauchyStress, *CauchyStressNP1, *defGrad, *angles, *temperature, *strain;
  // have to be checked if the additional effort is useful or not
  // deactivated for tests with implicit solver

  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&CauchyStress);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&CauchyStressNP1);

  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NP1)->ExtractView(&defGrad);

  double *vonMisesStress;
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&vonMisesStress);

  if (m_applyThermalStrains) dataManager.getData(m_temperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&temperature);
  dataManager.getData(m_strain, PeridigmField::STEP_NONE)->ExtractView(&strain);
  dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&angles);

  CORRESPONDENCE::getStrain(numOwnedPoints, defGrad, alpha, temperature, m_Tref, m_hencky, m_applyThermalStrains, strain);


  CORRESPONDENCE::updateElasticCauchyStressAnisotropic(strain,
                                                       CauchyStress,
                                                       CauchyStressNP1,
                                                       numOwnedPoints,
                                                       C,
                                                       angles,
                                                       m_type,
                                                       dt
                                                       );

  CORRESPONDENCE::getVonMisesStress(numOwnedPoints, CauchyStressNP1, vonMisesStress);
}
