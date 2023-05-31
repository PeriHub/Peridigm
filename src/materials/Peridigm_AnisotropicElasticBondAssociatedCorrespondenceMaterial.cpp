/*! \file Peridigm_AnisotropicElasticBondAssociatedCorrespondenceMaterial.cpp */

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

#include "Peridigm_AnisotropicElasticBondAssociatedCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::AnisotropicElasticBondAssociatedCorrespondenceMaterial::AnisotropicElasticBondAssociatedCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : BondAssociatedCorrespondenceMaterial(params),
    m_vonMisesStressFieldId(-1), 
    m_bondLevelVonMisesStressFieldId(-1),m_planeStrain(false),m_planeStress(false),m_type(0),m_applyThermalStrains(false)
{
  bool m_planeStrain = false, m_planeStress = false;


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

  getStiffnessmatrix(params, C, m_planeStrain, m_planeStress);
  m_applyThermalStrains = getThermalExpansionCoefficient(params,alpha);
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_vonMisesStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Von_Mises_Stress");
  m_deformationGradientFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Deformation_Gradient");

  m_modelAnglesId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Local_Angles");
  if (m_applyThermalStrains)
  {
    m_temperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature");
    m_fieldIds.push_back(m_temperatureFieldId);
  }
  m_strain = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Strain");
  m_fieldIds.push_back(m_modelAnglesId);
  m_fieldIds.push_back(m_deformationGradientFieldId);
  m_fieldIds.push_back(m_strain);
  m_fieldIds.push_back(m_vonMisesStressFieldId);

  m_bondLevelVonMisesStressFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Bond_Von_Mises_Stress");

  m_fieldIds.push_back(m_bondLevelVonMisesStressFieldId);
}

PeridigmNS::AnisotropicElasticBondAssociatedCorrespondenceMaterial::~AnisotropicElasticBondAssociatedCorrespondenceMaterial()
{
}

void
PeridigmNS::AnisotropicElasticBondAssociatedCorrespondenceMaterial::initialize(const double dt,
                                                                    const int numOwnedPoints,
                                                                    const int* ownedIDs,
                                                                    const int* neighborhoodList,
                                                                    PeridigmNS::DataManager& dataManager)
{

  PeridigmNS::BondAssociatedCorrespondenceMaterial::initialize(dt,
                                                               numOwnedPoints,
                                                               ownedIDs,
                                                               neighborhoodList,
                                                               dataManager);

  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelVonMisesStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
}

void
PeridigmNS::AnisotropicElasticBondAssociatedCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                                             const int numOwnedPoints,
                                                                             const int* neighborhoodList,
                                                                             PeridigmNS::DataManager& dataManager) const
  {

  // Compute the node-level stress values
  // This is only done for output (visualization) purposes
  double *unrotatedCauchyStressN;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressN);

  double *unrotatedCauchyStressNP1;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);

  double *unrotatedRateOfDeformation;
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);

  double *vonMisesStress;
  double  *defGrad, *angles, *temperature, *strain;
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&vonMisesStress);
  dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&angles);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NP1)->ExtractView(&defGrad);
  dataManager.getData(m_strain, PeridigmField::STEP_NONE)->ExtractView(&strain);
  CORRESPONDENCE::getStrain(numOwnedPoints, defGrad, alpha, temperature, false, m_applyThermalStrains, strain);


  CORRESPONDENCE::updateElasticCauchyStressAnisotropic(strain,
                                                       unrotatedCauchyStressN,
                                                       unrotatedCauchyStressNP1,
                                                       numOwnedPoints,
                                                       C,
                                                       angles,
                                                       m_type,
                                                       dt
                                                       );
  CORRESPONDENCE::updateElasticCauchyStress(unrotatedRateOfDeformation, 
                                            unrotatedCauchyStressN, 
                                            unrotatedCauchyStressNP1,
                                            vonMisesStress,
                                            numOwnedPoints,
                                            m_bulkModulus,
                                            m_shearModulus,
                                            dt);

  // Compute the bond-level stress values
  double *bondLevelUnrotatedCauchyStressXXN, *bondLevelUnrotatedCauchyStressXYN, *bondLevelUnrotatedCauchyStressXZN;
  double *bondLevelUnrotatedCauchyStressYXN, *bondLevelUnrotatedCauchyStressYYN, *bondLevelUnrotatedCauchyStressYZN;
  double *bondLevelUnrotatedCauchyStressZXN, *bondLevelUnrotatedCauchyStressZYN, *bondLevelUnrotatedCauchyStressZZN;
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressXXN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressXYN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressXZN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressYXN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressYYN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressYZN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressZXN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressZYN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressZZN);

  double *bondLevelUnrotatedCauchyStressXXNP1, *bondLevelUnrotatedCauchyStressXYNP1, *bondLevelUnrotatedCauchyStressXZNP1;
  double *bondLevelUnrotatedCauchyStressYXNP1, *bondLevelUnrotatedCauchyStressYYNP1, *bondLevelUnrotatedCauchyStressYZNP1;
  double *bondLevelUnrotatedCauchyStressZXNP1, *bondLevelUnrotatedCauchyStressZYNP1, *bondLevelUnrotatedCauchyStressZZNP1;
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXZNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYZNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZZNP1);

  double *bondLevelUnrotatedRateOfDeformationXX, *bondLevelUnrotatedRateOfDeformationXY, *bondLevelUnrotatedRateOfDeformationXZ;
  double *bondLevelUnrotatedRateOfDeformationYX, *bondLevelUnrotatedRateOfDeformationYY, *bondLevelUnrotatedRateOfDeformationYZ;
  double *bondLevelUnrotatedRateOfDeformationZX, *bondLevelUnrotatedRateOfDeformationZY, *bondLevelUnrotatedRateOfDeformationZZ;
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXZ);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYZ);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZZ);
  double *flyingPoint;
  double *bondLevelVonMisesStress;
  dataManager.getData(m_bondLevelVonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVonMisesStress);

  CORRESPONDENCE::updateBondLevelAnisotropicElasticCauchyStress(bondLevelUnrotatedRateOfDeformationXX,
                                                     bondLevelUnrotatedRateOfDeformationXY,
                                                     bondLevelUnrotatedRateOfDeformationXZ,
                                                     bondLevelUnrotatedRateOfDeformationYX,
                                                     bondLevelUnrotatedRateOfDeformationYY,
                                                     bondLevelUnrotatedRateOfDeformationYZ,
                                                     bondLevelUnrotatedRateOfDeformationZX,
                                                     bondLevelUnrotatedRateOfDeformationZY,
                                                     bondLevelUnrotatedRateOfDeformationZZ,
                                                     bondLevelUnrotatedCauchyStressXXN, 
                                                     bondLevelUnrotatedCauchyStressXYN, 
                                                     bondLevelUnrotatedCauchyStressXZN, 
                                                     bondLevelUnrotatedCauchyStressYXN, 
                                                     bondLevelUnrotatedCauchyStressYYN, 
                                                     bondLevelUnrotatedCauchyStressYZN, 
                                                     bondLevelUnrotatedCauchyStressZXN, 
                                                     bondLevelUnrotatedCauchyStressZYN, 
                                                     bondLevelUnrotatedCauchyStressZZN, 
                                                     bondLevelUnrotatedCauchyStressXXNP1, 
                                                     bondLevelUnrotatedCauchyStressXYNP1, 
                                                     bondLevelUnrotatedCauchyStressXZNP1, 
                                                     bondLevelUnrotatedCauchyStressYXNP1, 
                                                     bondLevelUnrotatedCauchyStressYYNP1, 
                                                     bondLevelUnrotatedCauchyStressYZNP1, 
                                                     bondLevelUnrotatedCauchyStressZXNP1, 
                                                     bondLevelUnrotatedCauchyStressZYNP1, 
                                                     bondLevelUnrotatedCauchyStressZZNP1, 
                                                     flyingPoint,
                                                     neighborhoodList,
                                                     numOwnedPoints,
                                                     C,
                                                     angles,
                                                     m_type,
                                                     dt);
}
