/*! \file Peridigm_UserCorrespondenceMaterial.cpp */

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
// funded by EMMA project reference
// Licence agreement
//
// Christian Willberg    christian.willberg@dlr.de
//@HEADER

#include "Peridigm_UserCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "user_material_interface_correspondence.h"
#include "correspondence.h"
#include <Teuchos_Assert.hpp>
#include <Teuchos_Exceptions.hpp>

//using namespace std;

PeridigmNS::UserCorrespondenceMaterial::UserCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_modelCoordinatesFieldId(-1),
    m_deformationGradientFieldId(-1),
    m_cauchyStressFieldId(-1),
    m_strainFieldId(-1),
    m_modelAnglesId(-1),
    m_flyingPointFlagFieldId(-1),
    m_rotationTensorFieldId(-1)    
{
  m_planeStrain = false;
  m_planeStress = false;
  m_type = 0;
  m_density = params.get<double>("Density");
  if (params.isParameter("Plane Strain")){
    m_planeStrain = params.get<bool>("Plane Strain");
    }
  if (params.isParameter("Plane Stress")){
    m_planeStress = params.get<bool>("Plane Stress");
    }
  if (m_planeStrain==true)m_type=1;
  if (m_planeStress==true)m_type=2;
 
  std::string var =params.name();
  std::string delimiter = "->";
  matName = var.substr(var.find_last_of(delimiter)+1,var.length());

  nprops = params.get<int>("Number of Properties");
  TEUCHOS_TEST_FOR_TERMINATION(nprops<1, 
     "****         The number properties must be greater than zero.\n");
  string prop = "Prop_";
  
  // https://docs.microsoft.com/de-de/cpp/cpp/delete-operator-cpp?view=msvc-170
  delete userProperties;
  userProperties = new double[nprops];

  for(int iID=1 ; iID<nprops+1 ; ++iID)
  {
    userProperties[iID-1] = params.get<double>(prop + std::to_string(iID));
  }
 
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  
  m_modelCoordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_cauchyStressFieldId                 = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_strainFieldId                       = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Strains");
  m_deformationGradientFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Deformation_Gradient");
  m_modelAnglesId                       = fieldManager.getFieldId(PeridigmField::NODE   , PeridigmField::VECTOR, PeridigmField::CONSTANT     , "Local_Angles");
  m_flyingPointFlagFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Flying_Point_Flag");
  m_rotationTensorFieldId               = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Rotation_Tensor");
  nstatev = 0;
  if (params.isParameter("Number of State Vars")){
    nstatev = params.get<int>("Number of State Vars");
    // FULL_TENSOR is used. Therefore, nine entries will be inialized at ones. 
    TEUCHOS_TEST_FOR_TERMINATION(nstatev<1, 
     "****         The number of state variables must be greater than zero.\n");

    int nstat = int(ceil(nstatev / 9.0));
    if (nstat > 0) {

      delete m_state;
      m_state = new int[nstat];
      prop = "State_Parameter_Field_";
      for(int iID=0 ; iID<nstat ; ++iID)
      {
        m_state[iID]                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, prop + std::to_string(iID+1));
        m_fieldIds.push_back(m_state[iID]);
      }
    }
  }

  m_modelAnglesId                       = fieldManager.getFieldId(PeridigmField::NODE   , PeridigmField::VECTOR, PeridigmField::CONSTANT     , "Local_Angles");
  
  m_fieldIds.push_back(m_strainFieldId);
  m_fieldIds.push_back(m_modelAnglesId);
  m_fieldIds.push_back(m_deformationGradientFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);
  m_fieldIds.push_back(m_flyingPointFlagFieldId);
  m_fieldIds.push_back(m_rotationTensorFieldId);

}

PeridigmNS::UserCorrespondenceMaterial::~UserCorrespondenceMaterial()
{
}

void
PeridigmNS::UserCorrespondenceMaterial::initialize(const double dt,
                                                             const int numOwnedPoints,
                                                             const int* ownedIDs,
                                                             const int* neighborhoodList,
                                                             PeridigmNS::DataManager& dataManager)
{
      PeridigmNS::CorrespondenceMaterial::initialize(dt,
                                                      numOwnedPoints,
                                                      ownedIDs,
                                                      neighborhoodList,
                                                      dataManager);
      dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_N)->PutScalar(-1.0);
      dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_NP1)->PutScalar(-1.0);
     
      double *angles;
      dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&angles);
      /*
       check if and where local coordinates exists; 
      */
      

}

void
PeridigmNS::UserCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                            const int numOwnedPoints,
                                                            PeridigmNS::DataManager& dataManager,
                                                            const double time) const
{
  ///////////////
  //PLACEHOLDER//
  ///////////////
  double *temperature = NULL, *dtemperature = NULL;
  // double time = 0.0;
  //////////////////////////////

  double *angles, *modelCoordinates;
  double *defGradN, *defGradNP1;
  double *GLStrainN, *GLStrainNP1;
  double *CauchyStressN, *CauchyStressNP1;
  double *RotationN, *RotationNP1;
  // have to be checked if the additional effort is useful or not
  // deactivated for tests with implicit solver
 
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&angles);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_N)->ExtractView(&defGradN);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NP1)->ExtractView(&defGradNP1);
  dataManager.getData(m_strainFieldId,              PeridigmField::STEP_N)->ExtractView(&GLStrainN);
  dataManager.getData(m_strainFieldId,              PeridigmField::STEP_NP1)->ExtractView(&GLStrainNP1);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&CauchyStressN);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&CauchyStressNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&RotationN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&RotationNP1);
 
  
  double *flyingPointFlag;
  dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_N)->ExtractView(&flyingPointFlag);

  double *statev = new double[nstatev*numOwnedPoints];
  // fill with data
  int nstat = int(ceil(nstatev / 9.0));
  if (nstat > 0) {
      double *stat;
      int tensorLen = 9;
      for(int iID=0 ; iID<nstat ; ++iID)
      { 
        dataManager.getData(m_state[iID], PeridigmField::STEP_NONE)->ExtractView(&stat);
        for(int jID=0 ; jID<tensorLen*numOwnedPoints; ++jID){
          statev[iID*tensorLen*numOwnedPoints + jID] = stat[jID];
        }
      }
  }

  // CORRESPONDENCE::computeGreenLagrangeStrain(defGradNP1,GLStrainNP1,flyingPointFlag,numOwnedPoints);
  // dataManager.getData(m_strainFieldId, PeridigmField::STEP_NP1)->ExtractView(&GLStrainNP1);
  double* props = new double[nprops]; //temp;
  for(int iID=0 ; iID<nprops ; ++iID){props[iID] = userProperties[iID];}
  CORRESPONDENCE::userMaterialInterface(modelCoordinates,
                                        defGradN, 
                                        defGradNP1, 
                                        GLStrainN,
                                        GLStrainNP1,
                                        CauchyStressN,
                                        CauchyStressNP1,
                                        numOwnedPoints,
                                        nstatev,
                                        statev,
                                        nprops,
                                        props,
                                        angles,
                                        time,
                                        dt,
                                        temperature,
                                        dtemperature,
                                        RotationN,
                                        RotationNP1,
                                        m_planeStress,
                                        m_planeStrain,
                                        matName);

   if (nstat > 0) {
      double *stat;
      int tensorLen = 9;
      for(int iID=0 ; iID<nstat ; ++iID)
      { 
        dataManager.getData(m_state[iID], PeridigmField::STEP_NONE)->ExtractView(&stat);
        for(int jID=0 ; jID<tensorLen*numOwnedPoints; ++jID){
          stat[jID] = statev[iID*tensorLen*numOwnedPoints + jID];
        }
      }
  }                                           
                                      
                                            
                                           
}
