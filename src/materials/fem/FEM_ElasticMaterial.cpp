/*! \file FEM_ElasticMaterial.cpp */

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
#include "FEM_ElasticMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_FEM.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::FEMElasticMaterial::FEMElasticMaterial(const Teuchos::ParameterList& params)
  : FEM(params),
    m_modelAnglesId(-1),
    m_deformationGradientFieldId(-1),
    m_cauchyStressFieldId(-1),
    m_displacementFieldId(-1),
    m_forceDensityFieldId(-1)
{
  bool m_planeStrain = false, m_planeStress = false;
  m_type = 0;
  m_density = params.get<double>("Density");
  order[0] = params.get<int>("Order");
  order[1] = params.get<int>("Order");
  order[2] = params.get<int>("Order");
  if (params.isParameter("Order_Y")) order[1] = params.get<int>("Order_Y");
  if (params.isParameter("Order_Z")) order[2] = params.get<int>("Order_Z");
  if (params.isParameter("Plane Strain")){
    m_planeStrain = params.get<bool>("Plane Strain");
    }
  if (params.isParameter("Plane Stress")){
    m_planeStress = params.get<bool>("Plane Stress");
    }
  if (m_planeStrain==true)m_type=1;
  if (m_planeStress==true)m_type=2;

 
  double C11=0.0, C44=0.0, C55=0.0, C66=0.0, C12=0.0, C13=0.0, C14=0.0, C15=0.0, C16=0.0, C22=0.0, C33=0.0;
  double  C23=0.0, C24=0.0, C25=0.0, C26=0.0, C34=0.0, C35=0.0, C36=0.0, C45=0.0, C46=0.0, C56=0.0;
  bool iso = false;
  if (params.isParameter("Material Symmetry")){
    if (params.get<string>("Material Symmetry")=="Isotropic"){
       iso = true;
       C11 = params.get<double>("C11");
       C44 = params.get<double>("C44");
       C55 = params.get<double>("C44");
       C66 = params.get<double>("C44");
       C12 = C11 - 2*C55;
       C13 = C12;
       C14 = 0.0;
       C15 = 0.0;
       C16 = 0.0;
       C22 = params.get<double>("C11");
       C33 = params.get<double>("C11");
       C23 = C12;
       C24 = 0.0;
       C25 = 0.0;
       C26 = 0.0;
       C34 = 0.0;
       C35 = 0.0;
       C36 = 0.0;
       C45 = 0.0;
       C46 = 0.0;
       C56 = 0.0; 
    }
   if (params.get<string>("Material Symmetry")=="Anisotropic"){
       C11 = params.get<double>("C11");
       C12 = params.get<double>("C12");
       C13 = params.get<double>("C13");
       C14 = params.get<double>("C14");
       C15 = params.get<double>("C15");
       C16 = params.get<double>("C16");
       C22 = params.get<double>("C22");
       C23 = params.get<double>("C23");
       C24 = params.get<double>("C24");
       C25 = params.get<double>("C25");
       C26 = params.get<double>("C26");
       C33 = params.get<double>("C33");
       C34 = params.get<double>("C34");
       C35 = params.get<double>("C35");
       C36 = params.get<double>("C36");
       C44 = params.get<double>("C44");
       C45 = params.get<double>("C45");
       C46 = params.get<double>("C46");
       C55 = params.get<double>("C55");
       C56 = params.get<double>("C56");
       C66 = params.get<double>("C66");
       //std::cout<<"??"<<std::endl;
   }
  }
  else{
         iso = true;
         m_bulkModulus = calculateBulkModulus(params);
         m_shearModulus = calculateShearModulus(params);
         double lam = m_bulkModulus - 2./3.*m_shearModulus; 
         C11 = 2.*m_shearModulus + lam;
         C44 = m_shearModulus;
         C55 = m_shearModulus;
         C66 = m_shearModulus;
         C12 = C11 - 2.*C55;
         C13 = C12;
         C14 = 0.0;
         C15 = 0.0;
         C16 = 0.0;
         C22 = 2.*m_shearModulus + lam;
         C33 = 2.*m_shearModulus + lam;
         C23 = C12;
         C24 = 0.0;
         C25 = 0.0;
         C26 = 0.0;
         C34 = 0.0;
         C35 = 0.0;
         C36 = 0.0;
         C45 = 0.0;
         C46 = 0.0;
         C56 = 0.0;
  }
// Equation (8) Dipasquale, D., Sarego, G., Zaccariotto, M., Galvanetto, U., A discussion on failure criteria
  // for ordinary state-based Peridynamics, Engineering Fracture Mechanics (2017), doi: https://doi.org/10.1016/
  // j.engfracmech.2017.10.011
  if (m_planeStrain==true)m_plane=true;
  if (m_planeStress==true)m_plane=true;
  if (m_plane==false){
   C[0][0] = C11;C[0][1] = C12;C[0][2]= C13; C[0][3] = C14; C[0][4] = C15; C[0][5]= C16;
   C[1][0] = C12;C[1][1] = C22;C[1][2]= C23; C[1][3] = C24; C[1][4] = C25; C[1][5]= C26;
   C[2][0] = C13;C[2][1] = C23;C[2][2]= C33; C[2][3] = C34; C[2][4] = C35; C[2][5]= C36;
   C[3][0] = C14;C[3][1] = C24;C[3][2]= C34; C[3][3] = C44; C[3][4] = C45; C[3][5]= C46;
   C[4][0] = C15;C[4][1] = C25;C[4][2]= C35; C[4][3] = C45; C[4][4] = C55; C[4][5]= C56;
   C[5][0] = C16;C[5][1] = C26;C[5][2]= C36; C[5][3] = C46; C[5][4] = C56; C[5][5]= C66;
  }
  else{
    if (m_planeStress==true && iso == true){
           //only transversal isotropic in the moment --> definition of iso missing
        C[0][0] = C11-C13*C13/C22;C[0][1] = C12-C13*C23/C22;C[0][2] = 0.0; C[0][3] = 0.0; C[0][4] = 0.0; C[0][5] = 0.0;
        C[1][0] = C12-C13*C23/C22;C[1][1] = C22-C13*C23/C22;C[1][2] = 0.0; C[1][3] = 0.0; C[1][4] = 0.0; C[1][5] = 0.0;
        C[2][0] = 0.0;            C[2][1] = 0.0;            C[2][2] = 0.0; C[2][3] = 0.0; C[2][4] = 0.0; C[2][5] = 0.0;
        C[3][0] = 0.0;            C[3][1] = 0.0;            C[3][2] = 0.0; C[3][3] = 0.0; C[3][4] = 0.0; C[3][5] = 0.0;
        C[4][0] = 0.0;            C[4][1] = 0.0;            C[4][2] = 0.0; C[4][3] = 0.0; C[4][4] = 0.0; C[4][5] = 0.0;
        C[5][0] = 0.0;            C[5][1] = 0.0;            C[5][2] = 0.0; C[5][3] = 0.0; C[5][4] = 0.0; C[5][5] = C66;
    }
    else{
        
         C[0][0] = C11;C[0][1] = C12;C[0][2] = 0.0;C[0][3] = 0.0;C[0][4] = 0.0;C[0][5] = C16;
         C[1][0] = C12;C[1][1] = C22;C[1][2] = 0.0;C[1][3] = 0.0;C[1][4] = 0.0;C[1][5] = C26;
         C[2][0] = 0.0;C[2][1] = 0.0;C[2][2] = 0.0;C[2][3] = 0.0;C[2][4] = 0.0;C[2][5] = 0.0;
         C[3][0] = 0.0;C[3][1] = 0.0;C[3][2] = 0.0;C[3][3] = 0.0;C[3][4] = 0.0;C[3][5] = 0.0;
         C[4][0] = 0.0;C[4][1] = 0.0;C[4][2] = 0.0;C[4][3] = 0.0;C[4][4] = 0.0;C[4][5] = 0.0;
         C[5][0] = C16;C[5][1] = C26;C[5][2] = 0.0;C[5][3] = 0.0;C[5][4] = 0.0;C[5][5] = C66;
        }
    }
    
  
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();

  m_cauchyStressFieldId                 = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_modelCoordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::CONSTANT, "Model_Coordinates");
  m_displacementFieldId                 = fieldManager.getFieldId(PeridigmField::NODE   , PeridigmField::VECTOR, PeridigmField::TWO_STEP     , "Displacements");
  m_coordinatesFieldId                  = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Coordinates");
  m_modelAnglesId                       = fieldManager.getFieldId(PeridigmField::NODE   , PeridigmField::VECTOR, PeridigmField::CONSTANT     , "Local_Angles");
  m_forceDensityFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_fieldIds.push_back(m_cauchyStressFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_displacementFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_modelAnglesId);
  m_fieldIds.push_back(m_forceDensityFieldId);

}

PeridigmNS::FEMElasticMaterial::~FEMElasticMaterial()
{
}

void
PeridigmNS::FEMElasticMaterial::initialize(const double dt,
                                                             const int numOwnedPoints,
                                                             const int* ownedIDs,
                                                             const int* neighborhoodList,
                                                             PeridigmNS::DataManager& dataManager)
{
      PeridigmNS::FEMElasticMaterial::initialize(dt,
                                                      numOwnedPoints,
                                                      ownedIDs,
                                                      neighborhoodList,
                                                      dataManager);

                             
}

void
PeridigmNS::FEMElasticMaterial::computeCauchyStress(const double dt,
                                                               const int numOwnedPoints,
                                                               const int* neighborhoodList,
                                                               PeridigmNS::DataManager& dataManager) const
{

  

  double *CauchyStressNP1, *coor, *angles, *displacements, *deformedCoor, *force;
 

  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&CauchyStressNP1);

  dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&angles);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&coor);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&deformedCoor);
  dataManager.getData(m_displacementFieldId, PeridigmField::STEP_NP1)->ExtractView(&displacements);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);
  FEM::getDisplacements(numOwnedPoints,coor,deformedCoor,displacements);
  int numElements;
  int *elementNodalList;
  FEM::elasticFEM(coor, 
                  displacements,
                  CauchyStressNP1,
                  numElements,
                  elementNodalList
                  C,
                  angles,
                  m_type,
                  dt,
                  order,
                  force); // --> ruft updateElasticCauchyStressAnisotropicCode auf
                                            
                                            
                                     
                                            
                                           
}
