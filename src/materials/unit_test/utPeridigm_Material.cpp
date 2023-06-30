/*! \file utPeridigm_Material.cpp */

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

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "Peridigm_ElasticLinearCorrespondenceMaterial.hpp"
#include "Peridigm_Material.hpp"
#include "Peridigm_SerialMatrix.hpp"
#include "Peridigm_Field.hpp"
#include <Epetra_SerialComm.h>
#include <iostream>


using namespace std;
using namespace PeridigmNS;
using namespace Teuchos;

TEUCHOS_UNIT_TEST(ElasticLinearCorrespondenceMaterial, testStateVariableAccessors) {

  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  params.set("C11", 10.0);
  params.set("C12", 10.0);
  params.set("C13", 10.0);
  params.set("C14", 10.0);
  params.set("C15", 10.0);
  params.set("C16", 10.0);
  params.set("C22", 10.0);
  params.set("C23", 10.0);
  params.set("C24", 10.0);
  params.set("C25", 10.0);
  params.set("C26", 10.0);
  params.set("C33", 10.0);
  params.set("C34", 10.0);
  params.set("C35", 10.0);
  params.set("C36", 10.0);
  params.set("C44", 10.0);
  params.set("C45", 10.0);
  params.set("C46", 10.0);
  params.set("C55", 10.0);
  params.set("C56", 10.0);
  params.set("C66", 10.0);

  ElasticLinearCorrespondenceMaterial mat(params);
  // \todo check field specs
}

TEUCHOS_UNIT_TEST(ElasticLinearCorrespondenceMaterial, getStiffnessmatrix) {
 // instantiate the material model
  double tolerance = 1e-8;
  ParameterList params;
  double m_bulkModulus = 130.0e9;
  double m_shearModulus = 78.0e9;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", m_bulkModulus);
  params.set("Shear Modulus", m_shearModulus);
  ElasticLinearCorrespondenceMaterial mat(params);

  double EMod = 9*m_bulkModulus * m_shearModulus / (3*m_bulkModulus + m_shearModulus);
  double nu = (3*m_bulkModulus - 2*m_shearModulus) / (6*m_bulkModulus + 2*m_shearModulus);

  double C[6][6];
  int i,j;
  double cTest[6][6];
  double cTestps[6][6];
  for(i=0; i<6; ++i){
    for (j=0; j<6; ++j){
        cTest[i][j] = 0.0;
        cTestps[i][j] = 0.0;
    }
  }
  double temp = 1 / ((1+nu) * (1-2*nu));
  cTest[0][0] = EMod  * (1-nu) * temp;
  cTest[1][1] = EMod  * (1-nu) * temp;
  cTest[2][2] = EMod  * (1-nu) * temp;
  cTest[0][1] = EMod * nu * temp;
  cTest[0][2] = EMod * nu * temp;
  cTest[1][0] = EMod * nu * temp;
  cTest[2][0] = EMod * nu * temp;
  cTest[2][1] = EMod * nu * temp;
  cTest[1][2] = EMod * nu * temp;
  cTest[3][3] = m_shearModulus;
  cTest[4][4] = m_shearModulus;
  cTest[5][5] = m_shearModulus;


  mat.getStiffnessmatrix(params,C,false,false);
  for(i=0; i<6; ++i){
    for (j=0; j<6; ++j){
        TEST_FLOATING_EQUALITY(C[i][j], cTest[i][j], tolerance);
    }
  }
  double Eps = cTest[5][5]*(3.0*cTest[0][0] - 4.0*cTest[5][5])/(cTest[0][0] - cTest[5][5]);     
  double nups = (Eps/2.0 - cTest[5][5])/cTest[5][5];
  
  //https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm

  cTestps[0][0] = Eps / (1-nups*nups);
  cTestps[0][1] = Eps*nups / (1-nups*nups);
  cTestps[1][0] = Eps*nups / (1-nups*nups);
  cTestps[1][1] = Eps / (1-nups*nups);
  cTestps[5][5] = cTest[5][5];

  mat.getStiffnessmatrix(params,C,false,true);
  for(i=0; i<6; ++i){
    for (j=0; j<6; ++j){
        TEST_FLOATING_EQUALITY(C[i][j], cTestps[i][j], tolerance);
    }
  }

  // test iso3D
  params.set("Material Symmetry", "Isotropic");
  mat.getStiffnessmatrix(params,C,false,false);
  for(i=0; i<6; ++i){
    for (j=0; j<6; ++j){
        TEST_FLOATING_EQUALITY(C[i][j], 0.0, tolerance);
    }
  }

  for(i=0; i<6; ++i){
    for (j=0; j<6; ++j){
      params.set("C" + std::to_string(i+1) + std::to_string(j+1), cTest[i][j]);
    }
  }
  mat.getStiffnessmatrix(params,C,false,false);
  for(i=0; i<6; ++i){
    for (j=0; j<6; ++j){
        TEST_FLOATING_EQUALITY(C[i][j], cTest[i][j], tolerance);
    }
  }
}

TEUCHOS_UNIT_TEST(ElasticLinearCorrespondenceMaterial, getThermalExpansionCoefficient) {
 // instantiate the material model
  double tolerance = 1e-8;
  ParameterList params;
  double m_bulkModulus = 130.0e9;
  double m_shearModulus = 78.0e9;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", m_bulkModulus);
  params.set("Shear Modulus", m_shearModulus);

  params.set("Apply Thermal Strain",false);

  ElasticLinearCorrespondenceMaterial mat(params);

  double alpha[3][3];
  double m_Tref = 10;
  int i,j;
  double alphaTest[3][3];
  double m_TrefTest = 10;

  for(i=0; i<3; ++i){
    alphaTest[i][i] = i/3 + 2;
    for (j=0; j<3; ++j){
        alphaTest[i][j] = 0.0;
    }
  }

  TEST_ASSERT(mat.getThermalExpansionCoefficient(params,alpha,m_Tref)==false);
  params.set("Apply Thermal Strain",true);
  params.set("Thermal Expansion Coefficient", alphaTest[0][0]);
  params.set("Thermal Expansion Reference Termperature", m_TrefTest);
  TEST_ASSERT(mat.getThermalExpansionCoefficient(params,alpha,m_Tref)==true);
  for(i=0; i<3; ++i){
    TEST_FLOATING_EQUALITY(alpha[i][i], alphaTest[0][0], tolerance);
    for (j=i+1; j<3; ++j){
        TEST_FLOATING_EQUALITY(alpha[i][j], alphaTest[i][j], tolerance);
        TEST_FLOATING_EQUALITY(alpha[j][i], alphaTest[j][i], tolerance);
    }
  }
  params.set("Thermal Expansion Coefficient Y", alphaTest[1][1]);
  params.set("Thermal Expansion Coefficient Z", alphaTest[2][2]);
  for(i=0; i<3; ++i){
    for (j=0; j<3; ++j){
        TEST_FLOATING_EQUALITY(alpha[i][j], alphaTest[i][j], tolerance);
        
    }
  }
  TEST_FLOATING_EQUALITY(m_Tref, m_TrefTest, tolerance);

}
int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
