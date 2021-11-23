/*! \file utPeridigm_Correspondence.cpp */

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
#include "Peridigm_ElasticMaterial.hpp"
#include "Peridigm_SerialMatrix.hpp"
#include "Peridigm_Field.hpp"
#include <Epetra_SerialComm.h>
#include <iostream>


using namespace std;
using namespace PeridigmNS;
using namespace Teuchos;

//! Tests state variable count and name accessor functions.

TEUCHOS_UNIT_TEST(CorrespondenceMaterial, testStateVariableAccessors) {

  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  params.set("Horizon", 10.0);
  ElasticMaterial mat(params);

  // \todo check field specs
}

//! Tests two-point system under compression against hand calculations.

TEUCHOS_UNIT_TEST(ElasticMaterial, testTwoPts) {

  // instantiate the material model
  ParameterList params;
  params.set("Density", 7800.0);
  params.set("Bulk Modulus", 130.0e9);
  params.set("Shear Modulus", 78.0e9);
  params.set("Horizon", 10.0);
  ElasticMaterial mat(params);

  // arguments for calls to material model
  Epetra_SerialComm comm;
  Epetra_Map nodeMap(2, 0, comm);
  Epetra_Map unknownMap(6, 0, comm);
  Epetra_Map bondMap(2, 0, comm);
  double dt = 1.0;
  int numOwnedPoints;
  int* ownedIDs;
  int* neighborhoodList;
  // \todo check field specs

  // set up discretization
  numOwnedPoints = 2;
  ownedIDs = new int[numOwnedPoints];
  ownedIDs[0] = 0;
  ownedIDs[1] = 1;
  neighborhoodList = new int[4];
  neighborhoodList[0] = 1;
  neighborhoodList[1] = 1;
  neighborhoodList[2] = 1;
  neighborhoodList[3] = 0;

  // create the data manager
  // in serial, the overlap and non-overlap maps are the same
  PeridigmNS::DataManager dataManager;
  dataManager.setMaps(Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&nodeMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&unknownMap, false),
                      Teuchos::rcp(&bondMap, false));
  dataManager.allocateData(mat.FieldIds());

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  int modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  int coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  int volumeFieldId = fieldManager.getFieldId("Volume");
  int weightedVolumeFieldId = fieldManager.getFieldId("Weighted_Volume");
  int dilatationFieldId = fieldManager.getFieldId("Dilatation");
  int bondDamageFieldId = fieldManager.getFieldId("Bond_Damage");
  int forceDensityFieldId = fieldManager.getFieldId("Force_Density");

  Epetra_Vector& x = *dataManager.getData(modelCoordinatesFieldId, PeridigmField::STEP_NONE);
  Epetra_Vector& y = *dataManager.getData(coordinatesFieldId, PeridigmField::STEP_NP1);
  Epetra_Vector& cellVolume = *dataManager.getData(volumeFieldId, PeridigmField::STEP_NONE);

  x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
  x[3] = 1.0; x[4] = 0.0; x[5] = 0.0;
  y[0] = 0.0; y[1] = 0.0; y[2] = 0.0;
  y[3] = 2.0; y[4] = 0.0; y[5] = 0.0;
  for(int i=0; i<cellVolume.MyLength(); ++i){
	cellVolume[i] = 1.0;
  }

  mat.initialize(dt, 
                 numOwnedPoints,
                 ownedIDs,
                 neighborhoodList,
                 dataManager);

  mat.computeForce(dt, 
				   numOwnedPoints,
				   ownedIDs,
				   neighborhoodList,
                   dataManager);

  double currentPositionX1 = y[0];
  TEST_COMPARE(currentPositionX1, <=, 1.0e-14);
  double currentPositionY1 = y[1];
  TEST_COMPARE(currentPositionY1, <=, 1.0e-14);
  double currentPositionZ1 = y[2];
  TEST_COMPARE(currentPositionZ1, <=, 1.0e-14);
  double currentPositionX2 = y[3];
  TEST_FLOATING_EQUALITY(currentPositionX2, 2.0, 1.0e-12);
  double currentPositionY2 = y[4];
  TEST_COMPARE(currentPositionY2, <=, 1.0e-14);
  double currentPositionZ2 = y[5];
  TEST_COMPARE(currentPositionZ2, <=,  1.0e-14);
  Epetra_Vector& weightedVolume = *dataManager.getData(weightedVolumeFieldId, PeridigmField::STEP_NONE);
  TEST_FLOATING_EQUALITY(weightedVolume[0], 1.0, 1.0e-15);
  TEST_FLOATING_EQUALITY(weightedVolume[1], 1.0, 1.0e-15);
  Epetra_Vector& dilatation = *dataManager.getData(dilatationFieldId, PeridigmField::STEP_NP1);
  TEST_FLOATING_EQUALITY(dilatation[0], 3.0, 1.0e-15);
  TEST_FLOATING_EQUALITY(dilatation[1], 3.0, 1.0e-15);
  Epetra_Vector& bondDamage = *dataManager.getData(bondDamageFieldId, PeridigmField::STEP_NP1);
  TEST_COMPARE(bondDamage[0], <=, 1.0e-15);
  TEST_COMPARE(bondDamage[1], <=, 1.0e-15);

  Epetra_Vector& force = *dataManager.getData(forceDensityFieldId, PeridigmField::STEP_NP1);
  TEST_FLOATING_EQUALITY(force[0], 2.34e+12, 1.0e-2);
  TEST_COMPARE(force[1], <=, 1.0e-14);
  TEST_COMPARE(force[2], <=, 1.0e-14);
  TEST_FLOATING_EQUALITY(force[3], -2.34e+12, 1.0e-2);
  TEST_COMPARE(force[4], <=, 1.0e-14);
  TEST_COMPARE(force[5], <=, 1.0e-14);

  delete[] ownedIDs;
  delete[] neighborhoodList;
}



int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
