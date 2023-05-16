/*! \file Peridigm_CorrespondenceMaterial.cpp */

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

#include "Peridigm_CorrespondenceMaterial.hpp"
#include "Peridigm_Logging.hpp"
#include "Peridigm_Field.hpp"
#include "correspondence.h"
#include "matrices.h"
#include "Peridigm_DegreesOfFreedomManager.hpp"
#include "Peridigm_Logging.hpp"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include "approximation.h"
using namespace std;

PeridigmNS::CorrespondenceMaterial::CorrespondenceMaterial(const Teuchos::ParameterList &params)
    : Material(params),
      m_density(0.0), m_C(0.0), m_hourglassCoefficient(0.0),
      m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
      m_horizonFieldId(-1), m_volumeFieldId(-1),
      m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_velocitiesFieldId(-1),
      m_hourglassForceDensityFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1),
      m_deformationGradientFieldId(-1),
      m_shapeTensorInverseFieldId(-1),
      m_cauchyStressFieldId(-1),
      m_leftStretchTensorFieldId(-1),
      m_rotationTensorFieldId(-1),
      m_unrotatedCauchyStressFieldId(-1),
      m_unrotatedRateOfDeformationFieldId(-1),
      m_unrotatedCauchyStressPlasticFieldId(-1),
      m_partialStressFieldId(-1),
      m_hourglassStiffId(-1),
      m_thermalFlowStateFieldId(-1),
      m_jacobianId(-1)
{

  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");

  m_stabilizationType = 1;
  m_bondbased = true;
  m_plane = false;
  m_approximation = false; 
  nonLin = true;
  linRateOfDeformation = true;
  m_plast = false;
  if (params.isParameter("Linear Elastic Correspondence"))
  {
    nonLin = false;
    linRateOfDeformation = false;
  }
  if (params.isParameter("Non linear"))
  {
    nonLin = params.get<bool>("Non linear");
  }

  bool m_planeStrain = false, m_planeStress = false;
  if (params.isParameter("Plane Strain"))
    m_planeStrain = params.get<bool>("Plane Strain");

  if (params.isParameter("Plane Stress"))
    m_planeStress = params.get<bool>("Plane Stress");
  if (m_planeStrain == true)
  {
    m_plane = true;
  }
  if (m_planeStress == true)
  {
    m_plane = true;
  }

  // if (params.get<std::string>("Material Model").find("Plastic")!=std::string::npos) {
  //     m_plast = true; // does not work. we have to check
  //
  // }

 if (params.isParameter("B-Spline Approximation")){
    m_approximation = params.get<bool>("B-Spline Approximation");
    if (m_approximation){
      m_num_control_points = params.get<int>("Control Points");
      m_degree = params.get<int>("Degree");
      if (m_num_control_points<m_degree+1)m_num_control_points = m_degree+1;
    }
  }
 if (params.isParameter("Thermal Bond Based"))
    m_bondbased = params.get<bool>("Thermal Bond Based");
  
  m_applyAutomaticDifferentiationJacobian = false;
  if (params.isParameter("Accumulated Plastic"))
    m_plast = params.get<bool>("Accumulated Plastic");

  if (params.isParameter("Apply Automatic Differentiation Jacobian"))
    m_applyAutomaticDifferentiationJacobian = params.get<bool>("Apply Automatic Differentiation Jacobian");
  if (params.isParameter("Stabilization Type"))
  {

    if (params.get<string>("Stabilization Type") == "Bond Based")
    {
      m_stabilizationType = 1;
      m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");
    }
    if (params.get<string>("Stabilization Type") == "State Based")
    {
      m_stabilizationType = 2;
      m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");
    }

    m_adaptHourGlass = false;
    if (params.get<string>("Stabilization Type") == "Global Stiffness")
    {
      getStiffnessmatrix(params, C, m_planeStrain, m_planeStress);

      m_stabilizationType = 3;
      m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");
      if (params.isParameter("Adapt Hourglass Stiffness"))
      {
        m_adaptHourGlass = params.get<bool>("Adapt Hourglass Stiffness");
      }
    }
  }

  m_applyThermalFlow = false;
  m_applyThermalPrintBedFlow = false;
  if (params.isParameter("Apply Thermal Flow")){
    m_applyThermalFlow = params.get<bool>("Apply Thermal Flow");
    if (m_applyThermalFlow){
      m_C = params.get<double>( "Specific Heat Capacity");
      m_lambda.push_back(params.get<double>("Thermal Conductivity"));
      m_lambda.push_back(0.0);
      m_lambda.push_back(0.0);
      m_lambda.push_back(0.0);
      if (params.isParameter("Thermal Conductivity 22"))m_lambda.push_back(params.get<double>("Thermal Conductivity 22"));
      else m_lambda.push_back(params.get<double>("Thermal Conductivity")) ;
      m_lambda.push_back(0.0);
      m_lambda.push_back(0.0);
      m_lambda.push_back(0.0);
      if (params.isParameter("Thermal Conductivity 33"))m_lambda.push_back(params.get<double>("Thermal Conductivity 33"));
      else m_lambda.push_back(params.get<double>("Thermal Conductivity")) ;
      
      materialProperties["Specific Heat Capacity"] = m_C;
      materialProperties["Thermal Conductivity 11"] = m_lambda[0];
      materialProperties["Thermal Conductivity 22"] = m_lambda[4];
      materialProperties["Thermal Conductivity 33"] = m_lambda[8];

      if (params.isParameter("Thermal Conductivity Print Bed") && params.isParameter("Print Bed Temperature")){
        m_applyThermalPrintBedFlow = true;
        m_lambdaBed = params.get<double>("Thermal Conductivity Print Bed");
        m_Tbed = params.get<double>("Print Bed Temperature");
      }
    }

    m_applyHeatTransfer = false;
    if (params.isParameter("Apply Heat Transfer")){
      m_applyHeatTransfer = params.get<bool>("Apply Heat Transfer");
      if (m_applyHeatTransfer){
        m_kappa = params.get<double>("Heat Transfer Coefficient");
        m_Tenv = params.get<double>("Environmental Temperature");
        m_factor = 1.0;
        m_surfaceCorrection = 1.0;
        m_limit = 0.8;
        if (params.isParameter("Volume Factor")) m_factor= params.get<double>("Volume Factor");
        if (params.isParameter("Surface Correction")) m_surfaceCorrection= params.get<double>("Surface Correction");
        if (params.isParameter("Volume Limit")) m_limit= params.get<double>("Volume Limit");
        
      }
    }
  }
  // TestForTermination(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the ElasticCorrespondence material model.\n");
  TestForTermination(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the ElasticCorrespondence material model.\n");

  PeridigmNS::FieldManager &fieldManager = PeridigmNS::FieldManager::self();
  m_horizonFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_volumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_modelAnglesId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Local_Angles");
  m_modelOrientationId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Orientations");
  m_coordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_velocitiesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity");
  m_forceDensityFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_hourglassForceDensityFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Hourglass_Force_Density");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_deformationGradientFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Deformation_Gradient");
  m_leftStretchTensorFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor");
  m_rotationTensorFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Rotation_Tensor");
  m_shapeTensorInverseFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Shape_Tensor_Inverse");
  m_unrotatedCauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_unrotatedCauchyStressPlasticFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Plastic_Cauchy_Stress");
  m_cauchyStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Cauchy_Stress");
  m_piolaStressTimesInvShapeTensorXId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorX");
  m_piolaStressTimesInvShapeTensorYId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorY");
  m_piolaStressTimesInvShapeTensorZId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorZ");
  m_partialStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Partial_Stress");
  m_detachedNodesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Detached_Nodes");
  m_specificVolumeFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Specific_Volume");
  
  m_hourglassStiffId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Hourglass_Stiffness");
  if (m_applyThermalFlow){
    
    m_temperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature");
    
    m_fieldIds.push_back(m_temperatureFieldId);
  }
  m_thermalFlowStateFieldId          = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Flux_Divergence"); // reuse of the field
  m_jacobianId = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::FULL_TENSOR,      PeridigmField::CONSTANT, "Jacobian_Matrix");
  m_fieldIds.push_back(m_thermalFlowStateFieldId); 
  m_fieldIds.push_back(m_horizonFieldId);
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_velocitiesFieldId);
  m_fieldIds.push_back(m_hourglassForceDensityFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_deformationGradientFieldId);

  m_fieldIds.push_back(m_leftStretchTensorFieldId);
  m_fieldIds.push_back(m_rotationTensorFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressPlasticFieldId);
  m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorXId);
  m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorYId);
  m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorZId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_partialStressFieldId);
  m_fieldIds.push_back(m_detachedNodesFieldId);
  m_fieldIds.push_back(m_specificVolumeFieldId);
  m_fieldIds.push_back(m_hourglassStiffId);
  m_fieldIds.push_back(m_modelAnglesId);
  m_fieldIds.push_back(m_modelOrientationId);
  m_fieldIds.push_back(m_jacobianId);

}

PeridigmNS::CorrespondenceMaterial::~CorrespondenceMaterial()
{
}

void PeridigmNS::CorrespondenceMaterial::initialize(const double dt,
                                                    const int numOwnedPoints,
                                                    const int *ownedIDs,
                                                    const int *neighborhoodList,
                                                    PeridigmNS::DataManager &dataManager)
{

  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_hourglassStiffId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressPlasticFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *leftStretchTensorN;
  double *leftStretchTensorNP1;
  double *rotationTensorN;
  double *rotationTensorNP1;
  double *detachedNodes;
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
  // Initialize the left stretch and rotation tenor to the identity matrix
  MATRICES::setOnesOnDiagonalFullTensor(leftStretchTensorN, numOwnedPoints);
  MATRICES::setOnesOnDiagonalFullTensor(leftStretchTensorNP1, numOwnedPoints);
  MATRICES::setOnesOnDiagonalFullTensor(rotationTensorN, numOwnedPoints);
  MATRICES::setOnesOnDiagonalFullTensor(rotationTensorNP1, numOwnedPoints);

  // Initialize the inverse of the shape tensor and the deformation gradient
  double *volume;
  double *horizon;
  double *modelCoordinates;
  double *coordinatesNP1;
  double *shapeTensorInverse;
  double *deformationGradient;
  double *bondDamageNP1;
  double *jacobian;
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);

  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradient);
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);
  dataManager.getData(m_jacobianId, PeridigmField::STEP_NONE)->ExtractView(&jacobian);

  int shapeTensorReturnCode = 0;
  shapeTensorReturnCode =
      CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                                 horizon,
                                                                                 modelCoordinates,
                                                                                 coordinatesNP1,
                                                                                 shapeTensorInverse,
                                                                                 deformationGradient,
                                                                                 bondDamageNP1,
                                                                                 neighborhoodList,
                                                                                 numOwnedPoints,
                                                                                 m_plane,
                                                                                 detachedNodes);

  string shapeTensorErrorMessage =
      "**** Error:  CorrespondenceMaterial::initialize() failed to compute shape tensor.\n";
  shapeTensorErrorMessage +=
      "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TestForTermination(shapeTensorReturnCode != 0, shapeTensorErrorMessage);
  
  
  if (m_approximation){
    int len = APPROXIMATION::get_field_size(numOwnedPoints, neighborhoodList, m_num_control_points, m_plane);
    delete approxMatrix;
    approxMatrix = new double[len];
    delete Bsplinegradient;
    Bsplinegradient = new double[len * PeridigmNS::dof()];

    if (m_plane) len = m_num_control_points * m_num_control_points;
    else len = m_num_control_points * m_num_control_points * m_num_control_points;
    double *contP;
    contP = new double[numOwnedPoints * len * PeridigmNS::dof()];
    delete gradient_function;
    gradient_function = new double[len * PeridigmNS::dof()];
    APPROXIMATION::get_approximation(numOwnedPoints, neighborhoodList,modelCoordinates,m_num_control_points,m_degree,m_plane,approxMatrix);
    APPROXIMATION::get_control_points(numOwnedPoints,neighborhoodList,m_num_control_points,modelCoordinates,approxMatrix,m_plane,contP);
    APPROXIMATION::get_gradient_functions(m_degree,m_num_control_points,m_plane,gradient_function);
    APPROXIMATION::get_jacobians(numOwnedPoints,contP,m_num_control_points,gradient_function,m_plane,jacobian);
    APPROXIMATION::get_B_spline_gradient(numOwnedPoints,neighborhoodList,m_num_control_points,approxMatrix,gradient_function,jacobian,
    m_plane,Bsplinegradient);

    
  }



  
  //APPROXIMATION::create_approximation(numOwnedPoints, neighborhoodList, modelCoordinates, num_control_points,degree,true,AMatrix);
  double *pointAngles;
  double *pointOrientations;
  dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&pointAngles);
  dataManager.getData(m_modelOrientationId, PeridigmField::STEP_NONE)->ExtractView(&pointOrientations);

  CORRESPONDENCE::getOrientations(numOwnedPoints,pointAngles,pointOrientations);

}

void PeridigmNS::CorrespondenceMaterial::computeForce(const double dt,
                                                      const int numOwnedPoints,
                                                      const int *ownedIDs,
                                                      const int *neighborhoodList,
                                                      PeridigmNS::DataManager &dataManager,
                                                      const double currentTime) const
{
  // Zero out the forces and partial stress

  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  double *horizon, *volume, *modelCoordinates, *coordinatesN, *coordinatesNP1, *shapeTensorInverse, *deformationGradient, *bondDamage, *bondDamageNP1, *pointAngles, *detachedNodes;
  // double *deformationGradientNonInc;

  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&pointAngles);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_N)->ExtractView(&coordinatesN);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
  dataManager.getData(m_shapeTensorInverseFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradient);

  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamage);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);

  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
  *(dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_N));
  double *unrotatedRateOfDeformation;
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);

  double *velocities, *leftStretchTensorN, *leftStretchTensorNP1, *rotationTensorN, *rotationTensorNP1;
  dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);
  // Compute the inverse of the shape tensor and the approximate deformation gradient
  // The approximate deformation gradient will be used by the derived class (specific correspondence material model)
  // to compute the Cauchy stress.
  // The inverse of the shape tensor is stored for later use after the Cauchy stress calculation
  
  
  if (m_approximation){
    double *jacobian;
    dataManager.getData(m_jacobianId, PeridigmField::STEP_NONE)->ExtractView(&jacobian);
    // int len;
    // if (m_plane) len = m_num_control_points * m_num_control_points;
    // else len = m_num_control_points * m_num_control_points * m_num_control_points;
    // double *contPu;
    // contPu = new double[numOwnedPoints * len * PeridigmNS::dof()];
    
    //APPROXIMATION::get_control_points(numOwnedPoints,neighborhoodList,m_num_control_points,coordinatesNP1,approxMatrix,m_plane,contPu);
    //APPROXIMATION::get_deformation_gradient(numOwnedPoints,contPu,m_num_control_points,gradient_function,m_plane,jacobian,deformationGradient);
    //Bsplinegradient
    APPROXIMATION::get_deformation_gradient_new(numOwnedPoints,neighborhoodList,coordinatesNP1,Bsplinegradient,deformationGradient);
    

  }
  else{
    int shapeTensorReturnCode = 0;
    shapeTensorReturnCode =
        CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                                  horizon,
                                                                                  modelCoordinates,
                                                                                  coordinatesNP1,
                                                                                  shapeTensorInverse,
                                                                                  deformationGradient,
                                                                                  bondDamageNP1,
                                                                                  neighborhoodList,
                                                                                  numOwnedPoints,
                                                                                  m_plane,
                                                                                  detachedNodes);
  
    string shapeTensorErrorMessage =
        "**** Error:  CorrespondenceMaterial::computeForce() failed to compute shape tensor.\n";
    shapeTensorErrorMessage +=
        "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
    TestForTermination(shapeTensorReturnCode != 0, shapeTensorErrorMessage);
  

    if (nonLin == true)
    {

      dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
      dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
      dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
      dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);

      // Compute left stretch tensor, rotation tensor, and unrotated rate-of-deformation.
      // Performs a polar decomposition via Flanagan & Taylor (1987) algorithm.
      //"A non-ordinary state-based peridynamic method to model solid material deformation and fracture"
      // Thomas L. Warren, Stewart A. Silling, Abe Askari, Olaf Weckner, Michael A. Epton, Jifeng Xu

      int rotationTensorReturnCode = 0;
      
      rotationTensorReturnCode = CORRESPONDENCE::computeUnrotatedRateOfDeformationAndRotationTensor(volume,
                                                                                                    horizon,
                                                                                                    modelCoordinates,
                                                                                                    velocities,
                                                                                                    deformationGradient,
                                                                                                    shapeTensorInverse,
                                                                                                    leftStretchTensorN,
                                                                                                    rotationTensorN,
                                                                                                    leftStretchTensorNP1,
                                                                                                    rotationTensorNP1,
                                                                                                    unrotatedRateOfDeformation,
                                                                                                    neighborhoodList,
                                                                                                    numOwnedPoints,
                                                                                                    dt,
                                                                                                    bondDamageNP1,
                                                                                                    m_plane,
                                                                                                    detachedNodes);

      string rotationTensorErrorMessage =
          "**** Error:  CorrespondenceMaterial::computeForce() failed to compute rotation tensor.\n";
      rotationTensorErrorMessage +=
          "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
      string rotationTensorErrorMessage2 =
          "**** Error:  CorrespondenceMaterial::computeForce() failed to invert deformation gradient tensor.\n";

      string rotationTensorErrorMessage3 =
          "**** Error:  CorrespondenceMaterial::computeForce() failed to invert temp.\n";

      TestForTermination(rotationTensorReturnCode == 1, rotationTensorErrorMessage);
      TestForTermination(rotationTensorReturnCode == 2, rotationTensorErrorMessage2);
      TestForTermination(rotationTensorReturnCode == 3, rotationTensorErrorMessage3);
    }
    else
    {

      // needed to create the rate of deformation for the step wise calculation
      // all rotations are excluded; this speeds up the calculations, but assumes no large rotations
      // is deactivated for linear elastic material, because the rate is not used in that case
      if (linRateOfDeformation)
      {
        CORRESPONDENCE::getLinearUnrotatedRateOfDeformation(volume,
                                                            horizon,
                                                            modelCoordinates,
                                                            velocities,
                                                            deformationGradient,
                                                            shapeTensorInverse,
                                                            unrotatedRateOfDeformation,
                                                            neighborhoodList,
                                                            numOwnedPoints,
                                                            bondDamageNP1,
                                                            m_plane,
                                                            detachedNodes);
      }
    }
  }


  if (m_applyThermalFlow){
    double *thermalFlow, *temperature; // quadratureWeights not included yet
    
    dataManager.getData(m_temperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&temperature);
    dataManager.getData(m_thermalFlowStateFieldId, PeridigmField::STEP_NP1)->ExtractView(&thermalFlow);
    dataManager.getData(m_thermalFlowStateFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    const double* lambda = &m_lambda[0];
    CORRESPONDENCE::computeHeatFlowState(
                                  modelCoordinates,
                                  numOwnedPoints,
                                  neighborhoodList,
                                  shapeTensorInverse,
                                  temperature,
                                  horizon,
                                  lambda,
                                  volume,
                                  detachedNodes,
                                  bondDamageNP1,
                                  pointAngles,
                                  m_plane,
                                  m_bondbased,
                                  m_applyThermalPrintBedFlow,
                                  m_lambdaBed,
                                  m_Tbed,
                                  thermalFlow);
    if (m_applyHeatTransfer){
      double *specificVolume;
      dataManager.getData(m_specificVolumeFieldId, PeridigmField::STEP_NP1)->ExtractView(&specificVolume);
  
      CORRESPONDENCE::computeHeatTransfer(
                                  modelCoordinates,
                                  numOwnedPoints,
                                  neighborhoodList,
                                  volume,
                                  temperature,
                                  horizon,
                                  detachedNodes,
                                  bondDamageNP1,
                                  m_plane,
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
  // Evaluate the Cauchy stress using the routine implemented in the derived class (specific correspondence material model)
  // The general idea is to compute the stress based on:
  //   1) The unrotated rate-of-deformation tensor
  //   2) The time step
  //   3) Whatever state variables are managed by the derived class
  //
  // computeCauchyStress() typically uses the following fields which are accessed via the DataManager:
  //   Input:  unrotated rate-of-deformation tensor
  //   Input:  unrotated Cauchy stress at step N
  //   Input:  internal state data (managed in the derived class)
  //   Output: unrotated Cauchy stress at step N+1

  // multiple Cauchy stresses will be provided over the datamanager --> Peridigm_ElasticLinearCorrespondence

  computeCauchyStress(dt, numOwnedPoints, dataManager, currentTime);

  // rotate back to the Eulerian frame
  double *unrotatedCauchyStressNP1, *cauchyStressNP1, *plasticStress;
  if (m_plast)
  {
    dataManager.getData(m_unrotatedCauchyStressPlasticFieldId, PeridigmField::STEP_NONE)->ExtractView(&plasticStress);
  }
  if (nonLin == true)
  {
    dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);
    dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1);

    CORRESPONDENCE::rotateCauchyStress(rotationTensorNP1,
                                       unrotatedCauchyStressNP1,
                                       cauchyStressNP1,
                                       numOwnedPoints);
    // Cauchy stress is now updated and in the rotated state.

    if (m_plast)
    {
      CORRESPONDENCE::rotateCauchyStress(rotationTensorNP1,
                                         plasticStress,
                                         plasticStress,
                                         numOwnedPoints);
    }
  }
  else
  {
    dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1);
    dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);
    *(dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1));
  }

  // }
  // std::cout<<*(cauchyStressNP1+1)<<" corr"<<std::endl;
  // std::cout<<*(cauchyStressNP1+1)<<" corrCau"<<std::endl;

  // Proceed with conversion to Piola-Kirchoff and force-vector states.
  //--------------------------------------------------------------------------
  double *forceDensity;
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  double *partialStress;
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&partialStress);
  dataManager.getData(m_partialStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *tempStressX, *tempStressY, *tempStressZ;
  dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->ExtractView(&tempStressX);
  dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->ExtractView(&tempStressY);
  dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->ExtractView(&tempStressZ);
  dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *hourglassStiff;

  dataManager.getData(m_hourglassStiffId, PeridigmField::STEP_NONE)->ExtractView(&hourglassStiff);
  
  CORRESPONDENCE::computeForcesAndStresses(
      numOwnedPoints,
      neighborhoodList,
      volume,
      horizon,
      modelCoordinates,
      coordinatesNP1,
      deformationGradient,
      cauchyStressNP1,
      plasticStress,
      shapeTensorInverse,
      bondDamageNP1,
      C,
      pointAngles,
      forceDensity,
      partialStress,
      tempStressX,
      tempStressY,
      tempStressZ,
      hourglassStiff,
      m_hourglassCoefficient,
      m_stabilizationType,
      m_plane,
      m_plast,
      m_adaptHourGlass,
      detachedNodes);

  //     std::cout<<numOwnedPoints<< " "<< *(deformationGradient)<<" "<<*(partialStress)<<std::endl;

  if (m_stabilizationType == 1)
  {
    CORRESPONDENCE::computeHourglassForce(volume,
                                          horizon,
                                          modelCoordinates,
                                          coordinatesNP1,
                                          deformationGradient,
                                          forceDensity,
                                          neighborhoodList,
                                          numOwnedPoints,
                                          m_bulkModulus,
                                          m_hourglassCoefficient,
                                          bondDamageNP1);
  }

  else if (m_stabilizationType == 2)
  {

    CORRESPONDENCE::computeCorrespondenceStabilityForce(volume,
                                                        horizon,
                                                        modelCoordinates,
                                                        coordinatesNP1,
                                                        deformationGradient,
                                                        forceDensity,
                                                        neighborhoodList,
                                                        numOwnedPoints,
                                                        m_bulkModulus,
                                                        m_hourglassCoefficient,
                                                        bondDamageNP1);
  }
}

void PeridigmNS::CorrespondenceMaterial::computeJacobian(const double dt,
                                                         const int numOwnedPoints,
                                                         const int *ownedIDs,
                                                         const int *neighborhoodList,
                                                         PeridigmNS::DataManager &dataManager,
                                                         PeridigmNS::SerialMatrix &jacobian,
                                                         PeridigmNS::Material::JacobianType jacobianType) const
{

  if (m_applyAutomaticDifferentiationJacobian)
  {
    // Compute the Jacobian via automatic differentiation
    computeAutomaticDifferentiationJacobian(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian, jacobianType);
  }
  else
  {
    //  // Call the base class function, which computes the Jacobian by finite difference
    computeJacobianFiniteDifference(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian, CENTRAL_DIFFERENCE, jacobianType);
  }
}

void PeridigmNS::CorrespondenceMaterial::computeAutomaticDifferentiationJacobian(const double dt,
                                                                                 const int numOwnedPoints,
                                                                                 const int *ownedIDs,
                                                                                 const int *neighborhoodList,
                                                                                 PeridigmNS::DataManager &dataManager,
                                                                                 PeridigmNS::SerialMatrix &jacobian,
                                                                                 PeridigmNS::Material::JacobianType jacobianType) const
{
  // Compute contributions to the tangent matrix on an element-by-element basis

  // To reduce memory re-allocation, use static variable to store Fad types for
  // current coordinates (independent variables).

  static vector<Sacado::Fad::DFad<double>> coordinatesNP1_AD;
  static vector<Sacado::Fad::DFad<double>> coordinates_AD;

  static Sacado::Fad::DFad<double> C_AD[6][6];

  // Loop over all points.
  int neighborhoodListIndex = 0;
  for (int iID = 0; iID < numOwnedPoints; ++iID)
  {

    // Create a temporary neighborhood consisting of a single point and its neighbors.
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    int numEntries = numNeighbors + 1;
    int numDof = 3 * numEntries;
    vector<int> tempMyGlobalIDs(numEntries);
    // Put the node at the center of the neighborhood at the beginning of the list.
    tempMyGlobalIDs[0] = dataManager.getOwnedScalarPointMap()->GID(iID);
    vector<int> tempNeighborhoodList(numEntries);
    tempNeighborhoodList[0] = numNeighbors;

    for (int iNID = 0; iNID < numNeighbors; ++iNID)
    {
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      tempMyGlobalIDs[iNID + 1] = dataManager.getOverlapScalarPointMap()->GID(neighborID);
      tempNeighborhoodList[iNID + 1] = iNID + 1;
    }

    Epetra_SerialComm serialComm;
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 3, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));

    // Create a temporary DataManager containing data for this point and its neighborhood.
    PeridigmNS::DataManager tempDataManager;
    tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
                            tempOneDimensionalMap,
                            Teuchos::RCP<const Epetra_BlockMap>(),
                            tempThreeDimensionalMap,
                            tempBondMap);

    // The temporary data manager will have the same field specs and data as the real data manager.
    vector<int> fieldIds = dataManager.getFieldIds();
    tempDataManager.allocateData(fieldIds);
    tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);

    // Set up numOwnedPoints and ownedIDs.
    // There is only one owned ID, and it has local ID zero in the tempDataManager.
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(tempNumOwnedPoints);
    tempOwnedIDs[0] = 0;

    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
    // Resize scratchMatrix if necessary
    if (scratchMatrix.Dimension() < numDof)
      scratchMatrix.Resize(numDof);

    // Create a list of global indices for the rows/columns in the scratch matrix.
    vector<int> globalIndices(numDof);
    for (int i = 0; i < numEntries; ++i)
    {
      int globalID = tempOneDimensionalMap->GID(i);
      for (int j = 0; j < 3; ++j)
        globalIndices[3 * i + j] = 3 * globalID + j;
    }

    // Extract pointers to the underlying data in the constitutive data array.

    double *modelCoordinates, *volume, *coordinatesNP1, *angles, *detachedNodes;
    tempDataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&angles);
    tempDataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    tempDataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
    tempDataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinatesNP1);
    dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
    double *horizon, *bondDamageNP1;
    tempDataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
    tempDataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);

    // Create arrays of Fad objects for the current coordinates, dilatation, and force density
    // Modify the existing vector of Fad objects for the current coordinates
    if ((int)coordinatesNP1_AD.size() < numDof)
    {
      coordinates_AD.resize(numDof);
      coordinatesNP1_AD.resize(numDof);
    }
    for (int i = 0; i < numDof; ++i)
    {
      coordinates_AD[i].diff(i, numDof);
      coordinatesNP1_AD[i].diff(i, numDof);
      coordinatesNP1_AD[i].val() = coordinatesNP1[i];
    }

    vector<Sacado::Fad::DFad<double>> shapeTensorInverse_AD;
    vector<Sacado::Fad::DFad<double>> deformationGradient_AD;
    vector<Sacado::Fad::DFad<double>> cauchyStress_AD;
    vector<Sacado::Fad::DFad<double>> cauchyStressNP1_AD;
    vector<Sacado::Fad::DFad<double>> cauchyStressPlastic_AD;
    vector<Sacado::Fad::DFad<double>> partialStress_AD;
    vector<Sacado::Fad::DFad<double>> tempStressX_AD;
    vector<Sacado::Fad::DFad<double>> tempStressY_AD;
    vector<Sacado::Fad::DFad<double>> tempStressZ_AD;
    vector<Sacado::Fad::DFad<double>> hourglassStiff_AD;

    partialStress_AD.resize(numDof * numDof);
    shapeTensorInverse_AD.resize(numDof * numDof);
    deformationGradient_AD.resize(numDof * numDof);
    cauchyStress_AD.resize(numDof * numDof);
    cauchyStressNP1_AD.resize(numDof * numDof);
    cauchyStressPlastic_AD.resize(numDof * numDof);
    tempStressX_AD.resize(numDof * numDof);
    tempStressY_AD.resize(numDof * numDof);
    tempStressZ_AD.resize(numDof * numDof);
    hourglassStiff_AD.resize(numDof * numDof);

    for (int i = 0; i < 6; ++i)
    {
      for (int j = 0; j < 6; ++j)
      {
        C_AD[i][j].val() = C[i][j];
      }
    }

    int shapeTensorReturnCode = 0;
    shapeTensorReturnCode =
        CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                                   horizon,
                                                                                   modelCoordinates,
                                                                                   &coordinatesNP1_AD[0],
                                                                                   &shapeTensorInverse_AD[0],
                                                                                   &deformationGradient_AD[0],
                                                                                   bondDamageNP1,
                                                                                   &tempNeighborhoodList[0],
                                                                                   tempNumOwnedPoints,
                                                                                   m_plane,
                                                                                   detachedNodes);

    string shapeTensorErrorMessage =
        "**** Error:  MATRICES::Invert2by2Matrix() failed to invert.\n";

    TestForTermination(shapeTensorReturnCode == 1, shapeTensorErrorMessage);

    double *cauchyStress, *cauchyStressNP1, *cauchyStressPlastic;
    tempDataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStress);
    tempDataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNP1);
    dataManager.getData(m_unrotatedCauchyStressPlasticFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressPlastic);

    for (int i = 0; i < numDof; ++i)
    {
      cauchyStress_AD[i].diff(i, numDof);
      cauchyStress_AD[i].val() = cauchyStress[i];
      cauchyStressNP1_AD[i].diff(i, numDof);
      cauchyStressNP1_AD[i].val() = cauchyStressNP1[i];
      cauchyStressPlastic_AD[i].diff(i, numDof);
      cauchyStressPlastic_AD[i].val() = cauchyStressPlastic[i];
    }
    // double *temperature;
    // tbd
    // CORRESPONDENCE::updateElasticCauchyStressSmallDef(&deformationGradient_AD[0],
    //                                      &cauchyStress_AD[0],
    //                                      &cauchyStressNP1_AD[0],
    //                                      tempNumOwnedPoints,
    //                                      &C_AD[0],
    //                                      angles,
    //                                      m_type,
    //                                      dt,
    //                                      m_alpha,
    //                                      temperature,
    //                                      m_applyThermalStrains,
    //                                      m_hencky);

    // define the sacadoFAD force vector
    vector<Sacado::Fad::DFad<double>> force_AD(numDof);
    // std::cout<<m_stabilizationType<<std::endl;
    CORRESPONDENCE::computeForcesAndStresses(
        tempNumOwnedPoints,
        &tempNeighborhoodList[0],
        volume,
        horizon,
        modelCoordinates,
        &coordinatesNP1_AD[0],
        &deformationGradient_AD[0],
        &cauchyStressNP1_AD[0],
        &cauchyStressPlastic_AD[0],
        &shapeTensorInverse_AD[0],
        bondDamageNP1,
        &C_AD[0],
        angles,
        &force_AD[0],
        &partialStress_AD[0],
        &tempStressX_AD[0],
        &tempStressY_AD[0],
        &tempStressZ_AD[0],
        &hourglassStiff_AD[0],
        m_hourglassCoefficient,
        m_stabilizationType,
        m_plane,
        m_plast,
        m_adaptHourGlass,
        detachedNodes);

    // Load derivative values into scratch matrix
    // Multiply by volume along the way to convert force density to force
    double value;
    for (int row = 0; row < numDof; ++row)
    {
      for (int col = 0; col < numDof; ++col)
      {
        value = force_AD[row].dx(col); //--> I think this must be it, because forces are already provided
                                       // value = force_AD[row].dx(col) * volume[row/3]; // given by peridigm org
        TestForTermination(!std::isfinite(value), "**** NaN detected in correspondence::computeAutomaticDifferentiationJacobian().\n");
        scratchMatrix(row, col) = value;
      }
    }

    // Sum the values into the global tangent matrix (this is expensive).
    if (jacobianType == PeridigmNS::Material::FULL_MATRIX)
      jacobian.addValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    else if (jacobianType == PeridigmNS::Material::BLOCK_DIAGONAL)
    {
      jacobian.addBlockDiagonalValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    }
    else // unknown jacobian type
      TestForTermination(true, "**** Unknown Jacobian Type\n");
  }
}

void PeridigmNS::CorrespondenceMaterial::computeJacobianFiniteDifference(const double dt,
                                                                         const int numOwnedPoints,
                                                                         const int *ownedIDs,
                                                                         const int *neighborhoodList,
                                                                         PeridigmNS::DataManager &dataManager,
                                                                         PeridigmNS::SerialMatrix &jacobian,
                                                                         FiniteDifferenceScheme finiteDifferenceScheme,
                                                                         PeridigmNS::Material::JacobianType jacobianType) const
{

  // The mechanics Jacobian is of the form:
  //
  // dF_0x/dx_0  dF_0x/dy_0  dF_0x/dz_0  dF_0x/dx_1  dF_0x/dy_1  dF_0x/dz_1  ...  dF_0x/dx_n  dF_0x/dy_n  dF_0x/dz_n
  // dF_0y/dx_0  dF_0y/dy_0  dF_0y/dz_0  dF_0y/dx_1  dF_0y/dy_1  dF_0y/dz_1  ...  dF_0y/dx_n  dF_0y/dy_n  dF_0y/dz_n
  // dF_0z/dx_0  dF_0z/dy_0  dF_0z/dz_0  dF_0z/dx_1  dF_0z/dy_1  dF_0z/dz_1  ...  dF_0z/dx_n  dF_0z/dy_n  dF_0z/dz_n
  // dF_1x/dx_0  dF_1x/dy_0  dF_1x/dz_0  dF_1x/dx_1  dF_1x/dy_1  dF_1x/dz_1  ...  dF_1x/dx_n  dF_1x/dy_n  dF_1x/dz_n
  // dF_1y/dx_0  dF_1y/dy_0  dF_1y/dz_0  dF_1y/dx_1  dF_1y/dy_1  dF_1y/dz_1  ...  dF_1y/dx_n  dF_1y/dy_n  dF_1y/dz_n
  // dF_1z/dx_0  dF_1z/dy_0  dF_1z/dz_0  dF_1z/dx_1  dF_1z/dy_1  dF_1z/dz_1  ...  dF_1z/dx_n  dF_1z/dy_n  dF_1z/dz_n
  //     .           .           .           .           .           .                .           .           .
  //     .           .           .           .           .           .                .           .           .
  //     .           .           .           .           .           .                .           .           .
  // dF_nx/dx_0  dF_nx/dy_0  dF_nx/dz_0  dF_nx/dx_1  dF_nx/dy_1  dF_nx/dz_1  ...  dF_nx/dx_n  dF_nx/dy_n  dF_nx/dz_n
  // dF_ny/dx_0  dF_ny/dy_0  dF_ny/dz_0  dF_ny/dx_1  dF_ny/dy_1  dF_ny/dz_1  ...  dF_ny/dx_n  dF_ny/dy_n  dF_ny/dz_n
  // dF_nz/dx_0  dF_nz/dy_0  dF_nz/dz_0  dF_nz/dx_1  dF_nz/dy_1  dF_nz/dz_1  ...  dF_nz/dx_n  dF_nz/dy_n  dF_nz/dz_n

  // Each entry is computed by finite difference:
  //
  // Forward difference:
  // dF_0x/dx_0 = ( F_0x(perturbed x_0) - F_0x(unperturbed) ) / epsilon
  //
  // Central difference:
  // dF_0x/dx_0 = ( F_0x(positive perturbed x_0) - F_0x(negative perturbed x_0) ) / ( 2.0*epsilon )

  TestForTermination(m_finiteDifferenceProbeLength == DBL_MAX, "**** Finite-difference Jacobian requires that the \"Finite Difference Probe Length\" parameter be set.\n");
  double epsilon = m_finiteDifferenceProbeLength;

  PeridigmNS::DegreesOfFreedomManager &dofManager = PeridigmNS::DegreesOfFreedomManager::self();
  bool solveForDisplacement = dofManager.displacementTreatedAsUnknown();
  // bool solveForTemperature = dofManager.temperatureTreatedAsUnknown();
  int numDof = dofManager.totalNumberOfDegreesOfFreedom();
  int numDisplacementDof = dofManager.numberOfDisplacementDegreesOfFreedom();
  // int numTemperatureDof = dofManager.numberOfTemperatureDegreesOfFreedom();
  int displacementDofOffset = dofManager.displacementDofOffset();
  // int temperatureDofOffset = dofManager.temperatureDofOffset();

  // Get field ids for all relevant data
  PeridigmNS::FieldManager &fieldManager = PeridigmNS::FieldManager::self();
  int volumeFId(-1), coordinatesFId(-1), velocityFId(-1), forceDensityFId(-1); //, temperatureFId(-1), fluxDivergenceFId(-1);
  volumeFId = fieldManager.getFieldId("Volume");
  if (solveForDisplacement)
  {
    coordinatesFId = fieldManager.getFieldId("Coordinates");
    velocityFId = fieldManager.getFieldId("Velocity");
    forceDensityFId = fieldManager.getFieldId("Force_Density");
  }


  int neighborhoodListIndex = 0;
  for (int iID = 0; iID < numOwnedPoints; ++iID)
  {

    // Create a temporary neighborhood consisting of a single point and its neighbors.
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    vector<int> tempMyGlobalIDs(numNeighbors + 1);
    // Put the node at the center of the neighborhood at the beginning of the list.
    tempMyGlobalIDs[0] = dataManager.getOwnedScalarPointMap()->GID(iID);
    vector<int> tempNeighborhoodList(numNeighbors + 1);
    tempNeighborhoodList[0] = numNeighbors;
    for (int iNID = 0; iNID < numNeighbors; ++iNID)
    {
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      tempMyGlobalIDs[iNID + 1] = dataManager.getOverlapScalarPointMap()->GID(neighborID);
      tempNeighborhoodList[iNID + 1] = iNID + 1;
    }

    Epetra_SerialComm serialComm;
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors + 1, numNeighbors + 1, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors + 1, numNeighbors + 1, &tempMyGlobalIDs[0], 3, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));

    // Create a temporary DataManager containing data for this point and its neighborhood.
    PeridigmNS::DataManager tempDataManager;
    tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
                            tempOneDimensionalMap,
                            Teuchos::RCP<const Epetra_BlockMap>(),
                            tempThreeDimensionalMap,
                            tempBondMap);

    // The temporary data manager will have the same fields and data as the real data manager.
    vector<int> fieldIds = dataManager.getFieldIds();
    tempDataManager.allocateData(fieldIds);
    tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);

    // Set up numOwnedPoints and ownedIDs.
    // There is only one owned ID, and it has local ID zero in the tempDataManager.
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(1);
    tempOwnedIDs[0] = 0;

    // Extract pointers to the underlying data.
    double *volume, *y, *v, *force; //, *temperature, *fluxDivergence;
    tempDataManager.getData(volumeFId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    if (solveForDisplacement)
    {
      tempDataManager.getData(coordinatesFId, PeridigmField::STEP_NP1)->ExtractView(&y);
      tempDataManager.getData(velocityFId, PeridigmField::STEP_NP1)->ExtractView(&v);
      tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1)->ExtractView(&force);
    }
    // if (solveForTemperature) {
    //   tempDataManager.getData(temperatureFId, PeridigmField::STEP_NP1)->ExtractView(&temperature);
    //   tempDataManager.getData(fluxDivergenceFId, PeridigmField::STEP_NP1)->ExtractView(&fluxDivergence);
    // }

    // Create a temporary vector for storing force and/or flux divergence.
    Teuchos::RCP<Epetra_Vector> forceVector, tempForceVector, fluxDivergenceVector, tempFluxDivergenceVector;
    double *tempForce; //, *tempFluxDivergence;
    if (solveForDisplacement)
    {
      forceVector = tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1);
      tempForceVector = Teuchos::rcp(new Epetra_Vector(*forceVector));
      tempForceVector->ExtractView(&tempForce);
    }
    // if (solveForTemperature) {
    //   fluxDivergenceVector = tempDataManager.getData(fluxDivergenceFId, PeridigmField::STEP_NP1);
    //   tempFluxDivergenceVector = Teuchos::rcp(new Epetra_Vector(*fluxDivergenceVector));
    //   tempFluxDivergenceVector->ExtractView(&tempFluxDivergence);
    // }

    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
    // Resize scratchMatrix if necessary
    if (scratchMatrix.Dimension() < numDof * (numNeighbors + 1))
      scratchMatrix.Resize(numDof * (numNeighbors + 1));

    // Create a list of global indices for the rows/columns in the scratch matrix.
    vector<int> globalIndices(numDof * (numNeighbors + 1));
    for (int i = 0; i < numNeighbors + 1; ++i)
    {
      int globalID = tempOneDimensionalMap->GID(i);
      for (int j = 0; j < numDof; ++j)
      {
        globalIndices[numDof * i + j] = numDof * globalID + j;
      }
    }

    if (finiteDifferenceScheme == FORWARD_DIFFERENCE)
    {
      if (solveForDisplacement)
      {
        // Compute and store the unperturbed force.
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);

        for (int i = 0; i < forceVector->MyLength(); ++i)
          tempForce[i] = force[i];
      }
      // if (solveForTemperature) {
      //   // Compute and store the unperturbed flux divergence.
      //   computeFluxDivergence(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
      //   for(int i=0 ; i<fluxDivergenceVector->MyLength() ; ++i)
      //     tempFluxDivergence[i] = fluxDivergence[i];
      // }
    }

    // Perturb one dof in the neighborhood at a time and compute the force and/or flux divergence.
    // The point itself plus each of its neighbors must be perturbed.
    for (int iNID = 0; iNID < numNeighbors + 1; ++iNID)
    {

      int perturbID;
      if (iNID < numNeighbors)
        perturbID = tempNeighborhoodList[iNID + 1];
      else
        perturbID = 0;

      // Displacement degrees of freedom
      for (int dof = 0; dof < numDisplacementDof; ++dof)
      {

        // Perturb a dof and compute the forces.
        double oldY = y[numDof * perturbID + dof];
        double oldV = v[numDof * perturbID + dof];

        if (finiteDifferenceScheme == CENTRAL_DIFFERENCE)
        {
          // Compute and store the negatively perturbed force.
          y[numDof * perturbID + dof] -= epsilon;
          v[numDof * perturbID + dof] -= epsilon / dt;
          computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
          y[numDof * perturbID + dof] = oldY;
          v[numDof * perturbID + dof] = oldV;
          for (int i = 0; i < forceVector->MyLength(); ++i)
            tempForce[i] = force[i];
        }

        // Compute the purturbed force.
        y[numDof * perturbID + dof] += epsilon;
        v[numDof * perturbID + dof] += epsilon / dt;
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        y[numDof * perturbID + dof] = oldY;
        v[numDof * perturbID + dof] = oldV;

        for (int i = 0; i < numNeighbors + 1; ++i)
        {
          int forceID;
          if (i < numNeighbors)
            forceID = tempNeighborhoodList[i + 1];
          else
            forceID = 0;

          for (int d = 0; d < numDof; ++d)
          {
            double value = (force[numDof * forceID + d] - tempForce[numDof * forceID + d]) / epsilon;
            if (finiteDifferenceScheme == CENTRAL_DIFFERENCE)
              value *= 0.5;
            scratchMatrix(numDof * forceID + displacementDofOffset + d, numDof * perturbID + displacementDofOffset + dof) = value;
          }
        }
      }

      // Temperature degrees of freedom
      // if(solveForTemperature){
      //
      //   // Perturb a temperature value and compute the flux divergence.
      //   double oldTemperature = temperature[perturbID];
      //
      //   if(finiteDifferenceScheme == CENTRAL_DIFFERENCE){
      //     // Compute and store the negatively perturbed flux divergence.
      //     temperature[perturbID] -= epsilon;
      //     computeFluxDivergence(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
      //     temperature[perturbID] = oldTemperature;
      //     for(int i=0 ; i<fluxDivergenceVector->MyLength() ; ++i)
      //       tempFluxDivergence[i] = fluxDivergence[i];
      //   }
      //
      //// Compute the purturbed flux divergence.
      //  temperature[perturbID] += epsilon;
      //  computeFluxDivergence(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
      //  temperature[perturbID] = oldTemperature;
      //
      //  for(int i=0 ; i<numNeighbors+1 ; ++i){
      //    int fluxDivergenceID;
      //    if(i < numNeighbors)
      //      fluxDivergenceID = tempNeighborhoodList[i+1];
      //    else
      //      fluxDivergenceID = 0;
      //
      //    double value = ( fluxDivergence[fluxDivergenceID] - tempFluxDivergence[fluxDivergenceID] ) / epsilon;
      //    if(finiteDifferenceScheme == CENTRAL_DIFFERENCE)
      //      value *= 0.5;
      //    scratchMatrix(numDof*fluxDivergenceID + temperatureDofOffset, numDof*perturbID + temperatureDofOffset) = value;
      //}
      //}
    }

    // Multiply by nodal volume
    for (unsigned int row = 0; row < globalIndices.size(); ++row)
    {
      for (unsigned int col = 0; col < globalIndices.size(); ++col)
      {
        scratchMatrix(row, col) *= volume[row / numDof];
      }
    }

    // Check for NaNs
    for (unsigned int row = 0; row < globalIndices.size(); ++row)
    {
      for (unsigned int col = 0; col < globalIndices.size(); ++col)
      {
        TestForTermination(!std::isfinite(scratchMatrix(row, col)), "**** NaN detected in finite-difference Jacobian.\n");
      }
    }

    // Sum the values into the global tangent matrix (this is expensive).
    if (jacobianType == PeridigmNS::Material::FULL_MATRIX)
      jacobian.addValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    else if (jacobianType == PeridigmNS::Material::BLOCK_DIAGONAL)
    {
      jacobian.addBlockDiagonalValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    }
    else // unknown jacobian type
      TestForTermination(true, "**** Unknown Jacobian Type\n");
  }
}



