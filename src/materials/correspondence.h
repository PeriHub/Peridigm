//! \file correspondence.h
//@HEADER
// ************************************************************************
//
//                             Peridigm
//               copyright (2011) Sandia Corporation
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
// Stewart A. Silling  Sasilli@sandia.gov
//
// ************************************************************************
//@HEADER
#ifndef CORRESPONDENCE_H
#define CORRESPONDENCE_H


namespace CORRESPONDENCE {

template<typename ScalarT>
int EigenVec2D
(
 const ScalarT* a,
 ScalarT* result
);

//! Hencky-Strain E = 0.5*ln(C).
template<typename ScalarT>
int computeLogStrain
(
 const ScalarT* defGrad,
 ScalarT* strain
);

template<typename ScalarT>
void rotateCauchyStress
(
const ScalarT* rotationTensor,
const ScalarT* unrotatedCauchyStress,
ScalarT* rotatedCauchyStress,
int numPoints
);

template<typename ScalarT>
int computeGradientWeights
(
    const double* horizon,
    const ScalarT* coordinates,
    const double* volume,
    const ScalarT* jacobianDeterminant,
    ScalarT* gradientWeight1,
    ScalarT* gradientWeight2,
    ScalarT* gradientWeight3,
    const int accuracyOrder,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints
);

template<typename ScalarT>
int computeShapeTensorInverseAndApproximateDeformationGradient
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const ScalarT* coordinatesNP1,
ScalarT* shapeTensorInverse,
ScalarT* deformationGradient,
const double* bondDamageNP1,
const int* neighborhoodList,
int numPoints,
const bool type,
double* detachedNodes
);

// Calculation of stretch rates following Flanagan & Taylor
template<typename ScalarT>
int computeUnrotatedRateOfDeformationAndRotationTensor(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const ScalarT* velocities,
ScalarT* deformationGradient,
const ScalarT* shapeTensorInverse,
ScalarT* leftStretchTensorN,
const ScalarT* rotationTensorN,
ScalarT* leftStretchTensorNP1,
ScalarT* rotationTensorNP1,
ScalarT* unrotatedRateOfDeformation,
const int* neighborhoodList,
int numPoints,
double dt,
const double* bondDamage,
const bool type,
double* detachedNodes
);

template<typename ScalarT>
void getStrain
(
const int numPoints,
const ScalarT* defGrad,
const double alpha[][3],
const ScalarT* temperature,
const bool hencky,
const bool applyThermalStrains,
ScalarT* strain
);

template<typename ScalarT>
void getVonMisesStress
(
const int numPoints,
const ScalarT* sigmaNP1,
ScalarT* vmStress
);

template<typename ScalarT>
void addTemperatureStrain
(
const double alpha[][3],
const ScalarT temperature,
ScalarT* strain
);

template<typename ScalarT>
void computeGreenLagrangeStrain
(
const ScalarT* defGrad,
ScalarT* strain
);

double FLAWFUNCTION(
const bool isFlaw,
const double yieldStress,
const double flawMagnitude, 
const double flawSize, 
const double* modelCoord, 
const double flawLocationX, 
const double flawLocationY, 
const double flawLocationZ, 
const int type
);
template<typename ScalarT>
void computeHourglassForce
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const ScalarT* coordinates,
const ScalarT* deformationGradient,
ScalarT* hourglassForceDensity,
const int* neighborhoodList,
int numPoints,
double bulkModulus,
double hourglassCoefficient,
const double* bondDamage
);

template<typename ScalarT>
void computeCorrespondenceStabilityForce
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const ScalarT* coordinates,
const ScalarT* deformationGradient,
ScalarT* hourglassForceDensity,
const int* neighborhoodList,
int numPoints,
double bulkModulus,
double hourglassCoefficient,
const double* bondDamage
);
template<typename ScalarT>
void computeCorrespondenceStabilityWanEtAl
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const ScalarT* coordinates,
const int* neighborhoodList,
int numPoints,
const ScalarT* deformationGradient,
const ScalarT* shapeTensorInverse,
const double* C,
const double* bondDamage,
ScalarT* hourglassForceDensity,
double m_hourglassCoefficient
);
template<typename ScalarT>
void createHourglassStiffness
(
const ScalarT C[][6],
const double alpha[],
const ScalarT* shapeTensorInverse,
ScalarT* hourglassStiff
);
template<typename ScalarT>
void computeCorrespondenceStabilityWanEtAlShort
(
const ScalarT FxsiX,
const ScalarT FxsiY,
const ScalarT FxsiZ,
const ScalarT deformedBondX,
const ScalarT deformedBondY,
const ScalarT deformedBondZ,
const ScalarT* hourglassStiff,
ScalarT* TS
);
template<typename ScalarT>
void setValuesForDetachedNodes
(
ScalarT* deformationGradient,
ScalarT* leftStretchTensor,
ScalarT* rotationTensor,
ScalarT* unrotatedRateOfDeformation,
ScalarT* shapeTensorInverse,
const double* detachedNodes,
const int numPoints
 );
template<typename ScalarT>
void MatMul
(
int n,
const ScalarT A[][6],
const ScalarT B[][6],
ScalarT C[][6],
bool transpose
);

template<typename ScalarT>
void setOnesOnDiagonalFullTensor(ScalarT* tensor, int numPoints);

void computeHeatFlowState(    
const double* modelCoord,
const int numOwnedPoints,
const int* neighborhoodList,
const double* shapeTensorInverse,
const double* temperature,
const double* horizon,
const double* lambda,
const double* volume,
const double* detachedNodes,
const double* bondDamage,
const double* angles,
const bool twoD,
const bool bondbased,
const bool applyThermalPrintBedFlow,
const double lambdaBed,
const double TBed,
double* heatFlowState
);

void computeHeatTransfer(  
    const double* modelCoordinates,  
    const int numOwnedPoints,
    const int* neighborhoodList,
    const double* volume,
    const double* temperature,
    const double* horizon,
    const double* detachedNodes,
    const double* bondDamage,
    const bool twoD,
    const double alpha,
    const double Tenv,
    const double factor,
    const double surfaceCorrection,
    const double limit,
    const bool applyThermalPrintBedFlow,
    double* specificVolume,
    double* heatFlowState
);

template<typename ScalarT>
void createRotatedStiff
(
const ScalarT C[][6],
const ScalarT* rotMat,
ScalarT Cnew[][6]
);

void getOrientations
  (
const int numOwnedPoints,
const double* angles,
double* orientations
);

template<typename ScalarT>
void computeForcesAndStresses
  (
const int numOwnedPoints,
const int* neighborhoodList,
const double* volume,
const double* horizon,
const double* modelCoordinates,
const ScalarT* coordinatesNP1,
const ScalarT* deformationGradient,
const ScalarT* cauchyStressNP1,
const ScalarT* cauchyStressPlastic,
const ScalarT* shapeTensorInverse,
const double* bondDamage,
const ScalarT C[][6],
const double* angles,
ScalarT* force,
ScalarT* partialStress,
ScalarT* tempStressX,
ScalarT* tempStressY,
ScalarT* tempStressZ,
ScalarT* hourglassStiff,
const double m_hourglassCoefficient,
const int m_stabilizationType,
const bool m_plane,
const bool m_plast,
const bool m_adaptHourGlass,
double* detachedNodes
    );
template<typename ScalarT>
void setOnesOnDiagonalFullTensor(ScalarT* tensor, int numPoints);

template<typename ScalarT>
void computeUndamagedWeightedVolume
(
const double* volume,
double* weightedVolume,
const ScalarT* jacobianDeterminant,
const double* horizon,
const ScalarT* coordinates,
const int* neighborhoodList,
int numPoints
);

template<typename ScalarT>
void updateGradientWeightEvaluationFlag
(
    const ScalarT* damageN,
    const ScalarT* damageNP1,
    ScalarT* gradientWeightEvaluationFlag,
    int numPoints
);

template<typename ScalarT>
int computeLagrangianGradientWeights
(
    const double* horizon,
    const double* modelCoordinates,
    const double* volume,
    ScalarT* gradientWeightX,
    ScalarT* gradientWeightY,
    ScalarT* gradientWeightZ,
    ScalarT* gradientWeightEvaluationFlag,
    const double* bondDamage,
    double* influenceState,
    const int accuracyOrder,
    const int* neighborhoodList,
    int numPoints
);

template<typename ScalarT>
void computeDeformationGradient
(
    const double* volume,
    const ScalarT* displacements,
    const ScalarT* velocities,
    const ScalarT* gradientWeightX,
    const ScalarT* gradientWeightY,
    const ScalarT* gradientWeightZ,
    ScalarT* deformationGradientX,
    ScalarT* deformationGradientY,
    ScalarT* deformationGradientZ,
    ScalarT* deformationGradientDotX,
    ScalarT* deformationGradientDotY,
    ScalarT* deformationGradientDotZ,
    const int* neighborhoodList,
    int numPoints
);

template<typename ScalarT>
void computeWeightedVolume
(
const double* volume,
double* weightedVolume,
const ScalarT* jacobianDeterminant,
const double* horizon,
const ScalarT* coordinates,
const double* flyingPointFlag,
const double* bondDamage,
const int* neighborhoodList,
int numPoints
);

template<typename ScalarT>
int computeShapeTensorInverseAndApproximateNodeLevelVelocityGradient
(
const double* volume,
const ScalarT* jacobianDeterminantN,
ScalarT* jacobianDeterminantNP1,
const double* horizon,
const ScalarT* coordinates,
const ScalarT* velocities,
ScalarT* shapeTensorInverse,
ScalarT* velocityGradient,
const double* flyingPointFlag,
const double* bondDamage,
const int* neighborhoodList,
int numPoints,
double dt
);

template<typename ScalarT>
void computeVelocityGradient
(
    const double* volume,
    const ScalarT* jacobianDeterminantN,
    ScalarT* jacobianDeterminantNP1,
    const ScalarT* velocities,
    const ScalarT* gradientWeight1,
    const ScalarT* gradientWeight2,
    const ScalarT* gradientWeight3,
    ScalarT* velocityGradient,
    ScalarT* velocityGradientX,
    ScalarT* velocityGradientY,
    ScalarT* velocityGradientZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints,
    double dt
);

template<typename ScalarT>
int computeShapeTensorInverseAndApproximateNodeLevelVelocityGradient
(
const double* volume,
const ScalarT* jacobianDeterminantN,
ScalarT* jacobianDeterminantNP1,
const double* horizon,
const ScalarT* coordinates,
const ScalarT* velocities,
ScalarT* shapeTensorInverse,
ScalarT* velocityGradient,
ScalarT* velocityGradientX,
ScalarT* velocityGradientY,
ScalarT* velocityGradientZ,
const double* flyingPointFlag,
const double* bondDamage,
const int* neighborhoodList,
int numPoints,
double dt
);


template<typename ScalarT>
void updateDeformationGradient
(
const ScalarT* velocityGradient,
const ScalarT* deformationGradientN,
ScalarT* deformationGradientNP1,
const double* flyingPointFlag,
int numPoints,
double dt
);

template<typename ScalarT>
void computeGreenLagrangeStrain
(
const ScalarT* deformationGradient,
ScalarT* greenLagrangeStrain,
const double* flyingPointFlag,
int numPoints
);

template<typename ScalarT>
void computeWeightedVolume
(
    const double* volume,
    ScalarT* weightedVolume,
    const double* influenceState,
    const int* neighborhoodList,
    int numPoints
);

template<typename ScalarT>
int computeNodeLevelUnrotatedRateOfDeformationAndRotationTensor(
const ScalarT* velocityGradient,
const ScalarT* leftStretchTensorN,
const ScalarT* rotationTensorN,
ScalarT* leftStretchTensorNP1,
ScalarT* rotationTensorNP1,
ScalarT* unrotatedRateOfDeformation,
const double* flyingPointFlag,
int numPoints,
double dt
);

template<typename ScalarT>
int computeNodeLevelUnrotatedRateOfDeformationAndRotationTensor
(
    const ScalarT* deformationGradientX,
    const ScalarT* deformationGradientY,
    const ScalarT* deformationGradientZ,
    const ScalarT* deformationGradientDotX,
    const ScalarT* deformationGradientDotY,
    const ScalarT* deformationGradientDotZ,
    const ScalarT* leftStretchTensorN,
    const ScalarT* rotationTensorN,
    ScalarT* leftStretchTensorNP1,
    ScalarT* rotationTensorNP1,
    ScalarT* unrotatedRateOfDeformation,
    int numPoints,
    double dt
);


template<typename ScalarT>
void rotateCauchyStress
(
const ScalarT* rotationTensor,
const ScalarT* unrotatedCauchyStress,
ScalarT* rotatedCauchyStress,
const double* flyingPointFlag,
int numPoints
);

template<typename ScalarT>
void updateGreenLagrangeStrain
(
    const ScalarT* deformationGradientX,
    const ScalarT* deformationGradientY,
    const ScalarT* deformationGradientZ,
    const ScalarT* deformationGradientDotX,
    const ScalarT* deformationGradientDotY,
    const ScalarT* deformationGradientDotZ,
    const ScalarT* greenLagrangeStrainN,
    ScalarT* greenLagrangeStrainNP1,
    int numPoints,
    double dt
);


template<typename ScalarT>
void getLinearUnrotatedRateOfDeformation
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const ScalarT* velocities,
ScalarT* deformationGradient,
const ScalarT* shapeTensorInverse,
ScalarT* unrotatedRateOfDeformation,
const int* neighborhoodList,
int numPoints,
const double* bondDamage,
const bool type,
double* detachedNodes
);


}

#endif // CORRESPONDENCE_H
