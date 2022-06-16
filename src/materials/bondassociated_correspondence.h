//! \file bondassociated_correspondence.h
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
void computeBondLevelVelocityGradient
(
const ScalarT* coordinates,
const ScalarT* velocities,
const ScalarT* velocityGradient,
ScalarT* bondLevelVelocityGradientXX,
ScalarT* bondLevelVelocityGradientXY,
ScalarT* bondLevelVelocityGradientXZ,
ScalarT* bondLevelVelocityGradientYX,
ScalarT* bondLevelVelocityGradientYY,
ScalarT* bondLevelVelocityGradientYZ,
ScalarT* bondLevelVelocityGradientZX,
ScalarT* bondLevelVelocityGradientZY,
ScalarT* bondLevelVelocityGradientZZ,
const double* flyingPointFlag,
const int* neighborhoodList,
int numPoints
);

template<typename ScalarT>
void computeBondLevelVelocityGradient
(
const ScalarT* coordinates,
const ScalarT* velocities,
const ScalarT* velocityGradientX,
const ScalarT* velocityGradientY,
const ScalarT* velocityGradientZ,
ScalarT* bondLevelVelocityGradientXX,
ScalarT* bondLevelVelocityGradientXY,
ScalarT* bondLevelVelocityGradientXZ,
ScalarT* bondLevelVelocityGradientYX,
ScalarT* bondLevelVelocityGradientYY,
ScalarT* bondLevelVelocityGradientYZ,
ScalarT* bondLevelVelocityGradientZX,
ScalarT* bondLevelVelocityGradientZY,
ScalarT* bondLevelVelocityGradientZZ,
const double* flyingPointFlag,
const int* neighborhoodList,
int numPoints
);


template<typename ScalarT>
int computeBondLevelUnrotatedRateOfDeformationAndRotationTensor(
const ScalarT* bondLevelVelocityGradientXX, 
const ScalarT* bondLevelVelocityGradientXY, 
const ScalarT* bondLevelVelocityGradientXZ,
const ScalarT* bondLevelVelocityGradientYX, 
const ScalarT* bondLevelVelocityGradientYY, 
const ScalarT* bondLevelVelocityGradientYZ, 
const ScalarT* bondLevelVelocityGradientZX,
const ScalarT* bondLevelVelocityGradientZY,
const ScalarT* bondLevelVelocityGradientZZ,
const ScalarT* bondLevelLeftStretchTensorXXN,
const ScalarT* bondLevelLeftStretchTensorXYN,
const ScalarT* bondLevelLeftStretchTensorXZN,
const ScalarT* bondLevelLeftStretchTensorYXN,
const ScalarT* bondLevelLeftStretchTensorYYN,
const ScalarT* bondLevelLeftStretchTensorYZN,
const ScalarT* bondLevelLeftStretchTensorZXN,
const ScalarT* bondLevelLeftStretchTensorZYN,
const ScalarT* bondLevelLeftStretchTensorZZN,
const ScalarT* bondLevelRotationTensorXXN, 
const ScalarT* bondLevelRotationTensorXYN, 
const ScalarT* bondLevelRotationTensorXZN, 
const ScalarT* bondLevelRotationTensorYXN, 
const ScalarT* bondLevelRotationTensorYYN, 
const ScalarT* bondLevelRotationTensorYZN, 
const ScalarT* bondLevelRotationTensorZXN, 
const ScalarT* bondLevelRotationTensorZYN, 
const ScalarT* bondLevelRotationTensorZZN, 
ScalarT* bondLevelLeftStretchTensorXXNP1,
ScalarT* bondLevelLeftStretchTensorXYNP1,
ScalarT* bondLevelLeftStretchTensorXZNP1,
ScalarT* bondLevelLeftStretchTensorYXNP1,
ScalarT* bondLevelLeftStretchTensorYYNP1,
ScalarT* bondLevelLeftStretchTensorYZNP1,
ScalarT* bondLevelLeftStretchTensorZXNP1,
ScalarT* bondLevelLeftStretchTensorZYNP1,
ScalarT* bondLevelLeftStretchTensorZZNP1,
ScalarT* bondLevelRotationTensorXXNP1,
ScalarT* bondLevelRotationTensorXYNP1,
ScalarT* bondLevelRotationTensorXZNP1,
ScalarT* bondLevelRotationTensorYXNP1,
ScalarT* bondLevelRotationTensorYYNP1,
ScalarT* bondLevelRotationTensorYZNP1,
ScalarT* bondLevelRotationTensorZXNP1,
ScalarT* bondLevelRotationTensorZYNP1,
ScalarT* bondLevelRotationTensorZZNP1,
ScalarT* bondLevelUnrotatedRateOfDeformationXX,
ScalarT* bondLevelUnrotatedRateOfDeformationXY,
ScalarT* bondLevelUnrotatedRateOfDeformationXZ,
ScalarT* bondLevelUnrotatedRateOfDeformationYX,
ScalarT* bondLevelUnrotatedRateOfDeformationYY,
ScalarT* bondLevelUnrotatedRateOfDeformationYZ,
ScalarT* bondLevelUnrotatedRateOfDeformationZX,
ScalarT* bondLevelUnrotatedRateOfDeformationZY,
ScalarT* bondLevelUnrotatedRateOfDeformationZZ,
const double* flyingPointFlag,
const int* neighborhoodList,
int numPoints,
double dt
);

template<typename ScalarT>
int computeBondLevelUnrotatedRateOfDeformationAndRotationTensor
(
    const double* modelCoordinates,
    const ScalarT* coordinates,
    const ScalarT* velocities,
    const ScalarT* deformationGradientX,
    const ScalarT* deformationGradientY,
    const ScalarT* deformationGradientZ,
    const ScalarT* deformationGradientDotX,
    const ScalarT* deformationGradientDotY,
    const ScalarT* deformationGradientDotZ,
    ScalarT* bondLevelDeformationGradientInvXX, 
    ScalarT* bondLevelDeformationGradientInvXY, 
    ScalarT* bondLevelDeformationGradientInvXZ,
    ScalarT* bondLevelDeformationGradientInvYX, 
    ScalarT* bondLevelDeformationGradientInvYY, 
    ScalarT* bondLevelDeformationGradientInvYZ, 
    ScalarT* bondLevelDeformationGradientInvZX,
    ScalarT* bondLevelDeformationGradientInvZY,
    ScalarT* bondLevelDeformationGradientInvZZ,
    ScalarT* bondLevelJacobianDeterminant,
    const ScalarT* bondLevelLeftStretchTensorXXN,
    const ScalarT* bondLevelLeftStretchTensorXYN,
    const ScalarT* bondLevelLeftStretchTensorXZN,
    const ScalarT* bondLevelLeftStretchTensorYXN,
    const ScalarT* bondLevelLeftStretchTensorYYN,
    const ScalarT* bondLevelLeftStretchTensorYZN,
    const ScalarT* bondLevelLeftStretchTensorZXN,
    const ScalarT* bondLevelLeftStretchTensorZYN,
    const ScalarT* bondLevelLeftStretchTensorZZN,
    const ScalarT* bondLevelRotationTensorXXN, 
    const ScalarT* bondLevelRotationTensorXYN, 
    const ScalarT* bondLevelRotationTensorXZN, 
    const ScalarT* bondLevelRotationTensorYXN, 
    const ScalarT* bondLevelRotationTensorYYN, 
    const ScalarT* bondLevelRotationTensorYZN, 
    const ScalarT* bondLevelRotationTensorZXN, 
    const ScalarT* bondLevelRotationTensorZYN, 
    const ScalarT* bondLevelRotationTensorZZN, 
    ScalarT* bondLevelLeftStretchTensorXXNP1,
    ScalarT* bondLevelLeftStretchTensorXYNP1,
    ScalarT* bondLevelLeftStretchTensorXZNP1,
    ScalarT* bondLevelLeftStretchTensorYXNP1,
    ScalarT* bondLevelLeftStretchTensorYYNP1,
    ScalarT* bondLevelLeftStretchTensorYZNP1,
    ScalarT* bondLevelLeftStretchTensorZXNP1,
    ScalarT* bondLevelLeftStretchTensorZYNP1,
    ScalarT* bondLevelLeftStretchTensorZZNP1,
    ScalarT* bondLevelRotationTensorXXNP1,
    ScalarT* bondLevelRotationTensorXYNP1,
    ScalarT* bondLevelRotationTensorXZNP1,
    ScalarT* bondLevelRotationTensorYXNP1,
    ScalarT* bondLevelRotationTensorYYNP1,
    ScalarT* bondLevelRotationTensorYZNP1,
    ScalarT* bondLevelRotationTensorZXNP1,
    ScalarT* bondLevelRotationTensorZYNP1,
    ScalarT* bondLevelRotationTensorZZNP1,
    ScalarT* bondLevelUnrotatedRateOfDeformationXX,
    ScalarT* bondLevelUnrotatedRateOfDeformationXY,
    ScalarT* bondLevelUnrotatedRateOfDeformationXZ,
    ScalarT* bondLevelUnrotatedRateOfDeformationYX,
    ScalarT* bondLevelUnrotatedRateOfDeformationYY,
    ScalarT* bondLevelUnrotatedRateOfDeformationYZ,
    ScalarT* bondLevelUnrotatedRateOfDeformationZX,
    ScalarT* bondLevelUnrotatedRateOfDeformationZY,
    ScalarT* bondLevelUnrotatedRateOfDeformationZZ,
    const double* influenceState,
    const int* neighborhoodList,
    int numPoints,
    double dt
);

template<typename ScalarT>
void rotateBondLevelCauchyStress(
const ScalarT* bondLevelRotationTensorXX,
const ScalarT* bondLevelRotationTensorXY,
const ScalarT* bondLevelRotationTensorXZ,
const ScalarT* bondLevelRotationTensorYX,
const ScalarT* bondLevelRotationTensorYY,
const ScalarT* bondLevelRotationTensorYZ,
const ScalarT* bondLevelRotationTensorZX,
const ScalarT* bondLevelRotationTensorZY,
const ScalarT* bondLevelRotationTensorZZ,
const ScalarT* bondLevelUnrotatedCauchyStressXX,
const ScalarT* bondLevelUnrotatedCauchyStressXY,
const ScalarT* bondLevelUnrotatedCauchyStressXZ,
const ScalarT* bondLevelUnrotatedCauchyStressYX,
const ScalarT* bondLevelUnrotatedCauchyStressYY,
const ScalarT* bondLevelUnrotatedCauchyStressYZ,
const ScalarT* bondLevelUnrotatedCauchyStressZX,
const ScalarT* bondLevelUnrotatedCauchyStressZY,
const ScalarT* bondLevelUnrotatedCauchyStressZZ,
ScalarT* bondLevelRotatedCauchyStressXX,
ScalarT* bondLevelRotatedCauchyStressXY,
ScalarT* bondLevelRotatedCauchyStressXZ,
ScalarT* bondLevelRotatedCauchyStressYX,
ScalarT* bondLevelRotatedCauchyStressYY,
ScalarT* bondLevelRotatedCauchyStressYZ,
ScalarT* bondLevelRotatedCauchyStressZX,
ScalarT* bondLevelRotatedCauchyStressZY,
ScalarT* bondLevelRotatedCauchyStressZZ,
const double* flyingPointFlag,
const int* neighborhoodList,
int numPoints
);

template<typename ScalarT>
void rotateBondLevelCauchyStress
(
    const ScalarT* bondLevelRotationTensorXX,
    const ScalarT* bondLevelRotationTensorXY,
    const ScalarT* bondLevelRotationTensorXZ,
    const ScalarT* bondLevelRotationTensorYX,
    const ScalarT* bondLevelRotationTensorYY,
    const ScalarT* bondLevelRotationTensorYZ,
    const ScalarT* bondLevelRotationTensorZX,
    const ScalarT* bondLevelRotationTensorZY,
    const ScalarT* bondLevelRotationTensorZZ,
    const ScalarT* bondLevelUnrotatedCauchyStressXX,
    const ScalarT* bondLevelUnrotatedCauchyStressXY,
    const ScalarT* bondLevelUnrotatedCauchyStressXZ,
    const ScalarT* bondLevelUnrotatedCauchyStressYX,
    const ScalarT* bondLevelUnrotatedCauchyStressYY,
    const ScalarT* bondLevelUnrotatedCauchyStressYZ,
    const ScalarT* bondLevelUnrotatedCauchyStressZX,
    const ScalarT* bondLevelUnrotatedCauchyStressZY,
    const ScalarT* bondLevelUnrotatedCauchyStressZZ,
    ScalarT* bondLevelRotatedCauchyStressXX,
    ScalarT* bondLevelRotatedCauchyStressXY,
    ScalarT* bondLevelRotatedCauchyStressXZ,
    ScalarT* bondLevelRotatedCauchyStressYX,
    ScalarT* bondLevelRotatedCauchyStressYY,
    ScalarT* bondLevelRotatedCauchyStressYZ,
    ScalarT* bondLevelRotatedCauchyStressZX,
    ScalarT* bondLevelRotatedCauchyStressZY,
    ScalarT* bondLevelRotatedCauchyStressZZ,
    const int* neighborhoodList,
    int numPoints
);


template<typename ScalarT>
void computeNonhomogeneityIntegral
(
const double* volume,
const double* weightedVolume,
const ScalarT* jacobianDeterminant,
const double* horizon,
const ScalarT* coordinates,
const ScalarT* bondLevelCauchyStressXX,
const ScalarT* bondLevelCauchyStressXY,
const ScalarT* bondLevelCauchyStressXZ,
const ScalarT* bondLevelCauchyStressYX,
const ScalarT* bondLevelCauchyStressYY,
const ScalarT* bondLevelCauchyStressYZ,
const ScalarT* bondLevelCauchyStressZX,
const ScalarT* bondLevelCauchyStressZY,
const ScalarT* bondLevelCauchyStressZZ,
ScalarT* nonhomogeneousIntegral,
const double* flyingPointFlag,
const double* bondDamage,
const int* neighborhoodList,
int numPoints
);


template<typename ScalarT>
void computeBondLevelPiolaStress
(
    const ScalarT* bondLevelJacobianDeterminant,
    const ScalarT* bondLevelCauchyStressXX,
    const ScalarT* bondLevelCauchyStressXY,
    const ScalarT* bondLevelCauchyStressXZ,
    const ScalarT* bondLevelCauchyStressYX,
    const ScalarT* bondLevelCauchyStressYY,
    const ScalarT* bondLevelCauchyStressYZ,
    const ScalarT* bondLevelCauchyStressZX,
    const ScalarT* bondLevelCauchyStressZY,
    const ScalarT* bondLevelCauchyStressZZ,
    const ScalarT* bondLevelDeformationGradientInvXX,
    const ScalarT* bondLevelDeformationGradientInvXY,
    const ScalarT* bondLevelDeformationGradientInvXZ,
    const ScalarT* bondLevelDeformationGradientInvYX,
    const ScalarT* bondLevelDeformationGradientInvYY,
    const ScalarT* bondLevelDeformationGradientInvYZ,
    const ScalarT* bondLevelDeformationGradientInvZX,
    const ScalarT* bondLevelDeformationGradientInvZY,
    const ScalarT* bondLevelDeformationGradientInvZZ,
    ScalarT* bondLevelPiolaStressXX,
    ScalarT* bondLevelPiolaStressXY,
    ScalarT* bondLevelPiolaStressXZ,
    ScalarT* bondLevelPiolaStressYX,
    ScalarT* bondLevelPiolaStressYY,
    ScalarT* bondLevelPiolaStressYZ,
    ScalarT* bondLevelPiolaStressZX,
    ScalarT* bondLevelPiolaStressZY,
    ScalarT* bondLevelPiolaStressZZ,
    const int* neighborhoodList,
    int numPoints
);

template<typename ScalarT>
void computeStressIntegral
(
    const double* volume,
    const double* weightedVolume,
    const double* modelCoordinates,
    const ScalarT* bondLevelStressXX,
    const ScalarT* bondLevelStressXY,
    const ScalarT* bondLevelStressXZ,
    const ScalarT* bondLevelStressYX,
    const ScalarT* bondLevelStressYY,
    const ScalarT* bondLevelStressYZ,
    const ScalarT* bondLevelStressZX,
    const ScalarT* bondLevelStressZY,
    const ScalarT* bondLevelStressZZ,
    ScalarT* stressIntegral,
    const double* influenceState,
    const int* neighborhoodList,
    int numPoints
);

}

#endif // CORRESPONDENCE_H
