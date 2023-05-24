//! \file bondassociated_correspondence.cxx

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

#include "bondassociated_correspondence.h"
#include "matrices.h"
#include "elastic_correspondence.h"
#include "material_utilities.h"
#include <Sacado.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <math.h>
#include <functional>
#include <cmath> 
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <vector> 
#include <string> 
namespace CORRESPONDENCE {



//This function computes bond-level velocity gradient
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
)
{
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  const ScalarT* velGrad = velocityGradient;
  const ScalarT* neighborVelGrad;

  ScalarT* bondLevelVelGradXX = bondLevelVelocityGradientXX;
  ScalarT* bondLevelVelGradXY = bondLevelVelocityGradientXY;
  ScalarT* bondLevelVelGradXZ = bondLevelVelocityGradientXZ;
  ScalarT* bondLevelVelGradYX = bondLevelVelocityGradientYX;
  ScalarT* bondLevelVelGradYY = bondLevelVelocityGradientYY;
  ScalarT* bondLevelVelGradYZ = bondLevelVelocityGradientYZ;
  ScalarT* bondLevelVelGradZX = bondLevelVelocityGradientZX;
  ScalarT* bondLevelVelGradZY = bondLevelVelocityGradientZY;
  ScalarT* bondLevelVelGradZZ = bondLevelVelocityGradientZZ;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength, deformedBondLengthSq;
  ScalarT velStateX, velStateY, velStateZ;
  ScalarT temp;

  std::vector<ScalarT> meanVelGradVector(9);
  ScalarT* meanVelGrad = &meanVelGradVector[0];

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, coord+=3, vel+=3, velGrad+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
          bondLevelVelGradXX++, bondLevelVelGradXY++, bondLevelVelGradXZ++, 
          bondLevelVelGradYX++, bondLevelVelGradYY++, bondLevelVelGradYZ++,
          bondLevelVelGradZX++, bondLevelVelGradZY++, bondLevelVelGradZZ++){

        neighborIndex = *neighborListPtr;
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVel = velocities + 3*neighborIndex;
        neighborVelGrad = velocityGradient + 9*neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);

        // The velState is the relative difference in velocities of the nodes at
        // each end of a bond. i.e., v_j - v_i
        velStateX = *(neighborVel)   - *(vel);
        velStateY = *(neighborVel+1) - *(vel+1);
        velStateZ = *(neighborVel+2) - *(vel+2);

        // average of the two points 
        for(int i=0; i<9; i++)
          *(meanVelGrad+i) = 0.5 * (*(velGrad+i) + *(neighborVelGrad+i));

        deformedBondLengthSq = deformedBondLength * deformedBondLength;

        temp = *(meanVelGrad+0) * deformedBondX + *(meanVelGrad+1) * deformedBondY + *(meanVelGrad+2) * deformedBondZ;
        *bondLevelVelGradXX = *(meanVelGrad+0) + (velStateX - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradXY = *(meanVelGrad+1) + (velStateX - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradXZ = *(meanVelGrad+2) + (velStateX - temp) * deformedBondZ/deformedBondLengthSq;

        temp = *(meanVelGrad+3) * deformedBondX + *(meanVelGrad+4) * deformedBondY + *(meanVelGrad+5) * deformedBondZ;
        *bondLevelVelGradYX = *(meanVelGrad+3) + (velStateY - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradYY = *(meanVelGrad+4) + (velStateY - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradYZ = *(meanVelGrad+5) + (velStateY - temp) * deformedBondZ/deformedBondLengthSq;

        temp = *(meanVelGrad+6) * deformedBondX + *(meanVelGrad+7) * deformedBondY + *(meanVelGrad+8) * deformedBondZ;
        *bondLevelVelGradZX = *(meanVelGrad+6) + (velStateZ - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradZY = *(meanVelGrad+7) + (velStateZ - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradZZ = *(meanVelGrad+8) + (velStateZ - temp) * deformedBondZ/deformedBondLengthSq;
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondLevelVelGradXX += numNeighbors; bondLevelVelGradXY += numNeighbors; bondLevelVelGradXZ += numNeighbors; 
      bondLevelVelGradYX += numNeighbors; bondLevelVelGradYY += numNeighbors; bondLevelVelGradYZ += numNeighbors;
      bondLevelVelGradZX += numNeighbors; bondLevelVelGradZY += numNeighbors; bondLevelVelGradZZ += numNeighbors;
    }
  }
}

//This function computes bond-level velocity gradient
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
)
{
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  const ScalarT* velGradX = velocityGradientX;
  const ScalarT* velGradY = velocityGradientY;
  const ScalarT* velGradZ = velocityGradientZ;
  const ScalarT* neighborVelGradX;
  const ScalarT* neighborVelGradY;
  const ScalarT* neighborVelGradZ;

  ScalarT* bondLevelVelGradXX = bondLevelVelocityGradientXX;
  ScalarT* bondLevelVelGradXY = bondLevelVelocityGradientXY;
  ScalarT* bondLevelVelGradXZ = bondLevelVelocityGradientXZ;
  ScalarT* bondLevelVelGradYX = bondLevelVelocityGradientYX;
  ScalarT* bondLevelVelGradYY = bondLevelVelocityGradientYY;
  ScalarT* bondLevelVelGradYZ = bondLevelVelocityGradientYZ;
  ScalarT* bondLevelVelGradZX = bondLevelVelocityGradientZX;
  ScalarT* bondLevelVelGradZY = bondLevelVelocityGradientZY;
  ScalarT* bondLevelVelGradZZ = bondLevelVelocityGradientZZ;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength, deformedBondLengthSq;
  ScalarT velStateX, velStateY, velStateZ;
  ScalarT temp;

  std::vector<ScalarT> meanVelGradVector(9);
  ScalarT* meanVelGrad = &meanVelGradVector[0];

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, coord+=3, vel+=3, 
      velGradX+=3,  velGradY+=3, velGradZ+=3, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
          bondLevelVelGradXX++, bondLevelVelGradXY++, bondLevelVelGradXZ++, 
          bondLevelVelGradYX++, bondLevelVelGradYY++, bondLevelVelGradYZ++,
          bondLevelVelGradZX++, bondLevelVelGradZY++, bondLevelVelGradZZ++){

        neighborIndex = *neighborListPtr;
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVel = velocities + 3*neighborIndex;
        neighborVelGradX = velocityGradientX + 3*neighborIndex;
        neighborVelGradY = velocityGradientY + 3*neighborIndex;
        neighborVelGradZ = velocityGradientZ + 3*neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);

        // The velState is the relative difference in velocities of the nodes at
        // each end of a bond. i.e., v_j - v_i
        velStateX = *(neighborVel)   - *(vel);
        velStateY = *(neighborVel+1) - *(vel+1);
        velStateZ = *(neighborVel+2) - *(vel+2);

        // average of the two points 
        for(int i=0; i<3; i++){
          *(meanVelGrad+i) = 0.5 * (*(velGradX+i) + *(neighborVelGradX+i));
          *(meanVelGrad+i+3) = 0.5 * (*(velGradY+i) + *(neighborVelGradY+i));
          *(meanVelGrad+i+6) = 0.5 * (*(velGradZ+i) + *(neighborVelGradZ+i));
        }

        deformedBondLengthSq = deformedBondLength * deformedBondLength;

        velStateX = *(neighborVel)   - *(vel);
        velStateY = *(neighborVel+1) - *(vel+1);
        velStateZ = *(neighborVel+2) - *(vel+2);

        temp = *(meanVelGrad+0) * deformedBondX + *(meanVelGrad+1) * deformedBondY + *(meanVelGrad+2) * deformedBondZ;
        *bondLevelVelGradXX = *(meanVelGrad+0) + (velStateX - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradXY = *(meanVelGrad+1) + (velStateX - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradXZ = *(meanVelGrad+2) + (velStateX - temp) * deformedBondZ/deformedBondLengthSq;

        temp = *(meanVelGrad+3) * deformedBondX + *(meanVelGrad+4) * deformedBondY + *(meanVelGrad+5) * deformedBondZ;
        *bondLevelVelGradYX = *(meanVelGrad+3) + (velStateY - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradYY = *(meanVelGrad+4) + (velStateY - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradYZ = *(meanVelGrad+5) + (velStateY - temp) * deformedBondZ/deformedBondLengthSq;

        temp = *(meanVelGrad+6) * deformedBondX + *(meanVelGrad+7) * deformedBondY + *(meanVelGrad+8) * deformedBondZ;
        *bondLevelVelGradZX = *(meanVelGrad+6) + (velStateZ - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradZY = *(meanVelGrad+7) + (velStateZ - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradZZ = *(meanVelGrad+8) + (velStateZ - temp) * deformedBondZ/deformedBondLengthSq;
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondLevelVelGradXX += numNeighbors; bondLevelVelGradXY += numNeighbors; bondLevelVelGradXZ += numNeighbors; 
      bondLevelVelGradYX += numNeighbors; bondLevelVelGradYY += numNeighbors; bondLevelVelGradYZ += numNeighbors;
      bondLevelVelGradZX += numNeighbors; bondLevelVelGradZY += numNeighbors; bondLevelVelGradZZ += numNeighbors;
    }
  }
}

//This function updates the node-level deformation gradient based on velocity gradient
template<typename ScalarT>
void updateDeformationGradient
(
    const ScalarT* velocityGradient,
    const ScalarT* deformationGradientN,
    ScalarT* deformationGradientNP1,
    const double* flyingPointFlag,
    int numPoints,
    double dt
)
{
  const ScalarT* velGrad = velocityGradient;
  const ScalarT* defGradN = deformationGradientN;
  ScalarT* defGradNP1 = deformationGradientNP1;

  const double* flyingPointFlg = flyingPointFlag;

  std::vector<ScalarT> FdotVector(9);
  ScalarT* Fdot = &FdotVector[0];

  for(int iID=0 ; iID<numPoints ; ++iID, velGrad+=9, defGradN+=9, defGradNP1+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // L = Fdot . inv(F)
      // Fdot = L . F
      // F_NP1 = F_N + Fdot . dt
      for(int i=0; i<9; i++)
        *(defGradNP1+i) = *(defGradN+i);

      MATRICES::MatrixMultiply(false, false, 1.0, velGrad, defGradN, Fdot);

      for(int i=0; i<9; i++)
        *(defGradNP1+i) += *(Fdot+i) * dt;
    }
  }
}







//Performs kinematic computations following Flanagan and Taylor (1987), returns
//unrotated rate-of-deformation and rotation tensors
//This function computes the node-based values
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
)
{
  int returnCode = 0;

  const ScalarT* velGradXX = bondLevelVelocityGradientXX;
  const ScalarT* velGradXY = bondLevelVelocityGradientXY;
  const ScalarT* velGradXZ = bondLevelVelocityGradientXZ;
  const ScalarT* velGradYX = bondLevelVelocityGradientYX;
  const ScalarT* velGradYY = bondLevelVelocityGradientYY;
  const ScalarT* velGradYZ = bondLevelVelocityGradientYZ;
  const ScalarT* velGradZX = bondLevelVelocityGradientZX;
  const ScalarT* velGradZY = bondLevelVelocityGradientZY;
  const ScalarT* velGradZZ = bondLevelVelocityGradientZZ;
  const ScalarT* leftStretchXXN = bondLevelLeftStretchTensorXXN;
  const ScalarT* leftStretchXYN = bondLevelLeftStretchTensorXYN;
  const ScalarT* leftStretchXZN = bondLevelLeftStretchTensorXZN;
  const ScalarT* leftStretchYXN = bondLevelLeftStretchTensorYXN;
  const ScalarT* leftStretchYYN = bondLevelLeftStretchTensorYYN;
  const ScalarT* leftStretchYZN = bondLevelLeftStretchTensorYZN;
  const ScalarT* leftStretchZXN = bondLevelLeftStretchTensorZXN;
  const ScalarT* leftStretchZYN = bondLevelLeftStretchTensorZYN;
  const ScalarT* leftStretchZZN = bondLevelLeftStretchTensorZZN;
  const ScalarT* rotTensorXXN = bondLevelRotationTensorXXN;
  const ScalarT* rotTensorXYN = bondLevelRotationTensorXYN;
  const ScalarT* rotTensorXZN = bondLevelRotationTensorXZN;
  const ScalarT* rotTensorYXN = bondLevelRotationTensorYXN;
  const ScalarT* rotTensorYYN = bondLevelRotationTensorYYN;
  const ScalarT* rotTensorYZN = bondLevelRotationTensorYZN;
  const ScalarT* rotTensorZXN = bondLevelRotationTensorZXN;
  const ScalarT* rotTensorZYN = bondLevelRotationTensorZYN;
  const ScalarT* rotTensorZZN = bondLevelRotationTensorZZN;

  ScalarT* leftStretchXXNP1 = bondLevelLeftStretchTensorXXNP1;
  ScalarT* leftStretchXYNP1 = bondLevelLeftStretchTensorXYNP1;
  ScalarT* leftStretchXZNP1 = bondLevelLeftStretchTensorXZNP1;
  ScalarT* leftStretchYXNP1 = bondLevelLeftStretchTensorYXNP1;
  ScalarT* leftStretchYYNP1 = bondLevelLeftStretchTensorYYNP1;
  ScalarT* leftStretchYZNP1 = bondLevelLeftStretchTensorYZNP1;
  ScalarT* leftStretchZXNP1 = bondLevelLeftStretchTensorZXNP1;
  ScalarT* leftStretchZYNP1 = bondLevelLeftStretchTensorZYNP1;
  ScalarT* leftStretchZZNP1 = bondLevelLeftStretchTensorZZNP1;
  ScalarT* rotTensorXXNP1 = bondLevelRotationTensorXXNP1;
  ScalarT* rotTensorXYNP1 = bondLevelRotationTensorXYNP1;
  ScalarT* rotTensorXZNP1 = bondLevelRotationTensorXZNP1;
  ScalarT* rotTensorYXNP1 = bondLevelRotationTensorYXNP1;
  ScalarT* rotTensorYYNP1 = bondLevelRotationTensorYYNP1;
  ScalarT* rotTensorYZNP1 = bondLevelRotationTensorYZNP1;
  ScalarT* rotTensorZXNP1 = bondLevelRotationTensorZXNP1;
  ScalarT* rotTensorZYNP1 = bondLevelRotationTensorZYNP1;
  ScalarT* rotTensorZZNP1 = bondLevelRotationTensorZZNP1;
  ScalarT* unrotRateOfDefXX = bondLevelUnrotatedRateOfDeformationXX;
  ScalarT* unrotRateOfDefXY = bondLevelUnrotatedRateOfDeformationXY;
  ScalarT* unrotRateOfDefXZ = bondLevelUnrotatedRateOfDeformationXZ;
  ScalarT* unrotRateOfDefYX = bondLevelUnrotatedRateOfDeformationYX;
  ScalarT* unrotRateOfDefYY = bondLevelUnrotatedRateOfDeformationYY;
  ScalarT* unrotRateOfDefYZ = bondLevelUnrotatedRateOfDeformationYZ;
  ScalarT* unrotRateOfDefZX = bondLevelUnrotatedRateOfDeformationZX;
  ScalarT* unrotRateOfDefZY = bondLevelUnrotatedRateOfDeformationZY;
  ScalarT* unrotRateOfDefZZ = bondLevelUnrotatedRateOfDeformationZZ;

  const double* flyingPointFlg = flyingPointFlag;

  std::vector<ScalarT> eulerianVelGradVector(9) ; ScalarT* eulerianVelGrad = &eulerianVelGradVector[0];
  std::vector<ScalarT> leftStretchNVector(9) ; ScalarT* leftStretchN = &leftStretchNVector[0];
  std::vector<ScalarT> leftStretchNP1Vector(9) ; ScalarT* leftStretchNP1 = &leftStretchNP1Vector[0];
  std::vector<ScalarT> rotTensorNVector(9) ; ScalarT* rotTensorN = &rotTensorNVector[0];
  std::vector<ScalarT> rotTensorNP1Vector(9) ; ScalarT* rotTensorNP1 = &rotTensorNP1Vector[0];
  std::vector<ScalarT> unrotRateOfDefVector(9) ; ScalarT* unrotRateOfDef = &unrotRateOfDefVector[0];

  std::vector<ScalarT> rateOfDefVector(9) ; ScalarT* rateOfDef = &rateOfDefVector[0];
  std::vector<ScalarT> spinVector(9) ; ScalarT* spin = &spinVector[0];
  std::vector<ScalarT> tempVector(9) ; ScalarT* temp = &tempVector[0];
  std::vector<ScalarT> tempInvVector(9) ; ScalarT* tempInv = &tempInvVector[0];
  std::vector<ScalarT> OmegaTensorVector(9) ; ScalarT* OmegaTensor = &OmegaTensorVector[0];
  std::vector<ScalarT> QMatrixVector(9) ; ScalarT* QMatrix = &QMatrixVector[0];
  std::vector<ScalarT> OmegaTensorSqVector(9) ; ScalarT* OmegaTensorSq = &OmegaTensorSqVector[0];
  std::vector<ScalarT> tempAVector(9) ; ScalarT* tempA = &tempAVector[0];
  std::vector<ScalarT> tempBVector(9) ; ScalarT* tempB = &tempBVector[0];
  std::vector<ScalarT> rateOfStretchVector(9) ; ScalarT* rateOfStretch = &rateOfStretchVector[0];

  ScalarT determinant;
  ScalarT omegaX, omegaY, omegaZ;
  ScalarT zX, zY, zZ;
  ScalarT wX, wY, wZ;
  ScalarT traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;
  int inversionReturnCode(0);
  std::string inversionErrorMessage = 
    "**** Error:  computeShapeTensorInverseAndApproximateVelocityGradient: Non-invertible matrix\n";

  int numNeighbors; //neighborIndex
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // All is bond level.
      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
          velGradXX++, velGradXY++, velGradXZ++, 
          velGradYX++, velGradYY++, velGradYZ++, 
          velGradZX++, velGradZY++, velGradZZ++,
          leftStretchXXN++, leftStretchXYN++, leftStretchXZN++, 
          leftStretchYXN++, leftStretchYYN++, leftStretchYZN++, 
          leftStretchZXN++, leftStretchZYN++, leftStretchZZN++,
          rotTensorXXN++, rotTensorXYN++, rotTensorXZN++,
          rotTensorYXN++, rotTensorYYN++, rotTensorYZN++,
          rotTensorZXN++, rotTensorZYN++, rotTensorZZN++,
          leftStretchXXNP1++, leftStretchXYNP1++, leftStretchXZNP1++, 
          leftStretchYXNP1++, leftStretchYYNP1++, leftStretchYZNP1++, 
          leftStretchZXNP1++, leftStretchZYNP1++, leftStretchZZNP1++,
          rotTensorXXNP1++, rotTensorXYNP1++, rotTensorXZNP1++,
          rotTensorYXNP1++, rotTensorYYNP1++, rotTensorYZNP1++,
          rotTensorZXNP1++, rotTensorZYNP1++, rotTensorZZNP1++,
          unrotRateOfDefXX++, unrotRateOfDefXY++, unrotRateOfDefXZ++,
          unrotRateOfDefYX++, unrotRateOfDefYY++, unrotRateOfDefYZ++,
          unrotRateOfDefZX++, unrotRateOfDefZY++, unrotRateOfDefZZ++){

        // neighborIndex = *neighborListPtr;

        // Store in a tensor form 
        *(eulerianVelGrad+0) = *velGradXX; *(eulerianVelGrad+1) = *velGradXY; *(eulerianVelGrad+2) = *velGradXZ;
        *(eulerianVelGrad+3) = *velGradYX; *(eulerianVelGrad+4) = *velGradYY; *(eulerianVelGrad+5) = *velGradYZ;
        *(eulerianVelGrad+6) = *velGradZX; *(eulerianVelGrad+7) = *velGradZY; *(eulerianVelGrad+8) = *velGradZZ;
        *(leftStretchN+0) = *leftStretchXXN; *(leftStretchN+1) = *leftStretchXYN; *(leftStretchN+2) = *leftStretchXZN;
        *(leftStretchN+3) = *leftStretchYXN; *(leftStretchN+4) = *leftStretchYYN; *(leftStretchN+5) = *leftStretchYZN;
        *(leftStretchN+6) = *leftStretchZXN; *(leftStretchN+7) = *leftStretchZYN; *(leftStretchN+8) = *leftStretchZZN;
        *(rotTensorN+0) = *rotTensorXXN; *(rotTensorN+1) = *rotTensorXYN; *(rotTensorN+2) = *rotTensorXZN;
        *(rotTensorN+3) = *rotTensorYXN; *(rotTensorN+4) = *rotTensorYYN; *(rotTensorN+5) = *rotTensorYZN;
        *(rotTensorN+6) = *rotTensorZXN; *(rotTensorN+7) = *rotTensorZYN; *(rotTensorN+8) = *rotTensorZZN;

        // Compute rate-of-deformation tensor, D = 1/2 * (L + Lt)
        *(rateOfDef)   = *(eulerianVelGrad);
        *(rateOfDef+1) = 0.5 * ( *(eulerianVelGrad+1) + *(eulerianVelGrad+3) );
        *(rateOfDef+2) = 0.5 * ( *(eulerianVelGrad+2) + *(eulerianVelGrad+6) );
        *(rateOfDef+3) = *(rateOfDef+1);
        *(rateOfDef+4) = *(eulerianVelGrad+4);
        *(rateOfDef+5) = 0.5 * ( *(eulerianVelGrad+5) + *(eulerianVelGrad+7) );
        *(rateOfDef+6) = *(rateOfDef+2);
        *(rateOfDef+7) = *(rateOfDef+5);
        *(rateOfDef+8) = *(eulerianVelGrad+8);

        // Compute spin tensor, W = 1/2 * (L - Lt)
        *(spin)   = 0.0;
        *(spin+1) = 0.5 * ( *(eulerianVelGrad+1) - *(eulerianVelGrad+3) );
        *(spin+2) = 0.5 * ( *(eulerianVelGrad+2) - *(eulerianVelGrad+6) );
        *(spin+3) = -1.0 * *(spin+1);
        *(spin+4) = 0.0;
        *(spin+5) = 0.5 * ( *(eulerianVelGrad+5) - *(eulerianVelGrad+7) );
        *(spin+6) = -1.0 * *(spin+2);
        *(spin+7) = -1.0 * *(spin+5);
        *(spin+8) = 0.0;
       
        //Following Flanagan & Taylor (T&F) 
        //
        //Find the vector z_i = \epsilon_{ikj} * D_{jm} * V_{mk} (T&F Eq. 13)
        //
        //where \epsilon_{ikj} is the alternator tensor.
        //
        //Components below copied from computer algebra solution to the expansion
        //above
        zX = - *(leftStretchN+2) *  *(rateOfDef+3) -  *(leftStretchN+5) *  *(rateOfDef+4) - 
               *(leftStretchN+8) *  *(rateOfDef+5) +  *(leftStretchN+1) *  *(rateOfDef+6) + 
               *(leftStretchN+4) *  *(rateOfDef+7) +  *(leftStretchN+7) *  *(rateOfDef+8);
        zY =   *(leftStretchN+2) *  *(rateOfDef)   +  *(leftStretchN+5) *  *(rateOfDef+1) + 
               *(leftStretchN+8) *  *(rateOfDef+2) -  *(leftStretchN)   *  *(rateOfDef+6) - 
               *(leftStretchN+3) *  *(rateOfDef+7) -  *(leftStretchN+6) *  *(rateOfDef+8);
        zZ = - *(leftStretchN+1) *  *(rateOfDef)   -  *(leftStretchN+4) *  *(rateOfDef+1) - 
               *(leftStretchN+7) *  *(rateOfDef+2) +  *(leftStretchN)   *  *(rateOfDef+3) + 
               *(leftStretchN+3) *  *(rateOfDef+4) +  *(leftStretchN+6) *  *(rateOfDef+5);

        //Find the vector w_i = -1/2 * \epsilon_{ijk} * W_{jk} (T&F Eq. 11)
        wX = 0.5 * ( *(spin+7) - *(spin+5) );
        wY = 0.5 * ( *(spin+2) - *(spin+6) );
        wZ = 0.5 * ( *(spin+3) - *(spin+1) );

        //Find trace(V)
        traceV = *(leftStretchN) + *(leftStretchN+4) + *(leftStretchN+8);

        // Compute (trace(V) * I - V) store in temp
        *(temp)   = traceV - *(leftStretchN);
        *(temp+1) = - *(leftStretchN+1);
        *(temp+2) = - *(leftStretchN+2);
        *(temp+3) = - *(leftStretchN+3);
        *(temp+4) = traceV - *(leftStretchN+4);
        *(temp+5) = - *(leftStretchN+5);
        *(temp+6) = - *(leftStretchN+6);
        *(temp+7) = - *(leftStretchN+7);
        *(temp+8) = traceV - *(leftStretchN+8);

        // Compute the inverse of the temp matrix
        MATRICES::Invert3by3Matrix(temp, determinant, tempInv);
        if(inversionReturnCode > 0){
          returnCode = inversionReturnCode;
          std::cout << inversionErrorMessage;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //Find omega vector, i.e. \omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
        omegaX =  wX + *(tempInv)   * zX + *(tempInv+1) * zY + *(tempInv+2) * zZ;
        omegaY =  wY + *(tempInv+3) * zX + *(tempInv+4) * zY + *(tempInv+5) * zZ;
        omegaZ =  wZ + *(tempInv+6) * zX + *(tempInv+7) * zY + *(tempInv+8) * zZ;

        //Find the tensor \Omega_{ij} = \epsilon_{ikj} * w_k (T&F Eq. 10)
        *(OmegaTensor) = 0.0;
        *(OmegaTensor+1) = -omegaZ;
        *(OmegaTensor+2) = omegaY;
        *(OmegaTensor+3) = omegaZ;
        *(OmegaTensor+4) = 0.0;
        *(OmegaTensor+5) = -omegaX;
        *(OmegaTensor+6) = -omegaY;
        *(OmegaTensor+7) = omegaX;
        *(OmegaTensor+8) = 0.0;

        //Increment R with (T&F Eq. 36 and 44) as opposed to solving (T&F 39) this
        //is desirable for accuracy in implicit solves and has no effect on
        //explicit solves (other than a slight decrease in speed).
        //
        // Compute Q with (T&F Eq. 44)
        //
        // Omega^2 = w_i * w_i (T&F Eq. 42)
        OmegaSq = omegaX*omegaX + omegaY*omegaY + omegaZ*omegaZ;
        // Omega = \sqrt{OmegaSq}
        Omega = sqrt(OmegaSq);

        // Avoid a potential divide-by-zero
        if(OmegaSq > 1.e-30){

          // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
          //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
          scaleFactor1 = sin(dt*Omega) / Omega;
          scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
          MATRICES::MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, OmegaTensorSq);
          *(QMatrix)   = 1.0 + scaleFactor1 * *(OmegaTensor)   + scaleFactor2 * *(OmegaTensorSq)   ;
          *(QMatrix+1) =       scaleFactor1 * *(OmegaTensor+1) + scaleFactor2 * *(OmegaTensorSq+1) ;
          *(QMatrix+2) =       scaleFactor1 * *(OmegaTensor+2) + scaleFactor2 * *(OmegaTensorSq+2) ;
          *(QMatrix+3) =       scaleFactor1 * *(OmegaTensor+3) + scaleFactor2 * *(OmegaTensorSq+3) ;
          *(QMatrix+4) = 1.0 + scaleFactor1 * *(OmegaTensor+4) + scaleFactor2 * *(OmegaTensorSq+4) ;
          *(QMatrix+5) =       scaleFactor1 * *(OmegaTensor+5) + scaleFactor2 * *(OmegaTensorSq+5) ;
          *(QMatrix+6) =       scaleFactor1 * *(OmegaTensor+6) + scaleFactor2 * *(OmegaTensorSq+6) ;
          *(QMatrix+7) =       scaleFactor1 * *(OmegaTensor+7) + scaleFactor2 * *(OmegaTensorSq+7) ;
          *(QMatrix+8) = 1.0 + scaleFactor1 * *(OmegaTensor+8) + scaleFactor2 * *(OmegaTensorSq+8) ;

        } else {
          *(QMatrix)   = 1.0 ; *(QMatrix+1) = 0.0 ; *(QMatrix+2) = 0.0 ;
          *(QMatrix+3) = 0.0 ; *(QMatrix+4) = 1.0 ; *(QMatrix+5) = 0.0 ;
          *(QMatrix+6) = 0.0 ; *(QMatrix+7) = 0.0 ; *(QMatrix+8) = 1.0 ;
        };

        // Compute R_STEP_NP1 = QMatrix * R_STEP_N (T&F Eq. 36)
        MATRICES::MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, rotTensorNP1);

        // Compute rate of stretch, Vdot = L*V - V*Omega
        // First tempA = L*V, 
        MATRICES::MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, tempA);

        // tempB = V*Omega
        MATRICES::MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, tempB);

        //Vdot = tempA - tempB
        for(int i=0 ; i<9 ; ++i)
          *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

        //V_STEP_NP1 = V_STEP_N + dt*Vdot
        for(int i=0 ; i<9 ; ++i)
          *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

        // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
        MATRICES::MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, temp);

        // d = Rt * temp
        MATRICES::MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, unrotRateOfDef);

        // Store back in element-wise format
        *leftStretchXXNP1 = *(leftStretchNP1+0); *leftStretchXYNP1 = *(leftStretchNP1+1); *leftStretchXZNP1 = *(leftStretchNP1+2);
        *leftStretchYXNP1 = *(leftStretchNP1+3); *leftStretchYYNP1 = *(leftStretchNP1+4); *leftStretchYZNP1 = *(leftStretchNP1+5);
        *leftStretchZXNP1 = *(leftStretchNP1+6); *leftStretchZYNP1 = *(leftStretchNP1+7); *leftStretchZZNP1 = *(leftStretchNP1+8);
        *rotTensorXXNP1 = *(rotTensorNP1+0); *rotTensorXYNP1 = *(rotTensorNP1+1); *rotTensorXZNP1 = *(rotTensorNP1+2);
        *rotTensorYXNP1 = *(rotTensorNP1+3); *rotTensorYYNP1 = *(rotTensorNP1+4); *rotTensorYZNP1 = *(rotTensorNP1+5);
        *rotTensorZXNP1 = *(rotTensorNP1+6); *rotTensorZYNP1 = *(rotTensorNP1+7); *rotTensorZZNP1 = *(rotTensorNP1+8);
        *unrotRateOfDefXX = *(unrotRateOfDef+0); *unrotRateOfDefXY = *(unrotRateOfDef+1); *unrotRateOfDefXZ = *(unrotRateOfDef+2);
        *unrotRateOfDefYX = *(unrotRateOfDef+3); *unrotRateOfDefYY = *(unrotRateOfDef+4); *unrotRateOfDefYZ = *(unrotRateOfDef+5);
        *unrotRateOfDefZX = *(unrotRateOfDef+6); *unrotRateOfDefZY = *(unrotRateOfDef+7); *unrotRateOfDefZZ = *(unrotRateOfDef+8);
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      velGradXX += numNeighbors; velGradXY += numNeighbors; velGradXZ += numNeighbors; 
      velGradYX += numNeighbors; velGradYY += numNeighbors; velGradYZ += numNeighbors; 
      velGradZX += numNeighbors; velGradZY += numNeighbors; velGradZZ += numNeighbors;
      leftStretchXXN += numNeighbors; leftStretchXYN += numNeighbors; leftStretchXZN += numNeighbors; 
      leftStretchYXN += numNeighbors; leftStretchYYN += numNeighbors; leftStretchYZN += numNeighbors; 
      leftStretchZXN += numNeighbors; leftStretchZYN += numNeighbors; leftStretchZZN += numNeighbors;
      rotTensorXXN += numNeighbors; rotTensorXYN += numNeighbors; rotTensorXZN += numNeighbors;
      rotTensorYXN += numNeighbors; rotTensorYYN += numNeighbors; rotTensorYZN += numNeighbors;
      rotTensorZXN += numNeighbors; rotTensorZYN += numNeighbors; rotTensorZZN += numNeighbors;
      leftStretchXXNP1 += numNeighbors; leftStretchXYNP1 += numNeighbors; leftStretchXZNP1 += numNeighbors; 
      leftStretchYXNP1 += numNeighbors; leftStretchYYNP1 += numNeighbors; leftStretchYZNP1 += numNeighbors; 
      leftStretchZXNP1 += numNeighbors; leftStretchZYNP1 += numNeighbors; leftStretchZZNP1 += numNeighbors;
      rotTensorXXNP1 += numNeighbors; rotTensorXYNP1 += numNeighbors; rotTensorXZNP1 += numNeighbors;
      rotTensorYXNP1 += numNeighbors; rotTensorYYNP1 += numNeighbors; rotTensorYZNP1 += numNeighbors;
      rotTensorZXNP1 += numNeighbors; rotTensorZYNP1 += numNeighbors; rotTensorZZNP1 += numNeighbors;
      unrotRateOfDefXX += numNeighbors; unrotRateOfDefXY += numNeighbors; unrotRateOfDefXZ += numNeighbors;
      unrotRateOfDefYX += numNeighbors; unrotRateOfDefYY += numNeighbors; unrotRateOfDefYZ += numNeighbors;
      unrotRateOfDefZX += numNeighbors; unrotRateOfDefZY += numNeighbors; unrotRateOfDefZZ += numNeighbors;
    }
  }

  return returnCode;
}


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
)
{
  const ScalarT* rotTensorXX = bondLevelRotationTensorXX;
  const ScalarT* rotTensorXY = bondLevelRotationTensorXY;
  const ScalarT* rotTensorXZ = bondLevelRotationTensorXZ;
  const ScalarT* rotTensorYX = bondLevelRotationTensorYX;
  const ScalarT* rotTensorYY = bondLevelRotationTensorYY;
  const ScalarT* rotTensorYZ = bondLevelRotationTensorYZ;
  const ScalarT* rotTensorZX = bondLevelRotationTensorZX;
  const ScalarT* rotTensorZY = bondLevelRotationTensorZY;
  const ScalarT* rotTensorZZ = bondLevelRotationTensorZZ;
  const ScalarT* unrotatedStressXX = bondLevelUnrotatedCauchyStressXX;
  const ScalarT* unrotatedStressXY = bondLevelUnrotatedCauchyStressXY;
  const ScalarT* unrotatedStressXZ = bondLevelUnrotatedCauchyStressXZ;
  const ScalarT* unrotatedStressYX = bondLevelUnrotatedCauchyStressYX;
  const ScalarT* unrotatedStressYY = bondLevelUnrotatedCauchyStressYY;
  const ScalarT* unrotatedStressYZ = bondLevelUnrotatedCauchyStressYZ;
  const ScalarT* unrotatedStressZX = bondLevelUnrotatedCauchyStressZX;
  const ScalarT* unrotatedStressZY = bondLevelUnrotatedCauchyStressZY;
  const ScalarT* unrotatedStressZZ = bondLevelUnrotatedCauchyStressZZ;
  ScalarT* rotatedStressXX = bondLevelRotatedCauchyStressXX;
  ScalarT* rotatedStressXY = bondLevelRotatedCauchyStressXY;
  ScalarT* rotatedStressXZ = bondLevelRotatedCauchyStressXZ;
  ScalarT* rotatedStressYX = bondLevelRotatedCauchyStressYX;
  ScalarT* rotatedStressYY = bondLevelRotatedCauchyStressYY;
  ScalarT* rotatedStressYZ = bondLevelRotatedCauchyStressYZ;
  ScalarT* rotatedStressZX = bondLevelRotatedCauchyStressZX;
  ScalarT* rotatedStressZY = bondLevelRotatedCauchyStressZY;
  ScalarT* rotatedStressZZ = bondLevelRotatedCauchyStressZZ;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT unrotatedStress[9], rotTensor[9], rotatedStress[9], temp[9];

  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // All is bond level.
      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++,
          rotTensorXX++, rotTensorXY++, rotTensorXZ++, 
          rotTensorYX++, rotTensorYY++, rotTensorYZ++, 
          rotTensorZX++, rotTensorZY++, rotTensorZZ++, 
          unrotatedStressXX++, unrotatedStressXY++, unrotatedStressXZ++, 
          unrotatedStressYX++, unrotatedStressYY++, unrotatedStressYZ++, 
          unrotatedStressZX++, unrotatedStressZY++, unrotatedStressZZ++, 
          rotatedStressXX++, rotatedStressXY++, rotatedStressXZ++, 
          rotatedStressYX++, rotatedStressYY++, rotatedStressYZ++, 
          rotatedStressZX++, rotatedStressZY++, rotatedStressZZ++){

        // write in matrix form 
        rotTensor[0] = *rotTensorXX; rotTensor[1] = *rotTensorXY; rotTensor[2] = *rotTensorXZ;
        rotTensor[3] = *rotTensorYX; rotTensor[4] = *rotTensorYY; rotTensor[5] = *rotTensorYZ;
        rotTensor[6] = *rotTensorZX; rotTensor[7] = *rotTensorZY; rotTensor[8] = *rotTensorZZ;
        unrotatedStress[0] = *unrotatedStressXX; unrotatedStress[1] = *unrotatedStressXY; unrotatedStress[2] = *unrotatedStressXZ;
        unrotatedStress[3] = *unrotatedStressYX; unrotatedStress[4] = *unrotatedStressYY; unrotatedStress[5] = *unrotatedStressYZ;
        unrotatedStress[6] = *unrotatedStressZX; unrotatedStress[7] = *unrotatedStressZY; unrotatedStress[8] = *unrotatedStressZZ;

        // temp = \sigma_unrot * Rt
         MATRICES::MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, temp);
        // \sigma_rot = R * temp
         MATRICES::MatrixMultiply(false, false, 1.0, rotTensor, temp, rotatedStress);

        // update bond-level data field
        *rotatedStressXX = rotatedStress[0]; *rotatedStressXY = rotatedStress[1]; *rotatedStressXZ = rotatedStress[2]; 
        *rotatedStressYX = rotatedStress[3]; *rotatedStressYY = rotatedStress[4]; *rotatedStressYZ = rotatedStress[5]; 
        *rotatedStressZX = rotatedStress[6]; *rotatedStressZY = rotatedStress[7]; *rotatedStressZZ = rotatedStress[8]; 
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      rotTensorXX += numNeighbors; rotTensorXY += numNeighbors; rotTensorXZ += numNeighbors; 
      rotTensorYX += numNeighbors; rotTensorYY += numNeighbors; rotTensorYZ += numNeighbors; 
      rotTensorZX += numNeighbors; rotTensorZY += numNeighbors; rotTensorZZ += numNeighbors; 
      unrotatedStressXX += numNeighbors; unrotatedStressXY += numNeighbors; unrotatedStressXZ += numNeighbors; 
      unrotatedStressYX += numNeighbors; unrotatedStressYY += numNeighbors; unrotatedStressYZ += numNeighbors; 
      unrotatedStressZX += numNeighbors; unrotatedStressZY += numNeighbors; unrotatedStressZZ += numNeighbors; 
      rotatedStressXX += numNeighbors; rotatedStressXY += numNeighbors; rotatedStressXZ += numNeighbors; 
      rotatedStressYX += numNeighbors; rotatedStressYY += numNeighbors; rotatedStressYZ += numNeighbors; 
      rotatedStressZX += numNeighbors; rotatedStressZY += numNeighbors; rotatedStressZZ += numNeighbors;
    }
  }
}

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
    int  numPoints
)
{
  const double* delta = horizon;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const double* w0 = weightedVolume;
  const double* neighborW0;
  const ScalarT* stressXX = bondLevelCauchyStressXX;
  const ScalarT* stressXY = bondLevelCauchyStressXY;
  const ScalarT* stressXZ = bondLevelCauchyStressXZ;
  const ScalarT* stressYX = bondLevelCauchyStressYX;
  const ScalarT* stressYY = bondLevelCauchyStressYY;
  const ScalarT* stressYZ = bondLevelCauchyStressYZ;
  const ScalarT* stressZX = bondLevelCauchyStressZX;
  const ScalarT* stressZY = bondLevelCauchyStressZY;
  const ScalarT* stressZZ = bondLevelCauchyStressZZ;
  ScalarT* integral = nonhomogeneousIntegral;

  const double* flyingPointFlg = flyingPointFlag;
  const double *bondDamagePtr = bondDamage;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength, deformedBondLengthSq;
  double neighborVolume, omega, scalarTemp;

  std::vector<ScalarT> tempVector(9);
  ScalarT* temp = &tempVector[0];

  std::vector<ScalarT> stressVector(9);
  ScalarT* stress = &stressVector[0];

  std::vector<ScalarT> integrandVector(9);
  ScalarT* integrand = &integrandVector[0];

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3, w0++, integral+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // Zero out data
      *(integral)   = 0.0 ; *(integral+1) = 0.0 ; *(integral+2) = 0.0 ;
      *(integral+3) = 0.0 ; *(integral+4) = 0.0 ; *(integral+5) = 0.0 ;
      *(integral+6) = 0.0 ; *(integral+7) = 0.0 ; *(integral+8) = 0.0 ;

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++,
          stressXX++, stressXY++, stressXZ++, 
          stressYX++, stressYY++, stressYZ++, 
          stressZX++, stressZY++, stressZZ++){

        neighborIndex = *neighborListPtr;
        neighborVolume = jacobianDeterminant[neighborIndex] * volume[neighborIndex];
        neighborCoord = coordinates + 3*neighborIndex;
        neighborW0 = weightedVolume + neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);
        deformedBondLengthSq = deformedBondX*deformedBondX +
                               deformedBondY*deformedBondY +
                               deformedBondZ*deformedBondZ;

        // write the stress in matrix form 
        stress[0] = *stressXX; stress[1] = *stressXY; stress[2] = *stressXZ; 
        stress[3] = *stressYX; stress[4] = *stressYY; stress[5] = *stressYZ; 
        stress[6] = *stressZX; stress[7] = *stressZY; stress[8] = *stressZZ; 

        // delta_jp - (y_j y_p)/|y|^2
        *(temp+0) = 1.0 - deformedBondX * deformedBondX / deformedBondLengthSq;
        *(temp+1) = - deformedBondX * deformedBondY / deformedBondLengthSq;
        *(temp+2) = - deformedBondX * deformedBondZ / deformedBondLengthSq;
        *(temp+3) = *(temp+1);
        *(temp+4) = 1.0 - deformedBondY * deformedBondY / deformedBondLengthSq;
        *(temp+5) = - deformedBondY * deformedBondZ / deformedBondLengthSq;
        *(temp+6) = *(temp+2);
        *(temp+7) = *(temp+5);
        *(temp+8) = 1.0 - deformedBondZ * deformedBondZ / deformedBondLengthSq;

        // Matrix multiply the stress and the second term to compute the integrand
        MATRICES::MatrixMultiply(false, false, 1.0, stress, temp, integrand);

        omega = (1.0 - *bondDamagePtr) * MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

        if(omega > 0.0){
          scalarTemp = omega * (0.5 / *w0 + 0.5 / *neighborW0) * neighborVolume;

          for(int i=0; i<9; i++)
            *(integral+i) += scalarTemp * *(integrand+i);
        }
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondDamagePtr += numNeighbors; 
      stressXX += numNeighbors; stressXY += numNeighbors; stressXZ += numNeighbors; 
      stressYX += numNeighbors; stressYY += numNeighbors; stressYZ += numNeighbors; 
      stressZX += numNeighbors; stressZY += numNeighbors; stressZZ += numNeighbors;
    }
  }
}



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
    double dt,
    bool m_plane
)
{
  int returnCode = 0;

  const ScalarT* defGradX = deformationGradientX;
  const ScalarT* defGradY = deformationGradientY;
  const ScalarT* defGradZ = deformationGradientZ;
  const ScalarT* defGradDotX = deformationGradientDotX;
  const ScalarT* defGradDotY = deformationGradientDotY;
  const ScalarT* defGradDotZ = deformationGradientDotZ;
  const ScalarT* leftStretchN = leftStretchTensorN;
  const ScalarT* rotTensorN = rotationTensorN;

  ScalarT* leftStretchNP1 = leftStretchTensorNP1;
  ScalarT* rotTensorNP1 = rotationTensorNP1;
  ScalarT* unrotRateOfDef = unrotatedRateOfDeformation;

  std::vector<ScalarT> defGradVector(9) ; ScalarT* defGrad = &defGradVector[0];
  std::vector<ScalarT> defGradDotVector(9) ; ScalarT* defGradDot = &defGradDotVector[0];
  std::vector<ScalarT> FinverseVector(9) ; ScalarT* Finverse = &FinverseVector[0];
  std::vector<ScalarT> eulerianVelGradVector(9) ; ScalarT* eulerianVelGrad = &eulerianVelGradVector[0];
  std::vector<ScalarT> rateOfDefVector(9) ; ScalarT* rateOfDef = &rateOfDefVector[0];
  std::vector<ScalarT> spinVector(9) ; ScalarT* spin = &spinVector[0];
  std::vector<ScalarT> tempVector(9) ; ScalarT* temp = &tempVector[0];
  std::vector<ScalarT> tempInvVector(9) ; ScalarT* tempInv = &tempInvVector[0];
  std::vector<ScalarT> OmegaTensorVector(9) ; ScalarT* OmegaTensor = &OmegaTensorVector[0];
  std::vector<ScalarT> QMatrixVector(9) ; ScalarT* QMatrix = &QMatrixVector[0];
  std::vector<ScalarT> OmegaTensorSqVector(9) ; ScalarT* OmegaTensorSq = &OmegaTensorSqVector[0];
  std::vector<ScalarT> tempAVector(9) ; ScalarT* tempA = &tempAVector[0];
  std::vector<ScalarT> tempBVector(9) ; ScalarT* tempB = &tempBVector[0];
  std::vector<ScalarT> rateOfStretchVector(9) ; ScalarT* rateOfStretch = &rateOfStretchVector[0];

  ScalarT determinant;
  ScalarT omegaX, omegaY, omegaZ;
  ScalarT zX, zY, zZ;
  ScalarT wX, wY, wZ;
  ScalarT traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;
  int inversionReturnCode(0);

  // placeholder for bond damage
  // double bondDamage = 0.0;

  for(int iID=0 ; iID<numPoints ; ++iID, defGradX+=3, defGradY+=3, defGradZ+=3, 
      defGradDotX+=3, defGradDotY+=3, defGradDotZ+=3, rotTensorN+=9, 
      rotTensorNP1+=9, leftStretchNP1+=9, leftStretchN+=9, unrotRateOfDef+=9){

    for(int i=0; i<3; i++){
      *(defGrad+i+0) = *(defGradX+i);
      *(defGrad+i+3) = *(defGradY+i);
      *(defGrad+i+6) = *(defGradZ+i);
      *(defGradDot+i+0) = *(defGradDotX+i);
      *(defGradDot+i+3) = *(defGradDotY+i);
      *(defGradDot+i+6) = *(defGradDotZ+i);
    }

    // Compute the inverse of the deformation gradient, Finverse
    inversionReturnCode = MATRICES::Invert3by3Matrix(defGrad, determinant, Finverse);
    if(inversionReturnCode > 0)
      returnCode = inversionReturnCode;

    // Compute the Eulerian velocity gradient L = Fdot * Finv
    MATRICES::MatrixMultiply(false, false, 1.0, defGradDot, Finverse, eulerianVelGrad);

    // Compute rate-of-deformation tensor, D = 1/2 * (L + Lt)
    *(rateOfDef)   = *(eulerianVelGrad);
    *(rateOfDef+1) = 0.5 * ( *(eulerianVelGrad+1) + *(eulerianVelGrad+3) );
    *(rateOfDef+2) = 0.5 * ( *(eulerianVelGrad+2) + *(eulerianVelGrad+6) );
    *(rateOfDef+3) = *(rateOfDef+1);
    *(rateOfDef+4) = *(eulerianVelGrad+4);
    *(rateOfDef+5) = 0.5 * ( *(eulerianVelGrad+5) + *(eulerianVelGrad+7) );
    *(rateOfDef+6) = *(rateOfDef+2);
    *(rateOfDef+7) = *(rateOfDef+5);
    *(rateOfDef+8) = *(eulerianVelGrad+8);

    // Compute spin tensor, W = 1/2 * (L - Lt)
    *(spin)   = 0.0;
    *(spin+1) = 0.5 * ( *(eulerianVelGrad+1) - *(eulerianVelGrad+3) );
    *(spin+2) = 0.5 * ( *(eulerianVelGrad+2) - *(eulerianVelGrad+6) );
    *(spin+3) = -1.0 * *(spin+1);
    *(spin+4) = 0.0;
    *(spin+5) = 0.5 * ( *(eulerianVelGrad+5) - *(eulerianVelGrad+7) );
    *(spin+6) = -1.0 * *(spin+2);
    *(spin+7) = -1.0 * *(spin+5);
    *(spin+8) = 0.0;
   
    //Following Flanagan & Taylor (T&F) 
    //
    //Find the vector z_i = \epsilon_{ikj} * D_{jm} * V_{mk} (T&F Eq. 13)
    //
    //where \epsilon_{ikj} is the alternator tensor.
    //
    //Components below copied from computer algebra solution to the expansion
    //above
    
    
    zX = - *(leftStretchN+2) *  *(rateOfDef+3) -  *(leftStretchN+5) *  *(rateOfDef+4) - 
           *(leftStretchN+8) *  *(rateOfDef+5) +  *(leftStretchN+1) *  *(rateOfDef+6) + 
           *(leftStretchN+4) *  *(rateOfDef+7) +  *(leftStretchN+7) *  *(rateOfDef+8);
    zY =   *(leftStretchN+2) *  *(rateOfDef)   +  *(leftStretchN+5) *  *(rateOfDef+1) + 
           *(leftStretchN+8) *  *(rateOfDef+2) -  *(leftStretchN)   *  *(rateOfDef+6) - 
           *(leftStretchN+3) *  *(rateOfDef+7) -  *(leftStretchN+6) *  *(rateOfDef+8);
    zZ = - *(leftStretchN+1) *  *(rateOfDef)   -  *(leftStretchN+4) *  *(rateOfDef+1) - 
           *(leftStretchN+7) *  *(rateOfDef+2) +  *(leftStretchN)   *  *(rateOfDef+3) + 
           *(leftStretchN+3) *  *(rateOfDef+4) +  *(leftStretchN+6) *  *(rateOfDef+5);

    //Find the vector w_i = -1/2 * \epsilon_{ijk} * W_{jk} (T&F Eq. 11)
    wX = 0.5 * ( *(spin+7) - *(spin+5) );
    wY = 0.5 * ( *(spin+2) - *(spin+6) );
    wZ = 0.5 * ( *(spin+3) - *(spin+1) );

    //Find trace(V)
    traceV = *(leftStretchN) + *(leftStretchN+4) + *(leftStretchN+8);

    // Compute (trace(V) * I - V) store in temp
    *(temp)   = traceV - *(leftStretchN);
    *(temp+1) = - *(leftStretchN+1);
    *(temp+2) = - *(leftStretchN+2);
    *(temp+3) = - *(leftStretchN+3);
    *(temp+4) = traceV - *(leftStretchN+4);
    *(temp+5) = - *(leftStretchN+5);
    *(temp+6) = - *(leftStretchN+6);
    *(temp+7) = - *(leftStretchN+7);
    *(temp+8) = traceV - *(leftStretchN+8);

    // Compute the inverse of the temp matrix
    MATRICES::Invert3by3Matrix(temp, determinant, tempInv);
    if(inversionReturnCode > 0)
      returnCode = inversionReturnCode;

    //Find omega vector, i.e. \omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
    omegaX =  wX + *(tempInv)   * zX + *(tempInv+1) * zY + *(tempInv+2) * zZ;
    omegaY =  wY + *(tempInv+3) * zX + *(tempInv+4) * zY + *(tempInv+5) * zZ;
    omegaZ =  wZ + *(tempInv+6) * zX + *(tempInv+7) * zY + *(tempInv+8) * zZ;

    //Find the tensor \Omega_{ij} = \epsilon_{ikj} * w_k (T&F Eq. 10)
    *(OmegaTensor) = 0.0;
    *(OmegaTensor+1) = -omegaZ;
    *(OmegaTensor+2) = omegaY;
    *(OmegaTensor+3) = omegaZ;
    *(OmegaTensor+4) = 0.0;
    *(OmegaTensor+5) = -omegaX;
    *(OmegaTensor+6) = -omegaY;
    *(OmegaTensor+7) = omegaX;
    *(OmegaTensor+8) = 0.0;

    //Increment R with (T&F Eq. 36 and 44) as opposed to solving (T&F 39) this
    //is desirable for accuracy in implicit solves and has no effect on
    //explicit solves (other than a slight decrease in speed).
    //
    // Compute Q with (T&F Eq. 44)
    //
    // Omega^2 = w_i * w_i (T&F Eq. 42)
    OmegaSq = omegaX*omegaX + omegaY*omegaY + omegaZ*omegaZ;
    // Omega = \sqrt{OmegaSq}
    Omega = sqrt(OmegaSq);

    // Avoid a potential divide-by-zero
    if( OmegaSq > 1.e-30){

      // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
      //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
      scaleFactor1 = sin(dt*Omega) / Omega;
      scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
      MATRICES::MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, OmegaTensorSq);
      *(QMatrix)   = 1.0 + scaleFactor1 * *(OmegaTensor)   + scaleFactor2 * *(OmegaTensorSq)   ;
      *(QMatrix+1) =       scaleFactor1 * *(OmegaTensor+1) + scaleFactor2 * *(OmegaTensorSq+1) ;
      *(QMatrix+2) =       scaleFactor1 * *(OmegaTensor+2) + scaleFactor2 * *(OmegaTensorSq+2) ;
      *(QMatrix+3) =       scaleFactor1 * *(OmegaTensor+3) + scaleFactor2 * *(OmegaTensorSq+3) ;
      *(QMatrix+4) = 1.0 + scaleFactor1 * *(OmegaTensor+4) + scaleFactor2 * *(OmegaTensorSq+4) ;
      *(QMatrix+5) =       scaleFactor1 * *(OmegaTensor+5) + scaleFactor2 * *(OmegaTensorSq+5) ;
      *(QMatrix+6) =       scaleFactor1 * *(OmegaTensor+6) + scaleFactor2 * *(OmegaTensorSq+6) ;
      *(QMatrix+7) =       scaleFactor1 * *(OmegaTensor+7) + scaleFactor2 * *(OmegaTensorSq+7) ;
      *(QMatrix+8) = 1.0 + scaleFactor1 * *(OmegaTensor+8) + scaleFactor2 * *(OmegaTensorSq+8) ;

    } else {
      *(QMatrix)   = 1.0 ; *(QMatrix+1) = 0.0 ; *(QMatrix+2) = 0.0 ;
      *(QMatrix+3) = 0.0 ; *(QMatrix+4) = 1.0 ; *(QMatrix+5) = 0.0 ;
      *(QMatrix+6) = 0.0 ; *(QMatrix+7) = 0.0 ; *(QMatrix+8) = 1.0 ;
    };

    // Compute R_STEP_NP1 = QMatrix * R_STEP_N (T&F Eq. 36)
    MATRICES::MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, rotTensorNP1);

    // Compute rate of stretch, Vdot = L*V - V*Omega
    // First tempA = L*V, 
    MATRICES::MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, tempA);

    // tempB = V*Omega
    MATRICES::MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, tempB);

    //Vdot = tempA - tempB
    for(int i=0 ; i<9 ; ++i)
      *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

    //V_STEP_NP1 = V_STEP_N + dt*Vdot
    for(int i=0 ; i<9 ; ++i)
      *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

    // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
    MATRICES::MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, temp);

    // d = Rt * temp
    MATRICES::MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, unrotRateOfDef);
  }

  return returnCode;
}


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
)
{
  int returnCode = 0;

  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  const ScalarT* defGradX = deformationGradientX;
  const ScalarT* neighborDefGradX;
  const ScalarT* defGradY = deformationGradientY;
  const ScalarT* neighborDefGradY;
  const ScalarT* defGradZ = deformationGradientZ;
  const ScalarT* neighborDefGradZ;
  const ScalarT* defGradDotX = deformationGradientDotX;
  const ScalarT* neighborDefGradDotX;
  const ScalarT* defGradDotY = deformationGradientDotY;
  const ScalarT* neighborDefGradDotY;
  const ScalarT* defGradDotZ = deformationGradientDotZ;
  const ScalarT* neighborDefGradDotZ;
  ScalarT* J = bondLevelJacobianDeterminant;
  ScalarT* defGradInvXX = bondLevelDeformationGradientInvXX;
  ScalarT* defGradInvXY = bondLevelDeformationGradientInvXY;
  ScalarT* defGradInvXZ = bondLevelDeformationGradientInvXZ;
  ScalarT* defGradInvYX = bondLevelDeformationGradientInvYX;
  ScalarT* defGradInvYY = bondLevelDeformationGradientInvYY;
  ScalarT* defGradInvYZ = bondLevelDeformationGradientInvYZ;
  ScalarT* defGradInvZX = bondLevelDeformationGradientInvZX;
  ScalarT* defGradInvZY = bondLevelDeformationGradientInvZY;
  ScalarT* defGradInvZZ = bondLevelDeformationGradientInvZZ;
  const ScalarT* leftStretchXXN = bondLevelLeftStretchTensorXXN;
  const ScalarT* leftStretchXYN = bondLevelLeftStretchTensorXYN;
  const ScalarT* leftStretchXZN = bondLevelLeftStretchTensorXZN;
  const ScalarT* leftStretchYXN = bondLevelLeftStretchTensorYXN;
  const ScalarT* leftStretchYYN = bondLevelLeftStretchTensorYYN;
  const ScalarT* leftStretchYZN = bondLevelLeftStretchTensorYZN;
  const ScalarT* leftStretchZXN = bondLevelLeftStretchTensorZXN;
  const ScalarT* leftStretchZYN = bondLevelLeftStretchTensorZYN;
  const ScalarT* leftStretchZZN = bondLevelLeftStretchTensorZZN;
  const ScalarT* rotTensorXXN = bondLevelRotationTensorXXN;
  const ScalarT* rotTensorXYN = bondLevelRotationTensorXYN;
  const ScalarT* rotTensorXZN = bondLevelRotationTensorXZN;
  const ScalarT* rotTensorYXN = bondLevelRotationTensorYXN;
  const ScalarT* rotTensorYYN = bondLevelRotationTensorYYN;
  const ScalarT* rotTensorYZN = bondLevelRotationTensorYZN;
  const ScalarT* rotTensorZXN = bondLevelRotationTensorZXN;
  const ScalarT* rotTensorZYN = bondLevelRotationTensorZYN;
  const ScalarT* rotTensorZZN = bondLevelRotationTensorZZN;
  const double *omega = influenceState;

  ScalarT* leftStretchXXNP1 = bondLevelLeftStretchTensorXXNP1;
  ScalarT* leftStretchXYNP1 = bondLevelLeftStretchTensorXYNP1;
  ScalarT* leftStretchXZNP1 = bondLevelLeftStretchTensorXZNP1;
  ScalarT* leftStretchYXNP1 = bondLevelLeftStretchTensorYXNP1;
  ScalarT* leftStretchYYNP1 = bondLevelLeftStretchTensorYYNP1;
  ScalarT* leftStretchYZNP1 = bondLevelLeftStretchTensorYZNP1;
  ScalarT* leftStretchZXNP1 = bondLevelLeftStretchTensorZXNP1;
  ScalarT* leftStretchZYNP1 = bondLevelLeftStretchTensorZYNP1;
  ScalarT* leftStretchZZNP1 = bondLevelLeftStretchTensorZZNP1;
  ScalarT* rotTensorXXNP1 = bondLevelRotationTensorXXNP1;
  ScalarT* rotTensorXYNP1 = bondLevelRotationTensorXYNP1;
  ScalarT* rotTensorXZNP1 = bondLevelRotationTensorXZNP1;
  ScalarT* rotTensorYXNP1 = bondLevelRotationTensorYXNP1;
  ScalarT* rotTensorYYNP1 = bondLevelRotationTensorYYNP1;
  ScalarT* rotTensorYZNP1 = bondLevelRotationTensorYZNP1;
  ScalarT* rotTensorZXNP1 = bondLevelRotationTensorZXNP1;
  ScalarT* rotTensorZYNP1 = bondLevelRotationTensorZYNP1;
  ScalarT* rotTensorZZNP1 = bondLevelRotationTensorZZNP1;
  ScalarT* unrotRateOfDefXX = bondLevelUnrotatedRateOfDeformationXX;
  ScalarT* unrotRateOfDefXY = bondLevelUnrotatedRateOfDeformationXY;
  ScalarT* unrotRateOfDefXZ = bondLevelUnrotatedRateOfDeformationXZ;
  ScalarT* unrotRateOfDefYX = bondLevelUnrotatedRateOfDeformationYX;
  ScalarT* unrotRateOfDefYY = bondLevelUnrotatedRateOfDeformationYY;
  ScalarT* unrotRateOfDefYZ = bondLevelUnrotatedRateOfDeformationYZ;
  ScalarT* unrotRateOfDefZX = bondLevelUnrotatedRateOfDeformationZX;
  ScalarT* unrotRateOfDefZY = bondLevelUnrotatedRateOfDeformationZY;
  ScalarT* unrotRateOfDefZZ = bondLevelUnrotatedRateOfDeformationZZ;

  std::vector<ScalarT> meanDefGradVector(9); ScalarT* meanDefGrad = &meanDefGradVector[0];
  std::vector<ScalarT> meanDefGradDotVector(9); ScalarT* meanDefGradDot = &meanDefGradDotVector[0];
  std::vector<ScalarT> defGradVector(9); ScalarT* defGrad = &defGradVector[0];
  std::vector<ScalarT> defGradDotVector(9); ScalarT* defGradDot = &defGradDotVector[0];
  std::vector<ScalarT> FinverseVector(9) ; ScalarT* Finverse = &FinverseVector[0];

  std::vector<ScalarT> eulerianVelGradVector(9) ; ScalarT* eulerianVelGrad = &eulerianVelGradVector[0];
  std::vector<ScalarT> leftStretchNVector(9) ; ScalarT* leftStretchN = &leftStretchNVector[0];
  std::vector<ScalarT> leftStretchNP1Vector(9) ; ScalarT* leftStretchNP1 = &leftStretchNP1Vector[0];
  std::vector<ScalarT> rotTensorNVector(9) ; ScalarT* rotTensorN = &rotTensorNVector[0];
  std::vector<ScalarT> rotTensorNP1Vector(9) ; ScalarT* rotTensorNP1 = &rotTensorNP1Vector[0];
  std::vector<ScalarT> unrotRateOfDefVector(9) ; ScalarT* unrotRateOfDef = &unrotRateOfDefVector[0];

  std::vector<ScalarT> rateOfDefVector(9) ; ScalarT* rateOfDef = &rateOfDefVector[0];
  std::vector<ScalarT> spinVector(9) ; ScalarT* spin = &spinVector[0];
  std::vector<ScalarT> tempVector(9) ; ScalarT* temp = &tempVector[0];
  std::vector<ScalarT> tempInvVector(9) ; ScalarT* tempInv = &tempInvVector[0];
  std::vector<ScalarT> OmegaTensorVector(9) ; ScalarT* OmegaTensor = &OmegaTensorVector[0];
  std::vector<ScalarT> QMatrixVector(9) ; ScalarT* QMatrix = &QMatrixVector[0];
  std::vector<ScalarT> OmegaTensorSqVector(9) ; ScalarT* OmegaTensorSq = &OmegaTensorSqVector[0];
  std::vector<ScalarT> tempAVector(9) ; ScalarT* tempA = &tempAVector[0];
  std::vector<ScalarT> tempBVector(9) ; ScalarT* tempB = &tempBVector[0];
  std::vector<ScalarT> rateOfStretchVector(9) ; ScalarT* rateOfStretch = &rateOfStretchVector[0];

  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLengthSq;
  ScalarT defStateX, defStateY, defStateZ;
  ScalarT velStateX, velStateY, velStateZ;

  ScalarT scalarTemp;
  ScalarT determinant;
  ScalarT omegaX, omegaY, omegaZ;
  ScalarT zX, zY, zZ;
  ScalarT wX, wY, wZ;
  ScalarT traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;
  int inversionReturnCode(0);
  std::string inversionErrorMessage = 
    "**** Error:  CORRESPONDENCE::computeBondLevelUnrotatedRateOfDeformationAndRotationTensor: Non-invertible matrix\n";

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, modelCoord+=3, coord+=3, defGradX+=3,
      defGradY+=3, defGradZ+=3, vel+=3, defGradDotX+=3, defGradDotY+=3, defGradDotZ+=3){

    // All is bond level.
    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, omega++, J++,
        defGradInvXX++, defGradInvXY++, defGradInvXZ++, 
        defGradInvYX++, defGradInvYY++, defGradInvYZ++, 
        defGradInvZX++, defGradInvZY++, defGradInvZZ++,
        leftStretchXXN++, leftStretchXYN++, leftStretchXZN++, 
        leftStretchYXN++, leftStretchYYN++, leftStretchYZN++, 
        leftStretchZXN++, leftStretchZYN++, leftStretchZZN++,
        rotTensorXXN++, rotTensorXYN++, rotTensorXZN++,
        rotTensorYXN++, rotTensorYYN++, rotTensorYZN++,
        rotTensorZXN++, rotTensorZYN++, rotTensorZZN++,
        leftStretchXXNP1++, leftStretchXYNP1++, leftStretchXZNP1++, 
        leftStretchYXNP1++, leftStretchYYNP1++, leftStretchYZNP1++, 
        leftStretchZXNP1++, leftStretchZYNP1++, leftStretchZZNP1++,
        rotTensorXXNP1++, rotTensorXYNP1++, rotTensorXZNP1++,
        rotTensorYXNP1++, rotTensorYYNP1++, rotTensorYZNP1++,
        rotTensorZXNP1++, rotTensorZYNP1++, rotTensorZZNP1++,
        unrotRateOfDefXX++, unrotRateOfDefXY++, unrotRateOfDefXZ++,
        unrotRateOfDefYX++, unrotRateOfDefYY++, unrotRateOfDefYZ++,
        unrotRateOfDefZX++, unrotRateOfDefZY++, unrotRateOfDefZZ++){

      if(*omega > 0.0){
        neighborIndex = *neighborListPtr;

        neighborModelCoord = modelCoordinates + 3*neighborIndex;
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVel = velocities + 3*neighborIndex;
        neighborDefGradX = deformationGradientX + 3*neighborIndex;
        neighborDefGradY = deformationGradientY + 3*neighborIndex;
        neighborDefGradZ = deformationGradientZ + 3*neighborIndex;
        neighborDefGradDotX = deformationGradientDotX + 3*neighborIndex;
        neighborDefGradDotY = deformationGradientDotY + 3*neighborIndex;
        neighborDefGradDotZ = deformationGradientDotZ + 3*neighborIndex;

        undeformedBondX = *(neighborModelCoord+0) - *(modelCoord+0);
        undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
        undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
        undeformedBondLengthSq = undeformedBondX*undeformedBondX +
                                 undeformedBondY*undeformedBondY +
                                 undeformedBondZ*undeformedBondZ;

        defStateX = *(neighborCoord+0) - *(coord+0);
        defStateY = *(neighborCoord+1) - *(coord+1);
        defStateZ = *(neighborCoord+2) - *(coord+2);

        velStateX = *(neighborVel+0) - *(vel+0);
        velStateY = *(neighborVel+1) - *(vel+1);
        velStateZ = *(neighborVel+2) - *(vel+2);

        // average of the two points 
        for(int i=0; i<3; i++){
          *(meanDefGrad+i+0) = 0.5 * (*(defGradX+i) + *(neighborDefGradX+i));
          *(meanDefGrad+i+3) = 0.5 * (*(defGradY+i) + *(neighborDefGradY+i));
          *(meanDefGrad+i+6) = 0.5 * (*(defGradZ+i) + *(neighborDefGradZ+i));
          *(meanDefGradDot+i+0) = 0.5 * (*(defGradDotX+i) + *(neighborDefGradDotX+i));
          *(meanDefGradDot+i+3) = 0.5 * (*(defGradDotY+i) + *(neighborDefGradDotY+i));
          *(meanDefGradDot+i+6) = 0.5 * (*(defGradDotZ+i) + *(neighborDefGradDotZ+i));
        }

        scalarTemp = *(meanDefGrad+0) * undeformedBondX + *(meanDefGrad+1) * undeformedBondY + *(meanDefGrad+2) * undeformedBondZ;
        *(defGrad+0) = *(meanDefGrad+0) + (defStateX - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
        *(defGrad+1) = *(meanDefGrad+1) + (defStateX - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
        *(defGrad+2) = *(meanDefGrad+2) + (defStateX - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

        scalarTemp = *(meanDefGrad+3) * undeformedBondX + *(meanDefGrad+4) * undeformedBondY + *(meanDefGrad+5) * undeformedBondZ;
        *(defGrad+3) = *(meanDefGrad+3) + (defStateY - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
        *(defGrad+4) = *(meanDefGrad+4) + (defStateY - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
        *(defGrad+5) = *(meanDefGrad+5) + (defStateY - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

        scalarTemp = *(meanDefGrad+6) * undeformedBondX + *(meanDefGrad+7) * undeformedBondY + *(meanDefGrad+8) * undeformedBondZ;
        *(defGrad+6) = *(meanDefGrad+6) + (defStateZ - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
        *(defGrad+7) = *(meanDefGrad+7) + (defStateZ - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
        *(defGrad+8) = *(meanDefGrad+8) + (defStateZ - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

        scalarTemp = *(meanDefGradDot+0) * undeformedBondX + *(meanDefGradDot+1) * undeformedBondY + *(meanDefGradDot+2) * undeformedBondZ;
        *(defGradDot+0) = *(meanDefGradDot+0) + (velStateX - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
        *(defGradDot+1) = *(meanDefGradDot+1) + (velStateX - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
        *(defGradDot+2) = *(meanDefGradDot+2) + (velStateX - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

        scalarTemp = *(meanDefGradDot+3) * undeformedBondX + *(meanDefGradDot+4) * undeformedBondY + *(meanDefGradDot+5) * undeformedBondZ;
        *(defGradDot+3) = *(meanDefGradDot+3) + (velStateY - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
        *(defGradDot+4) = *(meanDefGradDot+4) + (velStateY - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
        *(defGradDot+5) = *(meanDefGradDot+5) + (velStateY - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

        scalarTemp = *(meanDefGradDot+6) * undeformedBondX + *(meanDefGradDot+7) * undeformedBondY + *(meanDefGradDot+8) * undeformedBondZ;
        *(defGradDot+6) = *(meanDefGradDot+6) + (velStateZ - scalarTemp) * undeformedBondX/undeformedBondLengthSq;
        *(defGradDot+7) = *(meanDefGradDot+7) + (velStateZ - scalarTemp) * undeformedBondY/undeformedBondLengthSq;
        *(defGradDot+8) = *(meanDefGradDot+8) + (velStateZ - scalarTemp) * undeformedBondZ/undeformedBondLengthSq;

        // Compute the inverse of the deformation gradient, Finverse
        inversionReturnCode = MATRICES::Invert3by3Matrix(defGrad, determinant, Finverse);
        if(inversionReturnCode > 0){
          returnCode = inversionReturnCode;
          std::cout << inversionErrorMessage;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        *J = determinant;
        *(defGradInvXX) = *(Finverse+0); *(defGradInvXY) = *(Finverse+1); *(defGradInvXZ) = *(Finverse+2); 
        *(defGradInvYX) = *(Finverse+3); *(defGradInvYY) = *(Finverse+4); *(defGradInvYZ) = *(Finverse+5); 
        *(defGradInvZX) = *(Finverse+6); *(defGradInvZY) = *(Finverse+7); *(defGradInvZZ) = *(Finverse+8); 

        // Compute the Eulerian velocity gradient L = Fdot * Finv
        MATRICES::MatrixMultiply(false, false, 1.0, defGradDot, Finverse, eulerianVelGrad);

        // Store in a tensor form 
        *(leftStretchN+0) = *leftStretchXXN; *(leftStretchN+1) = *leftStretchXYN; *(leftStretchN+2) = *leftStretchXZN;
        *(leftStretchN+3) = *leftStretchYXN; *(leftStretchN+4) = *leftStretchYYN; *(leftStretchN+5) = *leftStretchYZN;
        *(leftStretchN+6) = *leftStretchZXN; *(leftStretchN+7) = *leftStretchZYN; *(leftStretchN+8) = *leftStretchZZN;
        *(rotTensorN+0) = *rotTensorXXN; *(rotTensorN+1) = *rotTensorXYN; *(rotTensorN+2) = *rotTensorXZN;
        *(rotTensorN+3) = *rotTensorYXN; *(rotTensorN+4) = *rotTensorYYN; *(rotTensorN+5) = *rotTensorYZN;
        *(rotTensorN+6) = *rotTensorZXN; *(rotTensorN+7) = *rotTensorZYN; *(rotTensorN+8) = *rotTensorZZN;

        // Compute rate-of-deformation tensor, D = 1/2 * (L + Lt)
        *(rateOfDef)   = *(eulerianVelGrad);
        *(rateOfDef+1) = 0.5 * ( *(eulerianVelGrad+1) + *(eulerianVelGrad+3) );
        *(rateOfDef+2) = 0.5 * ( *(eulerianVelGrad+2) + *(eulerianVelGrad+6) );
        *(rateOfDef+3) = *(rateOfDef+1);
        *(rateOfDef+4) = *(eulerianVelGrad+4);
        *(rateOfDef+5) = 0.5 * ( *(eulerianVelGrad+5) + *(eulerianVelGrad+7) );
        *(rateOfDef+6) = *(rateOfDef+2);
        *(rateOfDef+7) = *(rateOfDef+5);
        *(rateOfDef+8) = *(eulerianVelGrad+8);

        // Compute spin tensor, W = 1/2 * (L - Lt)
        *(spin)   = 0.0;
        *(spin+1) = 0.5 * ( *(eulerianVelGrad+1) - *(eulerianVelGrad+3) );
        *(spin+2) = 0.5 * ( *(eulerianVelGrad+2) - *(eulerianVelGrad+6) );
        *(spin+3) = -1.0 * *(spin+1);
        *(spin+4) = 0.0;
        *(spin+5) = 0.5 * ( *(eulerianVelGrad+5) - *(eulerianVelGrad+7) );
        *(spin+6) = -1.0 * *(spin+2);
        *(spin+7) = -1.0 * *(spin+5);
        *(spin+8) = 0.0;
       
        //Following Flanagan & Taylor (T&F) 
        //
        //Find the vector z_i = \epsilon_{ikj} * D_{jm} * V_{mk} (T&F Eq. 13)
        //
        //where \epsilon_{ikj} is the alternator tensor.
        //
        //Components below copied from computer algebra solution to the expansion
        //above
        zX = - *(leftStretchN+2) *  *(rateOfDef+3) -  *(leftStretchN+5) *  *(rateOfDef+4) - 
               *(leftStretchN+8) *  *(rateOfDef+5) +  *(leftStretchN+1) *  *(rateOfDef+6) + 
               *(leftStretchN+4) *  *(rateOfDef+7) +  *(leftStretchN+7) *  *(rateOfDef+8);
        zY =   *(leftStretchN+2) *  *(rateOfDef)   +  *(leftStretchN+5) *  *(rateOfDef+1) + 
               *(leftStretchN+8) *  *(rateOfDef+2) -  *(leftStretchN)   *  *(rateOfDef+6) - 
               *(leftStretchN+3) *  *(rateOfDef+7) -  *(leftStretchN+6) *  *(rateOfDef+8);
        zZ = - *(leftStretchN+1) *  *(rateOfDef)   -  *(leftStretchN+4) *  *(rateOfDef+1) - 
               *(leftStretchN+7) *  *(rateOfDef+2) +  *(leftStretchN)   *  *(rateOfDef+3) + 
               *(leftStretchN+3) *  *(rateOfDef+4) +  *(leftStretchN+6) *  *(rateOfDef+5);

        //Find the vector w_i = -1/2 * \epsilon_{ijk} * W_{jk} (T&F Eq. 11)
        wX = 0.5 * ( *(spin+7) - *(spin+5) );
        wY = 0.5 * ( *(spin+2) - *(spin+6) );
        wZ = 0.5 * ( *(spin+3) - *(spin+1) );

        //Find trace(V)
        traceV = *(leftStretchN) + *(leftStretchN+4) + *(leftStretchN+8);

        // Compute (trace(V) * I - V) store in temp
        *(temp)   = traceV - *(leftStretchN);
        *(temp+1) = - *(leftStretchN+1);
        *(temp+2) = - *(leftStretchN+2);
        *(temp+3) = - *(leftStretchN+3);
        *(temp+4) = traceV - *(leftStretchN+4);
        *(temp+5) = - *(leftStretchN+5);
        *(temp+6) = - *(leftStretchN+6);
        *(temp+7) = - *(leftStretchN+7);
        *(temp+8) = traceV - *(leftStretchN+8);

        // Compute the inverse of the temp matrix
        MATRICES::Invert3by3Matrix(temp, determinant, tempInv);
        if(inversionReturnCode > 0){
          returnCode = inversionReturnCode;
          std::cout << inversionErrorMessage;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //Find omega vector, i.e. \omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
        omegaX =  wX + *(tempInv)   * zX + *(tempInv+1) * zY + *(tempInv+2) * zZ;
        omegaY =  wY + *(tempInv+3) * zX + *(tempInv+4) * zY + *(tempInv+5) * zZ;
        omegaZ =  wZ + *(tempInv+6) * zX + *(tempInv+7) * zY + *(tempInv+8) * zZ;

        //Find the tensor \Omega_{ij} = \epsilon_{ikj} * w_k (T&F Eq. 10)
        *(OmegaTensor) = 0.0;
        *(OmegaTensor+1) = -omegaZ;
        *(OmegaTensor+2) = omegaY;
        *(OmegaTensor+3) = omegaZ;
        *(OmegaTensor+4) = 0.0;
        *(OmegaTensor+5) = -omegaX;
        *(OmegaTensor+6) = -omegaY;
        *(OmegaTensor+7) = omegaX;
        *(OmegaTensor+8) = 0.0;

        //Increment R with (T&F Eq. 36 and 44) as opposed to solving (T&F 39) this
        //is desirable for accuracy in implicit solves and has no effect on
        //explicit solves (other than a slight decrease in speed).
        //
        // Compute Q with (T&F Eq. 44)
        //
        // Omega^2 = w_i * w_i (T&F Eq. 42)
        OmegaSq = omegaX*omegaX + omegaY*omegaY + omegaZ*omegaZ;
        // Omega = \sqrt{OmegaSq}
        Omega = sqrt(OmegaSq);

        // Avoid a potential divide-by-zero
        if(OmegaSq > 1.e-30){

          // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
          //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
          scaleFactor1 = sin(dt*Omega) / Omega;
          scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
          MATRICES::MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, OmegaTensorSq);
          *(QMatrix)   = 1.0 + scaleFactor1 * *(OmegaTensor)   + scaleFactor2 * *(OmegaTensorSq)   ;
          *(QMatrix+1) =       scaleFactor1 * *(OmegaTensor+1) + scaleFactor2 * *(OmegaTensorSq+1) ;
          *(QMatrix+2) =       scaleFactor1 * *(OmegaTensor+2) + scaleFactor2 * *(OmegaTensorSq+2) ;
          *(QMatrix+3) =       scaleFactor1 * *(OmegaTensor+3) + scaleFactor2 * *(OmegaTensorSq+3) ;
          *(QMatrix+4) = 1.0 + scaleFactor1 * *(OmegaTensor+4) + scaleFactor2 * *(OmegaTensorSq+4) ;
          *(QMatrix+5) =       scaleFactor1 * *(OmegaTensor+5) + scaleFactor2 * *(OmegaTensorSq+5) ;
          *(QMatrix+6) =       scaleFactor1 * *(OmegaTensor+6) + scaleFactor2 * *(OmegaTensorSq+6) ;
          *(QMatrix+7) =       scaleFactor1 * *(OmegaTensor+7) + scaleFactor2 * *(OmegaTensorSq+7) ;
          *(QMatrix+8) = 1.0 + scaleFactor1 * *(OmegaTensor+8) + scaleFactor2 * *(OmegaTensorSq+8) ;

        } else {
          *(QMatrix)   = 1.0 ; *(QMatrix+1) = 0.0 ; *(QMatrix+2) = 0.0 ;
          *(QMatrix+3) = 0.0 ; *(QMatrix+4) = 1.0 ; *(QMatrix+5) = 0.0 ;
          *(QMatrix+6) = 0.0 ; *(QMatrix+7) = 0.0 ; *(QMatrix+8) = 1.0 ;
        };

        // Compute R_STEP_NP1 = QMatrix * R_STEP_N (T&F Eq. 36)
        MATRICES::MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, rotTensorNP1);

        // Compute rate of stretch, Vdot = L*V - V*Omega
        // First tempA = L*V, 
        MATRICES::MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, tempA);

        // tempB = V*Omega
        MATRICES::MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, tempB);

        //Vdot = tempA - tempB
        for(int i=0 ; i<9 ; ++i)
          *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

        //V_STEP_NP1 = V_STEP_N + dt*Vdot
        for(int i=0 ; i<9 ; ++i)
          *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

        // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
        MATRICES::MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, temp);

        // d = Rt * temp
        MATRICES::MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, unrotRateOfDef);

        // Store back in element-wise format
        *leftStretchXXNP1 = *(leftStretchNP1+0); *leftStretchXYNP1 = *(leftStretchNP1+1); *leftStretchXZNP1 = *(leftStretchNP1+2);
        *leftStretchYXNP1 = *(leftStretchNP1+3); *leftStretchYYNP1 = *(leftStretchNP1+4); *leftStretchYZNP1 = *(leftStretchNP1+5);
        *leftStretchZXNP1 = *(leftStretchNP1+6); *leftStretchZYNP1 = *(leftStretchNP1+7); *leftStretchZZNP1 = *(leftStretchNP1+8);
        *rotTensorXXNP1 = *(rotTensorNP1+0); *rotTensorXYNP1 = *(rotTensorNP1+1); *rotTensorXZNP1 = *(rotTensorNP1+2);
        *rotTensorYXNP1 = *(rotTensorNP1+3); *rotTensorYYNP1 = *(rotTensorNP1+4); *rotTensorYZNP1 = *(rotTensorNP1+5);
        *rotTensorZXNP1 = *(rotTensorNP1+6); *rotTensorZYNP1 = *(rotTensorNP1+7); *rotTensorZZNP1 = *(rotTensorNP1+8);
        *unrotRateOfDefXX = *(unrotRateOfDef+0); *unrotRateOfDefXY = *(unrotRateOfDef+1); *unrotRateOfDefXZ = *(unrotRateOfDef+2);
        *unrotRateOfDefYX = *(unrotRateOfDef+3); *unrotRateOfDefYY = *(unrotRateOfDef+4); *unrotRateOfDefYZ = *(unrotRateOfDef+5);
        *unrotRateOfDefZX = *(unrotRateOfDef+6); *unrotRateOfDefZY = *(unrotRateOfDef+7); *unrotRateOfDefZZ = *(unrotRateOfDef+8);
      }
      else{
        *leftStretchXXNP1 = *leftStretchXXN; *leftStretchXYNP1 = *leftStretchXYN; *leftStretchXZNP1 = *leftStretchXZN;
        *leftStretchYXNP1 = *leftStretchYXN; *leftStretchYYNP1 = *leftStretchYYN; *leftStretchYZNP1 = *leftStretchYZN;
        *leftStretchZXNP1 = *leftStretchZXN; *leftStretchZYNP1 = *leftStretchZYN; *leftStretchZZNP1 = *leftStretchZZN;
        *rotTensorXXNP1 = *rotTensorXXN; *rotTensorXYNP1 = *rotTensorXYN; *rotTensorXZNP1 = *rotTensorXZN;
        *rotTensorYXNP1 = *rotTensorYXN; *rotTensorYYNP1 = *rotTensorYYN; *rotTensorYZNP1 = *rotTensorYZN;
        *rotTensorZXNP1 = *rotTensorZXN; *rotTensorZYNP1 = *rotTensorZYN; *rotTensorZZNP1 = *rotTensorZZN;
        *unrotRateOfDefXX = 0.0; *unrotRateOfDefXY = 0.0; *unrotRateOfDefXZ = 0.0;
        *unrotRateOfDefYX = 0.0; *unrotRateOfDefYY = 0.0; *unrotRateOfDefYZ = 0.0;
        *unrotRateOfDefZX = 0.0; *unrotRateOfDefZY = 0.0; *unrotRateOfDefZZ = 0.0;
      }
    }
  }
  return returnCode;
}




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
)
{
  const ScalarT* rotTensorXX = bondLevelRotationTensorXX;
  const ScalarT* rotTensorXY = bondLevelRotationTensorXY;
  const ScalarT* rotTensorXZ = bondLevelRotationTensorXZ;
  const ScalarT* rotTensorYX = bondLevelRotationTensorYX;
  const ScalarT* rotTensorYY = bondLevelRotationTensorYY;
  const ScalarT* rotTensorYZ = bondLevelRotationTensorYZ;
  const ScalarT* rotTensorZX = bondLevelRotationTensorZX;
  const ScalarT* rotTensorZY = bondLevelRotationTensorZY;
  const ScalarT* rotTensorZZ = bondLevelRotationTensorZZ;
  const ScalarT* unrotatedStressXX = bondLevelUnrotatedCauchyStressXX;
  const ScalarT* unrotatedStressXY = bondLevelUnrotatedCauchyStressXY;
  const ScalarT* unrotatedStressXZ = bondLevelUnrotatedCauchyStressXZ;
  const ScalarT* unrotatedStressYX = bondLevelUnrotatedCauchyStressYX;
  const ScalarT* unrotatedStressYY = bondLevelUnrotatedCauchyStressYY;
  const ScalarT* unrotatedStressYZ = bondLevelUnrotatedCauchyStressYZ;
  const ScalarT* unrotatedStressZX = bondLevelUnrotatedCauchyStressZX;
  const ScalarT* unrotatedStressZY = bondLevelUnrotatedCauchyStressZY;
  const ScalarT* unrotatedStressZZ = bondLevelUnrotatedCauchyStressZZ;
  ScalarT* rotatedStressXX = bondLevelRotatedCauchyStressXX;
  ScalarT* rotatedStressXY = bondLevelRotatedCauchyStressXY;
  ScalarT* rotatedStressXZ = bondLevelRotatedCauchyStressXZ;
  ScalarT* rotatedStressYX = bondLevelRotatedCauchyStressYX;
  ScalarT* rotatedStressYY = bondLevelRotatedCauchyStressYY;
  ScalarT* rotatedStressYZ = bondLevelRotatedCauchyStressYZ;
  ScalarT* rotatedStressZX = bondLevelRotatedCauchyStressZX;
  ScalarT* rotatedStressZY = bondLevelRotatedCauchyStressZY;
  ScalarT* rotatedStressZZ = bondLevelRotatedCauchyStressZZ;

  ScalarT unrotatedStress[9], rotTensor[9], rotatedStress[9], temp[9];

  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID){

    // All is bond level.
    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++,
        rotTensorXX++, rotTensorXY++, rotTensorXZ++, 
        rotTensorYX++, rotTensorYY++, rotTensorYZ++, 
        rotTensorZX++, rotTensorZY++, rotTensorZZ++, 
        unrotatedStressXX++, unrotatedStressXY++, unrotatedStressXZ++, 
        unrotatedStressYX++, unrotatedStressYY++, unrotatedStressYZ++, 
        unrotatedStressZX++, unrotatedStressZY++, unrotatedStressZZ++, 
        rotatedStressXX++, rotatedStressXY++, rotatedStressXZ++, 
        rotatedStressYX++, rotatedStressYY++, rotatedStressYZ++, 
        rotatedStressZX++, rotatedStressZY++, rotatedStressZZ++){

      // write in matrix form 
      rotTensor[0] = *rotTensorXX; rotTensor[1] = *rotTensorXY; rotTensor[2] = *rotTensorXZ;
      rotTensor[3] = *rotTensorYX; rotTensor[4] = *rotTensorYY; rotTensor[5] = *rotTensorYZ;
      rotTensor[6] = *rotTensorZX; rotTensor[7] = *rotTensorZY; rotTensor[8] = *rotTensorZZ;
      unrotatedStress[0] = *unrotatedStressXX; unrotatedStress[1] = *unrotatedStressXY; unrotatedStress[2] = *unrotatedStressXZ;
      unrotatedStress[3] = *unrotatedStressYX; unrotatedStress[4] = *unrotatedStressYY; unrotatedStress[5] = *unrotatedStressYZ;
      unrotatedStress[6] = *unrotatedStressZX; unrotatedStress[7] = *unrotatedStressZY; unrotatedStress[8] = *unrotatedStressZZ;

      // temp = \sigma_unrot * Rt
      MATRICES::MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, temp);
      // \sigma_rot = R * temp
      MATRICES::MatrixMultiply(false, false, 1.0, rotTensor, temp, rotatedStress);

      // update bond-level data field
      *rotatedStressXX = rotatedStress[0]; *rotatedStressXY = rotatedStress[1]; *rotatedStressXZ = rotatedStress[2]; 
      *rotatedStressYX = rotatedStress[3]; *rotatedStressYY = rotatedStress[4]; *rotatedStressYZ = rotatedStress[5]; 
      *rotatedStressZX = rotatedStress[6]; *rotatedStressZY = rotatedStress[7]; *rotatedStressZZ = rotatedStress[8]; 
    }
  }
}


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
)
{
  const ScalarT* J = bondLevelJacobianDeterminant;
  const ScalarT* cauchyStressXX = bondLevelCauchyStressXX;
  const ScalarT* cauchyStressXY = bondLevelCauchyStressXY;
  const ScalarT* cauchyStressXZ = bondLevelCauchyStressXZ;
  const ScalarT* cauchyStressYX = bondLevelCauchyStressYX;
  const ScalarT* cauchyStressYY = bondLevelCauchyStressYY;
  const ScalarT* cauchyStressYZ = bondLevelCauchyStressYZ;
  const ScalarT* cauchyStressZX = bondLevelCauchyStressZX;
  const ScalarT* cauchyStressZY = bondLevelCauchyStressZY;
  const ScalarT* cauchyStressZZ = bondLevelCauchyStressZZ;
  const ScalarT* defGradInvXX = bondLevelDeformationGradientInvXX;
  const ScalarT* defGradInvXY = bondLevelDeformationGradientInvXY;
  const ScalarT* defGradInvXZ = bondLevelDeformationGradientInvXZ;
  const ScalarT* defGradInvYX = bondLevelDeformationGradientInvYX;
  const ScalarT* defGradInvYY = bondLevelDeformationGradientInvYY;
  const ScalarT* defGradInvYZ = bondLevelDeformationGradientInvYZ;
  const ScalarT* defGradInvZX = bondLevelDeformationGradientInvZX;
  const ScalarT* defGradInvZY = bondLevelDeformationGradientInvZY;
  const ScalarT* defGradInvZZ = bondLevelDeformationGradientInvZZ;
  ScalarT* piolaStressXX = bondLevelPiolaStressXX;
  ScalarT* piolaStressXY = bondLevelPiolaStressXY;
  ScalarT* piolaStressXZ = bondLevelPiolaStressXZ;
  ScalarT* piolaStressYX = bondLevelPiolaStressYX;
  ScalarT* piolaStressYY = bondLevelPiolaStressYY;
  ScalarT* piolaStressYZ = bondLevelPiolaStressYZ;
  ScalarT* piolaStressZX = bondLevelPiolaStressZX;
  ScalarT* piolaStressZY = bondLevelPiolaStressZY;
  ScalarT* piolaStressZZ = bondLevelPiolaStressZZ;

  std::vector<ScalarT> cauchyStressVector(9);
  ScalarT* cauchyStress = &cauchyStressVector[0];

  std::vector<ScalarT> defGradInvVector(9);
  ScalarT* defGradInv = &defGradInvVector[0];

  std::vector<ScalarT> piolaStressVector(9);
  ScalarT* piolaStress = &piolaStressVector[0];

  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID){

    // All is bond level.
    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, J++,
        cauchyStressXX++, cauchyStressXY++, cauchyStressXZ++, 
        cauchyStressYX++, cauchyStressYY++, cauchyStressYZ++, 
        cauchyStressZX++, cauchyStressZY++, cauchyStressZZ++, 
        defGradInvXX++, defGradInvXY++, defGradInvXZ++, 
        defGradInvYX++, defGradInvYY++, defGradInvYZ++, 
        defGradInvZX++, defGradInvZY++, defGradInvZZ++, 
        piolaStressXX++, piolaStressXY++, piolaStressXZ++, 
        piolaStressYX++, piolaStressYY++, piolaStressYZ++, 
        piolaStressZX++, piolaStressZY++, piolaStressZZ++){

      cauchyStress[0] = *cauchyStressXX;
      cauchyStress[1] = *cauchyStressXY;
      cauchyStress[2] = *cauchyStressXZ;
      cauchyStress[3] = *cauchyStressYX;
      cauchyStress[4] = *cauchyStressYY;
      cauchyStress[5] = *cauchyStressYZ;
      cauchyStress[6] = *cauchyStressZX;
      cauchyStress[7] = *cauchyStressZY;
      cauchyStress[8] = *cauchyStressZZ;

      defGradInv[0] = *defGradInvXX;
      defGradInv[1] = *defGradInvXY;
      defGradInv[2] = *defGradInvXZ;
      defGradInv[3] = *defGradInvYX;
      defGradInv[4] = *defGradInvYY;
      defGradInv[5] = *defGradInvYZ;
      defGradInv[6] = *defGradInvZX;
      defGradInv[7] = *defGradInvZY;
      defGradInv[8] = *defGradInvZZ;

      //P = J * \sigma * F^(-T)
      MATRICES::MatrixMultiply(false, true, *J, cauchyStress, defGradInv, piolaStress);

      *piolaStressXX = piolaStress[0];
      *piolaStressXY = piolaStress[1];
      *piolaStressXZ = piolaStress[2];
      *piolaStressYX = piolaStress[3];
      *piolaStressYY = piolaStress[4];
      *piolaStressYZ = piolaStress[5];
      *piolaStressZX = piolaStress[6];
      *piolaStressZY = piolaStress[7];
      *piolaStressZZ = piolaStress[8];
    }
  }
}


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
)
{
  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const double* w0 = weightedVolume;
  const double* neighborW0;
  const ScalarT* stressXX = bondLevelStressXX;
  const ScalarT* stressXY = bondLevelStressXY;
  const ScalarT* stressXZ = bondLevelStressXZ;
  const ScalarT* stressYX = bondLevelStressYX;
  const ScalarT* stressYY = bondLevelStressYY;
  const ScalarT* stressYZ = bondLevelStressYZ;
  const ScalarT* stressZX = bondLevelStressZX;
  const ScalarT* stressZY = bondLevelStressZY;
  const ScalarT* stressZZ = bondLevelStressZZ;
  ScalarT* integral = stressIntegral;

  const double *omega = influenceState;

  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLengthSq;
  double neighborVolume, scalarTemp;

  std::vector<ScalarT> tempVector(9);
  ScalarT* temp = &tempVector[0];

  std::vector<ScalarT> stressVector(9);
  ScalarT* stress = &stressVector[0];

  std::vector<ScalarT> integrandVector(9);
  ScalarT* integrand = &integrandVector[0];

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, modelCoord+=3, w0++, integral+=9){

    // Zero out data
    *(integral)   = 0.0 ; *(integral+1) = 0.0 ; *(integral+2) = 0.0 ;
    *(integral+3) = 0.0 ; *(integral+4) = 0.0 ; *(integral+5) = 0.0 ;
    *(integral+6) = 0.0 ; *(integral+7) = 0.0 ; *(integral+8) = 0.0 ;

    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, omega++,
        stressXX++, stressXY++, stressXZ++, 
        stressYX++, stressYY++, stressYZ++, 
        stressZX++, stressZY++, stressZZ++){

      if(*omega > 0.0){

        neighborIndex = *neighborListPtr;
        neighborVolume = volume[neighborIndex];
        neighborModelCoord = modelCoordinates + 3*neighborIndex;
        neighborW0 = weightedVolume + neighborIndex;

        undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
        undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
        undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
        undeformedBondLengthSq = undeformedBondX*undeformedBondX +
                               undeformedBondY*undeformedBondY +
                               undeformedBondZ*undeformedBondZ;

        // write the stress in matrix form 
        stress[0] = *stressXX; stress[1] = *stressXY; stress[2] = *stressXZ; 
        stress[3] = *stressYX; stress[4] = *stressYY; stress[5] = *stressYZ; 
        stress[6] = *stressZX; stress[7] = *stressZY; stress[8] = *stressZZ; 

        // delta_jp - (y_j y_p)/|y|^2
        *(temp+0) = 1.0 - undeformedBondX * undeformedBondX / undeformedBondLengthSq;
        *(temp+1) = - undeformedBondX * undeformedBondY / undeformedBondLengthSq;
        *(temp+2) = - undeformedBondX * undeformedBondZ / undeformedBondLengthSq;
        *(temp+3) = *(temp+1);
        *(temp+4) = 1.0 - undeformedBondY * undeformedBondY / undeformedBondLengthSq;
        *(temp+5) = - undeformedBondY * undeformedBondZ / undeformedBondLengthSq;
        *(temp+6) = *(temp+2);
        *(temp+7) = *(temp+5);
        *(temp+8) = 1.0 - undeformedBondZ * undeformedBondZ / undeformedBondLengthSq;

        // Matrix multiply the stress and the second term to compute the integrand
        MATRICES::MatrixMultiply(false, false, 1.0, stress, temp, integrand);

        scalarTemp = *omega * (0.5 / *w0 + 0.5 / *neighborW0) * neighborVolume;

        for(int i=0; i<9; i++)
          *(integral+i) += scalarTemp * *(integrand+i);
      }
    }
  }
}

/** Explicit template instantiation for double. */

template void computeBondLevelVelocityGradient<double>
(
    const double* coordinates,
    const double* velocities,
    const double* velocityGradient,
    double* bondLevelVelocityGradientXX,
    double* bondLevelVelocityGradientXY,
    double* bondLevelVelocityGradientXZ,
    double* bondLevelVelocityGradientYX,
    double* bondLevelVelocityGradientYY,
    double* bondLevelVelocityGradientYZ,
    double* bondLevelVelocityGradientZX,
    double* bondLevelVelocityGradientZY,
    double* bondLevelVelocityGradientZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints
);

template void computeBondLevelVelocityGradient<double>
(
    const double* coordinates,
    const double* velocities,
    const double* velocityGradientX,
    const double* velocityGradientY,
    const double* velocityGradientZ,
    double* bondLevelVelocityGradientXX,
    double* bondLevelVelocityGradientXY,
    double* bondLevelVelocityGradientXZ,
    double* bondLevelVelocityGradientYX,
    double* bondLevelVelocityGradientYY,
    double* bondLevelVelocityGradientYZ,
    double* bondLevelVelocityGradientZX,
    double* bondLevelVelocityGradientZY,
    double* bondLevelVelocityGradientZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints
);


template int computeBondLevelUnrotatedRateOfDeformationAndRotationTensor<double>
(
    const double* bondLevelVelocityGradientXX, 
    const double* bondLevelVelocityGradientXY, 
    const double* bondLevelVelocityGradientXZ,
    const double* bondLevelVelocityGradientYX, 
    const double* bondLevelVelocityGradientYY, 
    const double* bondLevelVelocityGradientYZ, 
    const double* bondLevelVelocityGradientZX,
    const double* bondLevelVelocityGradientZY,
    const double* bondLevelVelocityGradientZZ,
    const double* bondLevelLeftStretchTensorXXN,
    const double* bondLevelLeftStretchTensorXYN,
    const double* bondLevelLeftStretchTensorXZN,
    const double* bondLevelLeftStretchTensorYXN,
    const double* bondLevelLeftStretchTensorYYN,
    const double* bondLevelLeftStretchTensorYZN,
    const double* bondLevelLeftStretchTensorZXN,
    const double* bondLevelLeftStretchTensorZYN,
    const double* bondLevelLeftStretchTensorZZN,
    const double* bondLevelRotationTensorXXN, 
    const double* bondLevelRotationTensorXYN, 
    const double* bondLevelRotationTensorXZN, 
    const double* bondLevelRotationTensorYXN, 
    const double* bondLevelRotationTensorYYN, 
    const double* bondLevelRotationTensorYZN, 
    const double* bondLevelRotationTensorZXN, 
    const double* bondLevelRotationTensorZYN, 
    const double* bondLevelRotationTensorZZN, 
    double* bondLevelLeftStretchTensorXXNP1,
    double* bondLevelLeftStretchTensorXYNP1,
    double* bondLevelLeftStretchTensorXZNP1,
    double* bondLevelLeftStretchTensorYXNP1,
    double* bondLevelLeftStretchTensorYYNP1,
    double* bondLevelLeftStretchTensorYZNP1,
    double* bondLevelLeftStretchTensorZXNP1,
    double* bondLevelLeftStretchTensorZYNP1,
    double* bondLevelLeftStretchTensorZZNP1,
    double* bondLevelRotationTensorXXNP1,
    double* bondLevelRotationTensorXYNP1,
    double* bondLevelRotationTensorXZNP1,
    double* bondLevelRotationTensorYXNP1,
    double* bondLevelRotationTensorYYNP1,
    double* bondLevelRotationTensorYZNP1,
    double* bondLevelRotationTensorZXNP1,
    double* bondLevelRotationTensorZYNP1,
    double* bondLevelRotationTensorZZNP1,
    double* bondLevelUnrotatedRateOfDeformationXX,
    double* bondLevelUnrotatedRateOfDeformationXY,
    double* bondLevelUnrotatedRateOfDeformationXZ,
    double* bondLevelUnrotatedRateOfDeformationYX,
    double* bondLevelUnrotatedRateOfDeformationYY,
    double* bondLevelUnrotatedRateOfDeformationYZ,
    double* bondLevelUnrotatedRateOfDeformationZX,
    double* bondLevelUnrotatedRateOfDeformationZY,
    double* bondLevelUnrotatedRateOfDeformationZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints,
    double dt
);


template void rotateBondLevelCauchyStress(
    const double* bondLevelRotationTensorXX,
    const double* bondLevelRotationTensorXY,
    const double* bondLevelRotationTensorXZ,
    const double* bondLevelRotationTensorYX,
    const double* bondLevelRotationTensorYY,
    const double* bondLevelRotationTensorYZ,
    const double* bondLevelRotationTensorZX,
    const double* bondLevelRotationTensorZY,
    const double* bondLevelRotationTensorZZ,
    const double* bondLevelUnrotatedCauchyStressXX,
    const double* bondLevelUnrotatedCauchyStressXY,
    const double* bondLevelUnrotatedCauchyStressXZ,
    const double* bondLevelUnrotatedCauchyStressYX,
    const double* bondLevelUnrotatedCauchyStressYY,
    const double* bondLevelUnrotatedCauchyStressYZ,
    const double* bondLevelUnrotatedCauchyStressZX,
    const double* bondLevelUnrotatedCauchyStressZY,
    const double* bondLevelUnrotatedCauchyStressZZ,
    double* bondLevelRotatedCauchyStressXX,
    double* bondLevelRotatedCauchyStressXY,
    double* bondLevelRotatedCauchyStressXZ,
    double* bondLevelRotatedCauchyStressYX,
    double* bondLevelRotatedCauchyStressYY,
    double* bondLevelRotatedCauchyStressYZ,
    double* bondLevelRotatedCauchyStressZX,
    double* bondLevelRotatedCauchyStressZY,
    double* bondLevelRotatedCauchyStressZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints
);

template void computeNonhomogeneityIntegral<double>
(
    const double* volume,
    const double* weightedVolume,
    const double* jacobianDeterminant,
    const double* horizon,
    const double* coordinates,
    const double* bondLevelCauchyStressXX,
    const double* bondLevelCauchyStressXY,
    const double* bondLevelCauchyStressXZ,
    const double* bondLevelCauchyStressYX,
    const double* bondLevelCauchyStressYY,
    const double* bondLevelCauchyStressYZ,
    const double* bondLevelCauchyStressZX,
    const double* bondLevelCauchyStressZY,
    const double* bondLevelCauchyStressZZ,
    double* nonhomogeneousIntegral,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints
);


template int computeBondLevelUnrotatedRateOfDeformationAndRotationTensor<double>
(
    const double* modelCoordinates,
    const double* coordinates,
    const double* velocities,
    const double* deformationGradientX,
    const double* deformationGradientY,
    const double* deformationGradientZ,
    const double* deformationGradientDotX,
    const double* deformationGradientDotY,
    const double* deformationGradientDotZ,
    double* bondLevelDeformationGradientInvXX, 
    double* bondLevelDeformationGradientInvXY, 
    double* bondLevelDeformationGradientInvXZ,
    double* bondLevelDeformationGradientInvYX, 
    double* bondLevelDeformationGradientInvYY, 
    double* bondLevelDeformationGradientInvYZ, 
    double* bondLevelDeformationGradientInvZX,
    double* bondLevelDeformationGradientInvZY,
    double* bondLevelDeformationGradientInvZZ,
    double* bondLevelJacobianDeterminant,
    const double* bondLevelLeftStretchTensorXXN,
    const double* bondLevelLeftStretchTensorXYN,
    const double* bondLevelLeftStretchTensorXZN,
    const double* bondLevelLeftStretchTensorYXN,
    const double* bondLevelLeftStretchTensorYYN,
    const double* bondLevelLeftStretchTensorYZN,
    const double* bondLevelLeftStretchTensorZXN,
    const double* bondLevelLeftStretchTensorZYN,
    const double* bondLevelLeftStretchTensorZZN,
    const double* bondLevelRotationTensorXXN, 
    const double* bondLevelRotationTensorXYN, 
    const double* bondLevelRotationTensorXZN, 
    const double* bondLevelRotationTensorYXN, 
    const double* bondLevelRotationTensorYYN, 
    const double* bondLevelRotationTensorYZN, 
    const double* bondLevelRotationTensorZXN, 
    const double* bondLevelRotationTensorZYN, 
    const double* bondLevelRotationTensorZZN, 
    double* bondLevelLeftStretchTensorXXNP1,
    double* bondLevelLeftStretchTensorXYNP1,
    double* bondLevelLeftStretchTensorXZNP1,
    double* bondLevelLeftStretchTensorYXNP1,
    double* bondLevelLeftStretchTensorYYNP1,
    double* bondLevelLeftStretchTensorYZNP1,
    double* bondLevelLeftStretchTensorZXNP1,
    double* bondLevelLeftStretchTensorZYNP1,
    double* bondLevelLeftStretchTensorZZNP1,
    double* bondLevelRotationTensorXXNP1,
    double* bondLevelRotationTensorXYNP1,
    double* bondLevelRotationTensorXZNP1,
    double* bondLevelRotationTensorYXNP1,
    double* bondLevelRotationTensorYYNP1,
    double* bondLevelRotationTensorYZNP1,
    double* bondLevelRotationTensorZXNP1,
    double* bondLevelRotationTensorZYNP1,
    double* bondLevelRotationTensorZZNP1,
    double* bondLevelUnrotatedRateOfDeformationXX,
    double* bondLevelUnrotatedRateOfDeformationXY,
    double* bondLevelUnrotatedRateOfDeformationXZ,
    double* bondLevelUnrotatedRateOfDeformationYX,
    double* bondLevelUnrotatedRateOfDeformationYY,
    double* bondLevelUnrotatedRateOfDeformationYZ,
    double* bondLevelUnrotatedRateOfDeformationZX,
    double* bondLevelUnrotatedRateOfDeformationZY,
    double* bondLevelUnrotatedRateOfDeformationZZ,
    const double* influenceState,
    const int* neighborhoodList,
    int numPoints,
    double dt
);


template void rotateBondLevelCauchyStress
(
    const double* bondLevelRotationTensorXX,
    const double* bondLevelRotationTensorXY,
    const double* bondLevelRotationTensorXZ,
    const double* bondLevelRotationTensorYX,
    const double* bondLevelRotationTensorYY,
    const double* bondLevelRotationTensorYZ,
    const double* bondLevelRotationTensorZX,
    const double* bondLevelRotationTensorZY,
    const double* bondLevelRotationTensorZZ,
    const double* bondLevelUnrotatedCauchyStressXX,
    const double* bondLevelUnrotatedCauchyStressXY,
    const double* bondLevelUnrotatedCauchyStressXZ,
    const double* bondLevelUnrotatedCauchyStressYX,
    const double* bondLevelUnrotatedCauchyStressYY,
    const double* bondLevelUnrotatedCauchyStressYZ,
    const double* bondLevelUnrotatedCauchyStressZX,
    const double* bondLevelUnrotatedCauchyStressZY,
    const double* bondLevelUnrotatedCauchyStressZZ,
    double* bondLevelRotatedCauchyStressXX,
    double* bondLevelRotatedCauchyStressXY,
    double* bondLevelRotatedCauchyStressXZ,
    double* bondLevelRotatedCauchyStressYX,
    double* bondLevelRotatedCauchyStressYY,
    double* bondLevelRotatedCauchyStressYZ,
    double* bondLevelRotatedCauchyStressZX,
    double* bondLevelRotatedCauchyStressZY,
    double* bondLevelRotatedCauchyStressZZ,
    const int* neighborhoodList,
    int numPoints
);

template void computeBondLevelPiolaStress<double>
(
    const double* bondLevelJacobianDeterminant,
    const double* bondLevelCauchyStressXX,
    const double* bondLevelCauchyStressXY,
    const double* bondLevelCauchyStressXZ,
    const double* bondLevelCauchyStressYX,
    const double* bondLevelCauchyStressYY,
    const double* bondLevelCauchyStressYZ,
    const double* bondLevelCauchyStressZX,
    const double* bondLevelCauchyStressZY,
    const double* bondLevelCauchyStressZZ,
    const double* bondLevelDeformationGradientInvXX,
    const double* bondLevelDeformationGradientInvXY,
    const double* bondLevelDeformationGradientInvXZ,
    const double* bondLevelDeformationGradientInvYX,
    const double* bondLevelDeformationGradientInvYY,
    const double* bondLevelDeformationGradientInvYZ,
    const double* bondLevelDeformationGradientInvZX,
    const double* bondLevelDeformationGradientInvZY,
    const double* bondLevelDeformationGradientInvZZ,
    double* bondLevelPiolaStressXX,
    double* bondLevelPiolaStressXY,
    double* bondLevelPiolaStressXZ,
    double* bondLevelPiolaStressYX,
    double* bondLevelPiolaStressYY,
    double* bondLevelPiolaStressYZ,
    double* bondLevelPiolaStressZX,
    double* bondLevelPiolaStressZY,
    double* bondLevelPiolaStressZZ,
    const int* neighborhoodList,
    int numPoints
);

template void computeStressIntegral<double>
(
    const double* volume,
    const double* weightedVolume,
    const double* modelCoordinates,
    const double* bondLevelStressXX,
    const double* bondLevelStressXY,
    const double* bondLevelStressXZ,
    const double* bondLevelStressYX,
    const double* bondLevelStressYY,
    const double* bondLevelStressYZ,
    const double* bondLevelStressZX,
    const double* bondLevelStressZY,
    const double* bondLevelStressZZ,
    double* stressIntegral,
    const double* influenceState,
    const int* neighborhoodList,
    int numPoints
);


}
