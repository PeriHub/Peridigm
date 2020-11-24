//! \file correspondence.cxx

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

#include "correspondence.h"
#include "material_utilities.h"
#include <Sacado.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <math.h>
#include <functional>
#include <boost/math/constants/constants.hpp>
#include <vector>

namespace CORRESPONDENCE {

template<typename ScalarT>
void setOnesOnDiagonalFullTensor(ScalarT* tensor, int numPoints){
 
  ScalarT *tens = tensor;

  for(int iID=0; iID<numPoints; ++iID, tens+=9){
      *(tens) = 1.0;
      *(tens+4) = 1.0;
      *(tens+8) = 1.0;
  }  
}

template<typename ScalarT>
int Invert3by3Matrix
(
 const ScalarT* matrix,
 ScalarT& determinant,
 ScalarT* inverse
)
{
  int returnCode(0);

  ScalarT minor0 =  *(matrix+4) * *(matrix+8) - *(matrix+5) * *(matrix+7);
  ScalarT minor1 =  *(matrix+3) * *(matrix+8) - *(matrix+5) * *(matrix+6);
  ScalarT minor2 =  *(matrix+3) * *(matrix+7) - *(matrix+4) * *(matrix+6);
  ScalarT minor3 =  *(matrix+1) * *(matrix+8) - *(matrix+2) * *(matrix+7);
  ScalarT minor4 =  *(matrix)   * *(matrix+8) - *(matrix+6) * *(matrix+2);
  ScalarT minor5 =  *(matrix)   * *(matrix+7) - *(matrix+1) * *(matrix+6);
  ScalarT minor6 =  *(matrix+1) * *(matrix+5) - *(matrix+2) * *(matrix+4);
  ScalarT minor7 =  *(matrix)   * *(matrix+5) - *(matrix+2) * *(matrix+3);
  ScalarT minor8 =  *(matrix)   * *(matrix+4) - *(matrix+1) * *(matrix+3);
  determinant = *(matrix) * minor0 - *(matrix+1) * minor1 + *(matrix+2) * minor2;

  if(determinant == ScalarT(0.0)){
    returnCode = 1;
    *(inverse) = 0.0;
    *(inverse+1) = 0.0;
    *(inverse+2) = 0.0;
    *(inverse+3) = 0.0;
    *(inverse+4) = 0.0;
    *(inverse+5) = 0.0;
    *(inverse+6) = 0.0;
    *(inverse+7) = 0.0;
    *(inverse+8) = 0.0;
  }
  else{
    *(inverse) = minor0/determinant;
    *(inverse+1) = -1.0*minor3/determinant;
    *(inverse+2) = minor6/determinant;
    *(inverse+3) = -1.0*minor1/determinant;
    *(inverse+4) = minor4/determinant;
    *(inverse+5) = -1.0*minor7/determinant;
    *(inverse+6) = minor2/determinant;
    *(inverse+7) = -1.0*minor5/determinant;
    *(inverse+8) = minor8/determinant;
  }

  return returnCode;
}
  
//inversionReturnCode = Invert2by2Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);
template<typename ScalarT>
int Invert2by2Matrix
(
 const ScalarT* matrix,
 ScalarT& determinant,
 ScalarT* inverse
)
{
  int returnCode(0);
  ScalarT a =  *(matrix);
  ScalarT b =  *(matrix+1);
  ScalarT c =  *(matrix+3);
  ScalarT d =  *(matrix+4);
  
  
  determinant = a*d - b*c ;

  if(determinant == ScalarT(0.0)){
      returnCode = 1;
      *(inverse)   = 0.0;
      *(inverse+1) = 0.0;
      *(inverse+2) = 0.0;
      *(inverse+3) = 0.0;
      *(inverse+4) = 0.0;
      *(inverse+5) = 0.0;
      *(inverse+6) = 0.0;
      *(inverse+7) = 0.0;
      *(inverse+8) = 0.0;
  }
  else{
      *(inverse)   = d/determinant;
      *(inverse+1) = -1.0 * b/determinant;
      *(inverse+2) = 0.0;
      *(inverse+3) = -1.0 * c/determinant;
      *(inverse+4) = a/determinant;
      *(inverse+5) = 0.0;
      *(inverse+6) = 0.0;
      *(inverse+7) = 0.0;
      *(inverse+8) = 0.0;
  }
    
  return returnCode;
}


template<typename ScalarT>
void TransposeMatrix
(
 const ScalarT* matrix,
 ScalarT* transpose
)
{
  // Store some values so that the matrix and transpose can be the
  // same matrix (i.e., transpose in place)
  ScalarT temp_xy( *(matrix+1) );
  ScalarT temp_xz( *(matrix+2) );
  ScalarT temp_yz( *(matrix+5) );

  *(transpose)   = *(matrix);
  *(transpose+1) = *(matrix+3);
  *(transpose+2) = *(matrix+6);
  *(transpose+3) = temp_xy;
  *(transpose+4) = *(matrix+4);
  *(transpose+5) = *(matrix+7);
  *(transpose+6) = temp_xz;
  *(transpose+7) = temp_yz;
  *(transpose+8) = *(matrix+8);
}

template<typename ScalarT>
void MatrixMultiply
(
 bool transA,
 bool transB,
 ScalarT alpha,
 const ScalarT* a,
 const ScalarT* b,
 ScalarT* result
)
{
  // This function computes result = alpha * a * b
  // where alpha is a scalar and a and b are 3x3 matrices
  // The arguments transA and transB denote whether or not
  // to use the transpose of a and b, respectively.

  // The default ordering is row-major:
  //
  // XX(0) XY(1) XZ(2)
  // YX(3) YY(4) YZ(5)
  // ZX(6) ZY(7) ZZ(8)

  if(!transA && !transB){
    *(result+0) = *(a+0) * *(b+0) + *(a+1) * *(b+3) + *(a+2) * *(b+6);
    *(result+1) = *(a+0) * *(b+1) + *(a+1) * *(b+4) + *(a+2) * *(b+7);
    *(result+2) = *(a+0) * *(b+2) + *(a+1) * *(b+5) + *(a+2) * *(b+8);
    *(result+3) = *(a+3) * *(b+0) + *(a+4) * *(b+3) + *(a+5) * *(b+6);
    *(result+4) = *(a+3) * *(b+1) + *(a+4) * *(b+4) + *(a+5) * *(b+7);
    *(result+5) = *(a+3) * *(b+2) + *(a+4) * *(b+5) + *(a+5) * *(b+8);
    *(result+6) = *(a+6) * *(b+0) + *(a+7) * *(b+3) + *(a+8) * *(b+6);
    *(result+7) = *(a+6) * *(b+1) + *(a+7) * *(b+4) + *(a+8) * *(b+7);
    *(result+8) = *(a+6) * *(b+2) + *(a+7) * *(b+5) + *(a+8) * *(b+8);
  }
  else if(transA && !transB){
    *(result+0) = *(a+0) * *(b+0) + *(a+3) * *(b+3) + *(a+6) * *(b+6);
    *(result+1) = *(a+0) * *(b+1) + *(a+3) * *(b+4) + *(a+6) * *(b+7);
    *(result+2) = *(a+0) * *(b+2) + *(a+3) * *(b+5) + *(a+6) * *(b+8);
    *(result+3) = *(a+1) * *(b+0) + *(a+4) * *(b+3) + *(a+7) * *(b+6);
    *(result+4) = *(a+1) * *(b+1) + *(a+4) * *(b+4) + *(a+7) * *(b+7);
    *(result+5) = *(a+1) * *(b+2) + *(a+4) * *(b+5) + *(a+7) * *(b+8);
    *(result+6) = *(a+2) * *(b+0) + *(a+5) * *(b+3) + *(a+8) * *(b+6);
    *(result+7) = *(a+2) * *(b+1) + *(a+5) * *(b+4) + *(a+8) * *(b+7);
    *(result+8) = *(a+2) * *(b+2) + *(a+5) * *(b+5) + *(a+8) * *(b+8);
  }
  else if(!transA && transB){
    *(result+0) = *(a+0) * *(b+0) + *(a+1) * *(b+1) + *(a+2) * *(b+2);
    *(result+1) = *(a+0) * *(b+3) + *(a+1) * *(b+4) + *(a+2) * *(b+5);
    *(result+2) = *(a+0) * *(b+6) + *(a+1) * *(b+7) + *(a+2) * *(b+8);
    *(result+3) = *(a+3) * *(b+0) + *(a+4) * *(b+1) + *(a+5) * *(b+2);
    *(result+4) = *(a+3) * *(b+3) + *(a+4) * *(b+4) + *(a+5) * *(b+5);
    *(result+5) = *(a+3) * *(b+6) + *(a+4) * *(b+7) + *(a+5) * *(b+8);
    *(result+6) = *(a+6) * *(b+0) + *(a+7) * *(b+1) + *(a+8) * *(b+2);
    *(result+7) = *(a+6) * *(b+3) + *(a+7) * *(b+4) + *(a+8) * *(b+5);
    *(result+8) = *(a+6) * *(b+6) + *(a+7) * *(b+7) + *(a+8) * *(b+8);
  }
  else{
    *(result+0) = *(a+0) * *(b+0) + *(a+3) * *(b+1) + *(a+6) * *(b+2);
    *(result+1) = *(a+0) * *(b+3) + *(a+3) * *(b+4) + *(a+6) * *(b+5);
    *(result+2) = *(a+0) * *(b+6) + *(a+3) * *(b+7) + *(a+6) * *(b+8);
    *(result+3) = *(a+1) * *(b+0) + *(a+4) * *(b+1) + *(a+7) * *(b+2);
    *(result+4) = *(a+1) * *(b+3) + *(a+4) * *(b+4) + *(a+7) * *(b+5);
    *(result+5) = *(a+1) * *(b+6) + *(a+4) * *(b+7) + *(a+7) * *(b+8);
    *(result+6) = *(a+2) * *(b+0) + *(a+5) * *(b+1) + *(a+8) * *(b+2);
    *(result+7) = *(a+2) * *(b+3) + *(a+5) * *(b+4) + *(a+8) * *(b+5);
    *(result+8) = *(a+2) * *(b+6) + *(a+5) * *(b+7) + *(a+8) * *(b+8);
  }

  if(alpha != 1.0){
    for(int i=0 ; i<9 ; ++i)
      *(result+i) *= alpha;
  }
}

template<typename ScalarT>
int computeShapeTensorInverseAndApproximateDeformationGradient
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const ScalarT* coordinates,
ScalarT* shapeTensorInverse,
ScalarT* deformationGradient,
const double* bondDamage,
const int* neighborhoodList,
int numPoints,
const double* detachedNodes,
const bool type
)
{
  int returnCode = 0;

  const double* delta = horizon;
  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinates;
  
  const ScalarT* neighborCoord;
  ScalarT* shapeTensorInv = shapeTensorInverse;
  ScalarT* defGrad = deformationGradient;
  
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  ScalarT deformedBondX, deformedBondY, deformedBondZ;
  double neighborVolume, omega, temp;

  std::vector<ScalarT> shapeTensorVector(9);
  ScalarT* shapeTensor = &shapeTensorVector[0];
  ScalarT shapeTensorDeterminant;

  std::vector<ScalarT> defGradFirstTermVector(9);
  ScalarT* defGradFirstTerm = &defGradFirstTermVector[0];

  // placeholder for bond damage
  //double bondDamage = 0.0;

  int inversionReturnCode(0);

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, coord+=3,
        shapeTensorInv+=9, defGrad+=9){

    // Zero out data
    *(shapeTensor)   = 0.0 ; *(shapeTensor+1) = 0.0 ; *(shapeTensor+2) = 0.0 ;
    *(shapeTensor+3) = 0.0 ; *(shapeTensor+4) = 0.0 ; *(shapeTensor+5) = 0.0 ;
    *(shapeTensor+6) = 0.0 ; *(shapeTensor+7) = 0.0 ; *(shapeTensor+8) = 0.0 ;
    *(defGradFirstTerm)   = 0.0 ; *(defGradFirstTerm+1) = 0.0 ; *(defGradFirstTerm+2) = 0.0 ;
    *(defGradFirstTerm+3) = 0.0 ; *(defGradFirstTerm+4) = 0.0 ; *(defGradFirstTerm+5) = 0.0 ;
    *(defGradFirstTerm+6) = 0.0 ; *(defGradFirstTerm+7) = 0.0 ; *(defGradFirstTerm+8) = 0.0 ;

    numNeighbors = *neighborListPtr; neighborListPtr++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamage++){
      
      neighborIndex = *neighborListPtr;
      // if (*(detachedNodes+neighborIndex) !=0) continue; // if completely synchronized it is not needed --> one core 

      neighborVolume = volume[neighborIndex];
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborCoord = coordinates + 3*neighborIndex;
      if (*(detachedNodes+iID)!=0) continue;
      if (*(detachedNodes+neighborIndex)!=0) continue;
      undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);
      
      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);
      if (deformedBondX == 0 and deformedBondY == 0 and deformedBondZ == 0) continue;
      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta);
      
      
      temp = (1.0 - *bondDamage) * omega * neighborVolume;
      *(shapeTensor)   += temp * undeformedBondX * undeformedBondX;
      *(shapeTensor+1) += temp * undeformedBondX * undeformedBondY;
      
      *(shapeTensor+2) += temp * undeformedBondX * undeformedBondZ;
      
      *(shapeTensor+3) += temp * undeformedBondY * undeformedBondX;
      *(shapeTensor+4) += temp * undeformedBondY * undeformedBondY;
      
      *(shapeTensor+5) += temp * undeformedBondY * undeformedBondZ;
      *(shapeTensor+6) += temp * undeformedBondZ * undeformedBondX;
      *(shapeTensor+7) += temp * undeformedBondZ * undeformedBondY;
      *(shapeTensor+8) += temp * undeformedBondZ * undeformedBondZ;
      

      *(defGradFirstTerm)   += temp * deformedBondX * undeformedBondX;
      *(defGradFirstTerm+1) += temp * deformedBondX * undeformedBondY;
      *(defGradFirstTerm+2) += temp * deformedBondX * undeformedBondZ;
      *(defGradFirstTerm+3) += temp * deformedBondY * undeformedBondX;
      *(defGradFirstTerm+4) += temp * deformedBondY * undeformedBondY;
      *(defGradFirstTerm+5) += temp * deformedBondY * undeformedBondZ;
      *(defGradFirstTerm+6) += temp * deformedBondZ * undeformedBondX;
      *(defGradFirstTerm+7) += temp * deformedBondZ * undeformedBondY;
      *(defGradFirstTerm+8) += temp * deformedBondZ * undeformedBondZ;

    }
    
    if (*(detachedNodes+iID) == 0) {
        
        if (type==true){
            inversionReturnCode = Invert2by2Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);

        }
        else{
            inversionReturnCode = Invert3by3Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);
        }
            
        MatrixMultiply(false, false, 1.0, defGradFirstTerm, shapeTensorInv, defGrad);

        if(inversionReturnCode > 0){
           returnCode = inversionReturnCode;
           
           }
        }
    // Matrix multiply the first term and the shape tensor inverse to compute
    // the deformation gradient
    
  }

  return returnCode;
}

//Performs kinematic computations following Flanagan and Taylor (1987), returns
//unrotated rate-of-deformation and rotation tensors
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
const double* detachedNodes,
const bool type
)
{
  int returnCode = 0;

  const double* delta = horizon;
  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  ScalarT* defGrad = deformationGradient;
  const ScalarT* shapeTensorInv = shapeTensorInverse;
  ScalarT* leftStretchN = leftStretchTensorN;
  const ScalarT* rotTensorN = rotationTensorN;

  ScalarT* leftStretchNP1 = leftStretchTensorNP1;
  ScalarT* rotTensorNP1 = rotationTensorNP1;
  ScalarT* unrotRateOfDef = unrotatedRateOfDeformation;

  std::vector<ScalarT> FdotFirstTermVector(9) ; ScalarT* FdotFirstTerm = &FdotFirstTermVector[0];
  std::vector<ScalarT> FdotVector(9) ; ScalarT* Fdot = &FdotVector[0];
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
  ScalarT velStateX, velStateY, velStateZ;
  ScalarT traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double neighborVolume, omega, scalarTemp; 
  int inversionReturnCode(0);
  
  // placeholder for bond damage
 // double bondDamage = 0.0;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, vel+=3,
        shapeTensorInv+=9, rotTensorN+=9, rotTensorNP1+=9, leftStretchNP1+=9, leftStretchN+=9,
        unrotRateOfDef+=9, defGrad+=9){

    // Initialize data
    *(FdotFirstTerm)   = 0.0 ; *(FdotFirstTerm+1) = 0.0 ;  *(FdotFirstTerm+2) = 0.0;
    *(FdotFirstTerm+3) = 0.0 ; *(FdotFirstTerm+4) = 0.0 ;  *(FdotFirstTerm+5) = 0.0;
    *(FdotFirstTerm+6) = 0.0 ; *(FdotFirstTerm+7) = 0.0 ;  *(FdotFirstTerm+8) = 0.0;
    
    //Compute Fdot
    numNeighbors = *neighborListPtr; neighborListPtr++;
    
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){
      neighborIndex = *neighborListPtr;    
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborVel = velocities + 3*neighborIndex;
      if (*(detachedNodes+iID)!=0) continue;
      if (*(detachedNodes+neighborIndex)!=0) continue;
      neighborVolume = volume[neighborIndex];
      undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      // The velState is the relative difference in velocities of the nodes at
      // each end of a bond. i.e., v_j - v_i
      velStateX = *(neighborVel)   - *(vel);
      velStateY = *(neighborVel+1) - *(vel+1);
      velStateZ = *(neighborVel+2) - *(vel+2);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta);

      scalarTemp = (1.0 - *bondDamage) * omega * neighborVolume;
     // scalarTemp =  omega * neighborVolume;
      *(FdotFirstTerm)   += scalarTemp * velStateX * undeformedBondX;
      *(FdotFirstTerm+1) += scalarTemp * velStateX * undeformedBondY;
      *(FdotFirstTerm+2) += scalarTemp * velStateX * undeformedBondZ;
      *(FdotFirstTerm+3) += scalarTemp * velStateY * undeformedBondX;
      *(FdotFirstTerm+4) += scalarTemp * velStateY * undeformedBondY;
      *(FdotFirstTerm+5) += scalarTemp * velStateY * undeformedBondZ;
      *(FdotFirstTerm+6) += scalarTemp * velStateZ * undeformedBondX;
      *(FdotFirstTerm+7) += scalarTemp * velStateZ * undeformedBondY;
      *(FdotFirstTerm+8) += scalarTemp * velStateZ * undeformedBondZ;
    }
    if (*(detachedNodes+iID)!=0) continue;
    // Compute Fdot
    MatrixMultiply(false, false, 1.0, FdotFirstTerm, shapeTensorInv, Fdot);

    // Compute the inverse of the deformation gradient, Finverse
    if (type==true){
        inversionReturnCode = Invert2by2Matrix(defGrad, determinant, Finverse);
    }
    else{
        inversionReturnCode = Invert3by3Matrix(defGrad, determinant, Finverse);
        }
    if(inversionReturnCode > 0){
        returnCode = 2;
        return returnCode;
    }
    // Compute the Eulerian velocity gradient L = Fdot * Finv
    MatrixMultiply(false, false, 1.0, Fdot, Finverse, eulerianVelGrad);

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
    if (type==true){
        inversionReturnCode = Invert2by2Matrix(temp, determinant, tempInv);
    }
    else{
        inversionReturnCode = Invert3by3Matrix(temp, determinant, tempInv);
        }
    if(inversionReturnCode > 0){
      returnCode = inversionReturnCode;
      return returnCode;}

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
    if ( OmegaSq > 1.e-30){

      // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
      //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
      scaleFactor1 = sin(dt*Omega) / Omega;
      scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
      MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, OmegaTensorSq);
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
    MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, rotTensorNP1);

    // Compute rate of stretch, Vdot = L*V - V*Omega
    // First tempA = L*V, 
    MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, tempA);

    // tempB = V*Omega
    MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, tempB);

    //Vdot = tempA - tempB
    for(int i=0 ; i<9 ; ++i)
      *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

    //V_STEP_NP1 = V_STEP_N + dt*Vdot
    for(int i=0 ; i<9 ; ++i)
      *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

    // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
    MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, temp);

    // d = Rt * temp
    MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, unrotRateOfDef);
  }

  return returnCode;
}




template<typename ScalarT>
void computeGreenLagrangeStrain
(
  const ScalarT* deformationGradientXX,
  const ScalarT* deformationGradientXY,
  const ScalarT* deformationGradientXZ,
  const ScalarT* deformationGradientYX,
  const ScalarT* deformationGradientYY,
  const ScalarT* deformationGradientYZ,
  const ScalarT* deformationGradientZX,
  const ScalarT* deformationGradientZY,
  const ScalarT* deformationGradientZZ,
  ScalarT* greenLagrangeStrainXX,
  ScalarT* greenLagrangeStrainXY,
  ScalarT* greenLagrangeStrainXZ,
  ScalarT* greenLagrangeStrainYX,
  ScalarT* greenLagrangeStrainYY,
  ScalarT* greenLagrangeStrainYZ,
  ScalarT* greenLagrangeStrainZX,
  ScalarT* greenLagrangeStrainZY,
  ScalarT* greenLagrangeStrainZZ,
  int numPoints
)
{
  // Green-Lagrange Strain E = 0.5*(F^T F - I)

  const ScalarT* defGradXX = deformationGradientXX;
  const ScalarT* defGradXY = deformationGradientXY;
  const ScalarT* defGradXZ = deformationGradientXZ;
  const ScalarT* defGradYX = deformationGradientYX;
  const ScalarT* defGradYY = deformationGradientYY;
  const ScalarT* defGradYZ = deformationGradientYZ;
  const ScalarT* defGradZX = deformationGradientZX;
  const ScalarT* defGradZY = deformationGradientZY;
  const ScalarT* defGradZZ = deformationGradientZZ;
  ScalarT* strainXX = greenLagrangeStrainXX;
  ScalarT* strainXY = greenLagrangeStrainXY;
  ScalarT* strainXZ = greenLagrangeStrainXZ;
  ScalarT* strainYX = greenLagrangeStrainYX;
  ScalarT* strainYY = greenLagrangeStrainYY;
  ScalarT* strainYZ = greenLagrangeStrainYZ;
  ScalarT* strainZX = greenLagrangeStrainZX;
  ScalarT* strainZY = greenLagrangeStrainZY;
  ScalarT* strainZZ = greenLagrangeStrainZZ;

  for(int iID=0 ; iID<numPoints ; ++iID, 
        ++defGradXX, ++defGradXY, ++defGradXZ,
        ++defGradYX, ++defGradYY, ++defGradYZ,
        ++defGradZX, ++defGradZY, ++defGradZZ,
        ++strainXX, ++strainXY, ++strainXZ,
        ++strainYX, ++strainYY, ++strainYZ,
        ++strainZX, ++strainZY, ++strainZZ){

    *strainXX = 0.5 * ( *(defGradXX) * *(defGradXX) + *(defGradYX) * *(defGradYX) + *(defGradZX) * *(defGradZX) - 1.0 );
    *strainXY = 0.5 * ( *(defGradXX) * *(defGradXY) + *(defGradYX) * *(defGradYY) + *(defGradZX) * *(defGradZY) );
    *strainXZ = 0.5 * ( *(defGradXX) * *(defGradXZ) + *(defGradYX) * *(defGradYZ) + *(defGradZX) * *(defGradZZ) );
    *strainYX = 0.5 * ( *(defGradXY) * *(defGradXX) + *(defGradYY) * *(defGradYX) + *(defGradZY) * *(defGradZX) );
    *strainYY = 0.5 * ( *(defGradXY) * *(defGradXY) + *(defGradYY) * *(defGradYY) + *(defGradZY) * *(defGradZY) - 1.0 );
    *strainYZ = 0.5 * ( *(defGradXY) * *(defGradXZ) + *(defGradYY) * *(defGradYZ) + *(defGradZY) * *(defGradZZ) );
    *strainZX = 0.5 * ( *(defGradXZ) * *(defGradXX) + *(defGradYZ) * *(defGradYX) + *(defGradZZ) * *(defGradZX) );
    *strainZY = 0.5 * ( *(defGradXZ) * *(defGradXY) + *(defGradYZ) * *(defGradYY) + *(defGradZZ) * *(defGradZY) );
    *strainZZ = 0.5 * ( *(defGradXZ) * *(defGradXZ) + *(defGradYZ) * *(defGradYZ) + *(defGradZZ) * *(defGradZZ) - 1.0 );
  }
}

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
const double* detachedNodes,
const double* bondDamage
)
{
  double vol, neighborVol;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;
  ScalarT expectedNeighborLocationX, expectedNeighborLocationY, expectedNeighborLocationZ;
  ScalarT hourglassVectorX, hourglassVectorY, hourglassVectorZ;
  ScalarT dot, magnitude;
  int neighborIndex, numNeighbors;
  const ScalarT* defGrad = deformationGradient;
  const double* delta = horizon;
  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  ScalarT* hourglassForceDensityPtr = hourglassForceDensity;
  ScalarT* neighborHourglassForceDensityPtr;

  // placeholder for inclusion of bond damage
  //double bondDamage = 0.0;

  const double pi = boost::math::constants::pi<double>();
  double firstPartOfConstant = 18.0*hourglassCoefficient*bulkModulus/pi;
  double constant;

  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, coord+=3,
        defGrad+=9, hourglassForceDensityPtr+=3){

    constant = firstPartOfConstant/( (*delta)*(*delta)*(*delta)*(*delta) );

    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamage++){
      neighborIndex = *neighborListPtr;
      
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborCoord = coordinates + 3*neighborIndex;
      if (detachedNodes[neighborIndex]!=0) continue;
      
      undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);
      deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                deformedBondY*deformedBondY +
                                deformedBondZ*deformedBondZ);

      expectedNeighborLocationX = *(coord) +
        *(defGrad) * undeformedBondX +
        *(defGrad+1) * undeformedBondY +
        *(defGrad+2) * undeformedBondZ;
      expectedNeighborLocationY = *(coord+1) +
        *(defGrad+3) * undeformedBondX +
        *(defGrad+4) * undeformedBondY +
        *(defGrad+5) * undeformedBondZ;
      expectedNeighborLocationZ = *(coord+2) +
        *(defGrad+6) * undeformedBondX +
        *(defGrad+7) * undeformedBondY +
        *(defGrad+8) * undeformedBondZ;

      hourglassVectorX = expectedNeighborLocationX - *(neighborCoord);
      hourglassVectorY = expectedNeighborLocationY - *(neighborCoord+1);
      hourglassVectorZ = expectedNeighborLocationZ - *(neighborCoord+2);

      dot = hourglassVectorX*deformedBondX + hourglassVectorY*deformedBondY + hourglassVectorZ*deformedBondZ;
      dot *= -1.0;

      magnitude = (1.0-*bondDamage) * constant * (dot/undeformedBondLength) * (1.0/deformedBondLength);

      vol = volume[iID];
      neighborVol = volume[neighborIndex];
      neighborHourglassForceDensityPtr = hourglassForceDensity + 3*neighborIndex;

      *(hourglassForceDensityPtr)   += magnitude * deformedBondX * neighborVol;
      *(hourglassForceDensityPtr+1) += magnitude * deformedBondY * neighborVol;
      *(hourglassForceDensityPtr+2) += magnitude * deformedBondZ * neighborVol;

      *(neighborHourglassForceDensityPtr)   -= magnitude * deformedBondX * vol;
      *(neighborHourglassForceDensityPtr+1) -= magnitude * deformedBondY * vol;
      *(neighborHourglassForceDensityPtr+2) -= magnitude * deformedBondZ * vol;

    }
  }
}



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
)
{
  // S.A. Silling, "Stability of peridynamic correspondence material models and their
 // particle discretizations" in Comput. Methods Appl. Mech. Engrg. 322 (2017) 42–57,http://dx.doi.org/10.1016/j.cma.2017.03.043

  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  ScalarT deformedBondX, deformedBondY, deformedBondZ;
  ScalarT  magnitude;
  int neighborIndex, numNeighbors;
  const ScalarT* defGrad = deformationGradient;
  const double* delta1 = horizon;
  const double* delta2 = horizon;
  const double* modelCoord1 = modelCoordinates;
  const double* modelCoord2 = modelCoordinates;
  const double* neighborModelCoord;
  //const ScalarT* coord1 = coordinates;
  const ScalarT* coord2 = coordinates;
  const double* bondDamage1 = bondDamage;
  const double* bondDamage2 = bondDamage;
  double nonUniformDeformState[3], m_omega;
  ScalarT FxsiX, FxsiY, FxsiZ;
  const ScalarT* neighborCoord;
  double vol, neighborVol;
  ScalarT* hourglassForceDensityPtr = hourglassForceDensity;
  ScalarT* neighborHourglassForceDensityPtr;

  // placeholder for inclusion of bond damage
  //double bondDamage = 0.0;

  const double pi = boost::math::constants::pi<double>();
  double firstPartOfConstant = 18.0*hourglassCoefficient*bulkModulus/pi;
  double constant;
  double omega0 = 0.0;
  const int *neighborListPtr1 = neighborhoodList;
  const int *neighborListPtr2 = neighborhoodList;
  
  for(int iID=0 ; iID<numPoints ; ++iID, delta1++, modelCoord1+=3){

    numNeighbors = *neighborListPtr1; neighborListPtr1++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr1++, bondDamage1++){
      neighborIndex = *neighborListPtr1;
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      
      undeformedBondX = *(neighborModelCoord)   - *(modelCoord1);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord1+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord1+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);
      m_omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta1);
      neighborVol = volume[neighborIndex];
      omega0 += (1-*bondDamage1)*m_omega;
      
    }
  }
omega0 = 1;
  
  for(int iID=0 ; iID<numPoints ; ++iID, delta2++, modelCoord2+=3, coord2+=3,
        defGrad+=9, hourglassForceDensityPtr+=3){
    constant = firstPartOfConstant/( omega0 * *(delta2)* *(delta2)* *(delta2)* *(delta2)* *(delta2) );
    numNeighbors = *neighborListPtr2; neighborListPtr2++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr2++, bondDamage2++){
      neighborIndex = *neighborListPtr2;
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborCoord = coordinates + 3*neighborIndex;

      undeformedBondX = *(neighborModelCoord)   - *(modelCoord2);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord2+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord2+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);
      m_omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta2);
      
      
      
      deformedBondX = *(neighborCoord)   - *(coord2);
      deformedBondY = *(neighborCoord+1) - *(coord2+1);
      deformedBondZ = *(neighborCoord+2) - *(coord2+2);
      //deformedBondLength = sqrt(deformedBondX*deformedBondX +
      //                        deformedBondY*deformedBondY +
      //                        deformedBondZ*deformedBondZ);

      FxsiX = 
        *(defGrad)   * undeformedBondX +
        *(defGrad+1) * undeformedBondY +
        *(defGrad+2) * undeformedBondZ;
      FxsiY = 
        *(defGrad+3) * undeformedBondX +
        *(defGrad+4) * undeformedBondY +
        *(defGrad+5) * undeformedBondZ;
      FxsiZ = 
        *(defGrad+6) * undeformedBondX +
        *(defGrad+7) * undeformedBondY +
        *(defGrad+8) * undeformedBondZ;

      nonUniformDeformState[0] = deformedBondX - FxsiX;
      nonUniformDeformState[1] = deformedBondY - FxsiY;
      nonUniformDeformState[2] = deformedBondZ - FxsiZ;

      magnitude = (1.0-*bondDamage2) * constant * m_omega;

      vol = volume[iID];
      neighborVol = volume[neighborIndex];
      neighborHourglassForceDensityPtr = hourglassForceDensity + 3*neighborIndex;

      *(hourglassForceDensityPtr)   += magnitude * nonUniformDeformState[0] * neighborVol;
      *(hourglassForceDensityPtr+1) += magnitude * nonUniformDeformState[1] * neighborVol;
      *(hourglassForceDensityPtr+2) += magnitude * nonUniformDeformState[2] * neighborVol;

      *(neighborHourglassForceDensityPtr)   -= magnitude * nonUniformDeformState[0] * vol;
      *(neighborHourglassForceDensityPtr+1) -= magnitude * nonUniformDeformState[1] * vol;
      *(neighborHourglassForceDensityPtr+2) -= magnitude * nonUniformDeformState[2] * vol;

    }
  }
}



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
ScalarT* shapeTensorInverse,
const double* C,
const double* bondDamage,
const double* detachedNodes,
ScalarT* hourglassForceDensity,
double hourglassCoefficient
)
{
  // J. Wan et al., "Improved method for zero-energy mode suppression in peridynamic correspondence model" 
  // in Acta Mechanica Sinica https://doi.org/10.1007/s10409-019-00873-y

  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  ScalarT deformedBondX, deformedBondY, deformedBondZ;
  int neighborIndex, numNeighbors;
  const ScalarT* defGrad = deformationGradient;
  const double* delta = horizon;
  const double* modelCoord = modelCoordinates;
  const int *neighborListPtr = neighborhoodList;
  const double* neighborModelCoord;
  //const ScalarT* coord1 = coordinates;
  const ScalarT* coord = coordinates;
  const ScalarT* shapeTensorInv = shapeTensorInverse;
  double nonUniformDeformState[3], m_omega;
  double hourglassStiff[9];
  ScalarT TSx, TSy, TSz;
  ScalarT FxsiX, FxsiY, FxsiZ;
  const ScalarT* neighborCoord;
  double vol, neighborVol, factor;
  ScalarT* hourglassForceDensityPtr = hourglassForceDensity;
  ScalarT* neighborHourglassForceDensityPtr;

  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, coord+=3,
     defGrad+=9, shapeTensorInv+=9, hourglassForceDensityPtr+=3){
    //constant = firstPartOfConstant/( omega0 * *(delta2)* *(delta2)* *(delta2)* *(delta2)* *(delta2) );
    numNeighbors = *neighborListPtr; neighborListPtr++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamage++){
      
      neighborIndex = *neighborListPtr;
      
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborCoord = coordinates + 3*neighborIndex;
      
      if (*(detachedNodes+iID)!=0) continue;
      if (*(detachedNodes+neighborIndex)!=0) continue;
      undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);
      m_omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta);
      
      
      
      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);
      //deformedBondLength = sqrt(deformedBondX*deformedBondX +
      //                        deformedBondY*deformedBondY +
      //                        deformedBondZ*deformedBondZ);

      FxsiX = 
        *(defGrad)   * undeformedBondX +
        *(defGrad+1) * undeformedBondY +
        *(defGrad+2) * undeformedBondZ;
      FxsiY = 
        *(defGrad+3) * undeformedBondX +
        *(defGrad+4) * undeformedBondY +
        *(defGrad+5) * undeformedBondZ;
      FxsiZ = 
        *(defGrad+6) * undeformedBondX +
        *(defGrad+7) * undeformedBondY +
        *(defGrad+8) * undeformedBondZ;

      nonUniformDeformState[0] = deformedBondX - FxsiX;
      nonUniformDeformState[1] = deformedBondY - FxsiY;
      nonUniformDeformState[2] = deformedBondZ - FxsiZ;
      
      hourglassStiff[0] = C[0] * *(shapeTensorInv)  + C[1] * *(shapeTensorInv+4)  + C[2] * *(shapeTensorInv+8)  + C[3] *(*(shapeTensorInv+5) + *(shapeTensorInv+7) ) + C[4] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[5] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[1] = C[30]* *(shapeTensorInv)  + C[31]* *(shapeTensorInv+4)  + C[32]* *(shapeTensorInv+8)  + C[33]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[34] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[35]*(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[2] = C[24]* *(shapeTensorInv)  + C[25]* *(shapeTensorInv+4)  + C[26]* *(shapeTensorInv+8)  + C[27]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[28] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[29]*(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[3] = hourglassStiff[1];
      hourglassStiff[4] = C[6] * *(shapeTensorInv)  + C[7] * *(shapeTensorInv+4)  + C[8] * *(shapeTensorInv+8)  + C[9] *(*(shapeTensorInv+5) + *(shapeTensorInv+7) ) + C[10]*(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[11] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[5] = C[18]* *(shapeTensorInv)  + C[19]* *(shapeTensorInv+4)  + C[20]* *(shapeTensorInv+8)  + C[21]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[22] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[23] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[6] = hourglassStiff[2];
      hourglassStiff[7] = hourglassStiff[5];
      hourglassStiff[8] = C[12]* *(shapeTensorInv)  + C[13]* *(shapeTensorInv+4)  + C[14]* *(shapeTensorInv+8)  + C[15]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[16] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[17] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      
      factor = hourglassCoefficient*(1.0-*bondDamage) * m_omega;
      
      TSx = factor*(hourglassStiff[0] * nonUniformDeformState[0] + hourglassStiff[1] * nonUniformDeformState[1] + hourglassStiff[2] * nonUniformDeformState[2]);
      TSy = factor*(hourglassStiff[3] * nonUniformDeformState[0] + hourglassStiff[4] * nonUniformDeformState[1] + hourglassStiff[5] * nonUniformDeformState[2]);
      TSz = factor*(hourglassStiff[6] * nonUniformDeformState[0] + hourglassStiff[7] * nonUniformDeformState[1] + hourglassStiff[8] * nonUniformDeformState[2]);
      
      vol = volume[iID];
      neighborVol = volume[neighborIndex];
      neighborHourglassForceDensityPtr = hourglassForceDensity + 3*neighborIndex;

      *(hourglassForceDensityPtr)   += TSx * neighborVol;
      *(hourglassForceDensityPtr+1) += TSy * neighborVol;
      *(hourglassForceDensityPtr+2) += TSz * neighborVol;

      *(neighborHourglassForceDensityPtr)   -= TSx * vol;
      *(neighborHourglassForceDensityPtr+1) -= TSy * vol;
      *(neighborHourglassForceDensityPtr+2) -= TSz * vol;

    }
  }
  
}

template<typename ScalarT>
void computeCorrespondenceStabilityWanEtAlShort
(
const ScalarT FxsiX,
const ScalarT FxsiY,
const ScalarT FxsiZ,
const double deformedBondX,
const double deformedBondY,
const double deformedBondZ,
const double* hourglassStiff,
ScalarT TSx,
ScalarT TSy,
ScalarT TSz
)
{
  // J. Wan et al., "Improved method for zero-energy mode suppression in peridynamic correspondence model" 
  // in Acta Mechanica Sinica https://doi.org/10.1007/s10409-019-00873-y
    double nonUniformDeformState[3];

    nonUniformDeformState[0] = deformedBondX - FxsiX;
    nonUniformDeformState[1] = deformedBondY - FxsiY;
    nonUniformDeformState[2] = deformedBondZ - FxsiZ;
    
    TSx = (*(hourglassStiff)   * nonUniformDeformState[0] + *(hourglassStiff+1) * nonUniformDeformState[1] + *(hourglassStiff+2) * nonUniformDeformState[2]);
    TSy = (*(hourglassStiff+3) * nonUniformDeformState[0] + *(hourglassStiff+4) * nonUniformDeformState[1] + *(hourglassStiff+5) * nonUniformDeformState[2]);
    TSz = (*(hourglassStiff+6) * nonUniformDeformState[0] + *(hourglassStiff+7) * nonUniformDeformState[1] + *(hourglassStiff+8) * nonUniformDeformState[2]);
    
}

template<typename ScalarT>
void createHourglassStiffness
(
const double* C,
ScalarT* shapeTensorInverse,
double* hourglassStiff
)
{
    const ScalarT* shapeTensorInv = shapeTensorInverse;
    
    *(hourglassStiff)   = C[0] * *(shapeTensorInv)  + C[1] * *(shapeTensorInv+4)  + C[2] * *(shapeTensorInv+8)  + C[3] *(*(shapeTensorInv+5) + *(shapeTensorInv+7) ) + C[4] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[5] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
    *(hourglassStiff+1) = C[30]* *(shapeTensorInv)  + C[31]* *(shapeTensorInv+4)  + C[32]* *(shapeTensorInv+8)  + C[33]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[34] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[35]*(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
    *(hourglassStiff+2) = C[24]* *(shapeTensorInv)  + C[25]* *(shapeTensorInv+4)  + C[26]* *(shapeTensorInv+8)  + C[27]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[28] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[29]*(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
    *(hourglassStiff+3) = hourglassStiff[1];
    *(hourglassStiff+4) = C[6] * *(shapeTensorInv)  + C[7] * *(shapeTensorInv+4)  + C[8] * *(shapeTensorInv+8)  + C[9] *(*(shapeTensorInv+5) + *(shapeTensorInv+7) ) + C[10]*(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[11] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
    *(hourglassStiff+5) = C[18]* *(shapeTensorInv)  + C[19]* *(shapeTensorInv+4)  + C[20]* *(shapeTensorInv+8)  + C[21]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[22] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[23] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
    *(hourglassStiff+6) = hourglassStiff[2];
    *(hourglassStiff+7) = hourglassStiff[5];
    *(hourglassStiff+8) = C[12]* *(shapeTensorInv)  + C[13]* *(shapeTensorInv+4)  + C[14]* *(shapeTensorInv+8)  + C[15]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[16] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[17] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );

}
template<typename ScalarT>
void rotateCauchyStress
(
 const ScalarT* rotationTensor,
 const ScalarT* unrotatedCauchyStress,
 ScalarT* rotatedCauchyStress,
 int numPoints
)
{
  const ScalarT* rotTensor = rotationTensor;
  const ScalarT* unrotatedStress = unrotatedCauchyStress;
  ScalarT* rotatedStress = rotatedCauchyStress;
  ScalarT temp[9];

  for(int iID=0 ; iID<numPoints ; ++iID, 
        rotTensor+=9, unrotatedStress+=9, rotatedStress+=9){ 

      // temp = \sigma_unrot * Rt
      CORRESPONDENCE::MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, temp);
      // \sigma_rot = R * temp
      CORRESPONDENCE::MatrixMultiply(false, false, 1.0, rotTensor, temp, rotatedStress);
  }
}

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
 // set all matrices to values which does not make any problems later on
 // the detached nodes are not considerred
 )
{
    for(int iID=0 ; iID<numPoints ; ++iID,  unrotatedRateOfDeformation+=9,  deformationGradient+=9, leftStretchTensor+=9, rotationTensor+=9, ++detachedNodes){ 
       
        if (*detachedNodes!=0){
            
          for (int i=0; i<9; ++i){
            *(unrotatedRateOfDeformation+i) = 0.0;
            *(deformationGradient+i) = 0.0;
            *(leftStretchTensor+i) = 0.0;
            *(rotationTensor+i) = 0.0;
            *(shapeTensorInverse+i)=0.0;

            }
            *(deformationGradient) = 1.0;
            *(deformationGradient+4) = 1.0;
            *(deformationGradient+8) = 1.0;
            *(leftStretchTensor) = 1.0;
            *(leftStretchTensor+4) = 1.0;
            *(leftStretchTensor+8) = 1.0;
            *(rotationTensor) = 1.0;
            *(rotationTensor+4) = 1.0;
            *(rotationTensor+8) = 1.0;
            *(shapeTensorInverse) = 1.0;
            *(shapeTensorInverse+4) = 1.0;
            *(shapeTensorInverse+8) = 1.0;

        }
  }

}

/** Explicit template instantiation for double. */

template void TransposeMatrix<double>
(
 const double* matrix,
 double* transpose
);

template void MatrixMultiply<double>
(
 bool transA,
 bool transB,
 double alpha,
 const double* a,
 const double* b,
 double* result
);

template void rotateCauchyStress<double>
(
 const double* rotationTensor,
 const double* unrotatedCauchyStress,
 double* rotatedCauchyStress,
 int numPoints
 );

template int Invert3by3Matrix<double>
(
 const double* matrix,
 double& determinant,
 double* inverse
);
template int Invert2by2Matrix<double>
(
 const double* matrix,
 double& determinant,
 double* inverse
);
template int computeShapeTensorInverseAndApproximateDeformationGradient<double>
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const double* coordinates,
double* shapeTensorInverse,
double* deformationGradient,
const double* bondDamage,
const int* neighborhoodList,
int numPoints,
const double* detachedNodes,
const bool type
);

template int computeUnrotatedRateOfDeformationAndRotationTensor<double>
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const double* velocities,
double* deformationGradient,
const double* shapeTensorInverse,
double* leftStretchTensorN,
const double* rotationTensorN,
double* leftStretchTensorNP1,
double* rotationTensorNP1,
double* unrotatedRateOfDeformation,
const int* neighborhoodList,
int numPoints,
double dt,
const double* bondDamage,
const double* detachedNodes,
const bool type
);


template void computeGreenLagrangeStrain<double>
(
  const double* deformationGradientXX,
  const double* deformationGradientXY,
  const double* deformationGradientXZ,
  const double* deformationGradientYX,
  const double* deformationGradientYY,
  const double* deformationGradientYZ,
  const double* deformationGradientZX,
  const double* deformationGradientZY,
  const double* deformationGradientZZ,
  double* greenLagrangeStrainXX,
  double* greenLagrangeStrainXY,
  double* greenLagrangeStrainXZ,
  double* greenLagrangeStrainYX,
  double* greenLagrangeStrainYY,
  double* greenLagrangeStrainYZ,
  double* greenLagrangeStrainZX,
  double* greenLagrangeStrainZY,
  double* greenLagrangeStrainZZ,
  int numPoints
);

template void computeHourglassForce<double>
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const double* coordinates,
const double* deformationGradient,
double* hourglassForceDensity,
const int* neighborhoodList,
int numPoints,
double bulkModulus,
double hourglassCoefficient,
const double* detachedNodes,
const double* bondDamage
);

template void computeCorrespondenceStabilityForce<double>
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const double* coordinates,
const double* deformationGradient,
double* hourglassForceDensity,
const int* neighborhoodList,
int numPoints,
double bulkModulus,
double hourglassCoefficient,
const double* bondDamage
);

template void computeCorrespondenceStabilityWanEtAl<double>
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const double* coordinates,
const int* neighborhoodList,
int numPoints,
const double* deformationGradient,
double* shapeTensorInverse,
const double* C,
const double* bondDamage,
const double* detachedNodes,
double* hourglassForceDensity,
double hourglassCoefficient
);
template void createHourglassStiffness<double>
(
const double* C,
double* shapeTensorInverse,
double* hourglassStiff
);
template void computeCorrespondenceStabilityWanEtAlShort<double>
(
const double FxsiX,
const double FxsiY,
const double FxsiZ,
const double deformedBondX,
const double deformedBondY,
const double deformedBondZ,
const double* hourglassStiff,
double TSx,
double TSy,
double TSz
);
template void setValuesForDetachedNodes<double>
(
 double* deformationGradient,
 double* leftStretchTensor,
 double* rotationTensor,
 double* unrotatedRateOfDeformation,
 double* shapeTensorInverse,
 const double* detachedNodes,
 const int numPoints
 );


template void setOnesOnDiagonalFullTensor<double>
(
 double* tensor,
 int numPoints
);



/** Explicit template instantiation for Sacado::Fad::DFad<double>. */

template void TransposeMatrix<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* matrix,
 Sacado::Fad::DFad<double>* transpose
);

template void MatrixMultiply<Sacado::Fad::DFad<double> >
(
 bool transA,
 bool transB,
 Sacado::Fad::DFad<double> alpha,
 const Sacado::Fad::DFad<double>* a,
 const Sacado::Fad::DFad<double>* b,
 Sacado::Fad::DFad<double>* result
);

template int Invert3by3Matrix<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* matrix,
 Sacado::Fad::DFad<double>& determinant,
 Sacado::Fad::DFad<double>* inverse
);

template int Invert2by2Matrix<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* matrix,
 Sacado::Fad::DFad<double>& determinant,
 Sacado::Fad::DFad<double>* inverse
);

template void computeGreenLagrangeStrain<Sacado::Fad::DFad<double> >
(
  const Sacado::Fad::DFad<double>* deformationGradientXX,
  const Sacado::Fad::DFad<double>* deformationGradientXY,
  const Sacado::Fad::DFad<double>* deformationGradientXZ,
  const Sacado::Fad::DFad<double>* deformationGradientYX,
  const Sacado::Fad::DFad<double>* deformationGradientYY,
  const Sacado::Fad::DFad<double>* deformationGradientYZ,
  const Sacado::Fad::DFad<double>* deformationGradientZX,
  const Sacado::Fad::DFad<double>* deformationGradientZY,
  const Sacado::Fad::DFad<double>* deformationGradientZZ,
  Sacado::Fad::DFad<double>* greenLagrangeStrainXX,
  Sacado::Fad::DFad<double>* greenLagrangeStrainXY,
  Sacado::Fad::DFad<double>* greenLagrangeStrainXZ,
  Sacado::Fad::DFad<double>* greenLagrangeStrainYX,
  Sacado::Fad::DFad<double>* greenLagrangeStrainYY,
  Sacado::Fad::DFad<double>* greenLagrangeStrainYZ,
  Sacado::Fad::DFad<double>* greenLagrangeStrainZX,
  Sacado::Fad::DFad<double>* greenLagrangeStrainZY,
  Sacado::Fad::DFad<double>* greenLagrangeStrainZZ,
  int numPoints
);

}
