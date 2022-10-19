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
#include "temperature_diffusion.h"
#include "matrices.h"
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

template<typename ScalarT>
int EigenVec2D
(
 const ScalarT* a,
 ScalarT* result
)
{
  return MATRICES::EigenVec2D(a, result);
}



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
)
{
    const double* delta = horizon;
    const double* modelCoord = modelCoordinates;
    const double* neighborModelCoord;
    const ScalarT* vel = velocities;
    const ScalarT* neighborVel;
    ScalarT* defGrad = deformationGradient;
    const ScalarT* shapeTensorInv = shapeTensorInverse;

    ScalarT* rateOfDef = unrotatedRateOfDeformation;

    std::vector<ScalarT> FdotFirstTermVector(9) ; ScalarT* FdotFirstTerm = &FdotFirstTermVector[0];
    std::vector<ScalarT> FdotVector(9) ; ScalarT* Fdot = &FdotVector[0];
    std::vector<ScalarT> FinverseVector(9) ; ScalarT* Finverse = &FinverseVector[0];
    std::vector<ScalarT> eulerianVelGradVector(9) ; ScalarT* eulerianVelGrad = &eulerianVelGradVector[0];
    ScalarT determinant;
    ScalarT One = 1.0;
    ScalarT velStateX, velStateY, velStateZ;
    double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
    double neighborVolume, omega, scalarTemp; 
    int inversionReturnCode(0);

    int neighborIndex, numNeighbors;
    const int *neighborListPtr = neighborhoodList;
    //int returnCode;
    for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, vel+=3,
        shapeTensorInv+=9, rateOfDef+=9, defGrad+=9){

        // Initialize data
        *(FdotFirstTerm)   = 0.0 ; *(FdotFirstTerm+1) = 0.0 ;  *(FdotFirstTerm+2) = 0.0;
        *(FdotFirstTerm+3) = 0.0 ; *(FdotFirstTerm+4) = 0.0 ;  *(FdotFirstTerm+5) = 0.0;
        *(FdotFirstTerm+6) = 0.0 ; *(FdotFirstTerm+7) = 0.0 ;  *(FdotFirstTerm+8) = 0.0;
        
        //Compute Fdot
        numNeighbors = *neighborListPtr; neighborListPtr++;
        
        for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamage++){
          neighborIndex = *neighborListPtr;    
          neighborModelCoord = modelCoordinates + 3*neighborIndex;
          neighborVel = velocities + 3*neighborIndex;
          
         if (detachedNodes[iID]!=0) continue;
         if (detachedNodes[neighborIndex]!=0) continue;
          
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
        
        // Compute Fdot
        MATRICES::MatrixMultiply(false, false, One, FdotFirstTerm, shapeTensorInv, Fdot);

        // Compute the inverse of the deformation gradient, Finverse
        if (type==true){
            inversionReturnCode = MATRICES::Invert2by2Matrix(defGrad, determinant, Finverse);
        }
        else{
            inversionReturnCode = MATRICES::Invert3by3Matrix(defGrad, determinant, Finverse);
            }

        std::string matrixInversionErrorMessage =
          "**** Error:  Correspondence ::getLinearUnrotatedRateOfDeformation() failed to invert deformation gradient.\n";
        TEUCHOS_TEST_FOR_TERMINATION(inversionReturnCode != 0, matrixInversionErrorMessage);

        // Compute the Eulerian velocity gradient L = Fdot * Finv
        MATRICES::MatrixMultiply(false, false, One, Fdot, Finverse, eulerianVelGrad);

        
        
        
        
        
        // Versuch der nicht geklappt hat
        //if (type==true){
        //    inversionReturnCode = MATRICES::Invert2by2Matrix(defGradNP1, determinant, Finverse);
        //}
        //else{
        //    inversionReturnCode = MATRICES::Invert3by3Matrix(defGradNP1, determinant, Finverse);
        //}
        //for(int i = 0; i < 9; i++){
        //  *(Fdot+i) = (*(defGradNP1+i)-*(defGrad+i))/ dt;
        //}
        //MatrixMultiply(false, false, One, Fdot, Finverse, eulerianVelGrad);

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
        
     }

}





template<typename ScalarT>
int computeLogStrain
(
 const ScalarT* defGrad,
 ScalarT* strain
)
{
  // Hencky-Strain E = 0.5*ln(C)
  // C = F^T*F

  int returnCode(0);

  std::vector<ScalarT> defGradTempVector(9);
  ScalarT* defGradTemp = &defGradTempVector[0];
  std::vector<ScalarT> VVector(9);
  ScalarT* V = &VVector[0];
  std::vector<ScalarT> V_invVector(9);
  ScalarT* V_inv = &V_invVector[0];
  std::vector<ScalarT> A1_tempVector(9);
  ScalarT* A1_temp = &A1_tempVector[0];
  std::vector<ScalarT> A1Vector(9);
  ScalarT* A1 = &A1Vector[0];
  std::vector<ScalarT> A1_D_logVector(9);
  ScalarT* A1_D_log = &A1_D_logVector[0];
  std::vector<ScalarT> A_log_tempVector(9);
  ScalarT* A_log_temp = &A_log_tempVector[0];
  std::vector<ScalarT> strainTempVector(9);
  ScalarT* strainTemp = &strainTempVector[0];
  ScalarT determinant;
  int eigenVecReturnCode(1);
  int inversionReturnCode(1);
  ScalarT One = 1.0;
  ScalarT OneHalf = 0.5;
  
  MATRICES::MatrixMultiply(true,false,One,defGrad, defGrad, defGradTemp);
  
  eigenVecReturnCode = EigenVec2D(defGradTemp,V);
  
  if(eigenVecReturnCode !=0)
  {
    returnCode=1;
    return returnCode;
  }

  inversionReturnCode = MATRICES::Invert2by2Matrix(V,determinant,V_inv);

  if(inversionReturnCode !=0)
  {
    returnCode=2;
    return returnCode;
  }

  MATRICES::MatrixMultiply(false,false,One,V_inv, defGradTemp, A1_temp);

  MATRICES::MatrixMultiply(false,false,One,A1_temp, V, A1);
  
  *(A1_D_log+0) = log(*(A1+0));
  *(A1_D_log+1) = *(A1+1);
  *(A1_D_log+2) = 0.0 ;
  *(A1_D_log+3) = *(A1+3);
  *(A1_D_log+4) = log(*(A1+4));
  *(A1_D_log+5) = 0.0 ;
  *(A1_D_log+6) = 0.0 ;
  *(A1_D_log+7) = 0.0 ;
  *(A1_D_log+8) = 0.0 ;
  
  MATRICES::MatrixMultiply(false,false,One,V, A1_D_log, A_log_temp);
  
  for (int i = 0; i < 9; i++)
  {
    if(*(A_log_temp+i)!=*(A_log_temp+i) || *(V_inv+i)!=*(V_inv+i))
    {
      returnCode=3;
      return returnCode;
    }
  }
  
  MATRICES::MatrixMultiply(false,false,OneHalf,A_log_temp, V_inv, strainTemp);

  *(strain)   = *(strainTemp+0);
  *(strain+1) = *(strainTemp+1);
  *(strain+3) = *(strainTemp+3);
  *(strain+4) = *(strainTemp+4);

  return returnCode;
}


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
)
{
  //Using RK Implicit Gradient (Same as PD Gradient Operator) to construct
  //these shapes on the parametric space (2D gradients).
  //GMLS or other methods can be incorporated later.

  int returnCode = 0;

  const double* delta = horizon;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;

  ScalarT* phi1 = gradientWeight1;
  ScalarT* phi2 = gradientWeight2;
  ScalarT* phi3 = gradientWeight3;

  const double* flyingPointFlg = flyingPointFlag;
  const double *bondDamagePtr = bondDamage;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;

  double neighborVolume, omega, temp;

  int Qdim;

  int counter, thisOrder, p1, p2, p3;
  int i, j;

  // calculate dimension of Q vector
  Qdim = 0;
  for(thisOrder=1; thisOrder<=accuracyOrder; thisOrder++){
    for(p1=thisOrder; p1>=0; p1--){ // x-power
      for(p2=thisOrder-p1; p2>=0; p2--){ // y-power
        p3=thisOrder-p1-p2; //z-power
        Qdim++;
      }
    }
  }

  std::vector<ScalarT> QVector(Qdim);
  ScalarT* Q = &QVector[0];

  double phi1const = 1.0; int phi1ind = 0;
  double phi2const = 1.0; int phi2ind = 1;
  double phi3const = 1.0; int phi3ind = 2;

  std::vector<ScalarT> MVector(Qdim*Qdim);
  ScalarT* M = &MVector[0];

  std::vector<ScalarT> MinvVector(Qdim*Qdim);
  ScalarT* Minv = &MinvVector[0];

  double thresVal;

  int inversionReturnCode(0);

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3, flyingPointFlg++){

    if(*flyingPointFlg < 0.5){

      // Zero out data
      for(j=0; j<Qdim; j++)
        for(i=0; i<Qdim; i++)
          *(M+Qdim*i+j) = 0.0;

      // Calculate Moment matrix
      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++){

        neighborIndex = *neighborListPtr;
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVolume = jacobianDeterminant[neighborIndex] * volume[neighborIndex];

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);

        // Calculate Q for this bond
        counter = 0;
        for(thisOrder=1; thisOrder<=accuracyOrder; thisOrder++){
          for(p1=thisOrder; p1>=0; p1--){ // x-power
            for(p2=thisOrder-p1; p2>=0; p2--){ // y-power
              p3=thisOrder-p1-p2; //z-power

              Q[counter] = 1.0;
              for(i=0; i<p1; i++)
                Q[counter] *= deformedBondX / *delta;
              for(i=0; i<p2; i++)
                Q[counter] *= deformedBondY / *delta;
              for(i=0; i<p3; i++)
                Q[counter] *= deformedBondZ / *delta;

              counter++;
            }
          }
        }

        omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);
        omega *= (1.0 - *bondDamagePtr);

        temp = omega * neighborVolume;
        for(j=0; j<Qdim; j++)
          for(i=0; i<Qdim; i++)
            *(M+Qdim*i+j) += temp * Q[i] * Q[j];
      }

      thresVal = 1.0e-6 * *delta * *delta * *delta; // should be area (2d) or vol (3d) multiplied by a threshold value

      // calculate the inverse of the moment matrix (this must be a full rank matrix)
      //inversionReturnCode = computeSymmetrixMatrixInverse(M, Qdim, Minv);
      inversionReturnCode = MATRICES::invertAndCond(M, Minv, Qdim, thresVal);

      if(inversionReturnCode > 0)
        returnCode = inversionReturnCode;

      // Re-iterate over the neighbor set and compute Phis
      // Return the neighbor pointers to the beginning of set
      neighborListPtr -= numNeighbors; 
      bondDamagePtr -= numNeighbors; 
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, phi1++, phi2++, phi3++, bondDamagePtr++){

        neighborIndex = *neighborListPtr;
        neighborCoord = coordinates + 3*neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);

        // Calculate Q for this bond
        counter = 0;
        for(thisOrder=1; thisOrder<=accuracyOrder; thisOrder++){
          for(p1=thisOrder; p1>=0; p1--){ // x-power
            for(p2=thisOrder-p1; p2>=0; p2--){ // y-power
              p3=thisOrder-p1-p2; //z-power

              Q[counter] = 1.0;
              for(i=0; i<p1; i++)
                Q[counter] *= deformedBondX / *delta;
              for(i=0; i<p2; i++)
                Q[counter] *= deformedBondY / *delta;
              for(i=0; i<p3; i++)
                Q[counter] *= deformedBondZ / *delta;

              counter++;
            }
          }
        }

        omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);
        omega *= (1.0 - *bondDamagePtr);

        // caluclate phi1
        *phi1 = 0.0;
        temp = phi1const * omega;  
        for(j=0; j<Qdim; j++)
          *phi1 += temp * *(Minv+Qdim*phi1ind+j) * Q[j] / *delta;

        // caluclate phi2
        *phi2 = 0.0;
        temp = phi2const * omega;  
        for(j=0; j<Qdim; j++)
          *phi2 += temp * *(Minv+Qdim*phi2ind+j) * Q[j] / *delta;

        // caluclate phi3
        *phi3 = 0.0;
        temp = phi3const * omega;  
        for(j=0; j<Qdim; j++)
          *phi3 += temp * *(Minv+Qdim*phi3ind+j) * Q[j] / *delta;
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      phi1 += numNeighbors; phi2 += numNeighbors; phi3 += numNeighbors; 
      bondDamagePtr += numNeighbors; 
    }
  }

  return returnCode;
}


template int computeShapeTensorInverseAndApproximateDeformationGradient<Sacado::Fad::DFad<double> >
(
 const double* volume,
 const double* horizon,
 const double* modelCoordinates,
 const Sacado::Fad::DFad<double>* coordinatesNP1,
 Sacado::Fad::DFad<double>* shapeTensorInverse,
 Sacado::Fad::DFad<double>* deformationGradient,
 const double* bondDamageNP1,
 const int* neighborhoodList,
 int numPoints,
 const bool type,
 double* detachedNodes

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
)
{
  int returnCode = 0;

  const double* delta = horizon;
  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinatesNP1;
  //const ScalarT* neighborCoord;
  const ScalarT* neighborCoordNP1;
  ScalarT* shapeTensorInv = shapeTensorInverse;
  ScalarT* defGrad = deformationGradient;
  
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  ScalarT deformedBondX, deformedBondY, deformedBondZ;
  double neighborVolume, omega, temp;
  ScalarT One = 1.0;
  std::vector<ScalarT> shapeTensorVector(9);
  ScalarT* shapeTensor = &shapeTensorVector[0];
  ScalarT shapeTensorDeterminant;

  std::vector<ScalarT> defGradFirstTermVector(9);
  ScalarT* defGradFirstTerm = &defGradFirstTermVector[0];



  int inversionReturnCode(0);

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, coord+=3,
        shapeTensorInv+=9, defGrad+=9){
  
    //double bondCheck(0.0), bondCheckNP1(0.0);
    *(shapeTensor)   = 0.0 ; *(shapeTensor+1) = 0.0 ; *(shapeTensor+2) = 0.0 ;
    *(shapeTensor+3) = 0.0 ; *(shapeTensor+4) = 0.0 ; *(shapeTensor+5) = 0.0 ;
    *(shapeTensor+6) = 0.0 ; *(shapeTensor+7) = 0.0 ; *(shapeTensor+8) = 0.0 ;
    *(defGradFirstTerm)   = 0.0 ; *(defGradFirstTerm+1) = 0.0 ; *(defGradFirstTerm+2) = 0.0 ;
    *(defGradFirstTerm+3) = 0.0 ; *(defGradFirstTerm+4) = 0.0 ; *(defGradFirstTerm+5) = 0.0 ;
    *(defGradFirstTerm+6) = 0.0 ; *(defGradFirstTerm+7) = 0.0 ; *(defGradFirstTerm+8) = 0.0 ;

    numNeighbors = *neighborListPtr; neighborListPtr++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamageNP1++){
      
      neighborIndex = *neighborListPtr;

      neighborVolume = volume[neighborIndex];
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      //neighborCoord      = coordinates      + 3*neighborIndex;
      neighborCoordNP1   = coordinatesNP1   + 3*neighborIndex;
      if (*(detachedNodes+iID)!=0) continue;
      if (*(detachedNodes+neighborIndex)!=0) continue;

      undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
      
      // its increment independent to avoid problems with the influence function
      // currently the horizon does not realy deform
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);
      
      deformedBondX = *(neighborCoordNP1)   - *(coord);
      deformedBondY = *(neighborCoordNP1+1) - *(coord+1);
      deformedBondZ = *(neighborCoordNP1+2) - *(coord+2);
      
      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta);


      temp = (1.0 - *bondDamageNP1) * omega * neighborVolume;
      
      *(shapeTensor)   += temp * undeformedBondX * undeformedBondX;
      *(shapeTensor+1) += temp * undeformedBondX * undeformedBondY;
      
      *(shapeTensor+2) += temp * undeformedBondX * undeformedBondZ;
      
      *(shapeTensor+3) += temp * undeformedBondY * undeformedBondX;
      *(shapeTensor+4) += temp * undeformedBondY * undeformedBondY;
      
      *(shapeTensor+5) += temp * undeformedBondY * undeformedBondZ;
      *(shapeTensor+6) += temp * undeformedBondZ * undeformedBondX;
      *(shapeTensor+7) += temp * undeformedBondZ * undeformedBondY;
      *(shapeTensor+8) += temp * undeformedBondZ * undeformedBondZ;
      
      //if (deformedBondX-undeformedBondX != 0 and deformedBondY-undeformedBondY != 0)std::cout<<"es passiert was"<<std::endl;
          
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
                inversionReturnCode = MATRICES::Invert2by2Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);

            }
            else{
                inversionReturnCode = MATRICES::Invert3by3Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);
            }
            
        
        MATRICES::MatrixMultiply(false, false, One, defGradFirstTerm, shapeTensorInv, defGrad);
        }
        
        
        if(*(shapeTensor) == 0 or inversionReturnCode > 0){
           
          returnCode = inversionReturnCode;
          *(shapeTensorInv)   = 0.0 ; *(shapeTensorInv+1) = 0.0 ; *(shapeTensorInv+2) = 0.0 ;
          *(shapeTensorInv+3) = 0.0 ; *(shapeTensorInv+4) = 0.0 ; *(shapeTensorInv+5) = 0.0 ;
          *(shapeTensorInv+6) = 0.0 ; *(shapeTensorInv+7) = 0.0 ; *(shapeTensorInv+8) = 0.0 ;
          *(defGrad)   = 1.0 ;     *(defGrad+1) = 0.0 ;     *(defGrad+2) = 0.0 ;
          *(defGrad+3) = 0.0 ;     *(defGrad+4) = 1.0 ;     *(defGrad+5) = 0.0 ;
          *(defGrad+6) = 0.0 ;     *(defGrad+7) = 0.0 ;     *(defGrad+8) = 1.0 ;
        }

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
const bool type,
double* detachedNodes
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
  ScalarT One = 1.0;
  ScalarT omegaX, omegaY, omegaZ;
  ScalarT zX, zY, zZ;
  ScalarT wX, wY, wZ;
  ScalarT velStateX, velStateY, velStateZ;
  ScalarT traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double neighborVolume, omega, scalarTemp; 
  int inversionReturnCode(0);
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
    
    // Compute Fdot
    MATRICES::MatrixMultiply(false, false, One, FdotFirstTerm, shapeTensorInv, Fdot);

    // Compute the inverse of the deformation gradient, Finverse
    if (type==true){
        inversionReturnCode = MATRICES::Invert2by2Matrix(defGrad, determinant, Finverse);
    }
    else{
        inversionReturnCode = MATRICES::Invert3by3Matrix(defGrad, determinant, Finverse);
        }
    if(inversionReturnCode > 0){
        returnCode = 2;
        return returnCode;
    }
    // Compute the Eulerian velocity gradient L = Fdot * Finv
    MATRICES::MatrixMultiply(false, false, One, Fdot, Finverse, eulerianVelGrad);

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
        inversionReturnCode = MATRICES::Invert2by2Matrix(temp, determinant, tempInv);
    }
    else{
        inversionReturnCode = MATRICES::Invert3by3Matrix(temp, determinant, tempInv);
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
    //for(int i=0 ; i<9 ; ++i)
    //  *(unrotRateOfDef+i) = *(unrotRateOfDef+i) * dt;
  }

  return returnCode;
}

void getOrientations
  (
    const int numOwnedPoints,
    const double* angles,
    double* orientations
    )
{
  const double *pointAnglePtr;
  double *pointOrientationPtr;
  const double* pointAngles = angles;
  double* pointOrientation = orientations;

  std::vector<double> rotMatVector(9);
  double* rotMat = &rotMatVector[0];

  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    pointAnglePtr = pointAngles + 3*iID;
    pointOrientationPtr = pointOrientation + 3*iID;

    double alpha[3];

    alpha[0] = *(pointAnglePtr);
    alpha[1] = *(pointAnglePtr+1);
    alpha[2] = *(pointAnglePtr+2);

    MATRICES::createRotationMatrix(alpha, rotMat);

    *(pointOrientationPtr) = rotMat[0];
    *(pointOrientationPtr+1) = rotMat[3];
    *(pointOrientationPtr+2) = rotMat[6];
    
  }
}

template void computeForcesAndStresses<Sacado::Fad::DFad<double> >
  (
    const int numOwnedPoints,
    const int* neighborhoodList,
    const double* volume,
    const double* horizon,
    const double* modelCoordinates,
    const Sacado::Fad::DFad<double>* coordinatesNP1,
    const Sacado::Fad::DFad<double>* deformationGradient,
    const Sacado::Fad::DFad<double>* cauchyStressNP1,
    const Sacado::Fad::DFad<double>* cauchyStressPlastic,
    const Sacado::Fad::DFad<double>* shapeTensorInverse,
    const double* bondDamage,
    const Sacado::Fad::DFad<double> C[][6],
    const double* angles,
    Sacado::Fad::DFad<double>* force,
    Sacado::Fad::DFad<double>* partialStress,
    Sacado::Fad::DFad<double>* tempStressX,
    Sacado::Fad::DFad<double>* tempStressY,
    Sacado::Fad::DFad<double>* tempStressZ,
    Sacado::Fad::DFad<double>* hourglassStiff,
    const double m_hourglassCoefficient,
    const int m_stabilizationType,
    const bool m_plane,
    const bool plast,
    const bool adaptHourGlass,
    double* detachedNodes
);

template void computeForcesAndStresses<double>
  (
    const int numOwnedPoints,
    const int* neighborhoodList,
    const double* volume,
    const double* horizon,
    const double* modelCoordinates,
    const double* coordinatesNP1,
    const double* deformationGradient,
    const double* cauchyStressNP1,
    const double* cauchyStressPlastic,
    const double* shapeTensorInverse,
    const double* bondDamage,
    const double C[][6],
    const double* angles,
    double* force,
    double* partialStress,
    double* tempStressX,
    double* tempStressY,
    double* tempStressZ,
    double* hourglassStiff,
    const double m_hourglassCoefficient,
    const int m_stabilizationType,
    const bool m_plane,
    const bool plast,
    const bool adaptHourGlass,
    double* detachedNodes
);


template<typename ScalarT>
void computeForcesAndStresses
  (
    const int numOwnedPoints,
    const int* neighborhoodList,
    const double* volume,
    const double* horizon,
    const double* modelCoor,
    const ScalarT* coordinatesNP1,
    const ScalarT* deformationGradient,
    const ScalarT* cauchyStressNP1,
    const ScalarT* cauchyStressPlastic,
    const ScalarT* shapeTensorInverse,
    const double* bondDamage,
    const ScalarT C[][6],
    const double* angles,
    ScalarT* force,
    ScalarT* partialStressValues,
    ScalarT* StressX,
    ScalarT* StressY,
    ScalarT* StressZ,
    ScalarT* hourglassStiffValues,
    const double m_hourglassCoefficient,
    const int m_stabilizationType,
    const bool m_plane,
    const bool plast,
    const bool adaptHourGlass,
    double* detachedNodes
    )
{
    
  const double* delta = horizon;
  const double* modelCoordinates = modelCoor;
  const ScalarT* coorNP1 = coordinatesNP1;
  ScalarT* forceDensity = force;
  ScalarT* partialStress = partialStressValues;
  ScalarT* tempStressX = StressX;
  ScalarT* tempStressY = StressY;
  ScalarT* tempStressZ = StressZ;
  ScalarT* hourglassStiff = hourglassStiffValues;
  double hourglassScale = m_hourglassCoefficient;
  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  const double *modelCoordinatesPtr, *neighborModelCoordinatesPtr;
  ScalarT *partialStressPtr;
  const ScalarT* stress = cauchyStressNP1;
  const ScalarT* plasticStress = cauchyStressPlastic;
  const ScalarT* shapeTensorInv = shapeTensorInverse;
  const ScalarT* defGrad = deformationGradient;
  const ScalarT *deformedCoordinatesNP1Ptr, *neighborDeformedCoordinatesNP1Ptr;
  double X_dx, X_dy, X_dz, undeformedBondLength;
  ScalarT TX, TY, TZ;
  ScalarT *forceDensityPtr;
  double omega, vol, neighborVol;
  const double *pointAnglePtr;
  const double* pointAngles = angles;
  
  
  std::string matrixInversionErrorMessage =
    "**** Error:  CorrespondenceMaterialconst ::computeForce() failed to invert deformation gradient.\n";
  matrixInversionErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  matrixInversionErrorMessage +=
    "****         Force calculation";

  ScalarT jacobianDeterminant;

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities

  int bondIndex = 0;
  int matrixInversionReturnCode(0);
  
  std::vector<ScalarT> piolaStressVector(9), tempVector(9), tempPlastVector(9), defGradInvVector(9), hourglassStiffVector(9), 
  TSvector(3), scalVector(9);
  ScalarT* TS = &TSvector[0];
  ScalarT* temp = &tempVector[0];
  ScalarT* tempPlast = &tempPlastVector[0];
  ScalarT* defGradInv = &defGradInvVector[0];
  ScalarT* piolaStress = &piolaStressVector[0];

  ScalarT* hourglassStiffVal = &hourglassStiffVector[0];
  ScalarT* scal = &scalVector[0]; 
  ScalarT  One = 1.0;
  for(int i=0; i<9; i++){
    scal[i] = 1;
    }

  for(int iID=0 ; iID<numOwnedPoints ; ++iID, 
          ++delta, defGrad+=9, stress+=9, plasticStress+=9, shapeTensorInv+=9, hourglassStiff+=9){ //, defGradNonInc+=9
           
    // first Piola-Kirchhoff stress = J * cauchyStress * defGrad^-T

    // Invert the deformation gradient and store the determinant
   
    // Loop over the neighbors and compute contribution to force densities
    modelCoordinatesPtr = modelCoordinates + 3*iID;
    pointAnglePtr       = pointAngles      + 3*iID;
    //deformedCoordinatesPtr    = coordinates   + 3*iID;
    deformedCoordinatesNP1Ptr = coorNP1 + 3*iID;
    numNeighbors = *neighborListPtr; neighborListPtr++;
    forceDensityPtr = forceDensity + 3*iID;
    if (detachedNodes[iID]==0){
        if (m_plane==true){
            matrixInversionReturnCode =
            MATRICES::Invert2by2Matrix(defGrad, jacobianDeterminant, defGradInv);
            }
        else{
            matrixInversionReturnCode =
             MATRICES::Invert3by3Matrix(defGrad, jacobianDeterminant, defGradInv);
            }
    }
    TEUCHOS_TEST_FOR_TERMINATION(matrixInversionReturnCode != 0, matrixInversionErrorMessage);
    
    //P = J * \sigma * F^(-T)
     MATRICES::MatrixMultiply(false, true, jacobianDeterminant, stress, defGradInv, piolaStress);

    // Inner product of Piola stress and the inverse of the shape tensor
     MATRICES::MatrixMultiply(false, false, One, piolaStress, shapeTensorInv, temp);

    if (plast){
        
         MATRICES::MatrixMultiply(false, true, jacobianDeterminant, plasticStress, defGradInv, piolaStress);
         MATRICES::MatrixMultiply(false, false, One, piolaStress, shapeTensorInv, tempPlast);
        if (adaptHourGlass){// kann wahrscheinlich weg
            hourglassScale = 1;
            for(int i=0; i<9; i++){
                if (*(stress+i) != *(plasticStress+i)) {//scal[i] = 1 - *(plasticStress+i)/(*(stress+i));
                //scal[i] =0;
                hourglassScale = 0;
                break;}
                //else {scal[i]=1;}
            }
        }
    }
    if (matrixInversionReturnCode != 0){
        *(temp)   = 0;*(temp+1) = 0;*(temp+2) = 0;
        *(temp+3) = 0;*(temp+4) = 0;*(temp+5) = 0;
        *(temp+6) = 0;*(temp+7) = 0;*(temp+8) = 0;
        // as a last resort. might stabilize sometimes
        detachedNodes[iID] = 1;
    }

    if (m_stabilizationType == 3){
        double alpha[3];
        alpha[0] = *(pointAnglePtr);
        alpha[1] = *(pointAnglePtr+1);
        alpha[2] = *(pointAnglePtr+2);
        
        CORRESPONDENCE::createHourglassStiffness(C, alpha, shapeTensorInv, hourglassStiffVal);
        // scale the stiffness dependend to the plastic / elastic stress relation. If no elastic stresses are there, no additional stiffness is needed, because the bond is gone anyway

        *(hourglassStiff  ) = scal[0]*hourglassStiffVal[0]; *(hourglassStiff+1) = scal[1]*hourglassStiffVal[1]; *(hourglassStiff+2) = scal[2]*hourglassStiffVal[2];
        *(hourglassStiff+3) = scal[3]*hourglassStiffVal[3]; *(hourglassStiff+4) = scal[4]*hourglassStiffVal[4]; *(hourglassStiff+5) = scal[5]*hourglassStiffVal[5];
        *(hourglassStiff+6) = scal[6]*hourglassStiffVal[6]; *(hourglassStiff+7) = scal[7]*hourglassStiffVal[7]; *(hourglassStiff+8) = scal[8]*hourglassStiffVal[8];
        // checken ob Null
    }
    

    //int countNeighbors = 0;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){

      neighborIndex = *neighborListPtr;
      neighborModelCoordinatesPtr = modelCoordinates + 3*neighborIndex;
      neighborDeformedCoordinatesNP1Ptr = coorNP1 + 3*neighborIndex;

      if (detachedNodes[iID]==0 && detachedNodes[neighborIndex]==0){
          X_dx = *(neighborModelCoordinatesPtr)   - *(modelCoordinatesPtr);
          X_dy = *(neighborModelCoordinatesPtr+1) - *(modelCoordinatesPtr+1);
          X_dz = *(neighborModelCoordinatesPtr+2) - *(modelCoordinatesPtr+2);

          if (m_plane == true) {X_dz = 0;}
          
          undeformedBondLength = sqrt(X_dx*X_dx +
                                      X_dy*X_dy +
                                      X_dz*X_dz);
        
          omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta);

          if (m_stabilizationType == 3){
              //double hourglassStiff[9];
              ScalarT FxsiX, FxsiY, FxsiZ;

              ScalarT Y_dx = *(neighborDeformedCoordinatesNP1Ptr)   - *(deformedCoordinatesNP1Ptr);
              ScalarT Y_dy = *(neighborDeformedCoordinatesNP1Ptr+1) - *(deformedCoordinatesNP1Ptr+1);
              ScalarT Y_dz = *(neighborDeformedCoordinatesNP1Ptr+2) - *(deformedCoordinatesNP1Ptr+2);

              FxsiX = *(defGrad)   * X_dx + *(defGrad+1) * X_dy + *(defGrad+2) * X_dz;
              FxsiY = *(defGrad+3) * X_dx + *(defGrad+4) * X_dy + *(defGrad+5) * X_dz;
              FxsiZ = *(defGrad+6) * X_dx + *(defGrad+7) * X_dy + *(defGrad+8) * X_dz;
           
              
              CORRESPONDENCE::computeCorrespondenceStabilityWanEtAlShort(FxsiX,FxsiY,FxsiZ,Y_dx,Y_dy,Y_dz,hourglassStiffVal,TS);
              
              
          }
            
          TX =  (1-bondDamage[bondIndex]) * omega * ( *(temp)   * X_dx + *(temp+1) * X_dy + *(temp+2) * X_dz+ hourglassScale*TS[0]);
          TY =  (1-bondDamage[bondIndex]) * omega * ( *(temp+3) * X_dx + *(temp+4) * X_dy + *(temp+5) * X_dz+ hourglassScale*TS[1]);
          TZ =  (1-bondDamage[bondIndex]) * omega * ( *(temp+6) * X_dx + *(temp+7) * X_dy + *(temp+8) * X_dz+ hourglassScale*TS[2]);
          
         
          neighborVol = volume[neighborIndex];
          vol = volume[iID];

          MATERIAL_EVALUATION::setForces(TX, TY, TZ, vol, neighborVol, forceDensityPtr, &force[3*neighborIndex]);
          partialStressPtr = partialStress + 9*iID;
          MATERIAL_EVALUATION::setPartialStresses(TX, TY, TZ, X_dx, X_dy, X_dz, neighborVol, partialStressPtr);
      }

      //  countNeighbors += bondDamage[bondIndex];
      
      bondIndex += 1;

    }

    if (detachedNodes[iID]==0){
        if (plast){
            
            tempStressX[3*iID  ] = *(tempPlast);
            tempStressX[3*iID+1] = *(tempPlast+1);
            tempStressX[3*iID+2] = *(tempPlast+2);
            tempStressY[3*iID  ] = *(tempPlast+3);
            tempStressY[3*iID+1] = *(tempPlast+4);
            tempStressY[3*iID+2] = *(tempPlast+5);
            tempStressZ[3*iID  ] = *(tempPlast+6);
            tempStressZ[3*iID+1] = *(tempPlast+7);
            tempStressZ[3*iID+2] = *(tempPlast+8);
        }
        else {
            tempStressX[3*iID  ] = *(temp);
            tempStressX[3*iID+1] = *(temp+1);
            tempStressX[3*iID+2] = *(temp+2);
            tempStressY[3*iID  ] = *(temp+3);
            tempStressY[3*iID+1] = *(temp+4);
            tempStressY[3*iID+2] = *(temp+5);
            tempStressZ[3*iID  ] = *(temp+6);
            tempStressZ[3*iID+1] = *(temp+7);
            tempStressZ[3*iID+2] = *(temp+8);
        }
    }
    
    
    
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

  const double pi = PeridigmNS::value_of_pi();
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
 // particle discretizations" in Comput. Methods Appl. Mech. Engrg. 322 (2017) 42�57,http://dx.doi.org/10.1016/j.cma.2017.03.043

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

  const double pi = PeridigmNS::value_of_pi();
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
const double C[][6],
const double* bondDamage,
ScalarT* hourglassForceDensity,
double hourglassCoefficient//,
//double* hourglassStiff
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
  ScalarT hourglassStiff[9];
  ScalarT TSx, TSy, TSz;
  ScalarT FxsiX, FxsiY, FxsiZ;
  const ScalarT* neighborCoord;
  double vol, neighborVol, factor;
  ScalarT* hourglassForceDensityPtr = hourglassForceDensity;
  ScalarT* neighborHourglassForceDensityPtr;

  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, coord+=3,
     defGrad+=9, shapeTensorInv+=9, hourglassForceDensityPtr+=3){ //, hourglassStiff+=9){
    //constant = firstPartOfConstant/( omega0 * *(delta2)* *(delta2)* *(delta2)* *(delta2)* *(delta2) );
    numNeighbors = *neighborListPtr; neighborListPtr++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamage++){
      
      neighborIndex = *neighborListPtr;
      
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborCoord = coordinates + 3*neighborIndex;
                 
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
      //2d to be checked!!
      hourglassStiff[0] = C[0][0] * *(shapeTensorInv)  + C[0][1] * *(shapeTensorInv+4)  + C[0][2] * *(shapeTensorInv+8)  + C[0][3] *(*(shapeTensorInv+5) + *(shapeTensorInv+7) ) + C[0][4] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[0][5] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[1] = C[5][0]* *(shapeTensorInv)  + C[5][1]* *(shapeTensorInv+4)  + C[5][2]* *(shapeTensorInv+8)  + C[5][3]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[5][4] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[5][5]*(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[2] = C[4][0]* *(shapeTensorInv)  + C[4][1]* *(shapeTensorInv+4)  + C[4][2]* *(shapeTensorInv+8)  + C[4][3]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[4][4] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[4][5]*(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[3] = hourglassStiff[1];
      hourglassStiff[4] = C[1][0] * *(shapeTensorInv)  + C[1][1] * *(shapeTensorInv+4)  + C[1][2] * *(shapeTensorInv+8)  + C[1][3] *(*(shapeTensorInv+5) + *(shapeTensorInv+7) ) + C[1][4]*(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[1][5] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[5] = C[3][0]* *(shapeTensorInv)  + C[3][1]* *(shapeTensorInv+4)  + C[3][2]* *(shapeTensorInv+8)  + C[3][3]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[3][4] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[3][5] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      hourglassStiff[6] = hourglassStiff[2];
      hourglassStiff[7] = hourglassStiff[5];
      hourglassStiff[8] = C[2][0]* *(shapeTensorInv)  + C[2][1]* *(shapeTensorInv+4)  + C[2][2]* *(shapeTensorInv+8)  + C[2][3]*(*(shapeTensorInv+5) + *(shapeTensorInv+7) )+ C[2][4] *(*(shapeTensorInv+2) + *(shapeTensorInv+6) ) + C[2][5] *(*(shapeTensorInv+1) + *(shapeTensorInv+3) );
      
      factor = hourglassCoefficient*(1.0-*bondDamage) * m_omega;
      
      TSx = factor*(hourglassStiff[0] * nonUniformDeformState[0] + hourglassStiff[1] * nonUniformDeformState[1] + hourglassStiff[2] * nonUniformDeformState[2]);
      TSy = factor*(hourglassStiff[3] * nonUniformDeformState[0] + hourglassStiff[4] * nonUniformDeformState[1] + hourglassStiff[5] * nonUniformDeformState[2]);
      TSz = factor*(hourglassStiff[6] * nonUniformDeformState[0] + hourglassStiff[7] * nonUniformDeformState[1] + hourglassStiff[8] * nonUniformDeformState[2]);
      //TSx = factor*(*(hourglassStiff)   * nonUniformDeformState[0] + *(hourglassStiff+1) * nonUniformDeformState[1] + *(hourglassStiff+2) * nonUniformDeformState[2]);
      //TSy = factor*(*(hourglassStiff+3) * nonUniformDeformState[0] + *(hourglassStiff+4) * nonUniformDeformState[1] + *(hourglassStiff+5) * nonUniformDeformState[2]);
      //TSz = factor*(*(hourglassStiff+6) * nonUniformDeformState[0] + *(hourglassStiff+7) * nonUniformDeformState[1] + *(hourglassStiff+8) * nonUniformDeformState[2]);

      
      
      vol = volume[iID];
      neighborVol = volume[neighborIndex];
      neighborHourglassForceDensityPtr = hourglassForceDensity + 3*neighborIndex;

      MATERIAL_EVALUATION::setForces(TSx, TSy, TSz, vol, neighborVol, hourglassForceDensityPtr, neighborHourglassForceDensityPtr);


    }
  }
  
}

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
)
{
  // J. Wan et al., "Improved method for zero-energy mode suppression in peridynamic correspondence model" 
  // in Acta Mechanica Sinica https://doi.org/10.1007/s10409-019-00873-y
    ScalarT nonUniformDeformState[3];

    nonUniformDeformState[0] = deformedBondX - FxsiX;
    nonUniformDeformState[1] = deformedBondY - FxsiY;
    nonUniformDeformState[2] = deformedBondZ - FxsiZ;
    TS[0] = (hourglassStiff[0] * nonUniformDeformState[0] + hourglassStiff[1] * nonUniformDeformState[1] + hourglassStiff[2] * nonUniformDeformState[2]);
    TS[1] = (hourglassStiff[3] * nonUniformDeformState[0] + hourglassStiff[4] * nonUniformDeformState[1] + hourglassStiff[5] * nonUniformDeformState[2]);
    TS[2] = (hourglassStiff[6] * nonUniformDeformState[0] + hourglassStiff[7] * nonUniformDeformState[1] + hourglassStiff[8] * nonUniformDeformState[2]);
 
}

template<typename ScalarT>
void createHourglassStiffness
(
const ScalarT Cstiff[][6],
const double alpha[],
const ScalarT* shapeTensorInverse,
ScalarT* hourglassStiff
)
{
    const ScalarT* shapeTensorInv = shapeTensorInverse;
    ScalarT C[6][6];
    std::vector<ScalarT> rotMatVector(9);
    ScalarT* rotMat = &rotMatVector[0];
    MATRICES::createRotationMatrix(alpha,rotMat);
    CORRESPONDENCE::createRotatedStiff(Cstiff,rotMat,C);
    //CORRESPONDENCE::createRotatedPythonBasedStiff(Cstiff,alpha,C); 
    // to be checked for 2D
    hourglassStiff[0] =  *(shapeTensorInv)*C[0][0] + *(shapeTensorInv+1)*C[0][5] + *(shapeTensorInv+2)*C[0][4] + *(shapeTensorInv+3)*C[0][5] + *(shapeTensorInv+4)*C[0][1] + *(shapeTensorInv+5)*C[0][3] + *(shapeTensorInv+6)*C[0][4] + *(shapeTensorInv+7)*C[0][3] + *(shapeTensorInv+8)*C[0][2];
    hourglassStiff[1] =  *(shapeTensorInv)*C[0][5] + *(shapeTensorInv+1)*C[5][5] + *(shapeTensorInv+2)*C[3][5] + *(shapeTensorInv+3)*C[5][5] + *(shapeTensorInv+4)*C[1][5] + *(shapeTensorInv+5)*C[3][4] + *(shapeTensorInv+6)*C[3][5] + *(shapeTensorInv+7)*C[3][4] + *(shapeTensorInv+8)*C[2][3];
    hourglassStiff[2] =  *(shapeTensorInv)*C[0][4] + *(shapeTensorInv+1)*C[3][5] + *(shapeTensorInv+2)*C[4][4] + *(shapeTensorInv+3)*C[3][5] + *(shapeTensorInv+4)*C[1][4] + *(shapeTensorInv+5)*C[4][5] + *(shapeTensorInv+6)*C[4][4] + *(shapeTensorInv+7)*C[4][5] + *(shapeTensorInv+8)*C[2][4];
    hourglassStiff[3] =  *(shapeTensorInv)*C[0][5] + *(shapeTensorInv+1)*C[5][5] + *(shapeTensorInv+2)*C[3][5] + *(shapeTensorInv+3)*C[5][5] + *(shapeTensorInv+4)*C[1][5] + *(shapeTensorInv+5)*C[3][4] + *(shapeTensorInv+6)*C[3][5] + *(shapeTensorInv+7)*C[3][4] + *(shapeTensorInv+8)*C[2][3];
    hourglassStiff[4] =  *(shapeTensorInv)*C[0][1] + *(shapeTensorInv+1)*C[1][5] + *(shapeTensorInv+2)*C[1][4] + *(shapeTensorInv+3)*C[1][5] + *(shapeTensorInv+4)*C[1][1] + *(shapeTensorInv+5)*C[1][3] + *(shapeTensorInv+6)*C[1][4] + *(shapeTensorInv+7)*C[1][3] + *(shapeTensorInv+8)*C[1][2];
    hourglassStiff[5] =  *(shapeTensorInv)*C[0][3] + *(shapeTensorInv+1)*C[3][4] + *(shapeTensorInv+2)*C[4][5] + *(shapeTensorInv+3)*C[3][4] + *(shapeTensorInv+4)*C[1][3] + *(shapeTensorInv+5)*C[3][3] + *(shapeTensorInv+6)*C[4][5] + *(shapeTensorInv+7)*C[3][3] + *(shapeTensorInv+8)*C[2][3];
    hourglassStiff[6] =  *(shapeTensorInv)*C[0][4] + *(shapeTensorInv+1)*C[3][5] + *(shapeTensorInv+2)*C[4][4] + *(shapeTensorInv+3)*C[3][5] + *(shapeTensorInv+4)*C[1][4] + *(shapeTensorInv+5)*C[4][5] + *(shapeTensorInv+6)*C[4][4] + *(shapeTensorInv+7)*C[4][5] + *(shapeTensorInv+8)*C[2][4];
    hourglassStiff[7] =  *(shapeTensorInv)*C[0][3] + *(shapeTensorInv+1)*C[3][4] + *(shapeTensorInv+2)*C[4][5] + *(shapeTensorInv+3)*C[3][4] + *(shapeTensorInv+4)*C[1][3] + *(shapeTensorInv+5)*C[3][3] + *(shapeTensorInv+6)*C[4][5] + *(shapeTensorInv+7)*C[3][3] + *(shapeTensorInv+8)*C[2][3];
    hourglassStiff[8] =  *(shapeTensorInv)*C[0][2] + *(shapeTensorInv+1)*C[2][3] + *(shapeTensorInv+2)*C[2][4] + *(shapeTensorInv+3)*C[2][3] + *(shapeTensorInv+4)*C[1][2] + *(shapeTensorInv+5)*C[2][3] + *(shapeTensorInv+6)*C[2][4] + *(shapeTensorInv+7)*C[2][3] + *(shapeTensorInv+8)*C[2][2];

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
       MATRICES::MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, temp);
      // \sigma_rot = R * temp
       MATRICES::MatrixMultiply(false, false, 1.0, rotTensor, temp, rotatedStress);
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
)
{
  const double* delta = horizon;
  double* w0 = weightedVolume;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;
  double neighborVolume, omega;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3, w0++){

    // Zero out the weighted volume
    *w0 = 0.0;

    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){

      neighborIndex = *neighborListPtr;
      neighborVolume = jacobianDeterminant[neighborIndex] * volume[neighborIndex];
      neighborCoord = coordinates + 3*neighborIndex;

      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);
      deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                deformedBondY*deformedBondY +
                                deformedBondZ*deformedBondZ);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

      *w0 += omega * neighborVolume;
    }
  }
}

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
)
{
  const double* delta = horizon;
  double* w0 = weightedVolume;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;

  const double* flyingPointFlg = flyingPointFlag;
  const double *bondDamagePtr = bondDamage;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;
  double neighborVolume, omega;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3, w0++, flyingPointFlg++){

    // Zero out the weighted volume
    *w0 = 0.0;

    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++){

      neighborIndex = *neighborListPtr;
      neighborVolume = jacobianDeterminant[neighborIndex] * volume[neighborIndex];
      neighborCoord = coordinates + 3*neighborIndex;

      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);
      deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                deformedBondY*deformedBondY +
                                deformedBondZ*deformedBondZ);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

      *w0 += (1.0 - *bondDamagePtr) * omega * neighborVolume;
    }
  }
}

//This function computes the node-level velocity gradient
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
)
{
  int returnCode = 0;

  const double* delta = horizon;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  ScalarT* shapeTensorInv = shapeTensorInverse;
  ScalarT* velGrad = velocityGradient;
  ScalarT velGradTr;

  const double* flyingPointFlg = flyingPointFlag;
  const double *bondDamagePtr = bondDamage;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;
  ScalarT velStateX, velStateY, velStateZ;
  double neighborVolume, omega, temp;

  std::vector<ScalarT> shapeTensorVector(9);
  ScalarT* shapeTensor = &shapeTensorVector[0];
  ScalarT shapeTensorDeterminant;

  std::vector<ScalarT> velGradFirstTermVector(9);
  ScalarT* velGradFirstTerm = &velGradFirstTermVector[0];

  int inversionReturnCode(0);
  std::string inversionErrorMessage = 
    "**** Error:  computeShapeTensorInverseAndApproximateVelocityGradient: Non-invertible matrix\n";

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3,
      vel+=3, shapeTensorInv+=9, velGrad+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // Zero out data
      *(shapeTensor)   = 0.0 ; *(shapeTensor+1) = 0.0 ; *(shapeTensor+2) = 0.0 ;
      *(shapeTensor+3) = 0.0 ; *(shapeTensor+4) = 0.0 ; *(shapeTensor+5) = 0.0 ;
      *(shapeTensor+6) = 0.0 ; *(shapeTensor+7) = 0.0 ; *(shapeTensor+8) = 0.0 ;
      *(velGradFirstTerm)   = 0.0 ; *(velGradFirstTerm+1) = 0.0 ; *(velGradFirstTerm+2) = 0.0 ;
      *(velGradFirstTerm+3) = 0.0 ; *(velGradFirstTerm+4) = 0.0 ; *(velGradFirstTerm+5) = 0.0 ;
      *(velGradFirstTerm+6) = 0.0 ; *(velGradFirstTerm+7) = 0.0 ; *(velGradFirstTerm+8) = 0.0 ;

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++){

        neighborIndex = *neighborListPtr;
        neighborVolume = jacobianDeterminantN[neighborIndex] * volume[neighborIndex];
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVel = velocities + 3*neighborIndex;

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

        omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

        temp = (1.0 - *bondDamagePtr) * omega * neighborVolume / (deformedBondLength * deformedBondLength);

        *(shapeTensor)   += temp * deformedBondX * deformedBondX;
        *(shapeTensor+1) += temp * deformedBondX * deformedBondY;
        *(shapeTensor+2) += temp * deformedBondX * deformedBondZ;
        *(shapeTensor+3) += temp * deformedBondY * deformedBondX;
        *(shapeTensor+4) += temp * deformedBondY * deformedBondY;
        *(shapeTensor+5) += temp * deformedBondY * deformedBondZ;
        *(shapeTensor+6) += temp * deformedBondZ * deformedBondX;
        *(shapeTensor+7) += temp * deformedBondZ * deformedBondY;
        *(shapeTensor+8) += temp * deformedBondZ * deformedBondZ;

        *(velGradFirstTerm)   += temp * velStateX * deformedBondX;
        *(velGradFirstTerm+1) += temp * velStateX * deformedBondY;
        *(velGradFirstTerm+2) += temp * velStateX * deformedBondZ;
        *(velGradFirstTerm+3) += temp * velStateY * deformedBondX;
        *(velGradFirstTerm+4) += temp * velStateY * deformedBondY;
        *(velGradFirstTerm+5) += temp * velStateY * deformedBondZ;
        *(velGradFirstTerm+6) += temp * velStateZ * deformedBondX;
        *(velGradFirstTerm+7) += temp * velStateZ * deformedBondY;
        *(velGradFirstTerm+8) += temp * velStateZ * deformedBondZ;
      }
      
      inversionReturnCode = MATRICES::Invert3by3Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);
      if(inversionReturnCode > 0){
        returnCode = inversionReturnCode;
        std::cout << inversionErrorMessage;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // Matrix multiply the first term and the shape tensor inverse to compute
      // the velocity gradient
      MATRICES::MatrixMultiply(false, false, 1.0, velGradFirstTerm, shapeTensorInv, velGrad);

      // update volume
      // J dot = J . tr(L)
      velGradTr = *velGrad + *(velGrad+4) + *(velGrad+8);
      jacobianDeterminantNP1[iID] = jacobianDeterminantN[iID] * (1.0 + velGradTr*dt);
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondDamagePtr += numNeighbors; 
    }
  }

  return returnCode;
}

//This function computes the node-level velocity gradient
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
)
{
  int returnCode = 0;

  const double* delta = horizon;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  ScalarT* shapeTensorInv = shapeTensorInverse;
  ScalarT* velGrad = velocityGradient;
  ScalarT* velGradX = velocityGradientX;
  ScalarT* velGradY = velocityGradientY;
  ScalarT* velGradZ = velocityGradientZ;
  ScalarT velGradTr;

  const double* flyingPointFlg = flyingPointFlag;
  const double *bondDamagePtr = bondDamage;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;
  ScalarT velStateX, velStateY, velStateZ;
  double neighborVolume, omega, temp;

  std::vector<ScalarT> shapeTensorVector(9);
  ScalarT* shapeTensor = &shapeTensorVector[0];
  ScalarT shapeTensorDeterminant;

  std::vector<ScalarT> velGradFirstTermVector(9);
  ScalarT* velGradFirstTerm = &velGradFirstTermVector[0];

  int inversionReturnCode(0);
  std::string inversionErrorMessage = 
    "**** Error:  computeShapeTensorInverseAndApproximateVelocityGradient: Non-invertible matrix\n";

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3, vel+=3, shapeTensorInv+=9, 
      velGrad+=9, velGradX+=3, velGradY+=3, velGradZ+=3, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // Zero out data
      *(shapeTensor)   = 0.0 ; *(shapeTensor+1) = 0.0 ; *(shapeTensor+2) = 0.0 ;
      *(shapeTensor+3) = 0.0 ; *(shapeTensor+4) = 0.0 ; *(shapeTensor+5) = 0.0 ;
      *(shapeTensor+6) = 0.0 ; *(shapeTensor+7) = 0.0 ; *(shapeTensor+8) = 0.0 ;
      *(velGradFirstTerm)   = 0.0 ; *(velGradFirstTerm+1) = 0.0 ; *(velGradFirstTerm+2) = 0.0 ;
      *(velGradFirstTerm+3) = 0.0 ; *(velGradFirstTerm+4) = 0.0 ; *(velGradFirstTerm+5) = 0.0 ;
      *(velGradFirstTerm+6) = 0.0 ; *(velGradFirstTerm+7) = 0.0 ; *(velGradFirstTerm+8) = 0.0 ;

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++){

        neighborIndex = *neighborListPtr;
        neighborVolume = jacobianDeterminantN[neighborIndex] * volume[neighborIndex];
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVel = velocities + 3*neighborIndex;

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

        omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

        temp = (1.0 - *bondDamagePtr) * omega * neighborVolume / (deformedBondLength * deformedBondLength);

        *(shapeTensor)   += temp * deformedBondX * deformedBondX;
        *(shapeTensor+1) += temp * deformedBondX * deformedBondY;
        *(shapeTensor+2) += temp * deformedBondX * deformedBondZ;
        *(shapeTensor+3) += temp * deformedBondY * deformedBondX;
        *(shapeTensor+4) += temp * deformedBondY * deformedBondY;
        *(shapeTensor+5) += temp * deformedBondY * deformedBondZ;
        *(shapeTensor+6) += temp * deformedBondZ * deformedBondX;
        *(shapeTensor+7) += temp * deformedBondZ * deformedBondY;
        *(shapeTensor+8) += temp * deformedBondZ * deformedBondZ;

        *(velGradFirstTerm)   += temp * velStateX * deformedBondX;
        *(velGradFirstTerm+1) += temp * velStateX * deformedBondY;
        *(velGradFirstTerm+2) += temp * velStateX * deformedBondZ;
        *(velGradFirstTerm+3) += temp * velStateY * deformedBondX;
        *(velGradFirstTerm+4) += temp * velStateY * deformedBondY;
        *(velGradFirstTerm+5) += temp * velStateY * deformedBondZ;
        *(velGradFirstTerm+6) += temp * velStateZ * deformedBondX;
        *(velGradFirstTerm+7) += temp * velStateZ * deformedBondY;
        *(velGradFirstTerm+8) += temp * velStateZ * deformedBondZ;
      }
      
      inversionReturnCode = MATRICES::Invert3by3Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);
      if(inversionReturnCode > 0){
        returnCode = inversionReturnCode;
        std::cout << inversionErrorMessage;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // Matrix multiply the first term and the shape tensor inverse to compute
      // the velocity gradient
      MATRICES::MatrixMultiply(false, false, 1.0, velGradFirstTerm, shapeTensorInv, velGrad);

      *(velGradX+0) = *(velGrad+0); *(velGradX+1) = *(velGrad+1); *(velGradX+2) = *(velGrad+2); 
      *(velGradY+0) = *(velGrad+3); *(velGradY+1) = *(velGrad+4); *(velGradY+2) = *(velGrad+5); 
      *(velGradZ+0) = *(velGrad+6); *(velGradZ+1) = *(velGrad+7); *(velGradZ+2) = *(velGrad+8); 

      // update volume
      // J dot = J . tr(L)
      velGradTr = *velGrad + *(velGrad+4) + *(velGrad+8);
      jacobianDeterminantNP1[iID] = jacobianDeterminantN[iID] * (1.0 + velGradTr*dt);
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondDamagePtr += numNeighbors; 
    }
  }

  return returnCode;
}

//This function computes the node-level velocity gradient
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
)
{
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  const ScalarT* phi1 = gradientWeight1;
  const ScalarT* phi2 = gradientWeight2;
  const ScalarT* phi3 = gradientWeight3;
  ScalarT* velGrad = velocityGradient;
  ScalarT* velGradX = velocityGradientX;
  ScalarT* velGradY = velocityGradientY;
  ScalarT* velGradZ = velocityGradientZ;
  ScalarT velGradTr;

  const double* flyingPointFlg = flyingPointFlag;

  std::vector<ScalarT> velStateVector(3) ; ScalarT* velState = &velStateVector[0];

  double neighborVolume; // temp, omega, 

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, vel+=3, velGrad+=9, 
      velGradX+=3, velGradY+=3, velGradZ+=3, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.5){

      // Zero out data
      *(velGrad+0) = 0.0 ; *(velGrad+1) = 0.0 ; *(velGrad+2) = 0.0 ;
      *(velGrad+3) = 0.0 ; *(velGrad+4) = 0.0 ; *(velGrad+5) = 0.0 ;
      *(velGrad+6) = 0.0 ; *(velGrad+7) = 0.0 ; *(velGrad+8) = 0.0 ;

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, phi1++, phi2++, phi3++){

        neighborIndex = *neighborListPtr;
        neighborVolume = jacobianDeterminantN[neighborIndex] * volume[neighborIndex];
        neighborVel = velocities + 3*neighborIndex;

        for(int i=0; i<3; i++){
          *(velState+i) = *(neighborVel+i) - *(vel+i);
        }

        for(int i=0; i<3; i++){
          *(velGrad+i*3+0) += *(velState+i) * *phi1 * neighborVolume;
          *(velGrad+i*3+1) += *(velState+i) * *phi2 * neighborVolume;
          *(velGrad+i*3+2) += *(velState+i) * *phi3 * neighborVolume;
        }
      }
      
      *(velGradX+0) = *(velGrad+0); *(velGradX+1) = *(velGrad+1); *(velGradX+2) = *(velGrad+2); 
      *(velGradY+0) = *(velGrad+3); *(velGradY+1) = *(velGrad+4); *(velGradY+2) = *(velGrad+5); 
      *(velGradZ+0) = *(velGrad+6); *(velGradZ+1) = *(velGrad+7); *(velGradZ+2) = *(velGrad+8); 

      // update volume
      // J dot = J . tr(L)
      velGradTr = *velGrad + *(velGrad+4) + *(velGrad+8);
      jacobianDeterminantNP1[iID] = jacobianDeterminantN[iID] * (1.0 + velGradTr*dt);
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      phi1 += numNeighbors; phi2 += numNeighbors; phi3 += numNeighbors; 
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
)
{

    double reducedYieldStress = yieldStress;
    if (isFlaw){
      if (type == 1){
        reducedYieldStress = yieldStress * (1.0 - flawMagnitude 
                          * exp( ((- ( *(modelCoord) - flawLocationX)) * (*(modelCoord) - flawLocationX) -
                          (( *(modelCoord+1) - flawLocationY)) * (*(modelCoord+1) - flawLocationY) -
                          (( *(modelCoord+2) - flawLocationZ)) * (*(modelCoord+2) - flawLocationZ)
                          ) / flawSize / flawSize
                          ));
      }
      else{
        std::cout<<"not implemented yet"<<std::endl;
      }
    }
    return reducedYieldStress;
}
template<typename ScalarT>
void addTemperatureStrain
(
  const double alpha[][3],
  const ScalarT temperature,
  ScalarT* strain
)
{

  *(strain)   -= alpha[0][0] * temperature;  
  *(strain+1) -= alpha[0][1] * temperature;
  *(strain+2) -= alpha[0][2] * temperature;
  *(strain+3) -= alpha[1][0] * temperature;
  *(strain+4) -= alpha[1][1] * temperature;
  *(strain+5) -= alpha[1][2] * temperature;
  *(strain+6) -= alpha[2][0] * temperature;
  *(strain+7) -= alpha[2][1] * temperature;
  *(strain+8) -= alpha[2][2] * temperature;

}

template void addTemperatureStrain<Sacado::Fad::DFad<double>>
(
const double alpha[][3],
const Sacado::Fad::DFad<double> temperature,
Sacado::Fad::DFad<double>* strain
);
template void addTemperatureStrain<double>
(
const double alpha[][3],
const double temperature,
double* strain
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
)
{   std::string logStrainErrorMessage = "**** Error:  elastic_correspondence ::updateElasticCauchyStressAnisotropic() failed to compute logStrain.\n";
  for(int iID=0 ; iID<numPoints ; ++iID, defGrad+=9, strain+=9, temperature+=1){
    if (hencky){
      int defGradLogReturnCode = CORRESPONDENCE::computeLogStrain(defGrad,strain);
      TEUCHOS_TEST_FOR_TERMINATION(defGradLogReturnCode != 0, logStrainErrorMessage);
      }
    else
      {CORRESPONDENCE::computeGreenLagrangeStrain(defGrad,strain);}
    if (applyThermalStrains){
      CORRESPONDENCE::addTemperatureStrain(alpha,*(temperature),strain);
      }

  }

}

template void getStrain<Sacado::Fad::DFad<double>>
(
const int numPoints,
const Sacado::Fad::DFad<double>* defGrad,
const double alpha[][3],
const Sacado::Fad::DFad<double>* temperature,
const bool hencky,
const bool applyThermalStrains,
Sacado::Fad::DFad<double>* strain
);
template void getStrain<double>
(
const int numPoints,
const double* defGrad,
const double alpha[][3],
const double* temperature,
const bool hencky,
const bool applyThermalStrains,
double* strain
);


template<typename ScalarT>
void computeGreenLagrangeStrain
(
  const ScalarT* defGrad,
  ScalarT* strain
)
{

  *(strain)   = 0.5 * ( *(defGrad)   * *(defGrad)   + *(defGrad+3) * *(defGrad+3) + *(defGrad+6) * *(defGrad+6)  - 1.0 );
  *(strain+1) = 0.5 * ( *(defGrad)   * *(defGrad+1) + *(defGrad+3) * *(defGrad+4) + *(defGrad+6) * *(defGrad+7)  );
  *(strain+2) = 0.5 * ( *(defGrad)   * *(defGrad+2) + *(defGrad+3) * *(defGrad+5) + *(defGrad+6) * *(defGrad+8)  );
  *(strain+3) = 0.5 * ( *(defGrad)   * *(defGrad+1) + *(defGrad+3) * *(defGrad+4) + *(defGrad+6) * *(defGrad+7)  );
  *(strain+4) = 0.5 * ( *(defGrad+1) * *(defGrad+1) + *(defGrad+4) * *(defGrad+4) + *(defGrad+7) * *(defGrad+7)  - 1.0 );
  *(strain+5) = 0.5 * ( *(defGrad+1) * *(defGrad+2) + *(defGrad+4) * *(defGrad+5) + *(defGrad+7) * *(defGrad+8)  );
  *(strain+6) = 0.5 * ( *(defGrad)   * *(defGrad+2) + *(defGrad+3) * *(defGrad+5) + *(defGrad+6) * *(defGrad+8)  );
  *(strain+7) = 0.5 * ( *(defGrad+1) * *(defGrad+2) + *(defGrad+4) * *(defGrad+5) + *(defGrad+7) * *(defGrad+8)  );
  *(strain+8) = 0.5 * ( *(defGrad+2) * *(defGrad+2) + *(defGrad+5) * *(defGrad+5) + *(defGrad+8) * *(defGrad+8)  - 1.0 );
}


template void computeGreenLagrangeStrain<Sacado::Fad::DFad<double>>
(
const Sacado::Fad::DFad<double>* defGrad, 
Sacado::Fad::DFad<double>* strain

);
template void computeGreenLagrangeStrain<double>
(
const double* defGrad, 
double* strain

);







template<typename ScalarT>
void computeGreenLagrangeStrain
(
    const ScalarT* deformationGradient,
    ScalarT* greenLagrangeStrain,
    const double* flyingPointFlag,
    int numPoints
)
{
  // Green-Lagrange Strain E = 0.5*(F^T F - I)

  const ScalarT* defGrad = deformationGradient;
  ScalarT* strain = greenLagrangeStrain;

  const double* flyingPointFlg = flyingPointFlag;

  std::vector<ScalarT> tempVector(9);
  ScalarT* temp = &tempVector[0];

  for(int iID=0 ; iID<numPoints ; ++iID, defGrad+=9, strain+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      MATRICES::MatrixMultiply(true, false, 1.0, defGrad, defGrad, temp);

      for(int i=0; i<9; i++)
        *(strain+i) = 0.5 * *(temp+i);

      *(strain+0) -= 0.5;
      *(strain+4) -= 0.5;
      *(strain+8) -= 0.5;
    }
  }
}

//Performs kinematic computations following Flanagan and Taylor (1987), returns
//unrotated rate-of-deformation and rotation tensors
//This function computes the node-level values
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
)
{
  int returnCode = 0;

  const ScalarT* eulerianVelGrad = velocityGradient;
  const ScalarT* leftStretchN = leftStretchTensorN;
  const ScalarT* rotTensorN = rotationTensorN;

  ScalarT* leftStretchNP1 = leftStretchTensorNP1;
  ScalarT* rotTensorNP1 = rotationTensorNP1;
  ScalarT* unrotRateOfDef = unrotatedRateOfDeformation;

  const double* flyingPointFlg = flyingPointFlag;

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

  for(int iID=0 ; iID<numPoints ; ++iID, eulerianVelGrad+=9, rotTensorN+=9, 
      rotTensorNP1+=9, leftStretchNP1+=9, leftStretchN+=9, unrotRateOfDef+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

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
    }
  }

  return returnCode;
}

template<typename ScalarT>
void rotateCauchyStress
(
    const ScalarT* rotationTensor,
    const ScalarT* unrotatedCauchyStress,
    ScalarT* rotatedCauchyStress,
    const double* flyingPointFlag,
    int numPoints
)
{
  const ScalarT* rotTensor = rotationTensor;
  const ScalarT* unrotatedStress = unrotatedCauchyStress;
  ScalarT* rotatedStress = rotatedCauchyStress;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT temp[9];

  for(int iID=0 ; iID<numPoints ; ++iID, 
        rotTensor+=9, unrotatedStress+=9, rotatedStress+=9, flyingPointFlg++){ 

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // temp = \sigma_unrot * Rt
       MATRICES::MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, temp);
      // \sigma_rot = R * temp
       MATRICES::MatrixMultiply(false, false, 1.0, rotTensor, temp, rotatedStress);
    }
  }
}


template<typename ScalarT>
void updateGradientWeightEvaluationFlag
(
    const ScalarT* damageN,
    const ScalarT* damageNP1,
    ScalarT* gradientWeightEvaluationFlag,
    int numPoints
)
{
  const ScalarT* dmgN = damageN;
  const ScalarT* dmgNP1 = damageNP1;
  ScalarT* flag = gradientWeightEvaluationFlag;

  double tol = 1.0e-15;

  for(int iID=0 ; iID<numPoints ; ++iID, dmgN++, dmgNP1++, flag++){
    if(*dmgNP1 - *dmgN > tol)
      *flag = 0.0; // this will cause the weights to be recomputed
  }
}


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
)
{
  //Using RK Implicit Gradient (Same as PD Gradient Operator) to construct
  //these shapes on the parametric space (2D gradients).
  //GMLS or other methods can be incorporated later.

  int returnCode = 0;

  const double* delta = horizon;
  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;

  ScalarT* phi1 = gradientWeightX;
  ScalarT* phi2 = gradientWeightY;
  ScalarT* phi3 = gradientWeightZ;

  ScalarT* flag = gradientWeightEvaluationFlag;

  const double *bondDamagePtr = bondDamage;

  double *omega = influenceState;

  ScalarT undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;

  double neighborVolume, temp;

  int Qdim;

  int counter, thisOrder, p1, p2, p3;
  int i, j;

  // calculate dimension of Q vector
  Qdim = 0;
  for(thisOrder=1; thisOrder<=accuracyOrder; thisOrder++){
    for(p1=thisOrder; p1>=0; p1--){ // x-power
      for(p2=thisOrder-p1; p2>=0; p2--){ // y-power
        p3=thisOrder-p1-p2; //z-power
        Qdim++;
      }
    }
  }

  std::vector<ScalarT> QVector(Qdim);
  ScalarT* Q = &QVector[0];

  double phi1const = 1.0; int phi1ind = 0;
  double phi2const = 1.0; int phi2ind = 1;
  double phi3const = 1.0; int phi3ind = 2;

  std::vector<ScalarT> MVector(Qdim*Qdim);
  ScalarT* M = &MVector[0];

  std::vector<ScalarT> MinvVector(Qdim*Qdim);
  ScalarT* Minv = &MinvVector[0];

  double thresVal;

  int inversionReturnCode(0);

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, flag++){

    if(*flag < 0.5){

      // Zero out data
      for(j=0; j<Qdim; j++)
        for(i=0; i<Qdim; i++)
          *(M+Qdim*i+j) = 0.0;

      // Calculate Moment matrix
      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++, omega++){

        neighborIndex = *neighborListPtr;
        neighborModelCoord = modelCoordinates + 3*neighborIndex;
        neighborVolume = volume[neighborIndex];

        undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
        undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
        undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
        undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                    undeformedBondY*undeformedBondY +
                                    undeformedBondZ*undeformedBondZ);

        // Calculate Q for this bond
        counter = 0;
        for(thisOrder=1; thisOrder<=accuracyOrder; thisOrder++){
          for(p1=thisOrder; p1>=0; p1--){ // x-power
            for(p2=thisOrder-p1; p2>=0; p2--){ // y-power
              p3=thisOrder-p1-p2; //z-power

              Q[counter] = 1.0;
              for(i=0; i<p1; i++)
                Q[counter] *= undeformedBondX / *delta;
              for(i=0; i<p2; i++)
                Q[counter] *= undeformedBondY / *delta;
              for(i=0; i<p3; i++)
                Q[counter] *= undeformedBondZ / *delta;

              counter++;
            }
          }
        }

        *omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta);
        *omega *= (1.0 - *bondDamagePtr);

        temp = *omega * neighborVolume;
        for(j=0; j<Qdim; j++)
          for(i=0; i<Qdim; i++)
            *(M+Qdim*i+j) += temp * Q[i] * Q[j];
      }

      thresVal = 1.0e-6 * *delta * *delta * *delta; // should be area (2d) or vol (3d) multiplied by a threshold value

      // calculate the inverse of the moment matrix (this must be a full rank matrix)
      //inversionReturnCode = computeSymmetrixMatrixInverse(M, Qdim, Minv);
      inversionReturnCode = MATRICES::invertAndCond(M, Minv, Qdim, thresVal);

      if(inversionReturnCode > 0)
        returnCode = inversionReturnCode;

      // Re-iterate over the neighbor set and compute Phis
      // Return the neighbor pointers to the beginning of set
      neighborListPtr -= numNeighbors; 
      omega -= numNeighbors; 
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, phi1++, phi2++, phi3++, omega++){

        neighborIndex = *neighborListPtr;
        neighborModelCoord = modelCoordinates + 3*neighborIndex;

        undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
        undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
        undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
        undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                    undeformedBondY*undeformedBondY +
                                    undeformedBondZ*undeformedBondZ);

        // Calculate Q for this bond
        counter = 0;
        for(thisOrder=1; thisOrder<=accuracyOrder; thisOrder++){
          for(p1=thisOrder; p1>=0; p1--){ // x-power
            for(p2=thisOrder-p1; p2>=0; p2--){ // y-power
              p3=thisOrder-p1-p2; //z-power

              Q[counter] = 1.0;
              for(i=0; i<p1; i++)
                Q[counter] *= undeformedBondX / *delta;
              for(i=0; i<p2; i++)
                Q[counter] *= undeformedBondY / *delta;
              for(i=0; i<p3; i++)
                Q[counter] *= undeformedBondZ / *delta;

              counter++;
            }
          }
        }

        // caluclate phi1
        *phi1 = 0.0;
        temp = phi1const * *omega;  
        for(j=0; j<Qdim; j++)
          *phi1 += temp * *(Minv+Qdim*phi1ind+j) * Q[j] / *delta;

        // caluclate phi2
        *phi2 = 0.0;
        temp = phi2const * *omega;  
        for(j=0; j<Qdim; j++)
          *phi2 += temp * *(Minv+Qdim*phi2ind+j) * Q[j] / *delta;

        // caluclate phi3
        *phi3 = 0.0;
        temp = phi3const * *omega;  
        for(j=0; j<Qdim; j++)
          *phi3 += temp * *(Minv+Qdim*phi3ind+j) * Q[j] / *delta;
      }

      // indicate that the weights are evaluated (no recompute is needed unless
      // damage grows)
      *flag = 1.0;
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      phi1 += numNeighbors; phi2 += numNeighbors; phi3 += numNeighbors; 
      bondDamagePtr += numNeighbors; 
      omega += numNeighbors; 
    }
  }

  return returnCode;
}

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
)
{
  const ScalarT* disp = displacements;
  const ScalarT* neighborDisp;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  const ScalarT* phi1 = gradientWeightX;
  const ScalarT* phi2 = gradientWeightY;
  const ScalarT* phi3 = gradientWeightZ;
  ScalarT* defGradX = deformationGradientX;
  ScalarT* defGradY = deformationGradientY;
  ScalarT* defGradZ = deformationGradientZ;
  ScalarT* defGradDotX = deformationGradientDotX;
  ScalarT* defGradDotY = deformationGradientDotY;
  ScalarT* defGradDotZ = deformationGradientDotZ;

  std::vector<ScalarT> defGradVector(9) ; ScalarT* defGrad = &defGradVector[0];
  std::vector<ScalarT> defGradDotVector(9) ; ScalarT* defGradDot = &defGradDotVector[0];

  std::vector<ScalarT> dispStateVector(3) ; ScalarT* dispState = &dispStateVector[0];
  std::vector<ScalarT> velStateVector(3) ; ScalarT* velState = &velStateVector[0];

  double neighborVolume;//, temp;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, disp+=3, defGradX+=3, defGradY+=3, defGradZ+=3, 
      vel+=3, defGradDotX+=3, defGradDotY+=3, defGradDotZ+=3){

    // Initialize data
    *(defGrad+0) = 1.0 ; *(defGrad+1) = 0.0 ; *(defGrad+2) = 0.0 ;
    *(defGrad+3) = 0.0 ; *(defGrad+4) = 1.0 ; *(defGrad+5) = 0.0 ;
    *(defGrad+6) = 0.0 ; *(defGrad+7) = 0.0 ; *(defGrad+8) = 1.0 ;
    *(defGradDot+0) = 0.0 ; *(defGradDot+1) = 0.0 ; *(defGradDot+2) = 0.0 ;
    *(defGradDot+3) = 0.0 ; *(defGradDot+4) = 0.0 ; *(defGradDot+5) = 0.0 ;
    *(defGradDot+6) = 0.0 ; *(defGradDot+7) = 0.0 ; *(defGradDot+8) = 0.0 ;

    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, phi1++, phi2++, phi3++){

      neighborIndex = *neighborListPtr;
      neighborVolume = volume[neighborIndex];
      neighborDisp = displacements + 3*neighborIndex;
      neighborVel = velocities + 3*neighborIndex;

      for(int i=0; i<3; i++){
        *(dispState+i) = *(neighborDisp+i) - *(disp+i);
        *(velState+i) = *(neighborVel+i) - *(vel+i);
      }

      for(int i=0; i<3; i++){
        *(defGrad+i*3+0) += *(dispState+i) * *phi1 * neighborVolume;
        *(defGrad+i*3+1) += *(dispState+i) * *phi2 * neighborVolume;
        *(defGrad+i*3+2) += *(dispState+i) * *phi3 * neighborVolume;
        *(defGradDot+i*3+0) += *(velState+i) * *phi1 * neighborVolume;
        *(defGradDot+i*3+1) += *(velState+i) * *phi2 * neighborVolume;
        *(defGradDot+i*3+2) += *(velState+i) * *phi3 * neighborVolume;
      }
    }
    
    *(defGradX+0) = *(defGrad+0); *(defGradX+1) = *(defGrad+1); *(defGradX+2) = *(defGrad+2); 
    *(defGradY+0) = *(defGrad+3); *(defGradY+1) = *(defGrad+4); *(defGradY+2) = *(defGrad+5); 
    *(defGradZ+0) = *(defGrad+6); *(defGradZ+1) = *(defGrad+7); *(defGradZ+2) = *(defGrad+8); 
    *(defGradDotX+0) = *(defGradDot+0); *(defGradDotX+1) = *(defGradDot+1); *(defGradDotX+2) = *(defGradDot+2); 
    *(defGradDotY+0) = *(defGradDot+3); *(defGradDotY+1) = *(defGradDot+4); *(defGradDotY+2) = *(defGradDot+5); 
    *(defGradDotZ+0) = *(defGradDot+6); *(defGradDotZ+1) = *(defGradDot+7); *(defGradDotZ+2) = *(defGradDot+8); 
  }
}


template<typename ScalarT>
void computeWeightedVolume
(
    const double* volume,
    ScalarT* weightedVolume,
    const double* influenceState,
    const int* neighborhoodList,
    int numPoints
)
{
  ScalarT* w0 = weightedVolume;
  const double *omega = influenceState;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, w0++){

    // Zero out the weighted volume
    *w0 = 0.0;

    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, omega++){

      neighborIndex = *neighborListPtr;

      *w0 += *omega * volume[neighborIndex];
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
    double dt
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
)
{
  // Green-Lagrange Strain rate Edot = 0.5*(Fdot^T F + F^T Fdot)

  const ScalarT* defGradX = deformationGradientX;
  const ScalarT* defGradY = deformationGradientY;
  const ScalarT* defGradZ = deformationGradientZ;
  const ScalarT* defGradDotX = deformationGradientDotX;
  const ScalarT* defGradDotY = deformationGradientDotY;
  const ScalarT* defGradDotZ = deformationGradientDotZ;
  const ScalarT* strainN = greenLagrangeStrainN;
  ScalarT* strainNP1 = greenLagrangeStrainNP1;

  std::vector<ScalarT> tempVector(9);
  ScalarT* temp = &tempVector[0];

  std::vector<ScalarT> defGradVector(9) ; ScalarT* defGrad = &defGradVector[0];
  std::vector<ScalarT> defGradDotVector(9) ; ScalarT* defGradDot = &defGradDotVector[0];

  for(int iID=0 ; iID<numPoints ; ++iID, defGradX+=3, defGradY+=3, defGradZ+=3, 
      defGradDotX+=3, defGradDotY+=3, defGradDotZ+=3, strainN+=9, strainNP1+=9){

    for(int i=0; i<3; i++){
      *(defGrad+i+0) = *(defGradX+i);
      *(defGrad+i+3) = *(defGradY+i);
      *(defGrad+i+6) = *(defGradZ+i);
      *(defGradDot+i+0) = *(defGradDotX+i);
      *(defGradDot+i+3) = *(defGradDotY+i);
      *(defGradDot+i+6) = *(defGradDotZ+i);
    }

    MATRICES::MatrixMultiply(true, false, 1.0, defGradDot, defGrad, temp);

    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
        *(strainNP1+3*i+j) = *(strainN+3*i+j) + 0.5 * ( *(temp+3*i+j) + *(temp+3*j+i) ) * dt;
  }
}



/** Explicit template instantiation for double. */

template void rotateCauchyStress<double>
(
 const double* rotationTensor,
 const double* unrotatedCauchyStress,
 double* rotatedCauchyStress,
 int numPoints
 );



void computeHeatFlowState_correspondence(    
    const double* modelCoord,
    const int numOwnedPoints,
    const int* neighborhoodList,
    const double* shapeTensorInverse,
    const double* temperature,
    const double* horizon,
    const double* kappa,
    const double* volume,
    const double* detachedNodes,
    const double* bondDamage,
    const bool twoD,
    double* heatFlowState
    )

    {
     DIFFUSION::computeHeatFlowState_correspondence(
                                  modelCoord,
                                  numOwnedPoints,
                                  neighborhoodList,
                                  shapeTensorInverse,
                                  temperature,
                                  horizon,
                                  kappa,
                                  volume,
                                  detachedNodes,
                                  bondDamage,
                                  twoD,
                                  heatFlowState); 
    }

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template<typename ScalarT>
void createRotatedStiff
(
const ScalarT C[][6],
const ScalarT* rotMat,
ScalarT Cnew[][6]
){
    
     //  rM[0,0]**2     , rM[0,1]**2     , rM[0,2]**2     , 2.*rM[0,1]*rM[0,2]             , 2.*rM[0,0]*rM[0,2]             , 2.*rM[0,0]*rM[0,1]              ],
     //[ rM[1,0]**2     , rM[1,1]**2     , rM[1,2]**2     , 2.*rM[1,1]*rM[1,2]             , 2.*rM[1,0]*rM[1,2]             , 2.*rM[1,0]*rM[1,1]              ],
     //[ rM[2,0]**2     , rM[2,1]**2     , rM[2,2]**2     , 2.*rM[2,1]*rM[2,2]             , 2.*rM[2,0]*rM[2,2]             , 2.*rM[2,0]*rM[2,1]              ],
     //[ rM[1,0]*rM[2,0], rM[1,1]*rM[2,1], rM[1,2]*rM[2,2], rM[1,1]*rM[2,2]+rM[1,2]*rM[2,1], rM[1,0]*rM[2,2]+rM[1,2]*rM[2,0], rM[1,0]*rM[2,1]+rM[1,1]*rM[2,0] ],
     //[ rM[0,0]*rM[2,0], rM[0,1]*rM[2,1], rM[0,2]*rM[2,2], rM[0,1]*rM[2,2]+rM[0,2]*rM[2,1], rM[0,0]*rM[2,2]+rM[0,2]*rM[2,0], rM[0,0]*rM[2,1]+rM[0,1]*rM[2,0] ],
     //[ rM[0,0]*rM[1,0], rM[0,1]*rM[1,1], rM[0,2]*rM[1,2], rM[0,1]*rM[1,2]+rM[0,2]*rM[1,1], rM[0,0]*rM[1,2]+rM[0,2]*rM[1,0], rM[0,0]*rM[1,1]+rM[0,1]*rM[1,0] ],
     //])
    
    
    ScalarT tm[6][6];
    bool GlobToLoc = false; 
    if (GlobToLoc){
    
        tm[0][0] = *(rotMat)**(rotMat); tm[0][1] = *(rotMat+1)**(rotMat+1); tm[0][2] = *(rotMat+2)**(rotMat+2); tm[0][3] = 2**(rotMat+1)**(rotMat+2); tm[0][4] = 2**(rotMat)**(rotMat+2); tm[0][5] = 2**(rotMat)**(rotMat+1);
        tm[1][0] = *(rotMat+3)**(rotMat+3); tm[1][1] = *(rotMat+4)**(rotMat+4); tm[1][2] = *(rotMat+5)**(rotMat+5); tm[1][3] = 2**(rotMat+4)**(rotMat+5); tm[1][4] = 2**(rotMat+3)**(rotMat+5); tm[1][5] = 2**(rotMat+3)**(rotMat+4);  
        tm[2][0] = *(rotMat+6)**(rotMat+6); tm[2][1] = *(rotMat+7)**(rotMat+7); tm[2][2] = *(rotMat+8)**(rotMat+8); tm[2][3] = 2**(rotMat+7)**(rotMat+8); tm[2][4] = 2**(rotMat+6)**(rotMat+8); tm[2][5] = 2**(rotMat+6)**(rotMat+7);
        tm[3][0] = *(rotMat+3)**(rotMat+6); tm[3][1] = *(rotMat+4)**(rotMat+7); tm[3][2] = *(rotMat+5)**(rotMat+8); tm[3][3] = *(rotMat+4)**(rotMat+8)+*(rotMat+5)**(rotMat+7); tm[3][4] = *(rotMat+3)**(rotMat+8)+*(rotMat+5)**(rotMat+6); tm[3][5] =*(rotMat+3)**(rotMat+7)+*(rotMat+4)**(rotMat+6);
        tm[4][0] = *(rotMat)**(rotMat+6); tm[4][1] = *(rotMat+1)**(rotMat+7); tm[4][2] = *(rotMat+2)**(rotMat+8); tm[4][3] = *(rotMat+1)**(rotMat+8)+*(rotMat+2)**(rotMat+7); tm[4][4] = *(rotMat)**(rotMat+8)+*(rotMat+2)**(rotMat+6); tm[4][5] =*(rotMat)**(rotMat+7)+*(rotMat+1)**(rotMat+6);
        tm[5][0] = *(rotMat)**(rotMat+3); tm[5][1] = *(rotMat+1)**(rotMat+4); tm[5][2] = *(rotMat+2)**(rotMat+5); tm[5][3] = *(rotMat+1)**(rotMat+5)+*(rotMat+2)**(rotMat+4); tm[5][4] = *(rotMat)**(rotMat+5)+*(rotMat+2)**(rotMat+3); tm[5][5] =*(rotMat)**(rotMat+4)+*(rotMat+1)**(rotMat+3);
    }
    else{
        tm[0][0] = *(rotMat)**(rotMat);   tm[0][1] = *(rotMat+1)**(rotMat+1);   tm[0][2] = *(rotMat+2)**(rotMat+2);   tm[0][3] = *(rotMat+1)**(rotMat+2); tm[0][4] = *(rotMat)**(rotMat+2); tm[0][5] = *(rotMat)**(rotMat+1);
        tm[1][0] = *(rotMat+3)**(rotMat+3);   tm[1][1] = *(rotMat+4)**(rotMat+4);   tm[1][2] = *(rotMat+5)**(rotMat+5);   tm[1][3] = *(rotMat+4)**(rotMat+5); tm[1][4] = *(rotMat+3)**(rotMat+5); tm[1][5] = *(rotMat+3)**(rotMat+4);  
        tm[2][0] = *(rotMat+6)**(rotMat+6);   tm[2][1] = *(rotMat+7)**(rotMat+7);   tm[2][2] = *(rotMat+8)**(rotMat+8);   tm[2][3] = *(rotMat+7)**(rotMat+8); tm[2][4] = *(rotMat+6)**(rotMat+8); tm[2][5] = *(rotMat+6)**(rotMat+7);
        tm[3][0] = 2**(rotMat+3)**(rotMat+6); tm[3][1] = 2**(rotMat+4)**(rotMat+7); tm[3][2] = 2**(rotMat+5)**(rotMat+8); tm[3][3] = *(rotMat+4)**(rotMat+8)+*(rotMat+5)**(rotMat+7); tm[3][4] = *(rotMat+3)**(rotMat+8)+*(rotMat+5)**(rotMat+6); tm[3][5] =*(rotMat+3)**(rotMat+7)+*(rotMat+4)**(rotMat+6);
        tm[4][0] = 2**(rotMat)**(rotMat+6); tm[4][1] = 2**(rotMat+1)**(rotMat+7); tm[4][2] = 2**(rotMat+2)**(rotMat+8); tm[4][3] = *(rotMat+1)**(rotMat+8)+*(rotMat+2)**(rotMat+7); tm[4][4] = *(rotMat)**(rotMat+8)+*(rotMat+2)**(rotMat+6); tm[4][5] =*(rotMat)**(rotMat+7)+*(rotMat+1)**(rotMat+6);
        tm[5][0] = 2**(rotMat)**(rotMat+3); tm[5][1] = 2**(rotMat+1)**(rotMat+4); tm[5][2] = 2**(rotMat+2)**(rotMat+5); tm[5][3] = *(rotMat+1)**(rotMat+5)+*(rotMat+2)**(rotMat+4); tm[5][4] = *(rotMat)**(rotMat+5)+*(rotMat+2)**(rotMat+3); tm[5][5] =*(rotMat)**(rotMat+4)+*(rotMat+1)**(rotMat+3);
        
        
    }
    ScalarT Ctemp[6][6];
    bool transpose = true;
    MATRICES::MatMul(6,tm ,C ,Ctemp, transpose);
    transpose = false;
    MATRICES::MatMul(6,Ctemp ,tm ,Cnew, transpose);
}
template void createRotatedStiff<Sacado::Fad::DFad<double> >
(
const Sacado::Fad::DFad<double> C[][6],
const Sacado::Fad::DFad<double>* rotMat,
Sacado::Fad::DFad<double>  Cnew[][6]
);


// Explicit template instantiation for double
template void createRotatedStiff<double>
(
const double C[][6],
const double* rotMat,
double Cnew[][6]
);


template int computeLogStrain<double>
(
 const double* defGrad,
 double* strain
);
template int computeLogStrain<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* defGrad,
 Sacado::Fad::DFad<double>*  strain
);

template int computeShapeTensorInverseAndApproximateDeformationGradient<double>
(
const double* volume,
const double* horizon,
const double* modelCoordinates,
const double* coordinatesNP1,
double* shapeTensorInverse,
double* deformationGradient,
const double* bondDamageNP1,
const int* neighborhoodList,
int numPoints,
const bool type,
double* detachedNodes
);

template void computeVelocityGradient<double>
(
    const double* volume,
    const double* jacobianDeterminantN,
    double* jacobianDeterminantNP1,
    const double* velocities,
    const double* gradientWeight1,
    const double* gradientWeight2,
    const double* gradientWeight3,
    double* velocityGradient,
    double* velocityGradientX,
    double* velocityGradientY,
    double* velocityGradientZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints,
    double dt
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
const bool type,
double* detachedNodes
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

//template void computeCorrespondenceStabilityWanEtAl<double>
//(
//const double* volume,
//const double* horizon,
//const double* modelCoordinates,
//const double* coordinates,
//const int* neighborhoodList,
//int numPoints,
//const double* deformationGradient,
//double* shapeTensorInverse,
//const double C[][6],
//const double* bondDamage,
//const double* detachedNodes,
//double* hourglassForceDensity,
//double hourglassCoefficient
//);
template void createHourglassStiffness<double>
(
const double C[][6],
const double alpha[],
const double* shapeTensorInverse,
double* hourglassStiff
);
template void createHourglassStiffness<Sacado::Fad::DFad<double> >
(
const Sacado::Fad::DFad<double> C[][6],
const double alpha[],
const Sacado::Fad::DFad<double>* shapeTensorInverse,
Sacado::Fad::DFad<double>* hourglassStiff
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
double* TSx
);

template void computeCorrespondenceStabilityWanEtAlShort<Sacado::Fad::DFad<double> >
(
const Sacado::Fad::DFad<double> FxsiX,
const Sacado::Fad::DFad<double> FxsiY,
const Sacado::Fad::DFad<double> FxsiZ,
const Sacado::Fad::DFad<double> deformedBondX,
const Sacado::Fad::DFad<double> deformedBondY,
const Sacado::Fad::DFad<double> deformedBondZ,
const Sacado::Fad::DFad<double>* hourglassStiff,
Sacado::Fad::DFad<double>* TSx
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


/** Explicit template instantiation for Sacado::Fad::DFad<double>. */




template int computeGradientWeights<double>
(
    const double* horizon,
    const double* coordinates,
    const double* volume,
    const double* jacobianDeterminant,
    double* gradientWeight1,
    double* gradientWeight2,
    double* gradientWeight3,
    const int accuracyOrder,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints
);

template void computeUndamagedWeightedVolume<double>
(
    const double* volume,
    double* weightedVolume,
    const double* jacobianDeterminant,
    const double* horizon,
    const double* coordinates,
    const int* neighborhoodList,
    int numPoints
);

template void computeWeightedVolume<double>
(
    const double* volume,
    double* weightedVolume,
    const double* jacobianDeterminant,
    const double* horizon,
    const double* coordinates,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints
);

template int computeShapeTensorInverseAndApproximateNodeLevelVelocityGradient<double>
(
    const double* volume,
    const double* jacobianDeterminantN,
    double* jacobianDeterminantNP1,
    const double* horizon,
    const double* coordinates,
    const double* velocities,
    double* shapeTensorInverse,
    double* velocityGradient,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints,
    double dt
);

template int computeShapeTensorInverseAndApproximateNodeLevelVelocityGradient<double>
(
    const double* volume,
    const double* jacobianDeterminantN,
    double* jacobianDeterminantNP1,
    const double* horizon,
    const double* coordinates,
    const double* velocities,
    double* shapeTensorInverse,
    double* velocityGradient,
    double* velocityGradientX,
    double* velocityGradientY,
    double* velocityGradientZ,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints,
    double dt
);
template int EigenVec2D<double>(
    const double* a,
    double* result
);

template int EigenVec2D<Sacado::Fad::DFad<double>>(
    const  Sacado::Fad::DFad<double>* a,
     Sacado::Fad::DFad<double>* result
);

template int computeNodeLevelUnrotatedRateOfDeformationAndRotationTensor<double>
(
    const double* velocityGradient,
    const double* leftStretchTensorN,
    const double* rotationTensorN,
    double* leftStretchTensorNP1,
    double* rotationTensorNP1,
    double* unrotatedRateOfDeformation,
    const double* flyingPointFlag,
    int numPoints,
    double dt
);

template void updateDeformationGradient<double>
(
    const double* velocityGradient,
    const double* deformationGradientN,
    double* deformationGradientNP1,
    const double* flyingPointFlag,
    int numPoints,
    double dt
);

template void computeGreenLagrangeStrain<double>
(
    const double* deformationGradient,
    double* greenLagrangeStrain,
    const double* flyingPointFlag,
    int numPoints
);


template void rotateCauchyStress<double>
(
    const double* rotationTensor,
    const double* unrotatedCauchyStress,
    double* rotatedCauchyStress,
    const double* flyingPointFlag,
    int numPoints
);


template void updateGradientWeightEvaluationFlag<double>
(
    const double* damageN,
    const double* damageNP1,
    double* gradientWeightEvaluationFlag,
    int numPoints
);

template int computeLagrangianGradientWeights<double>
(
    const double* horizon,
    const double* modelCoordinates,
    const double* volume,
    double* gradientWeightX,
    double* gradientWeightY,
    double* gradientWeightZ,
    double* gradientWeightEvaluationFlag,
    const double* bondDamage,
    double* influenceState,
    const int accuracyOrder,
    const int* neighborhoodList,
    int numPoints
);

template void computeDeformationGradient<double>
(
    const double* volume,
    const double* displacements,
    const double* velocities,
    const double* gradientWeightX,
    const double* gradientWeightY,
    const double* gradientWeightZ,
    double* deformationGradientX,
    double* deformationGradientY,
    double* deformationGradientZ,
    double* deformationGradientDotX,
    double* deformationGradientDotY,
    double* deformationGradientDotZ,
    const int* neighborhoodList,
    int numPoints
);

template void computeWeightedVolume<double>
(
    const double* volume,
    double* weightedVolume,
    const double* influenceState,
    const int* neighborhoodList,
    int numPoints
);

template int computeNodeLevelUnrotatedRateOfDeformationAndRotationTensor<double>
(
    const double* deformationGradientX,
    const double* deformationGradientY,
    const double* deformationGradientZ,
    const double* deformationGradientDotX,
    const double* deformationGradientDotY,
    const double* deformationGradientDotZ,
    const double* leftStretchTensorN,
    const double* rotationTensorN,
    double* leftStretchTensorNP1,
    double* rotationTensorNP1,
    double* unrotatedRateOfDeformation,
    int numPoints,
    double dt
);

template void updateGreenLagrangeStrain<double>
(
    const double* deformationGradientX,
    const double* deformationGradientY,
    const double* deformationGradientZ,
    const double* deformationGradientDotX,
    const double* deformationGradientDotY,
    const double* deformationGradientDotZ,
    const double* greenLagrangeStrainN,
    double* greenLagrangeStrainNP1,
    int numPoints,
    double dt
);


template void getLinearUnrotatedRateOfDeformation<double>
(
    const double* volume,
    const double* horizon,
    const double* modelCoordinates,
    const double* velocities,
    double* deformationGradient,
    const double* shapeTensorInverse,
    double* unrotatedRateOfDeformation,
    const int* neighborhoodList,
    int numPoints,
    const double* bondDamage,
    const bool type,
    double* detachedNodes
);

}
