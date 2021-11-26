//! \file elastic_FEM.cxx

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

#include "elastic_correspondence.h"
#include "correspondence.h"
#include "matrices.h"
#include "material_utilities.h"
#include <Sacado.hpp>
#include <math.h>
#include <cmath> 

//#include <Teuchos_Assert.hpp>
//#include <Epetra_SerialComm.h>


using namespace std;
namespace FEM {

template<typename ScalarT>
void updateElasticCauchyStressFEM
(
ScalarT* modelCoordinates,
ScalarT* DeformationGradient, 
ScalarT* unrotatedCauchyStressN, 
ScalarT* unrotatedCauchyStressNP1, 
int numPoints, 
const ScalarT Cstiff[][6],
double* angles,
int type,
double dt,
int order[3]

)
{
  ScalarT* nodalCoor = modelCoordinates;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;
  ScalarT strain[3][3];
  ScalarT rotMat[3][3], rotMatT[3][3], temp[3][3];
  //int defGradLogReturnCode(0);
  bool rotation = false;
  std::vector<ScalarT> sigmaIntVector(9);
  ScalarT* sigmaInt = &sigmaIntVector[0];
  ScalarT* elNodalCoor;
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // for higher order it has to be adapted
  // it is a full integration
  // reduced integration is not included yet
  /////////////////////////////////////////////////////////////////////
  int numInt = (order[0]+1) * (order[1]+1) * (order[2]+1);
  
  std::vector<double> NxiVector(order[0]+1), NetaVector(order[1]+1), NpsiVector(order[2]+1);
  double* Nxi = &NxiVector[0], Neta = &NetaVector[0], Npsi = &NpsiVector[0];
  std::vector<double> BxiVector(order[0]+1), BetaVector(order[1]+1), BpsiVector(order[2]+1);
  double* Bxi = &BxiVector[0], Beta = &BetaVector[0], Bpsi = &BpsiVector[0];
  std::vector<double> elCoorxVector(order[0]+1), elCooryVector(order[1]+1), elCoorzVector(order[2]+1);
  double* elCoorx = &elCoorxVector[0], elCoory = &elCooryVector[0], elCoorz = &elCoorzVector[0];
  std::vector<double> weightxVector(order[0]+1), weightyVector(order[1]+1), weightzVector(order[2]+1);
  double* weightsx = &weightxVector[0], weightsy = &weightyVector[0], weightsz = &weightzVector[0];
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // irgendwie in die init?
  // das wäre der Austausch für verschiedene Elementtypen
  FEM::weightsAndIntegrationPoints(order[0], elCoorx, weightsx);
  FEM::weightsAndIntegrationPoints(order[1], elCoory, weightsy);
  FEM::weightsAndIntegrationPoints(order[2], elCoorz, weightsz);
 
 for (int jID=0 ; jID<numInt ; ++jID, Nxi+=order[0]+1, Neta+=order[1]+1, 
      Npsi+=order[2]+1, Bxi+=order[0]+1, Beta+=order[1]+1, Bpsi+=order[2]+1){
  FEM::getLagrangeElementData(order,elCoorx[jID],elCoory[jID],elCoorz[jID],Nxi,Neta,Npsi,Bxi,Beta,Bpsi);
  //FEM::getBMatrixAtIntegrationPoints
 }


  
  

  
  int dof = 3*(order[0]+1)*(order[1]+1)*(order[2]+1);  
  int numNeigh = *neighPtr; neighPtr++;
  for(int n=0;n<numNeigh;n++,neighPtr++){
              int localId = *neighPtr;
              elNodalCoor = &nodalCoor[localId]
              FEM::getJacobian(elNodalCoor,Nxi,Neta,Npsi,Bxi,Beta,Bpsi)
    for(int iID=0 ; iID<numEL ; ++iID, 
          defGrad+=9, sigmaNP1+=9, angles+=3){
            //for (int k=0...)
            //disp(k) = dispGlob(iID, k);

            //FEM::BMatrix(shape, modelCoorNeigh, BMatrix, detJ);
            for (int jID=0 ; jID<numInt ; ++jID){

              //
              //FEM::computeStrain(BMatrix,disp, dof,strain); --> Bu
              
              //https://www.continuummechanics.org/stressxforms.html
              // Q Q^T * sigma * Q Q^T = Q C Q^T epsilon Q Q^T
              if (rotation){  
                CORRESPONDENCE::createRotationMatrix(angles,rotMat);
                MATRICES::TransposeMatrix(rotMat,rotMatT);
                // geomNL
                MATRICES::MatrixMultiply3x3(rotMatT,strain, temp);
                MATRICES::MatrixMultiply3x3(temp,rotMat,strain);
              }

              CORRESPONDENCE::updateElasticCauchyStressAnisotropicCode(strain, sigmaInt, Cstiff, type);
              // rotation back
              if (rotation){  
                MATRICES::MatrixMultiply3x3fromVector(rotMat,sigmaInt, temp);
                MATRICES::MatrixMultiply3x3toVector(temp,rotMatT,sigmaInt);

              }

              //FEM::getNodelForce(BMatrix,n,m,sigmaInt, detJ, force)
              //force += Btranspose*sigmaInt*detJ;
              
              //sigmaNP1 += sigmaInt;

            }
            // gemittelte Spannungen
            // sigmaNP1 /= numInt;
            //globForce(topo) += force; ??
          }

}

template void updateElasticCauchyStressFEM<Sacado::Fad::DFad<double>>
(
Sacado::Fad::DFad<double>* DeformationGradient, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressN, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1, 
int numPoints, 
const Sacado::Fad::DFad<double> Cstiff[][6],
double* angles,
int type,
double dt,
bool incremental,
bool hencky
);
template void updateElasticCauchyStressFEM<double>
(
double* DeformationGradient, 
double* unrotatedCauchyStressN, 
double* unrotatedCauchyStressNP1, 
int numPoints, 
const double Cstiff[][6],
double* angles,
int type,
double dt,
bool incremental,
bool hencky
);







}


