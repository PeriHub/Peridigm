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
#include "elastic_FEM.h"
#include "FEM_routines.h"
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
  ScalarT* displacements, 
  ScalarT* unrotatedCauchyStressNP1, 
  int numElements, 
  int* elementNodalList,
  const ScalarT Cstiff[][6],
  double* angles,
  int type,
  double dt,
  int order[3]
  )
  { 

    ScalarT* nodalCoor = modelCoordinates;
    ScalarT* disp = displacements;
    ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;
    ScalarT strain[3][3];
    ScalarT rotMat[3][3], rotMatT[3][3], temp[3][3];
    //int defGradLogReturnCode(0);
    bool rotation = false;
    std::vector<ScalarT> sigmaIntVector(9);
    ScalarT* sigmaInt = &sigmaIntVector[0];
    ScalarT* elNodalCoor;
    int* neighPtr = neighborhoodList;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // for higher order it has to be adapted
    // it is a full integration
    // reduced integration is not included yet
    /////////////////////////////////////////////////////////////////////
    int numInt = (order[0]+1) * (order[1]+1) * (order[2]+1);
    int ndof = numInt * 3;
    std::vector<double> dispNodalVector(ndof);
    double* dispNodal = &dispNodalVector[0];
    std::vector<double> coorNodalVector(ndof);
    double* coorNodal = &coorNodalVector[0];
    
    std::vector<double> NxiVector(numInt), NetaVector(numInt), NpsiVector(numInt);
    double* Nxi = &NxiVector[0];
    double* Neta = &NetaVector[0];
    double* Npsi = &NpsiVector[0];
    std::vector<double> BxiVector(numInt), BetaVector(numInt), BpsiVector(numInt);
    double* Bxi = &BxiVector[0];
    double* Beta = &BetaVector[0];
    double* Bpsi = &BpsiVector[0];
    std::vector<double> elCoorxVector(order[0]+1), elCooryVector(order[1]+1), elCoorzVector(order[2]+1);
    double* elCoorx = &elCoorxVector[0];
    double* elCoory = &elCooryVector[0];
    double* elCoorz = &elCoorzVector[0];
    std::vector<double> weightxVector(order[0]+1), weightyVector(order[1]+1), weightzVector(order[2]+1);
    double* weightsx = &weightxVector[0];
    double* weightsy = &weightyVector[0];
    double* weightsz = &weightzVector[0];
    std::vector<double> BMatrixVector(ndof*numInt,6);
    double* &BMatrix;
    int topo[numInt][3];
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // irgendwie in die init?
    // das wäre der Austausch für verschiedene Elementtypen
    // falsch
    FEM::weightsAndIntegrationPoints(order[0], elCoorx, weightsx);
    FEM::weightsAndIntegrationPoints(order[1], elCoory, weightsy);
    FEM::weightsAndIntegrationPoints(order[2], elCoorz, weightsz);
    int count = 0;
    for (int kID=0 ; kID<order[2]+1 ; ++kID){
      for (int jID=0 ; jID<order[1]+1 ; ++jID){
        for (int iID=0 ; iID<order[0]+1 ; ++iID){
          FEM::getLagrangeElementData(order,elCoorx[iID],elCoory[jID],elCoorz[kID],Nxi,Neta,Npsi,Bxi,Beta,Bpsi);
          FEM::BMatrixData(order,Nxi,Neta,Npsi,Bxi,Beta,Bpsi,BMatrix);
          topo[count][0] = iID;
          topo[count][1] = jID;
          topo[count][2] = kID;
          count += 1;
          BMatrix+=ndof;
        }
      }
    }

    
  //for (int iID=0 ; iID<order[0]+1 ; ++iID){  
  // int dof = 3*(order[0]+1)*(order[1]+1)*(order[2]+1);  
  // int numNeigh = *neighPtr; neighPtr++;
  // for(int n=0;n<numNeigh;n++,neighPtr++){
  //             int localId = *neighPtr;
  //             elNodalCoor = &nodalCoor[localId];
  //             FEM::getJacobian(elNodalCoor,Nxi,Neta,Npsi,Bxi,Beta,Bpsi)}
  
  for(int iID=0 ; iID<numElements ; ++iID, sigmaNP1+=9, angles+=3){
        int numNeigh = *neighPtr; 
        for(int n=0 ; n<numNeigh ; ++n){
          neighPtr++;
          int localId = *neighPtr;
          coorNodal[3*n] =   nodalCoor[localId];
          coorNodal[3*n+1] = nodalCoor[localId+1];
          coorNodal[3*n+2] = nodalCoor[localId+2];
          dispNodal[3*n] =   disp[localId];
          dispNodal[3*n+1] = disp[localId+1];
          dispNodal[3*n+2] = disp[localId+2];
        }
        
              //FEM::BMatrix(shape, modelCoorNeigh, BMatrix, detJ);
        for (int jID=0 ; jID<numInt ; ++jID){
                FEM::getJacobian(elNodalCoor,Nxi,Neta,Npsi,Bxi,Beta,Bpsi)
                //
                FEM::computeStrain(BMatrix,disp, ndof,strain); 
                
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
    for(int n=0 ; n<numNeigh ; ++n){
          neighPtr++;
          int localId = *neighPtr;
          //forceGlob(iID, k)=force(k)
        }

              }
              // gemittelte Spannungen
              // sigmaNP1 /= numInt;
              //globForce(topo) += force; ??
            }

  }

  template void updateElasticCauchyStressFEM<Sacado::Fad::DFad<double>>
  (
  Sacado::Fad::DFad<double>* modelCoordinates,
  Sacado::Fad::DFad<double>* deformedCoordinates, 
  Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1, 
  int numPoints, 
  int* neighborhoodList,
  const Sacado::Fad::DFad<double> Cstiff[][6],
  double* angles,
  int type,
  double dt,
  int order[3]
  );
  template void updateElasticCauchyStressFEM<double>
  (
  double* modelCoordinates,
  double* deformedCoordinates, 
  double* unrotatedCauchyStressNP1, 
  int numPoints, 
  int* neighborhoodList,
  const double Cstiff[][6],
  double* angles,
  int type,
  double dt,
  int order[3]
  );







}


