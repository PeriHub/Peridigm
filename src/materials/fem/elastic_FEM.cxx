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

  //template<typename ScalarT>
  void elasticFEM
  (
  const double* modelCoordinates,
  const double* displacements, 
  double* unrotatedCauchyStressNP1, 
  const int numElements, 
  const int* elementNodalList,
  const double Cstiff[][6],
  const double* angles,
  const int type,
  const double dt,
  const int order[3],
  double* globalForce 
  )
  { 
    bool rotation = false;
    const double* nodalCoor = modelCoordinates;
    const double* disp = displacements;
   // double* force = globalForce;
    double* sigmaNP1 = unrotatedCauchyStressNP1;
    double strain[3][3];
    double rotMat[3][3], rotMatT[3][3], temp[3][3];
    //int defGradLogReturnCode(0);

    std::vector<double> sigmaIntVector(9);
    double* sigmaInt = &sigmaIntVector[0];
    const int* elemNodalPtr = elementNodalList;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // for higher order it has to be adapted
    // it is a full integration
    // reduced integration is not included yet
    /////////////////////////////////////////////////////////////////////
  
    bool twoD = false;
    if (type!=0) {twoD = true;}

    int numInt = FEM::getNumberOfIntegrationPoints(twoD, order);

    int ndof = numInt * 3;
    std::vector<double> dispNodalVector(ndof);
    double* dispNodal = &dispNodalVector[0];
    std::vector<double> elNodalCoorVector(ndof);
    double* elNodalCoor = &elNodalCoorVector[0];
    std::vector<double> elNodalForceVector(ndof);
    double* elNodalForces = &elNodalForceVector[0];
    std::vector<double> NxiVector(order[0]+1), NetaVector(order[1]+1), NpsiVector(order[2]+1);
    double* Nxi = &NxiVector[0];
    double* Neta = &NetaVector[0];
    double* Npsi = &NpsiVector[0];
    std::vector<double> BxiVector(order[0]+1), BetaVector(order[1]+1), BpsiVector(order[2]+1);
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
    int topo[numInt][3];
    std::vector<double> JMat(9);
    double* J = &JMat[0];
    std::vector<double> JinvMat(9);
    double* Jinv = &JinvMat[0];
    double detJ = 0.0;
    int localId, numNeigh;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // irgendwie in die init?
    // das wäre der Austausch für verschiedene Elementtypen
    // falsch
    
    FEM::weightsAndIntegrationPoints(order[0], elCoorx, weightsx);
    FEM::weightsAndIntegrationPoints(order[1], elCoory, weightsy);
    if (type == 0) FEM::weightsAndIntegrationPoints(order[2], elCoorz, weightsz);
    FEM::getElementTopo(twoD, order, topo);
    
    for (int iID=0 ; iID<order[0]+1 ; ++iID){
      FEM::getLagrangeElementData(order[0],elCoorx[iID],Nxi,Bxi);
    }  
    for (int iID=0 ; iID<order[1]+1 ; ++iID){
      FEM::getLagrangeElementData(order[1],elCoory[iID],Neta,Beta);
    }  
    if (twoD == false){
      for (int iID=0 ; iID<order[2]+1 ; ++iID){
        FEM::getLagrangeElementData(order[2],elCoorz[iID],Npsi,Bpsi);
      }  
    }
    for(int iID=0 ; iID<numElements ; ++iID, sigmaNP1+=9, angles+=3){
      numNeigh = *elemNodalPtr;
      elemNodalPtr++;
      for(int n=0 ; n<numNeigh ; ++n){

        localId = *elemNodalPtr;
        elNodalCoor[3*n]   = nodalCoor[localId];
        elNodalCoor[3*n+1] = nodalCoor[localId+1];
        elNodalCoor[3*n+2] = nodalCoor[localId+2];
        dispNodal[3*n]     = disp[localId];
        dispNodal[3*n+1]   = disp[localId+1];
        dispNodal[3*n+2]   = disp[localId+2];
        elNodalForces[3*n]     = 0.0;
        elNodalForces[3*n+1]   = 0.0;
        elNodalForces[3*n+2]   = 0.0;
        elemNodalPtr++;
      }

      elemNodalPtr -= numNeigh; // needed for the second loop to write forces in the global vector
      for (int jID=0 ; jID<numInt ; ++jID){
        // only if nodes and integration points are equal the topology is suitable here.
        
        FEM::getJacobian(Nxi,Neta,Npsi,Bxi,Beta,Bpsi,ndof,topo,elNodalCoor, twoD, J, detJ, Jinv);
        //
        FEM::computeStrain(Nxi,Neta,Npsi,Bxi,Beta,Bpsi,topo,dispNodal, ndof, Jinv, twoD, strain); 
        
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
        FEM::getNodelForce(Nxi,Neta,Npsi,Bxi,Beta,Bpsi,topo, sigmaInt, ndof, abs(detJ), Jinv, twoD, elNodalForces);

              // gemittelte Spannungen
              // sigmaNP1 /= numInt;
              //globForce(topo) += force; ??
      }

      for(int n=0 ; n<numNeigh ; ++n){
        localId = *elemNodalPtr;
        globalForce[localId]   += elNodalForces[3*n];      
        globalForce[localId+1] += elNodalForces[3*n+1];
        globalForce[localId+2] += elNodalForces[3*n+2];
        elemNodalPtr++;
      }



    }
 
  }



}


  