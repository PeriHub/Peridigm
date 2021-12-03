//! \file user_material_interface_correspondence.cxx

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
// Christian Willberg christian.willberg@dlr.de
//
// ************************************************************************
//
// funded by EMMA project reference
// Licence agreement
//
// Christian Willberg    christian.willberg@dlr.de
//@HEADER

#include "user_material_interface_correspondence.h"
#include "user_material.h"
#include "correspondence.h"
#include <Sacado.hpp>
#include <string>
#include "matrices.h"
namespace CORRESPONDENCE {



template<typename ScalarT>
void userMaterialInterface
(
const double* coords,
const ScalarT* DeformationGradientN, 
ScalarT* DeformationGradientNP1, 
const ScalarT* strainN, 
ScalarT* strainNP1, 
const ScalarT* unrotatedCauchyStressN, 
ScalarT* unrotatedCauchyStressNP1, 
const int numPoints, 
const int nstatev,
ScalarT* statev,
const int nprops,
const ScalarT* props,
const double* angles,
const double* flyingPointFlag,
const double time,
const double dtime,
const double* temp,
const double* dtemp,
const double* RotationN,
const double* RotationNP1,
const bool plane_stress,
const bool plane_strain,
const std::string matname
)
{
  // Hooke's law
  const ScalarT* defGradN = DeformationGradientN;
  ScalarT* defGradNP1 = DeformationGradientNP1;
  const ScalarT* GLStrainN = unrotatedCauchyStressN;
  ScalarT* GLStrainNP1 = strainNP1;
  const ScalarT* sigmaN = strainN;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;
  //const std::string matname;
  ScalarT deps[6], drot[9], strain[6], dstrain[9], stress[6];

  CORRESPONDENCE::computeGreenLagrangeStrain(defGradNP1,GLStrainNP1,flyingPointFlag ,numPoints);

  int nshr = 3;
  int nnormal = 3;
  int nstresscomp = 6;
  if (plane_stress) {nstresscomp = 3; nshr = 1; nnormal = 2;}
  if (plane_strain) {nstresscomp = 4; nshr = 1; nnormal = 3;}
  
  int NOEL;

  // not supported
  double PREDEF = -1, DPRED = -1, PNEWDT = -1, CELENT = -1;
  int NPT = -1, KSLAY = -1, KSPT = -1, JSTEP = -1, KINC = -1;
  double DSDDE[6][6];
  double DDSDDT[6], DRPLDE[6], DRPLDT;
  double SSE,SPD,SCD,RPL;
  //
  //  2D re-sort??
  //
  bool rotation = false;
  double* rotMat, rotMatT;
  for(int iID=0 ; iID<numPoints ; ++iID, 
        coords+=3, defGradN+=9, defGradNP1+=9, sigmaN+=9, GLStrainN+=9,GLStrainNP1+=9,sigmaNP1+=9, angles+=3){
          NOEL = iID;
          CORRESPONDENCE::DIFFTENSOR(GLStrainN, GLStrainNP1, dstrain);
          CORRESPONDENCE::DIFFTENSOR(RotationN, RotationNP1, drot);

        //if (rotation){  
        //    CORRESPONDENCE::createRotationMatrix(angles,rotMat);
        //    MATRICES::TransposeMatrix(rotMat,rotMatT);
        //    // geomNL
        //    MATRICES::MatrixMultiply3x3(rotMatT,strain, temp);
        //    MATRICES::MatrixMultiply3x3(temp,rotMat,strain);
        //  }
            CORRESPONDENCE::GETVOIGTNOTATION(GLStrainN,strain);
            CORRESPONDENCE::GETVOIGTNOTATION(dstrain,deps);
          /*
          Rotationstransformation

          CORRESPONDENCE::UMAT(stress,statev,DSDDE,SSE,SPD,SCD,RPL,
          DDSDDT, DRPLDE,DRPLDT,strain,deps,time,dtime,temp,dtemp,
          PREDEF,DPRED,matname,nnormal,nshr,nstresscomp,nstatev,props,
          nprops,coords,drot,PNEWDT,CELENT,defGradN,defGradNP1,
          NOEL,NPT,KSLAY,KSPT,JSTEP,KINC)
           
          Rotationstransformation 
            */  
           CORRESPONDENCE::GETTENSORNOTATION(stress,sigmaNP1);                        

        }

          // rotation back
          // if (rotation){  
          //   MATRICES::MatrixMultiply3x3fromVector(rotMat,sigmaNP1, temp);
          //   MATRICES::MatrixMultiply3x3toVector(temp,rotMatT,sigmaNP1);
          //
          // }

}
void GETVOIGTNOTATION
(
const double* TENSOR,
double VOIGT[6]
)
{
  VOIGT[0] = *(TENSOR);
  VOIGT[1] = *(TENSOR+4);
  VOIGT[2] = *(TENSOR+8);
  VOIGT[3] = 0.5*(*(TENSOR+5) + *(TENSOR+7));
  VOIGT[4] = 0.5*(*(TENSOR+2) + *(TENSOR+6));
  VOIGT[5] = 0.5*(*(TENSOR+1) + *(TENSOR+3));


}

void GETTENSORNOTATION
(
const double VOIGT[6],
double* TENSOR
)
{
  *(TENSOR)   = VOIGT[0];*(TENSOR+1) = VOIGT[5];*(TENSOR+2) = VOIGT[4];
  *(TENSOR+3) = VOIGT[5];*(TENSOR+4) = VOIGT[1];*(TENSOR+5) = VOIGT[3];
  *(TENSOR+6) = VOIGT[4];*(TENSOR+7) = VOIGT[3];*(TENSOR+8) = VOIGT[2];
}

template<typename ScalarT>
void DIFFTENSOR
(
const ScalarT* TENSORN,
const ScalarT* TENSORNP1, 
ScalarT* DTENSOR
)
{
  for(int iID=0 ; iID<9 ; ++iID){
    *(DTENSOR+iID) = *(TENSORNP1+iID)-*(TENSORN+iID);
  }
}

template void DIFFTENSOR<double>
(
const double* TENSORN,
const double* TENSORNP1, 
double* DTENSOR
);

/* 
template void userMaterialInterface<Sacado::Fad::DFad<double>>
(
const Sacado::Fad::DFad<double>* DeformationGradient, 
const Sacado::Fad::DFad<double>* unrotatedCauchyStressN, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1, 
const int numPoints, 
const int nprops,
const Sacado::Fad::DFad<double> props[],
const double* angles,
const int type,
const double dt,
const bool hencky
);
*/
template void userMaterialInterface<double>
(
const double* coords,
const double* DeformationGradientN, 
double* DeformationGradientNP1, 
const double* strainN, 
double* strainNP1, 
const double* unrotatedCauchyStressN, 
double* unrotatedCauchyStressNP1, 
const int numPoints, 
const int nstatev,
double* statev,
const int nprops,
const double* props,
const double* angles,
const double* flyingPointFlag,
const double time,
const double dtime,
const double* temp,
const double* dtemp,
const double* RotationN,
const double* RotationNP1,
const bool plane_stress,
const bool plane_strain,
const std::string matname
);





}


