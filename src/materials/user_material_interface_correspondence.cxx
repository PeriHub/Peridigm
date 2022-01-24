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
//@HEADER

#include "user_material_interface_correspondence.h"
#include "user_material.h"
#include "correspondence.h"
#include "elastic_correspondence.h"
#include "matrices.h"
#include <Sacado.hpp>
#include <string>

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
const std::string matname,
bool hencky
)
{
  // Hooke's law
  const ScalarT* defGradN = DeformationGradientN;
  ScalarT* defGradNP1 = DeformationGradientNP1;
  const ScalarT* GLStrainN = unrotatedCauchyStressN;
  ScalarT* GLStrainNP1 = strainNP1;
  const ScalarT* sigmaN = strainN;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;
  ScalarT* sigmaNP1Voigt;
  ScalarT* GLStrainNVoigt;
  double* depsVoigt;
  char matnameArray[80];
  int nname;
  double deps[9], drot[9];
  // double* deps; 
  // double* drot;

  // CORRESPONDENCE::computeGreenLagrangeStrain(defGradNP1,GLStrainNP1,flyingPointFlag ,numPoints);

  int nshr = 3;
  int nnormal = 3;
  int nstresscomp = 6;
  if (plane_stress) {nstresscomp = 3; nshr = 1; nnormal = 2;}
  if (plane_strain) {nstresscomp = 4; nshr = 1; nnormal = 3;}
  
  int NOEL;

  // not supported
  double PREDEF = -1, DPRED = -1, PNEWDT = -1, CELENT = -1;
  int NPT = -1, KSLAY = -1, KSPT = -1, JSTEP = -1, KINC = -1;
  double DSDDE[6*6];
  double DDSDDT[6], DRPLDE[6];
  // double* DDSDDT;
  // double* DRPLDE;
  double DRPLDT = -1;
  double SSE = -1,SPD = -1,SCD = -1,RPL = -1;
  //
  
  int defGradLogReturnCode(0);
  bool rotation = true;
  ScalarT strain[3][3];
  ScalarT rotMat[3][3], rotMatT[3][3], tempA[3][3];
  
  for(int iID=0 ; iID<numPoints ; ++iID, 
        coords+=3, defGradN+=9, defGradNP1+=9, sigmaN+=9, GLStrainN+=9,GLStrainNP1+=9,sigmaNP1+=9, angles+=3){
          NOEL = iID;

          CORRESPONDENCE::computeGreenLagrangeStrain(defGradNP1,GLStrainNP1,flyingPointFlag ,numPoints);
          
          std::string logStrainErrorMessage =
            "**** Error:  CorrespondenceMaterialconst ::updateElasticCauchyStressAnisotropic() failed to compute LogStrain.\n";
          TEUCHOS_TEST_FOR_TERMINATION(defGradLogReturnCode != 0, logStrainErrorMessage);

          CORRESPONDENCE::DIFFTENSOR(GLStrainN, GLStrainNP1, deps);
          CORRESPONDENCE::DIFFTENSOR(RotationN, RotationNP1, drot);

          nname = matname.length();

          for(unsigned int i=0; i<matname.length(); ++i)
          {
            matnameArray[i] = matname[i];
          }

          // Rotationstransformation
          //https://www.continuummechanics.org/stressxforms.html
          // Q Q^T * sigma * Q Q^T = Q C Q^T epsilon Q Q^T
          if (rotation){  
            CORRESPONDENCE::createRotationMatrix(angles,rotMat);
            MATRICES::TransposeMatrix(rotMat,rotMatT);
            // geomNL
            MATRICES::MatrixMultiply3x3fromVector(rotMatT,GLStrainN, tempA);
            MATRICES::MatrixMultiply3x3toVector(tempA,rotMat,GLStrainN);
          }

          CORRESPONDENCE::GetVoigt(sigmaNP1, sigmaNP1Voigt);
          CORRESPONDENCE::GetVoigt(GLStrainN, GLStrainNVoigt);
          CORRESPONDENCE::GetVoigt(deps, depsVoigt);

          CORRESPONDENCE::UMATINT(sigmaNP1Voigt,statev,DSDDE,&SSE,&SPD,&SCD,&RPL,
          DDSDDT, DRPLDE,&DRPLDT,GLStrainNVoigt,depsVoigt,&time,&dtime,temp,dtemp,
          &PREDEF,&DPRED,matnameArray,&nnormal,&nshr,&nstresscomp,&nstatev,props,
          &nprops,coords,drot,&PNEWDT,&CELENT,defGradN,defGradNP1,
          &NOEL,&NPT,&KSLAY,&KSPT,&JSTEP,&KINC,&nname); 

          CORRESPONDENCE::GetTensorFromVoigt(sigmaNP1Voigt, sigmaNP1);

          // Rotationstransformation 
          if (rotation){  
            MATRICES::MatrixMultiply3x3fromVector(rotMat,sigmaNP1, tempA);
            MATRICES::MatrixMultiply3x3toVector(tempA,rotMatT,sigmaNP1);
          }

        }

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

template<typename ScalarT>
void GetVoigt
(
const ScalarT* TENSOR,
ScalarT* VOIGT
)
{
  *(VOIGT) = *(TENSOR);
  *(VOIGT+1) = *(TENSOR+4);
  *(VOIGT+2) = *(TENSOR+8);
  *(VOIGT+3) = *(TENSOR+5);
  *(VOIGT+4) = *(TENSOR+2);
  *(VOIGT+5) = *(TENSOR+1);
}

template void GetVoigt<double>
(
const double* TENSOR,
double* VOIGT
);

template<typename ScalarT>
void GetTensorFromVoigt
(
const ScalarT* VOIGT,
ScalarT* TENSOR
)
{
  *(TENSOR) = *(VOIGT);
  *(TENSOR+1) = *(VOIGT+5);
  *(TENSOR+2) = *(VOIGT+4);
  *(TENSOR+3) = *(VOIGT+5);
  *(TENSOR+4) = *(VOIGT+1);
  *(TENSOR+5) = *(VOIGT+3);
  *(TENSOR+6) = *(VOIGT+4);
  *(TENSOR+7) = *(VOIGT+3);
  *(TENSOR+8) = *(VOIGT+2);
}

template void GetTensorFromVoigt<double>
(
const double* VOIGT,
double* TENSOR
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
const std::string matname,
bool hencky
);





}


