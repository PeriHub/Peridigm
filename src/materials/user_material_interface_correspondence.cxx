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
#include "correspondence.h"
#include "matrices.h"
#include <Sacado.hpp>
#include <string>
#include "matrices.h"
namespace CORRESPONDENCE {



template<typename ScalarT>
void userMaterialInterface
(
const double* coords,
const ScalarT* DeformationGradientN, 
const ScalarT* DeformationGradientNP1, 
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
  const ScalarT* defGradNP1 = DeformationGradientNP1;
  const ScalarT* GLStrainN = strainN;
  ScalarT* GLStrainNP1 = strainNP1;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;
  const ScalarT* sigmaN = unrotatedCauchyStressN;
  
  bool rotation = false;
  int nshr = 3;
  int nnormal = 3;
  int nstresscomp = 6;
  if (plane_stress) {nstresscomp = 3; nshr = 1; nnormal = 2;}
  if (plane_strain) {nstresscomp = 4; nshr = 1; nnormal = 3;}
  
  int NOEL;
  const double timeArray[2] = {time,time};

  // not supported
  double PREDEF = -1, DPRED = -1, PNEWDT = -1, CELENT = -1;
  int NPT = -1, KSLAY = -1, KSPT = -1, JSTEP = -1, KINC = -1;
  double DDSDDE[6*6];
  double DDSDDT[6], DRPLDE[6];
  // double* DDSDDT;
  // double* DRPLDE;
  double DRPLDT = -1;
  double SSE = -1,SPD = -1,SCD = -1,RPL = -1;
  //

  std::vector<ScalarT> strainLocVector(9),depsVector(9), drotVector(9), depsLocVector(9);
  ScalarT* strainLoc = &strainLocVector[0];
  ScalarT* deps = &depsVector[0];
  ScalarT* drot = &drotVector[0];
   // Voigt Notation
  std::vector<ScalarT> strainLocVoigtVector(6),depsLocVoigtVector(6), sigmaNP1LocVoigtVector(6);
  ScalarT* strainLocVoigt = &strainLocVoigtVector[0];
  ScalarT* depsLocVoigt = &depsLocVoigtVector[0];
  ScalarT* sigmaNP1LocVoigt = &sigmaNP1LocVoigtVector[0];

  char matnameArray[80];
  int nname = matname.length();
  for(unsigned int i=0; i<matname.length(); ++i){
    matnameArray[i] = matname[i];
  }

  for(int iID=0 ; iID<numPoints ; ++iID, 
          statev+=nstatev, coords+=3, defGradN+=9, defGradNP1+=9, GLStrainN+=9,GLStrainNP1+=9,sigmaN+=9,sigmaNP1+=9, angles+=3){
          NOEL = iID;
    if (MATRICES::vectorNorm(angles, 3)!=0)rotation=true;
    else rotation = false;
    CORRESPONDENCE::computeGreenLagrangeStrain(defGradNP1,GLStrainNP1);

    CORRESPONDENCE::DIFFTENSOR(GLStrainN, GLStrainNP1, deps);
    CORRESPONDENCE::DIFFTENSOR(RotationN, RotationNP1, drot);
    // Transformation global -> local
    //https://www.continuummechanics.org/stressxforms.html
    // Q Q^T * sigma * Q Q^T = Q C Q^T epsilon Q Q^T
    if (rotation==true){  
      MATRICES::tensorRotation(angles,GLStrainN,true,strainLoc);
      MATRICES::tensorRotation(angles,deps,true,deps);
      MATRICES::tensorRotation(angles,sigmaN,true,sigmaNP1);
// old stresses N -> hier rein

    }
    else{
      for(int jID=0 ; jID<9 ; ++jID){
        *(strainLoc+jID)=*(GLStrainN+jID);
        *(sigmaNP1+jID)=*(sigmaN+jID);
      }
    }
    
    
    CORRESPONDENCE::GetVoigtNotation(strainLoc, strainLocVoigt);
    CORRESPONDENCE::GetVoigtNotation(deps, depsLocVoigt);
    CORRESPONDENCE::GetVoigtNotation(sigmaNP1LocVoigt, sigmaNP1);
    if (plane_stress){
      CORRESPONDENCE::ReduceComp(sigmaNP1LocVoigt, true);
      CORRESPONDENCE::ReduceComp(strainLocVoigt, true);
      CORRESPONDENCE::ReduceComp(depsLocVoigt, true);
    }
    if (plane_strain){
      CORRESPONDENCE::ReduceComp(sigmaNP1LocVoigt, false);
      CORRESPONDENCE::ReduceComp(strainLocVoigt, false);
      CORRESPONDENCE::ReduceComp(depsLocVoigt, false);
    }
    // std::cout<<"Run UMAT"<<std::endl;
    CORRESPONDENCE::UMATINT(sigmaNP1LocVoigt,statev,DDSDDE,&SSE,&SPD,&SCD,&RPL,
    DDSDDT, DRPLDE,&DRPLDT,strainLocVoigt,depsLocVoigt,timeArray,&dtime,temp,dtemp,
    &PREDEF,&DPRED,matnameArray,&nnormal,&nshr,&nstresscomp,&nstatev,props,
    &nprops,coords,drot,&PNEWDT,&CELENT,defGradN,defGradNP1,
    &NOEL,&NPT,&KSLAY,&KSPT,&JSTEP,&KINC,&nname);
    // std::cout<<"Umat finished"<<std::endl;

    if (plane_stress){
      CORRESPONDENCE::ExtendToSixComp(sigmaNP1LocVoigt, true);
      CORRESPONDENCE::ExtendToSixComp(strainLocVoigt, true);
      CORRESPONDENCE::ExtendToSixComp(depsLocVoigt, true);
    }
    if (plane_strain){
      CORRESPONDENCE::ExtendToSixComp(sigmaNP1LocVoigt, false);
      CORRESPONDENCE::ExtendToSixComp(strainLocVoigt, false);
      CORRESPONDENCE::ExtendToSixComp(depsLocVoigt, false);
    }

    CORRESPONDENCE::GetTensorFromVoigtNotation(sigmaNP1LocVoigt, sigmaNP1);
    // back transformation local -> global 
    if (rotation==true){  
      MATRICES::tensorRotation(angles,sigmaNP1,false,sigmaNP1);
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
  MATRICES::DIFFTENSOR(TENSORN,TENSORNP1,DTENSOR);
}

template void DIFFTENSOR<double>
(
const double* TENSORN,
const double* TENSORNP1, 
double* DTENSOR
);


template<typename ScalarT>
void GetVoigtNotation
(
const ScalarT* TENSOR,
ScalarT* VOIGT
)
{
  *(VOIGT) = *(TENSOR);
  *(VOIGT+1) = *(TENSOR+4);
  *(VOIGT+2) = *(TENSOR+8);
  *(VOIGT+3) = 0.5*(*(TENSOR+5) + *(TENSOR+7));
  *(VOIGT+4) = 0.5*(*(TENSOR+2) + *(TENSOR+6));
  *(VOIGT+5) = 0.5*(*(TENSOR+1) + *(TENSOR+3));
}

template void GetVoigtNotation<double>
(
const double* TENSOR,
double* VOIGT
);

template<typename ScalarT>
void GetTensorFromVoigtNotation
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

template void GetTensorFromVoigtNotation<double>
(
const double* VOIGT,
double* TENSOR
);

template<typename ScalarT>
void ReduceComp
(
ScalarT* VOIGT,
bool threeOrFour
)
{
  if (threeOrFour){
    *(VOIGT) = *(VOIGT);
    *(VOIGT+1) = *(VOIGT+1);
    *(VOIGT+2) = *(VOIGT+5);
  }else{
    *(VOIGT) = *(VOIGT);
    *(VOIGT+1) = *(VOIGT+1);
    *(VOIGT+2) = *(VOIGT+2);
    *(VOIGT+3) = *(VOIGT+5);
  }
}

template void ReduceComp<double>
(
double* VOIGT,
bool threeOrFour
);

template<typename ScalarT>
void ExtendToSixComp
(
ScalarT* VOIGT,
bool threeOrFour
)
{
  if (threeOrFour){
    *(VOIGT) = *(VOIGT);
    *(VOIGT+1) = *(VOIGT+1);
    *(VOIGT+5) = *(VOIGT+2);
    *(VOIGT+2) = 0.0;
    *(VOIGT+3) = 0.0;
    *(VOIGT+4) = 0.0;
  }else{
    *(VOIGT) = *(VOIGT);
    *(VOIGT+1) = *(VOIGT+1);
    *(VOIGT+2) = *(VOIGT+2);
    *(VOIGT+5) = *(VOIGT+3);
    *(VOIGT+3) = 0.0;
    *(VOIGT+4) = 0.0;
  }
}

template void ExtendToSixComp<double>
(
double* VOIGT,
bool threeOrFour
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
const double* DeformationGradientNP1, 
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


