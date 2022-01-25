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
#include "correspondence.h"
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
const std::string matname,
const bool* coordinateTrafo
)
{
  // Hooke's law
  const ScalarT* defGradN = DeformationGradientN;
  ScalarT* defGradNP1 = DeformationGradientNP1;
  const ScalarT* GLStrainN = strainN;
  ScalarT* GLStrainNP1 = strainNP1;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;
  char matnameArray[80];
  int nname;


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

  std::vector<ScalarT> strainLocVector(9),depsVector(9), drotVector(9), depsLocVector(9), sigmaNP1LocVector(9);
  ScalarT* strainLoc = &strainLocVector[0];
  ScalarT* deps = &depsVector[0];
  ScalarT* depsLoc = &depsLocVector[0];
  ScalarT* drot = &drotVector[0];
  ScalarT* sigmaNP1Loc = &sigmaNP1LocVector[0];
   // Voigt Notation
  std::vector<ScalarT> strainLocVoigtVector(6),depsLocVoigtVector(6), sigmaNP1LocVoigtVector(6);
  ScalarT* strainLocVoigt = &strainLocVoigtVector[0];
  ScalarT* depsLocVoigt = &depsLocVoigtVector[0];
  ScalarT* sigmaNP1LocVoigt = &sigmaNP1LocVoigtVector[0];

  for(int iID=0 ; iID<numPoints ; ++iID, 
        coords+=3, defGradN+=9, defGradNP1+=9, GLStrainN+=9,GLStrainNP1+=9,sigmaNP1+=9, angles+=3){
          NOEL = iID;
    
          CORRESPONDENCE::DIFFTENSOR(GLStrainN, GLStrainNP1, deps);
          CORRESPONDENCE::DIFFTENSOR(RotationN, RotationNP1, drot);
          //CORRESPONDENCE::StoreAsMatrix(drotV, drot);
          nname = matname.length();

          for(unsigned int i=0; i<matname.length(); ++i){
            matnameArray[i] = matname[i];
          }

          // Rotationstransformation
          //https://www.continuummechanics.org/stressxforms.html
          // Q Q^T * sigma * Q Q^T = Q C Q^T epsilon Q Q^T
          if (coordinateTrafo[iID]==true){  
            MATRICES::tensorRotation(angles,GLStrainN,true,strainLoc);
            MATRICES::tensorRotation(angles,deps,true,depsLoc);
          }

          CORRESPONDENCE::GetVoigtNotation(strainLoc, strainLocVoigt);
          CORRESPONDENCE::GetVoigtNotation(depsLoc, depsLocVoigt);
          
          *(drot) = 1.0;
          *(drot+1) = 2.0;
          *(drot+2) = 3.0;
          *(drot+3) = 4.0;
          *(drot+4) = 5.0;
          *(drot+5) = 6.0;
          *(drot+6) = 7.0;
          *(drot+7) = 8.0;
          *(drot+8) = 9.0;
          
          std::cout<<"*(drot)"<<*(drot)<<std::endl;
          std::cout<<"*(drot+1)"<<*(drot+1)<<std::endl;
          std::cout<<"*(drot+2)"<<*(drot+2)<<std::endl;
          std::cout<<"*(drot+3)"<<*(drot+3)<<std::endl;
          std::cout<<"*(drot+4)"<<*(drot+4)<<std::endl;
          std::cout<<"*(drot+5)"<<*(drot+5)<<std::endl;
          std::cout<<"*(drot+6)"<<*(drot+6)<<std::endl;
          std::cout<<"*(drot+7)"<<*(drot+7)<<std::endl;
          std::cout<<"*(drot+8)"<<*(drot+8)<<std::endl;

          CORRESPONDENCE::UMATINT(sigmaNP1LocVoigt,statev,DSDDE,&SSE,&SPD,&SCD,&RPL,
          DDSDDT, DRPLDE,&DRPLDT,strainLocVoigt,depsLocVoigt,&time,&dtime,temp,dtemp,
          &PREDEF,&DPRED,matnameArray,&nnormal,&nshr,&nstresscomp,&nstatev,props,
          &nprops,coords,drot,&PNEWDT,&CELENT,defGradN,defGradNP1,
          &NOEL,&NPT,&KSLAY,&KSPT,&JSTEP,&KINC,&nname); 

          std::cout<<"*(sigmaNP1LocVoigt)"<<*(sigmaNP1LocVoigt)<<std::endl;
          std::cout<<"*(sigmaNP1LocVoigt+1)"<<*(sigmaNP1LocVoigt+1)<<std::endl;
          std::cout<<"*(sigmaNP1LocVoigt+2)"<<*(sigmaNP1LocVoigt+2)<<std::endl;
          std::cout<<"*(sigmaNP1LocVoigt+3)"<<*(sigmaNP1LocVoigt+3)<<std::endl;
          std::cout<<"*(sigmaNP1LocVoigt+4)"<<*(sigmaNP1LocVoigt+4)<<std::endl;
          std::cout<<"*(sigmaNP1LocVoigt+5)"<<*(sigmaNP1LocVoigt+5)<<std::endl;
          std::cout<<"*(sigmaNP1LocVoigt+6)"<<*(sigmaNP1LocVoigt+6)<<std::endl;
          std::cout<<"*(sigmaNP1LocVoigt+7)"<<*(sigmaNP1LocVoigt+7)<<std::endl;
          std::cout<<"*(sigmaNP1LocVoigt+8)"<<*(sigmaNP1LocVoigt+8)<<std::endl;

          CORRESPONDENCE::GetTensorFromVoigtNotation(sigmaNP1LocVoigt, sigmaNP1Loc);

          // Rotationstransformation 
          if (coordinateTrafo[iID]==true){  
            MATRICES::tensorRotation(angles,sigmaNP1Loc,false,sigmaNP1);
          }

        }

}
//template<typename double>
void StoreAsMatrix
(
  const double* vector,
  double matrix[][3]
)
{
  matrix[0][0] = *(vector);
  matrix[0][1] = *(vector+1);
  matrix[0][2] = *(vector+2);
  matrix[1][0] = *(vector+3);
  matrix[1][1] = *(vector+4);
  matrix[1][2] = *(vector+5);
  matrix[2][0] = *(vector+6);
  matrix[2][1] = *(vector+7);
  matrix[2][2] = *(vector+8);
}


//template void StoreAsMatrix<double>
//(
//  const double* vector,
//  double matrix[][3]
//);

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
const std::string matname,
const bool* coordinateTrafo
);





}


