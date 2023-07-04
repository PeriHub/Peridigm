//! \file user_material_interface_correspondence.h

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
// funded by EMMA project reference
// Licence agreement
//
// Christian Willberg    christian.willberg@dlr.de
//@HEADER
#ifndef USER_MATERIAL_INTERFACE_CORRESPONDENCE_H
#define USER_MATERIAL_INTERFACE_CORRESPONDENCE_H
#include <string>
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
int* step,
const int numPoints,
const int nstatev,
ScalarT* statevVector,
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
const bool testing
);

template<typename ScalarT>
void DIFFTENSOR
(
const ScalarT* TENSORN,
const ScalarT* TENSORNP1, 
ScalarT* DTENSOR
);

template<typename ScalarT>
void GetVoigtNotation
(
const ScalarT* TENSOR,
ScalarT* VOIGT
);

template<typename ScalarT>
void GetTensorFromVoigtNotation
(
const ScalarT* VOIGT,
ScalarT* TENSOR
);

template<typename ScalarT>
void ReduceComp
(
ScalarT* VOIGT,
bool threeOrFour
);

template<typename ScalarT>
void ExtendToSixComp
(
ScalarT* VOIGT,
bool threeOrFour
);

extern "C" void UPERMAT(double stressnew[], double statenew[], 
						double stressold[], double stateold[], 
						double straininc[], const double props[],
						const double *steptime, const double *totaltime, const double *dt, 
						int *NodeID, const int *MatID, 
						const int *NPOINT, const int *NDI, const int *NSHR, 
						const int *NSTATV, const int *NPROPS);

/*
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
*/   

/*
extern "C" void UMAT(double *sigmaNP1, double *statev, double DSDDE[6][6], double SSE, double SPD, double SCD, double RPL,
          double DDSDDT[6], double DRPLDE[6], double DRPLDT, const double *GLStrainN, double deps[9], const double time, const double dtime, const double* temp, const double* dtemp,
          double PREDEF, double DPRED, int nnormal, int nshr, int nstresscomp, const int nstatev, const double *props,
          const int nprops, const double* coords, double drot[3][3], double PNEWDT, double CELENT, const double* defGradN, const double* defGradNP1,
          int NOEL, int NPT, int KSLAY, int KSPT, int JSTEP, int KINC); //const std::string matname,
*/  

// extern "C" void UMAT(double *sigmaNP1, int *NSTRESS); //const std::string matname,

/*
extern "C" void UMAT(double *stress, double *statev, double *ddsdde, double *sse, double *spd,
		double *scd, double *rpl, double *ddsddt, double *drplde, double *drpldt,
		const double *GLStrainN, double *deps, const double *time, const double *dtime, const double *temp, const double *dtemp,
          double *PREDEF, double *DPRED, int *nnormal, int *nshr, int *nstresscomp, const int *nstatev, const double *props,
          const int *nprops, const double *coords, double *drot, double *PNEWDT, double *CELENT, const double *defGradN, const double *defGradNP1,
          int *NOEL, int *NPT, int *KSLAY, int *KSPT, int *JSTEP, int *KINC);
*/ 
extern "C" void UMATINT(double *stress, double *statev, double *ddsdde, double *sse, double *spd,
		double *scd, double *rpl, double *ddsddt, double *drplde, double *drpldt,
		const double *GLStrainN, double *deps, const double *time, const double *dtime, const double *temp, const double *dtemp,
          double *PREDEF, double *DPRED, char *matnameArray, int *nnormal, int *nshr, int *nstresscomp, const int *nstatev, const double *props,
          const int *nprops, const double *coords, double *drot, double *PNEWDT, double *CELENT, const double *defGradN, const double *defGradNP1,
          int *NOEL, int *NPT, int *KSLAY, int *KSPT, int *JSTEP, int *KINC, int *nname);

extern "C" void UMATINTTEST(double *stress, double *statev, double *ddsdde, double *sse, double *spd,
		double *scd, double *rpl, double *ddsddt, double *drplde, double *drpldt,
		const double *GLStrainN, double *deps, const double *time, const double *dtime, const double *temp, const double *dtemp,
          double *PREDEF, double *DPRED, char *matnameArray, int *nnormal, int *nshr, int *nstresscomp, const int *nstatev, const double *props,
          const int *nprops, const double *coords, double *drot, double *PNEWDT, double *CELENT, const double *defGradN, const double *defGradNP1,
          int *NOEL, int *NPT, int *KSLAY, int *KSPT, int *JSTEP, int *KINC, int *nname);
}

#endif // USER_MATERIAL_INTERFACE_CORRESPONDENCE_H
