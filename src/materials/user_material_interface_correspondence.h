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
);
template<typename ScalarT>
void DIFFTENSOR
(
const ScalarT* TENSORN,
const ScalarT* TENSORNP1, 
ScalarT* DTENSOR
);
void GETVOIGTNOTATION
(
const double* TENSOR,
double VOIGT[6]
);
void GETTENSORNOTATION
(
const double VOIGT[6],
double* TENSOR
);
}

#endif // USER_MATERIAL_INTERFACE_CORRESPONDENCE_H
