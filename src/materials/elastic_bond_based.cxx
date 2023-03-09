//! \file elastic_bond_based.cxx

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

#include <cmath>
#include <Sacado.hpp>
#include "elastic_bond_based.h"
#include "material_utilities.h"

namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void computeInternalForceElasticBondBased
(
    const double* xOverlap,
    const ScalarT* yOverlap,
    const double* volumeOverlap,
    const double* bondDamage,
    ScalarT* fInternalOverlap,
    const int* localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double horizon
)
{
  const double *xOwned = xOverlap;
  const ScalarT *yOwned = yOverlap;
  double volume, neighborVolume, damageOnBond;
  ScalarT stretch, t, zeta, dY;
  int neighborhoodIndex(0), bondDamageIndex(0), neighborId;

  const double pi = PeridigmNS::value_of_pi();
  double constant = 18.0*BULK_MODULUS/(pi*horizon*horizon*horizon*horizon);
  const int dof = PeridigmNS::dof();
  std::vector<double> X_dxVector(dof)  ; double*  X_dx = &X_dxVector[0];
  std::vector<ScalarT> Y_dxVector(dof) ; ScalarT*  Y_dx = &Y_dxVector[0];
  std::vector<ScalarT> f_Vector(dof) ; ScalarT*  f = &f_Vector[0];
  for(int p=0 ; p<numOwnedPoints ; p++){

    const double *X = xOwned;
    const ScalarT *Y = yOwned;
    volume = volumeOverlap[p];

    int numNeighbors = localNeighborList[neighborhoodIndex++];
    for(int n=0; n<numNeighbors; n++, xOwned +=3, yOwned +=3){

      neighborId = localNeighborList[neighborhoodIndex++];
      const double *XP = &xOverlap[3*neighborId];
      const ScalarT *YP = &yOverlap[3*neighborId];
      neighborVolume = volumeOverlap[neighborId];
      
      zeta = MATERIAL_EVALUATION::getDiffAndLen(X,XP,dof,X_dx);
      dY = MATERIAL_EVALUATION::getDiffAndLen(Y,YP,dof,Y_dx);
      stretch = (dY - zeta)/zeta;

      damageOnBond = bondDamage[bondDamageIndex++];

      t = 0.5*(1.0 - damageOnBond)*stretch*constant;

      MATERIAL_EVALUATION::getProjectedForces(t,Y_dx,dY,dof,f);
      MATERIAL_EVALUATION::setForces(f[0], f[1], f[2], volume, neighborVolume, &fInternalOverlap[3*p], &fInternalOverlap[3*neighborId]);

    }
  }
}

/** Explicit template instantiation for double. */
template void computeInternalForceElasticBondBased<double>
(
    const double* xOverlap,
    const double* yOverlap,
    const double* volumeOverlap,
    const double* bondDamage,
    double* fInternalOverlap,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double horizon
 );

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceElasticBondBased<Sacado::Fad::DFad<double> >
(
    const double* xOverlap,
    const Sacado::Fad::DFad<double>* yOverlap,
    const double* volumeOverlap,
    const double* bondDamage,
    Sacado::Fad::DFad<double>* fInternalOverlap,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double horizon
);

}
