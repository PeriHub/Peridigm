//! \file elastic.cxx

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
#include "elastic.h"
#include "material_utilities.h"

namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void computeInternalForceLinearElastic
(
        const double* xOverlap,
        const ScalarT* yOverlap,
        const double* mOwned,
        const double* volumeOverlap,
        const ScalarT* dilatationOwned,
        const double* bondDamage,
        ScalarT* fInternalOverlap,
        ScalarT* partialStressOverlap,
        const int*  localNeighborList,
        int numOwnedPoints,
        double BULK_MODULUS,
        double SHEAR_MODULUS,
        double horizon,
        double thermalExpansionCoefficient,
        const double* deltaTemperature,
        const bool planeStrain,
        const bool planeStress
)
{

    /*
     * Compute processor local contribution to internal force
     */
    double K = BULK_MODULUS;
  double MU = SHEAR_MODULUS;

  const double *xOwned = xOverlap;
  const ScalarT *yOwned = yOverlap;
  const double *deltaT = deltaTemperature;
  const double *m = mOwned;
  const double *v = volumeOverlap;
  const ScalarT *theta = dilatationOwned;
  ScalarT *fOwned = fInternalOverlap;
  ScalarT *psOwned = partialStressOverlap;

  const int *neighPtr = localNeighborList;
  double cellVolume, alpha, zeta, omega;
  ScalarT dY, t, e, c1;
  const int dof = PeridigmNS::dof();
  std::vector<double> X_dxVector(dof)  ; double*  X_dx = &X_dxVector[0];
  std::vector<ScalarT> Y_dxVector(dof) ; ScalarT*  Y_dx = &Y_dxVector[0];
  std::vector<ScalarT> f_Vector(dof) ; ScalarT*  f = &f_Vector[0];


  for(int p=0;p<numOwnedPoints;p++, xOwned+=dof, yOwned+=dof, fOwned+=dof, psOwned+=9, deltaT++, m++, theta++){

    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    const ScalarT *Y = yOwned;
    alpha = 15.0*MU/(*m);
    double selfCellVolume = v[p];
    for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
      int localId = *neighPtr;
      cellVolume = v[localId];
      const double *XP = &xOverlap[dof*localId];
      const ScalarT *YP = &yOverlap[dof*localId];

      zeta = MATERIAL_EVALUATION::getDiffAndLen(X,XP,dof,X_dx);
      dY = MATERIAL_EVALUATION::getDiffAndLen(Y,YP,dof,Y_dx);
      e = MATERIAL_EVALUATION::getStretch(zeta,dY);
      if(deltaTemperature)
        e -= thermalExpansionCoefficient*(*deltaT)*zeta;
      omega = scalarInfluenceFunction(zeta,horizon);
      // c1 = omega*(*theta)*(9.0*K-15.0*MU)/(3.0*(*m));
      c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0);
      t = (1.0-*bondDamage)*(c1 * zeta + (1.0-*bondDamage) * omega * alpha * e);

      MATERIAL_EVALUATION::getProjectedForces(t,Y_dx,dY,dof,f);

      MATERIAL_EVALUATION::setForces(f[0], f[1], f[2], selfCellVolume, cellVolume, fOwned, &fInternalOverlap[dof * localId]);
      if(partialStressOverlap != 0){
        MATERIAL_EVALUATION::setPartialStresses(f[0], f[1], f[2], X_dx[0], X_dx[1], X_dx[2], cellVolume, psOwned);
      }
    }
    
  }
}

/** Explicit template instantiation for double. */
template void computeInternalForceLinearElastic<double>
(
        const double* xOverlap,
        const double* yOverlap,
        const double* mOwned,
        const double* volumeOverlap,
        const double* dilatationOwned,
        const double* bondDamage,
        double* fInternalOverlap,
        double* partialStressOverlap,
        const int*  localNeighborList,
        int numOwnedPoints,
        double BULK_MODULUS,
        double SHEAR_MODULUS,
        double horizon,
        double thermalExpansionCoefficient,
        const double* deltaTemperature,
        const bool planeStrain,
        const bool planeStress
 );

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceLinearElastic<Sacado::Fad::DFad<double> >
(
        const double* xOverlap,
        const Sacado::Fad::DFad<double>* yOverlap,
        const double* mOwned,
        const double* volumeOverlap,
        const Sacado::Fad::DFad<double>* dilatationOwned,
        const double* bondDamage,
        Sacado::Fad::DFad<double>* fInternalOverlap,
        Sacado::Fad::DFad<double>* partialStressOverlap,
        const int*  localNeighborList,
        int numOwnedPoints,
        double BULK_MODULUS,
        double SHEAR_MODULUS,
        double horizon,
        double thermalExpansionCoefficient,
        const double* deltaTemperature,
        const bool planeStrain,
        const bool planeStress
);

}
