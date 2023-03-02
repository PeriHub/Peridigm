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
    ScalarT ed;
    double gamma=0;

    const int *neighPtr = localNeighborList;
    double cellVolume, alpha, X_dx[3], zeta = 0.0, omega;
    const int dof = 3;
    ScalarT Y_dx[3], dY = 0.0, t=0.0, fx, fy, fz, e, c1=0;
    for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, psOwned+=9, deltaT++, m++, theta++){

        int numNeigh = *neighPtr; neighPtr++;
        const double *X = xOwned;
        const ScalarT *Y = yOwned;
        alpha = 15.0*MU/(*m);
        double selfCellVolume = v[p];
        // c1 = omega*(*theta)*(9.0*K-15.0*MU)/(3.0*(*m));
        // Bobaru et al. "Handbook of Peridynamics Modeling; page 152 ff.;
        if (planeStress==true){
            c1 = 4.0*K*MU/(3.0*K+4.0*MU) * (*theta)/(*m);
            gamma = 4.0*MU/(3.0*K+4.0*MU);
        }
        // Bobaru et al. "Handbook of Peridynamics Modeling; page 152 ff.;
        if (planeStrain==true){
            c1 = (12.0*K-4.0*MU) / 9.0 * (*theta)/(*m);
            gamma = 2/3;
        }
        for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
            int localId = *neighPtr;
            cellVolume = v[localId];
            const double *XP = &xOverlap[3*localId];
            const ScalarT *YP = &yOverlap[3*localId];
            
            MATERIAL_EVALUATION::getDiffAndLen(XP,X,dof,X_dx,zeta);
            MATERIAL_EVALUATION::getDiffAndLen(YP,Y,dof,Y_dx,dY);
            
            e = dY - zeta;
            
            
            if(deltaTemperature)
              e -= thermalExpansionCoefficient*(*deltaT)*zeta;
            omega = scalarInfluenceFunction(zeta,horizon);
            //if (e == 0) continue;
            if (planeStrain or planeStress){            
                ed = e - gamma * (*theta)* zeta / 3;
                t = (1.0-*bondDamage)*omega*(c1 * zeta +  8.0 * MU * ed / (*m));
            }
            else {
                c1 = (*theta)*(3.0*K/(*m)-alpha/3.0);
                t = (1.0-*bondDamage)*omega*(c1 * zeta + alpha * e);
            }
            fx = t * Y_dx[0] / dY;
            fy = t * Y_dx[1] / dY;
            fz = t * Y_dx[2] / dY;

            MATERIAL_EVALUATION::setForces(fx, fy, fz, selfCellVolume, cellVolume, fOwned, &fInternalOverlap[3 * localId]);
            if(partialStressOverlap != 0){
              MATERIAL_EVALUATION::setPartialStresses(fx, fy, fz, X_dx[0], X_dx[1], X_dx[2], cellVolume, psOwned);
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
