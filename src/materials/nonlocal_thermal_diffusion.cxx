//! \file nonlocal_thermal_diffusion.cxx

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
#include <algorithm>
#include <vector>
#include <Sacado.hpp>
#include "nonlocal_thermal_diffusion.h"
#include "material_utilities.h"

// ----------------------------------HEAT FLOW---------------------------------
namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void computeHeatFlow
(
	const double*  xOverlap,
	const ScalarT* yOverlap,
	const double* volumeOverlap,
	const double* bondDamage,
	ScalarT* heatFlowOverlap,
	const int*  localNeighborList,
	int numOwnedPoints,
	double thermalConductivity,
	double specificHeat,
	double horizon,
	ScalarT* deltaTemperatureOverlap
)
{
	/*
	 * Compute processor local contribution to internal heat flux
	 */
	double K_T = thermalConductivity;
	const double PI_G = PeridigmNS::value_of_pi();
// 	const double *xOwned = xOverlap;
	const ScalarT *yOwned = yOverlap;
	const ScalarT *deltaTemperatureOwned = deltaTemperatureOverlap;
	const double *v = volumeOverlap;
	ScalarT *heatFlowOwned = heatFlowOverlap;
	const int *neighPtr = localNeighborList;
	double cellVolume;
// 	double X_dx, X_dy, X_dz, zeta, omega;
	ScalarT Y_dx, Y_dy, Y_dz, dY, dT, q1;
	double microConductivity = 6 * K_T /( PI_G * horizon*horizon*horizon*horizon);

// 	loop over all the nodes
	for(int p=0;p<numOwnedPoints;p++, deltaTemperatureOwned++, yOwned +=3, heatFlowOwned++){
		int numNeigh = *neighPtr; neighPtr++;
// 		const double *X = xOwned;
		const ScalarT *Y = yOwned;
// 		double selfCellVolume = v[p];
		const ScalarT *deltaT = deltaTemperatureOwned;
// 		loop over the horizon region
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
// 			const double *XP = &xOverlap[3*localId];
			const ScalarT *deltaTP = &deltaTemperatureOverlap[localId];
			const ScalarT *YP = &yOverlap[3*localId];
// 			X_dx = XP[0]-X[0];
// 			X_dy = XP[1]-X[1];
// 			X_dz = XP[2]-X[2];
// 			zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
			Y_dx = YP[0]-Y[0];
			Y_dy = YP[1]-Y[1];
			Y_dz = YP[2]-Y[2];
			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
// 			omega = scalarInfluenceFunction(zeta,horizon);
			dT = *deltaTP - *deltaT;
			q1 = (1-*bondDamage)*microConductivity*dT/dY;//*omega;
			*heatFlowOwned += q1*cellVolume*horizon;
// 			heatFlowOverlap[localId] -= q1*selfCellVolume*horizon;
		}
	}
}

/** Explicit template instantiation for double. */
template void computeHeatFlow<double>
(
	const double*  xOverlap,
	const double* yOverlap,
	const double* volumeOverlap,
	const double* bondDamage,
	double* heatFlowOverlap,
	const int*  localNeighborList,
	int numOwnedPoints,
// 	std::vector<int> neighPtrVector,
	double thermalConductivity,
	double specificHeat,
	double horizon,
	double* deltaTemperatureOverlap
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeHeatFlow<Sacado::Fad::DFad<double> >
(
	const double*  xOverlap,
	const Sacado::Fad::DFad<double>* yOverlap,
	const double* volumeOverlap,
	const double* bondDamage,
	Sacado::Fad::DFad<double>* heatFlowOverlap,
	const int*  localNeighborList,
	int numOwnedPoints,
	double thermalConductivity,
	double specificHeat,
	double horizon,
	Sacado::Fad::DFad<double>* deltaTemperatureOverlap
);

// --------------------------------INTERNAL FORCE-------------------------------

//! Computes contributions to the internal force resulting from owned points.
template<typename ScalarT>
void computeInternalForceLinearElasticCoupled
(
	const double* xOverlap,
	const ScalarT* yOverlap,
	const double* mOwned,
	const double* volumeOverlap,
	const ScalarT* dilatationOwned,
	const double* bondDamage,
	const double* dsfOwned,
	ScalarT* fInternalOverlap,
	const int*  localNeighborList,
	int numOwnedPoints,
	double BULK_MODULUS,
	double SHEAR_MODULUS,
	double horizon,
	double thermalExpansionCoefficient,
	ScalarT* deltaTemperature
)
{

// 	Compute processor local contribution to internal force

	double K = BULK_MODULUS;
	double MU = SHEAR_MODULUS;
	const double *xOwned = xOverlap;
	const ScalarT *yOwned = yOverlap;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	const double *dsf = dsfOwned;
	const ScalarT *theta = dilatationOwned;
	ScalarT *fOwned = fInternalOverlap;

	const int *neighPtr = localNeighborList;
	double cellVolume, alpha, X_dx, X_dy, X_dz, zeta, omega;
	ScalarT Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz, e, c1;
	for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, deltaTemperature++, m++, theta++, dsf++){

		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const ScalarT *Y = yOwned;
		alpha = 15.0*MU/(*m);
		alpha *= (*dsf);
		double selfCellVolume = v[p];
		ScalarT* deltaT = deltaTemperature;
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const ScalarT *YP = &yOverlap[3*localId];
			X_dx = XP[0]-X[0];
			X_dy = XP[1]-X[1];
			X_dz = XP[2]-X[2];
			zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
			Y_dx = YP[0]-Y[0];
			Y_dy = YP[1]-Y[1];
			Y_dz = YP[2]-Y[2];
			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
			e = dY - zeta - thermalExpansionCoefficient*(*deltaT)*zeta;

			omega = scalarInfluenceFunction(zeta,horizon);
			// c1 = omega*(*theta)*(9.0*K-15.0*MU)/(3.0*(*m));
			// NOTE: set pressure effect to maximum, this is not to be a permanent change
			c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0);
			t = (1.0-*bondDamage)*(c1 * zeta + (1.0-*bondDamage) * omega * alpha * e);
			fx = t * Y_dx / dY;
			fy = t * Y_dy / dY;
			fz = t * Y_dz / dY;

			*(fOwned+0) += fx*cellVolume;
			*(fOwned+1) += fy*cellVolume;
			*(fOwned+2) += fz*cellVolume;
			fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
			fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
			fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
		}
	}
}

// Explicit template instantiation for double.
template void computeInternalForceLinearElasticCoupled<double>
(
	const double* xOverlap,
	const double* yOverlap,
// 	const double* temperatureYOverlap,
	const double* mOwned,
	const double* volumeOverlap,
	const double* dilatationOwned,
	const double* bondDamage,
	const double* dsfOwned,
	double* fInternalOverlap,
	const int*  localNeighborList,
	int numOwnedPoints,
	double BULK_MODULUS,
	double SHEAR_MODULUS,
	double horizon,
	double thermalExpansionCoefficient,
	double* deltaTemperature
);

// Explicit template instantiation for Sacado::Fad::DFad<double>.
template void computeInternalForceLinearElasticCoupled<Sacado::Fad::DFad<double> >
(
	const double* xOverlap,
	const Sacado::Fad::DFad<double>* yOverlap,
// 	const Sacado::Fad::DFad<double>* temperatureYOverlap,
	const double* mOwned,
	const double* volumeOverlap,
	const Sacado::Fad::DFad<double>* dilatationOwned,
	const double* bondDamage,
	const double* dsfOwned,
	Sacado::Fad::DFad<double>* fInternalOverlap,
	const int*  localNeighborList,
	int numOwnedPoints,
	double BULK_MODULUS,
	double SHEAR_MODULUS,
	double horizon,
	double thermalExpansionCoefficient,
	Sacado::Fad::DFad<double>* deltaTemperature
);
}