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
#include <Eigen/Core>

namespace MATERIAL_EVALUATION
{

  template <typename ScalarT>
  void computeInternalForceElasticBondBased(
      const double *xOverlap,
      const ScalarT *yOverlap,
      const double *volumeOverlap,
      const double *bondDamage,
      ScalarT *fInternalOverlap,
      const int *localNeighborList,
      int numOwnedPoints,
      double BULK_MODULUS,
      double horizon)
  {
    double volume, neighborVolume, X[3], neighborX[3], initialBondLength, damageOnBond;
    ScalarT Y[3], neighborY[3], currentBondLength, stretch, t, fx, fy, fz;
    int neighborhoodIndex(0), bondDamageIndex(0), neighborId;

    const double pi = PeridigmNS::value_of_pi();
    double constant = 18.0 * BULK_MODULUS / (pi * horizon * horizon * horizon * horizon);

    int bondId = 0;
    for (int p = 0; p < numOwnedPoints; p++)
    {

      X[0] = xOverlap[p * 3];
      X[1] = xOverlap[p * 3 + 1];
      X[2] = xOverlap[p * 3 + 2];
      Y[0] = yOverlap[p * 3];
      Y[1] = yOverlap[p * 3 + 1];
      Y[2] = yOverlap[p * 3 + 2];
      volume = volumeOverlap[p];

      int numNeighbors = localNeighborList[neighborhoodIndex++];
      bool print = true;
      for (int n = 0; n < numNeighbors; n++)
      {

        neighborId = localNeighborList[neighborhoodIndex++];
        neighborX[0] = xOverlap[neighborId * 3];
        neighborX[1] = xOverlap[neighborId * 3 + 1];
        neighborX[2] = xOverlap[neighborId * 3 + 2];
        neighborY[0] = yOverlap[neighborId * 3];
        neighborY[1] = yOverlap[neighborId * 3 + 1];
        neighborY[2] = yOverlap[neighborId * 3 + 2];
        neighborVolume = volumeOverlap[neighborId];

        initialBondLength = std::sqrt((neighborX[0] - X[0]) * (neighborX[0] - X[0]) + (neighborX[1] - X[1]) * (neighborX[1] - X[1]) + (neighborX[2] - X[2]) * (neighborX[2] - X[2]));
        currentBondLength = std::sqrt((neighborY[0] - Y[0]) * (neighborY[0] - Y[0]) + (neighborY[1] - Y[1]) * (neighborY[1] - Y[1]) + (neighborY[2] - Y[2]) * (neighborY[2] - Y[2]));
        stretch = (currentBondLength - initialBondLength) / initialBondLength;

        damageOnBond = bondDamage[bondDamageIndex++];

        t = 0.5 * (1.0 - damageOnBond) * stretch * constant;

        fx = t * (neighborY[0] - Y[0]) / currentBondLength;
        fy = t * (neighborY[1] - Y[1]) / currentBondLength;
        fz = t * (neighborY[2] - Y[2]) / currentBondLength;

        fInternalOverlap[3 * p + 0] += fx * neighborVolume;
        fInternalOverlap[3 * p + 1] += fy * neighborVolume;
        fInternalOverlap[3 * p + 2] += fz * neighborVolume;
        fInternalOverlap[3 * neighborId + 0] -= fx * volume;
        fInternalOverlap[3 * neighborId + 1] -= fy * volume;
        fInternalOverlap[3 * neighborId + 2] -= fz * volume;
        if (X[0]== 1.5 & X[1]==-9){
          bondId = p;
          std::cout << "X[0]: " << X[0] << std::endl;
          std::cout << "X[1]: " << X[1] << std::endl;
          std::cout << "neighborX[0]: " << neighborX[0] << std::endl;
          std::cout << "neighborX[1]: " << neighborX[1] << std::endl;
          std::cout << "Un1[DPN * n + 0]: " << (neighborY[0] - neighborX[0]) << std::endl;
          std::cout << "Un1[DPN * n + 1]: " << (neighborY[1] - neighborX[1]) << std::endl;
          std::cout << "fx: " << fx << std::endl;
          std::cout << "fy: " << fy << std::endl;
        }
      }
      
    }
    if(bondId!=0){
      std::cout << "fInternalOverlap[3 * bondId + 0]: " << fInternalOverlap[3 * bondId + 0] << std::endl;
      std::cout << "fInternalOverlap[3 * bondId + 1]: " << fInternalOverlap[3 * bondId + 1] << std::endl;
    }
  }

  template <typename ScalarT>
  void computeInternalForceElasticBondBasedCollocation(
      const double *xOverlap,
      const ScalarT *yOverlap,
      const double *volumeOverlap,
      const double *bondDamage,
      ScalarT *fInternalOverlap,
      const int *localNeighborList,
      const int *localCollocationNeighborList,
      bool *useCollocationNodes,
      const Eigen::MatrixXd *UpdateMat,
      int numOwnedPoints,
      double BULK_MODULUS,
      double horizon,
      const double m_criticalStretch)
  {
    double volume, neighborVolume, X[3], neighborX[3], initialBondLength, damageOnBond;
    ScalarT Y[3], neighborY[3], currentBondLength, stretch, t, fx, fy, fz;
    int neighborhoodIndex(0), collocationNeighborhoodIndex(0), bondDamageIndex(0), neighborId;

    const double pi = PeridigmNS::value_of_pi();
    double constant = 18.0 * BULK_MODULUS / (pi * horizon * horizon * horizon * horizon);
    int bondId = 0;
    for (int p = 0; p < numOwnedPoints; p++)
    {

      X[0] = xOverlap[p * 3];
      X[1] = xOverlap[p * 3 + 1];
      X[2] = xOverlap[p * 3 + 2];
      Y[0] = yOverlap[p * 3];
      Y[1] = yOverlap[p * 3 + 1];
      Y[2] = yOverlap[p * 3 + 2];
      volume = volumeOverlap[p];

      int numNeighbors = localNeighborList[neighborhoodIndex++];
      int numCollocationNeighbors = localCollocationNeighborList[collocationNeighborhoodIndex++];

      if (useCollocationNodes[p] == true)
      {
        neighborhoodIndex += numNeighbors;
        for (int n = 0; n < numCollocationNeighbors; n++)
        {

          neighborId = localCollocationNeighborList[collocationNeighborhoodIndex++];

          if (neighborId == -1)
					{
						continue;
					}

          neighborX[0] = xOverlap[neighborId * 3];
          neighborX[1] = xOverlap[neighborId * 3 + 1];
          neighborX[2] = xOverlap[neighborId * 3 + 2];
          neighborY[0] = yOverlap[neighborId * 3];
          neighborY[1] = yOverlap[neighborId * 3 + 1];
          neighborY[2] = yOverlap[neighborId * 3 + 2];
          neighborVolume = volumeOverlap[neighborId];

          initialBondLength = std::sqrt((neighborX[0] - X[0]) * (neighborX[0] - X[0]) + (neighborX[1] - X[1]) * (neighborX[1] - X[1]) + (neighborX[2] - X[2]) * (neighborX[2] - X[2]));
          currentBondLength = std::sqrt((neighborY[0] - Y[0]) * (neighborY[0] - Y[0]) + (neighborY[1] - Y[1]) * (neighborY[1] - Y[1]) + (neighborY[2] - Y[2]) * (neighborY[2] - Y[2]));
          stretch = (currentBondLength - initialBondLength) / initialBondLength;

          if (stretch > 0.95 * m_criticalStretch)
          {
            useCollocationNodes[p] = false;
            useCollocationNodes[neighborId] = false;
          }
          
          // ScalarT ux = neighborY[0] - neighborX[0];
          // ScalarT uy = neighborY[1] - neighborX[1];
          ScalarT ux = neighborY[0] - neighborX[0];
          ScalarT uy = neighborY[1] - neighborX[1];
          // std::cout << "ux: " << ux << std::endl;
          // std::cout << "uy: " << uy << std::endl;
          fx = ux * UpdateMat[p].row(0)(n) + uy * UpdateMat[p].row(1)(n);
          fy = ux * UpdateMat[p].row(2)(n) + uy * UpdateMat[p].row(3)(n);

          // fz = neighborY[0] * UpdateMat[p].row(2)(j) + neighborY[1] * UpdateMat[p].row(3)(j);
          fz = 0.0;

          // damageOnBond = bondDamage[bondDamageIndex++];
          
          // for (int k = 0; k < HorizonLengths[i]; k++)
          // {
          //   damageOnBond = bondDamage[bondDamageIndex++];
          // }

          // t = 0.5 * (1.0 - damageOnBond) * stretch * constant;

          // fx = t * (neighborY[0] - Y[0]) / currentBondLength;
          // fy = t * (neighborY[1] - Y[1]) / currentBondLength;
          // fz = t * (neighborY[2] - Y[2]) / currentBondLength;

          fInternalOverlap[3 * p + 0] += fx;
          fInternalOverlap[3 * p + 1] += fy;
          fInternalOverlap[3 * p + 2] += fz;
          fInternalOverlap[3 * neighborId + 0] -= fx;
          fInternalOverlap[3 * neighborId + 1] -= fy;
          fInternalOverlap[3 * neighborId + 2] -= fz;

          if (X[0]== 1.5 & X[1]==-9){
            bondId = p;
            std::cout << "X[0]: " << X[0] << std::endl;
            std::cout << "X[1]: " << X[1] << std::endl;
            std::cout << "neighborX[0]: " << neighborX[0] << std::endl;
            std::cout << "neighborX[1]: " << neighborX[1] << std::endl;
            std::cout << "Un1[DPN * n + 0]: " << (neighborY[0] - neighborX[0]) << std::endl;
            std::cout << "Un1[DPN * n + 1]: " << (neighborY[1] - neighborX[1]) << std::endl;
            std::cout << "fx: " << fx << std::endl;
            std::cout << "fy: " << fy << std::endl;
          }
          // fInternalOverlap[3 * neighborId + 0] -= fx * volume * (1.0 - damageOnBond);
          // fInternalOverlap[3 * neighborId + 1] -= fy * volume * (1.0 - damageOnBond);
          // fInternalOverlap[3 * neighborId + 2] -= fz * volume * (1.0 - damageOnBond);
          // std::cout << " fInternalOverlap[3 * p + 0]: " <<  fInternalOverlap[3 * p + 0] << std::endl;
          // std::cout << " fInternalOverlap[3 * p + 1]: " <<  fInternalOverlap[3 * p + 1] << std::endl;
          // std::cout << " fInternalOverlap[3 * p + 2]: " <<  fInternalOverlap[3 * p + 2] << std::endl;
          // std::cout << " fInternalOverlap[3 * neighborId + 0]: " <<  fInternalOverlap[3 * neighborId + 0] << std::endl;
          // std::cout << " fInternalOverlap[3 * neighborId + 1]: " <<  fInternalOverlap[3 * neighborId + 1] << std::endl;
          // std::cout << " fInternalOverlap[3 * neighborId + 2]: " <<  fInternalOverlap[3 * neighborId + 2] << std::endl;
          // if (neighborId == p){
          //   std::cout << " p: " <<  p << std::endl;
          //   std::cout << " neighborId: " <<  neighborId << std::endl;
          //   std::cout << " fInternalOverlap[3 * p + 0]: " <<  fInternalOverlap[3 * p + 0] << std::endl;
          //   std::cout << " fx: " <<  fx << std::endl;
          //   std::cout << " neighborY[0]: " <<  neighborY[0] << std::endl;
          //   std::cout << " Y[0]: " <<  Y[0] << std::endl;
          //   std::cout << " currentBondLength: " <<  currentBondLength << std::endl;
          //   std::cout << " neighborVolume: " <<  neighborVolume << std::endl;
          // }
        }
      }
      else
      {
        collocationNeighborhoodIndex += numCollocationNeighbors;

        for (int n = 0; n < numNeighbors; n++)
        {
        //  TODO: IDS sind identisch !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          neighborId = localNeighborList[neighborhoodIndex++];
          neighborX[0] = xOverlap[neighborId * 3];
          neighborX[1] = xOverlap[neighborId * 3 + 1];
          neighborX[2] = xOverlap[neighborId * 3 + 2];
          neighborY[0] = yOverlap[neighborId * 3];
          neighborY[1] = yOverlap[neighborId * 3 + 1];
          neighborY[2] = yOverlap[neighborId * 3 + 2];
          neighborVolume = volumeOverlap[neighborId];

          initialBondLength = std::sqrt((neighborX[0] - X[0]) * (neighborX[0] - X[0]) + (neighborX[1] - X[1]) * (neighborX[1] - X[1]) + (neighborX[2] - X[2]) * (neighborX[2] - X[2]));
          currentBondLength = std::sqrt((neighborY[0] - Y[0]) * (neighborY[0] - Y[0]) + (neighborY[1] - Y[1]) * (neighborY[1] - Y[1]) + (neighborY[2] - Y[2]) * (neighborY[2] - Y[2]));
          stretch = (currentBondLength - initialBondLength) / initialBondLength;

          damageOnBond = bondDamage[bondDamageIndex++];

          t = 0.5 * (1.0 - damageOnBond) * stretch * constant;

          fx = t * (neighborY[0] - Y[0]) / currentBondLength;
          fy = t * (neighborY[1] - Y[1]) / currentBondLength;
          fz = t * (neighborY[2] - Y[2]) / currentBondLength;

          fInternalOverlap[3 * p + 0] += fx * neighborVolume;
          fInternalOverlap[3 * p + 1] += fy * neighborVolume;
          fInternalOverlap[3 * p + 2] += fz * neighborVolume;
          fInternalOverlap[3 * neighborId + 0] -= fx * volume;
          fInternalOverlap[3 * neighborId + 1] -= fy * volume;
          fInternalOverlap[3 * neighborId + 2] -= fz * volume;
        
        }
      }
    }
    if(bondId!=0){
      std::cout << "fInternalOverlap[3 * bondId + 0]: " << fInternalOverlap[3 * bondId + 0] << std::endl;
      std::cout << "fInternalOverlap[3 * bondId + 1]: " << fInternalOverlap[3 * bondId + 1] << std::endl;
    }
  }

  /** Explicit template instantiation for double. */
  template void computeInternalForceElasticBondBased<double>(
      const double *xOverlap,
      const double *yOverlap,
      const double *volumeOverlap,
      const double *bondDamage,
      double *fInternalOverlap,
      const int *localNeighborList,
      int numOwnedPoints,
      double BULK_MODULUS,
      double horizon);

  /** Explicit template instantiation for Sacado::Fad::DFad<double>. */
  template void computeInternalForceElasticBondBased<Sacado::Fad::DFad<double>>(
      const double *xOverlap,
      const Sacado::Fad::DFad<double> *yOverlap,
      const double *volumeOverlap,
      const double *bondDamage,
      Sacado::Fad::DFad<double> *fInternalOverlap,
      const int *localNeighborList,
      int numOwnedPoints,
      double BULK_MODULUS,
      double horizon);

  /** Explicit template instantiation for double. */
  template void computeInternalForceElasticBondBasedCollocation<double>(
      const double *xOverlap,
      const double *yOverlap,
      const double *volumeOverlap,
      const double *bondDamage,
      double *fInternalOverlap,
      const int *localNeighborList,
      const int *localCollocationNeighborList,
      bool *useCollocationNodes,
      const Eigen::MatrixXd *UpdateMat,
      int numOwnedPoints,
      double BULK_MODULUS,
      double horizon,
      const double m_criticalStretch);

  /** Explicit template instantiation for Sacado::Fad::DFad<double>. */
  template void computeInternalForceElasticBondBasedCollocation<Sacado::Fad::DFad<double>>(
      const double *xOverlap,
      const Sacado::Fad::DFad<double> *yOverlap,
      const double *volumeOverlap,
      const double *bondDamage,
      Sacado::Fad::DFad<double> *fInternalOverlap,
      const int *localNeighborList,
      const int *localCollocationNeighborList,
      bool *useCollocationNodes,
      const Eigen::MatrixXd *UpdateMat,
      int numOwnedPoints,
      double BULK_MODULUS,
      double horizon,
      const double m_criticalStretch);

}
