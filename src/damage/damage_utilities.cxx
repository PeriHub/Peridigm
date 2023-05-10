//! \file damage_utilities.cxx

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
// Christian Willberg    christian.willberg@dlr.de
//
// ************************************************************************
//@HEADER

#include "damage_utilities.h"
#include <cmath>
#include <vector>


namespace DAMAGE_UTILITIES {

void calculateDamageIndex
(
    const int numOwnedPoints,
    const int* ownedIDs,
    const double* vol,
    const int* neighborhoodList,
    const double* bondDamageNP1,
    double* damage
){

/*  void calculateDamageIndex()
    Function calculates the damage index. 
    The damage index is defined as damaged volume in relation the neighborhood volume.
    damageIndex = sum_i (brokenBonds_i * volume_i) / volumeNeighborhood
*/
    int neighborhoodListIndex = 0;
    int bondIndex = 0;
    double volume;
    int numNeighbors;
    int nodeId, neighborId;
    double totalDamage;
    for (int iID = 0; iID < numOwnedPoints; ++iID) {
        nodeId = ownedIDs[iID];
        numNeighbors = neighborhoodList[neighborhoodListIndex++];

        totalDamage = 0.0;
        volume = 0.0;
        for (int iNID = 0; iNID < numNeighbors; ++iNID) {
            
            neighborId = neighborhoodList[neighborhoodListIndex++];
            // must be zero to avoid synchronization errors
            totalDamage += bondDamageNP1[bondIndex++]*vol[neighborId];
            volume += vol[neighborId];
        }

        damage[nodeId] = totalDamage/volume;

    }
}


}
