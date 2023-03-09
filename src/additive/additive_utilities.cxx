//! \file additive_utilities.cxx

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

#include "additive_utilities.h"
#include <cmath>
#include <vector>


namespace ADDITIVE_UTILITIES {

void deleteAllBonds
(
    const int numOwnedPoints,
    const int* ownedIDs,
    const int* neighborhoodList,
    //double* bondDamageNP1,
    double* detachedNodes
){

/*  void calculateDamageIndex()
    Function delete all bonds in the init step. 
*/
    // int neighborhoodListIndex = 0;
    // int bondIndex = 0;
    // int numNeighbors;
    int nodeId;

    for (int iID = 0; iID < numOwnedPoints; ++iID) {
        nodeId = ownedIDs[iID];
        // numNeighbors = neighborhoodList[neighborhoodListIndex++];
        detachedNodes[nodeId] = 1;
        //for (int iNID = 0; iNID < numNeighbors; ++iNID) {       
        //    neighborhoodListIndex++;
        //    bondDamageNP1[bondIndex++] = 1;
        //}
    }


    }
}
