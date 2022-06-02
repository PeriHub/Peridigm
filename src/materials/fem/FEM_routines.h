//! \file FEM_routines.h

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
// ************************************************************************
//
// funded by dfg project reference
// Licence agreement
//
// Christian Willberg    christian.willberg@dlr.de
//@HEADER

#ifndef FEMROUTINES_H
#define FEMROUTINES_H

namespace FEM {

template<typename ScalarT>
void updateElasticCauchyStressFEM
(
ScalarT* DeformationGradient, 
ScalarT* unrotatedCauchyStressN, 
ScalarT* unrotatedCauchyStressNP1, 
int numPoints, 
const ScalarT Cstiff[][6],
double* angles,
int type,
double dt,
bool incremental,
bool hencky
);
void shapeFunctionsLagrangeRecursive
(
double* N, 
const int order, 
const double* xi,
const double elCoor
);
void getTopology(
    const int numOwnedPoints,
    const int *neighborhoodList,
    const double *nodeType,
    int numElements,
    std::vector<int> topology);

void derivativeShapeFunctionsLagrangeRecursive
(
double* B, 
const double* N,
const int order,
const double* xi,
const double elCoor
);
void getLagrangeElementData
(
const int order, 
const double elCoor,
double* Nxi,
double* Bxi
);

int getNumberOfIntegrationPoints
(
const bool twoD, 
const int numIntDir[3]
);
void setElementMatrices
(
const bool twoD,
const int offset,
const int order[3],
const double* Nxi,
const double* Neta,
const double* Npsi,
const double* Bxi,
const double* Beta,
const double* Bpsi,
double* Bx,
double* By,
double* Bz
);
void nodalForce
(
const double* Bx,
const double* By,
const double* Bz,
const int intPointPtr,
const int nnode,
const double detJ,
const double* Jinv, 
const bool twoD,
const double* sigmaInt, 
double* nforces
);

void tensorRotation
(
const double* angles,
const double* tensorIn,
const bool globToLoc,
double* tensorOut
);

void defineLagrangianGridSpace
(
const int order,
double* xi
);
void getElementTopo
(
const bool twoD, 
const int order[3], 
int *topo
);
void shapeFunctionsLagrange
(
double* Nmatrix, 
const int order[3], 
const double elCoor[3]
);

void weightsAndIntegrationPoints
(
const int order, 
double* elCoor,
double* weights
);



void computeStrain
(
const double* Bx,
const double* By,
const double* Bz,
const double* u, 
const int intPointPtr,
const int nnode,
const double* Jinv,
const bool twoD,
double* strain
);

double getJacobian
(
const double* Bx,
const double* By,
const double* Bz,
const int nnode, 
const int intPointPtr,
const double* coor,
const double weight,
const bool twoD,
double* J, 
double* Jinv
);

void setWeights
(
const int numIntDir[3],
const bool twoD, 
const double* weightsx,
const double* weightsy,
const double* weightsz,
double* weights
);
void setNodalStresses(
const int nnode,
const int elementID,
int topoPtr,
const int* topology,
double* sigmaNP1
);

void setToZero(
double* A,
int len
);
template<typename ScalarT>
void tensorRotation(
const double* angles,
const ScalarT* tensorIn,
const bool globToLoc,
ScalarT* tensorOut
);
void setGlobalStresses(
const int elementID,
const double* sigmaInt,
double* sigmaNP1
);
void setGlobalForces
(
const int nnode,
const int elementID,
int topoPtr,
const int* topology,
const double* elNodalForces,
const double* volume,
double* force
);

void getDisplacements
(
const int numOwnedPoints,
const double* modelCoordinates,
const double* coordinatesNP1,
double* displacements
);
int getNnode
(
const int order[3], 
const bool twoD
);
}

#endif // FEMROUTINES_H
