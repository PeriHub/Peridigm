//! \file matrices.h
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
// Author of this Routine
// Jan-Timo Hesse   Jan-Timo.Hesse@dlr.de
// German Aerospace Center
//@HEADER
#ifndef MATRICES_H
#define MATRICES_H

#include "Peridigm_Constants.hpp"

namespace MATRICES {

//set double fields to zero
template<typename ScalarT>
void setToZero(
    ScalarT* A,
    const int len
);

//! Invert a single 2-by-2 matrix; returns zero of successful, one if not successful (e.g., singular matrix).
template<typename ScalarT>
int Invert2by2Matrix
(
const ScalarT* matrix,
ScalarT& determinant,
ScalarT* inverse
);
// create a 3x3 rotation matrix with 3 input angles
template<typename ScalarT>
void createRotationMatrix
(
const double* alpha,
ScalarT* rotMat
);
//! Invert a single 3-by-3 matrix; returns zero of successful, one if not successful (e.g., singular matrix).
template<typename ScalarT>
int Invert3by3Matrix
(
const ScalarT* matrix,
ScalarT& determinant,
ScalarT* inverse
);

template<typename ScalarT>
int EigenVec2D
(
 const ScalarT* a,
 ScalarT* result
);
// rotates a second order tensor in a new configuration or back
template<typename ScalarT>
void tensorRotation
(
const double* angles,
const ScalarT* tensorIn,
const bool globToLoc,
ScalarT* tensorOut
);

double distance(
    const double a1, 
    const double a2, 
    const double a3,
    const double b1, 
    const double b2, 
    const double b3
    );

double vectorNorm
(
const double* vector,
const int len
);

//! Invert a single N-by-N symmetric matrix; returns zero of successful, one if not successful (e.g., singular matrix).

//! Inner product of two 3-by-3 matrices.
template<typename ScalarT>
void MatrixMultiply
(
bool transA,
bool transB,
ScalarT alpha,
const ScalarT* a,
const ScalarT* b,
ScalarT* result
);

template<typename ScalarT>
void MatrixMultiply3x3
(
const ScalarT A[][3],
const ScalarT B[][3],
ScalarT C[][3]
);

template<typename ScalarT>
void MatrixMultiply3x3fromVector
(
 const ScalarT  A[][3],
 const ScalarT* B,
 ScalarT C[][3]
);


template<typename ScalarT>
void MatMul
(
int n,
const ScalarT A[][6],
const ScalarT B[][6],
ScalarT C[][6],
bool transpose
);

template<typename ScalarT>
void setOnesOnDiagonalFullTensor(ScalarT* tensor, int numPoints);

template<typename ScalarT>
void setOnesOnDiagonalFullTensor(ScalarT* tensor, int numPoints);

//! Transpose matrix; if both arguments are the same pointer then the matrix is transposed in place.
template<typename ScalarT>
void TransposeMatrix
(
const ScalarT* matrix,
ScalarT* transpose
);
template<typename ScalarT>
int invertAndCond
(
const ScalarT* Min,
ScalarT *Mout,
const int size,
const double thresVal
);

template<typename ScalarT>
void DIFFTENSOR
(
const ScalarT* TENSORN,
const ScalarT* TENSORNP1, 
ScalarT* DTENSOR
);
}
#endif // MATRICES_H