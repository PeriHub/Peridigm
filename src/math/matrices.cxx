//! \file matrices.cxx

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

#include "matrices.h"
#include <Sacado.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <math.h>
#include <functional>
#include <cmath> 
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <vector> 
#include <string> 

namespace MATRICES {

template<typename ScalarT>
void MatMul
(
 int n,
 const ScalarT A[][6],
 const ScalarT B[][6],
 ScalarT C[][6],
 bool transpose
)
{
  // This function computes result = alpha * a * b
  // where alpha is a scalar and a and b are 3x3 matrices
  // The arguments transA and transB denote whether or not
  // to use the transpose of a and b, respectively.

  // The default ordering is row-major:
  //
  // XX(0) XY(1) XZ(2)
  // YX(3) YY(4) YZ(5)
  // ZX(6) ZY(7) ZZ(8)
  if (transpose==false){
    for (int iID = 0; iID < n; ++iID){
    /* For each column j of B */
        for(int jID=0 ; jID<n ; ++jID){
      /* Compute C(i,j) */
            C[iID][jID] = 0;
            for(int kID=0 ; kID<n ; ++kID){
                // transponiertes tm
                C[iID][jID] += A[iID][kID]*B[kID][jID];
                }
        }
    }
  }
  
  else {
    for (int iID = 0; iID < n; ++iID){
    /* For each column j of B */
        for(int jID=0 ; jID<n ; ++jID){
      /* Compute C(i,j) */
            C[iID][jID] = 0;
            for(int kID=0 ; kID<n ; ++kID){
                // transponiertes tm
                C[iID][jID] += A[kID][iID]*B[kID][jID];
                }
        }
    }
  }
}

template void MatMul<double>
(
 int n,
 const double A[][6],
 const double B[][6],
 double C[][6],
 bool transpose
);
template void MatMul<Sacado::Fad::DFad<double> >
(
 int n,
 const Sacado::Fad::DFad<double> A[][6],
 const Sacado::Fad::DFad<double> B[][6],
 Sacado::Fad::DFad<double> C[][6],
 bool transpose
);
}