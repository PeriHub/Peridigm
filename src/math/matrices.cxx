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

void setToZero(
    double* A,
    int len
)
{
    for (int i=0 ; i<len ; ++i){
        A[i] = 0.0;
    } 
}

template<typename ScalarT>
int Invert2by2Matrix
(
 const ScalarT* matrix,
 ScalarT& determinant,
 ScalarT* inverse
)
{
  int returnCode(0);
  ScalarT a =  *(matrix);
  ScalarT b =  *(matrix+1);
  ScalarT c =  *(matrix+3);
  ScalarT d =  *(matrix+4);
  
  
  determinant = a*d - b*c ;

  if(determinant == ScalarT(0.0)){
      returnCode = 1;
      *(inverse)   = 0.0;
      *(inverse+1) = 0.0;
      *(inverse+2) = 0.0;
      *(inverse+3) = 0.0;
      *(inverse+4) = 0.0;
      *(inverse+5) = 0.0;
      *(inverse+6) = 0.0;
      *(inverse+7) = 0.0;
      *(inverse+8) = 0.0;
  }
  else{
      *(inverse)   = d/determinant;
      *(inverse+1) = -1.0 * b/determinant;
      *(inverse+2) = 0.0;
      *(inverse+3) = -1.0 * c/determinant;
      *(inverse+4) = a/determinant;
      *(inverse+5) = 0.0;
      *(inverse+6) = 0.0;
      *(inverse+7) = 0.0;
      *(inverse+8) = 0.0;
  }
    
  return returnCode;
}

template<typename ScalarT>
int Invert3by3Matrix
(
 const ScalarT* matrix,
 ScalarT& determinant,
 ScalarT* inverse
)
{
  int returnCode(0);

  ScalarT minor0 =  *(matrix+4) * *(matrix+8) - *(matrix+5) * *(matrix+7);
  ScalarT minor1 =  *(matrix+3) * *(matrix+8) - *(matrix+5) * *(matrix+6);
  ScalarT minor2 =  *(matrix+3) * *(matrix+7) - *(matrix+4) * *(matrix+6);
  ScalarT minor3 =  *(matrix+1) * *(matrix+8) - *(matrix+2) * *(matrix+7);
  ScalarT minor4 =  *(matrix)   * *(matrix+8) - *(matrix+6) * *(matrix+2);
  ScalarT minor5 =  *(matrix)   * *(matrix+7) - *(matrix+1) * *(matrix+6);
  ScalarT minor6 =  *(matrix+1) * *(matrix+5) - *(matrix+2) * *(matrix+4);
  ScalarT minor7 =  *(matrix)   * *(matrix+5) - *(matrix+2) * *(matrix+3);
  ScalarT minor8 =  *(matrix)   * *(matrix+4) - *(matrix+1) * *(matrix+3);
  determinant = *(matrix) * minor0 - *(matrix+1) * minor1 + *(matrix+2) * minor2;

  if(determinant == ScalarT(0.0)){
    returnCode = 1;
    *(inverse) = 0.0;
    *(inverse+1) = 0.0;
    *(inverse+2) = 0.0;
    *(inverse+3) = 0.0;
    *(inverse+4) = 0.0;
    *(inverse+5) = 0.0;
    *(inverse+6) = 0.0;
    *(inverse+7) = 0.0;
    *(inverse+8) = 0.0;
  }
  else{
    *(inverse) = minor0/determinant;
    *(inverse+1) = -1.0*minor3/determinant;
    *(inverse+2) = minor6/determinant;
    *(inverse+3) = -1.0*minor1/determinant;
    *(inverse+4) = minor4/determinant;
    *(inverse+5) = -1.0*minor7/determinant;
    *(inverse+6) = minor2/determinant;
    *(inverse+7) = -1.0*minor5/determinant;
    *(inverse+8) = minor8/determinant;
  }

  return returnCode;
}

template<typename ScalarT>
void MatrixMultiply
(
 bool transA,
 bool transB,
 ScalarT alpha,
 const ScalarT* a,
 const ScalarT* b,
 ScalarT* result
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

  if(!transA && !transB){
    *(result+0) = *(a+0) * *(b+0) + *(a+1) * *(b+3) + *(a+2) * *(b+6);
    *(result+1) = *(a+0) * *(b+1) + *(a+1) * *(b+4) + *(a+2) * *(b+7);
    *(result+2) = *(a+0) * *(b+2) + *(a+1) * *(b+5) + *(a+2) * *(b+8);
    *(result+3) = *(a+3) * *(b+0) + *(a+4) * *(b+3) + *(a+5) * *(b+6);
    *(result+4) = *(a+3) * *(b+1) + *(a+4) * *(b+4) + *(a+5) * *(b+7);
    *(result+5) = *(a+3) * *(b+2) + *(a+4) * *(b+5) + *(a+5) * *(b+8);
    *(result+6) = *(a+6) * *(b+0) + *(a+7) * *(b+3) + *(a+8) * *(b+6);
    *(result+7) = *(a+6) * *(b+1) + *(a+7) * *(b+4) + *(a+8) * *(b+7);
    *(result+8) = *(a+6) * *(b+2) + *(a+7) * *(b+5) + *(a+8) * *(b+8);
  }
  else if(transA && !transB){
    *(result+0) = *(a+0) * *(b+0) + *(a+3) * *(b+3) + *(a+6) * *(b+6);
    *(result+1) = *(a+0) * *(b+1) + *(a+3) * *(b+4) + *(a+6) * *(b+7);
    *(result+2) = *(a+0) * *(b+2) + *(a+3) * *(b+5) + *(a+6) * *(b+8);
    *(result+3) = *(a+1) * *(b+0) + *(a+4) * *(b+3) + *(a+7) * *(b+6);
    *(result+4) = *(a+1) * *(b+1) + *(a+4) * *(b+4) + *(a+7) * *(b+7);
    *(result+5) = *(a+1) * *(b+2) + *(a+4) * *(b+5) + *(a+7) * *(b+8);
    *(result+6) = *(a+2) * *(b+0) + *(a+5) * *(b+3) + *(a+8) * *(b+6);
    *(result+7) = *(a+2) * *(b+1) + *(a+5) * *(b+4) + *(a+8) * *(b+7);
    *(result+8) = *(a+2) * *(b+2) + *(a+5) * *(b+5) + *(a+8) * *(b+8);
  }
  else if(!transA && transB){
    *(result+0) = *(a+0) * *(b+0) + *(a+1) * *(b+1) + *(a+2) * *(b+2);
    *(result+1) = *(a+0) * *(b+3) + *(a+1) * *(b+4) + *(a+2) * *(b+5);
    *(result+2) = *(a+0) * *(b+6) + *(a+1) * *(b+7) + *(a+2) * *(b+8);
    *(result+3) = *(a+3) * *(b+0) + *(a+4) * *(b+1) + *(a+5) * *(b+2);
    *(result+4) = *(a+3) * *(b+3) + *(a+4) * *(b+4) + *(a+5) * *(b+5);
    *(result+5) = *(a+3) * *(b+6) + *(a+4) * *(b+7) + *(a+5) * *(b+8);
    *(result+6) = *(a+6) * *(b+0) + *(a+7) * *(b+1) + *(a+8) * *(b+2);
    *(result+7) = *(a+6) * *(b+3) + *(a+7) * *(b+4) + *(a+8) * *(b+5);
    *(result+8) = *(a+6) * *(b+6) + *(a+7) * *(b+7) + *(a+8) * *(b+8);
  }
  else{
    *(result+0) = *(a+0) * *(b+0) + *(a+3) * *(b+1) + *(a+6) * *(b+2);
    *(result+1) = *(a+0) * *(b+3) + *(a+3) * *(b+4) + *(a+6) * *(b+5);
    *(result+2) = *(a+0) * *(b+6) + *(a+3) * *(b+7) + *(a+6) * *(b+8);
    *(result+3) = *(a+1) * *(b+0) + *(a+4) * *(b+1) + *(a+7) * *(b+2);
    *(result+4) = *(a+1) * *(b+3) + *(a+4) * *(b+4) + *(a+7) * *(b+5);
    *(result+5) = *(a+1) * *(b+6) + *(a+4) * *(b+7) + *(a+7) * *(b+8);
    *(result+6) = *(a+2) * *(b+0) + *(a+5) * *(b+1) + *(a+8) * *(b+2);
    *(result+7) = *(a+2) * *(b+3) + *(a+5) * *(b+4) + *(a+8) * *(b+5);
    *(result+8) = *(a+2) * *(b+6) + *(a+5) * *(b+7) + *(a+8) * *(b+8);
  }

  if(alpha != 1.0){
    for(int i=0 ; i<9 ; ++i)
      *(result+i) *= alpha;
  }
}

template<typename ScalarT>
void MatrixMultiply3x3
(
 const ScalarT A[][3],
 const ScalarT B[][3],
 ScalarT C[][3]
)
{
    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
    C[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
    C[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];
    C[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
    C[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
    C[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
}

template<typename ScalarT>
void MatrixMultiply3x3fromVector
(
 const ScalarT  A[][3],
 const ScalarT* B,
 ScalarT C[][3]
)
{
    C[0][0] = A[0][0] * *(B)   + A[0][1] * *(B+3) + A[0][2] * *(B+6);
    C[0][1] = A[0][0] * *(B+1) + A[0][1] * *(B+4) + A[0][2] * *(B+7);
    C[0][2] = A[0][0] * *(B+2) + A[0][1] * *(B+5) + A[0][2] * *(B+8);
    C[1][0] = A[1][0] * *(B)   + A[1][1] * *(B+3) + A[1][2] * *(B+6);
    C[1][1] = A[1][0] * *(B+1) + A[1][1] * *(B+4) + A[1][2] * *(B+7);
    C[1][2] = A[1][0] * *(B+2) + A[1][1] * *(B+5) + A[1][2] * *(B+8);
    C[2][0] = A[2][0] * *(B)   + A[2][1] * *(B+3) + A[2][2] * *(B+6);
    C[2][1] = A[2][0] * *(B+1) + A[2][1] * *(B+4) + A[2][2] * *(B+7);
    C[2][2] = A[2][0] * *(B+2) + A[2][1] * *(B+5) + A[2][2] * *(B+8);
}

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

template<typename ScalarT>
void setOnesOnDiagonalFullTensor(ScalarT* tensor, int numPoints){
 
  ScalarT *tens = tensor;

  for(int iID=0; iID<numPoints; ++iID, tens+=9){
      *(tens) = 1.0;
      *(tens+4) = 1.0;
      *(tens+8) = 1.0;
  }  
}
template<typename ScalarT>
void tensorRotation
(
    const double* angles,
    const ScalarT* tensorIn,
    const bool globToLoc,
    ScalarT* tensorOut
)
{   
    std::vector<ScalarT> rotMatVec(9);
    ScalarT* rotMat = &rotMatVec[0];
    std::vector<ScalarT> tempVec(9);
    ScalarT* temp = &tempVec[0];
    ScalarT A = 1;
    MATRICES::createRotationMatrix(angles,rotMat);
    //MATRICES::TransposeMatrix(rotMat,rotMatT);
    // geomNL
    if (globToLoc){
        MATRICES::MatrixMultiply(true,  false, A, rotMat,tensorIn, temp);
        MATRICES::MatrixMultiply(false, false, A, temp,rotMat,tensorOut);
    }
    else{
        MATRICES::MatrixMultiply(false, false, A, rotMat,tensorIn, temp);
        MATRICES::MatrixMultiply(false, true,  A, temp,rotMat,tensorOut);
    }
}

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template<typename ScalarT>
void createRotationMatrix
(
const double* alpha,
ScalarT* rotMat
){
    const double PI  = PeridigmNS::value_of_pi();
    double rad[3];
    std::vector<ScalarT> rotMatXVec(9);
    ScalarT* rotMatX = &rotMatXVec[0];
    std::vector<ScalarT> rotMatYVec(9);
    ScalarT* rotMatY = &rotMatYVec[0];
    std::vector<ScalarT> rotMatZVec(9);
    ScalarT* rotMatZ = &rotMatZVec[0];    
    std::vector<ScalarT> tempVec(9);
    ScalarT* temp = &tempVec[0];
    ScalarT A = 1.0;

    rad[0] = alpha[0]*(PI)/180.;
    rad[1] = alpha[1]*(PI)/180.;
    rad[2] = alpha[2]*(PI)/180.;
 

    // x - direction
    *(rotMatX) = 1;   *(rotMatX+1) = 0;           *(rotMatX+2) = 0;
    *(rotMatX+3) = 0; *(rotMatX+4) = cos(rad[0]); *(rotMatX+5) = -sin(rad[0]);
    *(rotMatX+6) = 0; *(rotMatX+7) = sin(rad[0]); *(rotMatX+8) =  cos(rad[0]);
    // y - direction
    *(rotMatY)   =  cos(rad[1]); *(rotMatY+1) = 0; *(rotMatY+2) = sin(rad[1]);
    *(rotMatY+3) = 0;            *(rotMatY+4) = 1; *(rotMatY+5) = 0;
    *(rotMatY+6) = -sin(rad[1]); *(rotMatY+7) = 0; *(rotMatY+8) = cos(rad[1]);
    // z - direction
    *(rotMatZ) = cos(rad[2]);   *(rotMatZ+1) = -sin(rad[2]); *(rotMatZ+2) = 0;
    *(rotMatZ+3) = sin(rad[2]); *(rotMatZ+4) =  cos(rad[2]); *(rotMatZ+5) = 0;
    *(rotMatZ+6) = 0;           *(rotMatZ+7) = 0;            *(rotMatZ+8) = 1;
    
            
    MATRICES::MatrixMultiply(false, false, A, rotMatX, rotMatY, temp);
    MATRICES::MatrixMultiply(false, false, A, temp, rotMatZ, rotMat);
}


//// Explicit template instantiation for double
template void createRotationMatrix<double>
(
const double* alpha,
double* rotMat
);
template void createRotationMatrix<Sacado::Fad::DFad<double> >
(
const double* alpha,
Sacado::Fad::DFad<double>* rotMat
);



template<typename ScalarT>
void TransposeMatrix
(
 const ScalarT* matrix,
 ScalarT* transpose
)
{
  // Store some values so that the matrix and transpose can be the
  // same matrix (i.e., transpose in place)
  ScalarT temp_xy( *(matrix+1) );
  ScalarT temp_xz( *(matrix+2) );
  ScalarT temp_yz( *(matrix+5) );

  *(transpose)   = *(matrix);
  *(transpose+1) = *(matrix+3);
  *(transpose+2) = *(matrix+6);
  *(transpose+3) = temp_xy;
  *(transpose+4) = *(matrix+4);
  *(transpose+5) = *(matrix+7);
  *(transpose+6) = temp_xz;
  *(transpose+7) = temp_yz;
  *(transpose+8) = *(matrix+8);
}


template int Invert2by2Matrix<double>
(
 const double* matrix,
 double& determinant,
 double* inverse
);

template int Invert2by2Matrix<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* matrix,
 Sacado::Fad::DFad<double>& determinant,
 Sacado::Fad::DFad<double>* inverse
);

template int Invert3by3Matrix<double>
(
 const double* matrix,
 double& determinant,
 double* inverse
);



template int Invert3by3Matrix<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* matrix,
 Sacado::Fad::DFad<double>& determinant,
 Sacado::Fad::DFad<double>* inverse
);

template void MatrixMultiply<double>
(
 bool transA,
 bool transB,
 double alpha,
 const double* a,
 const double* b,
 double* result
);

template void MatrixMultiply<Sacado::Fad::DFad<double> >
(
 bool transA,
 bool transB,
 Sacado::Fad::DFad<double> alpha,
 const Sacado::Fad::DFad<double>* a,
 const Sacado::Fad::DFad<double>* b,
 Sacado::Fad::DFad<double>* result
);

template void MatrixMultiply3x3<double>
(
 const double A[][3],
 const double B[][3],
 double C[][3]
);

template void MatrixMultiply3x3<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double> A[][3],
 const Sacado::Fad::DFad<double> B[][3],
 Sacado::Fad::DFad<double> C[][3]
);




template void MatrixMultiply3x3fromVector<double>
(
 const double  A[][3],
 const double* B,
 double C[][3]
);

template void MatrixMultiply3x3fromVector<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>  A[][3],
 const Sacado::Fad::DFad<double>* B,
 Sacado::Fad::DFad<double> C[][3]
);

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

template void setOnesOnDiagonalFullTensor<double>
(
 double* tensor,
 int numPoints
);

template void TransposeMatrix<double>
(
 const double* matrix,
 double* transpose
);

template void TransposeMatrix<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* matrix,
 Sacado::Fad::DFad<double>* transpose
);


#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double PYTHAG(double a, 
                     double b)
{
  double at = fabs(a), bt = fabs(b), ct, result;

  if (at > bt){ 
    ct = bt / at; 
    result = at * sqrt(1.0 + ct * ct); 
  } else if(bt > 0.0){
    ct = at / bt; 
    result = bt * sqrt(1.0 + ct * ct); 
  } else 
    result = 0.0;

  return(result);
}

int dsvd(double *a,
         int m,
         int n,
         double *w,
         double *v)
{
  int flag, i, its, j, jj, k;
  int l = 0;
  int nm = 0;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;
  double *rv1;

  if (m < n){
    fprintf(stderr, "#rows must be > #cols \n");
    return(0);
  }

  rv1 = (double *)malloc((unsigned int) n*sizeof(double));

  /* Householder reduction to bidiagonal form */
  for (i = 0; i < n; i++){
    /* left-hand reduction */
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;

    if (i < m){
      for (k = i; k < m; k++)
        scale += fabs((double)a[k*n+i]);

      if (scale){
        for (k = i; k < m; k++){
          a[k*n+i] = (double)((double)a[k*n+i]/scale);
          s += ((double)a[k*n+i] * (double)a[k*n+i]);
        }

        f = (double)a[i*n+i];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        a[i*n+i] = (double)(f - g);

        if (i != n - 1){
          for (j = l; j < n; j++){
            for (s = 0.0, k = i; k < m; k++)
              s += ((double)a[k*n+i] * (double)a[k*n+j]);

            f = s / h;

            for (k = i; k < m; k++)
              a[k*n+j] += (double)(f * (double)a[k*n+i]);
          }
        }

        for (k = i; k < m; k++)
          a[k*n+i] = (double)((double)a[k*n+i]*scale);
      }
    }
    w[i] = (double)(scale * g);

    /* right-hand reduction */
    g = s = scale = 0.0;

    if (i < m && i != n - 1){
      for (k = l; k < n; k++)
        scale += fabs((double)a[i*n+k]);

      if (scale){
        for (k = l; k < n; k++){
          a[i*n+k] = (double)((double)a[i*n+k]/scale);
          s += ((double)a[i*n+k] * (double)a[i*n+k]);
        }

        f = (double)a[i*n+l];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        a[i*n+l] = (double)(f - g);

        for (k = l; k < n; k++)
          rv1[k] = (double)a[i*n+k] / h;

        if (i != m - 1){
          for (j = l; j < m; j++){
            for (s = 0.0, k = l; k < n; k++)
              s += ((double)a[j*n+k] * (double)a[i*n+k]);

            for (k = l; k < n; k++)
              a[j*n+k] += (double)(s * rv1[k]);
          }
        }

        for (k = l; k < n; k++)
        a[i*n+k] = (double)((double)a[i*n+k]*scale);
      }
    }

    anorm = fmax(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
  }

  /* accumulate the right-hand transformation */
  for (i = n - 1; i >= 0; i--){
    if (i < n - 1){
      if (g){
        for (j = l; j < n; j++)
          v[j*n+i] = (double)(((double)a[i*n+j] / (double)a[i*n+l]) / g);

        /* double division to avoid underflow */
        for (j = l; j < n; j++){
          for (s = 0.0, k = l; k < n; k++)
            s += ((double)a[i*n+k] * (double)v[k*n+j]);

          for (k = l; k < n; k++)
            v[k*n+j] += (double)(s * (double)v[k*n+i]);
        }
      }
      for (j = l; j < n; j++)
        v[i*n+j] = v[j*n+i] = 0.0;
    }
    v[i*n+i] = 1.0;
    g = rv1[i];
    l = i;
  }

  /* accumulate the left-hand transformation */
  for (i = n - 1; i >= 0; i--){
    l = i + 1;
    g = (double)w[i];

    if (i < n - 1)
      for (j = l; j < n; j++)
        a[i*n+j] = 0.0;

    if (g){
      g = 1.0 / g;

      if (i != n - 1){
        for (j = l; j < n; j++){
          for (s = 0.0, k = l; k < m; k++)
            s += ((double)a[k*n+i] * (double)a[k*n+j]);

          f = (s / (double)a[i*n+i]) * g;

          for (k = i; k < m; k++)
            a[k*n+j] += (double)(f * (double)a[k*n+i]);
        }
      }

      for (j = i; j < m; j++)
      a[j*n+i] = (double)((double)a[j*n+i]*g);
    }
    else
    {
      for (j = i; j < m; j++)
        a[j*n+i] = 0.0;
    }
    ++a[i*n+i];
  }

  /* diagonalize the bidiagonal form */
  for (k = n - 1; k >= 0; k--){                             /* loop over singular values */
    for (its = 0; its < 30; its++){                         /* loop over allowed iterations */
      flag = 1;

      for (l = k; l >= 0; l--){                     /* test for splitting */
        nm = l - 1;

        if (fabs(rv1[l]) + anorm == anorm){
          flag = 0;
          break;
        }

        if (fabs((double)w[nm]) + anorm == anorm)
          break;
      }

      if (flag){
        s = 1.0;

        for (i = l; i <= k; i++){
          f = s * rv1[i];

          if (fabs(f) + anorm != anorm){
            g = (double)w[i];
            h = PYTHAG(f, g);
            w[i] = (double)h;
            h = 1.0 / h;
            c = g * h;
            s = (- f * h);
            for (j = 0; j < m; j++){
              y = (double)a[j*n+nm];
              z = (double)a[j*n+i];
              a[j*n+nm] = (double)(y * c + z * s);
              a[j*n+i] = (double)(z * c - y * s);
            }
          }
        }
      }

      z = (double)w[k];
      if (l == k){                  /* convergence */
        if (z < 0.0){              /* make singular value nonnegative */
          w[k] = (double)(-z);
          for (j = 0; j < n; j++)
            v[j*n+k] = (-v[j*n+k]);
        }
        break;
      }

      if (its >= 30) {
        free (rv1);

        fprintf(stderr, "No convergence after 30,000! iterations \n");
        return(0);
      }

      /* shift from bottom 2 x 2 minor */
      x = (double)w[l];
      nm = k - 1;
      y = (double)w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = PYTHAG(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

      /* next QR transformation */
      c = s = 1.0;

      for (j = l; j <= nm; j++){
        i = j + 1;
        g = rv1[i];
        y = (double)w[i];
        h = s * g;
        g = c * g;
        z = PYTHAG(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < n; jj++){
          x = (double)v[jj*n+j];
          z = (double)v[jj*n+i];
          v[jj*n+j] = (double)(x * c + z * s);
          v[jj*n+i] = (double)(z * c - x * s);
        }

        z = PYTHAG(f, h);
        w[j] = (double)z;

        if (z){
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }

        f = (c * g) + (s * y);
        x = (c * y) - (s * g);

        for (jj = 0; jj < m; jj++){
          y = (double)a[jj*n+j];
          z = (double)a[jj*n+i];
          a[jj*n+j] = (double)(y * c + z * s);
          a[jj*n+i] = (double)(z * c - y * s);
        }
      }

      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = (double)x;
    }
  }

  free (rv1);

  return(1);
}


//! Invert a single N-by-N symmetric matrix; returns zero of successful, one if not successful (e.g., singular matrix).
template<typename ScalarT>
int invertAndCond(const ScalarT* Min,
                  ScalarT *Mout,
                  const int size,
                  const double thresVal)
{
  
  double *v = new double[size*size];
  double *w = new double[size];
  double *u = new double[size*size];  
  

  int i,j,k;

  for(i=0;i<size*size;i++){
    u[i]=Min[i];
    Mout[i]=0.0;
  }
  dsvd(u, size, size, w, v);

  for(i=0 ; i < size ; i++ ){
    for(j=0 ; j < size ; j++ ){
      for( k=0 ; k < size ; k++ ){
        if(w[k] > thresVal){ // pseudo inverse approach (ignore zero eigs)
          Mout[i*size+j] += v[i*size+k]*1.0/w[k]*u[j*size+k];
        }
      }
    }
  }

  //free(u);
  //free(w);
  //free(v);

  return 0;
}

double distance(    
    const double a1, 
    const double a2, 
    const double a3,
    const double b1, 
    const double b2, 
    const double b3) 
{
    return ( sqrt( (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) + (a3-b3)*(a3-b3) ) );
}

double vectorNorm
(
    const double* vector,
    const int len
){
    double accum = 0.;
    for (int i = 0; i < len; ++i) {
        accum += vector[i] * vector[i];
    }
    return sqrt(accum);
}
    double dist(
    const double a1, 
    const double a2, 
    const double a3,
    const double b1, 
    const double b2, 
    const double b3);

template int invertAndCond<double>
(
    const double* Min,
    double *Mout,
    const int size,
    const double thresVal
);
template void tensorRotation<double>
(
  const double* angles,
  const double* tensorIn,
  const bool globToLoc,
  double* tensorOut
);

template void tensorRotation<Sacado::Fad::DFad<double> >
(
  const double* angles,
  const Sacado::Fad::DFad<double>* tensorIn,
  const bool globToLoc,
  Sacado::Fad::DFad<double>* tensorOut
);

template<typename ScalarT>
int EigenVec2D
(
 const ScalarT* a,
 ScalarT* result
)
{
  // This function computes the Eigenvector of a 2D 3x3 matrix
  int returnCode(0);

  *(result+0) = -2 * *(a+1) / ( *(a+0) - *(a+4) + sqrt(pow( *(a+0), 2) - 2 * *(a+0) * *(a+4) + 4 * pow( *(a+1), 2) + pow( *(a+4), 2)));
  *(result+1) = -2 * *(a+1) / ( *(a+0) - *(a+4) - sqrt(pow( *(a+0), 2) - 2 * *(a+0) * *(a+4) + 4 * pow( *(a+1), 2) + pow( *(a+4), 2)));
  *(result+2) = 0.0;
  *(result+3) = 1;
  *(result+4) = 1;
  *(result+5) = 0.0;
  *(result+6) = 0.0;
  *(result+7) = 0.0;
  *(result+8) = 0.0;

  if(*(result+0)!=*(result+0) || *(result+1)!=*(result+1))
  {
    returnCode=1;
  }

  return returnCode;
}
template int EigenVec2D<double>(
    const double* a,
    double* result
);

template int EigenVec2D<Sacado::Fad::DFad<double>>(
    const  Sacado::Fad::DFad<double>* a,
     Sacado::Fad::DFad<double>* result
);
}