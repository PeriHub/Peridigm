/*! \file ut_umats.cpp */

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

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "user_material_interface_correspondence.h"
#include "user_material.h"
using namespace std;
using namespace Teuchos;


TEUCHOS_UNIT_TEST(correspondence, DIFFTENSOR) {
    const double tolerance = 1.0e-15;
    std::vector<double> AVector(9);
    double* A = &AVector[0];
    std::vector<double> BVector(9);
    double* B = &BVector[0];
    std::vector<double> CVector(9);
    double* C = &CVector[0];
    std::vector<double> testVector(9);
    double* Ctest = &testVector[0];
    int num = 9;
    for (int n=0; n<9; n++){
        A[n] = num*num-2;
        B[n] = num;
        C[n] = B[n]-A[n];
        num++;
        }
      
    CORRESPONDENCE::DIFFTENSOR(A, B, Ctest);
    
     for (int n=0; n<9; n++){
        TEST_FLOATING_EQUALITY(C[n],Ctest[n],tolerance);
        }

}
TEUCHOS_UNIT_TEST(correspondence, GetVoigtNotation) {

    const double tolerance = 1.0e-15;
    std::vector<double> AVector(9);
    double* A = &AVector[0];
    std::vector<double> BVector(6);
    double* B = &BVector[0];
    std::vector<double> testVector(6);
    double* Btest = &testVector[0];
    int num = 9;
    for (int n=0; n<9; n++){
        A[n] = num*num-2;
        num++;
        }
    B[0] = A[0];
    B[1] = A[4];
    B[2] = A[8];
    B[3] = 0.5*(A[5]+A[7]);
    B[4] = 0.5*(A[2]+A[6]);
    B[5] = 0.5*(A[1]+A[3]);
    CORRESPONDENCE::GetVoigtNotation(A, Btest);
    
     for (int n=0; n<6; n++){
        TEST_FLOATING_EQUALITY(B[n],Btest[n],tolerance);
        }

}
TEUCHOS_UNIT_TEST(correspondence, GetTensorFromVoigtNotation) {
    const double tolerance = 1.0e-15;
    std::vector<double> AVector(6);
    double* A = &AVector[0];
    std::vector<double> BVector(9);
    double* B = &BVector[0];
    std::vector<double> testVector(9);
    double* Btest = &testVector[0];
    int num = 9;
    for (int n=0; n<6; n++){
        A[n] = num*num-2;
        num++;
        }
    B[0]= A[0];
    B[1]= A[5];
    B[2]= A[4];
    B[3]= A[5];
    B[4]= A[1];
    B[5]= A[3];
    B[6]= A[4];
    B[7]= A[3];
    B[8]= A[2];
    CORRESPONDENCE::GetTensorFromVoigtNotation(A, Btest);
    
    for (int n=0; n<9; n++){
        TEST_FLOATING_EQUALITY(B[n],B[n],tolerance);
        }
}

int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
