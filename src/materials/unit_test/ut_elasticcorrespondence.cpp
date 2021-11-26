/*! \file ut_correspondence.cpp */

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
#include "elastic_correspondence.h"

using namespace std;
using namespace Teuchos;

TEUCHOS_UNIT_TEST(correspondence, updateElasticCauchyStressAnisotropicCode) {
    const double tolerance = 1.0e-15;
    double A[6][6] ;

     int m, n, num;
     num = 1;
     for (m=0; m<6; m++){
             for (n=0; n<6; n++){
             A[m][n] = num;
             num++;
             }
     }   

    std::vector<double> BVector(9);
    double* B = &BVector[0];
    std::vector<double> testVector(9);
    double* Ctest = &testVector[0];

    num = 9;
    for (m=0; m<9; m++){
        B[m] = num;
        num++;
             
     }   
     //updateElasticCauchyStressAnisotropicCode
     //   const ScalarT strain[][3],
     //   ScalarT* sigmaNP1,
     //   const ScalarT C[][6],
     //   int type
    CORRESPONDENCE::updateElasticCauchyStressAnisotropicCode(B, Ctest, A, 0);
    
    double C[] =      { 468.0  ,3978.0, 3276.0,
                        3978.0 ,1170.0, 2574.0,
                        3276.0 ,2574.0, 1872.0};

    for (n=0; n<9; n++){

            TEST_FLOATING_EQUALITY(C[n],Ctest[n],tolerance);
        
    }
    double C2D[] = {167, 1487, 0, 
                    1487, 431., 0, 
                    0,   0,   0};   
    
    CORRESPONDENCE::updateElasticCauchyStressAnisotropicCode(B, Ctest, A, 1);
    
    for (n=0; n<9; n++){
        TEST_FLOATING_EQUALITY(C2D[n],Ctest[n],tolerance);
        }

}



int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
