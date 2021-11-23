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

    double B[3][3];
    std::vector<double> testVector(9);
    double* Ctest = &testVector[0];

    num = 9;
    for (m=0; m<3; m++){
             for (n=0; n<3; n++){
             B[m][n] = num;
             num++;
             }
     }   
     //updateElasticCauchyStressAnisotropicCode
     //   const ScalarT strain[][3],
     //   ScalarT* sigmaNP1,
     //   const ScalarT C[][6],
     //   int type
    CORRESPONDENCE::updateElasticCauchyStressAnisotropicCode(B, Ctest, A, 0);
    
    double C[] = {484,  2590, 3292, 
                  2590, 1186, 3994, 
                  3292, 2590, 1888};

    for (n=0; n<9; n++){

            TEST_FLOATING_EQUALITY(C[n],Ctest[n],tolerance);
        
    }
    double C2D[] = {123, 915, 0, 
                    915, 387, 0, 
                    0,   0,   0};   
    
    CORRESPONDENCE::updateElasticCauchyStressAnisotropicCode(B, Ctest, A, 1);
    
    for (n=0; n<9; n++){
        TEST_FLOATING_EQUALITY(C2D[n],Ctest[n],tolerance);
        }

}

TEUCHOS_UNIT_TEST(correspondence, createRotationMatrix) {
    const double tolerance = 1.0e-15;
        
    int m,n;
    double Btest[3][3];
    double A[3];
    // x-rotation
    A[0] = 14;A[1] = 0;A[2] = 0;
    CORRESPONDENCE::createRotationMatrix(A,Btest);
    
    double Bx[3][3] = {{1,  0, 0}, 
                  {0, 0.97029573, -0.2419219}, 
                  {0, 0.2419219, 0.97029573}};

    for (n=0; n<3; n++){
            for (m=0; m<3; m++){
                TEST_FLOATING_EQUALITY(Bx[n][m],Btest[n][m],tolerance); 
            }
    }
    // x-rotation
    A[0] = 0;A[1] = 70;A[2] = 0;
    CORRESPONDENCE::createRotationMatrix(A,Btest);
    
    double By[3][3] = {{0.34202014,  0, 0.93969262}, 
         {   0, 1, 0}, 
         {  -0.93969262, 0, 0.34202014}};

    for (n=0; n<3; n++){
            for (m=0; m<3; m++){
                TEST_FLOATING_EQUALITY(By[n][m], Btest[n][m],tolerance); 
            }
    }
    // z-rotation
    A[0] = 0;A[1] = 0;A[2] = 30;
    CORRESPONDENCE::createRotationMatrix(A,Btest);
    
    double Bz[3][3] = {{0.8660254,  -0.5, 0}, 
         {  0.5, 0.8660254, 0}, 
         {  0, 0, 1}};

    for (n=0; n<3; n++){
            for (m=0; m<3; m++){
                TEST_FLOATING_EQUALITY(Bz[n][m],Btest[n][m],tolerance); 
            }
    }
}

int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
