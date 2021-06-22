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
#include "correspondence.h"

using namespace std;
using namespace Teuchos;

TEUCHOS_UNIT_TEST(correspondence, MatMul) {

    double A[6][6];
    double B[6][6];
    double C[6][6] = {{441, 462, 483, 504, 525, 546},
    {1017, 1074, 1131, 1188, 1245, 1302},
    {1593, 1686, 1779, 1872, 1965, 2058},
    {2169, 2298, 2427, 2556, 2685, 2814},
    {2745, 2910, 3075, 3240, 3405, 3570},
    {3321, 3522, 3723, 3924, 4125, 4326}};

    double Ctrans[6][6] = {{2166, 2262, 2358, 2454, 2550, 2646},
    {2262, 2364, 2466, 2568, 2670, 2772},
    {2358, 2466, 2574, 2682, 2790, 2898},
    {2454, 2568, 2682, 2796, 2910, 3024},
    {2550, 2670, 2790, 2910, 3030, 3150},
    {2646, 2772, 2898, 3024, 3150, 3276}};
        
    double Ctest[6][6];
    double Ctranstest[6][6];
    int n, m, x;
    x=1;

    for (n=0; n<6; n++)
        for (m=0; m<6; m++)
        {
            A[n][m]=x;
            B[n][m]=x;
            x++;
        }

    CORRESPONDENCE::MatMul<double>(6,A,B,Ctest,false);
    CORRESPONDENCE::MatMul<double>(6,A,B,Ctranstest,true);
    
    for (n=0; n<6; n++)
        for (m=0; m<6; m++)
        {
            TEST_FLOATING_EQUALITY(Ctest[n][m],C[n][m],0);
            TEST_FLOATING_EQUALITY(Ctranstest[n][m],Ctrans[n][m],0);
        }
}

TEUCHOS_UNIT_TEST(correspondence, Invert3by3Matrix) {

    std::vector<double> AVector(9);
    double* A = &AVector[0];
    double determinant;
    std::vector<double> CVector(9);
    double* C = &CVector[0];
    std::vector<double> CtestVector(9);
    double* Ctest = &CtestVector[0];
    int n;

    *(A)=1;
    *(A+1)=2;
    *(A+2)=3;
    *(A+3)=2;
    *(A+4)=3;
    *(A+5)=4;
    *(A+6)=1;
    *(A+7)=2;
    *(A+8)=4;
    
    *(C)=-4;
    *(C+1)=2;
    *(C+2)=1;
    *(C+3)=4;
    *(C+4)=-1;
    *(C+5)=-2;
    *(C+6)=-1;
    *(C+7)=0;
    *(C+8)=1;

    CORRESPONDENCE::Invert3by3Matrix(A,determinant,Ctest);
    for (n=0; n<9; n++)
    {
        TEST_FLOATING_EQUALITY(*(Ctest+n),*(C+n),0);
    }
}

TEUCHOS_UNIT_TEST(correspondence, Invert2by2Matrix) {

    std::vector<double> AVector(9);
    double* A = &AVector[0];
    double determinant;
    std::vector<double> CVector(9);
    double* C = &CVector[0];
    std::vector<double> CtestVector(9);
    double* Ctest = &CtestVector[0];
    int n;

    *(A)=1;
    *(A+1)=2;
    *(A+3)=3;
    *(A+4)=4;
    
    *(C)=-2;
    *(C+1)=1;
    *(C+3)=1.5;
    *(C+4)=-0.5;

    CORRESPONDENCE::Invert2by2Matrix(A,determinant,Ctest);
    for (n=0; n<4; n++)
    {
        TEST_FLOATING_EQUALITY(*(Ctest+n),*(C+n),0);
    }
}

TEUCHOS_UNIT_TEST(correspondence, MatrixMultiply3x3) {

    double A[3][3];
    double B[3][3];
    double C[3][3] = {{30, 36, 42},
    {66, 81, 96},
    {102, 126, 150}};
    double Ctest[3][3];
    int n, m, x;
    x=1;

    for (n=0; n<3; n++)
        for (m=0; m<3; m++)
        {
            A[n][m]=x;
            B[n][m]=x;
            x++;
        }

    CORRESPONDENCE::MatrixMultiply3x3(A,B,Ctest);
    for (n=0; n<3; n++)
        for (m=0; m<3; m++)
        {
            TEST_FLOATING_EQUALITY(Ctest[n][m],C[n][m],0);
        }
}

int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
