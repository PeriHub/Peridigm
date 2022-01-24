/*! \file ut_matrices.cpp */

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
#include "matrices.h"

using namespace std;
using namespace Teuchos;
TEUCHOS_UNIT_TEST(matrices, createRotationMatrix) {
    const double tolerance = 2.0e-8;
        
    int m,n;
    double Btest[3][3];
    double A[3];
    // x-rotation
    A[0] = 14;A[1] = 0;A[2] = 0;
    MATRICES::createRotationMatrix(A,Btest);
    
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
    MATRICES::createRotationMatrix(A,Btest);
    
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
    MATRICES::createRotationMatrix(A,Btest);
    
    double Bz[3][3] = {{0.8660254,  -0.5, 0}, 
         {  0.5, 0.8660254, 0}, 
         {  0, 0, 1}};

    for (n=0; n<3; n++){
            for (m=0; m<3; m++){
                TEST_FLOATING_EQUALITY(Bz[n][m],Btest[n][m],tolerance); 
            }
    }
}
TEUCHOS_UNIT_TEST(matrices, MatMul) {

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
    int n, m, x;
    x=1;

    for (n=0; n<6; n++)
        for (m=0; m<6; m++)
        {
            A[n][m]=x;
            B[n][m]=x;
            x++;
        }

    MATRICES::MatMul<double>(6,A,B,Ctest,false);
    
    for (n=0; n<6; n++)
        for (m=0; m<6; m++)
        {
            TEST_FLOATING_EQUALITY(Ctest[n][m],C[n][m],0);
        }
    
    MATRICES::MatMul<double>(6,A,B,Ctest,true);

    for (n=0; n<6; n++)
        for (m=0; m<6; m++)
        {
            TEST_FLOATING_EQUALITY(Ctest[n][m],Ctrans[n][m],0);
        }
    
}

TEUCHOS_UNIT_TEST(matrices, Invert3by3Matrix) {

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

    MATRICES::Invert3by3Matrix(A,determinant,Ctest);
    for (n=0; n<9; n++)
    {
        TEST_FLOATING_EQUALITY(*(Ctest+n),*(C+n),0);
    }
}

TEUCHOS_UNIT_TEST(matrices, Invert2by2Matrix) {

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

    MATRICES::Invert2by2Matrix(A,determinant,Ctest);
    for (n=0; n<4; n++)
    {
        TEST_FLOATING_EQUALITY(*(Ctest+n),*(C+n),0);
    }
}

TEUCHOS_UNIT_TEST(matrices, TransposeMatrix) {

    std::vector<double> AVector(9);
    double* A = &AVector[0];
    std::vector<double> CVector(9);
    double* C = &CVector[0];
    std::vector<double> CtestVector(9);
    double* Ctest = &CtestVector[0];
    int n;

    for (n=0; n<9; n++)
    {
        *(A+n)=n+1;
    }
    
    *(C)=1;
    *(C+1)=4;
    *(C+2)=7;
    *(C+3)=2;
    *(C+4)=5;
    *(C+5)=8;
    *(C+6)=3;
    *(C+7)=6;
    *(C+8)=9;

    MATRICES::TransposeMatrix(A,Ctest);
    for (n=0; n<9; n++)
    {
        TEST_FLOATING_EQUALITY(*(Ctest+n),*(C+n),0);
    }
}

TEUCHOS_UNIT_TEST(matrices, MatrixMultiply3x3) {

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

    MATRICES::MatrixMultiply3x3(A,B,Ctest);
    for (n=0; n<3; n++)
        for (m=0; m<3; m++)
        {
            TEST_FLOATING_EQUALITY(Ctest[n][m],C[n][m],0);
        }
}

TEUCHOS_UNIT_TEST(matrices, MatrixMultiply3x3toVector) {

    double A[3][3];
    double B[3][3];
    double C[9] = {30, 36, 42, 66, 81, 96,102, 126, 150};
    std::vector<double> testVector(9);
    double* Ctest = &testVector[0];
    int n, m, x;
    x=1;

    for (n=0; n<3; n++)
        for (m=0; m<3; m++)
        {
            A[n][m]=x;
            B[n][m]=x;
            x++;
        }

    MATRICES::MatrixMultiply3x3toVector(A,B,Ctest);
    for (n=0; n<9; n++){
            TEST_FLOATING_EQUALITY(Ctest[n],C[n],0);
        }
}

TEUCHOS_UNIT_TEST(matrices, MatrixMultiply3x3fromVector) {

    double A[3][3];
    double C[3][3] = {{30, 36, 42},
                     {66, 81, 96},
                     {102, 126, 150}};
    std::vector<double> BVector(9);
    double* B = &BVector[0];
    double Ctest[3][3];
    int n, m, x;
    x=1;

    for (n=0; n<3; n++){
        for (m=0; m<3; m++)
        {
            A[n][m]=x;
            B[x-1]=x;
            x++;
        }
    }
    MATRICES::MatrixMultiply3x3fromVector(A,B,Ctest);
    for (n=0; n<3; n++){
        for (m=0; m<3; m++){
            TEST_FLOATING_EQUALITY(Ctest[n][m],C[n][m],0);
        }
        }
}

TEUCHOS_UNIT_TEST(matrices, MatrixMultiply) {

    std::vector<double> AVector(9);
    double* A = &AVector[0];
    std::vector<double> BVector(9);
    double* B = &BVector[0];
    std::vector<double> CVector(9);
    double* C = &CVector[0];
    std::vector<double> CTransAVector(9);
    double* CTransA = &CTransAVector[0];
    std::vector<double> CTransBVector(9);
    double* CTransB = &CTransBVector[0];
    std::vector<double> CTransABVector(9);
    double* CTransAB = &CTransABVector[0];
    std::vector<double> CtestVector(9);
    double* Ctest = &CtestVector[0];
    int n;

    for (n=0; n<9; n++)
    {
        *(A+n)=n+1;
        *(B+n)=n+1;
    }
    
    *(C)=-60;
    *(C+1)=-72;
    *(C+2)=-84;
    *(C+3)=-132;
    *(C+4)=-162;
    *(C+5)=-192;
    *(C+6)=-204;
    *(C+7)=-252;
    *(C+8)=-300;
    
    *(CTransA)=-132;
    *(CTransA+1)=-156;
    *(CTransA+2)=-180;
    *(CTransA+3)=-156;
    *(CTransA+4)=-186;
    *(CTransA+5)=-216;
    *(CTransA+6)=-180;
    *(CTransA+7)=-216;
    *(CTransA+8)=-252;
    
    *(CTransB)=-28;
    *(CTransB+1)=-64;
    *(CTransB+2)=-100;
    *(CTransB+3)=-64;
    *(CTransB+4)=-154;
    *(CTransB+5)=-244;
    *(CTransB+6)=-100;
    *(CTransB+7)=-244;
    *(CTransB+8)=-388;
    
    *(CTransAB)=-60;
    *(CTransAB+1)=-132;
    *(CTransAB+2)=-204;
    *(CTransAB+3)=-72;
    *(CTransAB+4)=-162;
    *(CTransAB+5)=-252;
    *(CTransAB+6)=-84;
    *(CTransAB+7)=-192;
    *(CTransAB+8)=-300;

    MATRICES::MatrixMultiply(false,false,-2.0,A,B,Ctest);

    for (n=0; n<9; n++)
    {
        TEST_FLOATING_EQUALITY(*(Ctest+n),*(C+n),0);
    }
    
    MATRICES::MatrixMultiply(true,false,-2.0,A,B,Ctest);

    for (n=0; n<9; n++)
    {
        TEST_FLOATING_EQUALITY(*(Ctest+n),*(CTransA+n),0);
    }
    
    MATRICES::MatrixMultiply(false,true,-2.0,A,B,Ctest);

    for (n=0; n<9; n++)
    {
        TEST_FLOATING_EQUALITY(*(Ctest+n),*(CTransB+n),0);
    }
    
    MATRICES::MatrixMultiply(true,true,-2.0,A,B,Ctest);

    for (n=0; n<9; n++)
    {
        TEST_FLOATING_EQUALITY(*(Ctest+n),*(CTransAB+n),0);
    }
}

int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
