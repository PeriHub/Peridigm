/*! \file ut_FEM_routines.cpp */

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
//
// funded by dfg project reference
// Licence agreement
//
// Christian Willberg    christian.willberg@dlr.de
//@HEADER

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "FEM_routines.h"

using namespace std;
using namespace Teuchos;

TEUCHOS_UNIT_TEST(FEM, weightsAndIntegrationPoints) {
    const double tolerance = 1.0e-15;
    int order;
    std::vector<double> elCoorTestVector(3);
    double* elCoorTest = &elCoorTestVector[0];
    std::vector<double> weightsTestVector(3);   
    double* weightsTest = &weightsTestVector[0];
    double elCoor[3],weights[3];
    if (order == 0){
        elCoor[0] = 0;
        weights[0] = 2;
        FEM::weightsAndIntegrationPoints(order, elCoorTest, weightsTest);
        TEST_FLOATING_EQUALITY(elCoor[0],elCoorTest[0],tolerance);
        TEST_FLOATING_EQUALITY(weights[0],weightsTest[0],tolerance);

        }
    if (order == 1){
        elCoor[0] = -sqrt(1.0/3.0);
        weights[0] = 1;
        elCoor[1] =  sqrt(1.0/3.0);
        weights[1] = 1;
        FEM::weightsAndIntegrationPoints(order, elCoorTest, weightsTest);
        for (int n=0; n<2; n++){
                TEST_FLOATING_EQUALITY(elCoor[n],elCoorTest[n],tolerance);
                TEST_FLOATING_EQUALITY(weights[n],weightsTest[n],tolerance);
        }
    }                        
    else if (order==2){
        elCoor[0] = -sqrt(3.0/5.0);
        weights[0] = 5.0/9.0;
        elCoor[1] =  0.0;
        weights[1] = sqrt(1.0/3.0);
        elCoor[2] =  sqrt(3.0/5.0);
        weights[2] = 5.0/9.0;
        FEM::weightsAndIntegrationPoints(order, elCoorTest, weightsTest);
        for (int n=0; n<2; n++){
                TEST_FLOATING_EQUALITY(elCoor[n],elCoorTest[n],tolerance);
                TEST_FLOATING_EQUALITY(weights[n],weightsTest[n],tolerance);
        }
    } 

}

TEUCHOS_UNIT_TEST(FEM, setElementMatrices) {
    double tolerance = 1e-14;
    std::vector<int> orderVector(3);
    int* order = &orderVector[0];
    order[0] = 1;order[1] = 1;order[2] = 1;
    std::vector<double> AVector(4);
    double* A = &AVector[0];
    A[0] = 1;A[1] = 1;A[2] = 1; A[3] = -1.2;
    std::vector<double> BVector(4);
    double* B = &BVector[0];
    B[0] = 1;B[1] = 1;B[2] = 1;B[3] = 5;
    std::vector<double> CVector(4);
    double* C = &CVector[0];
    C[0] = 1;C[1] = 1;C[2] = 1;C[3] = 5;
    std::vector<double> BxTestVector(16), ByTestVector(16), BzTestVector(16);
    double* BxTest = &BxTestVector[0];
    double* ByTest = &ByTestVector[0];
    double* BzTest = &BzTestVector[0];

    bool twoD = true;
    //to be filled;
    double Bx[16] ={1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0}, By[16] ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, Bz[16] ={0,0,0,0,0,0,0,0},0,0,0,0,0,0,0,0;
    if (twoD){
        for (int iID=0 ; iID<2 ; ++iID){
            FEM::setElementMatrices(twoD, 4*iID, order, A, B, C, B, C, A, BxTest, BzTest, BzTest);
        }
        for (int n=0; n<8; n++){
            TEST_FLOATING_EQUALITY(BxTest[n],Bx[n],tolerance);
            TEST_FLOATING_EQUALITY(ByTest[n],By[n],tolerance);
        }
    }
    else
    {
        for (int iID=0 ; iID<2 ; ++iID){
            FEM::setElementMatrices(twoD, 8*iID, order, A, B, C, B, C, A, BxTest, BzTest, BzTest);
        }
        for (int n=0; n<8; n++){
            TEST_FLOATING_EQUALITY(BxTest[n],Bx[n],tolerance);
            TEST_FLOATING_EQUALITY(ByTest[n],By[n],tolerance);
        }

    }
    
}
TEUCHOS_UNIT_TEST(FEM, shapeFunctionsLagrange) {
        const double tolerance = 1.0e-15;
        std::vector<double> NmatrixTestVector(8);
        double* NmatrixTest = &NmatrixTestVector[0];
        double Nmatrix[9];
        double elCoor[3];
        int order[3];
        order[0]=1;
        order[1]=1;
        order[2]=1;
        if (order[0] == 1){
                Nmatrix[0]= 0.09806875 ;
                Nmatrix[1]= 0.13268124999999997 ;
                Nmatrix[2]= 0.07144374999999999 ;
                Nmatrix[3]= 0.05280624999999999 ;
                Nmatrix[4]= 0.17818125 ;
                Nmatrix[5]= 0.24106875 ;
                Nmatrix[6]= 0.12980624999999998 ;
                Nmatrix[7]= 0.09594375 ;
                elCoor[0] = 0.15;
                elCoor[1] = -0.3;
                elCoor[2] = 0.29;

                order[0] = 1;
                order[1] = 1;
                order[2] = 1;
                FEM::shapeFunctionsLagrange(NmatrixTest, order, elCoor);
                for (int n=0; n<2; n++){
                        TEST_FLOATING_EQUALITY(Nmatrix[n],NmatrixTest[n],tolerance);
                       
                }

        }
 

}

TEUCHOS_UNIT_TEST(FEM, defineLagrangianGridSpace) {
    const double tolerance = 1.0e-15;
    std::vector<double> BTestVector(6);
    double* BTest = &BTestVector[0];
    double B[6];
    int A;
    A = 1;
    FEM::defineLagrangianGridSpace(A, BTest);
    B[0] = -1;
    B[1] =  1;
    for (int n=0; n<A+1; n++){
        TEST_FLOATING_EQUALITY(B[n],BTest[n],tolerance);
    }
    A = 5;
    B[0] = -1;
    B[1] = -0.6;
    B[2] = -0.2;
    B[3] = 0.2;
    B[4] = 0.6;
    B[5] =  1;
    FEM::defineLagrangianGridSpace(A, BTest);
    for (int n=0; n<A+1; n++){
        TEST_FLOATING_EQUALITY(B[n],BTest[n],tolerance);
    }
}

TEUCHOS_UNIT_TEST(FEM, getDisplacements) {
    const double tolerance = 1.0e-15;
    int A = 2;

    double B[6] = {0,0.1,-0.3,100,20.25,11};
    double C[6] = {0,-0.1,-0.3,1,0,12};
    double D[6] = {0,-0.2,0.0,-99,-20.25,1};
    std::vector<double> DVectorTest(6);
    double* DTest = &DVectorTest[0];

    FEM::getDisplacements(A, B, C, D);
 
    for (int n=0; n<3*A+2; n++){
        TEST_FLOATING_EQUALITY(D[n],DTest[n],tolerance);
    }
    
}
TEUCHOS_UNIT_TEST(FEM, getNumberOfIntegrationPoints) {
    bool A = true;
    int B[3] = [1,1,1];

    int Ctest = FEM::getNumberOfIntegrationPoints(A, B);
    int C = 8
    TEST_FLOATING_EQUALITY(C, Ctest, 0); 

    A = false;   
    Ctest = FEM::getNumberOfIntegrationPoints(A, B);
    C = 4; 
    TEST_FLOATING_EQUALITY(C, Ctest, 0); 
}
int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
