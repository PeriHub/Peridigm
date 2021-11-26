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
TEUCHOS_UNIT_TEST(FEM, shapeFunctionsLagrange) {
        const double tolerance = 1.0e-15;
        std::vector<double> NmatrixTestVector(8);
        double* NmatrixTest = &NmatrixTestVector[0];
        double Nmatrix[9];
        const double elCoor[3]

        double elCoor[3],weights[3], order[3];

        if (order == 1){
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
                        TEST_FLOATING_EQUALITY(N[n],NmatrixTest[n],tolerance);
                        TEST_FLOATING_EQUALITY(weights[n],weightsTest[n],tolerance);
                }

        }
    

    

}

int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
