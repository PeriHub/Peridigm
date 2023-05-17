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

TEUCHOS_UNIT_TEST(correspondence, addTemperatureStrain) {
    const double tolerance = 1.0e-15;
    double A[3][3] ;
    double temperature = 1.5;
    int m, n, num;
    num = 1;
    for (m=0; m<3; m++){
             for (n=0; n<3; n++){
             A[m][n] = num;
             num++;
             }
     }   
    

    std::vector<double> BVector(9);
    double* B = &BVector[0];
    std::vector<double> testVector(9);

    num = 9;
    for (m=0; m<9; m++){
        B[m] = num;
        num++;        
     }   

    CORRESPONDENCE::addTemperatureStrain(A,temperature,B);
    num = 1;
    for (m=0; m<3; m++){
             for (n=0; n<3; n++){
             A[m][n] = num;
             num++;
             }
     }   


    double Btest[9];
    Btest[0] = 9  - 1  * temperature;  
    Btest[1] = 10 - 2  * temperature;  
    Btest[2] = 11 - 3  * temperature;
    Btest[3] = 12 - 4  * temperature;
    Btest[4] = 13 - 5  * temperature;
    Btest[5] = 14 - 6  * temperature;
    Btest[6] = 15 - 7  * temperature;
    Btest[7] = 16 - 8  * temperature;
    Btest[8] = 17 - 9  * temperature;


    for (n=0; n<9; n++){

            TEST_FLOATING_EQUALITY(B[n],Btest[n],tolerance);
        
    }


}
TEUCHOS_UNIT_TEST(correspondence, computeGreenLagrangeStrain) {
    const double tolerance = 1.0e-15;
    std::vector<double> defGradVector(9);
    double* defGrad = &defGradVector[0];
    std::vector<double> strainVector(9);
    double* strain = &strainVector[0];
    std::vector<double> strainVectorTest(9);
    double* strainTest = &strainVectorTest[0];
    defGrad[0] = 1.0;
    defGrad[4] = 1.0;
    defGrad[8] = 1.0;
   
    CORRESPONDENCE::computeGreenLagrangeStrain(defGrad,strain);

    for (int n=0; n<9; n++) TEST_FLOATING_EQUALITY(strain[n],0.0,tolerance);
        
    defGrad[0] = 2.0;defGrad[1] = 1.0;defGrad[2] = 2.0;
    defGrad[3] = 2.0;defGrad[4] = 1.0;defGrad[5] = 2.3;
    defGrad[6] = 2.0;defGrad[7] = -1.0;defGrad[8] = 3.0;
   
    strainTest[0] = 11.0;strainTest[1] = 2.0;strainTest[2] = 14.6;
    strainTest[3] = 2.0;strainTest[4] = 2.0;strainTest[5] = 1.3;
    strainTest[6] = 14.6;strainTest[7] = 1.3;strainTest[8] = 17.29;
    CORRESPONDENCE::computeGreenLagrangeStrain(defGrad,strain);

    for (int n=0; n<9; n++) TEST_FLOATING_EQUALITY(strain[n],0.5*strainTest[n],tolerance);

}


int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
