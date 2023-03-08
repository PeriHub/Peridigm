/*! \file ut_ealsticcorrespondence.cpp */

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
#include "material_utilities.h"

using namespace std;
using namespace Teuchos;

TEUCHOS_UNIT_TEST(material_utilities, getDiffAndLen) {
    const double tolerance = 1.0e-15;
    std::vector<double> AVector(6);
    double* A = &AVector[0];
    std::vector<double> BVector(6);
    double* B = &BVector[0];
    std::vector<double> testVector(6);
    double* Ctest = &testVector[0];
    double d;
    int num = 1;
    for (int m=0; m<num; m++){
        A[m] = -2;
        B[m] =  2;             
     }   
     //updateElasticCauchyStressAnisotropicCode
     //   const ScalarT strain[][3],
     //   ScalarT* sigmaNP1,
     //   const ScalarT C[][6],
     //   int type
    d = MATERIAL_EVALUATION::getDiffAndLen(A,B,num,Ctest);
    TEST_FLOATING_EQUALITY(d,4,tolerance);
    TEST_FLOATING_EQUALITY(4.0,Ctest[0],tolerance);
    num = 3;
    for (int m=0; m<num; m++){
        A[m] = 2*m;
        B[m] = 1+m;             
     }   
    double C[] =      {1,0,-1};
    d = MATERIAL_EVALUATION::getDiffAndLen(A,B,num,Ctest);
    for (int m=0; m<num; m++) TEST_FLOATING_EQUALITY(C[m],Ctest[m],tolerance);    
    TEST_FLOATING_EQUALITY(d,sqrt(2.0),tolerance);
    

}
TEUCHOS_UNIT_TEST(material_utilities, getStretch) {
    const double tolerance = 1.0e-15;
    double A;

    A = MATERIAL_EVALUATION::getStretch(5.5,5.5);
    TEST_FLOATING_EQUALITY(A,0.0,tolerance);
    A = MATERIAL_EVALUATION::getStretch(5.6,-5.6);
    TEST_FLOATING_EQUALITY(A,-11.2,tolerance);
    A = MATERIAL_EVALUATION::getStretch(-5.6,5.6);
    TEST_FLOATING_EQUALITY(A, 11.2,tolerance);
    A = MATERIAL_EVALUATION::getStretch(2.0,0.0);
    TEST_FLOATING_EQUALITY(A,-2.0,tolerance);
    A = MATERIAL_EVALUATION::getStretch(0.0,2.0);
    TEST_FLOATING_EQUALITY(A,2.0,tolerance);

}
TEUCHOS_UNIT_TEST(material_utilities, getProjectedForces) {
    const double tolerance = 1.0e-15;
    double t;
    std::vector<double> AVector(3);
    double* A = &AVector[0]; 
    int dof = 3;
    double B;
    std::vector<double> fVector(3);
    double* f = &fVector[0]; 
    std::vector<double> fTestVector(3);
    double* fTest = &fTestVector[0]; 

    A[0] = 1; A[1]=2; A[2]=-1;
    B = sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
    t = 10;
    MATERIAL_EVALUATION::getProjectedForces(t,A,B,dof,f);
    fTest[0] = 10/B; fTest[1] = 2*10/B; fTest[2] = -10/B;
    for (int m=0; m<dof; m++) TEST_FLOATING_EQUALITY(f[m],fTest[m],tolerance);    
    dof = 2;
    A[0] = -1; A[1]=5;
    B = sqrt(A[0]*A[0]+A[1]*A[1]);
    t = 10;
    MATERIAL_EVALUATION::getProjectedForces(t,A,B,dof,f);
    fTest[0] = -10/B; fTest[1] = 5*10/B; fTest[2] = 0.0;
    for (int m=0; m<dof; m++) TEST_FLOATING_EQUALITY(f[m],fTest[m],tolerance);    


}




int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
