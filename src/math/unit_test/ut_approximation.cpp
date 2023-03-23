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
#include "approximation.h"

using namespace std;
using namespace Teuchos;


TEUCHOS_UNIT_TEST(approximation, get_field_size) {
   
    std::vector<int> AVector(20);
    int* A = &AVector[0];
    int n_cont = 4;
    // x-rotation
    A[0] = 5;
    A[6] = 6;
    
    int check = APPROXIMATION::get_field_size(1,A,n_cont);
    TEST_EQUALITY(check,10); 
    check = APPROXIMATION::get_field_size(2,A,n_cont);
    TEST_EQUALITY(check,21);
    check = APPROXIMATION::get_field_size(2,A,n_cont+5);
    TEST_EQUALITY(check,31); 
    A[6] = 0;
    check = APPROXIMATION::get_field_size(2,A,n_cont);
    TEST_EQUALITY(check,15); 
}

TEUCHOS_UNIT_TEST(approximation, get_sample_weighted) {
    const double tolerance = 2.0e-8;
    TEST_FLOATING_EQUALITY(APPROXIMATION::get_sample_weighted(0,1,0),0,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::get_sample_weighted(1,2,0),0.5,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::get_sample_weighted(0,0,-2),1,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::get_sample_weighted(7.5,0.5,-2),3.8,tolerance);
}

TEUCHOS_UNIT_TEST(approximation, basis_func) {
    const double tolerance = 2.0e-7;
    
    int p = 2, n = 4;
    std::vector<double> kVector(p+n+1);
    double* k = &kVector[0];
    APPROXIMATION::knots(n, p, true, k);
    // caclulated with Python
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(0,p,k,0.2),0.36,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(1,p,k,0.2),0.56,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(2,p,k,0.2),0.08,tolerance);
    
    p = 4;
    n = 7;
    std::vector<double> kVector2(p+n+1);
    double* k2 = &kVector2[0];
    APPROXIMATION::knots(n, p, true, k2);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(0,p,k2,0.1),0.2401,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(1,p,k2,0.2),0.429,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(1,p,k2,0.4),0.0512,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(1,p,k2,1.0),0.0,tolerance);
}

TEUCHOS_UNIT_TEST(approximation, knots) {
    const double tolerance = 4.0e-8;
    int n = 3;
    int p = 2;
    std::vector<double> kVector(20);
    double* k = &kVector[0];
    int num_knots = int(n + p + 1);
    double ktest[6] = {0.,  0., 0., 1., 1., 1.};
    APPROXIMATION::knots(n,p,true,k);
    for (n=0; n<num_knots; n++) TEST_FLOATING_EQUALITY(k[n],ktest[n],tolerance);  
    
    double ktest2[6] = {0.,  1./(num_knots  - 1), 2./(num_knots  - 1), 3./(num_knots  - 1), 4./(num_knots  - 1), 5./(num_knots  - 1)};
    n = 3;
    p = 2;
    APPROXIMATION::knots(n,p,false,k);
    for (n=0; n<num_knots; n++) TEST_FLOATING_EQUALITY(k[n],ktest2[n],tolerance);  

    n = 4;
    p = 2;
    num_knots = int(n + p + 1);
    double ktest3[7] = {0.,  0., 0., 0.5, 1., 1., 1.};
    APPROXIMATION::knots(n,p,true,k);
    for (n=0; n<num_knots; n++) TEST_FLOATING_EQUALITY(k[n],ktest3[n],tolerance);  
    n = 4;
    p = 3;
    num_knots = int(n + p + 1);
    double ktest4[8] = {0,  0, 0, 0, 1, 1, 1, 1};
    APPROXIMATION::knots(n,p,true,k);
    for (n=0; n<num_knots; n++) TEST_FLOATING_EQUALITY(k[n],ktest4[n],tolerance);  
    n = 6;
    p = 2;
    num_knots = int(n + p + 1);
    double ktest5[9] = {0,  0, 0, 0.25, 0.5, 0.75, 1 , 1, 1};
    APPROXIMATION::knots(n,p,true,k);
    for (n=0; n<num_knots; n++) TEST_FLOATING_EQUALITY(k[n],ktest5[n],tolerance);  
}
int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}