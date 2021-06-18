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

TEUCHOS_UNIT_TEST(correspondence, SetTest) {

    double A[6][6];
    double B[6][6];
    double C[6][6] = {{441, 1017, 1593, 2169, 2745, 3321},
    {462, 1074, 1686, 2298, 2910, 3522},
    {483, 1131, 1779, 2427, 3075, 3723},
    {504, 1188, 1872, 2556, 3240, 3924},
    {525, 1245, 1965, 2685, 3405, 4125},
    {546, 1302, 2058, 2814, 3570, 4326}};
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

    CORRESPONDENCE::MatMul<double>(6,A,B,Ctest,true);
    std::cout<< C << Ctest << std::endl;
    for (n=0; n<6; n++)
        for (m=0; m<6; m++)
        {
            TEST_ASSERT(Ctest[n][m]==C[n][m]);
        }
}



int main( int argc, char* argv[] ) {
  
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
