/*! \file ut_user_material_interface.cpp */

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
#include "user_material_interface_correspondence.h"

using namespace std;
using namespace Teuchos;

TEUCHOS_UNIT_TEST(correspondence, ut_user_material_interface) {
    const double tolerance = 1.0e-15;
    
    std::vector<double> strainLocVoigtVector(6), statevVector(9), depsLocVoigtVector(6), sigmaNP1LocVoigtVector(6), sigmaNP1LocVoigtTestVector(6), coordsVector(6), drotVector(9), drotTestVector(9), defGradNVector(9), defGradNP1Vector(9);
    double* strainLocVoigt = &strainLocVoigtVector[0];
    double* depsLocVoigt = &depsLocVoigtVector[0];
    double* sigmaNP1LocVoigt = &sigmaNP1LocVoigtVector[0];
    double* sigmaNP1LocVoigtTest = &sigmaNP1LocVoigtTestVector[0];
    double* statev = &statevVector[0];
    double* coords = &coordsVector[0];
    double* drot = &drotVector[0];
    double* drotTest = &drotTestVector[0];
    double* defGradN = &defGradNVector[0];
    double* defGradNP1 = &defGradNP1Vector[0];

    const double time = 0.1, dtime = 0.1;
    const double* temp;
    const double* dtemp;
    double DSDDE[6*6], DSDDETest[6*6];
    double DDSDDT[6], DRPLDE[6];
    double SSE = -1,SPD = -1,SCD = -1,RPL = -1, DRPLDT = -1;
    double PREDEF = -1, DPRED = -1, PNEWDT = -1, CELENT = -1;
    const std::string matname = "Elastic";
    char matnameArray[80];
    int nname = 7;
    int nnormal = 3;
    int nshr = 3;
    int nstresscomp = 6; 
    int NOEL;
    int NPT = -1, KSLAY = -1, KSPT = -1, JSTEP = -1, KINC = -1;
    const int nstatev = 9;
    const int nprops = 1;
    const double* props;
    
    nname = matname.length();

    for(unsigned int i=0; i<matname.length(); ++i){
    matnameArray[i] = matname[i];
    }
    
    int m, n, num;
    for (m=0; m<6; m++){
        *(sigmaNP1LocVoigt+m) = m;
        *(sigmaNP1LocVoigtTest+m) = m*2;
    } 
    for (m=0; m<6; m++){
        *(drot+m) = m;
        *(drotTest+m) = m*2;
    } 
    for (m=0; m<6; m++){
        DSDDE[m] = m;
        DRPLDE[m] = m*2;
    }  
    num = 1;
    for (m=0; m<6; m++){
        for (n=0; n<6; n++){
            DSDDE[m*n] = num;
            DSDDETest[m*n] = num*2;
            num++;
        }
    }   

    CORRESPONDENCE::UMATINTTEST(sigmaNP1LocVoigt,statev,DSDDE,&SSE,&SPD,&SCD,&RPL,
    DDSDDT, DRPLDE,&DRPLDT,strainLocVoigt,depsLocVoigt,&time,&dtime,temp,dtemp,
    &PREDEF,&DPRED,matnameArray,&nnormal,&nshr,&nstresscomp,&nstatev,props,
    &nprops,coords,drot,&PNEWDT,&CELENT,defGradN,defGradNP1,
    &NOEL,&NPT,&KSLAY,&KSPT,&JSTEP,&KINC,&nname); 
    
    for (m=0; m<6; m++){
        TEST_FLOATING_EQUALITY(*(sigmaNP1LocVoigt+m),*(sigmaNP1LocVoigtTest+m),tolerance);
    } 
    for (m=0; m<9; m++){
        TEST_FLOATING_EQUALITY(*(drot+m),*(drotTest+m),tolerance);
    } 
    for (m=0; m<6; m++){
        for (n=0; n<6; n++){
            TEST_FLOATING_EQUALITY(DSDDE[m*n],DSDDETest[m*n],tolerance);
        }
    } 
}

int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
