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
#include "correspondence.h"

using namespace std;
using namespace Teuchos;

TEUCHOS_UNIT_TEST(correspondence, ut_user_material_interface) {
    const double tolerance = 1.0e-15;
    const int nstatev = 6;
    std::vector<double> sigmaNP1LocVector(6), strainLocVector(6), depsLocVector(6), statevVector(nstatev), coordsVector(3), drotVector(9), defGradNVector(9), defGradNP1Vector(9);
    std::vector<double> sigmaNP1LocTestVector(6), strainLocTestVector(6), depsLocTestVector(6), statevTestVector(6), coordsTestVector(3), drotTestVector(9), defGradNTestVector(9), defGradNP1TestVector(9);
    double* sigmaNP1Loc = &sigmaNP1LocVector[0];
    double* sigmaNP1LocTest = &sigmaNP1LocTestVector[0];
    double* strainLoc= &strainLocVector[0];
    double* strainLocTest = &strainLocTestVector[0];
    double* depsLoc = &depsLocVector[0];
    double* depsLocTest = &depsLocTestVector[0];
    double* statev = &statevVector[0];
    double* statevTest = &statevTestVector[0];
    double* coords = &coordsVector[0];
    double* coordsTest = &coordsTestVector[0];
    double* drot = &drotVector[0];
    double* drotTest = &drotTestVector[0];
    double* defGradN = &defGradNVector[0];
    double* defGradNP1 = &defGradNP1Vector[0];
    double* defGradNTest = &defGradNTestVector[0];
    double* defGradNP1Test = &defGradNP1TestVector[0];

    const double time = 0.1, dtime = 0.1;
    const double timeArray[2] = {time,time};
    const double temp = 200.0;
    const double dtemp = 50.0;
    double DDSDDE[6*6], DDSDDETest[6*6];
    double DDSDDT[6], DRPLDE[6], DDSDDTTest[6], DRPLDETest[6];
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

    const int nprops = 1;
    const double props = 3;
    
    nname = matname.length();

    for(unsigned int i=0; i<matname.length(); ++i){
        matnameArray[i] = matname[i];
    }

    int m, n, num;
    for (m=0; m<9; m++){
        *(drot+m) = m;
        *(defGradN+m) = m;
        *(defGradNP1+m) = m;
        *(drotTest+m) = m*props;
        *(defGradNTest+m) = m*props;
        *(defGradNP1Test+m) = m*props;
    } 
    for (m=0; m<6; m++){
        *(sigmaNP1Loc+m) = m;
        *(sigmaNP1LocTest+m) = m*props;
        *(statev+m) = m;
        *(statevTest+m) = m*props;
        *(strainLoc+m) = m;
        *(strainLocTest+m) = m*props;
        *(depsLoc+m) = m;
        *(depsLocTest+m) = m*props;
        DDSDDT[m] = m;
        DDSDDTTest[m] = m*props;
        DRPLDE[m] = m;
        DRPLDETest[m] = m*props;
    }
    for (m=0; m<3; m++){
        *(coords+m) = m;
        *(coordsTest+m) = m*props;
    }
    num = 0;
    for (m=0; m<6; m++){
        for (n=0; n<6; n++){
            DDSDDE[num] = num;
            DDSDDETest[num] = num*props;
            num++;
        } 
    }   

    CORRESPONDENCE::UMATINTTEST(sigmaNP1Loc,statev,DDSDDE,&SSE,&SPD,&SCD,&RPL,
    DDSDDT, DRPLDE,&DRPLDT,strainLoc,depsLoc,timeArray,&dtime,&temp,&dtemp,
    &PREDEF,&DPRED,matnameArray,&nnormal,&nshr,&nstresscomp,&nstatev,&props,
    &nprops,coords,drot,&PNEWDT,&CELENT,defGradN,defGradNP1,
    &NOEL,&NPT,&KSLAY,&KSPT,&JSTEP,&KINC,&nname); 

    for (m=0; m<3; m++){
        TEST_FLOATING_EQUALITY(*(coords+m),*(coordsTest+m),tolerance);
    } 
    for (m=0; m<6; m++){
        TEST_FLOATING_EQUALITY(*(sigmaNP1Loc+m),*(sigmaNP1LocTest+m),tolerance);
        TEST_FLOATING_EQUALITY(*(statev+m),*(statevTest+m),tolerance);
        TEST_FLOATING_EQUALITY(*(strainLoc+m),*(strainLocTest+m),tolerance);
        TEST_FLOATING_EQUALITY(*(depsLoc+m),*(depsLocTest+m),tolerance);
        TEST_FLOATING_EQUALITY(DDSDDT[m],DDSDDTTest[m],tolerance);
        TEST_FLOATING_EQUALITY(DRPLDE[m],DRPLDETest[m],tolerance);
    } 
    for (m=0; m<9; m++){
        TEST_FLOATING_EQUALITY(*(drot+m),*(drotTest+m),tolerance);
        TEST_FLOATING_EQUALITY(*(defGradN+m),*(defGradNTest+m),tolerance);
        TEST_FLOATING_EQUALITY(*(defGradNP1+m),*(defGradNP1Test+m),tolerance);
    } 
    num = 0;
    for (m=0; m<6; m++){
        for (n=0; n<6; n++){
            TEST_FLOATING_EQUALITY(DDSDDE[num],DDSDDETest[num],tolerance);
            num++;
        }
    } 
}

TEUCHOS_UNIT_TEST(correspondence, userMaterialInterface) {
    const double tolerance = 1.0e-15;
    int nnodes = 2;
    const int nstatev = 9;
    std::vector<double> sigmaNLocVector(nnodes*9), sigmaNP1LocVector(nnodes*9), strainLocNVector(nnodes*9), strainNLocVector(nnodes*9), strainNP1LocVector(nnodes*9), statevVector(nnodes*nstatev), coordsVector(nnodes*3), anglesVector(nnodes*3), defGradNVector(nnodes*9), defGradNP1Vector(nnodes*9), rotNVector(nnodes*9), rotNP1Vector(nnodes*9);
    std::vector<double> sigmaNP1LocTestVector(nnodes*9), strainLocTestVector(nnodes*9), statevTestVector(nnodes*nstatev);
    double* sigmaNLoc = &sigmaNLocVector[0];
    double* sigmaNP1Loc = &sigmaNP1LocVector[0];
    double* sigmaNP1LocTest = &sigmaNP1LocTestVector[0];
    double* strainNLoc= &strainNLocVector[0];
    double* strainNP1Loc= &strainNP1LocVector[0];
    double* strainLocTest = &strainLocTestVector[0];
    double* statev = &statevVector[0];
    double* statevTest = &statevTestVector[0];
    double* coords = &coordsVector[0];
    double* rotN = &rotNVector[0];
    double* rotNP1 = &rotNP1Vector[0];
    double* angles = &anglesVector[0];
    double* defGradN = &defGradNVector[0];
    double* defGradNP1 = &defGradNP1Vector[0];


    const double time = 0.1, dtime = 0.1;
    std::vector<double> tempVector(nnodes), dtempVector(nnodes);
    double* temp = &tempVector[0];
    double* dtemp = &dtempVector[0];
    temp[0] = 200.0; temp[1] = 200.0;
    dtemp[0] = 200.0; dtemp[1] = 200.0;
    const std::string matname = "Elastic";
    
    const int nprops = 1;
    const double propsValue = 3;
    std::vector<double> propsVector(nprops);
    double* props = &propsVector[0];
    props[0] = propsValue;

    int m;
    for (m=0; m<nnodes*9; m++){
        *(defGradN+m) = 0;
        *(defGradNP1+m) = m+1;
        *(sigmaNLoc+m) = m;
        *(sigmaNP1Loc+m) = 0;
        *(sigmaNP1LocTest+m) = m*propsValue;

        *(strainNLoc+m) = 0;
        *(strainNP1Loc+m) = 0;
    } 
    for (m=0; m<nnodes*nstatev; m++){
        *(statev+m) = m;
        *(statevTest+m) = m*propsValue;
    }
            
    for (m=0; m<nnodes*3; m++){
        *(angles+m) = m;
    }
  for (m=0; m<nnodes; m++){
    CORRESPONDENCE::computeGreenLagrangeStrain(&defGradNP1[9*m],&strainLocTest[9*m]);
  }
  CORRESPONDENCE::userMaterialInterface(coords,
                                        defGradN, 
                                        defGradNP1, 
                                        strainNLoc,
                                        strainNP1Loc, //
                                        sigmaNLoc,
                                        sigmaNP1Loc,//
                                        nnodes,
                                        nstatev,
                                        statev,//
                                        nprops,
                                        props,
                                        angles,
                                        time,
                                        dtime,
                                        temp,
                                        dtemp,
                                        rotN,
                                        rotNP1,
                                        false,
                                        false,
                                        matname,
                                        true);
   

//    for (m=0; m<nnodes*9; m++){
//        TEST_FLOATING_EQUALITY(*(sigmaNP1Loc+m),*(sigmaNP1LocTest+m),tolerance);
//    }
    for (m=0; m<nnodes*9; m++){
        TEST_FLOATING_EQUALITY(*(strainNP1Loc+m),*(strainLocTest+m),tolerance);
    } 
    for (m=0; m<nnodes*nstatev; m++){
        TEST_FLOATING_EQUALITY(*(statev+m),*(statevTest+m),tolerance);
    } 
        

}

int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
