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
    
    int check = APPROXIMATION::get_field_size(1,A,n_cont,true);
    TEST_EQUALITY(check,96); 
    check = APPROXIMATION::get_field_size(2,A,n_cont,true);
    TEST_EQUALITY(check,208);
    check = APPROXIMATION::get_field_size(2,A,n_cont+5,true);
    TEST_EQUALITY(check,1053); 
    A[6] = 0;
    check = APPROXIMATION::get_field_size(2,A,n_cont,true);
    TEST_EQUALITY(check,112); 
    int nlist[50] = {24, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,24, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
 
    check = APPROXIMATION::get_field_size(2, nlist, n_cont-1,true);
    TEST_EQUALITY(check,450);
    check = APPROXIMATION::get_field_size(2, nlist, n_cont-1,false);
    TEST_EQUALITY(check,1350);
    check = APPROXIMATION::get_field_size(2,A,n_cont,false);
    TEST_EQUALITY(check,448); 
}

TEUCHOS_UNIT_TEST(approximation, get_sample_weighted) {
    const double tolerance = 2.0e-8;
    TEST_FLOATING_EQUALITY(APPROXIMATION::get_sample_weighted(0,0,1),0,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::get_sample_weighted(1,0,2),0.5,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::get_sample_weighted(0,-2,0),1,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::get_sample_weighted(-1,-1,1),0,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::get_sample_weighted(0,-1,1),0.5,tolerance);
    
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

    p = 2; 
    n = 3;
    std::vector<double> kVector3(p+n+1);
    double* k3 = &kVector3[0];
    APPROXIMATION::knots(n, p, true, k3);
    // caclulated with Python
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(0,p,k3,0.0),0.0,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(1,p,k3,0.0),0.0,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(2,p,k3,0.0),0.0,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(0,p,k3,0.5),0.25,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(1,p,k3,0.5),0.5,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::basis_func(2,p,k3,0.5),0.25,tolerance);


}

TEUCHOS_UNIT_TEST(approximation, deriv_basis_func) {
    const double tolerance = 2.0e-7;
    
    int p = 2, n = 4;
    std::vector<double> kVector(p+n+1);
    double* k = &kVector[0];
    APPROXIMATION::knots(n, p, true, k);
    // caclulated with Python
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(0,p,k,0.2),-2.4,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(1,p,k,0.2),1.599999999999,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(2,p,k,0.2),0.8,tolerance);
    
    p = 4;
    n = 7;
    std::vector<double> kVector2(p+n+1);
    double* k2 = &kVector2[0];
    APPROXIMATION::knots(n, p, true, k2);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(0,p,k2,0.1),-4.115999999999999,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(1,p,k2,0.2),-2.58,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(3,p,k2,0.2),1.2240000000000002,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(1,p,k2,0.4),-0.7679999999999996,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(1,p,k2,0.9),0.0,tolerance);

    p = 2;
    n = 3;
    std::vector<double> kVector3(p+n+1);
    double* k3 = &kVector3[0];
    APPROXIMATION::knots(n, p, true, k3);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(0,p,k3,0.5),-1.0,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(1,p,k3,0.5),0,tolerance);
    TEST_FLOATING_EQUALITY(APPROXIMATION::deriv_basis_func(2,p,k3,0.5),1.0,tolerance);
    
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
TEUCHOS_UNIT_TEST(approximation, create_approximation) {
    /*  
    create_approximation
        const int node,
        const int nneighbors,
        const int* neighborhoodlist,
        const double* coordinates,
        const int num_control_points,
        const int degree,
        const bool twoD,
        double* AMatrix
    */

    const double tolerance = 4.0e-8;
    
    int nPoints = 25;
    int numNode = 12;
    int p = 2;
    // lengths are hard coded to avoid warnings
    int nlist[24] = {0,  1, 2, 3, 4, 5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24};
    int ncont = 3;
    
    std::vector<double> AVector(9*25);
    double* A = &AVector[0];
    double ATest[9*25];
    double coor[3*25] =   {-1,  -1,0,
                          -0.5,-1,0,
                          0,   -1,0,
                          0.5, -1,0,
                          1,   -1,0,
                          -1,  -0.5,0,
                          -0.5,-0.5,0,
                          0,   -0.5,0,
                          0.5, -0.5,0,
                          1,   -0.5,0,
                          -1,  0,0,
                          -0.5,0,0,
                          0,0,0,
                          0.5, 0,0,
                          1,   0,0,
                          -1,  0.5,0,
                          -0.5,0.5,0,
                          0,   0.5,0,
                          0.5, 0.5,0,
                          1,   0.5,0,
                          -1,  1,0,
                          -0.5,1,0,
                          0,   1,0,
                          0.5, 1,0,
                          1,   1,0};
    APPROXIMATION::create_approximation(numNode,nPoints-1,nlist,coor,ncont,p,true,A);
    // calculated with Python
    ATest[0]= 0.562499999999996 ;
    ATest[1]= 0.0 ;
    ATest[2]= 0.0 ;
    ATest[3]= 0.0 ;
    ATest[4]= 0.0 ;
    ATest[5]= 0.0 ;
    ATest[6]= 0.0 ;
    ATest[7]= 5.062499999999993 ;
    ATest[8]= -1.687500000000004 ;
    ATest[9]= -2.812500000000002 ;
    ATest[10]= 1.6874999999999971 ;
    ATest[11]= 0.0 ;
    ATest[12]= -1.6875000000000073 ;
    ATest[13]= 0.9374999999999978 ;
    ATest[14]= -0.5625000000000022 ;
    ATest[15]= 0.0 ;
    ATest[16]= -2.812500000000005 ;
    ATest[17]= 0.9374999999999952 ;
    ATest[18]= 1.5624999999999976 ;
    ATest[19]= -0.9375000000000004 ;
    ATest[20]= 0.0 ;
    ATest[21]= 1.6875000000000013 ;
    ATest[22]= -0.5625000000000053 ;
    ATest[23]= -0.9375000000000051 ;
    ATest[24]= 0.5625000000000024 ;
    ATest[25]= -1.162499999999995 ;
    ATest[26]= 0.0 ;
    ATest[27]= 0.0 ;
    ATest[28]= 0.0 ;
    ATest[29]= 0.0 ;
    ATest[30]= 0.0 ;
    ATest[31]= 0.0 ;
    ATest[32]= -1.9124999999999934 ;
    ATest[33]= 0.6375000000000043 ;
    ATest[34]= 1.0625000000000022 ;
    ATest[35]= -0.6374999999999988 ;
    ATest[36]= 0.0 ;
    ATest[37]= 3.4875000000000083 ;
    ATest[38]= -1.9374999999999978 ;
    ATest[39]= 1.1624999999999996 ;
    ATest[40]= 0.0 ;
    ATest[41]= 3.2625000000000037 ;
    ATest[42]= -1.087499999999995 ;
    ATest[43]= -1.8124999999999967 ;
    ATest[44]= 1.0874999999999986 ;
    ATest[45]= 0.0 ;
    ATest[46]= -2.5875000000000017 ;
    ATest[47]= 0.862500000000007 ;
    ATest[48]= 1.4375000000000067 ;
    ATest[49]= -0.862500000000002 ;
    ATest[50]= 0.11249999999999957 ;
    ATest[51]= 0.0 ;
    ATest[52]= 0.0 ;
    ATest[53]= 0.0 ;
    ATest[54]= 0.0 ;
    ATest[55]= 0.0 ;
    ATest[56]= 0.0 ;
    ATest[57]= 0.11249999999999913 ;
    ATest[58]= -0.03750000000000018 ;
    ATest[59]= -0.06250000000000017 ;
    ATest[60]= 0.037499999999999395 ;
    ATest[61]= 0.0 ;
    ATest[62]= -0.3375000000000009 ;
    ATest[63]= 0.18749999999999992 ;
    ATest[64]= -0.11249999999999993 ;
    ATest[65]= 0.0 ;
    ATest[66]= 0.33749999999999963 ;
    ATest[67]= -0.11250000000000077 ;
    ATest[68]= -0.18750000000000058 ;
    ATest[69]= 0.11250000000000043 ;
    ATest[70]= 0.0 ;
    ATest[71]= 2.137500000000001 ;
    ATest[72]= -0.7125000000000012 ;
    ATest[73]= -1.1875000000000013 ;
    ATest[74]= 0.7125000000000004 ;
    ATest[75]= -1.1624999999999968 ;
    ATest[76]= 0.0 ;
    ATest[77]= 0.0 ;
    ATest[78]= 0.0 ;
    ATest[79]= 0.0 ;
    ATest[80]= 0.0 ;
    ATest[81]= 0.0 ;
    ATest[82]= -1.912499999999993 ;
    ATest[83]= 3.4875000000000056 ;
    ATest[84]= 3.2625000000000046 ;
    ATest[85]= -2.587499999999996 ;
    ATest[86]= 0.0 ;
    ATest[87]= 0.6375000000000062 ;
    ATest[88]= -1.0874999999999977 ;
    ATest[89]= 0.8625000000000029 ;
    ATest[90]= 0.0 ;
    ATest[91]= 1.0625000000000022 ;
    ATest[92]= -1.9374999999999984 ;
    ATest[93]= -1.8124999999999987 ;
    ATest[94]= 1.4375000000000009 ;
    ATest[95]= 0.0 ;
    ATest[96]= -0.6375000000000041 ;
    ATest[97]= 1.1625000000000016 ;
    ATest[98]= 1.0875000000000024 ;
    ATest[99]= -0.8625000000000024 ;
    ATest[100]= 2.402499999999998 ;
    ATest[101]= 0.0 ;
    ATest[102]= 0.0 ;
    ATest[103]= 0.0 ;
    ATest[104]= 0.0 ;
    ATest[105]= 0.0 ;
    ATest[106]= 0.0 ;
    ATest[107]= 0.7224999999999951 ;
    ATest[108]= -1.3175000000000028 ;
    ATest[109]= -1.2325000000000017 ;
    ATest[110]= 0.9774999999999985 ;
    ATest[111]= 0.0 ;
    ATest[112]= -1.3175000000000037 ;
    ATest[113]= 2.2475 ;
    ATest[114]= -1.7824999999999998 ;
    ATest[115]= 0.0 ;
    ATest[116]= -1.2325000000000017 ;
    ATest[117]= 2.247499999999998 ;
    ATest[118]= 2.102499999999999 ;
    ATest[119]= -1.6674999999999986 ;
    ATest[120]= 0.0 ;
    ATest[121]= 0.977500000000004 ;
    ATest[122]= -1.7825000000000029 ;
    ATest[123]= -1.6675000000000035 ;
    ATest[124]= 1.3225000000000018 ;
    ATest[125]= -0.23249999999999982 ;
    ATest[126]= 0.0 ;
    ATest[127]= 0.0 ;
    ATest[128]= 0.0 ;
    ATest[129]= 0.0 ;
    ATest[130]= 0.0 ;
    ATest[131]= 0.0 ;
    ATest[132]= -0.04249999999999919 ;
    ATest[133]= 0.0775000000000002 ;
    ATest[134]= 0.07250000000000018 ;
    ATest[135]= -0.057499999999999385 ;
    ATest[136]= 0.0 ;
    ATest[137]= 0.1275000000000005 ;
    ATest[138]= -0.21749999999999992 ;
    ATest[139]= 0.17249999999999988 ;
    ATest[140]= 0.0 ;
    ATest[141]= -0.12750000000000006 ;
    ATest[142]= 0.2325000000000006 ;
    ATest[143]= 0.2175000000000003 ;
    ATest[144]= -0.17250000000000043 ;
    ATest[145]= 0.0 ;
    ATest[146]= -0.807500000000001 ;
    ATest[147]= 1.4725000000000008 ;
    ATest[148]= 1.3775000000000008 ;
    ATest[149]= -1.0925000000000002 ;
    ATest[150]= 0.11249999999999957 ;
    ATest[151]= 0.0 ;
    ATest[152]= 0.0 ;
    ATest[153]= 0.0 ;
    ATest[154]= 0.0 ;
    ATest[155]= 0.0 ;
    ATest[156]= 0.0 ;
    ATest[157]= 0.11249999999999909 ;
    ATest[158]= -0.33750000000000013 ;
    ATest[159]= 0.3374999999999996 ;
    ATest[160]= 2.137499999999999 ;
    ATest[161]= 0.0 ;
    ATest[162]= -0.03750000000000062 ;
    ATest[163]= -0.11250000000000045 ;
    ATest[164]= -0.7125000000000005 ;
    ATest[165]= 0.0 ;
    ATest[166]= -0.06250000000000035 ;
    ATest[167]= 0.18749999999999958 ;
    ATest[168]= -0.18750000000000036 ;
    ATest[169]= -1.1875 ;
    ATest[170]= 0.0 ;
    ATest[171]= 0.03749999999999993 ;
    ATest[172]= -0.11250000000000018 ;
    ATest[173]= 0.11249999999999985 ;
    ATest[174]= 0.7125 ;
    ATest[175]= -0.2324999999999994 ;
    ATest[176]= 0.0 ;
    ATest[177]= 0.0 ;
    ATest[178]= 0.0 ;
    ATest[179]= 0.0 ;
    ATest[180]= 0.0 ;
    ATest[181]= 0.0 ;
    ATest[182]= -0.04249999999999941 ;
    ATest[183]= 0.12750000000000006 ;
    ATest[184]= -0.1274999999999999 ;
    ATest[185]= -0.8074999999999994 ;
    ATest[186]= 0.0 ;
    ATest[187]= 0.07750000000000068 ;
    ATest[188]= 0.23250000000000046 ;
    ATest[189]= 1.4725000000000001 ;
    ATest[190]= 0.0 ;
    ATest[191]= 0.07250000000000038 ;
    ATest[192]= -0.21749999999999947 ;
    ATest[193]= 0.21750000000000033 ;
    ATest[194]= 1.3774999999999997 ;
    ATest[195]= 0.0 ;
    ATest[196]= -0.057499999999999885 ;
    ATest[197]= 0.1725000000000002 ;
    ATest[198]= -0.17249999999999974 ;
    ATest[199]= -1.0924999999999998 ;
    ATest[200]= 0.022499999999999964 ;
    ATest[201]= 0.0 ;
    ATest[202]= 0.0 ;
    ATest[203]= 0.0 ;
    ATest[204]= 0.0 ;
    ATest[205]= 0.0 ;
    ATest[206]= 0.0 ;
    ATest[207]= 0.002499999999999919 ;
    ATest[208]= -0.007499999999999979 ;
    ATest[209]= 0.007500000000000062 ;
    ATest[210]= 0.04749999999999988 ;
    ATest[211]= 0.0 ;
    ATest[212]= -0.007500000000000007 ;
    ATest[213]= -0.022500000000000075 ;
    ATest[214]= -0.1425 ;
    ATest[215]= 0.0 ;
    ATest[216]= 0.007500000000000007 ;
    ATest[217]= -0.02250000000000016 ;
    ATest[218]= 0.02249999999999991 ;
    ATest[219]= 0.14250000000000007 ;
    ATest[220]= 0.0 ;
    ATest[221]= 0.04750000000000004 ;
    ATest[222]= -0.1425000000000001 ;
    ATest[223]= 0.14249999999999996 ;
    ATest[224]= 0.9025 ;
    for(int i=0 ; i<ncont*ncont ; ++i){       
        for(int j=0 ; j<nPoints; ++j){ 
            TEST_FLOATING_EQUALITY(float(i*nPoints + j),float(i*nPoints + j),tolerance); 
            TEST_FLOATING_EQUALITY(A[i*nPoints + j],ATest[i*nPoints + j],tolerance);  
        }
    }

   
}
TEUCHOS_UNIT_TEST(approximation, get_approximation) {
    /*  
    create_approximation
        const int node,
        const int nneighbors,
        const int* neighborhoodlist,
        const double* coordinates,
        const int num_control_points,
        const int degree,
        const bool twoD,
        double* AMatrix
    */

    const double tolerance = 4.0e-8;
    
    int nPoints = 25;
    int p = 2;
    // lengths are hard coded to avoid warnings
    int nlist[50] = {24, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,24, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
    int ncont = 3;
    int len = APPROXIMATION::get_field_size(2, nlist, nPoints, true);
    
    std::vector<double> AVector(len);
    double* A = &AVector[0];
    
    double ATest[9*50];
    double coor[3*26] =   {0,0,0,
                           0,0,0,
                          -1,  -1,0,
                          -0.5,-1,0,
                          0,   -1,0,
                          0.5, -1,0,
                          1,   -1,0,
                          -1,  -0.5,0,
                          -0.5,-0.5,0,
                          0,   -0.5,0,
                          0.5, -0.5,0,
                          1,   -0.5,0,
                          -1,  0,0,
                          -0.5,0,0,
                          0.5, 0,0,
                          1,   0,0,
                          -1,  0.5,0,
                          -0.5,0.5,0,
                          0,   0.5,0,
                          0.5, 0.5,0,
                          1,   0.5,0,
                          -1,  1,0,
                          -0.5,1,0,
                          0,   1,0,
                          0.5, 1,0,
                          1,   1,0};
    //APPROXIMATION::create_approximation(numNode,nPoints,nlist,coor,ncont,p,true,A);
    APPROXIMATION::get_approximation(2, nlist,coor,ncont,p,true,A);

    // calculated with Python
    ATest[0]= 0.562499999999996 ;
    ATest[1]= 0.0 ;
    ATest[2]= 0.0 ;
    ATest[3]= 0.0 ;
    ATest[4]= 0.0 ;
    ATest[5]= 0.0 ;
    ATest[6]= 0.0 ;
    ATest[7]= 5.062499999999993 ;
    ATest[8]= -1.687500000000004 ;
    ATest[9]= -2.812500000000002 ;
    ATest[10]= 1.6874999999999971 ;
    ATest[11]= 0.0 ;
    ATest[12]= -1.6875000000000073 ;
    ATest[13]= 0.9374999999999978 ;
    ATest[14]= -0.5625000000000022 ;
    ATest[15]= 0.0 ;
    ATest[16]= -2.812500000000005 ;
    ATest[17]= 0.9374999999999952 ;
    ATest[18]= 1.5624999999999976 ;
    ATest[19]= -0.9375000000000004 ;
    ATest[20]= 0.0 ;
    ATest[21]= 1.6875000000000013 ;
    ATest[22]= -0.5625000000000053 ;
    ATest[23]= -0.9375000000000051 ;
    ATest[24]= 0.5625000000000024 ;
    ATest[25]= -1.162499999999995 ;
    ATest[26]= 0.0 ;
    ATest[27]= 0.0 ;
    ATest[28]= 0.0 ;
    ATest[29]= 0.0 ;
    ATest[30]= 0.0 ;
    ATest[31]= 0.0 ;
    ATest[32]= -1.9124999999999934 ;
    ATest[33]= 0.6375000000000043 ;
    ATest[34]= 1.0625000000000022 ;
    ATest[35]= -0.6374999999999988 ;
    ATest[36]= 0.0 ;
    ATest[37]= 3.4875000000000083 ;
    ATest[38]= -1.9374999999999978 ;
    ATest[39]= 1.1624999999999996 ;
    ATest[40]= 0.0 ;
    ATest[41]= 3.2625000000000037 ;
    ATest[42]= -1.087499999999995 ;
    ATest[43]= -1.8124999999999967 ;
    ATest[44]= 1.0874999999999986 ;
    ATest[45]= 0.0 ;
    ATest[46]= -2.5875000000000017 ;
    ATest[47]= 0.862500000000007 ;
    ATest[48]= 1.4375000000000067 ;
    ATest[49]= -0.862500000000002 ;
    ATest[50]= 0.11249999999999957 ;
    ATest[51]= 0.0 ;
    ATest[52]= 0.0 ;
    ATest[53]= 0.0 ;
    ATest[54]= 0.0 ;
    ATest[55]= 0.0 ;
    ATest[56]= 0.0 ;
    ATest[57]= 0.11249999999999913 ;
    ATest[58]= -0.03750000000000018 ;
    ATest[59]= -0.06250000000000017 ;
    ATest[60]= 0.037499999999999395 ;
    ATest[61]= 0.0 ;
    ATest[62]= -0.3375000000000009 ;
    ATest[63]= 0.18749999999999992 ;
    ATest[64]= -0.11249999999999993 ;
    ATest[65]= 0.0 ;
    ATest[66]= 0.33749999999999963 ;
    ATest[67]= -0.11250000000000077 ;
    ATest[68]= -0.18750000000000058 ;
    ATest[69]= 0.11250000000000043 ;
    ATest[70]= 0.0 ;
    ATest[71]= 2.137500000000001 ;
    ATest[72]= -0.7125000000000012 ;
    ATest[73]= -1.1875000000000013 ;
    ATest[74]= 0.7125000000000004 ;
    ATest[75]= -1.1624999999999968 ;
    ATest[76]= 0.0 ;
    ATest[77]= 0.0 ;
    ATest[78]= 0.0 ;
    ATest[79]= 0.0 ;
    ATest[80]= 0.0 ;
    ATest[81]= 0.0 ;
    ATest[82]= -1.912499999999993 ;
    ATest[83]= 3.4875000000000056 ;
    ATest[84]= 3.2625000000000046 ;
    ATest[85]= -2.587499999999996 ;
    ATest[86]= 0.0 ;
    ATest[87]= 0.6375000000000062 ;
    ATest[88]= -1.0874999999999977 ;
    ATest[89]= 0.8625000000000029 ;
    ATest[90]= 0.0 ;
    ATest[91]= 1.0625000000000022 ;
    ATest[92]= -1.9374999999999984 ;
    ATest[93]= -1.8124999999999987 ;
    ATest[94]= 1.4375000000000009 ;
    ATest[95]= 0.0 ;
    ATest[96]= -0.6375000000000041 ;
    ATest[97]= 1.1625000000000016 ;
    ATest[98]= 1.0875000000000024 ;
    ATest[99]= -0.8625000000000024 ;
    ATest[100]= 2.402499999999998 ;
    ATest[101]= 0.0 ;
    ATest[102]= 0.0 ;
    ATest[103]= 0.0 ;
    ATest[104]= 0.0 ;
    ATest[105]= 0.0 ;
    ATest[106]= 0.0 ;
    ATest[107]= 0.7224999999999951 ;
    ATest[108]= -1.3175000000000028 ;
    ATest[109]= -1.2325000000000017 ;
    ATest[110]= 0.9774999999999985 ;
    ATest[111]= 0.0 ;
    ATest[112]= -1.3175000000000037 ;
    ATest[113]= 2.2475 ;
    ATest[114]= -1.7824999999999998 ;
    ATest[115]= 0.0 ;
    ATest[116]= -1.2325000000000017 ;
    ATest[117]= 2.247499999999998 ;
    ATest[118]= 2.102499999999999 ;
    ATest[119]= -1.6674999999999986 ;
    ATest[120]= 0.0 ;
    ATest[121]= 0.977500000000004 ;
    ATest[122]= -1.7825000000000029 ;
    ATest[123]= -1.6675000000000035 ;
    ATest[124]= 1.3225000000000018 ;
    ATest[125]= -0.23249999999999982 ;
    ATest[126]= 0.0 ;
    ATest[127]= 0.0 ;
    ATest[128]= 0.0 ;
    ATest[129]= 0.0 ;
    ATest[130]= 0.0 ;
    ATest[131]= 0.0 ;
    ATest[132]= -0.04249999999999919 ;
    ATest[133]= 0.0775000000000002 ;
    ATest[134]= 0.07250000000000018 ;
    ATest[135]= -0.057499999999999385 ;
    ATest[136]= 0.0 ;
    ATest[137]= 0.1275000000000005 ;
    ATest[138]= -0.21749999999999992 ;
    ATest[139]= 0.17249999999999988 ;
    ATest[140]= 0.0 ;
    ATest[141]= -0.12750000000000006 ;
    ATest[142]= 0.2325000000000006 ;
    ATest[143]= 0.2175000000000003 ;
    ATest[144]= -0.17250000000000043 ;
    ATest[145]= 0.0 ;
    ATest[146]= -0.807500000000001 ;
    ATest[147]= 1.4725000000000008 ;
    ATest[148]= 1.3775000000000008 ;
    ATest[149]= -1.0925000000000002 ;
    ATest[150]= 0.11249999999999957 ;
    ATest[151]= 0.0 ;
    ATest[152]= 0.0 ;
    ATest[153]= 0.0 ;
    ATest[154]= 0.0 ;
    ATest[155]= 0.0 ;
    ATest[156]= 0.0 ;
    ATest[157]= 0.11249999999999909 ;
    ATest[158]= -0.33750000000000013 ;
    ATest[159]= 0.3374999999999996 ;
    ATest[160]= 2.137499999999999 ;
    ATest[161]= 0.0 ;
    ATest[162]= -0.03750000000000062 ;
    ATest[163]= -0.11250000000000045 ;
    ATest[164]= -0.7125000000000005 ;
    ATest[165]= 0.0 ;
    ATest[166]= -0.06250000000000035 ;
    ATest[167]= 0.18749999999999958 ;
    ATest[168]= -0.18750000000000036 ;
    ATest[169]= -1.1875 ;
    ATest[170]= 0.0 ;
    ATest[171]= 0.03749999999999993 ;
    ATest[172]= -0.11250000000000018 ;
    ATest[173]= 0.11249999999999985 ;
    ATest[174]= 0.7125 ;
    ATest[175]= -0.2324999999999994 ;
    ATest[176]= 0.0 ;
    ATest[177]= 0.0 ;
    ATest[178]= 0.0 ;
    ATest[179]= 0.0 ;
    ATest[180]= 0.0 ;
    ATest[181]= 0.0 ;
    ATest[182]= -0.04249999999999941 ;
    ATest[183]= 0.12750000000000006 ;
    ATest[184]= -0.1274999999999999 ;
    ATest[185]= -0.8074999999999994 ;
    ATest[186]= 0.0 ;
    ATest[187]= 0.07750000000000068 ;
    ATest[188]= 0.23250000000000046 ;
    ATest[189]= 1.4725000000000001 ;
    ATest[190]= 0.0 ;
    ATest[191]= 0.07250000000000038 ;
    ATest[192]= -0.21749999999999947 ;
    ATest[193]= 0.21750000000000033 ;
    ATest[194]= 1.3774999999999997 ;
    ATest[195]= 0.0 ;
    ATest[196]= -0.057499999999999885 ;
    ATest[197]= 0.1725000000000002 ;
    ATest[198]= -0.17249999999999974 ;
    ATest[199]= -1.0924999999999998 ;
    ATest[200]= 0.022499999999999964 ;
    ATest[201]= 0.0 ;
    ATest[202]= 0.0 ;
    ATest[203]= 0.0 ;
    ATest[204]= 0.0 ;
    ATest[205]= 0.0 ;
    ATest[206]= 0.0 ;
    ATest[207]= 0.002499999999999919 ;
    ATest[208]= -0.007499999999999979 ;
    ATest[209]= 0.007500000000000062 ;
    ATest[210]= 0.04749999999999988 ;
    ATest[211]= 0.0 ;
    ATest[212]= -0.007500000000000007 ;
    ATest[213]= -0.022500000000000075 ;
    ATest[214]= -0.1425 ;
    ATest[215]= 0.0 ;
    ATest[216]= 0.007500000000000007 ;
    ATest[217]= -0.02250000000000016 ;
    ATest[218]= 0.02249999999999991 ;
    ATest[219]= 0.14250000000000007 ;
    ATest[220]= 0.0 ;
    ATest[221]= 0.04750000000000004 ;
    ATest[222]= -0.1425000000000001 ;
    ATest[223]= 0.14249999999999996 ;
    ATest[224]= 0.9025 ;
    for(int iID=0 ; iID<2 ; ++iID){ 
    
        for(int i=0 ; i<ncont*ncont ; ++i){       
            for(int j=0 ; j<nPoints; ++j){ 
                TEST_FLOATING_EQUALITY(float(i*nPoints + j),float(i*nPoints + j),tolerance); 
                TEST_FLOATING_EQUALITY(A[i*nPoints + j + iID*(nPoints*ncont*ncont)],ATest[i*nPoints + j],tolerance);  
            }
        }
    }
   
}
TEUCHOS_UNIT_TEST(approximation, get_control_point) {
   /*  
    create_approximation
        const int node,
        const int nneighbors,
        const int* neighborhoodlist,
        const double* coordinates,
        const int num_control_points,
        const int degree,
        const bool twoD,
        double* AMatrix
    */

    const double tolerance = 4.0e-8;
    
    int nPoints = 24;
    int numNode = 0;
    int p = 2;
    // lengths are hard coded to avoid warnings
    int nlist[24] = {1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
    int ncont = 3;
    double *contP;
    contP = new double[ncont * ncont * PeridigmNS::dof()];
    std::vector<double> AVector(9*25);
    double* A = &AVector[0];
    double contPTest[27];
    double coor[3*25] =   {0,0,0,
                          -1,  -1,0,
                          -0.5,-1,0,
                          0,   -1,0,
                          0.5, -1,0,
                          1,   -1,0,
                          -1,  -0.5,0,
                          -0.5,-0.5,0,
                          0,   -0.5,0,
                          0.5, -0.5,0,
                          1,   -0.5,0,
                          -1,  0,0,
                          -0.5,0,0,
                          0.5, 0,0,
                          1,   0,0,
                          -1,  0.5,0,
                          -0.5,0.5,0,
                          0,   0.5,0,
                          0.5, 0.5,0,
                          1,   0.5,0,
                          -1,  1,0,
                          -0.5,1,0,
                          0,   1,0,
                          0.5, 1,0,
                          1,   1,0};
    contPTest[0]= -0.9999999999999998 ;
    contPTest[1]= -1.000000000000005 ;
    contPTest[2]= 0 ;
    contPTest[3]= -1.0000000000000036 ;
    contPTest[4]= 0.0 ;
    contPTest[5]= 0 ;
    contPTest[6]= -1.0000000000000002 ;
    contPTest[7]= 0.9999999999999991 ;
    contPTest[8]= 0 ;
    contPTest[9]= 0.0 ;
    contPTest[10]= -1.0000000000000102 ;
    contPTest[11]= 0 ;
    contPTest[12]= 0.0 ;
    contPTest[13]= 0 ;
    contPTest[14]= 0 ;
    contPTest[15]= 0 ;
    contPTest[16]= 1.0 ;
    contPTest[17]= 0 ;
    contPTest[18]= 0.9999999999999987 ;
    contPTest[19]= -0.9999999999999999 ;
    contPTest[20]= 0 ;
    contPTest[21]= 1.0000000000000004 ;
    contPTest[22]= 0 ;
    contPTest[23]= 0 ;
    contPTest[24]= 0.9999999999999999 ;
    contPTest[25]= 0.9999999999999999 ;
    contPTest[26]= 0 ;
    APPROXIMATION::create_approximation(numNode,nPoints,nlist,coor,ncont,p,true,A);
    APPROXIMATION::get_control_point(numNode,nPoints,nlist,ncont,coor,A,true,contP);
    // calculated with Python
    
    for(int i=0 ; i<ncont*ncont*PeridigmNS::dof() ; ++i){       
            TEST_ASSERT(abs(contP[i]-contPTest[i])<tolerance);  
        }
}




TEUCHOS_UNIT_TEST(approximation, get_gradient) {
   /*  
    create_approximation
        const int node,
        const int nneighbors,
        const int* neighborhoodlist,
        const double* coordinates,
        const int num_control_points,
        const int degree,
        const bool twoD,
        double* AMatrix
    */

    const double tolerance = 4.0e-10;

    int p = 2;
    int ncont = 3;
    double *contP;
    contP = new double[ncont * ncont * PeridigmNS::dof()];
    double *contP3D;
    contP3D = new double[ncont * ncont * ncont * PeridigmNS::dof()];

    double gradTest[9];
    contP[0]= -0.9999999999999998 ;    contP[1]= -1.000000000000005 ;    contP[2]= 0 ;
    contP[3]= -1.0000000000000036 ;    contP[4]= 0.0 ;    contP[5]= 0 ;
    contP[6]= -1.0000000000000002 ;    contP[7]= 0.9999999999999991 ;    contP[8]= 0 ;
    contP[9]= 0.0 ;    contP[10]= -1.0000000000000102 ;    contP[11]= 0 ;
    contP[12]= 0.0 ;    contP[13]= 0 ;    contP[14]= 0 ;
    contP[15]= 0 ;    contP[16]= 1.0 ;    contP[17]= 0 ;
    contP[18]= 0.9999999999999987 ;    contP[19]= -0.9999999999999999 ;    contP[20]= 0 ;
    contP[21]= 1.0000000000000004 ;    contP[22]= 0 ;    contP[23]= 0 ;
    contP[24]= 0.9999999999999999 ;    contP[25]= 0.9999999999999999 ;    contP[26]= 0 ;

    std::vector<double> UVector(ncont + p + 1);
    double* U = &UVector[0];
    std::vector<double> VVector(ncont + p + 1);
    double* V = &VVector[0]; 
    std::vector<double> WVector(ncont + p + 1);
    double* W = &WVector[0]; 
    double u = 0.5, v = 0.5, w = 0.5;
    APPROXIMATION::knots(ncont,p,true,U);
    APPROXIMATION::knots(ncont,p,true,V);
    APPROXIMATION::knots(ncont,p,true,W);

    Eigen::MatrixXd gradientMxM(2,2);

    APPROXIMATION::get_gradient(p,ncont,contP,U,u,V,v,W,w,true,gradientMxM);
    gradTest[0]= 2.0000000000000013 ;
    gradTest[1]= -1.970645868709653e-15 ;
    gradTest[2]= -2.4702462297909733e-15 ;
    gradTest[3]= 2.000000000000006 ;



    for(int i=0 ; i<2 ; ++i){    
        for(int j=0 ; j<2 ; ++j){        
            TEST_ASSERT(abs(gradientMxM(i,j)-gradTest[2*i+j])<tolerance);     
        }
    }
    APPROXIMATION::get_gradient(p,ncont,contP,U,u,V,v,W,w,true,gradientMxM);
    gradTest[0]= 2.0;
    gradTest[1]= 0 ;
    gradTest[2]= 0 ;
    gradTest[3]= 0.0 ;
    gradTest[4]= 2.0 ;
    gradTest[5]= 0 ;
    gradTest[6]= 0 ;
    gradTest[7]= 0 ;
    gradTest[8]= 2.0 ;
    
    

    contP3D[0]= -0.9999999999999998 ;    contP3D[1]= -1.000000000000005 ;    contP3D[2]= -1 ;
    contP3D[3]= -1.0000000000000036 ;    contP3D[4]= 0.0 ;                   contP3D[5]= -1 ;
    contP3D[6]= -1.0000000000000002 ;    contP3D[7]= 0.9999999999999991 ;    contP3D[8]= -1 ;
    contP3D[9]= 0.0 ;                    contP3D[10]= -1.0000000000000102 ;  contP3D[11]= -1 ;
    contP3D[12]= 0.0 ;                   contP3D[13]= 0 ;                    contP3D[14]= -1 ;
    contP3D[15]= 0 ;                     contP3D[16]= 1.0 ;                  contP3D[17]= -1 ;
    contP3D[18]= 0.9999999999999987 ;    contP3D[19]= -0.9999999999999999 ;  contP3D[20]= -1 ;
    contP3D[21]= 1.0000000000000004 ;    contP3D[22]= 0 ;                    contP3D[23]= -1 ;
    contP3D[24]= 0.9999999999999999 ;    contP3D[25]= 0.9999999999999999 ;   contP3D[26]= -1 ; 
    int offset = 27;
    contP3D[0+offset]= -0.9999999999999998 ;    contP3D[1+offset]= -1.000000000000005 ;    contP3D[2+offset] = 0 ;
    contP3D[3+offset]= -1.0000000000000036 ;    contP3D[4+offset]= 0.0 ;                   contP3D[5+offset] = 0 ;
    contP3D[6+offset]= -1.0000000000000002 ;    contP3D[7+offset]= 0.9999999999999991 ;    contP3D[8+offset] = 0 ;
    contP3D[9+offset]= 0.0 ;                    contP3D[10+offset]= -1.0000000000000102 ;  contP3D[11+offset]= 0 ;
    contP3D[12+offset]= 0.0 ;                   contP3D[13+offset]= 0 ;                    contP3D[14+offset]= 0 ;
    contP3D[15+offset]= 0 ;                     contP3D[16+offset]= 1.0 ;                  contP3D[17+offset]= 0 ;
    contP3D[18+offset]= 0.9999999999999987 ;    contP3D[19+offset]= -0.9999999999999999 ;  contP3D[20+offset]= 0 ;
    contP3D[21+offset]= 1.0000000000000004 ;    contP3D[22+offset]= 0 ;                    contP3D[23+offset]= 0 ;
    contP3D[24+offset]= 0.9999999999999999 ;    contP3D[25+offset]= 0.9999999999999999 ;   contP3D[26+offset]= 0 ;  
    offset += 27;    
    contP3D[0+offset]= -0.9999999999999998 ;    contP3D[1+offset]= -1.000000000000005 ;    contP3D[2+offset] = 1 ;
    contP3D[3+offset]= -1.0000000000000036 ;    contP3D[4+offset]= 0.0 ;                   contP3D[5+offset] = 1 ;
    contP3D[6+offset]= -1.0000000000000002 ;    contP3D[7+offset]= 0.9999999999999991 ;    contP3D[8+offset] = 1 ;
    contP3D[9+offset]= 0.0 ;                    contP3D[10+offset]= -1.0000000000000102 ;  contP3D[11+offset]= 1 ;
    contP3D[12+offset]= 0.0 ;                   contP3D[13+offset]= 0 ;                    contP3D[14+offset]= 1 ;
    contP3D[15+offset]= 0 ;                     contP3D[16+offset]= 1.0 ;                  contP3D[17+offset]= 1 ;
    contP3D[18+offset]= 0.9999999999999987 ;    contP3D[19+offset]= -0.9999999999999999 ;  contP3D[20+offset]= 1 ;
    contP3D[21+offset]= 1.0000000000000004 ;    contP3D[22+offset]= 0 ;                    contP3D[23+offset]= 1 ;
    contP3D[24+offset]= 0.9999999999999999 ;    contP3D[25+offset]= 0.9999999999999999 ;   contP3D[26+offset]= 1 ;  
    Eigen::MatrixXd gradient3DMxM(3,3);
    APPROXIMATION::get_gradient(p,ncont,contP3D,U,u,V,v,W,w,false,gradient3DMxM);
    gradTest[0]= 0.0 ;
    gradTest[1]= 0.0 ;
    gradTest[2]= 2.0 ;
    gradTest[3]= 2.0 ;
    gradTest[4]= 0.0 ;
    gradTest[5]= 0.0 ;
    gradTest[6]= 0.0 ;
    gradTest[7]= 2.0 ;
    gradTest[8]= 0.0 ;

    for(int i=0 ; i<3 ; ++i){    
        for(int j=0 ; j<3 ; ++j){        
            TEST_ASSERT(abs(gradient3DMxM(i,j)-gradTest[3*i+j])<tolerance);  
        }
    }

 
}

TEUCHOS_UNIT_TEST(approximation, get_jacobian) {
   /*  
    create_approximation
        const int node,
        const int nneighbors,
        const int* neighborhoodlist,
        const double* coordinates,
        const int num_control_points,
        const int degree,
        const bool twoD,
        double* AMatrix
    */

    const double tolerance = 4.0e-10;
   
    int p = 2;
    // lengths are hard coded to avoid warnings
    int ncont = 3;
    double *jacobian;
    jacobian = new double[2 * PeridigmNS::dof() * PeridigmNS::dof()];
    double jacobianTest[18];
    double *contP;
    contP = new double[ncont * ncont * PeridigmNS::dof()];

    contP[0]= -0.9999999999999998 ;
    contP[1]= -1.000000000000005 ;
    contP[2]= 0 ;
    contP[3]= -1.0000000000000036 ;
    contP[4]= 0.0 ;
    contP[5]= 0 ;
    contP[6]= -1.0000000000000002 ;
    contP[7]= 0.9999999999999991 ;
    contP[8]= 0 ;
    contP[9]= 0.0 ;
    contP[10]= -1.0000000000000102 ;
    contP[11]= 0 ;
    contP[12]= 0.0 ;
    contP[13]= 0 ;
    contP[14]= 0 ;
    contP[15]= 0 ;
    contP[16]= 1.0 ;
    contP[17]= 0 ;
    contP[18]= 0.9999999999999987 ;
    contP[19]= -0.9999999999999999 ;
    contP[20]= 0 ;
    contP[21]= 1.0000000000000004 ;
    contP[22]= 0 ;
    contP[23]= 0 ;
    contP[24]= 0.9999999999999999 ;
    contP[25]= 0.9999999999999999 ;
    contP[26]= 0 ;
    std::vector<double> UVector(ncont + p + 1);
    double* U = &UVector[0];
    std::vector<double> VVector(ncont + p + 1);
    double* V = &VVector[0]; 
    std::vector<double> WVector(ncont + p + 1);
    double* W = &WVector[0]; 
    double u = 0.5, v = 0.5, w = 0.5;
    APPROXIMATION::knots(ncont,p,true,U);
    APPROXIMATION::knots(ncont,p,true,V);
    APPROXIMATION::knots(ncont,p,true,W);

    APPROXIMATION::get_jacobian(p,ncont,contP,U,u,V,v,W,w,true,jacobian);
    jacobianTest[0]= 0.49999999999999967 ;
    jacobianTest[1]= 4.926614671774113e-16 ;
    jacobianTest[2]= 0 ;
    jacobianTest[3]= 6.175615574477411e-16 ;
    jacobianTest[4]= 0.49999999999999845 ;
    jacobianTest[5]= 0 ;
    jacobianTest[6]= 0 ;
    jacobianTest[7]= 0 ;
    jacobianTest[8]= 0 ;

    // calculated with Python
     
    for(int i=0 ; i<PeridigmNS::dof()*PeridigmNS::dof() ; ++i){       
            TEST_ASSERT(abs(jacobian[i]-jacobianTest[i])<tolerance);  
        }
    
}

TEUCHOS_UNIT_TEST(approximation, get_jacobians) {
   /*  
    create_approximation
        const int node,
        const int nneighbors,
        const int* neighborhoodlist,
        const double* coordinates,
        const int num_control_points,
        const int degree,
        const bool twoD,
        double* AMatrix
    */

    const double tolerance = 4.0e-8;

    int p = 2;
    // lengths are hard coded to avoid warnings
    int ncont = 3;
    
    double *contP;
    contP = new double[2 * ncont * ncont * PeridigmNS::dof()];
;
    double *jacobian;
    jacobian = new double[2 * PeridigmNS::dof() * PeridigmNS::dof()];
    double jacobianTest[18];
    

    contP[0]= -0.9999999999999998 ;
    contP[1]= -1.000000000000005 ;
    contP[2]= 0 ;
    contP[3]= -1.0000000000000036 ;
    contP[4]= 0.0 ;
    contP[5]= 0 ;
    contP[6]= -1.0000000000000002 ;
    contP[7]= 0.9999999999999991 ;
    contP[8]= 0 ;
    contP[9]= 0.0 ;
    contP[10]= -1.0000000000000102 ;
    contP[11]= 0 ;
    contP[12]= 0.0 ;
    contP[13]= 0 ;
    contP[14]= 0 ;
    contP[15]= 0 ;
    contP[16]= 1.0 ;
    contP[17]= 0 ;
    contP[18]= 0.9999999999999987 ;
    contP[19]= -0.9999999999999999 ;
    contP[20]= 0 ;
    contP[21]= 1.0000000000000004 ;
    contP[22]= 0 ;
    contP[23]= 0 ;
    contP[24]= 0.9999999999999999 ;
    contP[25]= 0.9999999999999999 ;
    contP[26]= 0 ;
    contP[0+27]= -0.9999999999999998 ;
    contP[1+27]= -1.000000000000005 ;
    contP[2+27]= 0 ;
    contP[3+27]= -1.0000000000000036 ;
    contP[4+27]= 0.0 ;
    contP[5+27]= 0 ;
    contP[6+27]= -1.0000000000000002 ;
    contP[7+27]= 0.9999999999999991 ;
    contP[8+27]= 0 ;
    contP[9+27]= 0.0 ;
    contP[10+27]= -1.0000000000000102 ;
    contP[11+27]= 0 ;
    contP[12+27]= 0.0 ;
    contP[13+27]= 0 ;
    contP[14+27]= 0 ;
    contP[15+27]= 0 ;
    contP[16+27]= 1.0 ;
    contP[17+27]= 0 ;
    contP[18+27]= 0.9999999999999987 ;
    contP[19+27]= -0.9999999999999999 ;
    contP[20+27]= 0 ;
    contP[21+27]= 1.0000000000000004 ;
    contP[22+27]= 0 ;
    contP[23+27]= 0 ;
    contP[24+27]= 0.9999999999999999 ;
    contP[25+27]= 0.9999999999999999 ;
    contP[26+27]= 0 ;


    APPROXIMATION::get_jacobians(2,contP,ncont,p,true,jacobian);
    jacobianTest[0]= 0.49999999999999967 ;
    jacobianTest[1]= 4.926614671774113e-16 ;
    jacobianTest[2]= 0 ;
    jacobianTest[3]= 6.175615574477411e-16 ;
    jacobianTest[4]= 0.49999999999999845 ;
    jacobianTest[5]= 0 ;
    jacobianTest[6]= 0 ;
    jacobianTest[7]= 0 ;
    jacobianTest[8]= 0 ;

    // calculated with Python
    int temp = PeridigmNS::dof()*PeridigmNS::dof();
    for(int iID=0 ; iID<2 ; ++iID){    
        for(int i=0 ; i<temp ; ++i){    
                TEST_ASSERT(abs(jacobian[i+iID*temp]-jacobianTest[i])<tolerance);  
            }
    }
}


int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}