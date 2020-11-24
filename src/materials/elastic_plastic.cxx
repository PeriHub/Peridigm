//! \file elastic_plastic.cxx

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
#include <cmath>
#include <Sacado.hpp>
#include "elastic_plastic.h"
#include "material_utilities.h"

namespace MATERIAL_EVALUATION {

template<typename ScalarT>
ScalarT computeDeviatoricForceStateNorm
(
        int numNeigh,
        ScalarT theta,
        const int *neighPtr,
        const double *bondDamage,
        const double *deviatoricPlasticExtensionState,
        const double *isotropicPlasticExtensionState,
        const double *X,
        const ScalarT *Y,
        const double *xOverlap,
        const ScalarT *yOverlap,
        const double *volumeOverlap,
        double gamma,
        double horizon, 
        double alpha,
        const ScalarT kappa,
        ScalarT& pressure,
        const bool pressureSensitive
)
{
    ScalarT norm=0.0;
    const double *v = volumeOverlap;
    double cellVolume, dx_X, dy_X, dz_X, zeta, edpN, eipN;
    ScalarT dx_Y, dy_Y, dz_Y, dY, ed, tdTrial;
    double omega;
    ScalarT ti, ei;
    pressure = 0.0;
    for(int n=0;n<numNeigh;n++, neighPtr++, bondDamage++, deviatoricPlasticExtensionState++, isotropicPlasticExtensionState++){
        int localId = *neighPtr;
        cellVolume = v[localId];
        const double *XP = &xOverlap[3*localId];
        const ScalarT *YP = &yOverlap[3*localId];
        dx_X = XP[0]-X[0];
        dy_X = XP[1]-X[1];
        dz_X = XP[2]-X[2];
        zeta = sqrt(dx_X*dx_X+dy_X*dy_X+dz_X*dz_X);
        dx_Y = YP[0]-Y[0];
        dy_Y = YP[1]-Y[1];
        dz_Y = YP[2]-Y[2];
        dY = sqrt(dx_Y*dx_Y+dy_Y*dy_Y+dz_Y*dz_Y);
        // 2D - 3D included via gamma
        
        /*
         * Deviatoric extension state
         */
        ei = theta * gamma * zeta / 3;
        ed = dY-zeta-ei;

        /*
         * Deviatoric plastic extension state from last step
         */
        edpN = *deviatoricPlasticExtensionState;
        
        omega = scalarInfluenceFunction(zeta,horizon);
            
        /*
         * Compute trial stress
         * NOTE: include damage
         */
        tdTrial = (1-*bondDamage) * omega * alpha * (ed - edpN);
        
        if (pressureSensitive){
            eipN = *isotropicPlasticExtensionState;
            ti = (1-*bondDamage) * omega * kappa * (ei - eipN); // passt noch nicht; der Druck muss in den Routinen kontinuierlich bestimmt werden;
            pressure -= (ti + tdTrial) * zeta / 3 * cellVolume; // Silling et al. Peridynamic States and Constitutive Modeling (2007)
        }
        /*
         * Accumulate norm
         */
        norm += tdTrial * tdTrial * cellVolume;
    }

    return sqrt(norm);
}

template<typename ScalarT>
void computeInternalForceIsotropicElasticPlastic
(
        const double* xOverlap,
        const ScalarT* yNP1Overlap,
        const double* mOwned,
        const double* volumeOverlap,
        const ScalarT* dilatationOwned,
        const double* bondDamage,
        const double* deviatoricPlasticExtensionStateN,
        ScalarT* deviatoricPlasticExtensionStateNp1,
        const double* isotropicPlasticExtensionStateN,
        ScalarT* isotropicPlasticExtensionStateNp1,
        const double* lambdaN,
        ScalarT* lambdaNP1,
        ScalarT* fInternalOverlap,
        const int*  localNeighborList,
        int numOwnedPoints,
        double BULK_MODULUS,
        double SHEAR_MODULUS,
        double horizon,
        const double yieldStress,
        const bool planeStrain,
        const bool planeStress,
        const bool pressureSensitive,
        const double thickness,
        const double yieldPsiDP,
        const double yieldBetaDP
)
{
    /*
     * Compute processor local contribution to internal force
     */
    double K = BULK_MODULUS;
    double MU = SHEAR_MODULUS;
    //double OMEGA=1.0;
    double omega;
    bool elastic;
    double psi = 0.0, beta = 0.0;
    double eipN = 0.0, t0 = 0.0;
    double factor;
    
    const double *xOwned = xOverlap;
    const ScalarT *yOwned = yNP1Overlap;
    const double *m = mOwned;
    const double *v = volumeOverlap;
    const ScalarT *theta = dilatationOwned;
    ScalarT *fOwned = fInternalOverlap;
    //ScalarT c = 0.0;
    ScalarT kappa = 0.0;
    double gamma = 0.0;
    const int *neighPtr = localNeighborList;
    double cellVolume, alpha = 0.0, dx_X, dy_X, dz_X, zeta, edpN;
    ScalarT pTrial, tiTrial;
    ScalarT dx_Y, dy_Y, dz_Y, dY, tdTrial, t, ti, td;
    ScalarT ed, ei, pN1 = 0.0;
    ScalarT flowFunction;
    ScalarT tdNorm = 0.0;
    for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, m++, theta++, lambdaN++, lambdaNP1++){

        int numNeigh = *neighPtr; neighPtr++;
        const double *X = xOwned;
        const ScalarT *Y = yOwned;
        double weightedVol = *m;
        
        double selfCellVolume = v[p];
        

        ScalarT deltaLambda=0.0;
        if (planeStress==true){
          //  c = 4.0*K*MU/(3.0*K+4.0*MU) * (*theta) / weightedVol;
            // kappa = c / (*theta) / gamma / 3 (?) (above eq (7) lammiCJ
            kappa = 3.0 * K / weightedVol;
            gamma = 4.0*MU/(3.0*K+4.0*MU);
            alpha = 8.0*MU/weightedVol;
        }
        if (planeStrain==true){
           // c = (12.0*K-4.0*MU) / 9.0 * (*theta) / weightedVol;
            // kappa = c / (*theta) / gamma / 3 (?) (above eq (7) lammiCJ
            //kappa = 3*(12.0*K-4.0*MU) / 18.0 * (*theta) / weightedVol;
            
            gamma = 2.0/3.0;
            kappa = (6.0*K - 2.0*MU) / weightedVol ;
            alpha = 8.0*MU / weightedVol;
        }
        if (planeStress==false and planeStrain==false){
          //  c = 3 * K * (*theta) / weightedVol;
            // kappa = c / (*theta) / gamma / 3 (?) (above eq (7) lammiCJ
            kappa = 9 * K / weightedVol;
            gamma = 1.0;
            alpha = 15.0 * MU / weightedVol;
        }
        /*
        * 2d or 3d variety of yield value (uniaxial stress)
        */
        if(planeStrain or planeStress){ // thickness is in the volume included (tbd to guarantee this)
            
        //yieldValue = 225.0 / 3. * yieldStress * yieldStress / 4 / M_PI / thickness / pow(horizon,4); // have to be checked
            //factor = 8.0 / weightedVol * sqrt(4.0 / 75.0 * M_PI * pow(horizon,5)); // / pow(horizon,4));
            //factor = sqrt(225.0 / 3. / 2 / M_PI / thickness / pow(horizon,4)); // have to be checked
            
            //factor = 16.0 / 5.0 / weightedVol * sqrt(M_PI / 3.0 * pow(horizon,5)); // Simplified with Mathematica
            factor = 3.2 / weightedVol * sqrt(M_PI * pow(horizon,5) / 3.0); // Simplified with Mathematica
        }
        else {
            //yieldValue = 25.0 * yieldStress * yieldStress / 8 / M_PI / pow(horizon,5);
            //factor = sqrt(75.0 / 4 / M_PI / pow(horizon,5));        // equation (51) MitchelJA_2011; 1/2 is not needed 
            //factor = 15.0 / weightedVol * sqrt(4.0 / 75.0 * M_PI * pow(horizon,5));        // equation (51) MitchelJA_2011; 1/2 is not needed 
            factor = 2.0 / weightedVol * sqrt(3.0 * M_PI * pow(horizon,5)); // Simplified with Mathematica
        }
      //  factor = alpha / MU * sqrt(4.0 / 75.0 * M_PI * pow(horizon,5));
        t0 = factor * yieldStress;
        //t0 = sqrt(225.0 / 3. * yieldStress * yieldStress / 4 / M_PI / thickness / pow(horizon,4));
        if (pressureSensitive) {                                                         // equation (22) LammiCJ_2014 also 75 instead of 25
            //betaFactor = yieldStress / (-yieldPressure)*6*sqrt(M_PI*pow(horizon,5)/3);
                beta = yieldBetaDP * factor;
                psi  = yieldPsiDP  * factor; 
               
        }
            else {
                beta = 0.0;
                psi = 0.0;
        }
      //  t0   = yieldStress * 6 / weightedVol * sqrt(M_PI*pow(horizon,5)/3);
      //  t0   = sqrt(75.0 * yieldStress * yieldStress / 4 / M_PI / pow(horizon,5));
      //  t0   = factor * yieldStress;
      //  //std::cout<<sqrt(4.0/75.0)<<"t0 "<< factor<< "alpha " << alpha<<std::endl;
        /*
         * Compute norm of trial stress
         */
        pTrial = 0.0;
        
        tdNorm  = computeDeviatoricForceStateNorm(numNeigh,*theta,neighPtr,bondDamage,deviatoricPlasticExtensionStateN,isotropicPlasticExtensionStateN,X,Y,xOverlap,yNP1Overlap,v,gamma,horizon,alpha, kappa, pTrial, pressureSensitive);
        //if (pTrial>0){
        //    pTrial = 0;
        //}
        //else{
        //    pTrial = pTrial;
        // }
        /*
         * Evaluate yield function
         */
        //double pointWiseYieldValue =  yieldValue;
        // 
        // gleichung 21 - 25 Lammi
        //t0 = -beta * yieldPressure;
        
        // p = -int t*x/3 --> equation (86) in silling 2007
        // P = -1/3*ti*weightedVol/omega*x
        // t0 ist wahrscheinlich yieldvalue (checken)
        // beta = 0.0 if no pressure sensitivity
        //f = tdNorm - beta * pTrial - t0;
       // beta = -beta;
        //f = tdNorm - beta * pTrial - t0;
        // what if pressure positive??
        elastic = true;
        flowFunction =  beta * pTrial + t0;
        if (flowFunction < 0) flowFunction = 0; // Null garantieren --> Checken ob dann der Rest passt
        if(tdNorm > flowFunction){
            //std::cout<<tdNorm <<" bt"<< flowFunction <<" t0"<< pTrial<<std::endl;
            /*
             * This step is incrementally plastic
             */
            elastic = false;
            if (pressureSensitive){
            // beta = 0.0 if no pressure sensitivity
                //deltaLambda=( tdNorm - beta * pTrial - t0) / (alpha + K * beta * psi);
                //deltaLambda=( tdNorm - beta * pTrial - t0) / (alpha + K * beta * psi);
                //deltaLambda=( tdNorm - t0 - beta * pTrial) / (alpha + K * beta * psi); // assume that omega = 1; if not it has to be recalcuated for every step
                deltaLambda=( tdNorm - flowFunction) / (alpha + K * beta * psi); // assume that omega = 1; if not it has to be recalcuated for every step
                
            }
            else{
                deltaLambda=( tdNorm / t0 - 1.0 ) / alpha;
            }
            
            *lambdaNP1 = *lambdaN + deltaLambda;
        } else 
        {
            *lambdaNP1 = *lambdaN;
        }

        for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++, deviatoricPlasticExtensionStateN++, deviatoricPlasticExtensionStateNp1++, isotropicPlasticExtensionStateN++, isotropicPlasticExtensionStateNp1++){
            int localId = *neighPtr;
            cellVolume = v[localId];
            const double *XP = &xOverlap[3*localId];
            const ScalarT *YP = &yNP1Overlap[3*localId];
            dx_X = XP[0]-X[0];
            dy_X = XP[1]-X[1];
            dz_X = XP[2]-X[2];
            zeta = sqrt(dx_X*dx_X+dy_X*dy_X+dz_X*dz_X);
            dx_Y = YP[0]-Y[0];
            dy_Y = YP[1]-Y[1];
            dz_Y = YP[2]-Y[2];
            dY = sqrt(dx_Y*dx_Y+dy_Y*dy_Y+dz_Y*dz_Y);
            
            /*
             * Deviatoric and isotropic plastic extension state from last step
             */
            edpN = *deviatoricPlasticExtensionStateN;
            eipN = *isotropicPlasticExtensionStateN; 
            /*
             * Deviatoric extension state
             */
            /*
             * Compute trial stress
             */
            omega = scalarInfluenceFunction(zeta,horizon);
            ei = *theta * gamma * zeta / 3;
            ed = dY-zeta-ei;

            tdTrial = (1.0-*bondDamage) * alpha * omega * (ed - edpN);
            
            if (pressureSensitive){

                tiTrial = (1.0-*bondDamage) * kappa * omega * (ei - eipN);

            }
            else
            {
                tiTrial = (1.0-*bondDamage) * kappa * omega * ei;
            }
            /*
             * Evaluate yield function
             */
            if(elastic){
                /*
                 * Elastic case
                 */
                td = tdTrial;
                ti = tiTrial;
                /*
                 * Therefore edpNp1 = edpN
                 */
                *deviatoricPlasticExtensionStateNp1 = *deviatoricPlasticExtensionStateN;
                if (pressureSensitive){
                    *isotropicPlasticExtensionStateNp1 = *isotropicPlasticExtensionStateN;
                }
            } else {
                /*
                 * Compute deviatoric force state Eq.(34)
                 */
                if (pressureSensitive){
                    
                    //pN1  = pTrial + deltaLambda * K * beta;
                    //double temp = alpha + K*beta*psi;
                    
                    pN1  = pTrial  + deltaLambda * K * psi;
                    if (pN1 * beta + t0 < 0) pN1 = -t0/beta;
                    ti   = tiTrial - omega * kappa * deltaLambda  * (psi* gamma * zeta / 3.0);
                    *isotropicPlasticExtensionStateNp1  = eipN + deltaLambda * (psi* gamma * zeta / 3.0);
                    
                    //td = tdTrial / (1 + (alpha * omega * deltaLambda / (pN1 * beta + t0)));
                    
                    
                    td = tdTrial * (pN1 * beta + t0) / (pN1 * beta + t0 + alpha * omega * deltaLambda );
                    
                    //*deviatoricPlasticExtensionStateNp1 = edpN + deltaLambda * td;
                    
                    //ti = 3*(K - K*K*beta*psi/temp)*zeta* *theta / weightedVol;
                    //*isotropicPlasticExtensionStateNp1  = eipN + deltaLambda * (psi * gamma * zeta / 3.0);
                    
                    
                    //td = tdTrial / (1 + (alpha * omega * deltaLambda / tdNorm));
                    
                    //*deviatoricPlasticExtensionStateNp1 = edpN + deltaLambda * td / (pN1 * beta + t0);
                    *deviatoricPlasticExtensionStateNp1 = edpN + deltaLambda * tdTrial / (pN1 * beta + t0 + alpha * omega * deltaLambda);
                }
                else
                {
                    ti = tiTrial;
                    pN1 = 0.0;
                    *isotropicPlasticExtensionStateNp1 = 0.0;
                    td = tdTrial / (1 + alpha * deltaLambda);
                     *deviatoricPlasticExtensionStateNp1 = edpN + td * deltaLambda;
                     
                }

            } 
            
            /*
             * Force state (with damage)
             */
            //std::cout<<td << " td" <<  std::endl;
            t = (1.0-*bondDamage) * (ti + td);
                
                
                /*
                 * Update deviatoric plastic deformation state
                 */ 

            /*
             * Assemble pair wise force function
             */
            ScalarT fx = t * dx_Y / dY;
            ScalarT fy = t * dy_Y / dY;
            ScalarT fz = t * dz_Y / dY;

            *(fOwned+0) += fx*cellVolume;
            *(fOwned+1) += fy*cellVolume;
            *(fOwned+2) += fz*cellVolume;
            fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
            fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
            fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
        }
        //std::cout<< edpN << " plast "<< elastic << std::endl;
    }
    
}

/** Explicit template instantiation for double. */
template double computeDeviatoricForceStateNorm<double>
(
        int numNeigh,
        double theta,
        const int *neighPtr,
        const double *bondDamage,
        const double *deviatoricPlasticExtensionState,
        const double *isotropicPlasticExtensionState,
        const double *X,
        const double *Y,
        const double *xOverlap,
        const double *yOverlap,
        const double *volumeOverlap,
        double gamma,
        double horizon,
        double alpha,
        double kappa,
        double& pressure,
        bool pressureSensitive
);

/** Explicit template instantiation for double. */
template void computeInternalForceIsotropicElasticPlastic<double>
(
        const double* xOverlap,
        const double* yNP1Overlap,
        const double* mOwned,
        const double* volumeOverlap,
        const double* dilatationOwned,
        const double* bondDamage,
        const double* deviatoricPlasticExtensionStateN,
        double* deviatoricPlasticExtensionStateNp1,
        const double* isotropicPlasticExtensionStateN,
        double* isotropicPlasticExtensionStateNp1,
        const double* lambdaN,
        double* lambdaNP1,
        double* fInternalOverlap,
        const int*  localNeighborList,
        int numOwnedPoints,
        double BULK_MODULUS,
        double SHEAR_MODULUS,
        double horizon,
        const double yieldStress,
        const bool planeStrain,
        const bool planeStress,
        const bool pressureSensitive,
        const double thickness,
        const double yieldPsiDP,
        const double yieldBetaDP
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template Sacado::Fad::DFad<double> computeDeviatoricForceStateNorm<Sacado::Fad::DFad<double> >
(
        int numNeigh,
        Sacado::Fad::DFad<double> theta,
        const int *neighPtr,
        const double *bondDamage,
        const double *deviatoricPlasticExtensionState,
        const double *isotropicPlasticExtensionState,
        const double *X,
        const Sacado::Fad::DFad<double> *Y,
        const double *xOverlap,
        const Sacado::Fad::DFad<double> *yOverlap,
        const double *volumeOverlap,
        double gamma,
        double horizon,
        double alpha,
        const Sacado::Fad::DFad<double> kappa,
        Sacado::Fad::DFad<double>& pressure,
        bool pressureSensitive
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceIsotropicElasticPlastic<Sacado::Fad::DFad<double> >
(
        const double* xOverlap,
        const Sacado::Fad::DFad<double>* yNP1Overlap,
        const double* mOwned,
        const double* volumeOverlap,
        const Sacado::Fad::DFad<double>* dilatationOwned,
        const double* bondDamage,
        const double* deviatoricPlasticExtensionStateN,
        Sacado::Fad::DFad<double>* deviatoricPlasticExtensionStateNp1,
        const double* isotropicPlasticExtensionStateN,
        Sacado::Fad::DFad<double>* isotropicPlasticExtensionStateNp1,
        const double* lambdaN,
        Sacado::Fad::DFad<double>* lambdaNP1,
        Sacado::Fad::DFad<double>* fInternalOverlap,
        const int*  localNeighborList,
        int numOwnedPoints,
        double BULK_MODULUS,
        double SHEAR_MODULUS,
        double horizon,
        const double yieldStress,
        const bool planeStrain,
        const bool planeStress,
        const bool pressureSensitive,
        const double thickness,
        const double yieldPsiDP,
        const double yieldBetaDP
);

}