/*! \file Peridigm_EnergyReleaseDamageCorrepondenceModel.cpp */

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
// Author of this Routine
// Christian Willberg   christian.willberg@dlr.de
// German Aerospace Center
//@HEADER

#include "Peridigm_EnergyReleaseDamageCorrepondenceModel.hpp"
#include "Peridigm_Field.hpp"
#include "material_utilities.h"
#include "correspondence.h"
#include "matrices.h"
#include <thread>
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include "damage_utilities.h"

using namespace std;

PeridigmNS::EnergyReleaseDamageCorrepondenceModel::EnergyReleaseDamageCorrepondenceModel(const Teuchos::ParameterList& params)
: DamageModel(params),
m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
m_applyThermalStrains(false),
m_modelCoordinatesFieldId(-1),
m_coordinatesFieldId(-1),
m_damageFieldId(-1),
m_bondDamageFieldId(-1),
m_bondDamageDiffFieldId(-1),
m_deltaTemperatureFieldId(-1),
m_dilatationFieldId(-1),
m_weightedVolumeFieldId(-1),
m_horizonFieldId(-1),
m_detachedNodesFieldId(-1),
m_piolaStressTimesInvShapeTensorXId(-1),
m_piolaStressTimesInvShapeTensorYId(-1),
m_piolaStressTimesInvShapeTensorZId(-1),
m_forceDensityFieldId(-1),
m_deformationGradientFieldId(-1),
m_hourglassStiffId(-1) {

    
    if (params.isParameter("Critical Energy")) {
        m_criticalEnergyTension = params.get<double>("Critical Energy");
          
    } 
    if (params.isParameter("Degradation Factor")){
        degradationFactor = 1;
        degradationFactor = params.get<double>("Degradation Factor");
        
    }else{
        degradationFactor = 1; 
    }
    m_incremental = false;
    if (params.isParameter("Plastic")){
       m_incremental = params.get<bool>("Plastic");
    }
    m_criticalEnergyTension = params.get<double>("Critical Energy");
    m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");
    m_planeStrain=false;
    m_planeStress=false;
    detachedNodesCheck = false;
    
    
    if (params.isParameter("Detached Nodes Check"))
        detachedNodesCheck = params.get<bool>("Detached Nodes Check");
    if (params.isParameter("Plane Strain"))
      m_planeStrain = params.get<bool>("Plane Strain");
    if (params.isParameter("Plane Stress"))
      m_planeStress = params.get<bool>("Plane Stress");
    m_plane = false;
   
    m_onlyTension = false;
    if(params.isParameter("Only Tension"))
        m_onlyTension = params.get<bool>("Only Tension");
    m_criticalEnergyInterBlock = m_criticalEnergyTension;
    if (params.isParameter("Interblock damage energy")){
        m_criticalEnergyInterBlock = params.get<double>("Interblock damage energy");
        
        for (int iID = 0; iID < 8; ++iID){ block[iID] = 0;}
        if (params.isParameter("Block_12")) {block[0] = params.get<int>("Block_12");}
        if (params.isParameter("Block_21")) {block[1] = params.get<int>("Block_21");}
        if (params.isParameter("Block_34")) {block[2] = params.get<int>("Block_34");}
        if (params.isParameter("Block_43")) {block[3] = params.get<int>("Block_43");}
        if (params.isParameter("Block_56")) {block[4] = params.get<int>("Block_56");}
        if (params.isParameter("Block_65")) {block[5] = params.get<int>("Block_65");}
        if (params.isParameter("Block_78")) {block[6] = params.get<int>("Block_78");}
        if (params.isParameter("Block_87")) {block[7] = params.get<int>("Block_87");}
        
    }
    m_bondDiffSt = 1e10;
    if (params.isParameter("Stable Bond Difference"))
        m_bondDiffSt  = params.get<int>("Stable Bond Difference");
  //************************************
  // wie komme ich an den Namen??
  //************************************
  //if (params.isParameter("Linear Elastic Correspondence")){
   
    //std::cout<<"Use Material: Linear Elastic Correspondence"<<std::endl;

    if (m_planeStrain==true){
        //m_plane=true;
        m_plane = true;
        m_Thickness = params.get<double>("Thickness");
        //std::cout<<"Method 2D Plane Strain"<<std::endl;
    }
    if (m_planeStress==true){
       // m_plane=true;
        m_plane = true;
        m_Thickness = params.get<double>("Thickness");
        //std::cout<<"WRN: Method 2D Plane Stress --> not fully implemented yet"<<std::endl;
    }
   
// Params anpassen. Irgendwo muss das befuellt werden, da bspw. die thermische Verformung bei Material
// und Schaden existiert
//
//

    m_pi = M_PI;

    if (params.isParameter("Thermal Expansion Coefficient")) {
        m_alpha = params.get<double>("Thermal Expansion Coefficient");
        m_applyThermalStrains = true;
    }

    PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
    blockIdFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Block_Id");
    
    m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT,"Model_Coordinates");
    m_coordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
    m_volumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
    m_weightedVolumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
    m_dilatationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation");
    m_damageFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
    m_bondDamageFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
    m_bondDamageDiffFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage_Diff");
    m_horizonFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
    m_piolaStressTimesInvShapeTensorXId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorX");
    m_piolaStressTimesInvShapeTensorYId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorY");
    m_piolaStressTimesInvShapeTensorZId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "PiolaStressTimesInvShapeTensorZ");
    m_detachedNodesFieldId              = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Detached_Nodes");
    m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
    m_deformationGradientFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Deformation_Gradient");
    m_hourglassStiffId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Hourglass_Stiffness");
    m_bondEnergyFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Energy");
    
    
 
    m_fieldIds.push_back(blockIdFieldId);
    m_fieldIds.push_back(m_volumeFieldId);
    m_fieldIds.push_back(m_modelCoordinatesFieldId);
    m_fieldIds.push_back(m_coordinatesFieldId);
    m_fieldIds.push_back(m_weightedVolumeFieldId);
    m_fieldIds.push_back(m_dilatationFieldId);
    m_fieldIds.push_back(m_volumeFieldId);
    m_fieldIds.push_back(m_damageFieldId);
    m_fieldIds.push_back(m_detachedNodesFieldId);
    m_fieldIds.push_back(m_bondDamageFieldId);
    m_fieldIds.push_back(m_bondDamageDiffFieldId);
    m_fieldIds.push_back(m_horizonFieldId);
    m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorXId);
    m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorYId);
    m_fieldIds.push_back(m_piolaStressTimesInvShapeTensorZId);
    m_fieldIds.push_back(m_forceDensityFieldId);
    m_fieldIds.push_back(m_hourglassStiffId);
    m_fieldIds.push_back(m_deformationGradientFieldId);
    m_fieldIds.push_back(m_bondEnergyFieldId);
    

}

PeridigmNS::EnergyReleaseDamageCorrepondenceModel::~EnergyReleaseDamageCorrepondenceModel() {
}

void
PeridigmNS::EnergyReleaseDamageCorrepondenceModel::initialize(const double dt,
        const int numOwnedPoints,
        const int* ownedIDs,
        const int* neighborhoodList,
        PeridigmNS::DataManager& dataManager) const {
    double *damage, *bondDamage;


    dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
    
    dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
    // Initialize damage to zero
    int neighborhoodListIndex(0);
    int bondIndex(0);
    int nodeId, numNeighbors;
    int iID, iNID;

    for (iID = 0; iID < numOwnedPoints; ++iID) {
        nodeId = ownedIDs[iID];
        damage[nodeId] = 0.0;
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        neighborhoodListIndex += numNeighbors;
        for (iNID = 0; iNID < numNeighbors; ++iNID) {
            bondDamage[bondIndex] = 0.0;
            bondIndex += 1;
        }
    }

}

void
PeridigmNS::EnergyReleaseDamageCorrepondenceModel::computeDamage(const double dt,
        const int numOwnedPoints,
        const int* ownedIDs,
        const int* neighborhoodList,
        PeridigmNS::DataManager& dataManager,
        int blockInterfaceId = -1) const {

    double *x, *y, *yN, *damage, *bondDamage, *bondDamageNP1, *bondDamageDiff, *horizon, *vol, *detachedNodes, *blockNumber;
    double *bondEnergyN, *bondEnergyNP1;
    
    
    double criticalEnergyTension(-1.0);
    // for temperature dependencies easy to extent
    //double *deltaTemperature = NULL;
    double *tempStressX, *tempStressY, *tempStressZ;
    dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);

    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);

    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
    if (m_incremental)dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_N)->ExtractView(&yN);
    dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&vol);

    dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamage);
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);
    dataManager.getData(m_bondDamageDiffFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageDiff);
    dataManager.getData(blockIdFieldId, PeridigmField::STEP_NONE)->ExtractView(&blockNumber);
    dataManager.getData(m_bondEnergyFieldId, PeridigmField::STEP_N)->ExtractView(&bondEnergyN);
    dataManager.getData(m_bondEnergyFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondEnergyNP1);

    int iID, nodeId;
    //bool rerun;
    //rerun = false;
    //for (iID = 0; iID < numOwnedPoints; ++iID){
    //    nodeId = ownedIDs[iID];
    //    if (bondDamageDiff[nodeId]>m_bondDiffSt)
    //    {
    //        rerun = true;  
    //        cout << "rerun2: " << bondDamageDiff[nodeId] << " nodeId: " << nodeId << endl;
    //        bondDamageDiff[nodeId] = 0;
    //        break;
    //    }
    //}

    
    dataManager.getData(m_piolaStressTimesInvShapeTensorXId, PeridigmField::STEP_NP1)->ExtractView(&tempStressX);
    dataManager.getData(m_piolaStressTimesInvShapeTensorYId, PeridigmField::STEP_NP1)->ExtractView(&tempStressY);
    dataManager.getData(m_piolaStressTimesInvShapeTensorZId, PeridigmField::STEP_NP1)->ExtractView(&tempStressZ);
    dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
    ////////////////////////////////////////////////////////
    // Hourglass control
    // not synchronized yet; Hourglass correction is used only from nodeId site
    //
    double *defGrad, *hourglassStiff;
    vector<double> TSvector(3), hourglassStiffVector(9);
    double* hStiff = &hourglassStiffVector[0];
    double* TS = &TSvector[0];
    double hourglassScaling;
    ///////////////////////////////////////////////////////
    dataManager.getData(m_hourglassStiffId, PeridigmField::STEP_NONE)->ExtractView(&hourglassStiff);
    dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&defGrad);
    //std::cout<< "heredam"<<std::endl;
    // Set the bond damage to the previous value --> needed for iteration in implicit time integration
    //if(!rerun)
    *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N));
    if (m_incremental){
        *(dataManager.getData(m_bondEnergyFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_bondEnergyFieldId, PeridigmField::STEP_N));
        hourglassScaling = 0.0;
    }
    else{
        dataManager.getData(m_bondEnergyFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
        hourglassScaling = m_hourglassCoefficient;
    }
    double *forceDensity;
    dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);
    
    ////////////////////////////////////////////////////
    double trialDamage(0.0);
    int neighborhoodListIndex(0), bondIndex(0), bondCheck(0);
    int numNeighbors, neighborID = 0, iNID;
    double nodeInitialX[3], nodeCurrentX[3], nodeCurrentXN[3];
    nodeCurrentXN[0] = 0.0;
    nodeCurrentXN[1] = 0.0;
    nodeCurrentXN[2] = 0.0;
    double omegaP1, omegaP2;
    double critIso;
    double quadhorizon;
    // bond force density state
    double TX, TY, TZ, TXN, TYN, TZN;
    // bond force density state projected to deformed bond
    double TPX, TPY, TPZ, TPXN, TPYN, TPZN;
    double factor, factorN;
    // deformed state
    double Y[3]; //_dx, Y_dy, Y_dz;
    // initial state
    double X[3]; //_dx, X_dy, X_dz;
    // displacement vector state and abolute value
    double eta[3], incEta[3], normEtaSq;

    double dX, dY;//, dYSq;

    //---------------------------
    // INITIALIZE PROCESS STEP t
    //---------------------------
    // Foster 2009 Journal for Multiscale Computational Engineering
    // m_criticalEnergy = 4 G / (pi delta^4)
    // what if two horizons meet? two criterions will hit the same bond..
    //double quadhorizon =  4 /( m_pi * m_horizon * m_horizon * m_horizon * m_horizon );
    
    if (m_criticalEnergyTension > 0.0)
        criticalEnergyTension = m_criticalEnergyTension;
 
    // Update the bond damage
    // Break bonds if the bond energy potential is greater than the critical bond energy potential
    //---------------------------
    // DAMAGE ANALYSIS
    //---------------------------
    bondIndex = 0;
    for (iID = 0; iID < numOwnedPoints; ++iID, defGrad+=9, hourglassStiff+=9) {
        numNeighbors = neighborhoodList[neighborhoodListIndex++];
        
        bondCheck = 0;

        nodeId = ownedIDs[iID];
        for (int i = 0; i < 3; ++i){
            nodeInitialX[i] = x[nodeId*3 + i];
            nodeCurrentX[i] = y[nodeId*3 + i];
            if (m_incremental){
                nodeCurrentXN[i] = yN[nodeId*3 + i];       
            }
        }
        for (iNID = 0; iNID < numNeighbors; ++iNID) {
            neighborID = neighborhoodList[neighborhoodListIndex++];
            bool modelActive = true;
            if (m_onlyTension == true){
                // if this option is active bond break only if they are streched
                dX = distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2], x[neighborID*3], x[neighborID*3+1], x[neighborID*3+2]);
                dY = distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2], y[neighborID*3], y[neighborID*3+1], y[neighborID*3+2]);
                if (dY-dX<0) modelActive = false;
            } 
            if (modelActive == true){
                if (detachedNodes[nodeId]!=0) continue;
                if (detachedNodes[neighborID]!=0) continue;
                normEtaSq = 0.0;
                for (int i = 0; i < 3; ++i){
                    X[i] = x[neighborID*3 + i]   - nodeInitialX[i];
                    Y[i] = y[neighborID*3 + i]   - nodeCurrentX[i];
                    if (m_incremental){
                        //projection on bond, but energy only with incremental change
                        eta[i]  = Y[i]-X[i];
                        incEta[i]  = Y[i]-(yN[neighborID*3 + i]   - nodeCurrentXN[i]);
                    }
                    else{
                        eta[i]  = Y[i]-X[i];
                        incEta[i] = eta[i];
                    }
                    normEtaSq+=eta[i]*eta[i];
                }
                /*
                if m_incremental is true an energy accumlation is calculated using the trapezoidal integration;
                therefore bondEnergyN has to set to zero, if m_incremental is false
                */
                if (m_incremental == false) bondEnergyN[bondIndex] = 0.0;

                if (normEtaSq>0){// this is to avoid numerical issues
                    
                    criticalEnergyTension = m_criticalEnergyTension;
                    
                    if (blockNumber[neighborID]==blockInterfaceId)criticalEnergyTension = m_criticalEnergyInterBlock;
                    
                    omegaP1 = MATERIAL_EVALUATION::scalarInfluenceFunction(dX, horizon[nodeId]); 
                    omegaP2 = MATERIAL_EVALUATION::scalarInfluenceFunction(-dX, horizon[neighborID]); 
                    // average Force has to be taken
                    // if not the case where 
                    if (m_plane == true) X[2] = 0;
                    
                    double FxsiX = *(defGrad)   * X[0] + *(defGrad+1) * X[1] + *(defGrad+2) * X[2];
                    double FxsiY = *(defGrad+3) * X[0] + *(defGrad+4) * X[1] + *(defGrad+5) * X[2];
                    double FxsiZ = *(defGrad+6) * X[0] + *(defGrad+7) * X[1] + *(defGrad+8) * X[2];
                    hStiff[0] = *(hourglassStiff  ); hStiff[1] = *(hourglassStiff+1);  hStiff[2] = *(hourglassStiff+2);
                    hStiff[3] = *(hourglassStiff+3); hStiff[4] = *(hourglassStiff+4);  hStiff[5] = *(hourglassStiff+5);
                    hStiff[6] = *(hourglassStiff+6); hStiff[7] = *(hourglassStiff+7);  hStiff[8] = *(hourglassStiff+8);
                    CORRESPONDENCE::computeCorrespondenceStabilityWanEtAlShort(FxsiX,FxsiY,FxsiZ,Y[0],Y[1],Y[2],hStiff,TS);
                    TX  =   omegaP1 * ( tempStressX[3*nodeId]     * X[0] + tempStressX[3*nodeId+1]     * X[1] + tempStressX[3*nodeId+2]     * X[2] + hourglassScaling*TS[0]);
                    TY  =   omegaP1 * ( tempStressY[3*nodeId]     * X[0] + tempStressY[3*nodeId+1]     * X[1] + tempStressY[3*nodeId+2]     * X[2] + hourglassScaling*TS[1]);
                    TZ  =   omegaP1 * ( tempStressZ[3*nodeId]     * X[0] + tempStressZ[3*nodeId+1]     * X[1] + tempStressZ[3*nodeId+2]     * X[2] + hourglassScaling*TS[2]);
                    // undeformedBondX, undeformedBondY, undeformedBondZ of bond 1-2 equal to -undeformedBondX, -undeformedBondY, -undeformedBondZ of bond 2-1
                    TXN =   omegaP2 * ( tempStressX[3*neighborID] * X[0] + tempStressX[3*neighborID+1] * X[1] + tempStressX[3*neighborID+2] * X[2] + hourglassScaling*TS[0]);
                    TYN =   omegaP2 * ( tempStressY[3*neighborID] * X[0] + tempStressY[3*neighborID+1] * X[1] + tempStressY[3*neighborID+2] * X[2] + hourglassScaling*TS[1]);
                    TZN =   omegaP2 * ( tempStressZ[3*neighborID] * X[0] + tempStressZ[3*neighborID+1] * X[1] + tempStressZ[3*neighborID+2] * X[2] + hourglassScaling*TS[2]);
                    //std::cout<< "here2"<<std::endl;
                    // orthogonal projection of T and TN to the relative displacement vector Foster et al. "An energy based .."
                    // --> die senkrecht zur Projektion stehenden Anteile entsprechen eventuell den Schubanteilen. D.h. man könnte das Kriterium hier splitten.
                    factor = (eta[0]*TX + eta[1]*TY + eta[2]*TZ)/normEtaSq;
                    TPX = factor*eta[0]; TPY = factor*eta[1]; TPZ = factor*eta[2];
                    
                    //factor = (Y_dx*TX + Y_dy*TY + Y_dz*TZ)/dYSq;
                    //TPX = factor*Y_dx; TPY = factor*Y_dy; TPZ = factor*Y_dz;
                    
                    factorN = (eta[0]*TXN + eta[1]*TYN + eta[2]*TZN)/normEtaSq;
                    TPXN = factorN*eta[0]; TPYN = factorN*eta[1]; TPZN = factorN*eta[2];
                    
                    //factorN = (Y_dx*TXN + Y_dy*TYN + Y_dz*TZN)/dYSq;
                    //TPXN = factorN*Y_dx; TPYN = factorN*Y_dy; TPZN = factorN*Y_dz;

                    // 0.25 oder 0.5
                    bondEnergyNP1[bondIndex] = bondEnergyN[bondIndex] + 0.25*(1-bondDamageNP1[bondIndex])*(abs(TPX*incEta[0])+abs(TPXN*incEta[0])+abs(TPY*incEta[1])+abs(TPYN*incEta[1])+abs(TPZ*incEta[2])+abs(TPZN*incEta[2]));
                    
                }
                else
                {
                    bondEnergyNP1[bondIndex] = 0;
                }
                
                double avgHorizon = 0.5*(horizon[nodeId]+horizon[neighborID]);
                if (bondEnergyNP1[bondIndex]<0){
                    std::cout<<TPX<<" "<< labs(TPXN)<<" BE "<<bondEnergyNP1[bondIndex]<<" BD  "<< (1-bondDamageNP1[bondIndex])<<std::endl;
                }
                ////////////////////////////////////////////////////////////////
                //--> to check, depth is not included yet. How to handle??
                ////////////////////////////////////////////////////////////////
                if (m_planeStrain==false&&m_planeStress==false){
                   quadhorizon =  4 /( m_pi * avgHorizon * avgHorizon * avgHorizon * avgHorizon );
                }
                else
                {
                   quadhorizon =  3 /( avgHorizon * avgHorizon * avgHorizon * m_Thickness );
                }

                critIso = bondEnergyNP1[bondIndex]/(criticalEnergyTension*quadhorizon);

                trialDamage = 0.0;
                if (criticalEnergyTension > 0.0 && critIso > 1.0) {
                    trialDamage = bondDamageNP1[bondIndex] + degradationFactor;
                }
                if (trialDamage > bondDamageNP1[bondIndex]) {
                    if (trialDamage>1)trialDamage = 1;
                    bondDamageNP1[bondIndex] = trialDamage;
                    bondCheck++;    
                }
            }
            bondDamageDiff[nodeId] = bondCheck;

            if(bondCheck>m_bondDiffSt)
                break;

            bondIndex += 1;

        }

        if(bondCheck>m_bondDiffSt)
            break;

    }
    //  Update the element damage (percent of bonds broken)
    if (detachedNodesCheck == true){
        int check = 1;  //set check = 1 to start the loop
        //std::cout<< "detached"<<std::endl;
        while (check != 0){
            check = checkDetachedNodes(numOwnedPoints, ownedIDs, neighborhoodList, dataManager);
        }
    }

    
    DAMAGE_UTILITIES::calculateDamageIndex(numOwnedPoints,ownedIDs,vol,neighborhoodList,bondDamageNP1, damage);
   
}

int PeridigmNS::EnergyReleaseDamageCorrepondenceModel::checkDetachedNodes(
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,                                                     
                                                      PeridigmNS::DataManager& dataManager
                                                      ) const
{
  double *bondDamageNP1, *detachedNodes, *x, *volume, *horizon;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);
  dataManager.getData(m_detachedNodesFieldId, PeridigmField::STEP_NP1)->ExtractView(&detachedNodes);
 
  int neighborhoodListIndex(0), bondIndex(0);
  int nodeId, numNeighbors, neighborID, iID, iNID;
  double nodeInitialX[3],  initialDistance;
  double checkShapeTensor[9], neighborVolume, checkShapeTensorInv[9];
  double omega;
  double determinant;
  int check = 0, matrixInversionReturnCode = 0;

  for(iID=0 ; iID<numOwnedPoints ; ++iID){
      nodeId = ownedIDs[iID];
      nodeInitialX[0] = x[nodeId*3];
      nodeInitialX[1] = x[nodeId*3+1];
      nodeInitialX[2] = x[nodeId*3+2];

      numNeighbors = neighborhoodList[neighborhoodListIndex++];
      
      checkShapeTensor[0] = 0.0;
      checkShapeTensor[1] = 0.0;
      checkShapeTensor[2] = 0.0;
      checkShapeTensor[3] = 0.0;
      checkShapeTensor[4] = 0.0;
      checkShapeTensor[5] = 0.0;
      checkShapeTensor[6] = 0.0;
      checkShapeTensor[7] = 0.0;
      checkShapeTensor[8] = 0.0;
  
      for(iNID=0 ; iNID<numNeighbors ; ++iNID){
          
          neighborID = neighborhoodList[neighborhoodListIndex++];
          
          if (detachedNodes[nodeId]==0){
          
              neighborVolume = volume[neighborID];
              initialDistance = distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2],
                                  x[neighborID*3], x[neighborID*3+1], x[neighborID*3+2]);
              omega = MATERIAL_EVALUATION::scalarInfluenceFunction(initialDistance, horizon[iID]);
              //double omega = 1.0;
              double temp = (1.0 - bondDamageNP1[bondIndex]) * omega * neighborVolume;
              
              double undeformedBondX =  x[neighborID*3]   - nodeInitialX[0];
              double undeformedBondY =  x[neighborID*3+1] - nodeInitialX[1];
              double undeformedBondZ =  x[neighborID*3+2] - nodeInitialX[2];
              
              checkShapeTensor[0]  += temp * undeformedBondX * undeformedBondX;
              checkShapeTensor[1]  += temp * undeformedBondX * undeformedBondY;
              checkShapeTensor[2]  += temp * undeformedBondX * undeformedBondZ;
              checkShapeTensor[3]  += temp * undeformedBondY * undeformedBondX;
              checkShapeTensor[4]  += temp * undeformedBondY * undeformedBondY;
              checkShapeTensor[5]  += temp * undeformedBondY * undeformedBondZ;
              checkShapeTensor[6]  += temp * undeformedBondZ * undeformedBondX;
              checkShapeTensor[7]  += temp * undeformedBondZ * undeformedBondY;
              checkShapeTensor[8]  += temp * undeformedBondZ * undeformedBondZ;
          }
          if (detachedNodes[neighborID]!=0&&bondDamageNP1[bondIndex] != 1.0){// bondDamage check to avoid infinite loop; all bonds then destroyed
              bondDamageNP1[bondIndex] = 1.0;
              check = 1;
          }
          
          bondIndex += 1;
      }
      
      if (detachedNodes[nodeId]!=0) continue;
      
      if (m_plane==true){
        matrixInversionReturnCode =
        MATRICES::Invert2by2Matrix(checkShapeTensor, determinant, checkShapeTensorInv);
        }
      else{
        matrixInversionReturnCode =
        MATRICES::Invert3by3Matrix(checkShapeTensor, determinant, checkShapeTensorInv);
        }
      
      if (matrixInversionReturnCode != 0){// to be checked
          check = 1;
          detachedNodes[nodeId]=1.;
          int bondIndex2 = bondIndex-numNeighbors; // set index back, to have the same correct bonds
          // delete all connected bonds
          for(iNID=0 ; iNID<numNeighbors ; ++iNID){
                  bondDamageNP1[bondIndex2] = 1.0;
                  bondIndex2 += 1;
          }   
        matrixInversionReturnCode = 0;
      }
  }

  return check;
}


