/*! \file Peridigm_CorrespondenceMaterial.cpp */

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

#include "Peridigm_FEM.hpp"
#include "Peridigm_Field.hpp"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include "FEM_routines.h"
#include "matrices.h"
using namespace std;

PeridigmNS::FEMMaterial::FEMMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_density(0.0)
    
{
     
  m_density = params.get<double>("Density");

  bool m_planeStrain = false, m_planeStress = false;
  m_type = 0;
  order[0] = params.get<int>("Order");
  order[1] = params.get<int>("Order");
  order[2] = params.get<int>("Order");
  if (params.isParameter("Order_Y")) order[1] = params.get<int>("Order_Y");
  if (params.isParameter("Order_Z")) order[2] = params.get<int>("Order_Z");

  delete NxiVector;
  delete NetaVector;
  delete NpsiVector;
  NxiVector  = new double[order[0]+1];
  NetaVector = new double[order[1]+1];
  NpsiVector = new double[order[2]+1];
  delete BxiVector;
  delete BetaVector;
  delete BpsiVector;
  BxiVector  = new double[order[0]+1];
  BetaVector = new double[order[1]+1];
  BpsiVector = new double[order[2]+1];
  delete elCoorxVector;
  delete elCooryVector;
  delete elCoorzVector;
  elCoorxVector = new double[order[0]+1];
  elCooryVector = new double[order[1]+1];
  elCoorzVector = new double[order[2]+1];
  delete weightxVector;
  delete weightyVector; 
  delete weightzVector;
  weightxVector = new double[order[0]+1];
  weightyVector = new double[order[1]+1];
  weightzVector = new double[order[2]+1];




  if (params.isParameter("Plane Strain")){
    m_planeStrain = params.get<bool>("Plane Strain");
    }
  if (params.isParameter("Plane Stress")){
    m_planeStress = params.get<bool>("Plane Stress");
    }
  if (m_planeStrain==true)m_type=1;
  if (m_planeStress==true)m_type=2;
  twoD = false;
  if (m_type != 0)twoD = true;
  numInt = FEM::getNumberOfIntegrationPoints(twoD, order);
  delete localELtopo;
  localELtopo = new int[3*numInt];

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_modelCoordinatesFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_modelAnglesId                     = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Local_Angles");
  m_coordinatesFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_displacementFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacements");
  m_unrotatedCauchyStressFieldId      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_cauchyStressFieldId               = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Cauchy_Stress");
  m_partialStressFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Partial_Stress");
  
  

  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_displacementFieldId);

  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);

  m_fieldIds.push_back(m_partialStressFieldId);
  m_fieldIds.push_back(m_modelAnglesId);
  
}

PeridigmNS::FEMMaterial::~FEMMaterial()
{
}

void
PeridigmNS::FEMMaterial::initialize(const double dt,
                            const int numOwnedPoints,
                            const int* ownedIDs,
                            const int* topology,
                            PeridigmNS::DataManager& dataManager,
                            const int numElements)
{
  std::cout<<"init FEM"<<std::endl;
 // FEM::createLumbedMassesSomeHow()
  double* Nxi  = &NxiVector[0];
  double* Neta = &NetaVector[0];
  double* Npsi = &NpsiVector[0];
  
  double* Bxi  = &BxiVector[0];
  double* Beta = &BetaVector[0];
  double* Bpsi = &BpsiVector[0];
  
  double* elCoorx = &elCoorxVector[0];
  double* elCoory = &elCooryVector[0];
  double* elCoorz = &elCoorzVector[0];
  
  double* weightsx = &weightxVector[0];
  double* weightsy = &weightyVector[0];
  double* weightsz = &weightzVector[0];
  int* elTopo = &localELtopo[0];

  FEM::weightsAndIntegrationPoints(order[0], elCoorx, weightsx);
  FEM::weightsAndIntegrationPoints(order[1], elCoory, weightsy);
 
  FEM::getElementTopo(twoD, order, elTopo);
    
  for (int iID=0 ; iID<order[0]+1 ; ++iID){
    FEM::getLagrangeElementData(order[0],elCoorx[iID],Nxi,Bxi);
  }  
  for (int iID=0 ; iID<order[1]+1 ; ++iID){
    FEM::getLagrangeElementData(order[1],elCoory[iID],Neta,Beta);
  }  
  if (twoD == false){
    FEM::weightsAndIntegrationPoints(order[2],elCoorz,weightsz);
    for (int iID=0 ; iID<order[2]+1 ; ++iID){
      FEM::getLagrangeElementData(order[2],elCoorz[iID],Npsi,Bpsi);
    }  
  }
  //ggf. BxiNetaNpsi vektoren erstellen
    
  nnode = FEM::getNnode(order, twoD);
}

void
PeridigmNS::FEMMaterial::computeForce(const double dt,
                              const int numOwnedPoints,
                              const int* ownedIDs,
                              const int* topologyTemporaryNotUsed,
                              PeridigmNS::DataManager& dataManager,
                              const int numElementsTemporaryNotUsed) const
{
  
    double *CauchyStressNP1, *modelCoordinates, *nodeAngles, *displacements, *deformedCoor, *globalForce;
  

    dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&CauchyStressNP1);

    dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&nodeAngles);
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&deformedCoor);
    dataManager.getData(m_displacementFieldId, PeridigmField::STEP_NP1)->ExtractView(&displacements);
    dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&globalForce);




    bool rotation = false;
    const double* nodalCoor = modelCoordinates;
    double* disp = displacements;
    double* force = globalForce;
    double* sigmaNP1 = CauchyStressNP1;
    std::vector<double> strainVector(9);
    
    double* strain = &strainVector[0];
    double angles[3];
    //int defGradLogReturnCode(0);

    std::vector<double> sigmaIntVector(9);
    double* sigmaInt = &sigmaIntVector[0];
    //int* topoPtr = &topology[0];
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // for higher order it has to be adapted
    // it is a full integration
    // reduced integration is not included yet
    /////////////////////////////////////////////////////////////////////
  
    bool twoD = false;
    if (m_type!=0) {twoD = true;}

    
    int ndof = 3*nnode;

    std::vector<double> dispNodalVector(ndof);
    double* dispNodal = &dispNodalVector[0];
    std::vector<double> elNodalCoorVector(ndof);
    double* elNodalCoor = &elNodalCoorVector[0];
    std::vector<double> elNodalForceVector(ndof);
    double* elNodalForces = &elNodalForceVector[0];
    
    double* Nxi  = &NxiVector[0];
    double* Neta = &NetaVector[0];
    double* Npsi = &NpsiVector[0];
    
    double* Bxi  = &BxiVector[0];
    double* Beta = &BetaVector[0];
    double* Bpsi = &BpsiVector[0];
    
    double* weightsx = &weightxVector[0];
    double* weightsy = &weightyVector[0];
    double* weightsz = &weightzVector[0];
    int* elTopo = &localELtopo[0];
    int topoPtr = 0;
    std::vector<double> JMat(9);
    double* J = &JMat[0];
    std::vector<double> JinvMat(9);
    double* Jinv = &JinvMat[0];
    double detJ = 0.0, detJw = 0.0;
    int localId, numElemNodes;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // irgendwie in die init?
    // das wäre der Austausch für verschiedene Elementtypen
    // falsch
    int numElements = 3;
    bool test = true;
    std::vector<int> topologyVector(15);
    int* topology = &topologyVector[0];
    if (test){
      topology[0] = 4;
      topology[1] = 0;
      topology[2] = 1;
      topology[3] = 4;
      topology[4] = 5;
      topology[5] = 4;
      topology[6] = 1;
      topology[7] = 2;
      topology[8] = 5;
      topology[9] = 6;
      topology[10] = 4;
      topology[11] = 2;
      topology[12] = 3;
      topology[13] = 6;
      topology[14] = 7;
    }


    for(int iID=0 ; iID<numElements ; ++iID){
      // for averaging the element number to which the node is connected has to be known

      numElemNodes = topology[topoPtr];
      topoPtr++;
      angles[0] = 0.0;angles[1] = 0.0;angles[2] = 0.0;
      FEM::getDisplacements(nnode,modelCoordinates, deformedCoor,displacements);
      for(int n=0 ; n<numElemNodes ; ++n){
        localId = topology[topoPtr + n];
        elNodalCoor[3*n]   = nodalCoor[3*localId];
        elNodalCoor[3*n+1] = nodalCoor[3*localId+1];
        elNodalCoor[3*n+2] = nodalCoor[3*localId+2];
        disp[3*localId] = deformedCoor[3*localId]- nodalCoor[3*localId];
        disp[3*localId+1] = deformedCoor[3*localId+1]- nodalCoor[3*localId+1];
        disp[3*localId+2] = deformedCoor[3*localId+2]- nodalCoor[3*localId+2];
        dispNodal[3*n]     = disp[3*localId];
        dispNodal[3*n+1]   = disp[3*localId+1];
        dispNodal[3*n+2]   = disp[3*localId+2];
       // std::cout<<"u1 "<< iID<< " "<<dispNodal[3*n]<< " "<<dispNodal[3*n+1]<< " "<<dispNodal[3*n+2] <<" localId "<<localId<<std::endl;
        sigmaNP1[9*localId  ] = 0.0;sigmaNP1[9*localId+1] = 0.0; sigmaNP1[9*localId+2] = 0.0;
        sigmaNP1[9*localId+3] = 0.0;sigmaNP1[9*localId+4] = 0.0; sigmaNP1[9*localId+5] = 0.0;
        sigmaNP1[9*localId+6] = 0.0;sigmaNP1[9*localId+7] = 0.0; sigmaNP1[9*localId+8] = 0.0;
        angles[0] += nodeAngles[3*localId]/numElemNodes;angles[1] += nodeAngles[3*localId+1]/numElemNodes;angles[2] += nodeAngles[3*localId+2]/numElemNodes;
      }

      for (int jID=0 ; jID<numInt ; ++jID){
        // only if nodes and integration points are equal the topology is suitable here.
        
        detJ=FEM::getJacobian(Nxi,Neta,Npsi,Bxi,Beta,Bpsi,ndof,localELtopo,elNodalCoor, twoD, J, Jinv);
        FEM::computeStrain(Nxi,Neta,Npsi,Bxi,Beta,Bpsi,localELtopo,dispNodal, ndof, Jinv, twoD, strain); 
        
        //https://www.continuummechanics.org/stressxforms.html
        // Q Q^T * sigma * Q Q^T = Q C Q^T epsilon Q Q^T
        if (rotation){  
          MATRICES::tensorRotation(angles,strain,true,strain);
        }
        //CORRESPONDENCE::updateElasticCauchyStressAnisotropicCode(strain, sigmaInt, Cstiff, type);
        computeCauchyStress(strain, sigmaInt);
        // rotation back
        if (rotation){  
          MATRICES::tensorRotation(angles,sigmaInt,false,sigmaInt);
        }
        detJw = FEM::addWeights(detJ, twoD, jID, elTopo, weightsx,weightsy,weightsz);
        FEM::getNodalForce(Nxi,Neta,Npsi,Bxi,Beta,Bpsi,elTopo, sigmaInt, ndof, detJw, Jinv, twoD, elNodalForces);
        // has to be done for each integration point
        // it adds up the different parts of each integration point resulting element force
        
        FEM::getGlobalForcesAndElementStresses(numElemNodes,numInt, topoPtr, topology, elNodalForces, sigmaInt, force, sigmaNP1);

        //topology -= numNeigh;
              // avarage stresses
              // sigmaNP1 /= numInt;
              //globForce(topo) += force; ??
 
      }

      topoPtr+=numElemNodes;
 
 

    }
  
 }

