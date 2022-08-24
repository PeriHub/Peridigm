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
using namespace std;

PeridigmNS::FEMMaterial::FEMMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_density(0.0)
    
{
     
  m_density = params.get<double>("Density");

  bool m_planeStrain = false, m_planeStress = false;
  m_type = 0;
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
 
  order[0] = params.get<int>("Order");
  order[1] = params.get<int>("Order");
  order[2] = params.get<int>("Order");
  if (params.isParameter("Order_Y")) order[1] = params.get<int>("Order_Y");
  if (params.isParameter("Order_Z")) order[2] = params.get<int>("Order_Z");
  numIntDir[0] = order[0] + 1;
  numIntDir[1] = order[1] + 1;
  numIntDir[2] = order[2] + 1;
  if (params.isParameter("Integration Points in X")) numIntDir[0] = params.get<int>("Integration Points in X");
  if (params.isParameter("Integration Points in Y")) numIntDir[1] = params.get<int>("Integration Points in Y");
  if (params.isParameter("Integration Points in Z")) numIntDir[2] = params.get<int>("Integration Points in Z");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(numIntDir[0]<1, "**** Error:  Number of integration points in x direction must be greater zero.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(numIntDir[1]<1, "**** Error:  Number of integration points in y direction must be greater zero.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(numIntDir[2]<1, "**** Error:  Number of integration points in z direction must be greater zero.\n");
  
  // number of element nodes; if under or over integrated this number is not identical with the number
  // of integration points
  nnode = FEM::getNnode(order, twoD);
   // Integrationpoint topology maps global int number to x-direction number --> safe
  numInt = FEM::getNumberOfIntegrationPoints(twoD, numIntDir);
  delete Bx;
  delete By; 
  delete Bz;
  if (twoD){
    Bx = new double[nnode*numInt];
    By = new double[nnode*numInt];
    Bz = new double[1];
  }
  else
  {
    Bx = new double[nnode*numInt];
    By = new double[nnode*numInt];
    Bz = new double[nnode*numInt];
  }


 
  
  delete weightVector;
  weightVector = new double[numInt];
 

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_modelAnglesId                     = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Local_Angles");
  m_coordinatesFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_displacementFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacements");
  m_unrotatedCauchyStressFieldId      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_cauchyStressFieldId               = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Cauchy_Stress");
  m_partialStressFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Partial_Stress");
  m_nodeTypeFieldId                   = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Node_Type");
  
  

  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_displacementFieldId);

  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);

  m_fieldIds.push_back(m_partialStressFieldId);
  m_fieldIds.push_back(m_modelAnglesId);
  m_fieldIds.push_back(m_nodeTypeFieldId);
}

PeridigmNS::FEMMaterial::~FEMMaterial()
{
}

void
PeridigmNS::FEMMaterial::initialize(const double dt,
                            const int numOwnedPoints,
                            const int* ownedIDs,
                            const int* neighborhoodList,
                            PeridigmNS::DataManager& dataManager)
{
  std::cout<<"init FEM"<<std::endl;

 // FEM::createLumbedMassesSomeHow()
  std::vector<double> NxiVector(order[0]+1), NetaVector(order[1]+1), NpsiVector(order[2]+1);
  double* Nxi  = &NxiVector[0];
  double* Neta = &NetaVector[0];
  double* Npsi = &NpsiVector[0];
  std::vector<double> BxiVector(order[0]+1), BetaVector(order[1]+1), BpsiVector(order[2]+1);
  double* Bxi  = &BxiVector[0];
  double* Beta = &BetaVector[0];
  double* Bpsi = &BpsiVector[0];
  std::vector<double> elCoorxVector(numIntDir[0]), elCooryVector(numIntDir[1]), elCoorzVector(numIntDir[2]);
  double* elCoorx = &elCoorxVector[0];
  double* elCoory = &elCooryVector[0];
  double* elCoorz = &elCoorzVector[0];
  std::vector<double> weightxVector(numIntDir[0]), weightyVector(numIntDir[1]), weightzVector(numIntDir[2]);
  double* weightsx = &weightxVector[0];
  double* weightsy = &weightyVector[0];
  double* weightsz = &weightzVector[0];
  double* weights  = &weightVector[0];
 
  int intPointPtr = 0;
  // temporary vector; the length varies depended on the direction, but its max is nnode
  //std::vector<double> NIntVector(nnode), BIntVector(nnode);
  //double* NInt = &NIntVector[0];
  //double* BInt = &BIntVector[0];
  FEM::weightsAndIntegrationPoints(order[0], elCoorx, weightsx);
  FEM::weightsAndIntegrationPoints(order[1], elCoory, weightsy);

  if (twoD){
    for (int jID=0 ; jID<numIntDir[1] ; ++jID){
      FEM::getLagrangeElementData(order[1],elCoory[jID],Neta,Beta);
      for (int iID=0 ; iID<numIntDir[0] ; ++iID){
        FEM::getLagrangeElementData(order[0],elCoorx[iID],Nxi,Bxi);
        FEM::setElementMatrices(twoD, intPointPtr, order, Nxi, Neta, Npsi, Bxi, Beta, Bpsi, Bx, By, Bz);
        
        intPointPtr += nnode;
      }
    }
  }
  else{
    FEM::weightsAndIntegrationPoints(order[2], elCoorz, weightsz); 
    for (int kID=0 ; kID<numIntDir[2] ; ++kID){
      FEM::getLagrangeElementData(order[2],elCoorz[kID],Npsi,Bpsi);
      for (int jID=0 ; jID<numIntDir[1] ; ++jID){
        FEM::getLagrangeElementData(order[1],elCoory[jID],Neta,Beta);
        for (int iID=0 ; iID<numIntDir[0] ; ++iID){
          FEM::getLagrangeElementData(order[0],elCoorx[iID],Nxi,Bxi);
          FEM::setElementMatrices(twoD, intPointPtr, order, Nxi, Neta, Npsi, Bxi, Beta, Bpsi, Bx, By, Bz);
          intPointPtr += nnode;
        }
      }
    }
  }
  FEM::setWeights(numIntDir,twoD,weightsx,weightsy,weightsz,weights);
  double *nodeType;
  dataManager.getData(m_nodeTypeFieldId, PeridigmField::STEP_NONE)->ExtractView(&nodeType);

  // why it does not work?
  //FEM::getTopology(numOwnedPoints, neighborhoodList, nodeType, numElements, topologyVector);
  
  numElements = 0;
  int topoPtr = 0;
  for(int i=0 ; i<numOwnedPoints ; ++i){
    
    int numNodes = neighborhoodList[topoPtr];
    topoPtr++;
    
    if (nodeType[i]==2){
      numElements += 1;
      topologyVector.push_back(i);
      topologyVector.push_back(numNodes);
      for(int j=0 ; j<numNodes ; ++j){
        topologyVector.push_back(neighborhoodList[topoPtr]);
        topoPtr++;
        }
    }
    else{topoPtr+=numNodes;}
  } 



}
void
PeridigmNS::FEMMaterial::computeForce(const double dt,
                              const int numOwnedPoints,
                              const int* ownedIDs,
                              const int* neighborhoodList,
                              PeridigmNS::DataManager& dataManager,
                              const double currentTime) const
{
  
    double *CauchyStressNP1, *undeformedCoor, *nodeAngles, *deformedCoor, *globalForce;
    //double *displacements;
    //m_unrotatedCauchyStressFieldId    
    //m_cauchyStressFieldId
    dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&CauchyStressNP1);

    dataManager.getData(m_modelAnglesId, PeridigmField::STEP_NONE)->ExtractView(&nodeAngles);
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&undeformedCoor);
    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&deformedCoor);
   // dataManager.getData(m_displacementFieldId, PeridigmField::STEP_NP1)->ExtractView(&displacements);
    dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&globalForce);
    double *volume;
    dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    bool rotation = false;
   // double* disp = displacements;
    double* force = globalForce;
    double* sigmaNP1 = CauchyStressNP1;
    std::vector<double> strainVector(9);
    
    double* strain = &strainVector[0];
    std::vector<double> anglesVector(3);
    double* angles = &anglesVector[0];;
    //int defGradLogReturnCode(0);

    std::vector<double> sigmaIntVector(9);
    double* sigmaInt = &sigmaIntVector[0];

    const int* topology = &topologyVector[0];

    /////////////////////////////////////////////////////////////////////
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
    
    double* weight = &weightVector[0];
    
    std::vector<double> JMat(9);
    double* J = &JMat[0];
    std::vector<double> JinvMat(9);
    double* Jinv = &JinvMat[0];
    double detJ = 0.0;
    int globalId, numElemNodes;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //FEM::getDisplacements(nnode,modelCoordinates, deformedCoor,dispNodal);
    int topoPtr = 0;
    int elementID = 0;
    for(int iID=0 ; iID<numElements ; ++iID){
      // for averaging the element number to which the node is connected has to be known
      elementID = topology[topoPtr];
      topoPtr++;
      numElemNodes = topology[topoPtr];
      topoPtr++;  
      for(int i=0 ; i<3 ; ++i){
        angles[i] = nodeAngles[3*elementID+i]; // element angle is already the average of all element nodes            
        for(int nID=0 ; nID<numElemNodes ; ++nID){
          globalId = topology[topoPtr + nID];
          elNodalCoor[3*nID+i]   = deformedCoor[3*globalId+i];       
          dispNodal[3*nID+i]     = deformedCoor[3*globalId+i] - undeformedCoor[3*globalId+i];
          }
      } 
      if (FEM::vectorNorm(angles, 3)!=0)rotation=true;
      else rotation = false;
      //FEM::setToZero(&sigmaNP1[9 * elementID], 9);
      FEM::setToZero(&elNodalForces[0],3*numElemNodes);

      
      
      for (int jID=0 ; jID<numInt ; ++jID){
        int intPointPtr = jID*numElemNodes;
        // only if nodes and integration points are equal the topology is suitable here.

        detJ=FEM::getJacobian(Bx, By, Bz, numElemNodes, intPointPtr, elNodalCoor, weight[jID], twoD, J, Jinv);
        
        FEM::computeStrain(Bx, By, Bz, dispNodal, intPointPtr, numElemNodes, Jinv, twoD, strain); 
        
        //https://www.continuummechanics.org/stressxforms.html
        // Q Q^T * sigma * Q Q^T = Q C Q^T epsilon Q Q^T
        if (rotation){  
          FEM::tensorRotation(angles,strain,true,strain);
        }
        
        computeCauchyStress(strain, sigmaInt);
        
        // rotation back
        if (rotation){  
          FEM::tensorRotation(angles,sigmaInt,false,sigmaInt);
        }
        FEM::nodalForce(Bx, By, Bz, intPointPtr, numElemNodes, detJ, Jinv, twoD, sigmaInt, elNodalForces);
        // has to be done for each integration point
        // it adds up the different parts of each integration point resulting element force
        FEM::setGlobalStresses(elementID, sigmaInt, sigmaNP1);  
      }
 
      FEM::setNodalStresses(numElemNodes, elementID, topoPtr, topology, sigmaNP1);
      FEM::setGlobalForces(numElemNodes, elementID, topoPtr, topology, elNodalForces, volume, force); 
       
       

      topoPtr+=numElemNodes;

    }
  
 }

