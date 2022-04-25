/*! \file Peridigm_ElasticBondBasedMaterial.cpp */

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

#include "Peridigm_ElasticBondBasedMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_bond_based.h"
#include <Teuchos_Assert.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

PeridigmNS::ElasticBondBasedMaterial::ElasticBondBasedMaterial(const Teuchos::ParameterList &params)
    : Material(params),
      m_bulkModulus(0.0), m_density(0.0), m_horizon(0.0), m_volumeFieldId(-1), m_damageFieldId(-1),
      m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = params.get<double>("Bulk Modulus");
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");
  m_criticalStretch = params.get<double>("Critical Stretch");
  // m_criticalStretch = 10.0;
  if (params.isParameter("Young's Modulus") || params.isParameter("Poisson's Ratio") || params.isParameter("Shear Modulus"))
  {
    TEUCHOS_TEST_FOR_TERMINATION(true, "**** Error:  The Elastic bond based material model supports only one elastic constant, the bulk modulus.");
  }

  m_useCollocationNodes = false;
  if (params.isParameter("Use Collocation Nodes"))
  {
    m_useCollocationNodes = params.get<bool>("Use Collocation Nodes");
  }

  PeridigmNS::FieldManager &fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_damageFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
}

PeridigmNS::ElasticBondBasedMaterial::~ElasticBondBasedMaterial()
{
}

void PeridigmNS::ElasticBondBasedMaterial::initialize(const double dt,
                                                      const int numOwnedPoints,
                                                      const int *ownedIDs,
                                                      const int *neighborhoodList,
                                                      PeridigmNS::DataManager &dataManager)
{
  

  if (m_useCollocationNodes)
  {

  // Extract pointers to the underlying data
  // double *bondDamage;
  
  // dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);

  const int AuxNodesNo = 29;
  double AuxNodes[AuxNodesNo][2];

  useCollocationNodes = new bool[numOwnedPoints];

  AuxNodes[0][0] = 0.;
  AuxNodes[0][1] = -1.;
  AuxNodes[1][0] = -0.5;
  AuxNodes[1][1] = -0.75;
  AuxNodes[2][0] = -0.25;
  AuxNodes[2][1] = -0.75;
  AuxNodes[3][0] = 0.;
  AuxNodes[3][1] = -0.75;
  AuxNodes[4][0] = 0.25;
  AuxNodes[4][1] = -0.75;
  AuxNodes[5][0] = 0.5;
  AuxNodes[5][1] = -0.75;
  AuxNodes[6][0] = -0.75;
  AuxNodes[6][1] = -0.25;
  AuxNodes[7][0] = -0.25;
  AuxNodes[7][1] = -0.25;
  AuxNodes[8][0] = 0.;
  AuxNodes[8][1] = -0.25;
  AuxNodes[9][0] = 0.25;
  AuxNodes[9][1] = -0.25;
  AuxNodes[10][0] = 0.75;
  AuxNodes[10][1] = -0.25;
  AuxNodes[11][0] = -1.;
  AuxNodes[11][1] = 0.;
  AuxNodes[12][0] = -0.75;
  AuxNodes[12][1] = 0.;
  AuxNodes[13][0] = -0.25;
  AuxNodes[13][1] = 0.;
  AuxNodes[14][0] = 0.;
  AuxNodes[14][1] = 0.;
  AuxNodes[15][0] = 0.25;
  AuxNodes[15][1] = 0.;
  AuxNodes[16][0] = 0.75;
  AuxNodes[16][1] = 0.;
  AuxNodes[17][0] = 1.;
  AuxNodes[17][1] = 0.;
  AuxNodes[18][0] = -0.75;
  AuxNodes[18][1] = 0.25;
  AuxNodes[19][0] = -0.25;
  AuxNodes[19][1] = 0.25;
  AuxNodes[20][0] = 0.;
  AuxNodes[20][1] = 0.25;
  AuxNodes[21][0] = 0.25;
  AuxNodes[21][1] = 0.25;
  AuxNodes[22][0] = 0.75;
  AuxNodes[22][1] = 0.25;
  AuxNodes[23][0] = -0.75;
  AuxNodes[23][1] = 0.5;
  AuxNodes[24][0] = -0.25;
  AuxNodes[24][1] = 0.75;
  AuxNodes[25][0] = 0.;
  AuxNodes[25][1] = 0.75;
  AuxNodes[26][0] = 0.25;
  AuxNodes[26][1] = 0.75;
  AuxNodes[27][0] = 0.75;
  AuxNodes[27][1] = 0.5;
  AuxNodes[28][0] = 0.;
  AuxNodes[28][1] = 1.;

  const int MAX_BASIS = 6;
  Eigen::MatrixXd Coeffs = Eigen::MatrixXd::Zero(4, MAX_BASIS);

  Coeffs(0, 0) = 0.;
  Coeffs(0, 1) = 0.;
  Coeffs(0, 2) = 0.;
  Coeffs(0, 3) = 0.7853981633974483;
  Coeffs(0, 4) = 0.2617993877991494;
  Coeffs(0, 5) = 0.;
  Coeffs(1, 0) = 0.;
  Coeffs(1, 1) = 0.;
  Coeffs(1, 2) = 0.;
  Coeffs(1, 3) = 0.;
  Coeffs(1, 4) = 0.;
  Coeffs(1, 5) = 0.2617993877991494;
  Coeffs(2, 0) = 0.;
  Coeffs(2, 1) = 0.;
  Coeffs(2, 2) = 0.;
  Coeffs(2, 3) = 0.;
  Coeffs(2, 4) = 0.;
  Coeffs(2, 5) = 0.2617993877991494;
  Coeffs(3, 0) = 0.;
  Coeffs(3, 1) = 0.;
  Coeffs(3, 2) = 0.;
  Coeffs(3, 3) = 0.2617993877991494;
  Coeffs(3, 4) = 0.7853981633974483;
  Coeffs(3, 5) = 0.;

  double R = 1.0;
  double E = 3 * m_bulkModulus * (1 - 2 * 0.25);

  // std::cout << "E: " << E << std::endl;
  for (int k = 0; k < 4; k++)
  {
    for (int l = 0; l < 6; l++)
    {
      Coeffs(k, l) = Coeffs(k, l) * ((48.0 * E) / (5 * M_PI * m_horizon * m_horizon * m_horizon ));
      
      // std::cout << "Coeffs(k, l): " << Coeffs(k, l) << std::endl;
      // Coeffs(k,l) = (48.00 * E) / (5.0 * M_PI * (R * R * R));
    }
  }

  double *x;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);

  double dx = 0.25;
  int cloudNum;

  int neighborhoodListIndex = 0;
  int collocationNum = 0;

  double **Clouds;
  Clouds = new double *[numOwnedPoints];
  for (int i = 0; i < numOwnedPoints; i++)
  {
    Clouds[i] = new double[AuxNodesNo];
  }

  int cloudNumArray[numOwnedPoints];

  int nodeID, iNID, iID;
  // int bondDamageIndex(0);
  // double damageOnBond;

  for (iID = 0; iID < numOwnedPoints; ++iID)
  {

    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    
    // if (iID == 183242 | iID == 183243 ){
    // if (iID == 146335){
      useCollocationNodes[iID] = m_useCollocationNodes;
    // }
    // else
    // {
    //   useCollocationNodes[iID] = false;
    // }

    nodeID = ownedIDs[iID];

    cloudNum = 0;


    for (int j = 0; j < AuxNodesNo; ++j)
    {
      bool collocationNodeFound = false;

      double xAuxNode = AuxNodes[j][0] + x[nodeID * 3];
      double yAuxNode = AuxNodes[j][1] + x[nodeID * 3 + 1];
      // double zAuxNode = AuxNodes[j][2] + x[nodeID*3+2];

      // if (x[nodeID * 3] == 0 & x[nodeID * 3 + 1] == -100)
      // {
        // std::cout << "numNeighbors: " << numNeighbors << std::endl;
      // }
      iNID = 0;
      for (iNID; iNID < numNeighbors; ++iNID)
      {
        // damageOnBond = bondDamage[bondDamageIndex++];
        // std::cout << "bondDamageIndex: " << bondDamageIndex << std::endl;

        // if(damageOnBond != 0){
        //   std::cout << "damageOnBond: " << damageOnBond << std::endl;
        //   std::cout << "x[nodeID * 3]: " << x[nodeID * 3] << std::endl;
        //   std::cout << "x[nodeID * 3+1]: " << x[nodeID * 3+1] << std::endl;
        //   continue;}

        int neighborID = neighborhoodList[neighborhoodListIndex + iNID];

        // double dist = distance(x[neighborID*3],x[neighborID*3+1],x[neighborID*3+2],xAuxNode,yAuxNode,zAuxNode);
        double dist = distance(x[neighborID * 3], x[neighborID * 3 + 1], 0, xAuxNode, yAuxNode, 0);

        if (dist < dx / 1000.0)
        {
          
          // if (x[nodeID * 3] == 0 & x[nodeID * 3 + 1] == -100)
          // {
          //   std::cout << "x[neighborID * 3]: " << x[neighborID * 3] << std::endl;
          //   std::cout << "x[neighborID * 3+1]: " << x[neighborID * 3 + 1] << std::endl;
          //   std::cout << "dist: " << dist << std::endl;
          //  std::cout << "neighborID: " << neighborID << std::endl;
          // }
          Clouds[iID][cloudNum++] = neighborID;
          collocationNodeFound = true;
          iNID++;
          break;
        }
      }
      // bondDamageIndex-=iNID;
        // std::cout << "bondDamageIndex: " << bondDamageIndex << std::endl;

      if(!collocationNodeFound){
        Clouds[iID][cloudNum++] = -1;
        // if (x[nodeID * 3] == 0 & x[nodeID * 3 + 1] == -100)
        // {
        //   std::cout << "neighborID: " << "-1" << std::endl;
        // }
      }

    }
    // bondDamageIndex += numNeighbors;
    // std::cout << "bondDamageIndex: " << bondDamageIndex << std::endl;

    neighborhoodListIndex += numNeighbors;

    collocationNum += cloudNum + 1;
    cloudNumArray[iID] = cloudNum;
  }

  collocationNeighborhoodList = new int[collocationNum];

  int idx = 0;

  for (int i = 0; i < numOwnedPoints; ++i)
  {
    collocationNeighborhoodList[idx++] = cloudNumArray[i];
    for (int j = 0; j < cloudNumArray[i]; ++j)
    {
      collocationNeighborhoodList[idx++] = Clouds[i][j];
    }
  }

  double rmax = R;
  double rm = 2.0 * rmax;
  double cm = 0.25 * rmax;

  // Eigen::MatrixXd* UpdateMat = new Eigen::MatrixXd[numOwnedPoints];
  UpdateMat = new Eigen::MatrixXd[numOwnedPoints];
  // UpdateMat.resize(numOwnedPoints);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(MAX_BASIS, MAX_BASIS);
  Eigen::VectorXd p = Eigen::VectorXd::Zero(MAX_BASIS);

  neighborhoodListIndex = 0;

  for (int iID = 0; iID < numOwnedPoints; ++iID)
  {

    A.setZero();
    p.setZero();

    nodeID = ownedIDs[iID];

    double xi = x[nodeID * 3];
    double yi = x[nodeID * 3 + 1];

    int numCollocationNeighbors = collocationNeighborhoodList[neighborhoodListIndex++];

    // std::cout << "numCollocationNeighbors: " << numCollocationNeighbors << std::endl;

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(MAX_BASIS, AuxNodesNo);
    Eigen::MatrixXd MMatrix = Eigen::MatrixXd::Zero(MAX_BASIS, AuxNodesNo);

    for (int j = 0; j < numCollocationNeighbors; j++)
    {
      double xj = 0;
      double yj = 0;

      if (j == 14)
      {
        xj = xi;
        yj = yi;
        neighborhoodListIndex++;
      }
      else
      {
        int neighborID = collocationNeighborhoodList[neighborhoodListIndex++];
        // if (xi == 0 & yi == -100)
        // {
        //   std::cout << "neighborID: " << neighborID << std::endl;
        // }
        if (neighborID == -1)
          continue;

        xj = x[neighborID * 3];
        yj = x[neighborID * 3 + 1];
      }

      // int n = Clouds(iID, j);
      // if (n == -1) continue;

      double r = std::sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj));
      // if (xi == 0 & yi == -100)
      // {
      //   std::cout << "numCollocationNeighbors: " << numCollocationNeighbors << std::endl;
      //   std::cout << "xj: " << xj << std::endl;
      //   std::cout << "yj: " << yj << std::endl;
      //   std::cout << "r: " << r << std::endl;
      // }

      double w = weight(r, cm, rm);
      // if (xi == 0 & yi == -100)
      // {
      //   std::cout << "weight: " << w << std::endl;
      // }

      p = Poly(xj - xi, yj - yi);

      // if (xi == 0 & yi == -100)
      // {
      //   for (int k = 0; k < 6; k++)
      //   {
      //     std::cout << "p[" << k << "]: " << p[k] << std::endl;
      //   }
      // }

      A += w * p * p.transpose();

      B.col(j) = w * p;
    }

    // if (xi == 0 & yi == -100)
    // {
    //   for (int k = 0; k < 6; k++)
    //   {
		// 		for (int l = 0; l < 6; l++)
		// 		{
		// 			std::cout << "A(" << k << "," << l << "): " << A(k, l) << std::endl;
		// 		}
		// 		for (int l = 0; l < numCollocationNeighbors; l++)
		// 		{
		// 			std::cout << "B(" << k << "," << l << "): " << B(k, l) << std::endl;
		// 		}
    //   }
    // }

    MMatrix = A.completeOrthogonalDecomposition().pseudoInverse() * B;

    UpdateMat[iID].resize(4, AuxNodesNo);

    UpdateMat[iID].row(0) = Coeffs.row(0) * MMatrix;
    UpdateMat[iID].row(1) = Coeffs.row(1) * MMatrix;
    UpdateMat[iID].row(2) = Coeffs.row(2) * MMatrix;
    UpdateMat[iID].row(3) = Coeffs.row(3) * MMatrix;

    // if (xi == 0 & yi == -100)
    // {
    //   for (int k = 0; k < 4; k++)
    //   {
    //     for (int l = 0; l < numCollocationNeighbors; l++)
    //     {
    //       // if (UpdateMat[iID](k, l) != 0)
    //       // {
    //       std::cout << "xi: " << xi << std::endl;
    //       std::cout << "yi: " << yi << std::endl;
    //       std::cout << "MMatrix(" << k << "," << l << "): " << MMatrix(k, l) << std::endl;
    //       std::cout << "UpdateMat[" << iID << "](" << k << "," << l << "): " << UpdateMat[iID](k, l) << std::endl;
    //       // }
    //     }
    //   }
    // }
  }
}
}

void PeridigmNS::ElasticBondBasedMaterial::computeForce(const double dt,
                                                        const int numOwnedPoints,
                                                        const int *ownedIDs,
                                                        const int *neighborhoodList,
                                                        PeridigmNS::DataManager &dataManager,
                                                        const double currentTime) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *cellVolume, *bondDamage, *force;

  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);


  if (m_useCollocationNodes)
  {
  MATERIAL_EVALUATION::computeInternalForceElasticBondBasedCollocation(x, y, cellVolume, bondDamage, force, neighborhoodList, collocationNeighborhoodList, useCollocationNodes, UpdateMat, numOwnedPoints, m_bulkModulus, m_horizon, m_criticalStretch);
  }else{
    
  MATERIAL_EVALUATION::computeInternalForceElasticBondBased(x, y, cellVolume, bondDamage, force, neighborhoodList, numOwnedPoints, m_bulkModulus, m_horizon);
  }
}
