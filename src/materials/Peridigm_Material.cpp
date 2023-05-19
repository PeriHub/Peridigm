/*! \file Peridigm_Material.cpp */

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

#include "Peridigm_Material.hpp"
#include "Peridigm_Field.hpp"
#include "Peridigm_DegreesOfFreedomManager.hpp"
#include "Peridigm_Logging.hpp"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <cmath>
#include <correspondence.h> // For the semi-Lagrangian (Hypoelastic) models

using namespace std;

void PeridigmNS::Material::computeJacobian(const double dt,
                                           const int numOwnedPoints,
                                           const int* ownedIDs,
                                           const int* neighborhoodList,
                                           PeridigmNS::DataManager& dataManager,
                                           PeridigmNS::SerialMatrix& jacobian,
                                           PeridigmNS::Material::JacobianType jacobianType) const
{
  // Compute a finite-difference Jacobian using either FORWARD_DIFFERENCE or CENTRAL_DIFFERENCE
  computeFiniteDifferenceJacobian(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian, CENTRAL_DIFFERENCE, jacobianType);
}

void PeridigmNS::Material::computeFiniteDifferenceJacobian(const double dt,
                                                           const int numOwnedPoints,
                                                           const int* ownedIDs,
                                                           const int* neighborhoodList,
                                                           PeridigmNS::DataManager& dataManager,
                                                           PeridigmNS::SerialMatrix& jacobian,
                                                           FiniteDifferenceScheme finiteDifferenceScheme,
                                                           PeridigmNS::Material::JacobianType jacobianType) const
{
  // The mechanics Jacobian is of the form:
  //
  // dF_0x/dx_0  dF_0x/dy_0  dF_0x/dz_0  dF_0x/dx_1  dF_0x/dy_1  dF_0x/dz_1  ...  dF_0x/dx_n  dF_0x/dy_n  dF_0x/dz_n
  // dF_0y/dx_0  dF_0y/dy_0  dF_0y/dz_0  dF_0y/dx_1  dF_0y/dy_1  dF_0y/dz_1  ...  dF_0y/dx_n  dF_0y/dy_n  dF_0y/dz_n
  // dF_0z/dx_0  dF_0z/dy_0  dF_0z/dz_0  dF_0z/dx_1  dF_0z/dy_1  dF_0z/dz_1  ...  dF_0z/dx_n  dF_0z/dy_n  dF_0z/dz_n
  // dF_1x/dx_0  dF_1x/dy_0  dF_1x/dz_0  dF_1x/dx_1  dF_1x/dy_1  dF_1x/dz_1  ...  dF_1x/dx_n  dF_1x/dy_n  dF_1x/dz_n
  // dF_1y/dx_0  dF_1y/dy_0  dF_1y/dz_0  dF_1y/dx_1  dF_1y/dy_1  dF_1y/dz_1  ...  dF_1y/dx_n  dF_1y/dy_n  dF_1y/dz_n
  // dF_1z/dx_0  dF_1z/dy_0  dF_1z/dz_0  dF_1z/dx_1  dF_1z/dy_1  dF_1z/dz_1  ...  dF_1z/dx_n  dF_1z/dy_n  dF_1z/dz_n
  //     .           .           .           .           .           .                .           .           .
  //     .           .           .           .           .           .                .           .           .
  //     .           .           .           .           .           .                .           .           .
  // dF_nx/dx_0  dF_nx/dy_0  dF_nx/dz_0  dF_nx/dx_1  dF_nx/dy_1  dF_nx/dz_1  ...  dF_nx/dx_n  dF_nx/dy_n  dF_nx/dz_n
  // dF_ny/dx_0  dF_ny/dy_0  dF_ny/dz_0  dF_ny/dx_1  dF_ny/dy_1  dF_ny/dz_1  ...  dF_ny/dx_n  dF_ny/dy_n  dF_ny/dz_n
  // dF_nz/dx_0  dF_nz/dy_0  dF_nz/dz_0  dF_nz/dx_1  dF_nz/dy_1  dF_nz/dz_1  ...  dF_nz/dx_n  dF_nz/dy_n  dF_nz/dz_n

  // Each entry is computed by finite difference:
  //
  // Forward difference:
  // dF_0x/dx_0 = ( F_0x(perturbed x_0) - F_0x(unperturbed) ) / epsilon
  //
  // Central difference:
  // dF_0x/dx_0 = ( F_0x(positive perturbed x_0) - F_0x(negative perturbed x_0) ) / ( 2.0*epsilon )

  TestForTermination(m_finiteDifferenceProbeLength == DBL_MAX, "**** Finite-difference Jacobian requires that the \"Finite Difference Probe Length\" parameter be set.\n");
  double epsilon = m_finiteDifferenceProbeLength;

  PeridigmNS::DegreesOfFreedomManager& dofManager = PeridigmNS::DegreesOfFreedomManager::self();
  bool solveForDisplacement = dofManager.displacementTreatedAsUnknown();
  bool solveForTemperature = dofManager.temperatureTreatedAsUnknown();
  int numDof = dofManager.totalNumberOfDegreesOfFreedom();
  int numDisplacementDof = dofManager.numberOfDisplacementDegreesOfFreedom();
  //int numTemperatureDof = dofManager.numberOfTemperatureDegreesOfFreedom();
  int displacementDofOffset = dofManager.displacementDofOffset();
  int temperatureDofOffset = dofManager.temperatureDofOffset();

  // Get field ids for all relevant data
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  int volumeFId(-1), coordinatesFId(-1), velocityFId(-1), forceDensityFId(-1), temperatureFId(-1), fluxDivergenceFId(-1);
  volumeFId = fieldManager.getFieldId("Volume");
  if (solveForDisplacement) {
    coordinatesFId = fieldManager.getFieldId("Coordinates");
    velocityFId = fieldManager.getFieldId("Velocity");
    forceDensityFId = fieldManager.getFieldId("Force_Density");
  }
  if (solveForTemperature) {
    temperatureFId = fieldManager.getFieldId("Temperature");
    fluxDivergenceFId = fieldManager.getFieldId("Flux_Divergence");
  }

  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    // Create a temporary neighborhood consisting of a single point and its neighbors.
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    vector<int> tempMyGlobalIDs(numNeighbors+1);
    // Put the node at the center of the neighborhood at the beginning of the list.
    tempMyGlobalIDs[0] = dataManager.getOwnedScalarPointMap()->GID(iID);
    vector<int> tempNeighborhoodList(numNeighbors+1);
    tempNeighborhoodList[0] = numNeighbors;
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      tempMyGlobalIDs[iNID+1] = dataManager.getOverlapScalarPointMap()->GID(neighborID);
      tempNeighborhoodList[iNID+1] = iNID+1;
    }

    Epetra_SerialComm serialComm;
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 3, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));

    // Create a temporary DataManager containing data for this point and its neighborhood.
    PeridigmNS::DataManager tempDataManager;
    tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
                            tempOneDimensionalMap,
                            Teuchos::RCP<const Epetra_BlockMap>(),
                            tempThreeDimensionalMap,
                            tempBondMap);

    // The temporary data manager will have the same fields and data as the real data manager.
    vector<int> fieldIds = dataManager.getFieldIds();
    tempDataManager.allocateData(fieldIds);
    tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);

    // Set up numOwnedPoints and ownedIDs.
    // There is only one owned ID, and it has local ID zero in the tempDataManager.
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(1);
    tempOwnedIDs[0] = 0;

    // Extract pointers to the underlying data.
    double *volume, *y, *v, *force, *temperature, *fluxDivergence;
    tempDataManager.getData(volumeFId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    if (solveForDisplacement) {
      tempDataManager.getData(coordinatesFId, PeridigmField::STEP_NP1)->ExtractView(&y);
      tempDataManager.getData(velocityFId, PeridigmField::STEP_NP1)->ExtractView(&v);
      tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1)->ExtractView(&force);
    }
    if (solveForTemperature) {
      tempDataManager.getData(temperatureFId, PeridigmField::STEP_NP1)->ExtractView(&temperature);
      tempDataManager.getData(fluxDivergenceFId, PeridigmField::STEP_NP1)->ExtractView(&fluxDivergence);
    }

    // Create a temporary vector for storing force and/or flux divergence.
    Teuchos::RCP<Epetra_Vector> forceVector, tempForceVector, fluxDivergenceVector, tempFluxDivergenceVector;
    double *tempForce, *tempFluxDivergence;
    if (solveForDisplacement) {
      forceVector = tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1);
      tempForceVector = Teuchos::rcp(new Epetra_Vector(*forceVector));
      tempForceVector->ExtractView(&tempForce);
    }
    if (solveForTemperature) {
      fluxDivergenceVector = tempDataManager.getData(fluxDivergenceFId, PeridigmField::STEP_NP1);
      tempFluxDivergenceVector = Teuchos::rcp(new Epetra_Vector(*fluxDivergenceVector));
      tempFluxDivergenceVector->ExtractView(&tempFluxDivergence);
    }

    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
    // Resize scratchMatrix if necessary
    if(scratchMatrix.Dimension() < numDof*(numNeighbors+1))
      scratchMatrix.Resize(numDof*(numNeighbors+1));

    // Create a list of global indices for the rows/columns in the scratch matrix.
    vector<int> globalIndices(numDof*(numNeighbors+1));
    for(int i=0 ; i<numNeighbors+1 ; ++i){
      int globalID = tempOneDimensionalMap->GID(i);
      for(int j=0 ; j<numDof ; ++j){
        globalIndices[numDof*i+j] = numDof*globalID+j;
      }
    }

    if(finiteDifferenceScheme == FORWARD_DIFFERENCE){
      if (solveForDisplacement) {
        // Compute and store the unperturbed force.
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        for(int i=0 ; i<forceVector->MyLength() ; ++i)
          tempForce[i] = force[i];
      }
      if (solveForTemperature) {
        // Compute and store the unperturbed flux divergence.
        computeFluxDivergence(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        for(int i=0 ; i<fluxDivergenceVector->MyLength() ; ++i)
          tempFluxDivergence[i] = fluxDivergence[i];
      }
    }

    // Perturb one dof in the neighborhood at a time and compute the force and/or flux divergence.
    // The point itself plus each of its neighbors must be perturbed.
    for(int iNID=0 ; iNID<numNeighbors+1 ; ++iNID){

      int perturbID;
      if(iNID < numNeighbors)
        perturbID = tempNeighborhoodList[iNID+1];
      else
        perturbID = 0;

      // Displacement degrees of freedom
      for(int dof=0 ; dof<numDisplacementDof ; ++dof){

        // Perturb a dof and compute the forces.
        double oldY = y[numDof*perturbID+dof];
        double oldV = v[numDof*perturbID+dof];

        if(finiteDifferenceScheme == CENTRAL_DIFFERENCE){
          // Compute and store the negatively perturbed force.
          y[numDof*perturbID+dof] -= epsilon;
          v[numDof*perturbID+dof] -= epsilon/dt;
          computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
          y[numDof*perturbID+dof] = oldY;
          v[numDof*perturbID+dof] = oldV;
          for(int i=0 ; i<forceVector->MyLength() ; ++i)
            tempForce[i] = force[i];
        }

        // Compute the purturbed force.
        y[numDof*perturbID+dof] += epsilon;
        v[numDof*perturbID+dof] += epsilon/dt;
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        y[numDof*perturbID+dof] = oldY;
        v[numDof*perturbID+dof] = oldV;

        for(int i=0 ; i<numNeighbors+1 ; ++i){
          int forceID;
          if(i < numNeighbors)
            forceID = tempNeighborhoodList[i+1];
          else
            forceID = 0;

          for(int d=0 ; d<numDof ; ++d){
            double value = ( force[numDof*forceID+d] - tempForce[numDof*forceID+d] ) / epsilon;
            if(finiteDifferenceScheme == CENTRAL_DIFFERENCE)
              value *= 0.5;
            scratchMatrix(numDof*forceID + displacementDofOffset + d, numDof*perturbID + displacementDofOffset + dof) = value;
          }
        }
      }

      // Temperature degrees of freedom
      if(solveForTemperature){

        // Perturb a temperature value and compute the flux divergence.
        double oldTemperature = temperature[perturbID];

        if(finiteDifferenceScheme == CENTRAL_DIFFERENCE){
          // Compute and store the negatively perturbed flux divergence.
          temperature[perturbID] -= epsilon;
          computeFluxDivergence(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
          temperature[perturbID] = oldTemperature;
          for(int i=0 ; i<fluxDivergenceVector->MyLength() ; ++i)
            tempFluxDivergence[i] = fluxDivergence[i];
        }

        // Compute the purturbed flux divergence.
        temperature[perturbID] += epsilon;
        computeFluxDivergence(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        temperature[perturbID] = oldTemperature;

        for(int i=0 ; i<numNeighbors+1 ; ++i){
          int fluxDivergenceID;
          if(i < numNeighbors)
            fluxDivergenceID = tempNeighborhoodList[i+1];
          else
            fluxDivergenceID = 0;

          double value = ( fluxDivergence[fluxDivergenceID] - tempFluxDivergence[fluxDivergenceID] ) / epsilon;
          if(finiteDifferenceScheme == CENTRAL_DIFFERENCE)
            value *= 0.5;
          scratchMatrix(numDof*fluxDivergenceID + temperatureDofOffset, numDof*perturbID + temperatureDofOffset) = value;
        }
      }
    }

    // Multiply by nodal volume
    for(unsigned int row=0 ; row<globalIndices.size() ; ++row){
      for(unsigned int col=0 ; col<globalIndices.size() ; ++col){
        scratchMatrix(row, col) *= volume[row/numDof];
      }
    }

    // Check for NaNs
    for(unsigned int row=0 ; row<globalIndices.size() ; ++row){
      for(unsigned int col=0 ; col<globalIndices.size() ; ++col){
        TestForTermination(!std::isfinite(scratchMatrix(row, col)), "**** NaN detected in finite-difference Jacobian.\n");
      }
    }

    // Sum the values into the global tangent matrix (this is expensive).
    if (jacobianType == PeridigmNS::Material::FULL_MATRIX)
      jacobian.addValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    else if (jacobianType == PeridigmNS::Material::BLOCK_DIAGONAL) {
      jacobian.addBlockDiagonalValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    }
    else // unknown jacobian type
      TestForTermination(true, "**** Unknown Jacobian Type\n");
  }
}

template<typename ScalarT>
ScalarT getData(const Teuchos::ParameterList & params,
  std::string paramName
)
{
  ScalarT parameter = 0;
  if (params.isParameter(paramName))
      parameter = params.get<ScalarT>(paramName);
  return parameter;
}

template double getData<double>
(const Teuchos::ParameterList & params,
  std::string paramName
);
template bool getData<bool>
(const Teuchos::ParameterList & params,
  std::string paramName
);
template int getData<int>
(const Teuchos::ParameterList & params,
  std::string paramName
);

double PeridigmNS::Material::calculateBulkModulus(const Teuchos::ParameterList & params) const
{
  bool bulkModulusDefined(false), shearModulusDefined(false), youngsModulusDefined(false), poissonsRatioDefined(false);
  double bulkModulus(0.0), shearModulus(0.0), youngsModulus(0.0), poissonsRatio(0.0);
  double computedValue = 0.0;

  if( params.isParameter("Bulk Modulus") ){
    bulkModulusDefined = true;
    bulkModulus = params.get<double>("Bulk Modulus");
  }
  if( params.isParameter("Shear Modulus") ){
    shearModulus = params.get<double>("Shear Modulus");
    shearModulusDefined = true;
  }
  if( params.isParameter("Young's Modulus") ){
    youngsModulus = params.get<double>("Young's Modulus");
    youngsModulusDefined = true;
  }
  if( params.isParameter("Poisson's Ratio") ){
    poissonsRatio = params.get<double>("Poisson's Ratio");
    poissonsRatioDefined = true;
  }

  int numDefinedConstants = static_cast<int>(bulkModulusDefined) +
    static_cast<int>(shearModulusDefined) +
    static_cast<int>(youngsModulusDefined) +
    static_cast<int>(poissonsRatioDefined);

  TestForTermination(numDefinedConstants != 2, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  if(bulkModulusDefined)
    computedValue = bulkModulus;
  else if(youngsModulusDefined && shearModulusDefined)
    computedValue = (youngsModulus * shearModulus) / (3.0*(3.0*shearModulus - youngsModulus));
  else if(youngsModulusDefined && poissonsRatioDefined)
    computedValue = youngsModulus / (3.0*(1.0 - 2.0*poissonsRatio));
  else if(shearModulusDefined && poissonsRatioDefined)
    computedValue = (2.0*shearModulus*(1.0 + poissonsRatio)) / (3.0*(1.0 - 2.0*poissonsRatio));
  else
    TestForTermination(true, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  return computedValue;
}


double PeridigmNS::Material::calculateShearModulus(const Teuchos::ParameterList & params) const
{
  bool bulkModulusDefined(false), shearModulusDefined(false), youngsModulusDefined(false), poissonsRatioDefined(false);
  double bulkModulus(0.0), shearModulus(0.0), youngsModulus(0.0), poissonsRatio(0.0);
  double computedValue;

  if( params.isParameter("Bulk Modulus") ){
    bulkModulusDefined = true;
    bulkModulus = params.get<double>("Bulk Modulus");
  }
  if( params.isParameter("Shear Modulus") ){
    shearModulus = params.get<double>("Shear Modulus");
    shearModulusDefined = true;
  }
  if( params.isParameter("Young's Modulus") ){
    youngsModulus = params.get<double>("Young's Modulus");
    youngsModulusDefined = true;
  }
  if( params.isParameter("Poisson's Ratio") ){
    poissonsRatio = params.get<double>("Poisson's Ratio");
    poissonsRatioDefined = true;
  }

  int numDefinedConstants = static_cast<int>(bulkModulusDefined) +
    static_cast<int>(shearModulusDefined) +
    static_cast<int>(youngsModulusDefined) +
    static_cast<int>(poissonsRatioDefined);

  TestForTermination(numDefinedConstants != 2, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  if(shearModulusDefined)
    computedValue = shearModulus;
  else if(bulkModulusDefined && youngsModulusDefined)
    computedValue = (3.0*bulkModulus*youngsModulus) / (9.0*bulkModulus - youngsModulus);
  else if(bulkModulusDefined & poissonsRatioDefined)
    computedValue = (3.0*bulkModulus*(1.0 - 2.0*poissonsRatio)) / (2.0*(1.0 + poissonsRatio));
  else if(youngsModulusDefined && poissonsRatioDefined)
    computedValue = youngsModulus / (2.0*(1.0 + poissonsRatio));
  else
    TestForTermination(true, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  return computedValue;
}

bool PeridigmNS::Material::getThermalExpansionCoefficient(const Teuchos::ParameterList & params, double m_alpha[][3]) const
{
  bool m_applyThermalStrains = false;

  if (params.isParameter("Apply Thermal Strain")){
    m_applyThermalStrains = params.get<bool>("Apply Thermal Strain"); 
    if (params.isParameter("Thermal Expansion Coefficient"))
    {
      m_alpha[0][0] = params.get<double>("Thermal Expansion Coefficient");
      m_alpha[1][1] = m_alpha[0][0];
      m_alpha[2][2] = m_alpha[0][0];
      m_alpha[1][0] = 0.0;
      m_alpha[1][2] = 0.0;
      m_alpha[2][0] = 0.0;
      m_alpha[2][1] = 0.0;
      m_alpha[0][1] = 0.0;
      m_alpha[0][2] = 0.0; 
      if (params.isParameter("Thermal Expansion Coefficient Y"))
        m_alpha[1][1] = params.get<double>("Thermal Expansion Coefficient Y");
      if (params.isParameter("Thermal Expansion Coefficient Z"))
        m_alpha[2][2] = params.get<double>("Thermal Expansion Coefficient Z");
    }
  }  
  
  return m_applyThermalStrains;
}

void PeridigmNS::Material::getStiffnessmatrix(const Teuchos::ParameterList& params, double C[][6], bool pstrain, bool pstress) const
{
  double ctemp[6][6];
          
  for(int iID=0 ; iID<6 ; ++iID)
    {
    for(int jID=0 ; jID<6 ; ++jID){
      C[iID][jID]=0.0;
      ctemp[iID][jID]=0.0;
    }
  }

  // bool iso = false;
  bool iso = true;

  if (params.isParameter("Material Symmetry")){
    if (params.get<string>("Material Symmetry")=="Anisotropic"){
            iso = false;
    }
    for(int iID=0 ; iID<6 ; ++iID)
    { // because C is symmetric only the upper triangle is read
      for(int jID=iID ; jID<6 ; ++jID){
        if (params.isParameter("C" + std::to_string(iID+1) + std::to_string(jID+1))){
          ctemp[iID][jID] = params.get<double>("C" + std::to_string(iID+1) + std::to_string(jID+1));
        }
        else
        {
          ctemp[iID][jID] = 0.0;
        }
      }
    }
  }
  else
  {
    iso = true;
    double m_bulkModulus = calculateBulkModulus(params);
    double m_shearModulus = calculateShearModulus(params);
    
    double EMod = 9*m_bulkModulus * m_shearModulus / (3*m_bulkModulus + m_shearModulus);
    double nu = (3*m_bulkModulus - 2*m_shearModulus) / (6*m_bulkModulus + 2*m_shearModulus);
 
    for(int iID=0 ; iID<6 ; ++iID)
    {
      for(int jID=0 ; jID<6 ; ++jID){
        ctemp[iID][jID] = 0.0;
      }
    }
    double temp = 1 / ((1+nu) * (1-2*nu));
    ctemp[0][0] = EMod  * (1-nu) * temp;
    ctemp[1][1] = EMod  * (1-nu) * temp;
    ctemp[2][2] = EMod  * (1-nu) * temp;
    ctemp[0][1] = EMod * nu * temp;
    ctemp[0][2] = EMod * nu * temp;
    ctemp[1][0] = EMod * nu * temp;
    ctemp[2][0] = EMod * nu * temp;
    ctemp[2][1] = EMod * nu * temp;
    ctemp[1][2] = EMod * nu * temp;

    ctemp[3][3] = m_shearModulus;
    ctemp[4][4] = m_shearModulus;
    ctemp[5][5] = m_shearModulus;
    


  }
    // Equation (8) Dipasquale, D., Sarego, G., Zaccariotto, M., Galvanetto, U., A discussion on failure criteria
    // for ordinary state-based Peridynamics, Engineering Fracture Mechanics (2017), doi: https://doi.org/10.1016/
    // j.engfracmech.2017.10.011
  
  if (pstress==false&&pstrain == false){
    for(int iID=0 ; iID<6 ; ++iID)
      {
      for(int jID=iID ; jID<6 ; ++jID){
        C[iID][jID] = ctemp[iID][jID];
        C[jID][iID] = ctemp[iID][jID];
      }
    }
  }
  else{
    if (pstress==true&&iso == true){
      //only transversal isotropic in the moment --> definition of iso missing
      double Eps = ctemp[5][5]*(3.0*ctemp[0][0] - 4.0*ctemp[5][5])/(ctemp[0][0] - ctemp[5][5]);     
      double nups = (Eps/2.0 - ctemp[5][5])/ctemp[5][5];
      //https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_plane_stress.cfm
      C[0][0] = Eps / (1-nups*nups);
      C[0][1] = Eps*nups / (1-nups*nups);
      C[1][0] = Eps*nups / (1-nups*nups);
      C[1][1] = Eps / (1-nups*nups);
      C[5][5] = ctemp[5][5];
    }
    else{
      C[0][0] = ctemp[0][0];C[0][1] = ctemp[0][1];C[0][5] = ctemp[0][5];
      C[1][0] = ctemp[0][1];C[1][1] = ctemp[1][1];C[1][5] = ctemp[1][5];
      C[5][0] = ctemp[0][5];ctemp[5][1] = ctemp[1][5];C[5][5] = ctemp[5][5];
    }
  }
}










      

void
PeridigmNS::Material::computeDeviatoricStateNorm(
        const int numOwnedPoints,
        const int* ownedIDs,
        const int* neighborhoodList,
        const double BM,
        const double SM,
        PeridigmNS::DataManager& dataManager) const {
            

    //////////////////////////////////////////////////////////////////////////////////////    
    // Routine developed for Peridigm by DLR.
    // Dilatation has to be updated before the damage routine. Therefore, this routine has been added
    // to the model evaluator. The dilatation is needed to calculate the energy at both bonds (12 and 21) 
    // for the correct time step, prior to the calculation of the internal forces. The bonds with critical
    // energy are deleted and the internal loads calculation will be done.
    //////////////////////////////////////////////////////////////////////////////////////    
    // Questions?
    // Christian Willberg  christian.willberg@dlr.de
    // Martin Raedel       martin.raedel@dlr.de
    //////////////////////////////////////////////////////////////////////////////////////   
     
    PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
    //double m_horizon = params.get<double>("Horizon");
    //////////////////////////////////////////////////////////////////////////////////////
    // ALL IDS will be initialized in the material routine "Peridigm_ElasticMaterial.cpp"
    // carefull if changed from Peridigm_ElasticMaterial.cpp to something else
    //////////////////////////////////////////////////////////////////////////////////////    
    
    int m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
    int m_coordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
    int m_volumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
    int m_weightedVolumeFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
    int m_bondDamageFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
    int m_horizonFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
    int m_damageModelFieldId = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Damage_Model_Data");
    int m_damagePlasticModelFieldId = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Damage_Plastic_Model_Data");
    int m_deviatoricPlasticExtensionFieldId = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deviatoric_Plastic_Extension");
  
    ///////////////////////////////////////////////////////////////////////////////
    PeridigmNS::InfluenceFunction::functionPointer scalarInfluenceFunction;
    scalarInfluenceFunction = PeridigmNS::InfluenceFunction::self().getInfluenceFunction();
    double omega;
    double *xOwned, *yOwned, *bondDamage, *horizon;
    double *Volume, *weightedVolume, *damageModel, *damagePlasticModel, *deviatoricPlasticExtensionState;
    ///////////////////////////////////////////////////////////////////////////////
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOwned);
    dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
    dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&Volume);
    dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&yOwned);   
    dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
    dataManager.getData(m_deviatoricPlasticExtensionFieldId, PeridigmField::STEP_N)->ExtractView(&deviatoricPlasticExtensionState);
    dataManager.getData(m_damageModelFieldId, PeridigmField::STEP_NP1)->ExtractView(&damageModel);
    dataManager.getData(m_damagePlasticModelFieldId, PeridigmField::STEP_NP1)->ExtractView(&damagePlasticModel);
    
    *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)) = *(dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N));
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
    //double m_alpha = 0.0; // to be adapted for temperature
    double *deltaTemperature = NULL; 
    double alpha;
    const double *deltaT = deltaTemperature;
    const double *m = weightedVolume;
    const double *v = Volume;
    double cellVolume, dx_X, dy_X, dz_X, zeta, edpN;
    const int *neighPtr = neighborhoodList;
    double theta, gamma, c, Kdm, MUdm;
    double dx_Y, dy_Y, dz_Y, dY, ed, tdTrial, ti;
    
    for(int p=0; p<numOwnedPoints;p++, xOwned+=3, yOwned+=3, horizon++, deltaT++, m++, v++){
        int numNeigh = *neighPtr; neighPtr++;
        double *X = xOwned;
        double *Y = yOwned;
        double norm=0.0;
        gamma = damagePlasticModel[3*p+2];
        theta = damageModel[3*p];
        Kdm = damageModel[3*p+1];
        MUdm = damageModel[3*p+2];
        if (gamma == 1) {
            alpha = 15.0*Kdm;
            c = 3 * Kdm * (theta);
        }
        else if (gamma == 2/3){
            alpha = 8.0*Kdm;
            c = (12.0*Kdm-4.0*MUdm) / 9.0 * theta;
        }
        else{
            alpha = 8.0*damageModel[3*p+1];
            c = 4.0*Kdm*MUdm/(3.0*Kdm+4.0*MUdm) * theta;
        }
       
        damagePlasticModel[3*p] = 0.0;
        for(int n=0;n<numNeigh;n++, neighPtr++, bondDamage++, deviatoricPlasticExtensionState++){
            int localId = *neighPtr;
            cellVolume = v[localId];
            double *XP = &X[3*localId];
            double *YP = &Y[3*localId];
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
             * gamma decides if it is a 2D or 3D case
             */
            ed = dY- zeta - gamma * theta*zeta/3;

            /*
             * Deviatoric plastic extension state from last step
             */
            edpN = *deviatoricPlasticExtensionState;
            omega = scalarInfluenceFunction(zeta,*horizon);
                
            /*
             * Compute trial stress
             */

            
            tdTrial = (1-*bondDamage) * omega * alpha * (ed - edpN);
            ti      = (1-*bondDamage) * omega * c * zeta;
            /*
             * Accumulate norm 
             */ 
            norm += tdTrial * tdTrial * cellVolume;
            // peridynamic pressure Eq. (86) Silling 2007
           
            damagePlasticModel[3*p] += (tdTrial + ti) * zeta / 3;
        }

        
        damagePlasticModel[3*p+1] = sqrt(norm);
        
    }
}        
