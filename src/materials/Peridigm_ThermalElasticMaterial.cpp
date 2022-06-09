/*! \file Peridigm_ThermalElasticMaterial.cpp */

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

#include "Peridigm_ThermalElasticMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic.h"
#include "nonlocal_thermal_diffusion.h"
#ifdef PERIDIGM_KOKKOS
	#include "elastic_kokkos.h"
#endif
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>

using namespace std;

PeridigmNS::ThermalElasticMaterial::ThermalElasticMaterial(const Teuchos::ParameterList& params)
  : ThermalMaterial(params),
	m_bulkModulus(0.0), m_shearModulus(0.0), m_density(0.0), m_alpha(0.0), m_horizon(0.0),m_thermalConductivity(0.0),m_specificHeat(0.0),m_convectionConstant(0.0),m_fluidTemperature(0.0),
	m_applyAutomaticDifferentiationJacobian(true),
	m_applySurfaceCorrectionFactor(false),
	m_applyThermalStrains(true),
	m_thermalShock(false),
	m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
	m_volumeFieldId(-1), m_damageFieldId(-1), m_weightedVolumeFieldId(-1), m_dilatationFieldId(-1), m_modelCoordinatesFieldId(-1),
	m_coordinatesFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1), m_surfaceCorrectionFactorFieldId(-1),
	m_heatFlowFieldId(-1),m_deltaTemperatureFieldId(-1)
{
  //! TODO Add meaningful asserts on material properties.
	m_bulkModulus = calculateBulkModulus(params);
	m_shearModulus = calculateShearModulus(params);
	m_density = params.get<double>("Density");
	m_horizon = params.get<double>("Horizon");
	m_alpha = params.get<double>("Thermal Expansion Coefficient");
	m_thermalConductivity = params.get<double>("Thermal Conductivity");
	m_specificHeat = params.get<double>("Specific Heat");

	if(params.isParameter("Apply Automatic Differentiation Jacobian"))
		m_applyAutomaticDifferentiationJacobian = params.get<bool>("Apply Automatic Differentiation Jacobian");

	if(params.isParameter("Apply Shear Correction Factor"))
		m_applySurfaceCorrectionFactor = params.get<bool>("Apply Shear Correction Factor");

	if(params.isParameter("Fluid Temperature"))
		m_thermalShock = true;

	if(m_thermalShock)
	{
		m_convectionConstant = params.get<double>("Convection Constant");
		m_fluidTemperature = params.get<double>("Fluid Temperature"); //Suppose that the specimen is quenched into a fluid kept at constant temperature (heat sink).
	}

	materialProperties["Density"] = m_density;
	materialProperties["Horizon"] = m_horizon;
	materialProperties["Thermal Expansion Coefficient"] = m_alpha;
	materialProperties["Thermal Conductivity"] = m_thermalConductivity;
	materialProperties["Specific Heat"] = m_specificHeat;
	if(m_thermalShock)
	{
		materialProperties["Convection Constant"] = m_convectionConstant;
		materialProperties["Fluid Temperature"] = m_fluidTemperature;
	}

	PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
	m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
	m_damageFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
	m_weightedVolumeFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
	m_dilatationFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation");
	m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
	m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
	m_forceDensityFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
	m_bondDamageFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
	m_surfaceCorrectionFactorFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Surface_Correction_Factor");
	m_heatFlowFieldId        		 = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Heat_Flow");
	m_deltaTemperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");

	m_fieldIds.push_back(m_volumeFieldId);
	m_fieldIds.push_back(m_damageFieldId);
	m_fieldIds.push_back(m_weightedVolumeFieldId);
	m_fieldIds.push_back(m_dilatationFieldId);
	m_fieldIds.push_back(m_modelCoordinatesFieldId);
	m_fieldIds.push_back(m_coordinatesFieldId);
	m_fieldIds.push_back(m_forceDensityFieldId);
	m_fieldIds.push_back(m_bondDamageFieldId);
	m_fieldIds.push_back(m_surfaceCorrectionFactorFieldId);
	m_fieldIds.push_back(m_heatFlowFieldId);
	m_fieldIds.push_back(m_deltaTemperatureFieldId);
}

PeridigmNS::ThermalElasticMaterial::~ThermalElasticMaterial()
{
}

double
PeridigmNS::ThermalElasticMaterial::lookupMaterialProperty(const std::string keyname) const
{
	std::map<std::string, double>::const_iterator search = materialProperties.find(keyname);
	if(search != materialProperties.end()){
		return search->second;
	}
	else
		TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error: requested material property is not in Thermal Elastic Material");
	// This is a fallthrough case to make the compiler happy.
	return 0.0;
}

void
PeridigmNS::ThermalElasticMaterial::initialize(const double dt,
											const int numOwnedPoints,
											const int* ownedIDs,
											const int* neighborhoodList,
											PeridigmNS::DataManager& dataManager)
{
	// Extract pointers to the underlying data
	double *xOverlap,  *cellVolumeOverlap, *weightedVolume;
	dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
	dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolumeOverlap);
	dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);

	MATERIAL_EVALUATION::computeWeightedVolume(xOverlap,cellVolumeOverlap,weightedVolume,numOwnedPoints,neighborhoodList,m_horizon);

	dataManager.getData(m_surfaceCorrectionFactorFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
	if(m_applySurfaceCorrectionFactor){
		Epetra_Vector temp(*dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1));
		int lengthYOverlap = temp.MyLength();
		double  *yOverlap,  *surfaceCorrectionFactor;
		temp.ExtractView(&yOverlap);
		dataManager.getData(m_surfaceCorrectionFactorFieldId, PeridigmField::STEP_NONE)->ExtractView(&surfaceCorrectionFactor);
		MATERIAL_EVALUATION::computeShearCorrectionFactor(numOwnedPoints,lengthYOverlap,xOverlap,yOverlap,cellVolumeOverlap,weightedVolume,neighborhoodList,m_horizon,surfaceCorrectionFactor);
	}
}

void
PeridigmNS::ThermalElasticMaterial::computeForce(	const double dt,
													const int numOwnedPoints,
													const int* ownedIDs,
													const int* neighborhoodList,
													PeridigmNS::DataManager& dataManager) const
{
	// Zero out the forces
	dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
// 	dataManager.getData(m_heatFlowFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

	// Extract pointers to the underlying data
	double *x, *y, *cellVolume, *weightedVolume, *dilatation, *bondDamage, *force, *scf, *deltaTemperature; //*heatFlow;

	dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
	dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
	dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
	dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
	dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatation);
	dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
	dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);
	dataManager.getData(m_surfaceCorrectionFactorFieldId, PeridigmField::STEP_NONE)->ExtractView(&scf);
// 	dataManager.getData(m_heatFlowFieldId, PeridigmField::STEP_NP1)->ExtractView(&heatFlow);
	dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
	MATERIAL_EVALUATION::computeDilatation(x,y,weightedVolume,cellVolume,bondDamage,dilatation,neighborhoodList,numOwnedPoints,m_horizon,m_OMEGA,m_alpha,deltaTemperature);
//#ifdef PERIDIGM_KOKKOS
//  MATERIAL_EVALUATION::computeInternalForceLinearElasticKokkos(x,y,weightedVolume,cellVolume,dilatation,bondDamage,scf,force,neighborhoodList,numOwnedPoints,m_bulkModulus,m_shearModulus,m_horizon,m_alpha,deltaTemperature);
	//#else
	MATERIAL_EVALUATION::computeInternalForceLinearElasticCoupled(	x,
																	y,
																	weightedVolume,
																	cellVolume,
																	dilatation,
																	bondDamage,
																	scf,
																	force,
																	neighborhoodList,
																	numOwnedPoints,
																	m_bulkModulus,
																	m_shearModulus,
																	m_horizon,
																	m_alpha,
																	deltaTemperature);

// 	MATERIAL_EVALUATION::computeHeatFlow(	x,
// 											y,
// 											cellVolume,
// 											bondDamage,
// 											heatFlow,
// 											neighborhoodList,
// 											numOwnedPoints,
// 											m_thermalConductivity,
// 											m_specificHeat,
// 											m_horizon,
// 											deltaTemperature);
 }

void
PeridigmNS::ThermalElasticMaterial::computeHeatFlow(	const double dt,
															const int numOwnedPoints,
															const int* ownedIDs,
															const int* neighborhoodList,
															PeridigmNS::DataManager& dataManager) const
{

// 			Zero out the heat flow
	dataManager.getData(m_heatFlowFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

// 	Extract pointers to the underlying data
	double *x, *y, *cellVolume, *bondDamage, *heatFlow, *deltaTemperature;

	dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
	dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
	dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
	dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
	dataManager.getData(m_heatFlowFieldId, PeridigmField::STEP_NP1)->ExtractView(&heatFlow);
	dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);

	MATERIAL_EVALUATION::computeHeatFlow(	x,
											y,
											cellVolume,
											bondDamage,
											heatFlow,
											neighborhoodList,
											numOwnedPoints,
											m_thermalConductivity,
											m_specificHeat,
											m_horizon,
											deltaTemperature);
}

void
PeridigmNS::ThermalElasticMaterial::computeStoredElasticEnergyDensity(
															const double dt,
															const int numOwnedPoints,
															const int* ownedIDs,
															const int* neighborhoodList,
															PeridigmNS::DataManager& dataManager) const
{
// 	This function is intended to be called from a compute class.
// 	The compute class should have already created the Stored_Elastic_Energy_Density field id.
	int storedElasticEnergyDensityFieldId = PeridigmNS::FieldManager::self().getFieldId("Stored_Elastic_Energy_Density");

	double *x, *y, *cellVolume, *weightedVolume, *dilatation, *storedElasticEnergyDensity, *bondDamage, *surfaceCorrectionFactor, *deltaTemperature;
	dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
	dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
	dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
	dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
	dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatation);
	dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
	dataManager.getData(storedElasticEnergyDensityFieldId, PeridigmField::STEP_NONE)->ExtractView(&storedElasticEnergyDensity);
	dataManager.getData(m_surfaceCorrectionFactorFieldId, PeridigmField::STEP_NONE)->ExtractView(&surfaceCorrectionFactor);
	dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);

	int iID, iNID, numNeighbors, nodeId, neighborId;
	double omega, nodeInitialX[3], nodeCurrentX[3];
	double initialDistance, currentDistance, deviatoricExtension, neighborBondDamage;
	double nodeDilatation, alpha, temp;

	int neighborhoodListIndex(0), bondIndex(0);
	for(iID=0 ; iID<numOwnedPoints ; ++iID){
		nodeId = ownedIDs[iID];
		nodeInitialX[0] = x[nodeId*3];
		nodeInitialX[1] = x[nodeId*3+1];
		nodeInitialX[2] = x[nodeId*3+2];
		nodeCurrentX[0] = y[nodeId*3];
		nodeCurrentX[1] = y[nodeId*3+1];
		nodeCurrentX[2] = y[nodeId*3+2];
		nodeDilatation = dilatation[nodeId];
		alpha = 15.0*m_shearModulus/weightedVolume[nodeId];
		alpha *= surfaceCorrectionFactor[nodeId];

		temp = 0.0;

		numNeighbors = neighborhoodList[neighborhoodListIndex++];
		for(iNID=0 ; iNID<numNeighbors ; ++iNID){
			neighborId = neighborhoodList[neighborhoodListIndex++];
			neighborBondDamage = bondDamage[bondIndex++];
			initialDistance =
			distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2],
					x[neighborId*3], x[neighborId*3+1], x[neighborId*3+2]);
			currentDistance =
			distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2],
					y[neighborId*3], y[neighborId*3+1], y[neighborId*3+2]) -
					m_alpha*deltaTemperature[nodeId]*initialDistance;
			deviatoricExtension = (currentDistance - initialDistance) - nodeDilatation*initialDistance/3.0;
			omega=m_OMEGA(initialDistance,m_horizon);
			temp += (1.0 - neighborBondDamage) * omega * deviatoricExtension * deviatoricExtension * cellVolume[neighborId];
		}
	storedElasticEnergyDensity[nodeId] = 0.5 * m_bulkModulus * nodeDilatation * nodeDilatation + 0.5 * alpha * temp;
	}
}

void
PeridigmNS::ThermalElasticMaterial::computeJacobian(const double dt,
											const int numOwnedPoints,
											const int* ownedIDs,
											const int* neighborhoodList,
											PeridigmNS::DataManager& dataManager,
											PeridigmNS::SerialMatrix& jacobian,
											PeridigmNS::Material::JacobianType jacobianType) const
{
  if(m_applyAutomaticDifferentiationJacobian){
    // Compute the Jacobian via automatic differentiation
    computeAutomaticDifferentiationJacobian(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian, jacobianType);
		//jacobian.getFECrsMatrix()->Print(std::cout);
  }
  else{
    // Call the base class function, which computes the Jacobian by finite difference
		TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error: Base class finite difference Jacobian is incompatible with thermal elastic material.\n");
    //PeridigmNS::Material::computeJacobian(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian, jacobianType);
  }
}


void
PeridigmNS::ThermalElasticMaterial::computeAutomaticDifferentiationJacobian(const double dt,
																			const int numOwnedPoints,
																			const int* ownedIDs,
																			const int* neighborhoodList,
																			PeridigmNS::DataManager& dataManager,
																			PeridigmNS::SerialMatrix& jacobian,
																			PeridigmNS::Material::JacobianType jacobianType) const
{
	// Compute contributions to the tangent matrix on an element-by-element basis
	// To reduce memory re-allocation, use static variable to store Fad types for
	// current coordinates (independent variables).
	static vector<Sacado::Fad::DFad<double> > y_AD;
	static vector<Sacado::Fad::DFad<double> > dTY_AD;

	// Loop over all points.
	int neighborhoodListIndex = 0;
	for(int iID=0 ; iID<numOwnedPoints ; ++iID){
		// Create a temporary neighborhood consisting of a single point and its neighbors.
		int numNeighbors = neighborhoodList[neighborhoodListIndex++];
		int numEntries = numNeighbors+1;
		int dofPerNode = 4;
		int numTotalNeighborhoodDof = dofPerNode*numEntries;
		vector<int> tempMyGlobalIDs(numEntries);
		// Put the node at the center of the neighborhood at the beginning of the list.
		tempMyGlobalIDs[0] = dataManager.getOwnedScalarPointMap()->GID(iID);
		vector<int> tempNeighborhoodList(numEntries);
		tempNeighborhoodList[0] = numNeighbors;
		for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
			int neighborID = neighborhoodList[neighborhoodListIndex++];
			tempMyGlobalIDs[iNID+1] = dataManager.getOverlapScalarPointMap()->GID(neighborID);
			tempNeighborhoodList[iNID+1] = iNID+1;
    }

	Epetra_SerialComm serialComm;
	Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 1, 0, serialComm));
	// This is not a mistake, the material is indeed defined in terms of one dimensional and three dimensional variables.
	Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 3, 0, serialComm));
	Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));

	// Create a temporary DataManager containing data for this point and its neighborhood.
	PeridigmNS::DataManager tempDataManager;
	tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
							tempOneDimensionalMap,
							Teuchos::RCP<const Epetra_BlockMap>(),
							tempThreeDimensionalMap,
							tempBondMap);

	// The temporary data manager will have the same field specs and data as the real data manager.
	vector<int> fieldIds = dataManager.getFieldIds();
	tempDataManager.allocateData(fieldIds);
	tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);

	// Set up numOwnedPoints and ownedIDs.
	// There is only one owned ID, and it has local ID zero in the tempDataManager.
	int tempNumOwnedPoints = 1;
	vector<int> tempOwnedIDs(tempNumOwnedPoints);
	tempOwnedIDs[0] = 0;

	// Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
	// Resize scratchMatrix if necessary
	if(scratchMatrix.Dimension() < dofPerNode*numEntries)
		scratchMatrix.Resize(dofPerNode*numEntries);

	// Create a list of global indices for the rows/columns in the scratch matrix.
	vector<int> globalIndices(numTotalNeighborhoodDof);
	for(int i=0 ; i<numEntries ; ++i){
		int globalID = tempOneDimensionalMap->GID(i);
		for(int j=0 ; j<dofPerNode ; ++j)
			globalIndices[dofPerNode*i+j] = dofPerNode*globalID+j;
	}

	// Extract pointers to the underlying data in the constitutiveData array.
	double *x, *y, *cellVolume, *weightedVolume, *damage, *bondDamage, *scf, *heatFlow, *deltaTemperature;
	tempDataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
	tempDataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
	tempDataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
	tempDataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
	tempDataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
	tempDataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
	tempDataManager.getData(m_surfaceCorrectionFactorFieldId, PeridigmField::STEP_NONE)->ExtractView(&scf);
	tempDataManager.getData(m_heatFlowFieldId, PeridigmField::STEP_NP1)->ExtractView(&heatFlow);
	tempDataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
	// Create arrays of Fad objects for the current coordinates, dilatation, and force density
	// current delta temperature and heat flux density
	// Modify the existing vector of Fad objects for the current coordinates
	if((int)y_AD.size() < (dofPerNode-1)*numEntries)
		y_AD.resize((dofPerNode-1)*numEntries);
	if((int)dTY_AD.size() < numEntries)
		dTY_AD.resize(numEntries);

		// We want to get derivatives with respect to y and deltaTemperature at the same time
		// so we must determine:
		// Out of the total columns which of these are
		// entries for solids and which are entries for temperature?
	for(int i=0 ; i<numTotalNeighborhoodDof ; i+=dofPerNode){
		// First three dof in a pack of dofPerNode are for solids
		for(int j=0 ; j<3 ; ++j){
			y_AD[i*3/dofPerNode+j].diff(i+j, numTotalNeighborhoodDof);
			// Convert index stride and store value
			y_AD[i*3/dofPerNode+j].val() = y[i*3/dofPerNode+j];
		}
		// Last dof in a pack of dofPerNode is always delta temperature
		dTY_AD[i/dofPerNode].diff(i+3,numTotalNeighborhoodDof);
		dTY_AD[i/dofPerNode].val() = deltaTemperature[i/dofPerNode];
	}
	// Create vectors of empty AD types for the dependent variables
	vector<Sacado::Fad::DFad<double> > dilatation_AD(numEntries);
	vector<Sacado::Fad::DFad<double> > force_AD((dofPerNode-1)*numEntries);
	vector<Sacado::Fad::DFad<double> > heatFlow_AD(numEntries);

	// Evaluate the constitutive model using the AD types
	MATERIAL_EVALUATION::computeDilatation(x,&y_AD[0],weightedVolume,cellVolume,bondDamage,&dilatation_AD[0],&tempNeighborhoodList[0],tempNumOwnedPoints,m_horizon,m_OMEGA,m_alpha,deltaTemperature);
	MATERIAL_EVALUATION::computeInternalForceLinearElasticCoupled(x,&y_AD[0],weightedVolume,cellVolume,&dilatation_AD[0],bondDamage,scf,&force_AD[0],&tempNeighborhoodList[0],tempNumOwnedPoints,m_bulkModulus,m_shearModulus,m_horizon,m_alpha,&dTY_AD[0]);

	MATERIAL_EVALUATION::computeHeatFlow(x,&y_AD[0],cellVolume,bondDamage,&heatFlow_AD[0],&tempNeighborhoodList[0],tempNumOwnedPoints,m_thermalConductivity,m_specificHeat,m_horizon,&dTY_AD[0]);

	// Load derivative values into scratch matrix
	// Multiply by volume along the way to convert force density to force
	double value;
	for(int row=0 ; row<numTotalNeighborhoodDof ; row+=dofPerNode){
		for(int col=0 ; col<numTotalNeighborhoodDof ; col+=dofPerNode){
				for(int subcol=0 ; subcol<dofPerNode ; ++subcol){
					for(int subrow=0 ; subrow<(dofPerNode-1) ; ++subrow){
							value = force_AD[row*3/dofPerNode + subrow].dx(col + subcol) * cellVolume[row/dofPerNode];
							TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(value), "**** NaN detected in ThermalElasticMaterial::computeAutomaticDifferentiationJacobian() (internal force).\n");
							scratchMatrix(row+subrow, col+subcol) = value;
					}
					value = heatFlow_AD[row/dofPerNode].dx(col + subcol) * cellVolume[row/dofPerNode];
					TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(value), "**** NaN detected in ThermalElasticMaterial::computeAutomaticDifferentiationJacobian() (heat flow).\n");
			scratchMatrix(row +dofPerNode -1, col+subcol) = value;
				}
			}
		}

	// Sum the values into the global tangent matrix (this is expensive).
	if (jacobianType == PeridigmNS::Material::FULL_MATRIX)
		jacobian.addValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
	else if (jacobianType == PeridigmNS::Material::BLOCK_DIAGONAL) {
		jacobian.addBlockDiagonalValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
	}
	else // unknown jacobian type
		TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Unknown Jacobian Type\n");
	}
	}


/*
void
PeridigmNS::ThermalElasticMaterial::computeAutomaticDifferentiationJacobian(const double dt,
                                                                     const int numOwnedPoints,
                                                                     const int* ownedIDs,
                                                                     const int* neighborhoodList,
                                                                     PeridigmNS::DataManager& dataManager,
                                                                     PeridigmNS::SerialMatrix& jacobian,
                                                                     PeridigmNS::Material::JacobianType jacobianType) const
{
  // Compute contributions to the tangent matrix on an element-by-element basis
  // To reduce memory re-allocation, use static variable to store Fad types for
  // current coordinates (independent variables).
  static vector<Sacado::Fad::DFad<double> > y_AD;
  // Loop over all points.
  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    // Create a temporary neighborhood consisting of a single point and its neighbors.
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    int numEntries = numNeighbors+1;
    int numDof = 3*numEntries;
    vector<int> tempMyGlobalIDs(numEntries);
    // Put the node at the center of the neighborhood at the beginning of the list.
    tempMyGlobalIDs[0] = dataManager.getOwnedScalarPointMap()->GID(iID);
    vector<int> tempNeighborhoodList(numEntries);
    tempNeighborhoodList[0] = numNeighbors;
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      tempMyGlobalIDs[iNID+1] = dataManager.getOverlapScalarPointMap()->GID(neighborID);
      tempNeighborhoodList[iNID+1] = iNID+1;
    }
    Epetra_SerialComm serialComm;
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 3, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));
    // Create a temporary DataManager containing data for this point and its neighborhood.
    PeridigmNS::DataManager tempDataManager;
    tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
                            tempOneDimensionalMap,
                            Teuchos::RCP<const Epetra_BlockMap>(),
                            tempThreeDimensionalMap,
                            tempBondMap);
    // The temporary data manager will have the same field specs and data as the real data manager.
    vector<int> fieldIds = dataManager.getFieldIds();
    tempDataManager.allocateData(fieldIds);
    tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);
    // Set up numOwnedPoints and ownedIDs.
    // There is only one owned ID, and it has local ID zero in the tempDataManager.
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(tempNumOwnedPoints);
    tempOwnedIDs[0] = 0;
    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
    // Resize scratchMatrix if necessary
    if(scratchMatrix.Dimension() < numDof)
      scratchMatrix.Resize(numEntries*4);
    // Create a list of global indices for the rows/columns in the scratch matrix.
    vector<int> globalIndices(numEntries*4);
    for(int i=0 ; i<numEntries ; ++i){
      int globalID = tempOneDimensionalMap->GID(i);
      for(int j=0 ; j<4 ; ++j)
        globalIndices[4*i+j] = 4*globalID+j;
    }
    // Extract pointers to the underlying data in the constitutiveData array.
    double *x, *y, *cellVolume, *weightedVolume, *damage, *bondDamage, *scf, *deltaTemperature;
    tempDataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
    tempDataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
    tempDataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
    tempDataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
    tempDataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
    tempDataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
    tempDataManager.getData(m_surfaceCorrectionFactorFieldId, PeridigmField::STEP_NONE)->ExtractView(&scf);
    deltaTemperature = NULL;
    if(m_applyThermalStrains)
      tempDataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
    // Create arrays of Fad objects for the current coordinates, dilatation, and force density
    // Modify the existing vector of Fad objects for the current coordinates
    if((int)y_AD.size() < numDof)
      y_AD.resize(numDof);
    for(int i=0 ; i<numDof ; ++i){
      y_AD[i].diff(i, numDof);
      y_AD[i].val() = y[i];
    }
    // Create vectors of empty AD types for the dependent variables
    vector<Sacado::Fad::DFad<double> > dilatation_AD(numEntries);
    vector<Sacado::Fad::DFad<double> > force_AD(numDof);
    // Evaluate the constitutive model using the AD types
    MATERIAL_EVALUATION::computeDilatation(x,&y_AD[0],weightedVolume,cellVolume,bondDamage,&dilatation_AD[0],&tempNeighborhoodList[0],tempNumOwnedPoints,m_horizon,m_OMEGA,m_alpha,deltaTemperature);
    MATERIAL_EVALUATION::computeInternalForceLinearElastic(x,&y_AD[0],weightedVolume,cellVolume,&dilatation_AD[0],bondDamage,scf,&force_AD[0],&tempNeighborhoodList[0],tempNumOwnedPoints,m_bulkModulus,m_shearModulus,m_horizon,m_alpha,deltaTemperature);
    // Load derivative values into scratch matrix
    // Multiply by volume along the way to convert force density to force
    double value;
    for(int row=0 ; row<numEntries*4 ; row+=4){
    	for(int col=0 ; col<numEntries*4 ; col+=4){
				for(int subcol=0 ; subcol<4 ; ++subcol){
	  			for(int subrow=0 ; subrow<4 ; ++subrow){
						if(not (subcol == 3 or subrow == 3)){
	     				value = force_AD[row*3/4 + subrow].dx(col*3/4 + subcol) * cellVolume[row/4];
	     				TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite(value), "**** NaN detected in ThermalElasticMaterial::computeAutomaticDifferentiationJacobian() (internal force).\n");
						}
						else if(subcol = subrow)
							value = 1.0;
						else
							value = 0.0;
	     			scratchMatrix(row+subrow, col+subcol) = value;
	  			}
			}
     }
    }
    // Sum the values into the global tangent matrix (this is expensive).
    if (jacobianType == PeridigmNS::Material::FULL_MATRIX)
      jacobian.addValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    else if (jacobianType == PeridigmNS::Material::BLOCK_DIAGONAL) {
      jacobian.addBlockDiagonalValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
    }
    else // unknown jacobian type
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Unknown Jacobian Type\n");
  }
}
*/