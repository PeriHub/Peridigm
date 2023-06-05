/*! \file Peridigm_TextFileDiscretization.cpp */

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

#include "Peridigm_TextFileDiscretization.hpp"
#include "Peridigm_HorizonManager.hpp"
#include "Peridigm_Enums.hpp"
#include "NeighborhoodList.h"
#include "PdZoltan.h"
#include "Peridigm_Logging.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Peridigm_Logging.hpp>

#include <sstream>
#include <fstream>

using namespace std;

PeridigmNS::TextFileDiscretization::TextFileDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                                           const Teuchos::RCP<Teuchos::ParameterList>& params) :
  Discretization(epetra_comm),
  minElementRadius(1.0e50),
  maxElementRadius(0.0),
  maxElementDimension(0.0),
  numBonds(0),
  maxNumBondsPerElem(0),
  myPID(epetra_comm->MyPID()),
  numPID(epetra_comm->NumProc()),
  bondFilterCommand("None")
{
  TestForTermination(params->get<string>("Type") != "Text File", "Invalid Type in TextFileDiscretization");

  string meshFileName = params->get<string>("Input Mesh File");
  string topologyFileName = "";

  if (params->isParameter("Input FEM Topology File"))
    topologyFileName = params->get<string>("Input FEM Topology File");
  if (params->isParameter("Omit Bonds Between Blocks"))
    bondFilterCommand = params->get<string>("Omit Bonds Between Blocks");

  // Set up bond filters
  createBondFilters(params);

  QUICKGRID::Data decomp = getDecomp(meshFileName, topologyFileName, params);

  // \todo Refactor; the createMaps() call is currently inside getDecomp() due to order-of-operations issues with tracking element blocks.
  //createMaps(decomp);
  createNeighborhoodData(decomp);

  // 3D only
  TestForTermination(decomp.dimension != 3, "Invalid dimension in decomposition (only 3D is supported)");

  // fill the x vector with the current positions (owned positions only)
  initialX = Teuchos::rcp(new Epetra_Vector(Copy, *threeDimensionalMap, decomp.myX.get()));
  // fill with local coordinate distribution
  pointAngle = Teuchos::rcp(new Epetra_Vector(Copy, *threeDimensionalMap, decomp.myAngle.get()));
  // fill with point or element separation
  nodeType = Teuchos::rcp(new Epetra_Vector(Copy, *oneDimensionalMap, decomp.myNodeType.get()));
  // Create the bondMap, a local map used for constitutive data stored on bonds.
  createBondMapAndCheckForZeroNeighbors(bondMap, oneDimensionalMap, neighborhoodData, numBonds, maxNumBondsPerElem);

  // fill cell volumes
  cellVolume = Teuchos::rcp(new Epetra_Vector(Copy,*oneDimensionalMap,decomp.cellVolume.get()) ); 
   // fill point time
  pointTime = Teuchos::rcp(new Epetra_Vector(Copy,*oneDimensionalMap,decomp.myPointTime.get()) );

  // find the minimum element radius
  for(int i=0 ; i<cellVolume->MyLength() ; ++i){
    double radius = pow(0.238732414637843*(*cellVolume)[i], 0.33333333333333333);
    if(radius < minElementRadius)
      minElementRadius = radius;
    if(radius > maxElementRadius)
      maxElementRadius = radius;
  }
  vector<double> localMin(1);
  vector<double> globalMin(1);
  localMin[0] = minElementRadius;
  epetra_comm->MinAll(&localMin[0], &globalMin[0], 1);
  minElementRadius = globalMin[0];
  localMin[0] = maxElementRadius;
  epetra_comm->MaxAll(&localMin[0], &globalMin[0], 1);
  maxElementRadius = globalMin[0];
}

PeridigmNS::TextFileDiscretization::~TextFileDiscretization() {}

void PeridigmNS::TextFileDiscretization::getDiscretization(const string& textFileName,
                                                           vector<double> &coordinates,
                                                           vector<int> &blockIds,
                                                           vector<double> &volumes,
                                                           vector<double> &angles,
                                                           vector<double> &pointTime,
                                                           vector<double> &nodeType

)
{
  if (myPID == 0)
  {
    ifstream inFile(textFileName.c_str());
    TestForTermination(!inFile.is_open(), "**** Error opening discretization text file.\n");
    while(inFile.good()){
      string str;
      getline(inFile, str);
      str = trim(str);
      // Ignore comment lines, otherwise parse
      if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
        istringstream iss(str);
        vector<double> data;
        copy(istream_iterator<double>(iss),
             istream_iterator<double>(),
             back_inserter<vector<double> >(data));
        // Check for obvious problems with the data
        // Adapt to coordinate system and without --> to check
        if (data.size() != 5 && data.size() != 6 && data.size() != 8 && data.size() != 9)
        {
          string msg = "\n**** Error parsing text file, invalid line: " + str + "\n";
          TEUCHOS_TEST_FOR_EXCEPT_MSG(data.size() != 5, msg);
          TEUCHOS_TEST_FOR_EXCEPT_MSG(data.size() != 6, msg);
          TEUCHOS_TEST_FOR_EXCEPT_MSG(data.size() != 8, msg);
          TEUCHOS_TEST_FOR_EXCEPT_MSG(data.size() != 9, msg);
        }
        bool anglesImport = false;
        if (data.size() == 8) anglesImport = true;
        bool timeImport  = false;
        if (data.size() == 6) timeImport  = true;
        if (data.size() == 9){
          timeImport  = true;
          anglesImport = true;
        }
        // Store the coordinates, block id, volumes and angles
        coordinates.push_back(data[0]);
        coordinates.push_back(data[1]);
        coordinates.push_back(data[2]);
        nodeType.push_back(1);
        blockIds.push_back(static_cast<int>(data[3]));
        volumes.push_back(data[4]);
        if (timeImport && anglesImport){
          pointTime.push_back(data[5]);
          angles.push_back(data[6]);
          angles.push_back(data[7]);
          angles.push_back(data[8]);
        }else{
          if (anglesImport){
            angles.push_back(data[5]);
            angles.push_back(data[6]);
            angles.push_back(data[7]);
          }
          else{
            angles.push_back(0.0);
            angles.push_back(0.0);
            angles.push_back(0.0);
          }        
          if (timeImport){
            pointTime.push_back(data[5]);
          }
          else{
            pointTime.push_back(0.0);
          }
        }

      }
    }
    inFile.close();
  }
}

QUICKGRID::Data PeridigmNS::TextFileDiscretization::getDecomp(const string& textFileName,
                                                              const string& topologyFileName,
                                                              const Teuchos::RCP<Teuchos::ParameterList>& params)
{

  // Read data from the text file
  vector<double> coordinates;
  vector<double> volumes;
  vector<int> blockIds;
  vector<double> angles;  
  vector<double> pointTime;
  vector<double> horizon_of_element;
  vector<int> elementTopo;
  vector<int> nodeTypeTemporary;
  vector<double> nodeType;
  getDiscretization(textFileName, coordinates, blockIds, volumes, angles, pointTime, nodeType);
  int numFE = 0;
  if (params->isParameter("Input FEM Topology File"))
  {
    LOG(LogLevel::INFO,"Read FEM Topology");
    getFETopology(topologyFileName, coordinates, blockIds, volumes, angles, horizon_of_element, elementTopo, nodeType, numFE);
    int numNodes = static_cast<int>(blockIds.size());
    int etz = static_cast<int>(elementTopo.size());
    //https://stackoverflow.com/questions/51408632/mpi-bcast-c-stl-vector
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&numFE, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numNodes, 1 , MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&etz, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
    if (myPID != 0){
      elementTopo.resize(etz,0);
      horizon_of_element.resize(numNodes,0.0);
    }
    nodeTypeTemporary.resize(numNodes,0);
    if (myPID == 0){ 
      for(unsigned int i=0 ; i<nodeTypeTemporary.size() ; i++){
        nodeTypeTemporary[i] = nodeType[i];
      }
    } 
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(elementTopo.data(),  elementTopo.size(), MPI_INT, 0, MPI_COMM_WORLD);
    // nodeType is double, because int does not exist in the datamanager
    MPI_Bcast(nodeTypeTemporary.data(),  nodeTypeTemporary.size(), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(horizon_of_element.data(),  horizon_of_element.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::shared_ptr<PdBondFilter::BondFilter> bondFilter(new PdBondFilter::PreDefinedTopologyFilter(numNodes, numFE, elementTopo));
    bondFilters.push_back(bondFilter);
 
    
    LOG(LogLevel::INFO,"End Input FEM Topology");
    
  }
  
  int numElements = static_cast<int>(blockIds.size());
  TestForTermination(myPID == 0 && numElements < 1, "**** Error reading discretization text file, no data found.\n");

  // Record the block ids on the root processor
  set<int> uniqueBlockIds;
  if(myPID == 0){
    for(unsigned int i=0 ; i<blockIds.size() ; ++i)
      uniqueBlockIds.insert(blockIds[i]);
  }

  // Broadcast necessary data from root processor
  
  Teuchos::RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  
  int numGlobalElements;
  reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, 1, &numElements, &numGlobalElements);
  
  // Broadcast the unique block ids so that all processors are aware of the full block list
  // This is necessary because if a processor does not have any elements for a given block, it will be unaware the
  // given block exists, which causes problems downstream
  int numLocalUniqueBlockIds = static_cast<int>( uniqueBlockIds.size() );
  int numGlobalUniqueBlockIds;
  reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, 1, &numLocalUniqueBlockIds, &numGlobalUniqueBlockIds);
  vector<int> uniqueLocalBlockIds(numGlobalUniqueBlockIds, 0);
  int index = 0;
  for(set<int>::const_iterator it = uniqueBlockIds.begin() ; it != uniqueBlockIds.end() ; it++)
    uniqueLocalBlockIds[index++] = *it;
  vector<int> uniqueGlobalBlockIds(numGlobalUniqueBlockIds);  
  reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, numGlobalUniqueBlockIds, &uniqueLocalBlockIds[0], &uniqueGlobalBlockIds[0]);

  // Create list of global ids
  vector<int> globalIds(numElements);
  for(unsigned int i=0 ; i<globalIds.size() ; ++i)
    globalIds[i] = i;

  // Copy data into a decomp object
  int dimension = 3;
  QUICKGRID::Data decomp = QUICKGRID::allocatePdGridData(numElements, dimension);

  decomp.globalNumPoints = numGlobalElements;
  memcpy(decomp.myGlobalIDs.get(), &globalIds[0], numElements*sizeof(int)); 
  memcpy(decomp.cellVolume.get(), &volumes[0], numElements*sizeof(double)); 
  memcpy(decomp.myX.get(), &coordinates[0], 3*numElements*sizeof(double));
  memcpy(decomp.myAngle.get(), &angles[0], 3*numElements*sizeof(double));
  memcpy(decomp.myPointTime.get(), &pointTime[0], numElements*sizeof(double));
  // double, because the datamanage won't allow int
  memcpy(decomp.myNodeType.get(), &nodeType[0], numElements*sizeof(double));
  
  // Create a blockID vector in the current configuration
  // That is, the configuration prior to load balancing
  Epetra_BlockMap tempOneDimensionalMap(decomp.globalNumPoints,
                                        decomp.numPoints,
                                        decomp.myGlobalIDs.get(),
                                        1,
                                        0,
                                        *comm);
  Epetra_Vector tempBlockID(tempOneDimensionalMap);
  double* tempBlockIDPtr;
  tempBlockID.ExtractView(&tempBlockIDPtr);
  for(unsigned int i=0 ; i<blockIds.size() ; ++i)
    tempBlockIDPtr[i] = blockIds[i];

  // call the rebalance function on the current-configuration decomp
  
  decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);

  // create a (throw-away) one-dimensional owned map in the rebalanced configuration
  Epetra_BlockMap rebalancedMap(decomp.globalNumPoints, decomp.numPoints, decomp.myGlobalIDs.get(), 1, 0, *comm);

  // Create a (throw-away) blockID vector corresponding to the load balanced decomposition
  Epetra_Vector rebalancedBlockID(rebalancedMap);
  Epetra_Import rebalancedImporter(rebalancedBlockID.Map(), tempBlockID.Map());
  rebalancedBlockID.Import(tempBlockID, rebalancedImporter, Insert);

  // Initialize the element list for each block
  // Force blocks with no on-processor elements to have an entry in the elementBlocks map
  for(unsigned int i=0 ; i<uniqueGlobalBlockIds.size() ; i++){
    stringstream blockName;
    blockName << "block_" << uniqueGlobalBlockIds[i];
    (*elementBlocks)[blockName.str()] = std::vector<int>();
  }

  // Create the element list for each block
  for(int i=0 ; i<rebalancedBlockID.MyLength() ; ++i){
    stringstream blockName;
    blockName << "block_" << rebalancedBlockID[i];
    TestForTermination(elementBlocks->find(blockName.str()) == elementBlocks->end(),
                                "\n**** Error in TextFileDiscretization::getDecomp(), invalid block id.\n");
    int globalID = rebalancedBlockID.Map().GID(i);
    (*elementBlocks)[blockName.str()].push_back(globalID);
  }

  // Record the horizon for each point
  PeridigmNS::HorizonManager& horizonManager = PeridigmNS::HorizonManager::self();
  Teuchos::RCP<Epetra_Vector> rebalancedHorizonForEachPoint = Teuchos::rcp(new Epetra_Vector(rebalancedMap));
  double* rebalancedX = decomp.myX.get();
  for(map<string, vector<int> >::const_iterator it = elementBlocks->begin() ; it != elementBlocks->end() ; it++){
    const string& blockName = it->first;
    const vector<int>& globalIds = it->second;

    bool hasConstantHorizon = horizonManager.blockHasConstantHorizon(blockName);
    double constantHorizonValue(0.0);
    if(hasConstantHorizon)
      constantHorizonValue = horizonManager.getBlockConstantHorizonValue(blockName);
    unsigned int numGlobalIDs = globalIds.size();
    

    for(unsigned int i=0 ; i<numGlobalIDs ; ++i){

      int localId = rebalancedMap.LID(globalIds[i]);
      if(hasConstantHorizon){
        (*rebalancedHorizonForEachPoint)[localId] = constantHorizonValue;
        if (params->isParameter("Input FEM Topology File")){ 
          if (nodeTypeTemporary[globalIds[i]] == 2) (*rebalancedHorizonForEachPoint)[localId] = horizon_of_element[globalIds[i]];
        }

      }
      else{
        double x = rebalancedX[localId*3];
        double y = rebalancedX[localId*3 + 1];
        double z = rebalancedX[localId*3 + 2];
        double horizon = horizonManager.evaluateHorizon(blockName, x, y, z);
        (*rebalancedHorizonForEachPoint)[localId] = horizon;
        if (params->isParameter("Input FEM Topology File")){ 
          if (nodeTypeTemporary[globalIds[i]] == 2) (*rebalancedHorizonForEachPoint)[localId] = horizon_of_element[globalIds[i]];
        }
      }
    }
  }
  // execute neighbor search and update the decomp to include resulting ghosts
  std::shared_ptr<const Epetra_Comm> commSp(comm.getRawPtr(), NonDeleter<const Epetra_Comm>());
  Teuchos::RCP<PDNEIGH::NeighborhoodList> list;
  if(bondFilters.size() == 0){
    list = Teuchos::rcp(new PDNEIGH::NeighborhoodList(commSp,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,rebalancedHorizonForEachPoint));
  }
  else{
    list = Teuchos::rcp(new PDNEIGH::NeighborhoodList(commSp,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,rebalancedHorizonForEachPoint,bondFilters));
  }
  decomp.neighborhood=list->get_neighborhood();
  decomp.sizeNeighborhoodList=list->get_size_neighborhood_list();
  decomp.neighborhoodPtr=list->get_neighborhood_ptr();

  // Create all the maps.
  createMaps(decomp);

  // Create the blockID vector corresponding to the load balanced decomposition
  blockID = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  Epetra_Import tempImporter(blockID->Map(), tempBlockID.Map());
  blockID->Import(tempBlockID, tempImporter, Insert);

  // Create the horizonForEachPonit vector corresponding to the load balanced decomposition
  horizonForEachPoint = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  Epetra_Import horizonImporter(horizonForEachPoint->Map(), rebalancedHorizonForEachPoint->Map());
  horizonForEachPoint->Import(*rebalancedHorizonForEachPoint, horizonImporter, Insert);

  return decomp;
}

void PeridigmNS::TextFileDiscretization::getFETopology(const string& fileName,
                                                       vector<double>& coordinates,
                                                       vector<int>& blockIds,
                                                       vector<double>& volumes,
                                                       vector<double>& angles,
                                                       vector<double>& horizon,
                                                       vector<int>& elementTopo,
                                                       vector<double>& nodeType,
                                                       int& numFE)
{

  if (myPID == 0)
  {
    numFE = 0;
    double coorAvg[3], angAvg[3], volAvg;
    int numOfFiniteElements;
    std::string testString = "";
    for(unsigned int i=0 ; i<blockIds.size() ; i++) horizon.push_back(0.0);
    if (fileName.compare(testString) != 0)
    {
      ifstream inFile(fileName.c_str());
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!inFile.is_open(), "**** Error opening topology text file.\n");
      while (inFile.good())
      {
        // defines points as points
        // this is needed to separate between virtual element points and real points
 
        numOfFiniteElements += 1;
        string str;
        getline(inFile, str);
        str = trim(str);
        // Ignore comment lines, otherwise parse
        if (!(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0))
        {
          istringstream iss(str);
          vector<int> topo;
          copy(istream_iterator<int>(iss),
               istream_iterator<int>(),
               back_inserter<vector<int>>(topo));

          blockIds.push_back(topo[0]);      
          elementTopo.push_back(static_cast<int>(topo.size() - 1));
          // lenDecompElementNodes += static_cast<int>(topo.size()) + 1;

          /*
          tbd if multiple element types could exist in one block
          */
          coorAvg[0] = 0;
          coorAvg[1] = 0;
          coorAvg[2] = 0;
          angAvg[0] = 0;
          angAvg[1] = 0;
          angAvg[2] = 0;
          volAvg = 0;
          for (unsigned int n = 1; n < topo.size(); n++)
          {
            elementTopo.push_back(static_cast<int>(topo[n]));
            coorAvg[0] += coordinates[3 * topo[n]];
            coorAvg[1] += coordinates[3 * topo[n] + 1];
            coorAvg[2] += coordinates[3 * topo[n] + 2];
            angAvg[0] += angles[3 * topo[n]];
            angAvg[1] += angles[3 * topo[n] + 1];
            angAvg[2] += angles[3 * topo[n] + 2];
            volAvg += volumes[topo[n]];
          }

          for (unsigned int i = 0; i < 3; i++)
          {
            coorAvg[i] /= (topo.size() - 1);
            angAvg[i] /= (topo.size() - 1);
          }
          coordinates.push_back(coorAvg[0]);
          coordinates.push_back(coorAvg[1]);
          coordinates.push_back(coorAvg[2]);
          nodeType.push_back(2); // defines FE Elements
          angles.push_back(angAvg[0]);
          angles.push_back(angAvg[1]);
          angles.push_back(angAvg[2]);
          horizon.push_back(get_max_dist(coordinates, coorAvg, topo));
          volumes.push_back(volAvg / (topo.size() - 1));
          numFE = numFE + 1;
        }
      }

      inFile.close();
    }
    // bondfilteraufruf
  }
}
double PeridigmNS::TextFileDiscretization::get_max_dist(const vector<double> &coordinates, const double coorAvg[3],
                                                        const vector<int> &topo)
{
  double horizon = 0.0;
  double dist = 0.0;
  for (int n = 1; n < topo[0] + 1; n++)
  {
    dist = abs(distance(coorAvg[0], coorAvg[1], coorAvg[2],
                        coordinates[3 * topo[n]], coordinates[3 * topo[n] + 1], coordinates[3 * topo[n] + 2]));
    if (horizon < dist)
      // adding on percent to avoid round off errors
      // because the neighborhoodlist is later adapted it does not matter if
      // there are some extra nodes
      horizon = 1.01 * dist;
  }

  return horizon;
}

void
PeridigmNS::TextFileDiscretization::createMaps(const QUICKGRID::Data& decomp)
{
  int dimension;

  // oneDimensionalMap
  // used for global IDs and scalar data
  dimension = 1;
  oneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOwnedMap(*comm, decomp, dimension)));

  // oneDimensionalOverlapMap
  // used for global IDs and scalar data, includes ghosts
  dimension = 1;
  oneDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOverlapMap(*comm, decomp, dimension)));

  // threeDimensionalMap
  // used for R3 vector data, e.g., u, v, etc.
  dimension = 3;
  threeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOwnedMap(*comm, decomp, dimension)));

  // threeDimensionalOverlapMap
  // used for R3 vector data, e.g., u, v, etc.,  includes ghosts
  dimension = 3;
  threeDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOverlapMap(*comm, decomp, dimension)));
}

void
PeridigmNS::TextFileDiscretization::createNeighborhoodData(const QUICKGRID::Data& decomp)
{
   neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
   neighborhoodData->SetNumOwned(decomp.numPoints);
   memcpy(neighborhoodData->OwnedIDs(), 
 		 Discretization::getLocalOwnedIds(decomp, *oneDimensionalOverlapMap).get(),
 		 decomp.numPoints*sizeof(int));
   memcpy(neighborhoodData->NeighborhoodPtr(), 
 		 decomp.neighborhoodPtr.get(),
 		 decomp.numPoints*sizeof(int));
   neighborhoodData->SetNeighborhoodListSize(decomp.sizeNeighborhoodList);
   memcpy(neighborhoodData->NeighborhoodList(),
 		 Discretization::getLocalNeighborList(decomp, *oneDimensionalOverlapMap).get(),
 		 decomp.sizeNeighborhoodList*sizeof(int));
   neighborhoodData = filterBonds(neighborhoodData);
}

Teuchos::RCP<PeridigmNS::NeighborhoodData>
PeridigmNS::TextFileDiscretization::filterBonds(Teuchos::RCP<PeridigmNS::NeighborhoodData> unfilteredNeighborhoodData)
{
  // Set up a block bonding matrix, which defines whether or not bonds should be formed across blocks
  int numBlocks = getNumBlocks();
  std::vector< std::vector<bool> > blockBondingMatrix(numBlocks);
  for(int i=0 ; i<numBlocks ; ++i){
    blockBondingMatrix[i].resize(numBlocks, true);
  }

  if(bondFilterCommand == "None"){
    // All blocks are bonded, the blockBondingMatrix is unchanged
    return unfilteredNeighborhoodData;
  }
  else if(bondFilterCommand == "All"){
    // No blocks are bonded, the blockBondingMatrix is the identity matrix
    for(int i=0 ; i<numBlocks ; ++i){
      for(int j=0 ; j<numBlocks ; ++j){
        if(i != j)
          blockBondingMatrix[i][j] = false;
      }
    }
  }
  else{
    string msg = "**** Error, unrecognized value for \"Omit Bonds Between Blocks\":  ";
    msg += bondFilterCommand + "\n";
    msg += "**** Valid options are:  All, None\n";
    TestForTermination(true, msg);
  }

  // Create an overlap vector containing the block IDs of each cell
  Teuchos::RCP<const Epetra_BlockMap> ownedMap = getGlobalOwnedMap(1);
  Teuchos::RCP<const Epetra_BlockMap> overlapMap = getGlobalOverlapMap(1);
  Epetra_Vector blockIDs(*overlapMap);
  Epetra_Import importer(*overlapMap, *ownedMap);
  Teuchos::RCP<Epetra_Vector> ownedBlockIDs = getBlockID();
  blockIDs.Import(*ownedBlockIDs, importer, Insert);

  // Apply the block bonding matrix and create a new NeighborhoodData
  Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
  neighborhoodData->SetNumOwned(unfilteredNeighborhoodData->NumOwnedPoints());
  memcpy(neighborhoodData->OwnedIDs(), unfilteredNeighborhoodData->OwnedIDs(), neighborhoodData->NumOwnedPoints()*sizeof(int));
  vector<int> neighborhoodListVec;
  neighborhoodListVec.reserve(unfilteredNeighborhoodData->NeighborhoodListSize());
  int* const neighborhoodPtr = neighborhoodData->NeighborhoodPtr();

  int numOwnedPoints = neighborhoodData->NumOwnedPoints();
  int* const unfilteredNeighborhoodList = unfilteredNeighborhoodData->NeighborhoodList();
  int unfilteredNeighborhoodListIndex(0);
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    int blockID = static_cast<int>(blockIDs[iID]);
	int numUnfilteredNeighbors = unfilteredNeighborhoodList[unfilteredNeighborhoodListIndex++];
    unsigned int numNeighborsIndex = neighborhoodListVec.size();
    neighborhoodListVec.push_back(-1); // placeholder for number of neighbors
    int numNeighbors = 0;
	for(int iNID=0 ; iNID<numUnfilteredNeighbors ; ++iNID){
      int unfilteredNeighborID = unfilteredNeighborhoodList[unfilteredNeighborhoodListIndex++];
      int unfilteredNeighborBlockID = static_cast<int>(blockIDs[unfilteredNeighborID]);
      if(blockBondingMatrix[blockID-1][unfilteredNeighborBlockID-1] == true){
        neighborhoodListVec.push_back(unfilteredNeighborID);
        numNeighbors += 1;
      }
    }
    neighborhoodListVec[numNeighborsIndex] = numNeighbors;
    neighborhoodPtr[iID] = numNeighborsIndex;
  }

  neighborhoodData->SetNeighborhoodListSize(neighborhoodListVec.size());
  memcpy(neighborhoodData->NeighborhoodList(), &neighborhoodListVec[0], neighborhoodListVec.size()*sizeof(int));

  return neighborhoodData;
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::TextFileDiscretization::getGlobalOwnedMap(int d) const
{
  switch (d) {
    case 1:
      return oneDimensionalMap;
      break;
    case 3:
      return threeDimensionalMap;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         std::endl << "TextFileDiscretization::getGlobalOwnedMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::TextFileDiscretization::getGlobalOverlapMap(int d) const
{
  switch (d) {
    case 1:
      return oneDimensionalOverlapMap;
      break;
    case 3:
      return threeDimensionalOverlapMap;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         std::endl << "TextFileDiscretization::getOverlapMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::TextFileDiscretization::getGlobalBondMap() const
{
  return bondMap;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getInitialX() const
{
  return initialX;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getPointAngle() const
{ 
  return pointAngle;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getPointTime() const
{
  return pointTime;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getNodeType() const
{
  return nodeType;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getHorizon() const
{
  return horizonForEachPoint;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getCellVolume() const
{
  return cellVolume;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getBlockID() const
{
  return blockID;
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> 
PeridigmNS::TextFileDiscretization::getNeighborhoodData() const
{
  return neighborhoodData;
}

unsigned int
PeridigmNS::TextFileDiscretization::getNumBonds() const
{
  return numBonds;
}

unsigned int
PeridigmNS::TextFileDiscretization::getMaxNumBondsPerElem() const
{
  return maxNumBondsPerElem;
}

