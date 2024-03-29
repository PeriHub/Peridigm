//! \file Peridigm_ElasticBondBasedMaterial.hpp

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

#ifndef PERIDIGM_ELASTICBONDBASEDMATERIAL_HPP
#define PERIDIGM_ELASTICBONDBASEDMATERIAL_HPP

#include "Peridigm_Material.hpp"

namespace PeridigmNS {

  //! Elastic bond-based peridynamic material model.
  class ElasticBondBasedMaterial : public Material{
  public:

    //! Constructor.
    ElasticBondBasedMaterial(const Teuchos::ParameterList & params);

    //! Destructor.
    virtual ~ElasticBondBasedMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("Elastic Bond Based"); }

    //! Returns the density of the material.
    virtual double Density() const { return m_density; }

    //! Returns the bulk modulus of the material.
    virtual double BulkModulus() const { return m_bulkModulus; }

    //! Returns the shear modulus of the material.
    virtual double ShearModulus() const {
      double nu = 0.25;
      double shearModulus = (3.0*m_bulkModulus*(1.0 - 2.0*nu))/(2.0*(1.0 + nu));
      return shearModulus;
    }

    //! Returns a vector of field IDs corresponding to the variables associated with the material.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Initialized data containers and computes weighted volume.
    virtual void
    initialize(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               PeridigmNS::DataManager& dataManager);

    //! Evaluate the internal force.
    virtual void
    computeForce(const double dt,
                 const int numOwnedPoints,
                 const int* ownedIDs,
                 const int* neighborhoodList,
                 PeridigmNS::DataManager& dataManager,
                 const double currentTime = 0.0) const;
  virtual double lookupMaterialProperty(const std::string keyname) const 
    {
      std::map<std::string, double>::const_iterator search = materialProperties.find(keyname);
      if(search != materialProperties.end())
        return search->second;
      else
        TestForTermination(true, "**** Error: requested material property " + keyname + " is not in Multiphysics Elastic Material");
      // This is a fallthrough case to make the compiler happy.
      return 0.0;
    }
  protected:

    //! Computes the distance between nodes (a1, a2, a3) and (b1, b2, b3).
    inline double distance(double a1, double a2, double a3,
                           double b1, double b2, double b3) const
      {
        return ( sqrt( (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) + (a3-b3)*(a3-b3) ) );
      }

    // material parameters
    std::map<std::string, double> materialProperties;
    double m_bulkModulus;
    double m_density;
    double m_horizon;
    std::vector<double> m_lambda;
    double m_lambdaBed;
    double m_C;
    double m_kappa;
    double m_Tenv;
    double m_Tref;
    double m_Tbed;
    double m_factor;
    double m_surfaceCorrection;
    double m_limit;
    double alpha[3][3];
    // field spec ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_volumeFieldId;
    int m_damageFieldId;
    int m_modelCoordinatesFieldId;
    int m_coordinatesFieldId;
    int m_forceDensityFieldId;
    int m_bondDamageFieldId;
    int m_partialStressFieldId;
    int m_temperatureFieldId;
    int m_thermalFlowStateFieldId;
    int m_detachedNodesFieldId;
    int m_specificVolumeFieldId;
    int m_horizonFieldId;
    bool m_applyHeatTransfer;
    bool m_applyThermalFlow;
    bool m_applyThermalPrintBedFlow;
    bool m_applyThermalStrains;
  };
}

#endif // PERIDIGM_ELASTICBONDBASEDMATERIAL_HPP
