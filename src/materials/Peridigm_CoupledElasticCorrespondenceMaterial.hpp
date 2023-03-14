//! \file Peridigm_CoupledElasticCorrespondenceMaterial.hpp

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

#ifndef PERIDIGM_COUPLEDELASTICCORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_COUPLEDELASTICCORRESPONDENCEMATERIAL_HPP

#include "Peridigm_CorrespondenceMaterial.hpp"
#include "Peridigm_FEM.hpp"

namespace PeridigmNS {

  class CoupledElasticCorrespondenceMaterial : public CorrespondenceMaterial, public FEMMaterial{
  public:

    //! Constructor.
    CoupledElasticCorrespondenceMaterial(const Teuchos::ParameterList & params);

    //! Destructor.
    virtual ~CoupledElasticCorrespondenceMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("Coupled Elastic Correspondence"); }

    //! Returns the density of the material.
    virtual double Density() const override { return CorrespondenceMaterial::Density(); }

    //! Returns the bulk modulus of the material.
    virtual double BulkModulus() const override { return CorrespondenceMaterial::BulkModulus(); }

    //! Returns the shear modulus of the material.
    virtual double ShearModulus() const override { return CorrespondenceMaterial::ShearModulus(); }

    //! Returns a vector of field IDs corresponding to the variables associated with the material.
    virtual std::vector<int> FieldIds() const override { return CorrespondenceMaterial::FieldIds(); }

    //! Evaluate the Cauchy stress. --> call from Peridigm_CorrespondenceMaterial.cpp
    virtual void initialize(const double dt, 
                            const int numOwnedPoints, 
                            const int* ownedIDs,
                            const int* topology,
                            PeridigmNS::DataManager& dataManager);
    
    //! Evaluate the PD Cauchy stress.
    virtual void computeCauchyStress(const double dt,
                                     const int numOwnedPoints,
                                     PeridigmNS::DataManager& dataManager,
                                     const double time = 0.0) const override;

    //! Evaluate the FE Cauchy stress.
    virtual void computeCauchyStress(const double* strain,
                                     double* stress) const override;

    // //! Evaluate the Cauchy stress.
    // void computeCauchyStress(const double dt,
    //                           const int numOwnedPoints,
    //                           PeridigmNS::DataManager& dataManager,
    //                           const double* strain,
    //                           double* stress,
    //                           const double time = 0.0) const { 
    //                             CorrespondenceMaterial::computeCauchyStress(dt,
    //                               numOwnedPoints,
    //                               dataManager,
    //                               time); 
    //                             FEMMaterial::computeCauchyStress(strain,stress); 
    //                           }

    //! Evaluate the internal force.
    virtual void computeForce(const double dt,
                              const int numOwnedPoints,
                              const int* ownedIDs,
                              const int* neighborhoodList,
                              PeridigmNS::DataManager& dataManager,
                              const double currentTime) const override { 
                                CorrespondenceMaterial::computeForce(dt,
                                  numOwnedPoints,
                                  ownedIDs,
                                  neighborhoodList,
                                  dataManager,
                                  currentTime); 
                              }

    // virtual void computeCauchyStressFem(const double* strain,                                                  
    //                                  double* sigmaInt) const;
    //! Returns the requested material property
    //! A dummy method here.
    virtual double lookupMaterialProperty(const std::string keyname) const {return 0.0;}


  protected:

    bool m_applyThermalStrains;
    double m_alpha;

    int m_temperatureFieldId;
    int m_deltaTemperatureFieldId;
    int m_unrotatedRateOfDeformationFieldId;
    int m_unrotatedCauchyStressFieldId;
    
    double C[6][6];
    int m_type;
  };
}

#endif // PERIDIGM_COUPLEDELASTICCORRESPONDENCEMATERIAL2_HPP
