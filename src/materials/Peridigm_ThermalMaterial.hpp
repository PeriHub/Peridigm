#ifndef PERIDIGM_THERMALMATERIAL_HPP
#define PERIDIGM_THERMALMATERIAL_HPP

#include "Peridigm_Material.hpp"

namespace PeridigmNS{

	class ThermalMaterial : public Material{
	public: 

// 		! Constructor
		ThermalMaterial(const Teuchos::ParameterList & params): Material(params) {};

// 		! Destructor
		virtual ~ThermalMaterial(){};

		virtual void 
		computeHeatFlow(const double dt,
						const int numOwnedPoints,
						const int* ownedIDs,
						const int* neighborhoodList, 
						PeridigmNS::DataManager& dataManager) const = 0;
	};

}
#endif // PERIDIGM_THERMALMATERIAL_HPP