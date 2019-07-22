#ifndef NIMBLE_REGION_TEST_H
#define NIMBLE_REGION_TEST_H

#include "Test.h"
#include "TwoPhaseModel.h"
#include "PropertyDatabase.h"
#include "NimbleRegion.h"
#include "Region.h"

namespace csmp {

class NimbleRegion_Test : public Test	{
	public:
		std::string mesh_name;
		void run();
  
	private:
    void computeTotalMobility( Region<2U>&, const PropertyDatabase<2U>&, TwoPhaseModel<2U>& );
		void computeTotalMobility( NimbleRegion<2U>&, const PropertyDatabase<2U>&, TwoPhaseModel<2U>& );
		void constraintPlumeBoundary( NimbleRegion<2U>& ,const PropertyDatabase<2U>&);
		void releasePlumeBoundary( NimbleRegion<2U>&, const PropertyDatabase<2U>&);
	};

} // csmp

#endif // TRANSIENT_PRESSURE_TEST
