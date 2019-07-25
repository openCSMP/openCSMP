#ifndef NIMBLE_REGION_TEST_H
#define NIMBLE_REGION_TEST_H

#include "Test.h"
#include "TwoPhaseModel.h"
#include "PropertyDatabase.h"
#include "NimbleRegion.h"
#include "Region.h"

namespace csmp {

template<size_t> class Model;

/**
     Models the propagation of a saturation front, updating pressure only in the region where
     saturation changes, including a one element-wide halo, using the NimbleRegion to keep track
     of this region.
*/
class NimbleRegion_Test : public Test	{
	public:
		void run();
  
	private:
    // auxiliary functions
    void ComputeTotalMobility( Region<2U>&, const PropertyDatabase<2U>&, TwoPhaseModel<2U>& );
		void ComputeTotalMobility( NimbleRegion<2U>&, const PropertyDatabase<2U>&, TwoPhaseModel<2U>& );
		void ConstrainPlumeBoundary( NimbleRegion<2U>& ,const PropertyDatabase<2U>&);
		void ReleasePlumeBoundary( NimbleRegion<2U>&, const PropertyDatabase<2U>&);
	};


void  compute_vt_AtBaryCenter( Model<2U>& model, const char* flow_domain );

} // csmp

#endif // NIMBLE_REGION_TEST_H
