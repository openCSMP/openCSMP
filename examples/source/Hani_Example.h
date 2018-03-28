#ifndef HANI_EXAMPLE_H
#define HANI_EXAMPLE_H

#include "Example.h"
#include "ANSYS_Model3D.h"
#include "TwoPhaseModel.h"

using namespace std;

namespace csmp {

class  Hani_Example : public Example {
public:
	Hani_Example() { cout << " This is Hani Constructor \n"; Run();};
	virtual void Run();
	void computeTotalMobility(ANSYS_Model3D& mdl, TwoPhaseModel<3U>& relperm);
	virtual void Specifications();
};

} // csmp

#endif // ECLIPSE_MESH_INTERFACE_EXAMPLE_H
