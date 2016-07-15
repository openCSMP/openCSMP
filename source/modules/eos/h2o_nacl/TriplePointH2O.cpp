#include "TriplePointH2O.h"

using namespace std;

namespace csmp
{

  TriplePointH2O::TriplePointH2O(){}
  TriplePointH2O::~TriplePointH2O(){}

  /** Output H2O triple point temperature in [C] */
  double64 TriplePointH2O::Temperature()  const { return 0.01e0; }


  /** Output H2O triple point pressure in [Pa] */
  double64 TriplePointH2O::Pressure()     const { return 0.00611731672023643346e5; }

}//csmp
