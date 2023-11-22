#ifndef MultipleVariableCommunicator_h__
#define MultipleVariableCommunicator_h__

#include "CSP_definitions.h"
#include "VSetConnectivity.h"
#include "SuperGroup.h"
#include "PropertyDatabase.h"
#include "FV_CommunicationData.h"

namespace csp {


template<typename fT, stl_index dim>
class MultipleVariableCommunicator {
 public:
     ~MultipleVariableCommunicator();
     MultipleVariableCommunicator( SuperGroup<fT,dim>& ,FV_CommunicationData<fT,dim>& ,std::vector<const char*>& );

     void UpdateHaloNodes();

 private:
     FV_CommunicationData<fT,dim>         FV_comm;
     MemoryManager<fT,dim>&               pmem;
     PropertyDatabase&                    pref;
     ScalarVariable<csp_float>            property;
     std::vector<std::vector<csp_float> > isnddata,
                                          irecdata;
     int                                  my_rank;
     std::vector<const char*>             properties;
     Index                                prop_key;
     int                                  tag;
 };

} // END CSP NAMESPACE 
#endif
