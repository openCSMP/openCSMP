#ifndef FV_CommunicationData_h__
#define FV_CommunicationData_h__

#include "CSP_definitions.h"
#include "VSetConnectivity.h"
#include "SuperGroup.h"
#include "PropertyDatabase.h"
#include "SAMGp_CommunicationData.h"

namespace csp {


template<typename fT, stl_index dim>
class FV_CommunicationData {
 public:
     ~FV_CommunicationData();
     FV_CommunicationData( const FV_CommunicationData& );
     FV_CommunicationData( SuperGroup<fT,dim>&, VSetConnectivity<fT,dim>&, int& );
     FV_CommunicationData( SuperGroup<fT,dim>& sg,
                           SAMGp_CommunicationData<fT,dim>& samgp_comm,
                           int& my_id );

     FV_CommunicationData& operator=( const FV_CommunicationData& );

     int nshalo() const { return isndlist.size(); } // inner halo of processor
     int npsnd()  const { return iranksnd.size(); } // number of neighboring processors
     int nrhalo() const { return ireclist.size(); }
     int nprec()  const { return irankrec.size(); } // neighboring processors which will receive data

     std::vector<int> iranksnd;  // processors to which information is send
     std::vector<int> ipts;      // needed to distinguish which data to be send to which neighbour
     std::vector<int> isndlist;  // all send variables to all neighbours
     std::vector<int> irankrec;  // processors from which information is received
     std::vector<int> iptr;      // needed to distinguish which data to be received from which neighbour
     std::vector<int> ireclist;  // all send variables to all neighbours

     void Out();

     void UpdateHaloNodes(const char* prop, bool blockin_send_recv=false, bool verbose=false);
     void SetOuterhaloDirich(const char* prop);
     void SetOuterhaloPlain(const char* prop);
     void ConstantValueOnOuterhalo(const char* prop, csp_float val );
     void ConstantValueOnInnerhalo(const char* prop, csp_float val );

 private:
     MemoryManager<fT,dim>&      pmem;
     PropertyDatabase&           pref;
     std::vector<csp_float>      isnddata;
     std::vector<csp_float>      irecdata;
     int tag;
     int my_rank;
     std::vector<MPI::Request>   mpi_req;
     ScalarVariable<csp_float>   property;
 };

} // END CSP NAMESPACE 
#endif
