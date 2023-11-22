#ifndef SAMGp_CommunicationData_h__
#define SAMGp_CommunicationData_h__

#include "CSP_definitions.h"
#include "VSetConnectivity.h"

namespace csp {

template<typename fT, stl_index dim>
struct SAMGp_CommunicationData {
     SAMGp_CommunicationData();
     ~SAMGp_CommunicationData();
     SAMGp_CommunicationData( const SAMGp_CommunicationData& );
     SAMGp_CommunicationData( VSetConnectivity<fT,dim>& VSetConn, int& my_rank );

     SAMGp_CommunicationData& operator=( const SAMGp_CommunicationData& );

     void Out();

     int nshalo() const { return isndlist.size(); } // inner halo of processor
     int npsnd() const { return iranksnd.size(); } // number of neighboring processors
     int nrhalo() const { return ireclist.size(); }
     int nprec() const { return irankrec.size(); } // neighboring processors which will receive data

     std::vector<int> iranksnd;  // processors to which information is send
     std::vector<int> ipts;      // needed to distinguish which data to be send to which neighbour
     std::vector<int> isndlist;  // all send variables to all neighbours
     std::vector<int> irankrec;  // processors from which information is received
     std::vector<int> iptr;      // needed to distinguish which data to be received from which neighbour
     std::vector<int> ireclist;  // all send variables to all neighbours

     void DumpToTextFile(const char* fname);
     bool InitializeFromTextFile(const char* fname, int my_rank);

  private:
     int my_rank_;
 };

} // END CSP NAMESPACE 
#endif
