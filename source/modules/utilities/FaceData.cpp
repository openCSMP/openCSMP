#include <climits>

#include "FaceData.h"
#include "Box.h"

using namespace std;

namespace csmp {

FaceData& FaceData::operator=( const FaceData& fd )
 {
    if ( &fd != this ) {
         key      = fd.key;
         nodes    = fd.nodes;
         efnumber = fd.efnumber;
         nbors    = fd.nbors;
         etype    = fd.etype;
      }
    return *this;
 }



void  FaceData::Reset()
 {
     key.erase( key.begin(), key.end() );
     nodes.erase( nodes.begin(), nodes.end() );
     efnumber.first  = ULONG_MAX;
     efnumber.second = ULONG_MAX;
     nbors.first     = IRREGULAR_OUTSIDE;
     nbors.second    = IRREGULAR_OUTSIDE;
     etype           = UNKNOWN;
 }

} // end namespace csmp
