#include "ElementGrid.h"

using namespace std;

namespace csmp {

void ElementGrid::Out(std::ostream& os) const
 {
    pair<int32,int32>  index;
    
    os <<"\nElementGrid::Out(): printing grid...\n" << endl;
    for ( cEGridIterator it=grid.begin(); it!=grid.end(); it++ )
      {
          index = (*it).first;
          os <<"ElementGrid["<< index.first <<","<< index.second <<"] ";
          os << (*it).second << endl;
      }
 }

} // end namespace csmp
