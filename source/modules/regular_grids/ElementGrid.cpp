#include "ElementGrid.h"

using namespace std;

namespace csmp {

void ElementGrid::Out() const
 {
    pair<int32,int32>  index;
    
    cout <<"\nElementGrid::Out(): printing grid...\n" << endl;
    for ( cEGridIterator it=grid.begin(); it!=grid.end(); it++ )
      {
          index = (*it).first;
          cout <<"ElementGrid["<< index.first <<","<< index.second <<"] ";
          cout << (*it).second << endl;
      }
 }

} // end namespace csmp
