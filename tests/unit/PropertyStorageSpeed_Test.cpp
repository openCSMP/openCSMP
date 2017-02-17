#include "PropertyStorageSpeed_Test.h"
#include "ANSYS_Model3D.h"
#include "BE_Time.h"

using namespace std;

namespace csmp {

PropertyStorageSpeed_Test::PropertyStorageSpeed_Test( ostream* osptr )
 : Test(osptr)
 {
 }
 

void PropertyStorageSpeed_Test::run()
 {
   cout <<"******************************************************************************************\n";
   cout << endl;
   cout <<" Property Storage I/O benchmark\n";
   cout << endl;
   cout <<"******************************************************************************************\n";
   // cout <<"\n\nmain: Enter name of input variables file: ";
   string  var_file("PropertyStorage_test1.txt");
   //cin >> var_file;
   
   ANSYS_Model3D  model( "PropertyStorage_test", var_file.c_str() );
   
   // 1. measuring the amount of memory allocated to build the model including 
   //    property storage (the storage requirements without property storage are established using 
   //    the variables file 'memory_benchmark_no_variables')
   cout <<"\nmain: Have a look at the TaskManager to compare total memory allocation at this point."<< endl;

   // 2. read write operations for all variables in the model using the CopyReplace function
   map<string,csmp::Index>  properties;
   model.Database().ListProperties( properties );
   
   BE_Time  stop_watch;
   // random numbers
   // cout <<"\nmain: Enter seed for random number: ";
   unsigned int seed(stop_watch.Second());
   // cin >> seed;
   srand(seed);
   double64            va;
   ScalarVariable      sc;  
   VectorVariable<3U>  vc;  
   TensorVariable<3U>  ts;  
   
   // timing the operations
   for ( int i=0; i<25000; i++ )
     {
        // initialising the input properties
        sc() = static_cast<double64>(rand()) / 3.;
        vc = static_cast<double64>(rand()) / 3.;
        ts = static_cast<double64>(rand()) / 3.;
        
        // writing properties
        for ( map<string,csmp::Index>::const_iterator it=properties.begin(); it!=properties.end(); it++ ) {
             if ( (*it).second.type == SCALAR )      model.InputPropertyValue( (*it).first.c_str(), sc );
             else if ( (*it).second.type == VECTOR ) model.InputPropertyValue( (*it).first.c_str(), vc );
             else if ( (*it).second.type == TENSOR ) model.InputPropertyValue( (*it).first.c_str(), ts );
          }
        // reading properties including some manipulations
        double64 pmin, pmax;  
        cout <<"\n\nmain: property ranges: ";
        for ( map<string,csmp::Index>::const_iterator it=properties.begin(); it!=properties.end(); it++ ) {
             model.MinMaxOf( (*it).first.c_str(), pmin, pmax );
             cout <<"\nmain: min/max of '"<< (*it).first <<"': "<< pmin <<" - "<< pmax;
          }
        cout <<"\n\nmain: property averages: ";
        for ( map<string,csmp::Index>::const_iterator it=properties.begin(); it!=properties.end(); it++ ) {
             va = model.Region("Model").Average( (*it).first.c_str() );
             cout <<"\nmain: average of '"<< (*it).first <<"': "<< va;
          }
     }
   cout <<"\nmain: Completed property read/write in "<< stop_watch.delta(BE_Time()) <<" secs."<< endl;
 
} // end run 

} // end csmp


