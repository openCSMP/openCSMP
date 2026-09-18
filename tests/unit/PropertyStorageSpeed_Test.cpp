#include "PropertyStorageSpeed_Test.h"
#include "ANSYS_Model3D.h"
#include "Region.h"
#include "BE_Time.h"

using namespace std;

namespace csmp {

PropertyStorageSpeed_Test::PropertyStorageSpeed_Test( ostream* osptr )
 : Test(osptr)
 {
 }
 
// TODO:  fix test for the different variable types; tensor does not work with random matrix entries
void PropertyStorageSpeed_Test::run()
 {
   cout <<"******************************************************************************************\n";
   cout << endl;
   cout <<" Property Storage I/O benchmark\n";
   cout << endl;
   cout <<"******************************************************************************************\n";
   string  var_file("PropertyStorage_test1.txt");
   
   ANSYS_Model3D  model( "PropertyStorage_test", var_file.c_str() );
   
   // 1. measuring the amount of memory allocated to build the model including 
   //    property storage (the storage requirements without property storage are established using 
   //    the variables file 'memory_benchmark_no_variables')
   cout <<"\nmain: Have a look at the TaskManager to compare total memory allocation at this point."<< endl;

   // 2. read write operations for all variables in the model using the CopyReplace function
   map<string,csmp::Index>  properties;
   model.Database().ListVariables( properties );
   
   BE_Time  stop_watch;
   // random numbers
   // cout <<"\nmain: Enter seed for random number: ";
   unsigned int seed(stop_watch.Second());
   // cin >> seed;
   srand(seed);
   double              va;
   ScalarVariable      sc;
   VectorVariable<3U>  vc;  
   TensorVariable<3U>  ts;  
   
   // timing the operations
   auto t0 = chrono::high_resolution_clock::now();
   for ( size_t i=0; i<100; i++ )
     {
        // initialising the input properties
        sc() = static_cast<double>(rand()) / 3.;
        vc   = static_cast<double>(rand()) / 3.;
        ts   = static_cast<double>(rand()) / 3.;
        
        // writing properties
        for ( map<string,csmp::Index>::const_iterator it=properties.begin(); it!=properties.end(); it++ ) {
             if ( (*it).second.type == SCALAR )      model.InputPropertyValue( (*it).first.c_str(), sc );
             else if ( (*it).second.type == VECTOR ) model.InputPropertyValue( (*it).first.c_str(), vc );
             else if ( (*it).second.type == TENSOR ) model.InputPropertyValue( (*it).first.c_str(), ts );
          }
        // reading properties including some manipulations
        double pmin, pmax;  
        if ( verbose_ ) cout <<"\n\nmain: property ranges: ";
        for ( map<string,csmp::Index>::const_iterator it=properties.begin(); it!=properties.end(); it++ ) {
             model.MinMaxOf( (*it).first.c_str(), pmin, pmax );
             if ( verbose_ ) cout <<"\nmain: min/max of '"<< (*it).first <<"': "<< pmin <<" - "<< pmax;
          }
        if ( verbose_ ) cout <<"\n\nmain: property averages: ";
        for ( map<string,csmp::Index>::const_iterator it=properties.begin(); it!=properties.end(); it++ ) {
             va = model.Region("Model").Average( (*it).first.c_str() );
             if ( verbose_ ) cout <<"\nmain: average of '"<< (*it).first <<"': "<< va;
          }
     }
   auto t1 = chrono::high_resolution_clock::now();
   cout <<"\n\n"<<"PropertyStorageSpeed_Test::run: Completed property read/write in ";
   cout << chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
   cout <<" milliseconds."<< endl;
 
} // end run 

} // end csmp


