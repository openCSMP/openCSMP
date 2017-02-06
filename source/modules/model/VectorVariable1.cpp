#include "VectorVariable1.h"

using namespace std;

namespace csmp {


void  VectorVariable<1U>::In()
 {
     string  status;
     cout <<"\nEnter status for x-component of variable: ";
     cin >> status;
     flag = parseStatus( status.c_str() );
     
     cout <<"\nEnter vector element: ";
     cin >> data;

 } // end In


/// Multiplies by negative unity vector
void  VectorVariable<1U>::Invert()
  {
    data *= -1.;
  }



void  VectorVariable<1U>::Out(std::ostream& os) const 
 {
     os <<"\nStatus: "<< parseStatus(flag) <<"\t\t";
     os << endl;
     os << data <<"\t\t";
     os << endl;
 } // end Out
 
 
 
} // end namespace csmp

