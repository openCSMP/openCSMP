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



void  VectorVariable<1U>::Out() const 
 {
     cout <<"\nStatus: "<< parseStatus(flag) <<"\t\t";
     cout << endl;
     cout << data <<"\t\t";
     cout << endl;

 } // end Out
 
 
 
} // end namespace csmp

