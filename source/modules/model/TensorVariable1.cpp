#include "TensorVariable1.h"

using namespace std;

namespace csmp {


void  TensorVariable<1U>::In()
 {
     string  status;
     
     cout <<"\nEnter ["<< 1U <<"]["<< 1U <<"] tensor variable status: ";
     cin >> status;
     flag = parseStatus( status.c_str() );
     
     cout <<"\nEnter single element: ";
     cin >> data;

 } // end In





/// @test tested: O.K.
void  TensorVariable<1U>::Out() const
 {
     cout <<"\nStatus: "<< parseStatus( flag ) << endl;
     cout <<"\nValue:  "<< data << endl;

 } // end Out




// stubs

bool  TensorVariable<1U>::EigenValues( VectorVariable<1U>& vecEigenvalues ) const
 {
    vecEigenvalues = data;
    return true;
 }



bool  TensorVariable<1U>::Eigen( VectorVariable<1U>& vvEigenvalues, 
                                 TensorVariable<1U>& tvEigenvectors,
                                 bool bNormalize ) const
 {
    vvEigenvalues  = data;
    tvEigenvectors = data;
    if ( bNormalize ) tvEigenvectors = 1.;
    return true;
 }

} // end namespace csmp

