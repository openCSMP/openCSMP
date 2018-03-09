#include "TensorVariable1.h"

using namespace std;

namespace csmp {

double64& TensorVariable<1U>::operator()( size_t, size_t ) 
 {
    return data; 
 }


const double64& TensorVariable<1U>::operator()( size_t, size_t ) const
 {
    return data; 
 }



double64 TensorVariable<1U>::Component( size_t ) const
 {
    return data; 
 }



void  TensorVariable<1U>::Component( size_t, double64 val ) 
 {
    data = val; 
 }



VARIABLE_FLAG& TensorVariable<1U>::Flag( size_t )
 { 
    return flag; 
 }



VARIABLE_FLAG  TensorVariable<1U>::Flag( size_t ) const 
 { 
    return flag; 
 }

size_t TensorVariable<1U>::Size() const
  {
    return 1U;
  }

void TensorVariable<1U>::Resize( size_t, double64 newValue )
  {
    data = newValue;
  }


TensorVariable<1U>::TensorVariable()
   : flag(ANY), data(std::numeric_limits<double64>::quiet_NaN())
  {
  }



TensorVariable<1U>&  TensorVariable<1U>::operator=( const TensorVariable<1U>& ts )
 {
    if ( &ts != this ) {
         flag = ts.flag;
         data = ts.data;
      }  
    return *this; 
 }




TensorVariable<1U>::TensorVariable( const TensorVariable<1U>& t )
 {
    *this = t;
 }




TensorVariable<1U>::TensorVariable( VARIABLE_FLAG f, double64 val )
 : flag(f), data(val)
 {
 } 
 
                                   
 

TensorVariable<1U>::~TensorVariable() {}




// here the flag of the lefthand tensor-variable is sustained

TensorVariable<1U>  TensorVariable<1U>::operator+( const TensorVariable<1U>& t ) const
 {
    return TensorVariable<1U>( flag, t.data+data );
 }



TensorVariable<1U>  TensorVariable<1U>::operator-( const TensorVariable<1U>& t ) const
 {
      return TensorVariable<1U>( flag, data-t.data );
 }


TensorVariable<1U>  TensorVariable<1U>::operator+( double64 val ) const
 {
      return TensorVariable<1U>( flag, data+val );
 }
 
 

TensorVariable<1U>  TensorVariable<1U>::operator-( double64 val ) const
 {
      return TensorVariable<1U>( flag, data-val );
 }
 
 

TensorVariable<1U>  TensorVariable<1U>::operator*( double64 val ) const
 {
      return TensorVariable<1U>( flag, data*val );
 }
 
 

TensorVariable<1U>  TensorVariable<1U>::operator/( double64 val ) const
 {
      return TensorVariable<1U>( flag, data/val );
 }


// matrix vector multiplication: v = M * v

VectorVariable<1U>  TensorVariable<1U>::operator*( const VectorVariable<1U>& vc ) const
 {
    return VectorVariable<1U>( flag, data * vc[0] );
 } 


Point<1U>  TensorVariable<1U>::operator*( const Point<1U>& v ) const
 {
    return Point<1U>( data * v[0] );
 } 




TensorVariable<1U> TensorVariable<1U>::Adjoint() const
 {
      return TensorVariable( flag, 1. );
 }



TensorVariable<1U>  TensorVariable<1U>::operator*( const TensorVariable<1U>& ts ) const 
 {
   return TensorVariable<1U>( flag, data * ts.data );
 } 




TensorVariable<1U>&  TensorVariable<1U>::operator+=( const ScalarVariable& sc )
 {
    data += sc();
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator-=( const ScalarVariable& sc )
 {
    data -= sc();
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator*=( const ScalarVariable& sc )
 {
    data *= sc();
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator/=( const ScalarVariable& sc )
 {
    data /= sc();
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator+=( const TensorVariable<1U>& ts )
 {
    data += ts.data;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator-=( const TensorVariable<1U>& ts )
 {
    data -= ts.data;
    return *this; 
 }


// element by element division

TensorVariable<1U>&  TensorVariable<1U>::operator/=( const TensorVariable<1U>& ts )
 {
    data /= ts.data;
    return *this; 
 }



// element by element division

TensorVariable<1U>  TensorVariable<1U>::operator/( const TensorVariable<1U>& ts ) const
 {
    return TensorVariable<1U>( flag, data / ts.data ); 
 }



// matrix multiplication

TensorVariable<1U>&  TensorVariable<1U>::operator*=( const TensorVariable<1U>& ts )
 {
    data *= ts.data;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator+=( double64 val )
 {
    data += val;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator-=( double64 val )
 {
    data -= val;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator*=( double64 val )
 {
    data *= val;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator/=( double64 val )
 {
    data /= val;
    return *this; 
 }

// --------------------
// ASSIGNMENT OPERATORS
// --------------------


TensorVariable<1U>&  TensorVariable<1U>::operator=( double64 val )
 {
    data = val;
    return *this; 
 }


// the flag is adopted from the scalar variable

TensorVariable<1U>&  TensorVariable<1U>::operator=( const ScalarVariable& sc )
 {
    flag = sc.Flag();
    data = sc();
    return *this; 
 }



// writes vector into the diagonal of the zero'd tensor

TensorVariable<1U>&  TensorVariable<1U>::operator=( const VectorVariable<1U>& vc )
 {
    flag = vc.Flag();
    data = vc[0];
    return *this; 
 }





bool  TensorVariable<1U>::operator==( const TensorVariable<1U>& ts ) const
 {
    return ( data == ts.data && flag == ts.flag );
 }

 
  

bool  TensorVariable<1U>::operator!=( const TensorVariable<1U>& t ) const
 {
     return !(*this == t);
 } 



bool  TensorVariable<1U>::operator<( const TensorVariable<1U>& t ) const
 {
     return (this < &t);
 } 


// -------
// METHODS
// -------



void TensorVariable<1U>::Identity()
 {
    data = static_cast<double64>(1.0);
 }



TensorVariable<1U>  TensorVariable<1U>::Transposed() const
 {
    return TensorVariable( flag, data );
 }



double64 TensorVariable<1U>::Determinant() const
 {
    return data;
 }
 
double64 TensorVariable<1U>::Trace() const
 {
    return data;
 }

TensorVariable<1U> TensorVariable<1U>::Inverse() const
 {
    return TensorVariable( flag, 1. / data );
 }
 
 


 double64  TensorVariable<1U>::MinElement() const
 {
    return data;
 }
  
  


 double64  TensorVariable<1U>::MaxElement() const
 {
    return data;
 }
 
 


bool  TensorVariable<1U>::IsWithinRange( double64 vmin, double64 vmax ) const
 {
    if ( data < vmin || data > vmax ) return false;
    return true;
 }
 
 
 
// vector-matrix multiplication: v^T = (v^T * A)^T = A^T v  

VectorVariable<1U>  operator*( const VectorVariable<1U>& vc, const TensorVariable<1U>& ts )
 {
    return VectorVariable<1U>( vc.Flag(), ts(0,0) * vc[0] );
 }



void TensorVariable<1U>::AssignToRow( size_t, VectorVariable<1U>& vc )
{
	flag = vc.Flag(0U); 
	data = vc[0U];
}


void TensorVariable<1U>::AssignToColumn( size_t, VectorVariable<1U>& vc )
{
	flag = vc.Flag(0U);
	data = vc[0U];
}




VectorVariable<1U> TensorVariable<1U>::Row( size_t ) const
{
	return VectorVariable<1U>( flag, data );
}


VectorVariable<1U> TensorVariable<1U>::Column( size_t ) const
{
	return VectorVariable<1U>( flag, data );
}

 bool TensorVariable<1>::Out( FILE* fp ) const
  {
  fwrite( (void*)this, sizeof(TensorVariable<1>), 1, fp );
  return true;
  }

 bool TensorVariable<1>::In( FILE* fp )
  {
  fread( (void*)this, sizeof(TensorVariable<1>), 1, fp );
  return true;
  }


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

