// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "TensorVariable1.h"

using namespace std;

namespace csmp {

double& TensorVariable<1U>::operator()( uint32_t, uint32_t ) noexcept
 {
    return data; 
 }


const double& TensorVariable<1U>::operator()( uint32_t, uint32_t ) const noexcept
 {
    return data; 
 }



double TensorVariable<1U>::Component( uint32_t ) const noexcept
 {
    return data; 
 }



void  TensorVariable<1U>::Component( uint32_t, double val ) noexcept
 {
    data = val; 
 }



VARIABLE_FLAG& TensorVariable<1U>::Flag( uint32_t ) noexcept
 { 
    return flag; 
 }



VARIABLE_FLAG  TensorVariable<1U>::Flag( uint32_t ) const noexcept
 { 
    return flag; 
 }



TensorVariable<1U>::TensorVariable() noexcept
   : flag(ANY), data(std::numeric_limits<double>::quiet_NaN())
  {
  }





TensorVariable<1U>::TensorVariable( VARIABLE_FLAG f, double val ) noexcept
 : flag(f), data(val)
 {
 } 
 




// here the flag of the lefthand tensor-variable is sustained

TensorVariable<1U>  TensorVariable<1U>::operator+( const TensorVariable<1U>& t ) const noexcept
 {
    return TensorVariable<1U>( flag, t.data+data );
 }



TensorVariable<1U>  TensorVariable<1U>::operator-( const TensorVariable<1U>& t ) const noexcept
 {
      return TensorVariable<1U>( flag, data-t.data );
 }


TensorVariable<1U>  TensorVariable<1U>::operator+( double val ) const noexcept
 {
      return TensorVariable<1U>( flag, data+val );
 }
 
 

TensorVariable<1U>  TensorVariable<1U>::operator-( double val ) const noexcept
 {
      return TensorVariable<1U>( flag, data-val );
 }
 
 

TensorVariable<1U>  TensorVariable<1U>::operator*( double val ) const noexcept
 {
      return TensorVariable<1U>( flag, data*val );
 }
 
 

TensorVariable<1U>  TensorVariable<1U>::operator/( double val ) const noexcept
 {
      return TensorVariable<1U>( flag, data/val );
 }


// matrix vector multiplication: v = M * v

VectorVariable<1U>  TensorVariable<1U>::operator*( const VectorVariable<1U>& vc ) const noexcept
 {
    return VectorVariable<1U>( flag, data * vc[0] );
 } 


Point<1U>  TensorVariable<1U>::operator*( const Point<1U>& v ) const noexcept
 {
    return Point<1U>( data * v[0] );
 } 




TensorVariable<1U> TensorVariable<1U>::Adjoint() const noexcept
 {
      return TensorVariable( flag, 1. );
 }



TensorVariable<1U>  TensorVariable<1U>::operator*( const TensorVariable<1U>& ts ) const noexcept 
 {
   return TensorVariable<1U>( flag, data * ts.data );
 } 




TensorVariable<1U>&  TensorVariable<1U>::operator+=( const ScalarVariable& sc ) noexcept
 {
    data += sc();
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator-=( const ScalarVariable& sc ) noexcept
 {
    data -= sc();
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator*=( const ScalarVariable& sc ) noexcept
 {
    data *= sc();
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator/=( const ScalarVariable& sc ) noexcept
 {
    data /= sc();
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator+=( const TensorVariable<1U>& ts ) noexcept
 {
    data += ts.data;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator-=( const TensorVariable<1U>& ts ) noexcept
 {
    data -= ts.data;
    return *this; 
 }


// element by element division

TensorVariable<1U>&  TensorVariable<1U>::operator/=( const TensorVariable<1U>& ts ) noexcept
 {
    data /= ts.data;
    return *this; 
 }



// element by element division

TensorVariable<1U>  TensorVariable<1U>::operator/( const TensorVariable<1U>& ts ) const noexcept
 {
    return TensorVariable<1U>( flag, data / ts.data ); 
 }



// matrix multiplication

TensorVariable<1U>&  TensorVariable<1U>::operator*=( const TensorVariable<1U>& ts ) noexcept
 {
    data *= ts.data;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator+=( double val ) noexcept
 {
    data += val;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator-=( double val ) noexcept
 {
    data -= val;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator*=( double val ) noexcept
 {
    data *= val;
    return *this; 
 }



TensorVariable<1U>&  TensorVariable<1U>::operator/=( double val ) noexcept
 {
    data /= val;
    return *this; 
 }

// --------------------
// ASSIGNMENT OPERATORS
// --------------------


TensorVariable<1U>&  TensorVariable<1U>::operator=( double val ) noexcept
 {
    data = val;
    return *this; 
 }


// the flag is adopted from the scalar variable

TensorVariable<1U>&  TensorVariable<1U>::operator=( const ScalarVariable& sc ) noexcept
 {
    flag = sc.Flag();
    data = sc();
    return *this; 
 }



// writes vector into the diagonal of the zero'd tensor

TensorVariable<1U>&  TensorVariable<1U>::operator=( const VectorVariable<1U>& vc ) noexcept
 {
    flag = vc.Flag();
    data = vc[0];
    return *this; 
 }





bool  TensorVariable<1U>::operator==( const TensorVariable<1U>& ts ) const noexcept
 {
    return ( data == ts.data && flag == ts.flag );
 }

 
  

bool  TensorVariable<1U>::operator!=( const TensorVariable<1U>& t ) const noexcept
 {
     return !(*this == t);
 } 



bool  TensorVariable<1U>::operator<( const TensorVariable<1U>& t ) const noexcept
 {
     return (this < &t);
 } 


// -------
// METHODS
// -------



void TensorVariable<1U>::Identity() noexcept
 {
    data = static_cast<double>(1.0);
 }



TensorVariable<1U>  TensorVariable<1U>::Transposed() const noexcept
 {
    return TensorVariable( flag, data );
 }



double TensorVariable<1U>::Determinant() const noexcept
 {
    return data;
 }
 
double TensorVariable<1U>::Trace() const noexcept
 {
    return data;
 }

TensorVariable<1U> TensorVariable<1U>::Inverse() const noexcept
 {
    return TensorVariable( flag, 1. / data );
 }
 
 


 double  TensorVariable<1U>::MinElement() const noexcept
 {
    return data;
 }
  
  


 double  TensorVariable<1U>::MaxElement() const noexcept
 {
    return data;
 }
 
 


bool  TensorVariable<1U>::IsWithinRange( double vmin, double vmax ) const noexcept
 {
    if ( data < vmin || data > vmax ) return false;
    return true;
 }
 
 
 
// vector-matrix multiplication: v^T = (v^T * A)^T = A^T v  

  VectorVariable<1U>  operator*( const VectorVariable<1U>& vc, const TensorVariable<1U>& ts ) noexcept
  {
    return VectorVariable<1U>( vc.Flag(), ts(0,0) * vc[0] );
  }
  
  // vector-matrix multiplication: v^T = (v^T * A)^T = A^T v

  Point<1U>  operator*( const Point<1U>& vc, const TensorVariable<1U>& ts ) noexcept
  {
    return Point<1U>( ts(0,0) * vc[0] );
  }
  


void TensorVariable<1U>::AssignToRow( uint32_t, VectorVariable<1U>& vc ) noexcept
{
	flag = vc.Flag(0U); 
	data = vc[0U];
}


void TensorVariable<1U>::AssignToColumn( uint32_t, VectorVariable<1U>& vc ) noexcept
{
	flag = vc.Flag(0U);
	data = vc[0U];
}




VectorVariable<1U> TensorVariable<1U>::Row( uint32_t ) const noexcept
{
	return VectorVariable<1U>( flag, data );
}


VectorVariable<1U> TensorVariable<1U>::Column( uint32_t ) const noexcept
{
	return VectorVariable<1U>( flag, data );
}

bool TensorVariable<1U>::Out( std::fstream& fp ) const
{
  const int32_t flag_0( flag );
  fp.write( (char*)&flag_0, sizeof( int32_t ) ); // VARIABLE_FLAG
  fp.write( (char*)&data, sizeof( double ) );
  return true;
}

bool TensorVariable<1U>::In( std::fstream& fp )
{
  fp.read( (char*)&flag, sizeof( int32_t ) ); // VARIABLE_FLAG
  fp.read( (char*)&data, sizeof( double ) );
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
void  TensorVariable<1U>::Out() const noexcept
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
    tvEigenvectors = bNormalize ? 1. : data;
    return true;
 }
 
 
// dummy function, normalises Eigen vectors
bool TensorVariable<1U>::EigenWeaklyNonSymmetric( VectorVariable<1U>& eigenVals,
                                                  TensorVariable<1U>& eigenVecs, double ) const
{
    eigenVals = data;
    eigenVecs = 1.;
    return true;
}

} // end namespace csmp

