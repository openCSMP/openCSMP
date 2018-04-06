#include "TensorVariable2.h"
#include <cassert>

using namespace std;

namespace csmp {

/**
    @todo SKM initialise array elements in the initialiser list
*/
TensorVariable<2U>::TensorVariable()
  : flag{{ANY,ANY}},
    data{numeric_limits<double64>::quiet_NaN(),numeric_limits<double64>::quiet_NaN(),
         numeric_limits<double64>::quiet_NaN(),numeric_limits<double64>::quiet_NaN()}
  {
  }



TensorVariable<2U>&  TensorVariable<2U>::operator=( const TensorVariable<2U>& ts )
 {
    if ( &ts != this ) 
      {
         flag = ts.flag;
         data = ts.data;
      }  
    return *this; 
 }



/**
    @todo SKM do not use assigment when constructing a new object
*/
TensorVariable<2U>::TensorVariable( const TensorVariable<2U>& t )
 : flag(t.flag),
   data(t.data)
 {
 }



/**
    initialises variable as diagonal isotropic tensor with flag and value
    SKM 17/6/2015
*/
TensorVariable<2U>::TensorVariable( VARIABLE_FLAG f, double64 val )
  : flag{{f,f}},
    data{val,0.,0.,val}
 {
 }




TensorVariable<2U>::TensorVariable( VARIABLE_FLAG f, 
                                    double64 v11, double64 v12,
                                    double64 v21, double64 v22 )
  : flag{{f,f}},
    data{v11,v12,v21,v22}
 {
 }
 
 
 
TensorVariable<2U>::TensorVariable( const VARIABLE_FLAG& f11, const VARIABLE_FLAG& f22,
                                    const double64&  v11, const double64&  v12,
                                    const double64&  v21, const double64&  v22 )
  : flag{{f11,f22}},
    data{v11,v12,v21,v22}
 {
 }





double64& TensorVariable<2U>::operator()( size_t i, size_t j ) 
 {
#ifndef NDEBUG 
    if ( i >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::operator(): row access violation, i="<< i << std::endl;
         return data[0][0];
      }
    if ( j >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::operator(): column access violation, j="<< j << std::endl;
         return data[0][0];
      }
#endif
    return data[i][j]; 
 }




const double64& TensorVariable<2U>::operator()( size_t i, size_t j ) const
 {
#ifndef NDEBUG 
    if ( i >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::operator(): row access violation, i="<< i << std::endl;
         return data[0][0];
      }
    if ( j >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::operator(): column access violation, j="<< j << std::endl;
         return data[0][0];
      }
#endif
    return data[i][j]; 
 }



void TensorVariable<2U>::Component( size_t i, double64 val )
 { 
    assert( i < Size() );
    // row by row
    if ( i == 0U )      data[0U][0U] = val;
    else if ( i == 1U ) data[0U][1U] = val;
    else if ( i == 2U ) data[1U][0U] = val;
    else                data[1U][1U] = val; // remaining case
 }



double64 TensorVariable<2U>::Component( size_t i ) const
 { 
    assert( i < Size() );
    // row by row
    if ( i == 0U ) return data[0U][0U];
    if ( i == 1U ) return data[0U][1U];
    if ( i == 2U ) return data[1U][0U];
    return data[1U][1U]; // remaining case
 }


 
 

VARIABLE_FLAG& TensorVariable<2U>::Flag( size_t i )      
 { 
#ifndef NDEBUG 
    if ( i >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::Flag(): diagonal access violation, i="<< i << std::endl;
         return flag[0];
      }
#endif
    return flag[i]; 
 }



VARIABLE_FLAG  TensorVariable<2U>::Flag( size_t i ) const 
 { 
#ifndef NDEBUG 
    if ( i >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::Flag(): diagonal access violation, i="<< i << std::endl;
         return flag[0];
      }
#endif
    return flag[i]; 
 }


size_t TensorVariable<2U>::Size() const
  {
    return 4U;
  }

void TensorVariable<2U>::Resize( size_t, double64 newValue )
  {
    data[0][0] = newValue;
    data[0][1] = newValue;
    data[1][0] = newValue;
    data[1][1] = newValue;
  }
 
                                   
 

TensorVariable<2U>::~TensorVariable() {}


// keep for storage of tensors in associative containers
bool  TensorVariable<2U>::operator==( const TensorVariable<2U>& ts ) const
 {
    return ( data == ts.data && flag == ts.flag );
 }


bool  TensorVariable<2U>::operator!=( const TensorVariable<2U>& t ) const
 {
     return !(*this == t);
 } 



bool  TensorVariable<2U>::operator<( const TensorVariable<2U>& t ) const
 {
     return (this < &t);
 } 


// -------
// METHODS
// -------



// re-tested: SKM 29-9-2001
// @test tested: O.K.
double64 TensorVariable<2U>::Determinant() const
 {
    return data[0][0]*data[1][1] - data[0][1]*data[1][0];
 }
 


double64 TensorVariable<2U>::Trace() const
 {
    return data[0][0]+data[1][1];
 }



void TensorVariable<2U>::AssignToRow( size_t iRow, VectorVariable<2U>& vc )
{
	if ( iRow == 0U )
	    flag[0U] = vc.Flag(0U);
	else
	    flag[1U] = vc.Flag(1U);
	data[iRow][0U] = vc[0U];
	data[iRow][1U] = vc[1U];
}


void TensorVariable<2U>::AssignToColumn( size_t iCol, VectorVariable<2U>& vc )
{
	if ( iCol == 0U )
	    flag[0U] = vc.Flag(0U);
	else
	    flag[1U] = vc.Flag(1U);
	data[0U][iCol] = vc[0U];
	data[1U][iCol] = vc[1U];
}




VectorVariable<2U> TensorVariable<2U>::Row( size_t iRow ) const
{
	return VectorVariable<2U>( flag[iRow], flag[iRow], 
	                           data[iRow][0U], data[iRow][1U] );
}


VectorVariable<2U> TensorVariable<2U>::Column( size_t iCol ) const
{
	return VectorVariable<2U>( flag[iCol], flag[iCol], 
	                           data[0U][iCol], data[1U][iCol] );
}


 bool TensorVariable<2>::Out( FILE* fp ) const
  {
  fwrite( (void*)this, sizeof(TensorVariable<2>), 1, fp );
  return true;
  }

 bool TensorVariable<2>::In( FILE* fp )
  {
  fread( (void*)this, sizeof(TensorVariable<2>), 1, fp );
  return true;
  }




// here the flag of the lefthand tensor-variable is sustained

TensorVariable<2U>  TensorVariable<2U>::operator+( const TensorVariable<2U>& t ) const
 {
      return TensorVariable( flag[0], flag[1],  
                             t.data[0][0]+data[0][0], t.data[0][1]+data[0][1],
                             t.data[1][0]+data[1][0], t.data[1][1]+data[1][1] );
 }



TensorVariable<2U>  TensorVariable<2U>::operator-( const TensorVariable<2U>& t ) const
 {
      return TensorVariable( flag[0], flag[1], 
                             data[0][0]-t.data[0][0], data[0][1]-t.data[0][1],
                             data[1][0]-t.data[1][0], data[1][1]-t.data[1][1] );
 }


TensorVariable<2U>  TensorVariable<2U>::operator+( double64 val ) const
 {
      return TensorVariable( flag[0], flag[1],  
                             data[0][0]+val, data[0][1]+val, 
                             data[1][0]+val, data[1][1]+val );
 }
 
 

TensorVariable<2U>  TensorVariable<2U>::operator-( double64 val ) const
 {
      return TensorVariable( flag[0], flag[1], 
                             data[0][0]-val, data[0][1]-val,
                             data[1][0]-val, data[1][1]-val );
 }
 
 

TensorVariable<2U>  TensorVariable<2U>::operator*( double64 val ) const
 {
      return TensorVariable( flag[0], flag[1],  
                             data[0][0]*val, data[0][1]*val,
                             data[1][0]*val, data[1][1]*val );
 }
 
 

TensorVariable<2U>  TensorVariable<2U>::operator/( double64 val ) const
 {
      return TensorVariable( flag[0], flag[1], 
                             data[0][0]/val, data[0][1]/val,
                             data[1][0]/val, data[1][1]/val );
 }


// matrix vector multiplication: v = M * v
// @test tested: O.K. SKM 29-9-2001

VectorVariable<2U>  TensorVariable<2U>::operator*( const VectorVariable<2U>& vc ) const
 {
    VectorVariable<2U> temp( vc.Flag(0), vc.Flag(1),
                             data[0][0] * vc[0] + data[0][1] * vc[1],
                             data[1][0] * vc[0] + data[1][1] * vc[1] );
    return temp;
 } 

Point<2U>  TensorVariable<2U>::operator*( const Point<2U>& v ) const
 {
    return Point<2U>(
        data[0][0] * v[0] + data[0][1] * v[1],
        data[1][0] * v[0] + data[1][1] * v[1] );
 } 



// re-tested: SKM 29-9-2001
// retested after partial specialization SKM9/12/03

TensorVariable<2U> TensorVariable<2U>::Adjoint() const
 {
    // forming B matrix
      return TensorVariable( flag[0], flag[1],
                             data[1][1], -data[1][0],
                            -data[0][1],  data[0][0] );
 }


// @test tested: O.K. 

TensorVariable<2U>  TensorVariable<2U>::operator*( const TensorVariable<2U>& ts ) const 
 {
   TensorVariable<2U> temp;

   // temp.data[i][j] += data[i][k]*ts.data[k][j];
   temp.data[0][0] = data[0][0] * ts.data[0][0] + data[0][1] * ts.data[1][0];
   temp.data[0][1] = data[0][0] * ts.data[0][1] + data[0][1] * ts.data[1][1];
   temp.data[1][0] = data[1][0] * ts.data[0][0] + data[1][1] * ts.data[1][0];
   temp.data[1][1] = data[1][0] * ts.data[0][1] + data[1][1] * ts.data[1][1];
   temp.flag[0] = flag[0];
   temp.flag[1] = flag[1];
   
   return temp;
 } 




TensorVariable<2U>&  TensorVariable<2U>::operator+=( const ScalarVariable& sc )
 {
    data[0][0] += sc();
    data[0][1] += sc();
    data[1][0] += sc();
    data[1][1] += sc();

    return *this; 
 }



TensorVariable<2U>&  TensorVariable<2U>::operator-=( const ScalarVariable& sc )
 {
    data[0][0] -= sc();
    data[0][1] -= sc();
    data[1][0] -= sc();
    data[1][1] -= sc();

    return *this; 
 }



TensorVariable<2U>&  TensorVariable<2U>::operator*=( const ScalarVariable& sc )
 {
    data[0][0] *= sc();
    data[0][1] *= sc();
    data[1][0] *= sc();
    data[1][1] *= sc();

    return *this; 
 }



TensorVariable<2U>&  TensorVariable<2U>::operator/=( const ScalarVariable& sc )
 {
    data[0][0] /= sc();
    data[0][1] /= sc();
    data[1][0] /= sc();
    data[1][1] /= sc();

    return *this; 
 }



TensorVariable<2U>&  TensorVariable<2U>::operator+=( const TensorVariable<2U>& ts )
 {
    data[0][0] += ts.data[0][0];
    data[0][1] += ts.data[0][1];
    data[1][0] += ts.data[1][0];
    data[1][1] += ts.data[1][1];

    return *this; 
 }



TensorVariable<2U>&  TensorVariable<2U>::operator-=( const TensorVariable<2U>& ts )
 {
    data[0][0] -= ts.data[0][0];
    data[0][1] -= ts.data[0][1];
    data[1][0] -= ts.data[1][0];
    data[1][1] -= ts.data[1][1];

    return *this; 
 }


// element by element division

TensorVariable<2U>&  TensorVariable<2U>::operator/=( const TensorVariable<2U>& ts )
 {
    data[0][0] /= ts.data[0][0];
    data[0][1] /= ts.data[0][1];
    data[1][0] /= ts.data[1][0];
    data[1][1] /= ts.data[1][1];

    return *this; 
 }



// element by element division

TensorVariable<2U>  TensorVariable<2U>::operator/( const TensorVariable<2U>& ts ) const
 {
    TensorVariable  temp;
 
    temp.data[0][0] = data[0][0] / ts.data[0][0];
    temp.data[0][1] = data[0][1] / ts.data[0][1];
    temp.data[1][0] = data[1][0] / ts.data[1][0];
    temp.data[1][1] = data[1][1] / ts.data[1][1];
    temp.flag[0] = flag[0];
    temp.flag[1] = flag[1];
             
    return temp; 
 }



// matrix multiplication

TensorVariable<2U>&  TensorVariable<2U>::operator*=( const TensorVariable<2U>& ts )
 {
   TensorVariable<2U> temp;

   temp.data[0][0] = data[0][0] * ts.data[0][0] + data[0][1] * ts.data[1][0];
   temp.data[0][1] = data[0][0] * ts.data[0][1] + data[0][1] * ts.data[1][1];
   temp.data[1][0] = data[1][0] * ts.data[0][0] + data[1][1] * ts.data[1][0];
   temp.data[1][1] = data[1][0] * ts.data[0][1] + data[1][1] * ts.data[1][1];
   temp.flag[0] = flag[0];
   temp.flag[1] = flag[1];
  
    return *this = temp; 
 }



TensorVariable<2U>&  TensorVariable<2U>::operator+=( double64 val )
 {
    data[0][0] += val;
    data[0][1] += val;
    data[1][0] += val;
    data[1][1] += val;

    return *this; 
 }



TensorVariable<2U>&  TensorVariable<2U>::operator-=( double64 val )
 {
    data[0][0] -= val;
    data[0][1] -= val;
    data[1][0] -= val;
    data[1][1] -= val;

    return *this; 
 }



TensorVariable<2U>&  TensorVariable<2U>::operator*=( double64 val )
 {
    data[0][0] *= val;
    data[0][1] *= val;
    data[1][0] *= val;
    data[1][1] *= val;

    return *this; 
 }



TensorVariable<2U>&  TensorVariable<2U>::operator/=( double64 val )
 {
    data[0][0] /= val;
    data[0][1] /= val;
    data[1][0] /= val;
    data[1][1] /= val;

    return *this; 
 }


// --------------------
// ASSIGNMENT OPERATORS
// --------------------


TensorVariable<2U>&  TensorVariable<2U>::operator=( double64 val )
 {
    data[0][0] = val;
    data[0][1] = val;
    data[1][0] = val;
    data[1][1] = val;

    return *this; 
 }


// the flag is adopted from the scalar variable

TensorVariable<2U>&  TensorVariable<2U>::operator=( const ScalarVariable& sc )
 {
    flag[0] = flag[1] = sc.Flag();
    data[0][0] = sc();
    data[0][1] = sc();
    data[1][0] = sc();
    data[1][1] = sc();

    return *this; 
 }


// writes vector into the diagonal of the zero'd tensor

TensorVariable<2U>&  TensorVariable<2U>::operator=( const VectorVariable<2U>& vc )
 {
    flag[0] = vc.Flag(0);
    flag[1] = vc.Flag(1);
    data[0][0] = vc(0);
    data[0][1] = static_cast<double64>(0.0);
    data[1][0] = static_cast<double64>(0.0);
    data[1][1] = vc(1);

    return *this; 
 }





// @test tested: O.K.

void TensorVariable<2U>::Identity()
 {
    data[0][0] = static_cast<double64>(1.0);
    data[0][1] = static_cast<double64>(0.0);
    data[1][0] = static_cast<double64>(0.0);
    data[1][1] = static_cast<double64>(1.0);
 }



void TensorVariable<2U>::DiagonalValues( double64 f_00, double64 f_11 )
{
    data[0][0] = f_00;
    data[1][1] = f_11;
}



void TensorVariable<2U>::DiagonalValues(const vector<double64>& vecDiags )
{
    data[0][0] = vecDiags[0];
    data[1][1] = vecDiags[1];
}



void TensorVariable<2U>::DiagonalValues(const VectorVariable<2U>& vecDiags )
{
    data[0][0] = vecDiags[0];
    data[1][1] = vecDiags[1];
    flag[0] = vecDiags.Flag(0);
    flag[1] = vecDiags.Flag(1);
}


/// @test tested: O.K.
TensorVariable<2U>  TensorVariable<2U>::Transposed() const
 {
      return TensorVariable( flag[0], flag[1], 
                             data[0][0], data[1][0],
                             data[0][1], data[1][1] );
 }


// re-tested: SKM 29-9-2001
// @test tested: O.K.
TensorVariable<2U> TensorVariable<2U>::Inverse() const
 {
    double64  det = Determinant();
    
    if ( det == static_cast<double64>(0.) ) {
         std::cerr <<"\nTensorVariable<2U>::Inverse: Determinant = 0" << std::endl;
         return TensorVariable<2U>();
      }
      
    det = 1. / det;  
      
      return TensorVariable( flag[0], flag[1], 
                             det * data[1][1], det * -data[0][1],
                             det * -data[1][0], det * data[0][0] );
 }
 
 
 
/// @test re-tested: SKM 29-9-2001
double64  TensorVariable<2U>::MinElement() const
 {
    double64 me = data[0][0];
     
    if ( data[0][1] < me ) me = data[0][1];
    if ( data[1][0] < me ) me = data[1][0];
    if ( data[1][1] < me ) me = data[1][1];

    return me;
 }
  
  
/// @test re-tested: SKM 29-9-2001

double64  TensorVariable<2U>::MaxElement() const
 {
    double64 me = data[0][0];

    if ( data[0][1] > me ) me = data[0][1];
    if ( data[1][0] > me ) me = data[1][0];
    if ( data[1][1] > me ) me = data[1][1];

    return me;
 }
 
 
 
/// @test re-tested: SKM 29-9-2001

bool  TensorVariable<2U>::IsWithinRange( double64 vmin, double64 vmax ) const
 {
    if ( data[0][0] < vmin || data[0][0] > vmax ) return false;
    if ( data[0][1] < vmin || data[0][1] > vmax ) return false;
    if ( data[1][0] < vmin || data[1][0] > vmax ) return false;
    if ( data[1][1] < vmin || data[1][1] > vmax ) return false; 

    return true;
 }
 
 
/// vector-matrix multiplication: v^T = (v^T * A)^T = A^T v
VectorVariable<2U>  operator*( const VectorVariable<2U>& vc, const TensorVariable<2U>& ts )
 {
    VectorVariable<2U> temp( vc.Flag(0), vc.Flag(1), 
                             ts(0,0) * vc[0] + ts(1,0) * vc[1],
                             ts(0,1) * vc[0] + ts(1,1) * vc[1] );
    return temp;
 }


 
// INPUT OUTPUT


void  TensorVariable<2U>::In()
 {
     cout.flush();
     cout <<"\nEnter ["<< 2U << "] tensor variable status: ";
     cout.flush();
     string  status;
     cin >> status;
     for ( size_t i=0; i<2U; i++ )
       flag[i] = parseStatus( status.c_str() );
     
     cout <<"\nEnter first row of elements : ";
     cout.flush();
     for ( size_t i=0; i<2U; i++ ) cin >> data[0][i];
     cout <<"Enter second row of elements: ";
     cout.flush();
     for ( size_t i=0; i<2U; i++ ) cin >> data[1][i];

 } // end In



/// @test tested: O.K.

void  TensorVariable<2U>::Out() const
 {
     size_t   i, j;
     
     cout <<"\nStatus: "<< endl;
    for ( j=0; j<2U; j++ ) cout << parseStatus(flag[j]) <<"  ";
    cout << endl;

     cout <<"\nValues: "<< endl;
     for ( i=0; i<2U; i++ )
       {
          for ( j=0; j<2U; j++ ) cout << data[i][j] <<"\t\t";
          cout << endl;
       }
 } // end Out



/**
 
Computes eigenvalues and eigenvectors assuming that the tensor variable
is symmetric. If not symmetric, off-diagonal elements are averaged.
No check of symmetry is performed.  

By default the length of the eigenvectors is equivalent to the eigenvalues
but it can be normalized to one.  

@section arguments Input Arguments

Flag to normalize length of eigenvectors.  

@return Eigenvalues and vectors are returned into the supplied Vector and
TensorVariables, respectively.

Will return false if the rank of the matrix is zero.  

@section implementation Implementation 

Just uses the normal binomial root formula.  
*/

bool TensorVariable<2U>::Eigen( VectorVariable<2U>& evals, 
                                TensorVariable<2U>& evecs, 
                                bool bNormalize ) const
{
	//calculate eigenvalues
	//set eigenvalues
	if(!EigenValues(evals)) return false;
		
	//suppose tensor is symmetric, by averaging the non diagonal elements
	const double64 f_data01((data[0][1]+data[1][0])/2.);
	if(f_data01 == 0) 
	{
	  evecs.Identity();
      return true;
	}
  
  VectorVariable<2U> v;	
	v(0)= -(data[1][1]-evals[0])/f_data01;
	v(1)=1.;
	
	if(bNormalize)
	  v.EuclideanNormalize();
	
	evecs.AssignToRow(0,v);
	
	v(0)= -(data[1][1]-evals[1])/f_data01;
	v(1)=1.;
	
	if(bNormalize)
	  v.EuclideanNormalize();
	
	evecs.AssignToRow(1,v);
	
	return true;
}



bool TensorVariable<2U>::EigenValues( VectorVariable<2U>& vecEigenvalues ) const
{
   // a0 and a1  are coefficients of quadratic equation " x^2 + a1 * x + a0 = 0". The roots of this equation are 
   // the Eigen values of the tensor
   const double64 a1 = -( data[ 0 ][ 0 ] + data[ 1 ][ 1 ] );
   const double64 a0 = ( data[ 0 ][ 0 ] * data[ 1 ][ 1 ] -
                         data[ 1 ][ 0 ] * data[ 0 ][ 1 ] );
   const double64 D = a1 * a1 - 4 * a0;
   
   if ( D >= 0 ) // The equation has two real roots
   {
      vecEigenvalues( 0 ) = ( -a1 + std::sqrt( D ) ) / 2.0;
      vecEigenvalues( 1 ) = ( -a1 - std::sqrt( D ) ) / 2.0;
   }
   
   if ( D < 0 ) // The equation has imaginary root
   {
      vecEigenvalues( 0 ) = -1 * std::sqrt( -1. );
      vecEigenvalues( 1 ) = -1 * std::sqrt( -1. );
   }
     
   return true;   
}




bool TensorVariable<2U>::EigenValues( vector<double64>& vecEigenvalues ) const
{
   // a0 and a1  are coefficients of quadratic equation " x^2 + a1 * x + a0 = 0". The roots of this equation are 
   // the Eigen values of the tensor
   const double64 a1 = -( data[ 0 ][ 0 ] + data[ 1 ][ 1 ] );
   const double64 a0 = ( data[ 0 ][ 0 ] * data[ 1 ][ 1 ] -
                         data[ 1 ][ 0 ] * data[ 0 ][ 1 ] );
   const double64 D = a1 * a1 - 4 * a0;
   
   if ( D >= 0 ) // The equation has two real roots
   {
      vecEigenvalues[ 0 ] = ( -a1 + std::sqrt( D ) ) / 2.0;
      vecEigenvalues[ 1 ] = ( -a1 - std::sqrt( D ) ) / 2.0;
   }
   
   if ( D < 0 ) // The equation has imaginary root
   {
      vecEigenvalues[ 0 ] = -1 * std::sqrt( -1. );
      vecEigenvalues[ 1 ] = -1 * std::sqrt( -1. );
   }
     
   return true;   
}


bool TensorVariable<2U>::EigenNonSymmetric( VectorVariable<2U>& eigenVals,
                                            TensorVariable<2U>& eigenVecs ) const
{
    const int n = 2;
    double V[n][n], d[n], e[n];

    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        V[i][j] = (*this)(i,j);

    // Symmetric Householder reduction to tridiagonal form.
    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.
    for (int j = 0; j < n; j++) {
        d[j] = V[n-1][j];
      }

    // Householder reduction to tridiagonal form.

    for (int i = n-1; i > 0; i--) {

      // Scale to avoid under/overflow.

      double scale = 0.0;
      double h = 0.0;
      for (int k = 0; k < i; k++) {
        scale = scale + fabs(d[k]);
      }
      if (scale == 0.0) {
        e[i] = d[i-1];
        for (int j = 0; j < i; j++) {
          d[j] = V[i-1][j];
          V[i][j] = 0.0;
          V[j][i] = 0.0;
        }
      } else {

        // Generate Householder vector.

        for (int k = 0; k < i; k++) {
          d[k] /= scale;
          h += d[k] * d[k];
        }
        double f = d[i-1];
        double g = sqrt(h);
        if (f > 0) {
          g = -g;
        }
        e[i] = scale * g;
        h = h - f * g;
        d[i-1] = f - g;
        for (int j = 0; j < i; j++) {
          e[j] = 0.0;
        }

        // Apply similarity transformation to remaining columns.

        for (int j = 0; j < i; j++) {
          f = d[j];
          V[j][i] = f;
          g = e[j] + V[j][j] * f;
          for (int k = j+1; k <= i-1; k++) {
            g += V[k][j] * d[k];
            e[k] += V[k][j] * f;
          }
          e[j] = g;
        }
        f = 0.0;
        for (int j = 0; j < i; j++) {
          e[j] /= h;
          f += e[j] * d[j];
        }
        double hh = f / (h + h);
        for (int j = 0; j < i; j++) {
          e[j] -= hh * d[j];
        }
        for (int j = 0; j < i; j++) {
          f = d[j];
          g = e[j];
          for (int k = j; k <= i-1; k++) {
            V[k][j] -= (f * e[k] + g * d[k]);
          }
          d[j] = V[i-1][j];
          V[i][j] = 0.0;
        }
      }
      d[i] = h;
    }

    // Accumulate transformations.

    for (int i = 0; i < n-1; i++) {
      V[n-1][i] = V[i][i];
      V[i][i] = 1.0;
      double h = d[i+1];
      if (h != 0.0) {
        for (int k = 0; k <= i; k++) {
          d[k] = V[k][i+1] / h;
        }
        for (int j = 0; j <= i; j++) {
          double g = 0.0;
          for (int k = 0; k <= i; k++) {
            g += V[k][i+1] * V[k][j];
          }
          for (int k = 0; k <= i; k++) {
            V[k][j] -= g * d[k];
          }
        }
      }
      for (int k = 0; k <= i; k++) {
        V[k][i+1] = 0.0;
      }
    }
    for (int j = 0; j < n; j++) {
      d[j] = V[n-1][j];
      V[n-1][j] = 0.0;
    }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;

  // Symmetric tridiagonal QL algorithm.
  //  This is derived from the Algol procedures tql2, by
  //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  //  Fortran subroutine in EISPACK.

    for (int i = 1; i < n; i++) {
      e[i-1] = e[i];
    }
    e[n-1] = 0.0;

    double f = 0.0;
    double tst1 = 0.0;
    double eps = pow(2.0,-52.0);
    for (int l = 0; l < n; l++) {

      // Find small subdiagonal element

      tst1 = std::max(tst1,fabs(d[l]) + fabs(e[l]));
      int m = l;
      while (m < n) {
        if (fabs(e[m]) <= eps*tst1) {
          break;
        }
        m++;
      }

      // If m == l, d[l] is an eigenvalue,
      // otherwise, iterate.

      if (m > l) {
        int iter = 0;
        do {
          iter = iter + 1;  // (Could check iteration count here.)

          // Compute implicit shift

          double g = d[l];
          double p = (d[l+1] - g) / (2.0 * e[l]);
          double r = sqrt(p*p + 1.0);
          if (p < 0) {
            r = -r;
          }
          d[l] = e[l] / (p + r);
          d[l+1] = e[l] * (p + r);
          double dl1 = d[l+1];
          double h = g - d[l];
          for (int i = l+2; i < n; i++) {
            d[i] -= h;
          }
          f = f + h;

          // Implicit QL transformation.

          p = d[m];
          double c = 1.0;
          double c2 = c;
          double c3 = c;
          double el1 = e[l+1];
          double s = 0.0;
          double s2 = 0.0;
          for (int i = m-1; i >= l; i--) {
            c3 = c2;
            c2 = c;
            s2 = s;
            g = c * e[i];
            h = c * p;
            r = sqrt(p*p + e[i]*e[i]);
            e[i+1] = s * r;
            s = e[i] / r;
            c = p / r;
            p = c * d[i] - s * g;
            d[i+1] = h + s * (c * g + s * d[i]);

            // Accumulate transformation.

            for (int k = 0; k < n; k++) {
              h = V[k][i+1];
              V[k][i+1] = s * V[k][i] + c * h;
              V[k][i] = c * V[k][i] - s * h;
            }
          }
          p = -s * s2 * c3 * el1 * e[l] / dl1;
          e[l] = s * p;
          d[l] = c * p;

          // Check for convergence.

        } while (fabs(e[l]) > eps*tst1);
      }
      d[l] = d[l] + f;
      e[l] = 0.0;
    }
    
    // Sort eigenvalues and corresponding vectors.

    for (int i = 0; i < n-1; i++) {
      int k = i;
      double p = d[i];
      for (int j = i+1; j < n; j++) {
        if (d[j] < p) {
          k = j;
          p = d[j];
        }
      }
      if (k != i) {
        d[k] = d[i];
        d[i] = p;
        for (int j = 0; j < n; j++) {
          p = V[j][i];
          V[j][i] = V[j][k];
          V[j][k] = p;
        }
      }
    }

    for (int i=0; i<n; i++)
    {
      for (int j=0; j<n; j++)
            eigenVecs(i,j) = V[i][j];
      eigenVals(i) = d[i];
    }

    return true;
}


} // end namespace csmp


