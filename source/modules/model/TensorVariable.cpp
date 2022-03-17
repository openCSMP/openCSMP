#include "TensorVariable.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {


TensorVariable<3U>::TensorVariable( const TensorVariable<3U>& t )
  : flag( t.flag ),
  data( t.data )
{
}


TensorVariable<3U>::TensorVariable()
  : flag{ { ANY,ANY,ANY } },
  data{ numeric_limits<double>::quiet_NaN(),numeric_limits<double>::quiet_NaN(),numeric_limits<double>::quiet_NaN(),
  numeric_limits<double>::quiet_NaN(),numeric_limits<double>::quiet_NaN(),numeric_limits<double>::quiet_NaN(),
  numeric_limits<double>::quiet_NaN(),numeric_limits<double>::quiet_NaN(),numeric_limits<double>::quiet_NaN() }
{
}



TensorVariable<3U>&  TensorVariable<3U>::operator=( const TensorVariable<3U>& ts )
{
  if ( &ts != this ) {
    flag = ts.flag;
    data = ts.data;
  }
  return *this;
}



/**
initialises variable as diagonal isotropic tensor with flag and value
SKM 17/6/2015
*/
TensorVariable<3U>::TensorVariable( VARIABLE_FLAG f, double val )
  : flag{ { f,f,f } },
  data{ { val,0.,0.,0.,val,0.,0.,0.,val } }
{
}

// ******************************************************************************************
//
//            INLINE METHODS START HERE
//
// ******************************************************************************************

TensorVariable<3U>::~TensorVariable() {}


double& TensorVariable<3U>::operator()( uint32_t i, uint32_t j )
{
#ifndef NDEBUG 
  if ( i >= 3U ) {
    std::cerr << "\nTensorVariable<3U>::operator(): row access violation, i=" << i << std::endl;
    return data[0][0];
  }
  if ( j >= 3U ) {
    std::cerr << "\nTensorVariable<3U>::operator(): column access violation, j=" << j << std::endl;
    return data[0][0];
  }
#endif
  return data[i][j];
}




const double& TensorVariable<3U>::operator()( uint32_t i, uint32_t j ) const
{
#ifndef NDEBUG 
  if ( i >= 3U ) {
    std::cerr << "\nTensorVariable<3U>::operator(): row access violation, i=" << i << std::endl;
    return data[0][0];
  }
  if ( j >= 3U ) {
    std::cerr << "\nTensorVariable<3U>::operator(): column access violation, j=" << j << std::endl;
    return data[0][0];
  }
#endif
  return data[i][j];
}



void TensorVariable<3U>::Component( uint32_t i, double val )
{
  assert( i < Size() );
  // row by row
  if ( i == 0U )      data[0U][0U] = val;
  else if ( i == 1U ) data[0U][1U] = val;
  else if ( i == 2U ) data[0U][2U] = val;

  else if ( i == 3U ) data[1U][0U] = val;
  else if ( i == 4U ) data[1U][1U] = val;
  else if ( i == 5U ) data[1U][2U] = val;

  else if ( i == 6U ) data[2U][0U] = val;
  else if ( i == 7U ) data[2U][1U] = val;

  else data[2U][2U] = val; // remaining case
}




double TensorVariable<3U>::Component( uint32_t i ) const
{
  assert( i < Size() );
  // row by row
  if ( i == 0U ) return data[0U][0U];
  if ( i == 1U ) return data[0U][1U];
  if ( i == 2U ) return data[0U][2U];

  if ( i == 3U ) return data[1U][0U];
  if ( i == 4U ) return data[1U][1U];
  if ( i == 5U ) return data[1U][2U];

  if ( i == 6U ) return data[2U][0U];
  if ( i == 7U ) return data[2U][1U];

  return data[2U][2U]; // remaining case
}




VARIABLE_FLAG& TensorVariable<3U>::Flag( uint32_t i )
{
#ifndef NDEBUG 
  if ( i >= 3U ) {
    std::cerr << "\nTensorVariable<3U>::Flag(): diagonal access violation, i=" << i << std::endl;
    return flag[0];
  }
#endif
  return flag[i];
}



VARIABLE_FLAG  TensorVariable<3U>::Flag( uint32_t i ) const
{
#ifndef NDEBUG 
  if ( i >= 3U ) {
    std::cerr << "\nTensorVariable<3U>::Flag(): diagonal access violation, i=" << i << std::endl;
    return flag[0];
  }
#endif
  return flag[i];
}




bool TensorVariable<3U>::EigenValues( VectorVariable<3U>& Ev ) const
{
  return EigenValuesPositiveDefiniteSymmetricMatrix( Ev( 0 ), Ev( 1 ), Ev( 2 ) );
}


bool TensorVariable<3U>::EigenValues( std::vector<double>& Ev ) const
{
  return EigenValuesPositiveDefiniteSymmetricMatrix( Ev[0], Ev[1], Ev[2] );
}






double TensorVariable<3U>::Trace() const
{
  return data[0][0] + data[1][1] + data[2][2];
}



void TensorVariable<3U>::AssignToRow( uint32_t iRow, VectorVariable<3U>& vc )
{
  if ( iRow == 0U )
    flag[0U] = vc.Flag( 0U );
  else if ( iRow == 1U )
    flag[1U] = vc.Flag( 1U );
  else
    flag[2U] = vc.Flag( 2U );

  data[iRow][0U] = vc[0U];
  data[iRow][1U] = vc[1U];
  data[iRow][2U] = vc[2U];
}


void TensorVariable<3U>::AssignToColumn( uint32_t iCol, VectorVariable<3U>& vc )
{
  if ( iCol == 0U )
    flag[0U] = vc.Flag( 0U );
  else if ( iCol == 1U )
    flag[1U] = vc.Flag( 1U );
  else
    flag[2U] = vc.Flag( 2U );

  data[0U][iCol] = vc[0U];
  data[1U][iCol] = vc[1U];
  data[2U][iCol] = vc[2U];
}




VectorVariable<3U> TensorVariable<3U>::Row( uint32_t iRow ) const
{
  return VectorVariable<3U>( flag[iRow], flag[iRow], flag[iRow],
                             data[iRow][0U], data[iRow][1U], data[iRow][2U] );
}


VectorVariable<3U> TensorVariable<3U>::Column( uint32_t iCol ) const
{
  return VectorVariable<3U>( flag[iCol], flag[iCol], flag[iCol],
                             data[0U][iCol], data[1U][iCol], data[2U][iCol] );
}


bool  TensorVariable<3U>::operator==( const TensorVariable<3U>& ts ) const
{
  return (data == ts.data && flag == ts.flag);
}


bool  TensorVariable<3U>::operator!=( const TensorVariable<3U>& t ) const
{
  return !(*this == t);
}


/// to achieve reasonable ordering in associative STL containers
bool  TensorVariable<3U>::operator<( const TensorVariable<3U>& t ) const
{
  return (this < &t);
}

bool TensorVariable<3U>::Out( std::fstream& fp ) const
{
  const int32_t flag_0( this->flag[0] );
  const int32_t flag_1( this->flag[1] );
  const int32_t flag_2( this->flag[2] );
  const size_t flag_size = sizeof( int32_t );
  fp.write( (char*)&flag_0, flag_size );
  fp.write( (char*)&flag_1, flag_size );
  fp.write( (char*)&flag_2, flag_size );
  const size_t data_size = sizeof( double );
  fp.write( (char*)&this->data[0][0], data_size );
  fp.write( (char*)&this->data[1][0], data_size );
  fp.write( (char*)&this->data[2][0], data_size );
  fp.write( (char*)&this->data[0][1], data_size );
  fp.write( (char*)&this->data[1][1], data_size );
  fp.write( (char*)&this->data[2][1], data_size );
  fp.write( (char*)&this->data[0][2], data_size );
  fp.write( (char*)&this->data[1][2], data_size );
  fp.write( (char*)&this->data[2][2], data_size );
  return true;
}

bool TensorVariable<3U>::In( std::fstream& fp )
{
  const size_t flag_size = sizeof( int32_t );
  fp.read( (char*)&this->flag[0], flag_size );
  fp.read( (char*)&this->flag[1], flag_size );
  fp.read( (char*)&this->flag[2], flag_size );
  const size_t data_size = sizeof( double );
  fp.read( (char*)&this->data[0][0], data_size );
  fp.read( (char*)&this->data[1][0], data_size );
  fp.read( (char*)&this->data[2][0], data_size );
  fp.read( (char*)&this->data[0][1], data_size );
  fp.read( (char*)&this->data[1][1], data_size );
  fp.read( (char*)&this->data[2][1], data_size );
  fp.read( (char*)&this->data[0][2], data_size );
  fp.read( (char*)&this->data[1][2], data_size );
  fp.read( (char*)&this->data[2][2], data_size );
  return true;
}


/// fastest way to insert a tensor into an STL container; tensor only has flags for diagonal elements
TensorVariable<2U> makeTensor( VARIABLE_FLAG f1, VARIABLE_FLAG f2,
                               double v11, double v12,
                               double v21, double v22 )
{
  return TensorVariable<2U>( f1, f2, v11, v12, v21, v22 );
}


/// fastest way to insert a tensor into an STL container; tensor only has flags for diagonal elements
TensorVariable<3U> makeTensor( VARIABLE_FLAG f1, VARIABLE_FLAG f2, VARIABLE_FLAG f3,
                               double v11, double v12, double v13,
                               double v21, double v22, double v23,
                               double v31, double v32, double v33 )
{
  return TensorVariable<3U>( f1, f2, f3, v11, v12, v13, v21, v22, v23, v31, v32, v33 );
}




TensorVariable<3U>::TensorVariable( VARIABLE_FLAG f,
                                    double v11, double v12, double v13,
                                    double v21, double v22, double v23,
                                    double v31, double v32, double v33 )
  : flag{ { f,f,f } },
  data{ { v11,v12,v13,v21,v22,v23,v31,v32,v33 } }
{
}



TensorVariable<3U>::TensorVariable( const VARIABLE_FLAG f11, const VARIABLE_FLAG f22, const VARIABLE_FLAG f33,
                                    const double  v11, const double  v12, const double  v13,
                                    const double  v21, const double  v22, const double  v23,
                                    const double  v31, const double  v32, const double  v33 )
  : flag{ { f11,f22,f33 } },
  data{ { v11,v12,v13,v21,v22,v23,v31,v32,v33 } }
{
}



// here the flag of the lefthand tensor-variable is sustained

TensorVariable<3U>  TensorVariable<3U>::operator+( const TensorVariable<3U>& t ) const
{
  return TensorVariable( flag[0], flag[1], flag[2],
                    t.data[0][0] + data[0][0], t.data[0][1] + data[0][1], t.data[0][2] + data[0][2],
                    t.data[1][0] + data[1][0], t.data[1][1] + data[1][1], t.data[1][2] + data[1][2],
                    t.data[2][0] + data[2][0], t.data[2][1] + data[2][1], t.data[2][2] + data[2][2] );
}



TensorVariable<3U>  TensorVariable<3U>::operator-( const TensorVariable<3U>& t ) const
{
  return TensorVariable( flag[0], flag[1], flag[2],
                    data[0][0] - t.data[0][0], data[0][1] - t.data[0][1], data[0][2] - t.data[0][2],
                    data[1][0] - t.data[1][0], data[1][1] - t.data[1][1], data[1][2] - t.data[1][2],
                    data[2][0] - t.data[2][0], data[2][1] - t.data[2][1], data[2][2] - t.data[2][2] );
}


TensorVariable<3U>  TensorVariable<3U>::operator+( double val ) const
{
  return TensorVariable( flag[0], flag[1], flag[2],
                    data[0][0] + val, data[0][1] + val, data[0][2] + val,
                    data[1][0] + val, data[1][1] + val, data[1][2] + val,
                    data[2][0] + val, data[2][1] + val, data[2][2] + val );
}



TensorVariable<3U>  TensorVariable<3U>::operator-( double val ) const
{
  return TensorVariable( flag[0], flag[1], flag[2],
                    data[0][0] - val, data[0][1] - val, data[0][2] - val,
                    data[1][0] - val, data[1][1] - val, data[1][2] - val,
                    data[2][0] - val, data[2][1] - val, data[2][2] - val );
}



TensorVariable<3U>  TensorVariable<3U>::operator*( double val ) const
{
  return TensorVariable( flag[0], flag[1], flag[2],
                    data[0][0] * val, data[0][1] * val, data[0][2] * val,
                    data[1][0] * val, data[1][1] * val, data[1][2] * val,
                    data[2][0] * val, data[2][1] * val, data[2][2] * val );
}



TensorVariable<3U>  TensorVariable<3U>::operator/( double val ) const
{
  return TensorVariable( flag[0], flag[1], flag[2],
                    data[0][0] / val, data[0][1] / val, data[0][2] / val,
                    data[1][0] / val, data[1][1] / val, data[1][2] / val,
                    data[2][0] / val, data[2][1] / val, data[2][2] / val );
}


// matrix vector multiplication: v = M * v
/// @test tested: O.K. SKM 29-9-2001

VectorVariable<3U>  TensorVariable<3U>::operator*( const VectorVariable<3U>& vc ) const
{
  VectorVariable<3U> temp( vc.Flag( 0 ), vc.Flag( 1 ), vc.Flag( 2 ),
                           data[0][0] * vc[0] + data[0][1] * vc[1] + data[0][2] * vc[2],
                           data[1][0] * vc[0] + data[1][1] * vc[1] + data[1][2] * vc[2],
                           data[2][0] * vc[0] + data[2][1] * vc[1] + data[2][2] * vc[2] );
  return temp;
}

Point<3U>  TensorVariable<3U>::operator*( const Point<3U>& v ) const
{
  return Point<3U>(
    data[0][0] * v[0] + data[0][1] * v[1] + data[0][2] * v[2],
    data[1][0] * v[0] + data[1][1] * v[1] + data[1][2] * v[2],
    data[2][0] * v[0] + data[2][1] * v[1] + data[2][2] * v[2] );
}



// re-tested: SKM 29-9-2001
TensorVariable<3U> TensorVariable<3U>::Adjoint() const
{
  return TensorVariable( flag[0], flag[1], flag[2],
                    data[1][1] * data[2][2] - data[1][2] * data[2][1],
                    -data[1][0] * data[2][2] + data[1][2] * data[2][0],
                    data[1][0] * data[2][1] - data[2][0] * data[1][1],
                    -data[0][1] * data[2][2] + data[0][2] * data[2][1],
                    data[0][0] * data[2][2] - data[0][2] * data[2][0],
                    -data[0][0] * data[2][1] + data[0][1] * data[2][0],
                    data[0][1] * data[1][2] - data[0][2] * data[1][1],
                    -data[0][0] * data[1][2] + data[0][2] * data[1][0],
                    data[0][0] * data[1][1] - data[0][1] * data[1][0] );
}



/// @test tested: O.K.
TensorVariable<3U>  TensorVariable<3U>::operator*( const TensorVariable<3U>& ts ) const
{
  TensorVariable<3U> temp;

  // Unrolled loop(i...3,j...3), expression: data[i][j] = ts.data[j][j]
  // row 0
  temp.data[0][0] = data[0][0] * ts.data[0][0] + data[0][1] * ts.data[1][0] + data[0][2] * ts.data[2][0];
  temp.data[0][1] = data[0][0] * ts.data[0][1] + data[0][1] * ts.data[1][1] + data[0][2] * ts.data[2][1];
  temp.data[0][2] = data[0][0] * ts.data[0][2] + data[0][1] * ts.data[1][2] + data[0][2] * ts.data[2][2];
  // row 1
  temp.data[1][0] = data[1][0] * ts.data[0][0] + data[1][1] * ts.data[1][0] + data[1][2] * ts.data[2][0];
  temp.data[1][1] = data[1][0] * ts.data[0][1] + data[1][1] * ts.data[1][1] + data[1][2] * ts.data[2][1];
  temp.data[1][2] = data[1][0] * ts.data[0][2] + data[1][1] * ts.data[1][2] + data[1][2] * ts.data[2][2];
  // row 2
  temp.data[2][0] = data[2][0] * ts.data[0][0] + data[2][1] * ts.data[1][0] + data[2][2] * ts.data[2][0];
  temp.data[2][1] = data[2][0] * ts.data[0][1] + data[2][1] * ts.data[1][1] + data[2][2] * ts.data[2][1];
  temp.data[2][2] = data[2][0] * ts.data[0][2] + data[2][1] * ts.data[1][2] + data[2][2] * ts.data[2][2];

  temp.flag[0] = flag[0];
  temp.flag[1] = flag[1];
  temp.flag[2] = flag[2];

  return temp;
}




TensorVariable<3U>&  TensorVariable<3U>::operator+=( const ScalarVariable& sc )
{
  data[0][0] += sc();
  data[0][1] += sc();
  data[1][0] += sc();
  data[1][1] += sc();
  data[0][2] += sc();
  data[1][2] += sc();
  data[2][0] += sc();
  data[2][1] += sc();
  data[2][2] += sc();

  return *this;
}



TensorVariable<3U>&  TensorVariable<3U>::operator-=( const ScalarVariable& sc )
{
  data[0][0] -= sc();
  data[0][1] -= sc();
  data[1][0] -= sc();
  data[1][1] -= sc();
  data[0][2] -= sc();
  data[1][2] -= sc();
  data[2][0] -= sc();
  data[2][1] -= sc();
  data[2][2] -= sc();

  return *this;
}



TensorVariable<3U>&  TensorVariable<3U>::operator*=( const ScalarVariable& sc )
{
  data[0][0] *= sc();
  data[0][1] *= sc();
  data[1][0] *= sc();
  data[1][1] *= sc();
  data[0][2] *= sc();
  data[1][2] *= sc();
  data[2][0] *= sc();
  data[2][1] *= sc();
  data[2][2] *= sc();

  return *this;
}



TensorVariable<3U>&  TensorVariable<3U>::operator/=( const ScalarVariable& sc )
{
  data[0][0] /= sc();
  data[0][1] /= sc();
  data[1][0] /= sc();
  data[1][1] /= sc();
  data[0][2] /= sc();
  data[1][2] /= sc();
  data[2][0] /= sc();
  data[2][1] /= sc();
  data[2][2] /= sc();

  return *this;
}



TensorVariable<3U>&  TensorVariable<3U>::operator+=( const TensorVariable<3U>& ts )
{
  data[0][0] += ts.data[0][0];
  data[0][1] += ts.data[0][1];
  data[1][0] += ts.data[1][0];
  data[1][1] += ts.data[1][1];
  data[0][2] += ts.data[0][2];
  data[1][2] += ts.data[1][2];
  data[2][0] += ts.data[2][0];
  data[2][1] += ts.data[2][1];
  data[2][2] += ts.data[2][2];

  return *this;
}



TensorVariable<3U>&  TensorVariable<3U>::operator-=( const TensorVariable<3U>& ts )
{
  data[0][0] -= ts.data[0][0];
  data[0][1] -= ts.data[0][1];
  data[1][0] -= ts.data[1][0];
  data[1][1] -= ts.data[1][1];
  data[0][2] -= ts.data[0][2];
  data[1][2] -= ts.data[1][2];
  data[2][0] -= ts.data[2][0];
  data[2][1] -= ts.data[2][1];
  data[2][2] -= ts.data[2][2];

  return *this;
}


/// element by element division
TensorVariable<3U>&  TensorVariable<3U>::operator/=( const TensorVariable<3U>& ts )
{
  data[0][0] /= ts.data[0][0];
  data[0][1] /= ts.data[0][1];
  data[1][0] /= ts.data[1][0];
  data[1][1] /= ts.data[1][1];
  data[0][2] /= ts.data[0][2];
  data[1][2] /= ts.data[1][2];
  data[2][0] /= ts.data[2][0];
  data[2][1] /= ts.data[2][1];
  data[2][2] /= ts.data[2][2];

  return *this;
}



/// element by element division
TensorVariable<3U>  TensorVariable<3U>::operator/( const TensorVariable<3U>& ts ) const
{
  TensorVariable  temp;

  temp.data[0][0] = data[0][0] / ts.data[0][0];
  temp.data[0][1] = data[0][1] / ts.data[0][1];
  temp.data[1][0] = data[1][0] / ts.data[1][0];
  temp.data[1][1] = data[1][1] / ts.data[1][1];
  temp.data[0][2] = data[0][2] / ts.data[0][2];
  temp.data[1][2] = data[1][2] / ts.data[1][2];
  temp.data[2][0] = data[2][0] / ts.data[2][0];
  temp.data[2][1] = data[2][1] / ts.data[2][1];
  temp.data[2][2] = data[2][2] / ts.data[2][2];

  temp.flag[0] = flag[0];
  temp.flag[1] = flag[1];
  temp.flag[2] = flag[2];

  return temp;
}



/// matrix multiplication
TensorVariable<3U>&  TensorVariable<3U>::operator*=( const TensorVariable<3U>& ts )
{
  TensorVariable<3U> temp;

  // Unrolled loop(i...3,j...3), expression: data[i][j] = ts.data[j][j]
  // row 0
  temp.data[0][0] = data[0][0] * ts.data[0][0] + data[0][1] * ts.data[1][0] + data[0][2] * ts.data[2][0];
  temp.data[0][1] = data[0][0] * ts.data[0][1] + data[0][1] * ts.data[1][1] + data[0][2] * ts.data[2][1];
  temp.data[0][2] = data[0][0] * ts.data[0][2] + data[0][1] * ts.data[1][2] + data[0][2] * ts.data[2][2];
  // row 1
  temp.data[1][0] = data[1][0] * ts.data[0][0] + data[1][1] * ts.data[1][0] + data[1][2] * ts.data[2][0];
  temp.data[1][1] = data[1][0] * ts.data[0][1] + data[1][1] * ts.data[1][1] + data[1][2] * ts.data[2][1];
  temp.data[1][2] = data[1][0] * ts.data[0][2] + data[1][1] * ts.data[1][2] + data[1][2] * ts.data[2][2];
  // row 2
  temp.data[2][0] = data[2][0] * ts.data[0][0] + data[2][1] * ts.data[1][0] + data[2][2] * ts.data[2][0];
  temp.data[2][1] = data[2][0] * ts.data[0][1] + data[2][1] * ts.data[1][1] + data[2][2] * ts.data[2][1];
  temp.data[2][2] = data[2][0] * ts.data[0][2] + data[2][1] * ts.data[1][2] + data[2][2] * ts.data[2][2];

  temp.flag[0] = flag[0];
  temp.flag[1] = flag[1];
  temp.flag[2] = flag[2];

  return *this = std::move( temp );
}




TensorVariable<3U>&  TensorVariable<3U>::operator+=( double val )
{
  data[0][0] += val;
  data[0][1] += val;
  data[1][0] += val;
  data[1][1] += val;
  data[0][2] += val;
  data[1][2] += val;
  data[2][0] += val;
  data[2][1] += val;
  data[2][2] += val;

  return *this;
}



TensorVariable<3U>&  TensorVariable<3U>::operator-=( double val )
{
  data[0][0] -= val;
  data[0][1] -= val;
  data[1][0] -= val;
  data[1][1] -= val;
  data[0][2] -= val;
  data[1][2] -= val;
  data[2][0] -= val;
  data[2][1] -= val;
  data[2][2] -= val;

  return *this;
}



TensorVariable<3U>&  TensorVariable<3U>::operator*=( double val )
{
  data[0][0] *= val;
  data[0][1] *= val;
  data[1][0] *= val;
  data[1][1] *= val;
  data[0][2] *= val;
  data[1][2] *= val;
  data[2][0] *= val;
  data[2][1] *= val;
  data[2][2] *= val;

  return *this;
}



TensorVariable<3U>&  TensorVariable<3U>::operator/=( double val )
{
  data[0][0] /= val;
  data[0][1] /= val;
  data[1][0] /= val;
  data[1][1] /= val;
  data[0][2] /= val;
  data[1][2] /= val;
  data[2][0] /= val;
  data[2][1] /= val;
  data[2][2] /= val;

  return *this;
}


// --------------------
// ASSIGNMENT OPERATORS
// --------------------


TensorVariable<3U>&  TensorVariable<3U>::operator=( double val )
{
  data[0][0] = val;
  data[0][1] = val;
  data[1][0] = val;
  data[1][1] = val;
  data[0][2] = val;
  data[1][2] = val;
  data[2][0] = val;
  data[2][1] = val;
  data[2][2] = val;

  return *this;
}


// the flag is adopted from the scalar variable

TensorVariable<3U>&  TensorVariable<3U>::operator=( const ScalarVariable& sc )
{
  flag[0] = flag[1] = flag[2] = sc.Flag();

  data[0][0] = sc();
  data[0][1] = sc();
  data[1][0] = sc();
  data[1][1] = sc();
  data[0][2] = sc();
  data[1][2] = sc();
  data[2][0] = sc();
  data[2][1] = sc();
  data[2][2] = sc();

  return *this;
}


/// writes vector into the diagonal of the zero'd tensor
TensorVariable<3U>&  TensorVariable<3U>::operator=( const VectorVariable<3U>& vc )
{
  flag[0] = vc.Flag( 0 );
  flag[1] = vc.Flag( 1 );
  flag[2] = vc.Flag( 2 );

  data[0][0] = vc( 0 );
  data[0][1] = static_cast<double>(0.0);
  data[1][0] = static_cast<double>(0.0);
  data[1][1] = vc( 1 );
  data[0][2] = static_cast<double>(0.0);
  data[1][2] = static_cast<double>(0.0);
  data[2][0] = static_cast<double>(0.0);
  data[2][1] = static_cast<double>(0.0);
  data[2][2] = vc( 2 );

  return *this;
}








// -------
// METHODS
// -------

/// @test tested: O.K.

void TensorVariable<3U>::Identity()
{
  data[0][0] = static_cast<double>(1.0);
  data[0][1] = static_cast<double>(0.0);
  data[0][2] = static_cast<double>(0.0);

  data[1][0] = static_cast<double>(0.0);
  data[1][1] = static_cast<double>(1.0);
  data[1][2] = static_cast<double>(0.0);

  data[2][0] = static_cast<double>(0.0);
  data[2][1] = static_cast<double>(0.0);
  data[2][2] = static_cast<double>(1.0);
}



void TensorVariable<3U>::DiagonalValues( double f_00, double f_11, double f_22 )
{
  data[0][0] = f_00;
  data[1][1] = f_11;
  data[2][2] = f_22;
}



void TensorVariable<3U>::DiagonalValues( const vector<double>& vecDiags )
{
  data[0][0] = vecDiags[0];
  data[1][1] = vecDiags[1];
  data[2][2] = vecDiags[2];
}



void TensorVariable<3U>::DiagonalValues( const VectorVariable<3U>& vecDiags )
{
  data[0][0] = vecDiags[0];
  data[1][1] = vecDiags[1];
  data[2][2] = vecDiags[2];
  flag[0] = vecDiags.Flag( 0 );
  flag[1] = vecDiags.Flag( 1 );
  flag[2] = vecDiags.Flag( 2 );
}


/// @test tested: O.K.

TensorVariable<3U>  TensorVariable<3U>::Transposed() const
{
  return ( TensorVariable( flag[0], flag[1], flag[2],
            data[0][0], data[1][0], data[2][0],
            data[0][1], data[1][1], data[2][1],
            data[0][2], data[1][2], data[2][2] ) );
}


// re-tested: SKM 29-9-2001
/// @test tested: O.K.

double TensorVariable<3U>::Determinant() const
{
  double det = data[0][0] * (data[1][1] * data[2][2] - data[2][1] * data[1][2]);
  det -= data[0][1] * (data[1][0] * data[2][2] - data[2][0] * data[1][2]);
  det += data[0][2] * (data[1][0] * data[2][1] - data[2][0] * data[1][1]);

  return det;
}


// re-tested: SKM 29-9-2001
/// @test tested: O.K.

TensorVariable<3U> TensorVariable<3U>::Inverse() const
{
  double  det = Determinant();

  if ( det == static_cast<double>(0.) ) {
    std::cerr << "\nTensorVariable<dim>::Inverse: Determinant = 0" << std::endl;
    return TensorVariable<3U>();
  }

  det = 1. / det;

  // A^-1 = (1 / det A) B^T
  return TensorVariable( flag[0], flag[1], flag[2],
                         // B00 
                         det * (data[1][1] * data[2][2] - data[1][2] * data[2][1]),
                         // B10
                         det * (-data[0][1] * data[2][2] + data[0][2] * data[2][1]),
                         // B20
                         det * (data[0][1] * data[1][2] - data[0][2] * data[1][1]),
                         // B01
                         det * (-data[1][0] * data[2][2] + data[1][2] * data[2][0]),
                         // B11
                         det * (data[0][0] * data[2][2] - data[0][2] * data[2][0]),
                         // B21
                         det * (-data[0][0] * data[1][2] + data[0][2] * data[1][0]),
                         // B02
                         det * (data[1][0] * data[2][1] - data[2][0] * data[1][1]),
                         // B12
                         det * (-data[0][0] * data[2][1] + data[0][1] * data[2][0]),
                         // B22
                         det * (data[0][0] * data[1][1] - data[0][1] * data[1][0]) );
}





/// @test re-tested: SKM 29-9-2001

double  TensorVariable<3U>::MinElement() const
{
  double me = data[0][0];

  if ( data[0][1] < me ) me = data[0][1];
  if ( data[1][0] < me ) me = data[1][0];
  if ( data[1][1] < me ) me = data[1][1];
  if ( data[0][2] < me ) me = data[0][2];
  if ( data[1][2] < me ) me = data[1][2];
  if ( data[2][0] < me ) me = data[2][0];
  if ( data[2][1] < me ) me = data[2][1];
  if ( data[2][2] < me ) me = data[2][2];

  return me;
}


/// @test re-tested: SKM 29-9-2001

double  TensorVariable<3U>::MaxElement() const
{
  double me = data[0][0];

  if ( data[0][1] > me ) me = data[0][1];
  if ( data[1][0] > me ) me = data[1][0];
  if ( data[1][1] > me ) me = data[1][1];
  if ( data[0][2] > me ) me = data[0][2];
  if ( data[1][2] > me ) me = data[1][2];
  if ( data[2][0] > me ) me = data[2][0];
  if ( data[2][1] > me ) me = data[2][1];
  if ( data[2][2] > me ) me = data[2][2];

  return me;
}



/// @test re-tested: SKM 29-9-2001

bool  TensorVariable<3U>::IsWithinRange( double vmin, double vmax ) const
{
  if ( data[0][0] < vmin || data[0][0] > vmax ) return false;
  //if ( data[0][1] < vmin || data[0][1] > vmax ) return false;
  //if ( data[0][2] < vmin || data[0][2] > vmax ) return false;
  //if ( data[1][0] < vmin || data[1][0] > vmax ) return false; 
  if ( data[1][1] < vmin || data[1][1] > vmax ) return false;
  //if ( data[1][2] < vmin || data[1][2] > vmax ) return false; 
  //if ( data[2][0] < vmin || data[2][0] > vmax ) return false;
  //if ( data[2][1] < vmin || data[2][1] > vmax ) return false;
  if ( data[2][2] < vmin || data[2][2] > vmax ) return false;

  return true;
}




/// vector-matrix multiplication: v^T = (v^T * A)^T = A^T v
VectorVariable<3U>  operator*( const VectorVariable<3U>& vc, const TensorVariable<3U>& ts )
{
  return VectorVariable<3U>( vc.Flag( 0 ), vc.Flag( 1 ), vc.Flag( 2 ),
                             ts( 0, 0 ) * vc[0] + ts( 1, 0 ) * vc[1] + ts( 2, 0 ) * vc[2],
                             ts( 0, 1 ) * vc[0] + ts( 1, 1 ) * vc[1] + ts( 2, 1 ) * vc[2],
                             ts( 0, 2 ) * vc[0] + ts( 1, 2 ) * vc[1] + ts( 2, 2 ) * vc[2] );
}

/// vector-matrix multiplication: v^T = (v^T * A)^T = A^T v
Point<3U>  operator*( const Point<3U>& vc, const TensorVariable<3U>& ts )
{
  return Point<3u>(
    ts( 0, 0 ) * vc[0] + ts( 1, 0 ) * vc[1] + ts( 2, 0 ) * vc[2],
    ts( 0, 1 ) * vc[0] + ts( 1, 1 ) * vc[1] + ts( 2, 1 ) * vc[2],
    ts( 0, 2 ) * vc[0] + ts( 1, 2 ) * vc[1] + ts( 2, 2 ) * vc[2] );
}



template<uint32_t dim>
ostream&  operator<<( ostream& stream, const TensorVariable<dim>& o )
{
  for ( auto i = 0U; i<dim; i++ )
  {
    stream << "flag" << i + 1U << ": " << parseStatus( o.Flag( i ) );
    stream << ", row" << i + 1U << ":";
    for ( uint32_t j = 0U; j<dim; j++ )
      stream <<" "<< o( i, j );
    stream <<"; ";

  }

  return stream;
}



void  TensorVariable<3U>::In()
{
  size_t  i;
  string  status;

  cout.flush();
  cout << "\nEnter [" << 3U << "] tensor variable status: ";
  cout.flush();
  cin >> status;
  for ( i = 0; i<3U; i++ )
    flag[i] = parseStatus( status.c_str() );

  cout << "\nEnter first row of elements : ";
  cout.flush();
  for ( i = 0; i<3U; i++ ) cin >> data[0][i];
  cout << "Enter second row of elements: ";
  cout.flush();
  for ( i = 0; i<3U; i++ ) cin >> data[1][i];
  cout << "Enter third row of elements : ";
  cout.flush();
  for ( i = 0; i<3U; i++ ) cin >> data[2][i];

} // end In





  /// @test tested: O.K.
void  TensorVariable<3U>::Out() const
{
  size_t   i, j;

  cout << "\nStatus: " << endl;
  for ( i = 0; i<3U; i++ )
    cout << parseStatus( flag[i] ) << "  ";

  cout << endl;

  cout << "\nValues: " << endl;
  for ( i = 0; i<3U; i++ )
  {
    for ( j = 0; j<3U; j++ ) cout << data[i][j] << "\t\t";
    cout << endl;
  }
} // end Out






  /**

  Computes eigenvalues and eigenvectors assuming that the tensor variable
  is symmetric. If not symmetric, off-diagonal elements are averaged.
  No check of symmetry is performed.

  @attention the Eigenvectors are returned row by row into the evecs
  tensor.

  By default the length of the eigenvectors is equivalent to the eigenvalues
  but it can be normalized to one.

  @section arguments Input Arguments

  Flag to normalize length of eigenvectors to 1.

  @return Eigenvalues and vectors are returned into the supplied Vector and
  TensorVariables, respectively. Method returns false if the rank of the matrix is zero.

  @section implementation Implementation

  Cubic root finding algorithm using Cardano's formula and the eigenvectors
  are calculated via LU-backsubstitution.

  @test SKM - observed, by comparison with Maple, that all the values are
  the same, but that the first element of the
  second Eigenvector always has the opposite sign (does this matter?)

  */
bool TensorVariable<3U>::Eigen( VectorVariable<3U>& evals, TensorVariable<3U>& evecs, bool bNormalize ) const
{
  //calculate eigenvalues
  //set eigenvalues
  if ( !EigenValues( evals ) ) return false;

  //suppose tensor is symmetric, by averaging the non diagonal elements
  const double f_0_1( (data[1][0] + data[0][1]) / 2. ); //=1_0
  const double f_0_2( (data[2][0] + data[0][2]) / 2. ); //=2_0
  const double f_1_2( (data[1][2] + data[2][1]) / 2. ); //=1_2

                                                          //EigenValueCalculation(Ev0_, Ev1_, Ev2_);
  TensorVariable<3U> Ev0_I( ANY, evals( 0 ) );
  TensorVariable<3U> Ev1_I( ANY, evals( 1 ) );
  TensorVariable<3U> Ev2_I( ANY, evals( 2 ) );

  TensorVariable<3U> A( ANY, data[0][0], f_0_1, f_0_2, f_0_1, data[1][1], f_1_2, f_0_2, f_1_2, data[2][2] );
  //A.Out();

  TensorVariable<3U> A_Ev0_I = A - Ev0_I;
  TensorVariable<3U> A_Ev1_I = A - Ev1_I;
  TensorVariable<3U> A_Ev2_I = A - Ev2_I;

  TensorVariable<3U> AEv0byAEv1 = A_Ev0_I * A_Ev1_I;
  TensorVariable<3U> AEv0byAEv2 = A_Ev0_I * A_Ev2_I;
  TensorVariable<3U> AEv1byAEv2 = A_Ev1_I * A_Ev2_I;

  VectorVariable<3U> LEvec( ANY, ANY, ANY,
                            sqrt( AEv1byAEv2( 0, 0 )*AEv1byAEv2( 0, 0 ) + AEv1byAEv2( 1, 0 )*AEv1byAEv2( 1, 0 ) + AEv1byAEv2( 2, 0 )*AEv1byAEv2( 2, 0 ) ),
                            sqrt( AEv0byAEv2( 0, 1 )*AEv0byAEv2( 0, 1 ) + AEv0byAEv2( 1, 1 )*AEv0byAEv2( 1, 1 ) + AEv0byAEv2( 2, 1 )*AEv0byAEv2( 2, 1 ) ),
                            sqrt( AEv0byAEv1( 0, 2 )*AEv0byAEv1( 0, 2 ) + AEv0byAEv1( 1, 2 )*AEv0byAEv1( 1, 2 ) + AEv0byAEv1( 2, 2 )*AEv0byAEv1( 2, 2 ) ) );

  uint32_t EvecLengthColID[3] = { 0, 1, 2 };

  for ( uint32_t i( 0 ); i < 3; i++ ) {
    uint32_t j = i;
    if ( LEvec( i ) < numeric_limits<double>::epsilon() ) {
      j++;
      if ( j > 2 ) { j = 0; }
      LEvec( i ) = sqrt( AEv1byAEv2( 0, j )*AEv1byAEv2( 0, j ) + AEv1byAEv2( 1, j )*AEv1byAEv2( 1, j ) + AEv1byAEv2( 2, j )*AEv1byAEv2( 2, j ) );
      if ( LEvec( i ) < numeric_limits<double>::epsilon() ) {
        j++;
        if ( j > 2 ) { j = 0; }
        LEvec( i ) = sqrt( AEv1byAEv2( 0, j )*AEv1byAEv2( 0, j ) + AEv1byAEv2( 1, j )*AEv1byAEv2( 1, j ) + AEv1byAEv2( 2, j )*AEv1byAEv2( 2, j ) );
      }
      EvecLengthColID[i] = j;
    }
  }

  // Normalised eigen vectors. Each column of the following matrix is the eigenvector for the corresponding eigenvalue.
  for ( uint32_t i( 0 ); i < 3; i++ ) {
    evecs( i, 0 ) = AEv1byAEv2( i, EvecLengthColID[0] ) / LEvec( EvecLengthColID[0] );
    evecs( i, 1 ) = AEv0byAEv2( i, EvecLengthColID[1] ) / LEvec( EvecLengthColID[1] );
    evecs( i, 2 ) = AEv0byAEv1( i, EvecLengthColID[2] ) / LEvec( EvecLengthColID[2] );
  }

  // Not Normalised eignvectors; i.e. vectors of principal stresses/strains.
  if ( !bNormalize ) {
    for ( uint32_t j( 0 ); j< 3; j++ ) {
      for ( uint32_t i( 0 ); i < 3; i++ ) {
        evecs( i, j ) *= evals( j );
      }
    }
  }

  // evecs inherits the flags from the original tensor variable
  evecs.flag[0] = flag[0];
  evecs.flag[1] = flag[1];
  evecs.flag[2] = flag[2];

  return true;

} // end Eigen





  /**

  Computes eigenvalues and eigenvectors for 3d TensorVariable (doesn't have to be symmetric)

  @section arguments Input Arguments

  eigenVectors, eigenValues

  @return Eigenvalues and vectors are returned into the supplied Vector and
  TensorVariables, respectively.

  @section implementation Implementation
  Uses:
  1) Symmetric Householder reduction to tridiagonal form,
  derived from the Algol procedures tred2 by
  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  Fortran subroutine in EISPACK.
  2) Symmetric tridiagonal QL algorithm,
  derived from the Algol procedures tql2, by
  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  Fortran subroutine in EISPACK.

  @author Lukas Mosser?

  */
bool TensorVariable<3U>::EigenNonSymmetric( VectorVariable<3U>& eigenVals,
                                            TensorVariable<3U>& eigenVecs ) const
{
  const int n = 3;
  double V[n][n], d[n], e[n];

  for ( int i = 0; i<n; i++ )
    for ( int j = 0; j<n; j++ )
      V[i][j] = (*this)(i, j);

  // Symmetric Householder reduction to tridiagonal form.
  //  This is derived from the Algol procedures tred2 by
  //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  //  Fortran subroutine in EISPACK.
  for ( int j = 0; j < n; j++ ) {
    d[j] = V[n - 1][j];
  }

  // Householder reduction to tridiagonal form.

  for ( int i = n - 1; i > 0; i-- ) {

    // Scale to avoid under/overflow.

    double scale = 0.0;
    double h = 0.0;
    for ( int k = 0; k < i; k++ ) {
      scale = scale + fabs( d[k] );
    }
    if ( scale == 0.0 ) {
      e[i] = d[i - 1];
      for ( int j = 0; j < i; j++ ) {
        d[j] = V[i - 1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    }
    else {

      // Generate Householder vector.

      for ( int k = 0; k < i; k++ ) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      double f = d[i - 1];
      double g = sqrt( h );
      if ( f > 0 ) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i - 1] = f - g;
      for ( int j = 0; j < i; j++ ) {
        e[j] = 0.0;
      }

      // Apply similarity transformation to remaining columns.

      for ( int j = 0; j < i; j++ ) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for ( int k = j + 1; k <= i - 1; k++ ) {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for ( int j = 0; j < i; j++ ) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      double hh = f / (h + h);
      for ( int j = 0; j < i; j++ ) {
        e[j] -= hh * d[j];
      }
      for ( int j = 0; j < i; j++ ) {
        f = d[j];
        g = e[j];
        for ( int k = j; k <= i - 1; k++ ) {
          V[k][j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[i - 1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.

  for ( int i = 0; i < n - 1; i++ ) {
    V[n - 1][i] = V[i][i];
    V[i][i] = 1.0;
    double h = d[i + 1];
    if ( h != 0.0 ) {
      for ( int k = 0; k <= i; k++ ) {
        d[k] = V[k][i + 1] / h;
      }
      for ( int j = 0; j <= i; j++ ) {
        double g = 0.0;
        for ( int k = 0; k <= i; k++ ) {
          g += V[k][i + 1] * V[k][j];
        }
        for ( int k = 0; k <= i; k++ ) {
          V[k][j] -= g * d[k];
        }
      }
    }
    for ( int k = 0; k <= i; k++ ) {
      V[k][i + 1] = 0.0;
    }
  }
  for ( int j = 0; j < n; j++ ) {
    d[j] = V[n - 1][j];
    V[n - 1][j] = 0.0;
  }
  V[n - 1][n - 1] = 1.0;
  e[0] = 0.0;

  // Symmetric tridiagonal QL algorithm.
  //  This is derived from the Algol procedures tql2, by
  //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  //  Fortran subroutine in EISPACK.

  for ( int i = 1; i < n; i++ ) {
    e[i - 1] = e[i];
  }
  e[n - 1] = 0.0;

  double f = 0.0;
  double tst1 = 0.0;
  double eps = pow( 2.0, -52.0 );
  for ( int l = 0; l < n; l++ ) {

    // Find small subdiagonal element

    tst1 = std::max( tst1, fabs( d[l] ) + fabs( e[l] ) );
    int m = l;
    while ( m < n ) {
      if ( fabs( e[m] ) <= eps*tst1 ) {
        break;
      }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if ( m > l ) {
      int iter = 0;
      do {
        iter = iter + 1;  // (Could check iteration count here.)

                          // Compute implicit shift

        double g = d[l];
        double p = (d[l + 1] - g) / (2.0 * e[l]);
        double r = sqrt( p*p + 1.0 );
        if ( p < 0 ) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l + 1] = e[l] * (p + r);
        double dl1 = d[l + 1];
        double h = g - d[l];
        for ( int i = l + 2; i < n; i++ ) {
          d[i] -= h;
        }
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        double c = 1.0;
        double c2 = c;
        double c3 = c;
        double el1 = e[l + 1];
        double s = 0.0;
        double s2 = 0.0;
        for ( int i = m - 1; i >= l; i-- ) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = sqrt( p*p + e[i] * e[i] );
          e[i + 1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for ( int k = 0; k < n; k++ ) {
            h = V[k][i + 1];
            V[k][i + 1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while ( fabs( e[l] ) > eps*tst1 );
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }

  // Sort eigenvalues and corresponding vectors.

  for ( int i = 0; i < n - 1; i++ ) {
    int k = i;
    double p = d[i];
    for ( int j = i + 1; j < n; j++ ) {
      if ( d[j] < p ) {
        k = j;
        p = d[j];
      }
    }
    if ( k != i ) {
      d[k] = d[i];
      d[i] = p;
      for ( int j = 0; j < n; j++ ) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }

  for ( int i = 0; i<n; i++ )
  {
    for ( int j = 0; j<n; j++ )
      eigenVecs( i, j ) = V[i][j];
    eigenVals( i ) = d[i];
  }

  return true;
}


/**

Code revision by Hossein Ageshlui (Nov. 2015), fixing several errors:

Previous code included two functions for eigenvalue calculation, using basically the same code,
one outputing eigenvalues as a vector and the other one as an object of the class VectorVariable<3U>&.

To shorten the code and to not to cause any issues with previously written code,
this new function has been added which calculates eigenvalues as scalar variables.

The old functions now only call this new function and output the eigenvalues as a vector or as an object of a VectorVariable.
The new code also sorts the eigenvalues.

Also all the pow (1./3.) calls are replaced with cbrt(). The pow (1./.3) does not produce a result when used with negative values and it is slow.
The source for solving the cubic equation is: http://mathworld.wolfram.com/CubicFormula.html

*/
bool TensorVariable<3U>::EigenValuesPositiveDefiniteSymmetricMatrix( double& eigenValue0, double& eigenValue1, double& eigenValue2 ) const
{
  const double PI( 3.14159265358979323846 );

  //HA: The tensor is expected to be symmetric. However, The symmetry of the matrix is ensured by averaging non diagonal elements.
  const double f_0_1( (data[1][0] + data[0][1]) / 2. ); //=1_0
  const double f_0_2( (data[2][0] + data[0][2]) / 2. ); //=2_0
  const double f_1_2( (data[1][2] + data[2][1]) / 2. ); //=1_2

                                                          // a0, a1 and a2 are coefficients of cubic equation " x^3 + a2 * x^2 + a1 * x + a0 = 0". The roots of this equation are 
                                                          // the Eigen values of the tensor
  const double a2 = -(data[0][0] + data[1][1] + data[2][2]);
  const double a1 = -(f_0_1 * f_0_1 + f_0_2 * f_0_2 +
                         f_1_2 * f_1_2 - data[0][0] * data[1][1] -
                         data[0][0] * data[2][2] - data[1][1] * data[2][2]);
  const double a0 = -(data[0][0] * data[1][1] * data[2][2] -
                         data[0][0] * f_1_2 * f_1_2 -
                         f_0_1 * f_0_1 * data[2][2] +
                         f_0_1 * f_0_2 * f_1_2 +
                         f_0_2  * f_0_1 * f_1_2 -
                         f_0_2 * f_0_2 * data[1][1]);
  const double R = (9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54.0;
  const double Q = (3 * a1 - a2 * a2) / 9.0;
  const double D = Q * Q * Q + R * R;

  if ( D > 0. ) // One real root and pair of complex conjugate roots
  {
    /*
    HA. This check will never be used since the stress tensor is symmetric. Hence:
    (1) A has exactly n(not necessarily distinct) eigenvalues.
    (2) Sets of n eigenvectors exist for each of eigenvalues, and they are mututally orthogonal.
    */
    double m = R + std::sqrt( D );
    double n = R - std::sqrt( D );
    m = (m > 0) ? std::cbrt( m ) : -1 * std::cbrt( -1 * m );
    n = (n > 0) ? std::cbrt( n ) : -1 * std::cbrt( -1 * n );
    eigenValue0 = m + n - a2 / 3.0;

    //csmp_error.notice( WARNING, "TensorVariable::EigenValues(double,double,double):",
    //                  "found complex conjugate roots when calculating the eigenvalues of a tensor.");
    return false;
  }

  if ( fabs( D ) <= numeric_limits<double>::epsilon() ) // The equation has three real roots, at least two of them are equal
  {

    eigenValue0 = 2 * std::cbrt( R ) - a2 / 3.0;
    eigenValue1 = eigenValue2 = -1 * std::cbrt( R ) - a2 / 3.0;

    if ( eigenValue0 < eigenValue1 ) // to return the roots in descending order
    {
      double temp = eigenValue0;
      eigenValue0 = eigenValue1;
      eigenValue2 = temp;
    }
  }

  if ( D < 0. ) // The equation has three real roots
  {
    double teta = std::acos( R / (std::abs( Q ) * std::sqrt( -1 * Q )) );
    eigenValue0 = 2 * std::sqrt( -1 * Q ) * std::cos( teta / 3.0 ) - a2 / 3.0;
    eigenValue1 = 2 * std::sqrt( -1 * Q ) * std::cos( teta / 3.0 + 120.0 / 180.0 * PI ) - a2 / 3.0;
    eigenValue2 = 2 * std::sqrt( -1 * Q ) * std::cos( teta / 3.0 + 240.0 / 180.0 * PI ) - a2 / 3.0;

    // returns the roots in descending order
    if ( eigenValue2 > eigenValue1 ) {
      double temp = eigenValue1;
      eigenValue1 = eigenValue2;
      eigenValue2 = temp;
    }
    else if ( eigenValue1 > eigenValue0 ) {
      double temp = eigenValue0;
      eigenValue0 = eigenValue1;
      eigenValue1 = temp;
    }
  }
  return true;
}



template ostream&  operator<<( ostream& stream, const TensorVariable<1U>& o );
template ostream&  operator<<( ostream& stream, const TensorVariable<2U>& o );
template ostream&  operator<<( ostream& stream, const TensorVariable<3U>& o );

} // end namespace csmp

