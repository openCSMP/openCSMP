#include "PropertyData.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "binaryReadWrite.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

// helper template functions

// classic way
size_t flagOffset( VARIABLE_TYPE type, size_t spatial_dimension, size_t array_length ) {
    if ( type == SCALAR or type == ARRAY ) return 1U;
    if ( type == VECTOR or type == TENSOR ) return spatial_dimension;
    // flagged array
    return array_length;
 }

size_t valueOffset( VARIABLE_TYPE type, size_t spatial_dimension, size_t array_length ) {
    if ( type == SCALAR ) return 1U;
    if ( type == VECTOR ) return spatial_dimension;
    if ( type == TENSOR ) return spatial_dimension * spatial_dimension;
    // array and flagged array
    return array_length;
 }

// for flags
/// for all variable types
template<VARIABLE_TYPE,size_t=1U> size_t flagStride( size_t array_length );

// full specialisations
template<> size_t flagStride<SCALAR>( size_t ) { return 1U; }

template<> size_t flagStride<ARRAY>( size_t ) { return 1U; }
template<> size_t flagStride<FLAGGEDARRAY>( size_t array_length ) { return array_length; }

template<> size_t flagStride<VECTOR>( size_t ) { return 1U; }
template<> size_t flagStride<VECTOR,2U>( size_t ) { return 2U; }
template<> size_t flagStride<VECTOR,3U>( size_t ) { return 3U; }

template<> size_t flagStride<TENSOR>( size_t ) { return 1U; }
template<> size_t flagStride<TENSOR,2U>( size_t ) { return 2U; }
template<> size_t flagStride<TENSOR,3U>( size_t ) { return 3U; }

// for data values
template<VARIABLE_TYPE,size_t=1U> size_t dataStride( size_t array_length );

// full specialisations
template<> size_t dataStride<SCALAR>( size_t ) { return 1U; }
template<> size_t dataStride<ARRAY>( size_t array_length ) { return array_length; }
template<> size_t dataStride<FLAGGEDARRAY>( size_t array_length ) { return array_length; }

template<> size_t dataStride<VECTOR>( size_t ) { return 1U; }
template<> size_t dataStride<VECTOR,2U>( size_t ) { return 2U; }
template<> size_t dataStride<VECTOR,3U>( size_t) { return 3U; }

template<> size_t dataStride<TENSOR>( size_t) { return 1U; }
template<> size_t dataStride<TENSOR,2U>( size_t ) { return 4U; }
template<> size_t dataStride<TENSOR,3U>( size_t ) { return 9U; }


/**
*/
PropertyData::PropertyData( PLACEMENT place, VARIABLE_TYPE type, size_t spatial_dimension, size_t array_length )
 : place_(place),
   type_(type),
   dim_(spatial_dimension),
   flag_stride_( flagOffset( type, spatial_dimension, array_length ) ),
   data_stride_( valueOffset( type, spatial_dimension, array_length ) )
 {
 } 


PropertyData::PropertyData( const PropertyData& pd )
 : place_(pd.place_),
   type_(pd.type_),
   dim_(pd.dim_),
   flags_(pd.flags_),
   data_(pd.data_),
   flag_stride_(pd.flag_stride_),
   data_stride_(pd.data_stride_)
 {
//    cerr <<"\nPropertyData: called copy constructor.\n";
 }


/// moves all the contents from pd over to this object
PropertyData::PropertyData( PropertyData&& pd )
 : place_{pd.place_},
   type_{pd.type_},
   dim_{pd.dim_},
   flags_{pd.flags_},
   data_{pd.data_},
   flag_stride_{pd.flag_stride_},
   data_stride_{pd.data_stride_}
 {
//    cerr <<"\nPropertyData: called move constructor.\n";
 }


PropertyData& PropertyData::operator=( const PropertyData& pd )
 {
    if ( &pd != this ) {
        assert( place_ == pd.place_ );
        assert( type_ == pd.type_ );
        assert( dim_ == pd.dim_ );
        flags_  = pd.flags_;
        data_   = pd.data_;
     }
   return *this;
 } 

// no dynamically allocated memory in this class
PropertyData::~PropertyData()
 {
 } 


// ACCESSORS AND MUTATORS

void PropertyData::PushBack( VARIABLE_FLAG flag ) { flags_.push_back( flag ); }

void PropertyData::PushBack( double64 val ) { data_.push_back( val ); }


/// scalars
VARIABLE_FLAG& PropertyData::Flag( size_t nth_value ) {
     assert( Type() == SCALAR );
     assert( !flags_.empty() );
     assert( nth_value < flags_.size() );
     return flags_[nth_value];
  }
  

/// scalars
VARIABLE_FLAG  PropertyData::Flag( size_t nth_value ) const {
     assert( Type() == SCALAR );
     assert( !flags_.empty() );
     assert( nth_value < flags_.size() );
     return flags_[nth_value];
  }
  
  
/// set flags: vectors and array variables
VARIABLE_FLAG& PropertyData::Flag( size_t nth_value, size_t ith_dim ) {
    assert( Type() == VECTOR or Type() == ARRAY or Type() == FLAGGEDARRAY );
    assert( !flags_.empty() );
    assert( nth_value < flags_.size() );
    assert( (Type() == VECTOR && ith_dim <= dim_) or (Type() == ARRAY && ith_dim == 0) or (Type() == FLAGGEDARRAY && ith_dim < flag_stride_) );
    return flags_[ nth_value * flag_stride_ + ith_dim ];
  }

/// read flags: vectors and array variables
VARIABLE_FLAG PropertyData::Flag( size_t nth_value, size_t ith_dim ) const {
    assert( Type() == VECTOR or Type() == ARRAY or Type() == FLAGGEDARRAY );
     assert( !flags_.empty() );
    assert( nth_value < flags_.size() );
    assert( (Type() == VECTOR && ith_dim <= dim_) or (Type() == ARRAY && ith_dim == 0) or (Type() == FLAGGEDARRAY && ith_dim < flag_stride_) );
    return flags_[ nth_value * flag_stride_ + ith_dim ];
 }
  
  
/// tensors
VARIABLE_FLAG& PropertyData::Flag( size_t nth_value, size_t ith_row, size_t jth_col ) {
     assert( Type() == TENSOR );
     assert( !flags_.empty() );
     assert( nth_value < flags_.size() );
     assert( ith_row <= dim_ );
     assert( jth_col <= dim_ );
     assert( jth_col == ith_row );
// since only 3 values are stored and i==j on the diagonal, only j is used
//     return flags_[ nth_value * flag_stride_ + ith_row * dim_ + jth_col ];
     return flags_[ nth_value * flag_stride_ + jth_col ];
 }

/// tensors (only their diagonal elements have flags)
VARIABLE_FLAG PropertyData::Flag( size_t nth_value, size_t ith_row, size_t jth_col ) const {
     assert( Type() == TENSOR );
     assert( !flags_.empty() );
     assert( nth_value < flags_.size() );
     assert( ith_row <= dim_ );
     assert( jth_col <= dim_ );
     assert( jth_col == ith_row );
// since only 3 values are stored and i==j on the diagonal, only j is used
//     return flags_[ nth_value * flag_stride_ + ith_row * dim_ + jth_col ];
     return flags_[ nth_value * flag_stride_ + jth_col ];
 }



/// scalars
double64& PropertyData::Value( size_t nth_value ) {
    assert( Type() == SCALAR );
    assert( !data_.empty() );
    assert( nth_value < data_.size() );
    return data_[nth_value];
 }
 
 
double64  PropertyData::Value( size_t nth_value ) const {
    assert( Type() == SCALAR );
    assert( !data_.empty() );
    assert( nth_value < data_.size() );
    return data_[nth_value];
 }
 
 
/// assign values: vectors and array variables
double64& PropertyData::Value( size_t nth_value, size_t ith_dim ) {
    assert( Type() == VECTOR or Type() == ARRAY or Type() == FLAGGEDARRAY );
    assert( !data_.empty() );
    assert( nth_value < data_.size() );
    assert( (Type() == VECTOR && ith_dim <= dim_) or (Type() == ARRAY && ith_dim <= data_stride_) or (Type() == FLAGGEDARRAY && ith_dim <= data_stride_) );
    return data_[ nth_value * data_stride_ + ith_dim ];
 }
 
/// retrieve values: vectors and array variables
double64  PropertyData::Value( size_t nth_value, size_t ith_dim ) const {
    assert( Type() == VECTOR or Type() == ARRAY or Type() == FLAGGEDARRAY );
    assert( !data_.empty() );
    assert( nth_value < data_.size() );
    assert( (Type() == VECTOR && ith_dim <= dim_) or (Type() == ARRAY && ith_dim <= data_stride_) or (Type() == FLAGGEDARRAY && ith_dim <= data_stride_) );
    return data_[ nth_value * data_stride_ + ith_dim ];
 }

 
/// tensors
double64& PropertyData::Value( size_t nth_value, size_t ith_row, size_t jth_col ) {
    assert( Type() == TENSOR );
    assert( !data_.empty() );
    assert( nth_value < data_.size() );
    assert( ith_row <= dim_ );
    assert( jth_col <= dim_ );
    return data_[ nth_value * data_stride_ + ith_row * dim_ + jth_col ];
 }
 
 
double64  PropertyData::Value( size_t nth_value, size_t ith_row, size_t jth_col ) const {
    assert( Type() == TENSOR );
    assert( !data_.empty() );
    assert( nth_value < data_.size() );
    assert( ith_row <= dim_ );
    assert( jth_col <= dim_ );
    return data_[ nth_value * data_stride_ + ith_row * dim_ + jth_col ];
 }



/**
    you may only reset using a key of the same  (i.e., ELEMENT )
*/
void PropertyData::Clear()
 {
    flags_.clear();
    data_.clear();
 } 


/**
     reserves storage for n-objects as opposed to determining the 
     size of the sheer storage which is what the Reserve() method
     with 2 arguments does.
*/
void PropertyData::Reserve( size_t n_objects )
 {
    flags_.reserve( n_objects * flag_stride_ );
    data_.reserve( n_objects * data_stride_ );
 }


/** 
    if extra values are created they will be initialised to ANY and not-a-number (NaN)
*/
void PropertyData::Reserve( size_t flag_capacity, size_t value_capacity )
 {
    flags_.reserve( flag_capacity );
    data_.reserve( value_capacity );
 }



/**
     reserves storage for n-objects as opposed to determining the 
     size of the sheer storage which is what the Reserve() method
     with 2 arguments does.
*/
void PropertyData::Resize( size_t n_objects )
 {
    flags_.resize( n_objects * flag_stride_, ANY );
    data_.resize( n_objects * data_stride_, std::numeric_limits<double64>::quiet_NaN() );
 }


/**
    if extra values are created they will be initialised to ANY and not-a-number (NaN)
*/
void PropertyData::Resize( size_t n_flags, size_t n_values )
 {
    flags_.resize( n_flags, ANY );
    data_.resize( n_values, std::numeric_limits<double64>::quiet_NaN() );
 }


/**
     keeps only those values with the indices stored in the map

     indices are assumed to range between 0..n-1
*/
void PropertyData::ReduceTo( const map<size_t,size_t>& o_n_elmt_ids )
 {
    if ( o_n_elmt_ids.empty() ) { 
         cerr <<"\nPropertyData::ReduceTo: Supplied map is empty. "<< std::endl;
         cerr <<" Therefore the PropertyData will be emptied."<< std::endl;
         Clear();
         return;
      }
   
    // complex operation that must proceed from the first to the last element
    // because the variables were written into one-dimensional array
   
    // values
    vector<VARIABLE_FLAG>  new_flags( o_n_elmt_ids.size() * flag_stride_ );
    vector<double64>       new_data( o_n_elmt_ids.size() * data_stride_ );
   
    if ( type_ == SCALAR or type_ == ARRAY ) {
         for ( map<size_t,size_t>::const_iterator
               it=o_n_elmt_ids.begin(); it!=o_n_elmt_ids.end(); it++ ) {
              new_flags[ (*it).second ] = flags_[ (*it).first ];
              new_data[ (*it).second ]  = data_[ (*it).first ];
           }
      }
    else if ( type_ == VECTOR or type_ == FLAGGEDARRAY ) {
             for ( map<size_t,size_t>::const_iterator
                  it=o_n_elmt_ids.begin(); it!=o_n_elmt_ids.end(); it++ ) {
                // flag and data stride are the same
                assert( flag_stride_ == data_stride_ );
                for ( size_t i=0U; i<data_stride_; ++i ) {
                     new_flags[ (*it).second + i ] = flags_[ (*it).first + i ];
                     new_data[ (*it).second + i ] = data_[ (*it).first + i ];
                  }
            }
       }
    else if ( type_ == TENSOR ) {
             for ( map<size_t,size_t>::const_iterator
                  it=o_n_elmt_ids.begin(); it!=o_n_elmt_ids.end(); it++ ) {
                // the flags have to be handled differently than the values
                for ( size_t i=0U; i<flag_stride_; ++i ) {
                     new_flags[ (*it).second + i ] = flags_[ (*it).first + i ];
                  }
                for ( size_t i=0U; i<data_stride_; ++i ) {
                     new_data[ (*it).second + i ] = data_[ (*it).first + i ];
                  }
            }
       }
   
    // resetting the containers
    flags_ = new_flags;
    data_  = new_data;
    // trimming excess storage
    vector<VARIABLE_FLAG>( flags_ ).swap( flags_ );
    vector<double64>( data_ ).swap( data_ );

 } // end ReduceTo
 
 


void  PropertyData::MinMaxOf( double64& tmin, double64& tmax ) const
{
    tmin = (*min_element( data_.begin(), data_.end() ));
    tmax = (*max_element( data_.begin(), data_.end() ));
}


// O.K.
void PropertyData::ScaleRangeTo( double64 tmin, double64 tmax )
 {
    double64  old_min, old_max;
    MinMaxOf( old_min, old_max );

    // test whether we are already O.K.
    if ( fabs(fabs(tmax) - fabs(old_max)) < numeric_limits<double64>::epsilon() &&
         fabs(fabs(tmin) - fabs(old_min)) < numeric_limits<double64>::epsilon() )
      return;
  
    const double64  old_range = old_max - old_min;
    const double64  new_range = tmax - tmin;

    transform( data_.begin(), data_.end(), data_.begin(),
              [tmin,tmax,old_min,old_range,new_range](double64 val)
               { return tmin + ((val - old_min)/old_range) * new_range; } );

    // testing for correctness
    double64 new_min, new_max;
    MinMaxOf( new_min, new_max );
    if ( new_min != tmin || new_max != tmax ) {
        std::cout <<"\nPropertyData<T>::ScaleRangeTo: Scaling failed."<< endl;
         Out(std::cout);
         throw range_error("PropertyData<T>::ScaleRangeTo");
      }

 } // end ScaleRangeTo 



/**
    Adds offset to all data values in the container
 
    TODO: make sure that this is done only where it makes sense!
*/
void PropertyData::OffsetRangeBy( double64 offset )
 {
    transform( data_.begin(), data_.end(), data_.begin(), [offset](double64 val){ return val + offset; } );

 } // end OffsetRangeBy 



/**
    Applies any kind of function with the syntax  double f(double) to
    the double64 entries of the stored data.
    
    TODO: make sure that this is done only where it makes sense!
*/
void PropertyData::TransformValues( double64 (*f)(double64) )
 {
     transform( data_.begin(), data_.end(), data_.begin(), (*f) );
 }




/**
 @fn  void PropertyData::OutBinary( FILE* fp ) const

 @brief Out binary. 
 
    const PLACEMENT             place_;        ///< placement of variable
    const VARIABLE_TYPE         type_;         ///< any of scalar..flagged array
    const size_t                dim_;          ///< spatial dimension 
    const size_t                flag_stride_;  ///< variable to variable offset (e.g. dim in a vector var)
    const size_t                data_stride_;  ///< variable to variable offset (e.g. dim in a vector var)
    std::vector<VARIABLE_FLAG>  flags_;        ///< variable flags
    std::vector<double64>       data_;         ///< variable values

 */
bool PropertyData::OutBinary( FILE* fp ) const
 {
     // writing the variable placement
     int32  var_placement = static_cast<int32>(place_);
     fwrite( (void*) &var_placement, sizeof(int32), 1, fp );

     // writing the variable type
     int32  var_type = static_cast<int32>(type_);
     fwrite( (void*) &var_type, sizeof(int32), 1, fp );

     // writing the spatial dimension
     int32  var_dim = static_cast<int32>(dim_);
     fwrite( (void*) &var_dim, sizeof(int32), 1, fp );

     // writing the flag stride
     int32  var_flag_stride = static_cast<int32>(flag_stride_);
     fwrite( (void*) &var_flag_stride, sizeof(int32), 1, fp );

     // writing the data stride
     int32  var_data_stride = static_cast<int32>(data_stride_);
     fwrite( (void*) &var_data_stride, sizeof(int32), 1, fp );
   
     // writing the number of records followed by flag values
     bool return_value = skm_C_fwrite( fp, flags_ );
   
     // writing the data values
     return_value = skm_C_fwrite( fp, data_ );
   
     return return_value;
 }
 

/**
    Reads PropertyData record written by PropertyData:OutBinary() and returns an initialised data object.

    Not a class member, this method constructs right-sized propery data container and returns it
    as an rvalue using the move constructor.
*/
PropertyData inBinaryPropertyData( FILE* fp )
 {
     // reading the variable placement
     int32  var_placement(UNSPECIFIED);
     fread( (void*) &var_placement, sizeof(int32), 1, fp );
     //assert( place_ == static_cast<PLACEMENT>(var_placement) );

     // reading the variable type
     int32  var_type(UNSPECIFIED);
     fread( (void*) &var_type, sizeof(int32), 1, fp );
     //assert( type_ == static_cast<VARIABLE_TYPE>(var_type) );

     // reading the spatial dimension
     int32  var_dim(UNSPECIFIED);
     fread( (void*) &var_dim, sizeof(int32), 1, fp );
     //assert( dim_ == static_cast<size_t>(var_dim) );

     // reading the flag stride
     int32  var_flag_stride(UNSPECIFIED);
     fread( (void*) &var_flag_stride, sizeof(int32), 1, fp );
     //assert( flag_stride_ == static_cast<size_t>(var_flag_stride) );

     // reading the data stride
     int32  var_data_stride(UNSPECIFIED);
     fread( (void*) &var_data_stride, sizeof(int32), 1, fp );
     //assert( data_stride_ == static_cast<size_t>(var_data_stride) );
   
     // calculating the array length
     size_t array_length = ( var_type == ARRAY or var_type == FLAGGEDARRAY ) ? array_length = var_data_stride : 0U;
   
     PropertyData data( static_cast<PLACEMENT>(var_placement), static_cast<VARIABLE_TYPE>(var_type), var_dim, array_length );

     // reading the number of records followed by flag values
     std::vector<VARIABLE_FLAG>  flags;
     skm_C_fread( fp, flags );
   
     // reading the data values
     std::vector<double64> values;
     skm_C_fread( fp, values );
   
     // pushing the data into Property record
     data.Reserve( flags.size(), values.size() );
     for ( auto it=flags.begin(); it!=flags.end(); ++it ) data.PushBack( (*it) );
     for ( auto it=values.begin(); it!=values.end(); ++it ) data.PushBack( (*it) );
   
     return data;
 }




void PropertyData::Out(std::ostream& os) const
 {
    os <<"\nPropertyData::Out: storage for "<< data_.size()/data_stride_ <<" "<< string(parseType(type_)) <<" objects."<< endl;
    os <<"Data placement: "<< string(parsePlacement(place_)) << endl;
    os <<"Spatial dimension: "<< dim_ << endl;
    os <<"Flag stride:       "<< flag_stride_ << endl;
    os <<"Data stride:       "<< data_stride_ << endl;
    os <<"\nflag values:\n";
    for ( size_t i=0U; i<flags_.size(); i++ )
      os << parseStatus(flags_[i]) <<" ";
    os << endl;

    os <<"\nvariable component values:\n";
    for ( size_t i=0U; i<data_.size(); i++ )
      os << data_[i] <<" ";
    os << endl;

 } // end Out 



bool  PropertyData::operator==( const PropertyData& d ) const
 {
    if ( !(place_       == d.place_) ) return false;
    if ( !(type_        == d.type_) ) return false;
    if ( !(dim_         == d.dim_) ) return false;
    if ( !(flag_stride_ == d.flag_stride_) ) return false;
    if ( !(data_stride_ == d.data_stride_) ) return false;
    if ( !(flags_       == d.flags_) ) return false;
    if ( !(data_        == d.data_) ) return false;
    return true;
 }

// PUSH BACK

void pushBack( PropertyData& data, const ScalarVariable& sc ) {
    assert( data.Type() == SCALAR );
    data.PushBack( sc.Flag() );
    data.PushBack( sc() );
 }
  
  
void pushBack( PropertyData& data, const ArrayVariable& ar ) {
    assert( data.Type() == ARRAY );
    data.PushBack( ar.Flag() );
    for ( size_t i=0U; i<ar.Size(); ++i )
      data.PushBack( ar[i] );
 }
 
 
void pushBack( PropertyData& data, const FlaggedArrayVariable& fa ) {
    assert( data.Type() == FLAGGEDARRAY );
    for ( size_t i=0U; i<fa.Size(); ++i ) {
         data.PushBack( fa.Flag(i) );
         data.PushBack( fa[i] );
      }
 }


template<size_t dim>
void pushBack( PropertyData& data, const VectorVariable<dim>& vc ) {
    assert( data.Type() == VECTOR );
    assert( data.Dim() == dim );
    for ( size_t i=0U; i<dim; ++i ) {
         data.PushBack( vc.Flag(i) );
         data.PushBack( vc[i] );
      }
 }
 
template void pushBack( PropertyData& data, const VectorVariable<1U>& vc );
template void pushBack( PropertyData& data, const VectorVariable<2U>& vc );
template void pushBack( PropertyData& data, const VectorVariable<3U>& vc );
 
template<size_t dim>
void pushBack( PropertyData& data, const TensorVariable<dim>& ts ) {
    assert( data.Type() == TENSOR );
    assert( data.Dim() == dim );
    for ( size_t i=0U; i<dim; ++i ) {
         data.PushBack( ts.Flag(i) );
         for ( size_t j=0U; j<dim; ++j )
           data.PushBack( ts(i,j) );
      }
 }

template void pushBack( PropertyData& data, const TensorVariable<1U>& vc );
template void pushBack( PropertyData& data, const TensorVariable<2U>& vc );
template void pushBack( PropertyData& data, const TensorVariable<3U>& vc );

// READ AND STORE FUNCTIONS

template<>
void store( PropertyData& data, size_t position, const ScalarVariable& sc ) {
    assert( data.Type() == SCALAR );
    data.Flag( position )  = sc.Flag();
    data.Value( position ) = sc();
 }

template<>
void store( PropertyData& data, size_t position, const ArrayVariable& ary ) {
    assert( data.Type() == ARRAY );
    data.Flag( position, 0U ) = ary.Flag();
    const size_t array_size(ary.Size());
    for ( size_t i=0U; i<array_size; ++i )
      data.Value( position * array_size + i ) = ary[i];
 }

template<>
void store( PropertyData& data, size_t position, const FlaggedArrayVariable& ary ) {
    assert( data.Type() == FLAGGEDARRAY );
    const size_t array_size(ary.Size());
    for ( size_t i=0U; i<array_size; ++i ) {
         data.Flag( position * array_size + i  ) = ary.Flag(i);
         data.Value( position * array_size + i ) = ary[i];
      }
 }

template<size_t dim>
void store( PropertyData& data, size_t position, const VectorVariable<dim>& vc ) {
    assert( data.Type() == VECTOR );
    assert( data.Dim() == dim );
    for ( size_t i=0U; i<dim; ++i ) {
         data.Flag( position, i )  = vc.Flag(i);
         data.Value( position, i ) = vc[i];
      }
 }

template void store( PropertyData&, size_t, const VectorVariable<1U>& );
template void store( PropertyData&, size_t, const VectorVariable<2U>& );
template void store( PropertyData&, size_t, const VectorVariable<3U>& );


template<size_t dim>
void store( PropertyData& data, size_t position, const TensorVariable<dim>& ts ) {
    assert( data.Type() == TENSOR );
    assert( data.Dim() == dim );
    for ( size_t i=0U; i<dim; ++i ) {
         data.Flag( position, i, i ) = ts.Flag(i);
         for ( size_t j=0U; j<dim; ++j )
           data.Value( position, i, j ) = ts(i,j);
      }
 }

template void store( PropertyData&, size_t, const TensorVariable<1U>& );
template void store( PropertyData&, size_t, const TensorVariable<2U>& );
template void store( PropertyData&, size_t, const TensorVariable<3U>& );


// READING PROPERTY DATA

template<>
void read( const PropertyData& data, size_t position, ScalarVariable& sc ) {
    assert( data.Type() == SCALAR );
    sc.Flag() = data.Flag( position );
    sc()      = data.Value( position );
 }

template<>
void read( const PropertyData& data, size_t position, ArrayVariable& ary ) {
    assert( data.Type() == ARRAY );
    ary.Flag() = data.Flag( position, 0U );
    const size_t array_size(data.Components());
    ary.Resize(array_size);
    for ( size_t i=0U; i<array_size; ++i )
      ary(i) = data.Value( position, i );
 }

template<>
void read( const PropertyData& data, size_t position, FlaggedArrayVariable& ary ) {
    assert( data.Type() == FLAGGEDARRAY );
    const size_t array_size(data.Components());
    ary.Resize(array_size);
    for ( size_t i=0U; i<array_size; ++i ) {
         ary.Flag(i) = data.Flag( position, i );
         ary(i) = data.Value( position, i );
      }
 }

template<size_t dim>
void read( const PropertyData& data, size_t position, VectorVariable<dim>& vc ) {
    assert( data.Type() == VECTOR );
    assert( data.Dim() == dim );
    for ( size_t i=0U; i<dim; ++i ) {
         vc.Flag(i) = data.Flag( position, i );
         vc(i) = data.Value( position, i );
      }
 }

template void read( const PropertyData&, size_t, VectorVariable<1U>& );
template void read( const PropertyData&, size_t, VectorVariable<2U>& );
template void read( const PropertyData&, size_t, VectorVariable<3U>& );


template<size_t dim>
void read( const PropertyData& data, size_t position, TensorVariable<dim>& ts ) {
    assert( data.Type() == TENSOR );
    assert( data.Dim() == dim );
    for ( size_t i=0U; i<dim; ++i ) {
         ts.Flag(i) = data.Flag( position, i, i );
         for ( size_t j=0U; j<dim; ++j )
           ts(i,j) = data.Value( position, i, j );
      }
 }

template void read( const PropertyData&, size_t, TensorVariable<1U>& );
template void read( const PropertyData&, size_t, TensorVariable<2U>& );
template void read( const PropertyData&, size_t, TensorVariable<3U>& );


} // end namespace csmp
