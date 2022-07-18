#ifndef CSMP_PROPERTY_DATA_H
#define CSMP_PROPERTY_DATA_H

#include "CSMP_definitions.h"

namespace csmp {

class ScalarVariable;
class ArrayVariable;
class FlaggedArrayVariable;
template<uint32_t> class VectorVariable;
template<uint32_t> class TensorVariable;
class PropertyData;

/// prototypes for appending variable values to container
void pushBack( PropertyData&, const ScalarVariable& );
void pushBack( PropertyData&, const ArrayVariable& );
void pushBack( PropertyData&, const FlaggedArrayVariable& );
template<uint32_t dim> void pushBack( PropertyData&, const VectorVariable<dim>& );
template<uint32_t dim> void pushBack( PropertyData&, const TensorVariable<dim>& );

/// storing the PropertyData data:  store( data, i, makeScalar(ANY,0.) );
template<typename csmp_var_type> void store( PropertyData&, size_t position, const csmp_var_type& );
template<uint32_t dim> void store( PropertyData&, size_t position, const VectorVariable<dim>& );
template<uint32_t dim> void store( PropertyData&, size_t position, const TensorVariable<dim>& );

/// reading PropertyData data: VectorVariable<3U> a;  read( data, i, a );
template<typename csmp_var_type> void read( const PropertyData&, size_t position, csmp_var_type& );
template<uint32_t dim> void read( const PropertyData&, size_t position, VectorVariable<dim>& );
template<uint32_t dim> void read( const PropertyData&, size_t position, TensorVariable<dim>& );

/**

  @brief flexible-size container for CSMP data types, providing range-check and conversion functions

  @author S.K. Matthai
  @date 8/4/2016

  @section motivation Motivation
  
  Efficient data transfer object with an efficient binary interface to 
  store variables to disk.

  @section design Design Intent

  Replace FEM_Data in VSet to streamline variable exchange and storage
  to disk; make it possible to read variables selectively from file.

  @section examples Application Examples

  Create contained for 1000 arrays, initialize it and write it to disk, for example:

  @code
  PropertyData data( NODE, ARRAY, DIM3, 21 );
  data.Reserve(1000);
  ArrayVariable  array(1000);
  ...
  for( a : 1000 )  data.PushBack( array );
  data.OutBinary( FILE );  
  @endcode

  @attention use read() and store() from above to retrieve entire csmp variables,
  there are no member accessors Read() or Store()

*/
class PropertyData {
  public:
    /// constructor (but no initializer) for all possible csmp variable types; array length gives the number of elements in array or flagged array variable
    PropertyData( PLACEMENT, VARIABLE_TYPE, uint32_t dim, uint32_t array_length=0U );
  
    PropertyData( const PropertyData& );
    PropertyData( PropertyData&& );
    ~PropertyData();
    PropertyData& operator=( const PropertyData& );
  
    /// compare to another dataset of this sort
    bool operator==( const PropertyData& ) const;
  
    PLACEMENT     Placement() const { return place_; }
    VARIABLE_TYPE Type() const { return type_; }
    uint32_t      Dim() const { return dim_; }
    size_t        Size() const { return data_.size(); }
    /// number of value entries stored in the variable (=size of array)
    uint32_t      Components() const { return data_stride_; };
  
    /// changing container size
    void Clear();

    /// for a user-defined number of objects of the type the property data was constructed for
    void Reserve( size_t n_objects );
  
    /// sheer storage not considering types of variables
    void Reserve( size_t flag_capacity, size_t value_capacity );

    /// for a user-defined number of objects of the type the property data was constructed for
    void Resize( size_t n_objects );

    /// sheer storage size not considering types of variables
    void Resize( size_t n_flags, size_t n_values );
  
    /// change the size if this required by parent container, e.g., VSet
    void ReduceTo( const std::map<size_t,size_t>& o_n_elmt_ids );
  
    // accessors / mutators for values (operator[] cannot be overloaded for this)
    /// @attention use read() and store() from above to retrieve entire csmp variables
  
    /// 'raw' inserter for flags (user must ascertain that stride is correct); use pushBack() for objects
    void PushBack( VARIABLE_FLAG );

    /// 'raw' inserter for values (user must ascertain that stride is correct); use pushBack() for objects
    void PushBack( double );
  
    /// inserter from another PropertyData
    void PushBackFrom( const PropertyData& data, size_t nth_value );
  
    /// accessors / mutators for flags
    /// scalars
    VARIABLE_FLAG& Flag( size_t nth_value );
    VARIABLE_FLAG  Flag( size_t nth_value ) const;
    /// vectors and array variables
    VARIABLE_FLAG& Flag( size_t nth_value, uint32_t ith_dim );
    VARIABLE_FLAG  Flag( size_t nth_value, uint32_t ith_dim ) const;
    /// tensors
    VARIABLE_FLAG& Flag( size_t nth_value, uint32_t ith_row, uint32_t jth_col );
    VARIABLE_FLAG  Flag( size_t nth_value, uint32_t ith_row, uint32_t jth_col ) const;
 
    /// scalars
    double&      Value( size_t nth_value );
    double       Value( size_t nth_value ) const;
    /// vectors and array variables
    double&      Value( size_t nth_value, uint32_t ith_dim );
    double       Value( size_t nth_value, uint32_t ith_dim ) const;
    /// tensors
    double&      Value( size_t nth_value, uint32_t ith_row, uint32_t jth_col );
    double       Value( size_t nth_value, uint32_t ith_row, uint32_t jth_col ) const;
  
    /// checks range; all values included
    void MinMaxOf( double& tmin, double& tmax ) const;
  
    /// scales the current variable range to the new one; all values are scaled
    void ScaleRangeTo( double tmin, double tmax );

    /// add supplied value to all entries (look at conventions for vectors and tensors in the doc of these classes)
    void OffsetRangeBy( double offset_value );
  
    /// takes a pointer to a typical math function like double sqrt(double) as argument; applied this function to all values
    void TransformValues( double (*f)(double) );

    /// writing stored flag and data values to file
    bool OutBinary( std::fstream& fp ) const;

    void Out() const;

  private:
    PropertyData() = delete;

  private:
    const PLACEMENT             place_;        ///< placement of variable
    const VARIABLE_TYPE         type_;         ///< any of scalar..flagged array
    const uint32_t              dim_;          ///< spatial dimension
    const uint32_t              flag_stride_;  ///< variable to variable offset (e.g. dim in a vector var)
    const uint32_t              data_stride_;  ///< variable to variable offset (e.g. dim in a vector var)
    std::vector<VARIABLE_FLAG>  flags_;        ///< variable flags
    std::vector<double>         data_;         ///< variable values
 };


/// reading stored flag and data vaues from file
PropertyData inBinaryPropertyData( std::fstream& fp );

/// calculating the distance bwetween consecutive data entries in the container
uint32_t flagOffset( VARIABLE_TYPE, uint32_t spatial_dimension, uint32_t array_length );
uint32_t valueOffset( VARIABLE_TYPE, uint32_t spatial_dimension, uint32_t array_length );

} // csmp

#endif





