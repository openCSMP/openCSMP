#ifndef CSMP_FEM_DATA_H
#define CSMP_FEM_DATA_H

#include "CSMP_definitions.h"
#include "ArrayVariable.h"
#include "binaryReadWrite.h"

namespace csmp {

/**
    flexible-size container for CSMP data types, providing range-check and conversion functions
    
    @todo SKM: since this class is used to store models to disk, padded CSMP types rather than
    short int's and doubles are stored. This makes models on disk very big.
*/
template<typename csp_type>
class FEM_Data {
  public:
    FEM_Data();
    FEM_Data( const FEM_Data<csp_type>& );
    FEM_Data( PLACEMENT, size_t i );
    FEM_Data( PLACEMENT, const std::vector<csp_type>& );
    ~FEM_Data();
    FEM_Data<csp_type>& operator=( const FEM_Data<csp_type>& );
    bool operator==( const FEM_Data<csp_type>& ) const;
    PLACEMENT Placement() const { return place; };
    size_t Size() const { return data.size(); };
  
    /// changing container size
    void Reset( const csmp::Index&, size_t size, const csp_type& );
    void Extend( size_t new_size );
    void ReduceTo( const std::map<size_t,size_t>& o_n_elmt_ids );
  
    /// access
    csp_type& operator[]( size_t i );
    const csp_type& operator[]( size_t i ) const;
    csp_type  operator()( size_t i ) const;
  
    /// checking range
    void MinMaxOf( csp_type& tmin, csp_type& tmax ) const;
    /// rescaling
    void ScaleRangeTo( const csp_type& tmin, const csp_type& tmax );
    /// add supplied value to all entries (look at conventions for vectors and tensors in the doc of these classes)
    void OffsetRangeBy( csp_type df );
  
    /// checks whether variable has already been taken the log of
    bool Logarithmitized() const { return logarithmitized; };
    /// convert values to their natural logarithm
    void LogarithmOfValues();
    /// convert values to base 10 logarithm
    void DecadicLogarithmOfValues();
    void SquareRootOfValues();

    bool OutBinary( std::FILE* fp ) const;
    void InBinary( std::FILE* fp );
    void Out(std::ostream& os) const;

protected:
    PLACEMENT              place;
    std::vector<csp_type>  data;
    bool                   logarithmitized;
 };


template<typename csp_type>
inline csp_type& FEM_Data<csp_type>::operator[]( size_t i ) 
 {
#ifndef NDEBUG
    if ( i >= data.size() || data.empty() ) {
         std::cerr <<"\nFEM_Data<csp_type>[]: Range exceeded: ";
         std::cerr << i <<" instead of: 0 up to "<< (data.size() - 1) << std::endl;
         throw std::range_error("FEM_Data<csp_type>::operator[]");
         return data[0];
    }
#endif
    return data[i]; 
 }


template<typename csp_type>
inline const csp_type& FEM_Data<csp_type>::operator[]( size_t i ) const
 {
#ifndef NDEBUG
    if ( i >= data.size() || data.empty() ) {
         std::cerr <<"\nFEM_Data<csp_type>[]: Range exceeded: ";
         std::cerr << i <<" instead of: 0 up to "<< (data.size() - 1) << std::endl;
         throw std::range_error("FEM_Data<csp_type>::operator[]");
         return data[0];
    }
#endif
    return data[i]; 
 }



template<typename csp_type>
inline csp_type  FEM_Data<csp_type>::operator()( size_t i ) const 
 {
#ifndef NDEBUG
    if ( i >= data.size() || data.empty() ) {
         std::cerr <<"\nFEM_Data<csp_type>(): Range exceeded: ";
         std::cerr << i <<" instead of: 0 up to "<< (data.size() - 1) << std::endl;
         throw std::range_error("FEM_Data<csp_type>::operator()");
         return data[0];
    }
#endif
    return data[i]; 
}


} // csmp

#endif





