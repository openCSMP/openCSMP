#include "FEM_Data.h"
#include "TensorVariable.h"
#include "binaryReadWrite.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

using namespace std;

namespace csmp {

template<typename csp_type>
FEM_Data<csp_type>::FEM_Data()
 : place(NODE), logarithmitized(false)
 {
 } 


template<typename csp_type>
FEM_Data<csp_type>::FEM_Data( PLACEMENT p, size_t i )
 : place(p), logarithmitized(false)
 {
    data.reserve(i);
    for ( size_t j=0U; j<i; j++ ) data.push_back( csp_type() );
 } 


template<typename csp_type>
FEM_Data<csp_type>::FEM_Data( PLACEMENT p, const vector<csp_type>& d )
 : place(p), data(d), logarithmitized(false)
 {
 } 


template<typename csp_type>
FEM_Data<csp_type>& FEM_Data<csp_type>::operator=( const FEM_Data<csp_type>& p )
 {
   if ( &p == this ) return *this;
   place           = p.place;
   data            = p.data; 
   logarithmitized = p.logarithmitized;
   return *this;
 } 


template<typename csp_type>
FEM_Data<csp_type>::FEM_Data( const FEM_Data<csp_type>& d )
 : place(d.place), data(d.data), logarithmitized(d.logarithmitized)
 {
 } 


template<typename csp_type>
FEM_Data<csp_type>::~FEM_Data()
 {
 } 


/// you may only reset using a key of the same <csp_type> (i.e., ELEMENT )
template<typename csp_type>
void FEM_Data<csp_type>::Reset( const csmp::Index& setting, size_t size, const csp_type& val )
 {
    place           = setting.place;
    logarithmitized = false;
    data.erase( data.begin(), data.end() );
    if ( size > 0 ) {
         data.reserve(size);
         for( size_t i=0U; i<size; i++ ) data.push_back( val );
      } 
 } 


template<typename csp_type>
void FEM_Data<csp_type>::Extend( size_t new_size )
 {
    if ( new_size <= data.size() )
      {
         cout <<"\nFEM_Data<csp_type>::Extend: New size is smaller than present size ! ";
         cout <<"Nothing is done."<< endl;
         return;
      }
    data.reserve( new_size );
    for ( size_t i=data.size(); i<new_size; i++ ) data.push_back( csp_type() );  
 }


// indices are assumed to range between 0..n-1
template<typename csp_type>
void FEM_Data<csp_type>::ReduceTo( const map<size_t,size_t>& o_n_elmt_ids )
 {
    if ( o_n_elmt_ids.empty() ) { 
         cout <<"\nFEM_Data<csp_type>::ReduceTo: Supplied list is empty. "<< std::endl;
         cout <<" Therefore the FEM data will also be emptied."<< std::endl;
      }
      
    vector<csp_type>  new_data(o_n_elmt_ids.size());
    for ( map<size_t,size_t>::const_iterator
          it=o_n_elmt_ids.begin(); it!=o_n_elmt_ids.end(); it++ ) 
      new_data[ (*it).second ] = data[ (*it).first ];
    
    data = new_data;
 }



template<typename csp_type>
void FEM_Data<csp_type>::ScaleRangeTo( const csp_type& tmin, const csp_type& tmax )
 {
    csp_type  old_min, old_max;

    MinMaxOf( old_min, old_max );

    // test whether we are already O.K.
    if ( tmax == old_max && tmin == old_min ) return;
  
    csp_type  old_range = old_max - old_min;
    csp_type  new_range = tmax - tmin;  

    for ( typename vector<csp_type>::iterator
          it=data.begin(); it!=data.end(); it++ )
       *it = tmin + ((*it - old_min)/old_range) * new_range;

    // testing for correctness
    csp_type new_min, new_max;
    MinMaxOf( new_min, new_max );
    if ( new_min != tmin || new_max != tmax ) {
         cout <<"\nFEM_Data<T>::ScaleRangeTo: Scaling failed."<< endl;
         Out(cout);
         throw range_error("FEM_Data<T>::ScaleRangeTo");
      }

 } // end ScaleRangeTo 



template<typename csp_type>
void FEM_Data<csp_type>::OffsetRangeBy( csp_type df )
 {
    for ( typename vector<csp_type>::iterator it=data.begin(); it!=data.end(); it++ ) *it += df;

 } // end OffsetRangeBy 


/*
template<typename csp_type>
void FEM_Data<csp_type>::LogarithmOfValues()
{
  typename vector<csp_type>::iterator
  it = min_element( data.begin(), data.end() );

  for ( it=data.begin(); it!=data.end(); it++ )
    {
       (*it).Ln();
    }

 logarithmitized = true;
}



template<typename csp_type>
void FEM_Data<csp_type>::DecadicLogarithmOfValues()
{
  typename vector<csp_type>::iterator
  it = min_element( data.begin(), data.end() );
  for ( it=data.begin(); it!=data.end(); it++ )
    {
       (*it).Log10();
    }

 logarithmitized = true;
}


template<typename csp_type>
void FEM_Data<csp_type>::SquareRootOfValues()
{
  typename vector<csp_type>::iterator
  it = min_element( data.begin(), data.end() );
  for ( it=data.begin(); it!=data.end(); it++ ) 
    {
       (*it).Sqrt();
    }
}
*/

template<typename csp_type>
void  FEM_Data<csp_type>::MinMaxOf( csp_type& tmin, csp_type& tmax ) const
{
    tmin = (*min_element( data.begin(), data.end() ));
    tmax = (*max_element( data.begin(), data.end() ));
}


namespace femDataBinaryDispatch {

  template<typename csp_type>
  bool outBinary( FILE* fp, const vector<csp_type>& cntr )
    {
      return skm_C_fwrite( fp, cntr );
    }

  template<typename csp_type>
  bool inBinary( FILE* fp, vector<csp_type>& cntr )
    {
      return skm_C_fread( fp, cntr );
    }

  template<>
  bool outBinary( FILE* fp, const vector<ArrayVariable>& arrayVariables )
    {      
      const size_t dataDepth( arrayVariables.size() );

      // size
      fwrite( (void*) &dataDepth, sizeof(size_t), 1, fp );

      // data
      for ( size_t i(0); i < dataDepth; ++i )
        if( !arrayVariables[i].Out(fp) )
          return false;

      return true;
    }

  template<>
  bool inBinary( FILE* fp, vector<ArrayVariable>& arrayVariables )
    {   
      // size
      size_t dataDepth(0);
      fread( (void*) &dataDepth, sizeof(size_t), 1, fp );
      arrayVariables = vector<ArrayVariable>( dataDepth, ArrayVariable() );

      // data
      for ( size_t i(0); i < dataDepth; ++i )
        if( !arrayVariables[i].In(fp) )
          return false;
      arrayVariables.swap(arrayVariables);
      return true;
    }

  template<>
  bool outBinary( FILE* fp, const vector<FlaggedArrayVariable>& flaggedArrayVariables )
    {
      const size_t dataDepth( flaggedArrayVariables.size() );

      // size
      fwrite( (void*) &dataDepth, sizeof(size_t), 1, fp );

      // data
      for ( size_t i(0); i < dataDepth; ++i )
        if( !flaggedArrayVariables[i].Out(fp) )
          return false;

      return true;
    }

  template<>
  bool inBinary( FILE* fp, vector<FlaggedArrayVariable>& flaggedArrayVariables )
    {
      // size
      size_t dataDepth(0);
      fread( (void*) &dataDepth, sizeof(size_t), 1, fp );
      flaggedArrayVariables = vector<FlaggedArrayVariable>( dataDepth, FlaggedArrayVariable() );

      // data
      for ( size_t i(0); i < dataDepth; ++i )
        if( !flaggedArrayVariables[i].In(fp) )
          return false;
      flaggedArrayVariables.swap(flaggedArrayVariables);
      return true;
    }

  } // femDataBinaryDispatch


/**
 @fn  template<typename csp_type> void FEM_Data<csp_type>::OutBinary( FILE* fp ) const

 @brief Out binary. 

 @tparam  csp_type  Type of the csp type.
 @param [in,out]  fp  If non-null, the fp.
 */
template<typename csp_type>
bool FEM_Data<csp_type>::OutBinary( FILE* fp ) const
 {
     // writing the variable placement
     int32  var_placement = static_cast<int32>(place);
     fwrite( (void*) &var_placement, sizeof(int32), 1, fp );
     
     // writing the dataset
     return femDataBinaryDispatch::outBinary( fp, data );
 }
 
 
template<typename csp_type>
void FEM_Data<csp_type>::InBinary( FILE* fp )
 {
     // reading the variable placement
     int32  var_placement(UNSPECIFIED);
     fread( (void*) &var_placement, sizeof(int32), 1, fp );
     place = static_cast<PLACEMENT>(var_placement);
     
     // reading the dataset
     femDataBinaryDispatch::inBinary( fp, data );
 }


template<typename csp_type>
void FEM_Data<csp_type>::Out(std::ostream& os) const
 {
    os <<"\nFEM_Data::Out:"<< endl;
    os <<"Data placement: "<< string(parsePlacement(place)) << endl;
    if ( logarithmitized ) os <<"\nLogarithmitized data: "<< endl;
    else                   os <<"\nData ("<< data.size()<<" vals): "<< endl;
 
    for ( size_t i=0U; i<data.size(); i++ )
      os << data[i] <<"\t";
    os << endl;

 } // end Out 


template<typename csp_type>
bool  FEM_Data<csp_type>::operator==( const FEM_Data<csp_type>& d ) const
 {
    if ( !(place == d.place) ) return false;
    if ( !(logarithmitized == d.logarithmitized) ) return false;
    if ( !(data == d.data) ) return false;
    return true;
 }


#ifdef SKM_TEMPLATE_TESTING
  void
main()
 {
    FEM_Data<ScalarVariable > test_data( NODE, 100 );
    double64 min, max;
    
    test_data.Range( min, max );
    test_data.ScaleRangeTo( 0.001, 7.2 );
 }
 
#endif  

template class FEM_Data<ScalarVariable>;
template class FEM_Data<VectorVariable<1U> >;
template class FEM_Data<TensorVariable<1U> >;
template class FEM_Data<VectorVariable<2U> >;
template class FEM_Data<TensorVariable<2U> >;
template class FEM_Data<VectorVariable<3U> >;
template class FEM_Data<TensorVariable<3U> >;
template class FEM_Data<ArrayVariable>;
template class FEM_Data<FlaggedArrayVariable>;

} // end namespace csmp
