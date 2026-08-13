#include <type_traits>
#include "LocalVariableStorage.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Edge.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Model.h"

using namespace std;

namespace csmp {

// ---------------------------------------------------------------
//       CONSTRUCTORS
// ---------------------------------------------------------------

/// Ctor for storees without IntegrationPointVariables
template<uint32_t dim, template<uint32_t> class STOREE>
LocalVariableStorage<dim, STOREE>::LocalVariableStorage( const LocalVariables& lv )
    : data_()
{
    ResizePropertyStorage(lv);
}


/// Goto overload for storees with IntegrationPointVariables @attention Will fail AT COMPILE TIME (very nice...) if used for an storee without integration points
template<uint32_t dim, template<uint32_t> class STOREE>
LocalVariableStorage<dim, STOREE>::LocalVariableStorage( const LocalVariables& lv, const IntegrationPointVariables& ipv )
    : data_()
{
    ResizePropertyStorage(lv, ipv);
}



template<uint32_t dim, template<uint32_t> class STOREE>
LocalVariableStorage<dim, STOREE>::LocalVariableStorage( const LocalVariableStorage<dim,STOREE>& lvs )
 : data_{lvs.data_}
 {}
 

template<uint32_t dim, template<uint32_t> class STOREE>
LocalVariableStorage<dim, STOREE>:: LocalVariableStorage( LocalVariableStorage<dim,STOREE>&& lvs )
 : data_{lvs.data_}
 {}



/// Goto overload for storees without IntegrationPointVariables
template<uint32_t dim, template<uint32_t> class STOREE>
bool LocalVariableStorage<dim,STOREE>::operator==( const LocalVariableStorage<dim,STOREE>& lvs ) const
  {
     if ( this == &lvs ) return true;
     if ( this->data_.flags != lvs.data_.flags ) return false;
     if ( this->data_.data  != lvs.data_.data ) return false;
     return true;
  }




/// Goto overload for storees without IntegrationPointVariables
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::ResizePropertyStorage( const LocalVariables& lv )
  {
    ResizePropertyStorage( lv.totalDataDepth, lv.totalFlagDepth );

// initialise auxiliary parameters for debugging
#ifndef NDEBUG
    StoreLocalState(lv);
#endif
  }
  


/// Goto overload for storees with IntegrationPointVariables
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::ResizePropertyStorage( const LocalVariables& lv, const IntegrationPointVariables& ipv )
  {
    const pair<uint32_t, uint32_t> newContainerSize = localVariableDispatch::containerNewSize( static_cast<const STOREE<dim>*>(this), lv, ipv );
    ResizePropertyStorage(newContainerSize.first, newContainerSize.second);

// initialise auxiliary parameters for debugging
#ifndef NDEBUG
    StoreLocalState(lv);
#endif
  }





/// Resizes the storage preserving original values if any (this is where all ResizePropertyStorage end up)
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::ResizePropertyStorage( int_type newDataComponentCount, int_type newFlagComponentCount )
  {
    // needed?
    if( data_.flags.size() == newFlagComponentCount && data_.data.size() == newDataComponentCount )
      return;

    // resizing
    data_.flags.resize( newFlagComponentCount, ANY );
    data_.data.resize( newDataComponentCount, numeric_limits<double>::quiet_NaN() );

    // trimming excess capacity
    vector<VARIABLE_FLAG>( data_.flags ).swap( data_.flags );
    vector<double>( data_.data ).swap( data_.data );
  }





#ifndef NDEBUG
/// We check in debug only
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim, STOREE>::AssertPlacement(const csmp::Index& idx) const
{
   assert( (parsePlacement<dim,STOREE>()) == idx.place );
}


/// We check in debug only
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim, STOREE>::AssertIntegrationPointPlacement(const csmp::Index& idx) const
{
    assert(idx.place == ELEMENT_INTEGRATION_POINT || idx.place == FACE_INTEGRATION_POINT || idx.place == INTER_FACE_INTEGRATION_POINT);
}

/// We check in debug only
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim, STOREE>::AssertFiniteVolumeIntegrationPointPlacement(const csmp::Index& idx) const
{
    assert(idx.place == SECTOR_INTEGRATION_POINT || idx.place == FACET_INTEGRATION_POINT ||
        idx.place == FACE_SECTOR_INTEGRATION_POINT || idx.place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
        idx.place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT);
}

/// We store in debug only
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim, STOREE>::StoreLocalState(const LocalVariables& lv)
{
    data_.scalars = lv.scalars;
    data_.vectors = lv.vectors;
    data_.tensors = lv.tensors;
    data_.arrays = lv.arrayCount;
    data_.arrayLength = lv.arrayLength;
    data_.flaggedArrays = lv.flaggedArrayCount;
    data_.flaggedArrayLength = lv.flaggedArrayLength;
}
#endif



/// Creates an empty space in the storage preserving the original variable values. LocalVariables in Index have to reflect NEW state (incl added)
/**
@todo (2-F) Integration point properties not supported
*/
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::AddProperty( const csmp::Index& prop_key )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    /// @todo (2-F) Asserts missing

    // Resize to new state (this already supports ip vars)
    const pair<uint32_t,uint32_t> newTotaDataDepth = localVariableDispatch::containerTotalDataDepth( storeePtr, prop_key );
    ResizePropertyStorage( newTotaDataDepth.first, newTotaDataDepth.second );
    const uint32_t dataSize = static_cast<uint32_t>(data_.data.size());
    const uint32_t flagSize = static_cast<uint32_t>(data_.flags.size());

    /// Roman, 2013: with ipvs support

    // Offset
    const pair<uint32_t,uint32_t> newOffset = localVariableDispatch::containerOffset( storeePtr, prop_key );
    const uint32_t dataOffset         ( newOffset.first       );
    const uint32_t flagOffset         ( newOffset.second      );

    // DataDepth
    const uint32_t dataDepth          ( prop_key.dataDepth    );
    const uint32_t flagDepth          ( prop_key.flagDepth    );

    // Cycles
    const pair<int_type,int_type> newCycles = localVariableDispatch::containerIPCycles( storeePtr, prop_key );
    const int_type nipCycles1            ( newCycles.first       );
    const int_type nipCycles2            ( newCycles.second      );

    const pair<int_type,int_type> newCycles1Data = localVariableDispatch::containerIPCycle1Offset( storeePtr, prop_key );
    const int_type ipCycle1DataOffset    ( newCycles1Data.first  );
    const int_type ipCycle1FlagOffset    ( newCycles1Data.second );

    const pair<int_type,int_type> newCycles2Data = localVariableDispatch::containerIPCycle2Offset( storeePtr, prop_key );
    const int_type ipCycle2DataOffset    ( newCycles2Data.first  );
    const int_type ipCycle2FlagOffset    ( newCycles2Data.second );

    // Data and Flag Bounds
    vector<vector<pair<int_type,int_type> > > dataBounds( nipCycles1, vector<pair<int_type,int_type> >(nipCycles2+1,pair<int_type,int_type>(0,0)) );
    for ( int_type cycle1=0; cycle1<nipCycles1; ++cycle1  )
    {
        for ( int_type cycle2=0; cycle2<nipCycles2; ++cycle2  )
        {
            dataBounds[cycle1][cycle2].first  = dataOffset+cycle1*ipCycle1DataOffset + cycle2*ipCycle2DataOffset;
            dataBounds[cycle1][cycle2].second = dataBounds[cycle1][cycle2].first + dataDepth;
        }
        dataBounds[cycle1][nipCycles2].first  = dataOffset+(cycle1+1)*ipCycle1DataOffset;
        dataBounds[cycle1][nipCycles2].second = dataOffset+(cycle1+1)*ipCycle1DataOffset;
    }
    assert( nipCycles1 > 0
        && "LocalVariableStorage::AddProperty: nipCycles1 is zero. "
           "This means containerIPCycles returned a zero first cycle count "
           "for the property being added. "
           "For SECTOR_INTEGRATION_POINT and FACET_INTEGRATION_POINT "
           "placements, the element sector/facet counts must be non-zero "
           "at the time CreateProperty is called. "
           "Check that the Index was constructed with correct ipFactorSector "
           "and offsetFactorSector values, and that the mesh has been fully "
           "initialised before creating sector or facet integration point "
           "properties." );

    dataBounds[nipCycles1-1][nipCycles2].first  = dataSize;
    dataBounds[nipCycles1-1][nipCycles2].second = dataSize;

    vector<vector<pair<int_type,int_type> > > flagBounds( nipCycles1, vector<pair<int_type,int_type> >(nipCycles2+1,pair<int_type,int_type>(0,0)) );
    for ( int_type cycle1=0; cycle1<nipCycles1; ++cycle1  )
    {
        for ( int_type cycle2=0; cycle2<nipCycles2; ++cycle2  )
        {
            flagBounds[cycle1][cycle2].first  = flagOffset+cycle1*ipCycle1FlagOffset + cycle2*ipCycle2FlagOffset;
            flagBounds[cycle1][cycle2].second = flagBounds[cycle1][cycle2].first + flagDepth;
        }
        flagBounds[cycle1][nipCycles2].first  = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
        flagBounds[cycle1][nipCycles2].second = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
    }

    assert( nipCycles1 > 0
        && "LocalVariableStorage::AddProperty: nipCycles1 is zero "
           "for flagBounds. Same root cause as dataBounds assert above." );

    flagBounds[nipCycles1-1][nipCycles2].first  = flagSize;
    flagBounds[nipCycles1-1][nipCycles2].second = flagSize;

    // Moving data (same for all types) note: decrementing unsigned int avoiding comparison with zero
    for ( auto cycle1=nipCycles1; cycle1-- > 0; ){
        for ( auto cycle2=nipCycles2; cycle2-- > 0; ){
            const auto dataStart = dataBounds[ cycle1 ][ cycle2     ].second;
            const auto dataEnd   = dataBounds[ cycle1 ][ cycle2 + 1 ].second;
            const int_type stride = (cycle1 + 1) * (cycle2 + 1) * dataDepth;
            for ( int_type i=dataEnd; i-- > dataStart; )
              data_.data[i] = data_.data[i-stride];
        }
    }
    // Moving flags (same for all types)
    for ( auto cycle1 = nipCycles1; cycle1-- > 0; ){
        for ( auto cycle2 = nipCycles2; cycle2-- > 0; ){
            const auto flagsStart = flagBounds[ cycle1 ][ cycle2     ].second;
            const auto flagsEnd   = flagBounds[ cycle1 ][ cycle2 + 1 ].second;
            const int_type stride = (cycle1 + 1) * (cycle2 + 1) * flagDepth;
            for ( int_type i = flagsEnd; i-- > flagsStart; )
                data_.flags[i] = data_.flags[i - stride];
        }
    }
    
    // Adding non-initialized data of new variable
    for ( int_type cycle1=0; cycle1<nipCycles1; ++cycle1  )
        for ( int_type cycle2=0; cycle2<nipCycles2; ++cycle2  )
            for ( int_type i=dataBounds[cycle1][cycle2].first; i<dataBounds[cycle1][cycle2].second; ++i )
                data_.data[i] = numeric_limits<double>::quiet_NaN();
                
    // Adding non-initialized flags of new variable
    for ( int_type cycle1=0; cycle1<nipCycles1; ++cycle1  )
        for ( int_type cycle2=0; cycle2<nipCycles2; ++cycle2  )
            for ( int_type i=flagBounds[cycle1][cycle2].first; i<flagBounds[cycle1][cycle2].second; ++i )
                data_.flags[i] = ANY;


#ifndef NDEBUG
    StoreLocalState(prop_key.localVariables);
#endif
  } // end AddProperty





/// Deletes data and flags of property LocalVariables in Index have to reflect CURRENT state (incl property to be removed)
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::DeleteProperty( const csmp::Index& prop_key )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);

    /// @todo (2-F) Asserts missing
    const auto flagSize = data_.flags.size();
    const auto dataSize = data_.data.size();

        /// Roman, 2013: with ipvs support

        // Offset
        const pair<int_type,int_type> newOffset = localVariableDispatch::containerOffset( storeePtr, prop_key );
        const int_type dataOffset         ( newOffset.first       );
        const int_type flagOffset         ( newOffset.second      );

        // DataDepth
        const int_type dataDepth          ( prop_key.dataDepth    );
        const int_type flagDepth          ( prop_key.flagDepth    );

        // Cycles
        const pair<int_type,int_type> newCycles = localVariableDispatch::containerIPCycles( storeePtr, prop_key );
        const int_type nipCycles1            ( newCycles.first       );
        const int_type nipCycles2            ( newCycles.second      );

        const pair<int_type,int_type> newCycles1Data = localVariableDispatch::containerIPCycle1Offset( storeePtr, prop_key );
        const int_type ipCycle1DataOffset    ( newCycles1Data.first  );
        const int_type ipCycle1FlagOffset    ( newCycles1Data.second );

        const pair<int_type,int_type> newCycles2Data = localVariableDispatch::containerIPCycle2Offset( storeePtr, prop_key );
        const int_type ipCycle2DataOffset    ( newCycles2Data.first  );
        const int_type ipCycle2FlagOffset    ( newCycles2Data.second );

        // Data and Flag Bounds
        vector<vector<pair<int_type,int_type> > > dataBounds( nipCycles1, vector<pair<int_type,int_type> >(nipCycles2+1,pair<int_type,int_type>(0,0)) );
        for ( int_type cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for ( int_type cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                dataBounds[cycle1][cycle2].first  = dataOffset+cycle1*ipCycle1DataOffset + cycle2*ipCycle2DataOffset;
                dataBounds[cycle1][cycle2].second = dataBounds[cycle1][cycle2].first + dataDepth;
            }
            dataBounds[cycle1][nipCycles2].first  = dataOffset+(cycle1+1)*ipCycle1DataOffset;
            dataBounds[cycle1][nipCycles2].second = dataOffset+(cycle1+1)*ipCycle1DataOffset;
        }
        dataBounds[nipCycles1-1][nipCycles2].first  = static_cast<int_type>(dataSize);
        dataBounds[nipCycles1-1][nipCycles2].second = static_cast<int_type>(dataSize);

        vector<vector<pair<int_type,int_type> > > flagBounds( nipCycles1, vector<pair<int_type,int_type> >(nipCycles2+1,pair<int_type,int_type>(0,0)) );
        for ( int_type cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for ( int_type cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                flagBounds[cycle1][cycle2].first  = flagOffset+cycle1*ipCycle1FlagOffset + cycle2*ipCycle2FlagOffset;
                flagBounds[cycle1][cycle2].second = flagBounds[cycle1][cycle2].first + flagDepth;
            }
            flagBounds[cycle1][nipCycles2].first  = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
            flagBounds[cycle1][nipCycles2].second = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
        }
        flagBounds[nipCycles1-1][nipCycles2].first  = static_cast<int_type>(flagSize);
        flagBounds[nipCycles1-1][nipCycles2].second = static_cast<int_type>(flagSize);


        // Shifting data (same for all types)
        for ( int_type cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for ( int_type cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                const int_type dataStart  ( dataBounds[ cycle1 ][ cycle2     ].first - (cycle1*nipCycles2+cycle2)*dataDepth );
                const int_type dataEnd    ( dataBounds[ cycle1 ][ cycle2 + 1 ].first - (cycle1*nipCycles2+cycle2+1)*dataDepth);
                for ( int_type i=dataStart; i<dataEnd; ++i )
                  data_.data[i] = data_.data[i+(cycle1*nipCycles2+cycle2+1)*dataDepth];
            }
        }
        // Shifting flags (same for all types)
        for ( int_type cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for ( int_type cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                const int_type flagsStart ( flagBounds[ cycle1 ][ cycle2     ].first - (cycle1*nipCycles2+cycle2)*flagDepth );
                const int_type flagsEnd   ( flagBounds[ cycle1 ][ cycle2 + 1 ].first - (cycle1*nipCycles2+cycle2+1)*flagDepth);
                for ( int_type i=flagsStart; i<flagsEnd; ++i )
                    data_.flags[i] = data_.flags[i+(cycle1*nipCycles2+cycle2+1)*flagDepth];
            }
        }

        // trim excessive
        const int_type newDataSize( static_cast<int_type>(dataSize) - nipCycles1*nipCycles2*dataDepth );
        const int_type newFlagSize( static_cast<int_type>(flagSize) - nipCycles1*nipCycles2*flagDepth );
        ResizePropertyStorage( newDataSize, newFlagSize );

     #ifndef NDEBUG
        if( ( prop_key.ipFactorCell + prop_key.ipFactorSector + prop_key.ipFactorFacet ) == 0 )
        {
            if( prop_key.type       == SCALAR )
                --data_.scalars;
            else if( prop_key.type  == VECTOR )
                --data_.vectors;
            else if( prop_key.type  == TENSOR )
                --data_.tensors;
            else if( prop_key.type  == ARRAY )
              {
                --data_.arrays;
                data_.arrayLength -= prop_key.dataDepth;
              }
            else if( prop_key.type  == FLAGGEDARRAY )
              {
                --data_.flaggedArrays;
                data_.flaggedArrayLength -= prop_key.dataDepth;
              }
        }
     #endif

  } // end DeleteProperty



template<uint32_t dim, template<uint32_t> class STOREE>
bool LocalVariableStorage<dim,STOREE>::EmptyLVS() const
  {
    return data_.data.empty();
  }


/**
      revised SKM 5/6/2020 - Michael Liem highlighted output problem; not tested for vectors and tensors yet
      
      @attention this method only outputs variables stored on the Model
*/
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::OutLVS() const
  {
    #ifndef NDEBUG
    const PLACEMENT varPlacement = parsePlacement<dim,STOREE>(); 
    cout <<"\nLocalVariableStorage<" << dim << ">::Out: ";
    cout <<"\n\tstored scalar variables: ";
    if ( data_.scalars > 0U ) {
      csmp::Index  idx(SCALAR,varPlacement,0U);
      while ( idx.index < data_.scalars ) {
        double sc = Read( idx );
        cout << endl <<"\t\t"<< sc;
        idx.index++;
        idx.dataOffset++;
        }
      }

    if ( data_.vectors > 0U ) {
      cout <<"\n\n\tstored vector variables: ";
      csmp::Index  idx(VECTOR,varPlacement,0U);
      VectorVariable<dim>  vc;
      while ( idx.index < data_.vectors ) {
        Read( idx, vc );
        cout << endl <<"\t\t"<< vc;
        idx.index++;
        idx.dataOffset += idx.dataDepth;
        }
      }

    if ( data_.tensors > 0U ) {
      cout <<"\n\n\tstored tensor variables: ";
      const uint32_t  tensors( data_.tensors );
      csmp::Index  idx(TENSOR,varPlacement,0U);
      TensorVariable<dim>  ts;
      while ( idx.index < tensors ) {
        Read( idx, ts );
        cout << endl <<"\t\t"<< ts;
        idx.index++;
        idx.dataOffset += idx.dataDepth;
        }
      }
      
    if ( data_.arrays > 0U ) {
        cout <<"\n\n\tstored array variables: ";
        const uint32_t  array_variables( data_.arrays );
        for ( auto i{0U}; i<array_variables; ++i )
          {
             csmp::Index    idx(ARRAY,varPlacement,i);
             ArrayVariable  ary( idx.dataDepth ); 
             Read( idx, ary );
             cout << endl <<"\t\t"<< ary;
          }
      }
      
    if ( data_.flaggedArrays > 0U ) {
        cout <<"\n\n\tstored array variables: ";
        const uint32_t  flaggged_array_variables( data_.flaggedArrays );
        for ( auto i{0U}; i<flaggged_array_variables; ++i )
          {
             csmp::Index           idx(FLAGGEDARRAY,varPlacement,i);
             FlaggedArrayVariable  ary( idx.dataDepth ); 
             Read( idx, ary );
             cout << endl <<"\t\t"<< ary;
          }
      }

    cout << endl;
    #endif
  }



// ===============
// LOCAL VARIABLES
// ===============

/// Scalar variable value
template<uint32_t dim, template<uint32_t> class STOREE>
double LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx ) const
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
 assert( idx.type == SCALAR );
 assert( idx.index < data_.scalars );
 assert( idx.dataOffset < data_.data.size() );
    return data_.data[idx.dataOffset];
 }
 
 
/// Scalar variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, ScalarVariable& sc ) const  
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
 assert( idx.type == SCALAR );
 assert( idx.index < data_.scalars );
 assert( idx.dataOffset < data_.data.size() );
 assert( idx.flagOffset < data_.flags.size() );
    sc.Flag() = data_.flags[idx.flagOffset];
    sc        = data_.data[idx.dataOffset];
 }


/// Scalar variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const ScalarVariable& sc )  
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
 assert( idx.type == SCALAR );
 assert( idx.index < data_.scalars );
 assert( idx.flagOffset < data_.flags.size() );
 assert( idx.dataOffset < data_.data.size() );

    data_.flags[idx.flagOffset] = sc.Flag();
    data_.data[idx.dataOffset]  = sc();
 }


/// Scalar & Array variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx ) const 
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
 assert( idx.type == SCALAR || idx.type == ARRAY );
 assert( idx.flagOffset < data_.flags.size() );

    return data_.flags[idx.flagOffset];
 }
 

/// Vector, Tensor, FlaggedArray variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx, int_type i ) const
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
  assert( (idx.flagOffset+i) < data_.flags.size() );

    return data_.flags[ idx.flagOffset + i ];
 }
 
 
/// Scalar & Array variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx, VARIABLE_FLAG flag ) 
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
 assert( idx.type == SCALAR  || idx.type == ARRAY );
 assert( idx.flagOffset < data_.flags.size() );
    data_.flags[idx.flagOffset] = flag;
 }

 
/// Vector,Tensor, FlaggedArray variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx, int_type i, VARIABLE_FLAG flag )
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
 assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
 assert( (idx.flagOffset+i) < data_.flags.size() );

    data_.flags[ idx.flagOffset + i ] = flag;
 }
 

/// Vector variable 
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const VectorVariable<dim>& vc )  
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
 assert( idx.type == VECTOR );
 assert( idx.index < data_.vectors );
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim-1) < data_.data.size() );

    for ( uint32_t i{0u}; i<dim; ++i ) {
         data_.flags[ idx.flagOffset+i ] = vc.Flag(i);
         data_.data[ idx.dataOffset+i ]  = vc[i];
      }
 }


/// Vector variable 
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, VectorVariable<dim>& vc ) const  
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
 assert( idx.type == VECTOR );
 assert( idx.index < data_.vectors );
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim-1) < data_.data.size() );

    for ( uint32_t i{0u}; i<dim; ++i ) {
         vc.Flag(i) = data_.flags[ idx.flagOffset+i ];
         vc(i)      = data_.data[ idx.dataOffset+i ];
      }
 }
 

/// Tensor variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const TensorVariable<dim>& ts )
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
  assert( idx.type == TENSOR );
  assert( idx.index < data_.tensors );
  assert( (idx.flagOffset+dim-1) < data_.flags.size() );
  assert( (idx.dataOffset+dim*dim-1) < data_.data.size() );

    const uint32_t dataOffset(idx.dataOffset);
    const uint32_t flagOffset(idx.flagOffset);
    for ( uint32_t i{0U}; i<dim; i++ )
      {
        data_.flags[ flagOffset+i ] = ts.Flag(i);
        for ( uint32_t j{0U}; j<dim; j++ )
          data_.data[ dataOffset+i*dim+j ] = ts(i,j);
      }
 }


/// Tensor variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, TensorVariable<dim>& ts ) const  
 {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
 assert( idx.type == TENSOR );
 assert( idx.index < data_.tensors );
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim*dim-1) < data_.data.size() );

   const uint32_t dataOffset(idx.dataOffset);
   const uint32_t flagOffset(idx.flagOffset);
   for ( uint32_t i{0U}; i<dim; i++ )
     {
       ts.Flag(i) = data_.flags[ flagOffset+i ] ;
       for ( uint32_t j{0U}; j<dim; j++ )
         ts(i,j) = data_.data[ dataOffset+i*dim+j ];
     }
 }


/// Array variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const ArrayVariable& av )
  {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
  assert( idx.type == ARRAY );
  assert( av.Size() == idx.dataDepth );
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset) < data_.flags.size() );

    // Use std::copy for better performance
    std::copy(av.Begin(), av.End(), &data_.data[idx.dataOffset]);
    data_.flags[idx.flagOffset] = av.Flag();
  }


/// Array variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, ArrayVariable& av ) const
  {
    av.Resize( idx.dataDepth );
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
  assert( idx.type == ARRAY );
  assert( av.Size() == idx.dataDepth );
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset) < data_.flags.size() );

    // Use a loop for assignment to avoid const_iterator issues
    for (uint32_t i = 0; i < av.Size(); ++i) {
        av(i) = data_.data[idx.dataOffset + i];
    }
    av.Flag() = data_.flags[idx.flagOffset];
  }


/// FlaggedArray variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const FlaggedArrayVariable& av )
  {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
  assert( idx.type == FLAGGEDARRAY );
  assert( av.Size() == idx.dataDepth ); 
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset+idx.dataDepth-1) < data_.flags.size() );

    // Use std::copy for better performance
    std::copy(av.Begin(), av.End(), &data_.data[idx.dataOffset]);
    // Store flags one by one
    for (uint32_t i = 0; i < av.Size(); ++i) {
        data_.flags[idx.flagOffset + i] = av.Flag(i);
    }
  }


/// FlaggedArray variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, FlaggedArrayVariable& av ) const
  {
  av.Resize( idx.dataDepth );
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
  assert( idx.type == FLAGGEDARRAY );
  assert( av.Size() == idx.dataDepth ); 
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset+idx.dataDepth-1) < data_.flags.size() );

    // Use a loop for assignment to avoid const_iterator issues
    for (uint32_t i = 0; i < av.Size(); ++i) {
        av(i) = data_.data[idx.dataOffset + i];
        av.Flag(i) = data_.flags[idx.flagOffset + i];
    }
  }




// ===================================
// ELEMENT INTEGRATION POINT VARIABLES
// ===================================

#include "LocalVariableStorageIndexArithmetic.h"


/// Scalar variable value at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
double LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::Index& idx ) const
  {
    const int_type offset(DATA_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );

    return data_.data[offset];
  }


/// Scalar variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::Index& idx, ScalarVariable& sc ) const
  {
    const int_type offset(DATA_OFFSET_IP);
    const int_type flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );

    sc.Flag() = data_.flags[flagOffset];
    sc        = data_.data[offset];
  }


/// Scalar variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::Index& idx, const ScalarVariable& sc )
  {
    const int_type offset(DATA_OFFSET_IP);
    const int_type flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );

    data_.flags[flagOffset] = sc.Flag();
    data_.data[offset] = sc();
  }


/// Scalar & Array variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::Index& idx ) const
  {
    const int_type flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == SCALAR || idx.type == ARRAY );
  assert( flagOffset < data_.flags.size() );

    return data_.flags[flagOffset];
  }


/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::Index& idx, int_type i ) const
  {
    const int_type flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
    assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
    assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
    assert( flagOffset+i < data_.flags.size() );

    return data_.flags[ flagOffset+i ];
  }


/// Scalar & Array variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::Index& idx, VARIABLE_FLAG flag )
  {
    const int_type flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == SCALAR  || idx.type == ARRAY );
  assert( flagOffset < data_.flags.size() );

    data_.flags[flagOffset] = flag;
  }


/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::Index& idx, int_type i, VARIABLE_FLAG flag )
  {
    const int_type flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
  assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
  assert( flagOffset+i < data_.flags.size() );

    data_.flags[ flagOffset+i ] = flag;
  }


/// Vector variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::Index& idx, const VectorVariable<dim>& vc )
  {
    const int_type offset(DATA_OFFSET_IP);
    const int_type flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );

  for ( uint32_t i(0); i<dim; ++i ) {
    data_.flags[ flagOffset+i ] = vc.Flag(i);
    data_.data[ offset+i ]  = vc[i];
    }
  }


/// Vector variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::Index& idx, VectorVariable<dim>& vc ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );

    for ( uint32_t i(0u); i<dim; ++i ) {
      vc.Flag(i) = data_.flags[ flagOffset+i ];
      vc(i)      = data_.data[ offset+i ];
      }
  }


/// Tensor variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::Index& idx, const TensorVariable<dim>& ts )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );

  for ( uint32_t i{0U}; i<dim; i++ )
    {
    data_.flags[ flagOffset+i ] = ts.Flag(i);
    for ( uint32_t j{0U}; j<dim; j++ )
      data_.data[ offset+i*dim+j ] = ts(i,j);
    }
  }


/// Tensor variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::Index& idx, TensorVariable<dim>& ts ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );

  for ( uint32_t i{0U}; i<dim; i++ )
    {
    ts.Flag(i) = data_.flags[flagOffset+i];
    for ( uint32_t j{0U}; j<dim; j++ )
      ts(i,j) = data_.data[ offset+i*dim+j ];
    }
  }

/// Array variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::Index& idx, const ArrayVariable& av )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);
    const uint32_t arraySize( idx.dataDepth );

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );

    for( uint32_t i(0); i < arraySize; ++i )
      data_.data[ offset     + i ] = av[i];
    data_.flags[flagOffset] = av.Flag();
  }


/// Array variables
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::Index& idx, ArrayVariable& av ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);
    const uint32_t arraySize( idx.dataDepth );

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset < data_.flags.size() );

    for( uint32_t i(0); i < arraySize; ++i )
      av(i) = data_.data[ offset+i ];

    av.Flag( data_.flags[flagOffset] );
  }


/// FlaggedArray variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::Index& idx, const FlaggedArrayVariable& av )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);
    const uint32_t arraySize( idx.dataDepth );

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );

    for( uint32_t i(0); i < arraySize; ++i )
    {
      data_.data[ offset     + i ] = av[i];
      data_.flags[flagOffset + i ] = av.Flag(i);
    }
  }


/// FlaggedArray variables
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::Index& idx, FlaggedArrayVariable& av ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);
    const uint32_t arraySize( idx.dataDepth );

#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );

    for( uint32_t i(0); i < arraySize; ++i )
    {
      av(i) = data_.data[ offset + i ];
      av.Flag( i, data_.flags[flagOffset+i] );
    }
  }



/**
    @todo replace this super wasteful method
*/
template<uint32_t dim, template<uint32_t> class STOREE>
bool LocalVariableStorage<dim,STOREE>::IsWithinRange( const csmp::Index& idx,
                                                      double vmin, double vmax ) const
  {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
 AssertPlacement(idx);
#endif
    if ( idx.type == SCALAR ) {
      const double val = Read( idx );
      return ( val >= vmin and val <= vmax ) ? true : false;
      }
    if ( idx.type == VECTOR ) {
      VectorVariable<dim> vc = ReadVector( idx );
      return vc.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == TENSOR ) {
      TensorVariable<dim>  ts = ReadTensor( idx );
      return ts.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == ARRAY ) {
      vector<double> av = ReadArray( idx );
      auto [min_it, max_it] = minmax_element(av.begin(), av.end());
      return *min_it >= vmin && *max_it <= vmax;
      }
    if ( idx.type == FLAGGEDARRAY ) {
      FlaggedArrayVariable  fa;
      Read( idx, fa );
      return fa.IsWithinRange( vmin, vmax );
      }

    std::cerr <<"\nLocalVariableStorage<dim,STOREE>::IsWithinRange: range check could not be performed."<< std::endl;
    return false;
  }


template<uint32_t dim, template<uint32_t> class STOREE>
VectorVariable<dim> LocalVariableStorage<dim,STOREE>::ReadVector( uint32_t ip, const csmp::Index& idx ) const
 {
    const int_type dataOffset(DATA_OFFSET_IP);
    const int_type flagOffset(FLAG_OFFSET_IP);

#if !defined(NDEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == VECTOR );
  assert( dataOffset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );

    if constexpr( dim == 3U )
      return VectorVariable<3U>( data_.flags[ flagOffset ], data_.flags[ flagOffset+1 ], data_.flags[ flagOffset+2 ],
                                 data_.data[ dataOffset ], data_.data[ dataOffset+1 ], data_.data[ dataOffset+2 ] );

    else if constexpr( dim == 2U )
      return VectorVariable<2U>( data_.flags[ flagOffset ], data_.flags[ flagOffset+1 ],
                                 data_.data[ dataOffset ], data_.data[ dataOffset+1 ] );

    else if constexpr( dim == 1U )
      return VectorVariable<1U>( data_.flags[ flagOffset ], data_.data[ dataOffset ] );
 }


template<uint32_t dim, template<uint32_t> class STOREE>
TensorVariable<dim> LocalVariableStorage<dim,STOREE>::ReadTensor( uint32_t ip, const csmp::Index& idx ) const
 {
    const int_type dataOffset(DATA_OFFSET_IP);
    const int_type flagOffset(FLAG_OFFSET_IP);

#if !defined(NDEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == TENSOR );
  assert( dataOffset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );

    if constexpr ( dim == 3U )
      return TensorVariable<3U>( data_.flags[ flagOffset ], data_.flags[ flagOffset+1 ], data_.flags[ flagOffset+2 ],
                                 data_.data[ dataOffset ], data_.data[ dataOffset+1 ], data_.data[ dataOffset+2 ],
                                 data_.data[ dataOffset+dim ], data_.data[ dataOffset+dim+1 ], data_.data[ dataOffset+dim+2 ],
                                 data_.data[ dataOffset+2*dim ], data_.data[ dataOffset+2*dim+1 ], data_.data[ dataOffset+2*dim+2 ] );

    else if constexpr ( dim == 2U )
      return TensorVariable<2U>( data_.flags[ flagOffset ], data_.flags[ flagOffset+1 ],
                                 data_.data[ dataOffset ], data_.data[ dataOffset+1 ],
                                 data_.data[ dataOffset+dim ], data_.data[ dataOffset+dim+1 ] );

    else if constexpr ( dim == 1U )
      return TensorVariable<1U>( data_.flags[ flagOffset ], data_.data[ dataOffset ] );
 }


template<uint32_t dim, template<uint32_t> class STOREE>
vector<double> LocalVariableStorage<dim,STOREE>::ReadArray( uint32_t ip, const csmp::Index& idx ) const
 {
    const int_type dataOffset(DATA_OFFSET_IP);
    const int_type arraySize( idx.dataDepth );

#if !defined(NDEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
  assert( idx.type == ARRAY );
  assert( dataOffset+arraySize <= data_.data.size() );

    return std::vector<double>( next(data_.data.begin(),dataOffset), next(data_.data.begin(),dataOffset+arraySize) );
 }



template<uint32_t dim, template<uint32_t> class STOREE>
bool LocalVariableStorage<dim,STOREE>::IsWithinRange( uint32_t ip, const csmp::Index& idx,
                                                      double vmin, double vmax ) const
  {
#if defined(DEBUG) && defined(CSMP_VARIABLE_STORAGE_DEBUG)
  AssertIntegrationPointPlacement(idx);
#endif
    if ( idx.type == SCALAR ) {
      const double val = Read( ip, idx );
      return ( val >= vmin and val <= vmax ) ? true : false;
      }
    if ( idx.type == VECTOR ) {
      VectorVariable<dim> vc = ReadVector( ip, idx );
      return vc.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == TENSOR ) {
      TensorVariable<dim>  ts = ReadTensor( ip, idx );
      return ts.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == ARRAY ) {
      vector<double> av = ReadArray( ip, idx );
      auto [min_it, max_it] = minmax_element(av.begin(), av.end());
      return *min_it >= vmin && *max_it <= vmax;
      }
    if ( idx.type == FLAGGEDARRAY ) {
      FlaggedArrayVariable  fa;
      Read( ip, idx, fa );
      return fa.IsWithinRange( vmin, vmax );
      }
    std::cout <<"\nLocalVariableStorage<dim,STOREE>::IsWithinRange: range check could not be performed."<< std::endl;
    return false;
  }







// =========================================
// FINITE VOLUME INTEGRATION POINT VARIABLES
// =========================================


/** 
    Scalar variable value at sector or facet integration points only.
*/
template<uint32_t dim, template<uint32_t> class STOREE>
double LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    // offset to first instance of idx variable in the data vector
    const pair<uint32_t, uint32_t> offsetData = localVariableDispatch::containerOffset( storeePtr, idx);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offsetData.first < data_.data.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth; // for facet integration points only
    ///                                                                                    ^^^^^^^^
    return data_.data[offsetData.first + sector_ip_offset];
  }



/// Scalar variable at facet or sector integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, ScalarVariable& sc ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<uint32_t, uint32_t> offsetData = localVariableDispatch::containerOffset( storeePtr, idx);
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif
    const uint32_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    sc.Flag() = data_.flags[flagOffset + sector_ip_flag_offset];
    sc        = data_.data[offset + sector_ip_offset];
  }




/// Scalar variable at facet or sector integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, const ScalarVariable& sc )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<uint32_t, uint32_t> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif
    const uint32_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    data_.flags[flagOffset + sector_ip_flag_offset] = sc.Flag();
    data_.data[offset + sector_ip_offset]           = sc();
  }





/// Scalar & Array variable flag at facet or sector integration point
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<uint32_t, uint32_t> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR || idx.type == ARRAY );
  assert( offsetData.second < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;
    
    return data_.flags[offsetData.second + sector_ip_offset];
  }




/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, int_type i ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type,int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx);
    const int_type flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
  assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
  assert( flagOffset+i < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;
    
    return data_.flags[ flagOffset + sector_ip_offset + i ];
  }




/// Scalar & Array variable flag at facet or sector integration point (note that the Array has only a single flag)
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, VARIABLE_FLAG flag )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR  || idx.type == ARRAY );
  assert( flagOffset < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    data_.flags[flagOffset + sector_ip_offset] = flag;
  }



/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, int_type i, VARIABLE_FLAG flag )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<uint32_t, uint32_t> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
  assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
  assert( flagOffset+i < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    data_.flags[ flagOffset + sector_ip_offset + i ] = flag;
  }




/// Vector variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, const VectorVariable<dim>& vc )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type offset(offsetData.first);
    const int_type flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const int_type sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const int_type sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i(0); i<dim; ++i ) {
         data_.flags[ flagOffset + sector_ip_flag_offset + i ] = vc.Flag(i);
         data_.data[ offset + sector_ip_offset + i ]            = vc[i];
      }
  }




/// Vector variable at facet or sector integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, VectorVariable<dim>& vc ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type offset(offsetData.first);
    const int_type flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const int_type sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const int_type sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i(0); i<dim; ++i ) {
         vc.Flag(i) = data_.flags[ flagOffset + sector_ip_flag_offset + i ];
         vc(i)      = data_.data[ offset + sector_ip_offset + i ];
      }
  }




/// Tensor variable at sector or facet integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, const TensorVariable<dim>& ts )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<uint32_t, uint32_t> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const uint32_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( auto i{0U}; i<dim; ++i ) {
          data_.flags[ flagOffset + sector_ip_flag_offset + i ] = ts.Flag(i);
          for ( uint32_t j{0U}; j<dim; ++j )
            data_.data[ offset + sector_ip_offset + i*dim + j ] = ts(i,j);
      }
 }


/// Tensor variable at facet or sector integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, TensorVariable<dim>& ts ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type offset(offsetData.first);
    const int_type flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const int_type sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const int_type sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i{0U}; i<dim; ++i ) {
          ts.Flag(i) = data_.flags[flagOffset + sector_ip_flag_offset + i];
          for ( uint32_t j{0U}; j<dim; ++j )
            ts(i,j) = data_.data[ offset + sector_ip_offset + i*dim + j ];
      }
  }


/// Array variable at sector or facet integration poin
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, const ArrayVariable& av )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type offset(offsetData.first);
    const int_type flagOffset(offsetData.second);

    const int_type arraySize( idx.dataDepth );

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif
    const int_type sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const int_type sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i(0u); i < arraySize; ++i )
      data_.data[ offset + sector_ip_offset + i ]     = av[i];
    data_.flags[ flagOffset + sector_ip_flag_offset ] = av.Flag();
  }





/// Array variables at sector or facet integration points
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, ArrayVariable& av ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type offset(offsetData.first);
    const int_type flagOffset(offsetData.second);

    const int_type arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif
    const int_type sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const int_type sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i(0); i < arraySize; ++i )
      av(i) = data_.data[ offset + sector_ip_offset + i ];
    av.Flag( data_.flags[flagOffset + sector_ip_flag_offset] );
  }



/// FlaggedArray variable at facet or sector integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, const FlaggedArrayVariable& av )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type offset(offsetData.first);
    const int_type flagOffset(offsetData.second);

    const int_type arraySize( idx.dataDepth );

#ifndef NDEBUG  
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif
    const int_type sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const int_type sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i(0u); i < arraySize; ++i ) {
         data_.data[ offset     + sector_ip_offset + i ]      = av[i];
         data_.flags[flagOffset + sector_ip_flag_offset + i ] = av.Flag(i);
      }
  }



/// FlaggedArray variables
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, FlaggedArrayVariable& av ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type offset(offsetData.first);
    const int_type flagOffset(offsetData.second);

    const int_type arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif
    const int_type sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const int_type sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i(0); i < arraySize; ++i ) {
         av(i) = data_.data[ offset + sector_ip_offset + i ];
         av.Flag( i, data_.flags[flagOffset + sector_ip_flag_offset + i] );
      }
  }


template<uint32_t dim, template<uint32_t> class STOREE>
VectorVariable<dim> LocalVariableStorage<dim,STOREE>::ReadVector( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx ) const
 {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type dataOffset(offsetData.first);
    const int_type flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( dataOffset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const int_type ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                     (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                     (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const int_type ip_data_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;
    if constexpr( dim == 3U )
      return VectorVariable<3U>( data_.flags[ flagOffset+ip_flag_offset ], data_.flags[ flagOffset+ip_flag_offset+1 ], data_.flags[ flagOffset+ip_flag_offset+2 ],
                                 data_.data[ dataOffset+ip_data_offset ], data_.data[ dataOffset+ip_data_offset+1 ], data_.data[ dataOffset+ip_data_offset+2 ] );

    else if constexpr( dim == 2U )
      return VectorVariable<2U>( data_.flags[ flagOffset+ip_flag_offset ], data_.flags[ flagOffset+ip_flag_offset+1 ],
                                 data_.data[ dataOffset+ip_data_offset ], data_.data[ dataOffset+ip_data_offset+1 ] );

    else if constexpr( dim == 1U )
      return VectorVariable<1U>( data_.flags[ flagOffset+ip_flag_offset ], data_.data[ dataOffset+ip_data_offset ] );
 }


template<uint32_t dim, template<uint32_t> class STOREE>
TensorVariable<dim> LocalVariableStorage<dim,STOREE>::ReadTensor( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx ) const
 {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type dataOffset(offsetData.first);
    const int_type flagOffset(offsetData.second);

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == TENSOR );
  assert( dataOffset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const int_type ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                     (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                     (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const int_type ip_data_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;
    if constexpr ( dim == 3U )
      return TensorVariable<3U>( data_.flags[ flagOffset+ip_flag_offset ], data_.flags[ flagOffset+ip_flag_offset+1 ], data_.flags[ flagOffset+ip_flag_offset+2 ],
                                 data_.data[ dataOffset+ip_data_offset ], data_.data[ dataOffset+ip_data_offset+1 ], data_.data[ dataOffset+ip_data_offset+2 ],
                                 data_.data[ dataOffset+ip_data_offset+dim ], data_.data[ dataOffset+ip_data_offset+dim+1 ], data_.data[ dataOffset+ip_data_offset+dim+2 ],
                                 data_.data[ dataOffset+ip_data_offset+2*dim ], data_.data[ dataOffset+ip_data_offset+2*dim+1 ], data_.data[ dataOffset+ip_data_offset+2*dim+2 ] );

    else if constexpr ( dim == 2U )
      return TensorVariable<2U>( data_.flags[ flagOffset+ip_flag_offset ], data_.flags[ flagOffset+ip_flag_offset+1 ],
                                 data_.data[ dataOffset+ip_data_offset ], data_.data[ dataOffset+ip_data_offset+1 ],
                                 data_.data[ dataOffset+ip_data_offset+dim ], data_.data[ dataOffset+ip_data_offset+dim+1 ] );

    else if constexpr ( dim == 1U )
      return TensorVariable<1U>( data_.flags[ flagOffset+ip_flag_offset ], data_.data[ dataOffset+ip_data_offset ] );
}


template<uint32_t dim, template<uint32_t> class STOREE>
vector<double> LocalVariableStorage<dim,STOREE>::ReadArray( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx ) const
 {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const pair<int_type, int_type> offsetData = localVariableDispatch::containerOffset( storeePtr, idx );
    const int_type dataOffset(offsetData.first);
    const int_type arraySize( idx.dataDepth );

#ifndef NDEBUG
  localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == ARRAY );
  assert( dataOffset+arraySize <= data_.data.size() );
#endif
    const int_type ip_data_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    return std::vector<double>( next(data_.data.begin(),dataOffset+ip_data_offset), next(data_.data.begin(),dataOffset+ip_data_offset+arraySize) );
 }




/**
    checks whether any value of the target variable is within the range supplied via the arguments vmin and vmax.
    
     @attention for a vector, tensor, or array this range check is performed on all their elements.
     
     @todo refactor: for vectors, tensors and array variables this method is terribly inefficient as it creates temporaries for the checking.
*/
template<uint32_t dim, template<uint32_t> class STOREE>
bool LocalVariableStorage<dim,STOREE>::IsWithinRange( uint32_t sector_or_facet, uint32_t ip,
                                                      const csmp::Index& idx,
                                                      double vmin, double vmax ) const
  {
#ifndef NDEBUG
   const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
   localVariableDispatch::assertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
   AssertFiniteVolumeIntegrationPointPlacement(idx);
#endif
    if ( idx.type == SCALAR ) {
      const double val = Read( sector_or_facet, ip, idx );
      return ( val >= vmin and val <= vmax ) ? true : false;
      }
    if ( idx.type == VECTOR ) {
      VectorVariable<dim> vc = ReadVector( sector_or_facet, ip, idx );
      return vc.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == TENSOR ) {
      TensorVariable<dim>  ts = ReadTensor( sector_or_facet, ip, idx );
      return ts.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == ARRAY ) {
      vector<double> av = ReadArray( sector_or_facet, ip, idx );
      auto [min_it, max_it] = minmax_element(av.begin(), av.end());
      return *min_it >= vmin && *max_it <= vmax;
      }
    if ( idx.type == FLAGGEDARRAY ) {
      FlaggedArrayVariable  fa;
      Read( sector_or_facet, ip, idx, fa );
      return fa.IsWithinRange( vmin, vmax );
      }

    cout <<"\nLocalVariableStorage<dim,STOREE>::IsWithinRange: range check could not be performed."<< endl;
    return false;
  }
  
  
  
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::StoreArrayEntry( const csmp::Index& idx, double value, int_type array_elmt )
  {
    assert( idx.type == ARRAY );
#ifndef NDEBUG
  assert( array_elmt <= idx.dataDepth );
#endif
    data_.data[ idx.dataOffset + array_elmt ] = value;
 }


    // SKM 1/11/23
/// reads single double element inside of the target array
template<uint32_t dim, template<uint32_t> class STOREE>
double LocalVariableStorage<dim,STOREE>::ReadArrayEntry( const csmp::Index& idx, int_type array_elmt ) const
 {
    assert( idx.type == ARRAY );
#ifndef NDEBUG
  assert( array_elmt <= idx.dataDepth );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
#endif
    return data_.data[ idx.dataOffset + array_elmt ];
 } // end Read (Array[array_elmt]





template class LocalVariableStorage<1U,Node>;
template class LocalVariableStorage<1U,Element>;
template class LocalVariableStorage<1U,Face>;
template class LocalVariableStorage<1U,InterFace>;
template class LocalVariableStorage<1U,Edge>;
template class LocalVariableStorage<1U,Region>;
template class LocalVariableStorage<1U,Boundary>;
template class LocalVariableStorage<1U,SplitBoundary>;
template class LocalVariableStorage<1U,Model>;

template class LocalVariableStorage<2U,Node>;
template class LocalVariableStorage<2U,Element>;
template class LocalVariableStorage<2U,Face>;
template class LocalVariableStorage<2U,InterFace>;
template class LocalVariableStorage<2U,Edge>;
template class LocalVariableStorage<2U,Region>;
template class LocalVariableStorage<2U,Boundary>;
template class LocalVariableStorage<2U,SplitBoundary>;
template class LocalVariableStorage<2U,Model>;

template class LocalVariableStorage<3U,Node>;
template class LocalVariableStorage<3U,Element>;
template class LocalVariableStorage<3U,Face>;
template class LocalVariableStorage<3U,InterFace>;
template class LocalVariableStorage<3U,Edge>;
template class LocalVariableStorage<3U,Region>;
template class LocalVariableStorage<3U,Boundary>;
template class LocalVariableStorage<3U,SplitBoundary>;
template class LocalVariableStorage<3U,Model>;

} // end csmp
