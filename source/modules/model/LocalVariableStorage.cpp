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

namespace csmp {


// ---------------------------------------------------------------
//       CONSTRUCTORS
// ---------------------------------------------------------------



/// Ctor for plain(empty) variable storage
template<uint32_t dim, template<uint32_t> class STOREE>
LocalVariableStorage<dim, STOREE>::LocalVariableStorage()
    : data_()
{
}


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
LocalVariableStorage<dim, STOREE>::~LocalVariableStorage()
{
}


template<uint32_t dim, template<uint32_t> class STOREE>
LocalVariableStorage<dim, STOREE>::LocalVariableStorage( const LocalVariableStorage<dim,STOREE>& ps )
    : data_(ps.data_)
{
}


template<uint32_t dim, template<uint32_t> class STOREE>
LocalVariableStorage<dim, STOREE>& LocalVariableStorage<dim, STOREE>::operator=( const LocalVariableStorage& ps )
{
    if ( &ps != this )
      data_.data = ps.data_.data;

    return *this;
}



/// Goto overload for storees without IntegrationPointVariables
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::ResizePropertyStorage( const LocalVariables& lv )
  {
    ResizePropertyStorage( lv.totalDataDepth , lv.totalFlagDepth );

// initialise auxiliary parameters for debugging
#ifndef NDEBUG
    StoreLocalState(lv);
#endif
  }


/// Goto overload for storees with IntegrationPointVariables
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::ResizePropertyStorage( const LocalVariables& lv, const IntegrationPointVariables& ipv )
  {
    const std::pair<uint32_t, uint32_t> newContainerSize = lvsCompileTimeDispatch::containerNewSize( static_cast<const STOREE<dim>*>(this), lv, ipv );
    ResizePropertyStorage(newContainerSize.first, newContainerSize.second);

// initialise auxiliary parameters for debugging
#ifndef NDEBUG
    StoreLocalState(lv);
#endif
  }


/// Resizes the storage preserving original values if any (this is where all ResizePropertyStorage end up)
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::ResizePropertyStorage( uint32_t newDataComponentCount, uint32_t newFlagComponentCount )
  {
    // needed?
    if( data_.flags.size() == newFlagComponentCount && data_.data.size() == newDataComponentCount )
      return;

    // resizing
    data_.flags.resize( newFlagComponentCount, ANY );
    data_.data.resize( newDataComponentCount, std::numeric_limits<double>::quiet_NaN() );

    // trimming excess capacity
    std::vector<VARIABLE_FLAG>( data_.flags ).swap( data_.flags );
    std::vector<double>( data_.data ).swap( data_.data );
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
    const std::pair<uint32_t,uint32_t> newTotaDataDepth = lvsCompileTimeDispatch::containerTotalDataDepth( storeePtr, prop_key );
    ResizePropertyStorage( newTotaDataDepth.first, newTotaDataDepth.second );
    const uint32_t dataSize = static_cast<uint32_t>(data_.data.size());
    const uint32_t flagSize = static_cast<uint32_t>(data_.flags.size());

    /// Roman, 2013: with ipvs support

    // Offset
    const std::pair<uint32_t,uint32_t> newOffset = lvsCompileTimeDispatch::containerOffset( storeePtr, prop_key );
    const uint32_t dataOffset         ( newOffset.first       );
    const uint32_t flagOffset         ( newOffset.second      );

    // DataDepth
    const uint32_t dataDepth          ( prop_key.dataDepth    );
    const uint32_t flagDepth          ( prop_key.flagDepth    );

    // Cycles
    const std::pair<int,int> newCycles = lvsCompileTimeDispatch::containerIPCycles( storeePtr, prop_key );
    const int nipCycles1            ( newCycles.first       );
    const int nipCycles2            ( newCycles.second      );

    const std::pair<int,int> newCycles1Data = lvsCompileTimeDispatch::containerIPCycle1Offset( storeePtr, prop_key );
    const int ipCycle1DataOffset    ( newCycles1Data.first  );
    const int ipCycle1FlagOffset    ( newCycles1Data.second );

    const std::pair<int,int> newCycles2Data = lvsCompileTimeDispatch::containerIPCycle2Offset( storeePtr, prop_key );
    const int ipCycle2DataOffset    ( newCycles2Data.first  );
    const int ipCycle2FlagOffset    ( newCycles2Data.second );

    // Data and Flag Bounds
    std::vector<std::vector<std::pair<uint32_t,uint32_t> > > dataBounds( nipCycles1, std::vector<std::pair<uint32_t,uint32_t> >(nipCycles2+1,std::pair<uint32_t,uint32_t>(0,0)) );
    for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
    {
        for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
        {
            dataBounds[cycle1][cycle2].first  = dataOffset+cycle1*ipCycle1DataOffset + cycle2*ipCycle2DataOffset;
            dataBounds[cycle1][cycle2].second = dataBounds[cycle1][cycle2].first + dataDepth;
        }
        dataBounds[cycle1][nipCycles2].first  = dataOffset+(cycle1+1)*ipCycle1DataOffset;
        dataBounds[cycle1][nipCycles2].second = dataOffset+(cycle1+1)*ipCycle1DataOffset;
    }
    dataBounds[nipCycles1-1][nipCycles2].first  = dataSize;
    dataBounds[nipCycles1-1][nipCycles2].second = dataSize;

    std::vector<std::vector<std::pair<uint32_t,uint32_t> > > flagBounds( nipCycles1, std::vector<std::pair<uint32_t,uint32_t> >(nipCycles2+1,std::pair<uint32_t,uint32_t>(0,0)) );
    for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
    {
        for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
        {
            flagBounds[cycle1][cycle2].first  = flagOffset+cycle1*ipCycle1FlagOffset + cycle2*ipCycle2FlagOffset;
            flagBounds[cycle1][cycle2].second = flagBounds[cycle1][cycle2].first + flagDepth;
        }
        flagBounds[cycle1][nipCycles2].first  = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
        flagBounds[cycle1][nipCycles2].second = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
    }
    flagBounds[nipCycles1-1][nipCycles2].first  = flagSize;
    flagBounds[nipCycles1-1][nipCycles2].second = flagSize;

    // Moving data (same for all types)
    for (int cycle1=nipCycles1-1U; cycle1>=0; --cycle1  ){
        for (int cycle2=nipCycles2-1U; cycle2>=0; --cycle2  ){
            const uint32_t dataStart  ( dataBounds[ cycle1 ][ cycle2     ].second - 1 );
            const uint32_t dataEnd    ( dataBounds[ cycle1 ][ cycle2 + 1 ].second - 1 );
            for ( uint32_t i=dataEnd; i>dataStart; --i )
              data_.data[i] = data_.data[i-(cycle1+1)*(cycle2+1)*dataDepth];
        }
    }
    // Moving flags (same for all types)
    for (int cycle1=nipCycles1-1U; cycle1>=0; --cycle1  ){
        for (int cycle2=nipCycles2-1U; cycle2>=0; --cycle2  ){
            const uint32_t flagsStart ( flagBounds[ cycle1 ][ cycle2     ].second - 1 );
            const uint32_t flagsEnd   ( flagBounds[ cycle1 ][ cycle2 + 1 ].second - 1 );
            for ( uint32_t i=flagsEnd; i>flagsStart; --i )
                data_.flags[i] = data_.flags[i-(cycle1+1)*(cycle2+1)*flagDepth];
        }
    }

    // Adding non-initialized data of new variable
    for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            for ( uint32_t i=dataBounds[cycle1][cycle2].first; i<dataBounds[cycle1][cycle2].second; ++i )
                data_.data[i] = std::numeric_limits<double>::quiet_NaN();
    // Adding non-initialized flags of new variable
    for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            for ( uint32_t i=flagBounds[cycle1][cycle2].first; i<flagBounds[cycle1][cycle2].second; ++i )
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
    const uint32_t flagSize = static_cast<uint32_t>(data_.flags.size());
    const uint32_t dataSize = static_cast<uint32_t>(data_.data.size());

        /// Roman, 2013: with ipvs support

        // Offset
        const std::pair<uint32_t,uint32_t> newOffset = lvsCompileTimeDispatch::containerOffset( storeePtr, prop_key );
        const uint32_t dataOffset         ( newOffset.first       );
        const uint32_t flagOffset         ( newOffset.second      );

        // DataDepth
        const uint32_t dataDepth          ( prop_key.dataDepth    );
        const uint32_t flagDepth          ( prop_key.flagDepth    );

        // Cycles
        const std::pair<int,int> newCycles = lvsCompileTimeDispatch::containerIPCycles( storeePtr, prop_key );
        const int nipCycles1            ( newCycles.first       );
        const int nipCycles2            ( newCycles.second      );

        const std::pair<int,int> newCycles1Data = lvsCompileTimeDispatch::containerIPCycle1Offset( storeePtr, prop_key );
        const int ipCycle1DataOffset    ( newCycles1Data.first  );
        const int ipCycle1FlagOffset    ( newCycles1Data.second );

        const std::pair<int,int> newCycles2Data = lvsCompileTimeDispatch::containerIPCycle2Offset( storeePtr, prop_key );
        const int ipCycle2DataOffset    ( newCycles2Data.first  );
        const int ipCycle2FlagOffset    ( newCycles2Data.second );

        // Data and Flag Bounds
        std::vector<std::vector<std::pair<uint32_t,uint32_t> > > dataBounds( nipCycles1, std::vector<std::pair<uint32_t,uint32_t> >(nipCycles2+1,std::pair<uint32_t,uint32_t>(0,0)) );
        for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                dataBounds[cycle1][cycle2].first  = dataOffset+cycle1*ipCycle1DataOffset + cycle2*ipCycle2DataOffset;
                dataBounds[cycle1][cycle2].second = dataBounds[cycle1][cycle2].first + dataDepth;
            }
            dataBounds[cycle1][nipCycles2].first  = dataOffset+(cycle1+1)*ipCycle1DataOffset;
            dataBounds[cycle1][nipCycles2].second = dataOffset+(cycle1+1)*ipCycle1DataOffset;
        }
        dataBounds[nipCycles1-1][nipCycles2].first  = dataSize;
        dataBounds[nipCycles1-1][nipCycles2].second = dataSize;

        std::vector<std::vector<std::pair<uint32_t,uint32_t> > > flagBounds( nipCycles1, std::vector<std::pair<uint32_t,uint32_t> >(nipCycles2+1,std::pair<uint32_t,uint32_t>(0,0)) );
        for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                flagBounds[cycle1][cycle2].first  = flagOffset+cycle1*ipCycle1FlagOffset + cycle2*ipCycle2FlagOffset;
                flagBounds[cycle1][cycle2].second = flagBounds[cycle1][cycle2].first + flagDepth;
            }
            flagBounds[cycle1][nipCycles2].first  = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
            flagBounds[cycle1][nipCycles2].second = flagOffset+(cycle1+1)*ipCycle1FlagOffset;
        }
        flagBounds[nipCycles1-1][nipCycles2].first  = flagSize;
        flagBounds[nipCycles1-1][nipCycles2].second = flagSize;


        // Shifting data (same for all types)
        for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                const uint32_t dataStart  ( dataBounds[ cycle1 ][ cycle2     ].first - (cycle1*nipCycles2+cycle2)*dataDepth );
                const uint32_t dataEnd    ( dataBounds[ cycle1 ][ cycle2 + 1 ].first - (cycle1*nipCycles2+cycle2+1)*dataDepth);
                for ( uint32_t i=dataStart; i<dataEnd; ++i )
                  data_.data[i] = data_.data[i+(cycle1*nipCycles2+cycle2+1)*dataDepth];
            }
        }
        // Shifting flags (same for all types)
        for (int cycle1=0; cycle1<nipCycles1; ++cycle1  )
        {
            for (int cycle2=0; cycle2<nipCycles2; ++cycle2  )
            {
                const uint32_t flagsStart ( flagBounds[ cycle1 ][ cycle2     ].first - (cycle1*nipCycles2+cycle2)*flagDepth );
                const uint32_t flagsEnd   ( flagBounds[ cycle1 ][ cycle2 + 1 ].first - (cycle1*nipCycles2+cycle2+1)*flagDepth);
                for ( uint32_t i=flagsStart; i<flagsEnd; ++i )
                    data_.flags[i] = data_.flags[i+(cycle1*nipCycles2+cycle2+1)*flagDepth];
            }
        }

        // trim excessive
        const uint32_t newDataSize ( dataSize - nipCycles1*nipCycles2*dataDepth );
        const uint32_t newFlagSize ( flagSize - nipCycles1*nipCycles2*flagDepth );
        ResizePropertyStorage( newDataSize, newFlagSize );

     #ifndef NDEBUG
        if( ( prop_key.ipFactorSimplex + prop_key.ipFactorSector + prop_key.ipFactorFacet ) == 0 )
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
    std::cout <<"\nLocalVariableStorage<" << dim << ">::Out: ";
    std::cout <<"\n\tstored scalar variables: ";
    if ( data_.scalars > 0U ) {
      csmp::Index  idx(SCALAR,varPlacement,0U);
      while ( idx.index < data_.scalars ) {
        double sc = Read( idx );
        std::cout << std::endl <<"\t\t"<< sc;
        idx.index++;
        idx.dataOffset++;
        }
      }

    if ( data_.vectors > 0U ) {
      std::cout <<"\n\n\tstored vector variables: ";
      csmp::Index  idx(VECTOR,varPlacement,0U);
      VectorVariable<dim>  vc;
      while ( idx.index < data_.vectors ) {
        Read( idx, vc );
        std::cout << std::endl <<"\t\t"<< vc;
        idx.index++;
        idx.dataOffset += idx.dataDepth;
        }
      }

    if ( data_.tensors > 0U ) {
      std::cout <<"\n\n\tstored tensor variables: ";
      const uint32_t  tensors( data_.tensors );
      csmp::Index  idx(TENSOR,varPlacement,0U);
      TensorVariable<dim>  ts;
      while ( idx.index < tensors ) {
        Read( idx, ts );
        std::cout << std::endl <<"\t\t"<< ts;
        idx.index++;
        idx.dataOffset += idx.dataDepth;
        }
      }
      
    if ( data_.arrays > 0U ) {
        std::cout <<"\n\n\tstored array variables: ";
        const uint32_t  array_variables( data_.arrays );
        for ( auto i{0U}; i<array_variables; ++i )
          {
             csmp::Index    idx(ARRAY,varPlacement,i);
             ArrayVariable  ary( idx.dataDepth ); 
             Read( idx, ary );
             std::cout << std::endl <<"\t\t"<< ary;
          }
      }
      
    if ( data_.flaggedArrays > 0U ) {
        std::cout <<"\n\n\tstored array variables: ";
        const uint32_t  flaggged_array_variables( data_.flaggedArrays );
        for ( auto i{0U}; i<flaggged_array_variables; ++i )
          {
             csmp::Index           idx(FLAGGEDARRAY,varPlacement,i);
             FlaggedArrayVariable  ary( idx.dataDepth ); 
             Read( idx, ary );
             std::cout << std::endl <<"\t\t"<< ary;
          }
      }

    std::cout << std::endl;
    #endif
  }






















// ===============
// LOCAL VARIABLES
// ===============

/// Scalar variable value
template<uint32_t dim, template<uint32_t> class STOREE>
double LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx ) const  
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR );
 assert( idx.index < data_.scalars );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.dataOffset < data_.data.size() );
#endif
    return data_.data[idx.dataOffset];
 }
 
 
/// Scalar variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, ScalarVariable& sc ) const  
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR );
 assert( idx.index < data_.scalars );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.dataOffset < data_.data.size() );
 assert( idx.flagOffset < data_.flags.size() );
#endif
    sc.Flag() = data_.flags[idx.flagOffset];
    sc        = data_.data[idx.dataOffset];
 }


/// Scalar variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const ScalarVariable& sc )  
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR );
 assert( idx.index < data_.scalars );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
 assert( idx.dataOffset < data_.data.size() );
#endif
    data_.flags[idx.flagOffset] = sc.Flag();
    data_.data[idx.dataOffset]  = sc();
 }


/// Scalar & Array variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx ) const 
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR || idx.type == ARRAY );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
#endif
    return data_.flags[idx.flagOffset];
 }
 

/// Vector, Tensor, FlaggedArray variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx, uint32_t i ) const
 {
 #ifndef NDEBUG
  AssertPlacement(idx);
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    return data_.flags[ idx.flagOffset + i ];
 }
 
 
/// Scalar & Array variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx, VARIABLE_FLAG flag ) 
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == SCALAR  || idx.type == ARRAY );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
#endif
    data_.flags[idx.flagOffset] = flag;
 }

 
/// Vector,Tensor, FlaggedArray variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( const csmp::Index& idx, uint32_t i, VARIABLE_FLAG flag )
 {
#ifndef NDEBUG
 AssertPlacement(idx);
 assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    data_.flags[ idx.flagOffset + i ] = flag;
 }
 

/// Vector variable 
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const VectorVariable<dim>& vc )  
 {
#ifndef NDEBUG
 AssertPlacement( idx );
 assert( idx.type == VECTOR );
 assert( idx.index < data_.vectors );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim-1) < data_.data.size() );
#endif
    for ( uint32_t i{0u}; i<dim; ++i ) {
         data_.flags[ idx.flagOffset+i ] = vc.Flag(i);
         data_.data[ idx.dataOffset+i ]  = vc[i];
      }
 }


/// Vector variable 
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, VectorVariable<dim>& vc ) const  
 {
#ifndef NDEBUG
 AssertPlacement( idx );
 assert( idx.type == VECTOR );
 assert( idx.index < data_.vectors );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim-1) < data_.data.size() );
#endif
    for ( uint32_t i{0u}; i<dim; ++i ) {
         vc.Flag(i) = data_.flags[ idx.flagOffset+i ];
         vc(i)      = data_.data[ idx.dataOffset+i ];
      }
 }
 

/// Tensor variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const TensorVariable<dim>& ts )
 {
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == TENSOR );
  assert( idx.index < data_.tensors );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.flagOffset+dim-1) < data_.flags.size() );
  assert( (idx.dataOffset+dim*dim-1) < data_.data.size() );
#endif
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
#ifndef NDEBUG
 AssertPlacement( idx );
 assert( idx.type == TENSOR );
 assert( idx.index < data_.tensors );
#endif  
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim*dim-1) < data_.data.size() );
#endif
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
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == ARRAY );
  assert( av.Size() == idx.dataDepth );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset) < data_.flags.size() );
#endif
    const uint32_t data_offset( idx.dataOffset );
    const uint32_t flags_offset( idx.flagOffset );
    const uint32_t arraySize( idx.dataDepth );
    for( uint32_t i{0u}; i < arraySize; ++i )
      data_.data[ data_offset   + i ] = av[i];
    data_.flags[ flags_offset] = av.Flag();
  }


/// Array variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, ArrayVariable& av ) const
  {
    av.Resize( idx.dataDepth );
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == ARRAY );
  assert( av.Size() == idx.dataDepth );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset) < data_.flags.size() );
#endif
    const uint32_t data_offset( idx.dataOffset );
    const uint32_t flags_offset( idx.flagOffset );
    const uint32_t arraySize( idx.dataDepth );
    for( uint32_t i{0u}; i < arraySize; ++i )
      av(i)     = data_.data[ data_offset +  i ];
    av.Flag()= data_.flags[ flags_offset];
}

/// FlaggedArray variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( const csmp::Index& idx, const FlaggedArrayVariable& av )
  {
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == FLAGGEDARRAY );
  assert( av.Size() == idx.dataDepth ); 
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset+idx.dataDepth-1) < data_.flags.size() );
#endif
    const uint32_t data_offset( idx.dataOffset );
    const uint32_t flags_offset( idx.flagOffset );
    const uint32_t arraySize( idx.dataDepth );
    for( uint32_t i(0); i < arraySize; ++i )
    {
      data_.data[ data_offset   + i ] = av[i];
      data_.flags[ flags_offset + i ] = av.Flag(i);
    }
  }


/// FlaggedArray variable
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( const csmp::Index& idx, FlaggedArrayVariable& av ) const
  {
  av.Resize( idx.dataDepth );
#ifndef NDEBUG
  AssertPlacement( idx );
  assert( idx.type == FLAGGEDARRAY );
  assert( av.Size() == idx.dataDepth ); 
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset+idx.dataDepth-1) < data_.flags.size() );
#endif
    const uint32_t data_offset( idx.dataOffset );
    const uint32_t flags_offset( idx.flagOffset );
    const uint32_t arraySize( idx.dataDepth );
    for( uint32_t i{0u}; i < arraySize; ++i )
    {
      av(i)     = data_.data[ data_offset +  i ];
      av.Flag(i)= data_.flags[ flags_offset + i ];
    }
}


/**
    @todo replace this super wasteful method
*/
template<uint32_t dim, template<uint32_t> class STOREE>
bool LocalVariableStorage<dim,STOREE>::IsWithinRange( const csmp::Index& idx, 
                                                      double vmin, double vmax ) const
  {
#ifndef NDEBUG
  AssertPlacement( idx );
#endif
    if ( idx.type == SCALAR ) {
      const double val = Read( idx );
      return ( val >= vmin and val <= vmax ) ? true : false;
      }
    if ( idx.type == VECTOR ) {
      VectorVariable<dim>  vc;
      Read( idx, vc );
      return vc.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == TENSOR ) {
      TensorVariable<dim>  ts;
      Read( idx, ts );
      return ts.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == ARRAY ) {
      ArrayVariable  av;
      Read( idx, av );
      return av.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == FLAGGEDARRAY ) {
      FlaggedArrayVariable  fa;
      Read( idx, fa );
      return fa.IsWithinRange( vmin, vmax );
      }

    std::cerr <<"\nLocalVariableStorage<dim,STOREE>::IsWithinRange: range check could not be performed."<< std::endl;
    return false;
  }





// ===================================
// ELEMENT INTEGRATION POINT VARIABLES
// ===================================

#include "LocalVariableStorageIndexArithmetic.h"


/// Scalar variable value at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
double LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::Index& idx ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
#endif

    return data_.data[offset];
  }


/// Scalar variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::Index& idx, ScalarVariable& sc ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif

    sc.Flag() = data_.flags[flagOffset];
    sc        = data_.data[offset];
  }


/// Scalar variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::Index& idx, const ScalarVariable& sc )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif

    data_.flags[flagOffset] = sc.Flag();
    data_.data[offset] = sc();
  }


/// Scalar & Array variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::Index& idx ) const
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR || idx.type == ARRAY );
  assert( flagOffset < data_.flags.size() );
#endif

    return data_.flags[flagOffset];
  }


/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::Index& idx, uint32_t i ) const
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
    AssertIntegrationPointPlacement(idx);
    assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
    assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
    assert( flagOffset+i < data_.flags.size() );
#endif

    return data_.flags[ flagOffset+i ];
  }


/// Scalar & Array variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::Index& idx, VARIABLE_FLAG flag )
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == SCALAR  || idx.type == ARRAY );
  assert( flagOffset < data_.flags.size() );
#endif

    data_.flags[flagOffset] = flag;
  }


/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::Index& idx, uint32_t i, VARIABLE_FLAG flag )
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR || idx.type == TENSOR || idx.type == FLAGGEDARRAY );
  assert( (i < dim)&&(idx.type != FLAGGEDARRAY) || (i < idx.dataDepth )&&(idx.type == FLAGGEDARRAY) );
  assert( flagOffset+i < data_.flags.size() );
#endif

    data_.flags[ flagOffset+i ] = flag;
  }


/// Vector variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::Index& idx, const VectorVariable<dim>& vc )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif

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

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif

    for ( uint32_t i(0); i<dim; ++i ) {
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

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif

  for ( auto i{0U}; i<dim; i++ )
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

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( idx.type == TENSOR );
  assert( offset+dim*dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif

  for ( auto i{0U}; i<dim; i++ )
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

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif

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

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif

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

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif

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

#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif

    for( uint32_t i(0); i < arraySize; ++i )
    {
      av(i) = data_.data[ offset+i ];
      av.Flag( i, data_.flags[flagOffset+i] );
    }
  }


template<uint32_t dim, template<uint32_t> class STOREE>
bool LocalVariableStorage<dim,STOREE>::IsWithinRange( uint32_t ip, const csmp::Index& idx,
                                                      double vmin, double vmax ) const
  {
#ifndef NDEBUG
  AssertIntegrationPointPlacement(idx);
#endif

    if ( idx.type == SCALAR ) {
      const double val = Read( ip, idx );
      return ( val >= vmin and val <= vmax ) ? true : false;
      }
    if ( idx.type == VECTOR ) {
      VectorVariable<dim>  vc;
      Read( ip, idx, vc );
      return vc.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == TENSOR ) {
      TensorVariable<dim>  ts;
      Read( ip, idx, ts );
      return ts.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == ARRAY ) {
      ArrayVariable  av;
      Read( ip, idx, av );
      return av.IsWithinRange( vmin, vmax );
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx);
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
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
VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, uint32_t i ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
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
void LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, uint32_t i, VARIABLE_FLAG flag )
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const uint32_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( idx.type == VECTOR );
  assert( offset+dim-1 < data_.data.size() );
  assert( flagOffset+dim-1 < data_.flags.size() );
#endif
    const uint32_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

    const uint32_t arraySize( idx.dataDepth );

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif
    const uint32_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i(0); i < arraySize; ++i )
      data_.data[ offset + sector_ip_offset + i ]     = av[i];
    data_.flags[ flagOffset + sector_ip_flag_offset ] = av.Flag();
  }





/// Array variables at sector or facet integration points
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, ArrayVariable& av ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

    const uint32_t arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == ARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif
    const uint32_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
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
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

    const uint32_t arraySize( idx.dataDepth );

#ifndef NDEBUG  
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif
    const uint32_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i(0); i < arraySize; ++i ) {
         data_.data[ offset     + sector_ip_offset + i ]      = av[i];
         data_.flags[flagOffset + sector_ip_flag_offset + i ] = av.Flag(i);
      }
  }



/// FlaggedArray variables
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& idx, FlaggedArrayVariable& av ) const
  {
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

    const uint32_t arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
  lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip );
  AssertFiniteVolumeIntegrationPointPlacement(idx);
  assert( av.Size() == idx.dataDepth );
  assert( idx.type == FLAGGEDARRAY );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif
    const uint32_t sector_ip_flag_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                         (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth;

    for ( uint32_t i(0); i < arraySize; ++i ) {
         av(i) = data_.data[ offset + sector_ip_offset + i ];
         av.Flag( i, data_.flags[flagOffset + sector_ip_flag_offset + i] );
      }
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
   lvsCompileTimeDispatch::AssertFiniteVolumeIntegrationPointIndex( storeePtr, sector_or_facet, ip);
   AssertFiniteVolumeIntegrationPointPlacement(idx);
#endif

    if ( idx.type == SCALAR ) {
         const double val = Read( sector_or_facet, ip, idx );
         return ( val >= vmin and val <= vmax ) ? true : false;
      }
    if ( idx.type == VECTOR ) {
         VectorVariable<dim>  vc;
         Read( sector_or_facet, ip, idx, vc );
         return vc.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == TENSOR ) {
         TensorVariable<dim>  ts;
         Read( sector_or_facet, ip, idx, ts );
         return ts.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == ARRAY ) {
         ArrayVariable  av;
         Read( sector_or_facet, ip, idx, av );
         return av.IsWithinRange( vmin, vmax );
      }
    if ( idx.type == FLAGGEDARRAY ) {
         FlaggedArrayVariable  fa;
         Read( sector_or_facet, ip, idx, fa );
         return fa.IsWithinRange( vmin, vmax );
      }

    std::cout <<"\nLocalVariableStorage<dim,STOREE>::IsWithinRange: range check could not be performed."<< std::endl;
    return false;
  }
  
  
  
template<uint32_t dim, template<uint32_t> class STOREE>
void LocalVariableStorage<dim,STOREE>::StoreArrayEntry( const csmp::Index& idx, double value, uint32_t array_elmt )
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
double LocalVariableStorage<dim,STOREE>::ReadArrayEntry( const csmp::Index& idx, uint32_t array_elmt ) const
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
