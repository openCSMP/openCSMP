#ifndef LOCAL_VARIABLE_STORAGE_H
#define LOCAL_VARIABLE_STORAGE_H

#include "LocalVariableStorageIndexArithmetic.h"
#include "Index.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include "enumTypeCompatibilityChecks.h"
#include "compileTimeVariableDispatch.h"

//#define VARIABLE_STORAGE_DEBUG

namespace csmp {

/** CSMP local/physical variable storage

@author P. Lang
@author S.K. Matthai
@date 2007-2012

@attention CSMP runs should be performed at least once in DEBUG mode with conceptual models to check validity of Read/Store calls. For performance, these checks are omitted in RELEASE.

@section storees Storee Concept

Storees have to follow:

 - providing a @code PLACEMENT Placement() const { return PLACEMENT; } @endcode function, in which PLACEMENT is i.e. ELEMENT, NODE etc. (global csmp enum)
 - in the assignment operator (and copy ctor), add: @code this->LVS( el.LVS() ); @endcode

@section design Design

We use CRT static polymorphism here to allow for non-virtual Placement() functions and exploit compile time type information.
All methods specific to storees with integration points (Element, Face, InterFace) will fail AT COMPILE TIME when
used for storees without integration points (Region, Boundary, SplitBoundary)

@note Ideally we would like to use type traits instead of crt, but this would require change of csmp interface and wouldn't provide any performance increase

@section implementation Implementation

We store always in order: all ScalarVariable data/flags, then all VectorVariable data/flags, then all TensorVariable data/flags and finally all ArrayVariable data/flags
At integration points this holds, whith all ELEMENT_INTEGRATION_POINT variables first (storing all such values for the first integration point, then for the second and so forth), then all SECTOR_INTEGRATION_POINT and last all FACET_INTEGRATION_POINT (for each IP)
We perform validity checks on Read/Store etc in debug mode only, that is with no 'NDEBUG' symbol. Additional checks are activated if 'VARIABLE_STORAGE_DEBUG' is set.

@attention The variable storage is implicitly coupled with the PropertyDatabase through the Index and its offset calculation performed there

@todo (3-D) We could now even go for a single Read/Write function templatized on the Variable type. Variable placement would then be needed to be available in variables.
@note It's probably not sensible to merge duplicated functionality here for performance reasons
@todo (2-D) DocMe (index arithmetic)
*/
template<uint32_t dim, template<uint32_t> class STOREE>
class LocalVariableStorage {
  public:
    // ctors, dtor and assignment
    LocalVariableStorage();
    ~LocalVariableStorage();
    LocalVariableStorage( const LocalVariables& lv );
    LocalVariableStorage( const LocalVariables& lv, const IntegrationPointVariables& iv );
    LocalVariableStorage( const LocalVariableStorage& );
    LocalVariableStorage( LocalVariableStorage&& ) = default;
    LocalVariableStorage& operator=( const LocalVariableStorage& );
    LocalVariableStorage& operator=( LocalVariableStorage&& ) = default;

    // size ops
    void            ResizePropertyStorage   ( const LocalVariables& lv );
    void            ResizePropertyStorage   ( const LocalVariables& lv, const IntegrationPointVariables& iv );
    void            ResizePropertyStorage   ( uint32_t dataComponents, uint32_t flagComponents );
    void            AddProperty             ( const csmp::Index& );
    void            DeleteProperty          ( const csmp::Index& );

    // local variables access
    bool            IsWithinRange( const csmp::Index&, double, double ) const;
    double          Read    ( const csmp::Index& )                          const;
    void            Read    ( const csmp::Index&, ScalarVariable& )         const;
    void            Read    ( const csmp::Index&, VectorVariable<dim>& )    const;
    void            Read    ( const csmp::Index&, TensorVariable<dim>& )    const;
    void            Read    ( const csmp::Index&, ArrayVariable& )          const;
    void            Read    ( const csmp::Index&, FlaggedArrayVariable& )   const;
    void            Store   ( const csmp::Index&, const ScalarVariable& );
    void            Store   ( const csmp::Index&, const VectorVariable<dim>& );
    void            Store   ( const csmp::Index&, const TensorVariable<dim>& );
    void            Store   ( const csmp::Index&, const ArrayVariable& );
    void            Store   ( const csmp::Index&, const FlaggedArrayVariable& );
    VARIABLE_FLAG   Status  ( const csmp::Index& ) const;                   // scalars & arrays
    VARIABLE_FLAG   Status  ( const csmp::Index&, uint32_t ) const;           // vectors & tensors & flagged arrays
    void            Status  ( const csmp::Index&, VARIABLE_FLAG );          // scalars & arrays
    void            Status  ( const csmp::Index&, uint32_t, VARIABLE_FLAG );  // vectors & tensors & flagged arrays

    // integration point variables (will fail at COMPILE TIME when used for storees without integration points)
    bool            IsWithinRange( uint32_t ip, const csmp::Index&, double, double )  const;
    double          Read    ( uint32_t ip, const csmp::Index& )                           const;
    void            Read    ( uint32_t ip, const csmp::Index&, ScalarVariable& )          const;
    void            Read    ( uint32_t ip, const csmp::Index&, VectorVariable<dim>& )     const;
    void            Read    ( uint32_t ip, const csmp::Index&, TensorVariable<dim>& )     const;
    void            Read    ( uint32_t ip, const csmp::Index&, ArrayVariable& )           const;
    void            Read    ( uint32_t ip, const csmp::Index&, FlaggedArrayVariable& )    const;
    void            Store   ( uint32_t ip, const csmp::Index&, const ScalarVariable& );
    void            Store   ( uint32_t ip, const csmp::Index&, const VectorVariable<dim>& );
    void            Store   ( uint32_t ip, const csmp::Index&, const TensorVariable<dim>& );
    void            Store   ( uint32_t ip, const csmp::Index&, const ArrayVariable& );
    void            Store   ( uint32_t ip, const csmp::Index&, const FlaggedArrayVariable& );
    VARIABLE_FLAG   Status  ( uint32_t ip, const csmp::Index& ) const;                   // scalars & arrays
    VARIABLE_FLAG   Status  ( uint32_t ip, const csmp::Index&, uint32_t ) const;           // vectors & tensors  & flagged arrays
    void            Status  ( uint32_t ip, const csmp::Index&, VARIABLE_FLAG );          // scalars & arrays
    void            Status  ( uint32_t ip, const csmp::Index&, uint32_t, VARIABLE_FLAG );  // vectors & tensors  & flagged arrays

    // finite volume integration point variables (will fail at COMPILE TIME when used for storees without fv integration points)
    bool            IsWithinRange( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, double, double )  const;
    double          Read    ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& )                           const;
    void            Read    ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, ScalarVariable& )          const;
    void            Read    ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, VectorVariable<dim>& )     const;
    void            Read    ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, TensorVariable<dim>& )     const;
    void            Read    ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, ArrayVariable& )           const;
    void            Read    ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, FlaggedArrayVariable& )    const;
    void            Store   ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, const ScalarVariable&      );
    void            Store   ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, const VectorVariable<dim>& );
    void            Store   ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, const TensorVariable<dim>& );
    void            Store   ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, const ArrayVariable&       );
    void            Store   ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, const FlaggedArrayVariable&   );
    VARIABLE_FLAG   Status  ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index& ) const;                   // scalars & arrays
    VARIABLE_FLAG   Status  ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, uint32_t ) const;           // vectors & tensors  & flagged arrays
    void            Status  ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, VARIABLE_FLAG );          // scalars & arrays
    void            Status  ( uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, uint32_t, VARIABLE_FLAG );  // vectors & tensors  & flagged arrays

    // -----------------------------------------------------------------------------------------------------------------------
    // METHODS that make use of static dispatching via the template<VARIABLE_TYPE ty,PLACEMENT pl> struct INDEX : public Index
    // -----------------------------------------------------------------------------------------------------------------------
    // SKM 20/6/2020
    // nodes, elements, faces, interfaces
    template<PLACEMENT place> double Read( const csmp::INDEX<SCALAR,place>& ) const;
    template<PLACEMENT place> void     Read( const csmp::INDEX<SCALAR,place>&, ScalarVariable& ) const;
    template<PLACEMENT place> void     Read( const csmp::INDEX<VECTOR,place>&, VectorVariable<dim>& ) const;
    template<PLACEMENT place> void     Read( const csmp::INDEX<TENSOR,place>&, TensorVariable<dim>& ) const;
    template<PLACEMENT place> void     Read( const csmp::INDEX<ARRAY,place>&,  ArrayVariable& ) const;
    template<PLACEMENT place> void     Read( const csmp::INDEX<FLAGGEDARRAY,place>&, FlaggedArrayVariable& ) const;
    template<PLACEMENT place> void     Store( const csmp::INDEX<SCALAR,place>&, const ScalarVariable& );
    template<PLACEMENT place> void     Store( const csmp::INDEX<VECTOR,place>&, const VectorVariable<dim>& );
    template<PLACEMENT place> void     Store( const csmp::INDEX<TENSOR,place>&, const TensorVariable<dim>& );
    template<PLACEMENT place> void     Store( const csmp::INDEX<ARRAY,place>&, const ArrayVariable& );
    template<PLACEMENT place> void     Store( const csmp::INDEX<FLAGGEDARRAY,place>&, const FlaggedArrayVariable& );
    template<PLACEMENT place> VARIABLE_FLAG  Status( const csmp::INDEX<SCALAR,place>& ) const;                   // scalars & arrays
    template<PLACEMENT place> VARIABLE_FLAG  Status( const csmp::INDEX<ARRAY,place>& ) const;                    // scalars & arrays
    template<PLACEMENT place> VARIABLE_FLAG  Status( const csmp::INDEX<VECTOR,place>&, uint32_t ) const;           // vectors & tensors & flagged arrays
    template<PLACEMENT place> VARIABLE_FLAG  Status( const csmp::INDEX<TENSOR,place>&, uint32_t ) const;           // vectors & tensors & flagged arrays
    template<PLACEMENT place> VARIABLE_FLAG  Status( const csmp::INDEX<FLAGGEDARRAY,place>&, uint32_t ) const;           // vectors & tensors & flagged arrays
    template<PLACEMENT place> void     Status( const csmp::INDEX<SCALAR,place>&, VARIABLE_FLAG );          // scalars & arrays
    template<PLACEMENT place> void     Status( const csmp::INDEX<ARRAY,place>&, VARIABLE_FLAG );          // scalars & arrays
    template<PLACEMENT place> void     Status( const csmp::INDEX<VECTOR,place>&, uint32_t, VARIABLE_FLAG );  // vectors & tensors & flagged arrays
    template<PLACEMENT place> void     Status( const csmp::INDEX<TENSOR,place>&, uint32_t, VARIABLE_FLAG );  // vectors & tensors & flagged arrays
    template<PLACEMENT place> void     Status( const csmp::INDEX<FLAGGEDARRAY,place>&, uint32_t, VARIABLE_FLAG );  // vectors & tensors & flagged arrays
    
    // TODO: create some methods that allow to read individual entries in ARRAY or FLAGGEDARRAY variable rather than the whole array
    // TODO: IsWithinRange() methods have not been adapted to INDEX yet

    // integration point variables (will fail at COMPILE TIME when used for storees without integration points)
    template<PLACEMENT place> double        Read  ( uint32_t ip, const csmp::INDEX<SCALAR,place>& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t ip, const csmp::INDEX<SCALAR,place>&, ScalarVariable& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t ip, const csmp::INDEX<VECTOR,place>&, VectorVariable<dim>& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t ip, const csmp::INDEX<TENSOR,place>&, TensorVariable<dim>& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t ip, const csmp::INDEX<ARRAY,place>&,  ArrayVariable& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>&, FlaggedArrayVariable& ) const;
    template<PLACEMENT place> void            Store ( uint32_t ip, const csmp::INDEX<SCALAR,place>&, const ScalarVariable& );
    template<PLACEMENT place> void            Store ( uint32_t ip, const csmp::INDEX<VECTOR,place>&, const VectorVariable<dim>& );
    template<PLACEMENT place> void            Store ( uint32_t ip, const csmp::INDEX<TENSOR,place>&, const TensorVariable<dim>& );
    template<PLACEMENT place> void            Store ( uint32_t ip, const csmp::INDEX<ARRAY,place>&,  const ArrayVariable& );
    template<PLACEMENT place> void            Store ( uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>&, const FlaggedArrayVariable& );
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t ip, const csmp::INDEX<SCALAR,place>& ) const;  // scalars & arrays
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t ip, const csmp::INDEX<ARRAY,place>& ) const;
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t ip, const csmp::INDEX<VECTOR,place>&, uint32_t ) const; // vectors & tensors & flagged arrays
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t ip, const csmp::INDEX<TENSOR,place>&, uint32_t ) const;
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>&, uint32_t ) const;
    template<PLACEMENT place> void            Status( uint32_t ip, const csmp::INDEX<SCALAR,place>&, VARIABLE_FLAG ); // scalars & arrays
    template<PLACEMENT place> void            Status( uint32_t ip, const csmp::INDEX<ARRAY,place>&, VARIABLE_FLAG );
    template<PLACEMENT place> void            Status( uint32_t ip, const csmp::INDEX<VECTOR,place>&, uint32_t, VARIABLE_FLAG ); // vectors & tensors & flagged arrays
    template<PLACEMENT place> void            Status( uint32_t ip, const csmp::INDEX<TENSOR,place>&, uint32_t, VARIABLE_FLAG );
    template<PLACEMENT place> void            Status( uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>&, uint32_t, VARIABLE_FLAG );

    // finite volume integration point variables (will fail at COMPILE TIME when used for storees without fv integration points)
    template<PLACEMENT place> double        Read  ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>&, ScalarVariable& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<VECTOR,place>&, VectorVariable<dim>& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<TENSOR,place>&, TensorVariable<dim>& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<ARRAY,place>&,  ArrayVariable& ) const;
    template<PLACEMENT place> void            Read  ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>&, FlaggedArrayVariable& ) const;
    template<PLACEMENT place> void            Store ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>&, const ScalarVariable& );
    template<PLACEMENT place> void            Store ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<VECTOR,place>&, const VectorVariable<dim>& );
    template<PLACEMENT place> void            Store ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<TENSOR,place>&, const TensorVariable<dim>& );
    template<PLACEMENT place> void            Store ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<ARRAY,place>&,  const ArrayVariable& );
    template<PLACEMENT place> void            Store ( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>&, const FlaggedArrayVariable& );
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>& ) const; // scalars & arrays
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<ARRAY,place>& ) const;
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<VECTOR,place>&, uint32_t ) const; // vectors & tensors & flagged arrays
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<TENSOR,place>&, uint32_t ) const;
    template<PLACEMENT place> VARIABLE_FLAG   Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>&, uint32_t ) const;
    template<PLACEMENT place> void            Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>&, VARIABLE_FLAG ); // scalars & arrays
    template<PLACEMENT place> void            Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<ARRAY,place>&, VARIABLE_FLAG );
    template<PLACEMENT place> void            Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<VECTOR,place>&, uint32_t, VARIABLE_FLAG ); // vectors & tensors & flagged arrays
    template<PLACEMENT place> void            Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<TENSOR,place>&, uint32_t, VARIABLE_FLAG );
    template<PLACEMENT place> void            Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>&, uint32_t, VARIABLE_FLAG );

    // output
    bool EmptyLVS() const;
    void OutLVS() const;

#ifndef NDEBUG
    void AssertPlacement( const csmp::Index& ) const;
    void AssertIntegrationPointPlacement( const csmp::Index& ) const;
    void AssertFiniteVolumeIntegrationPointPlacement( const csmp::Index& ) const;
    void StoreLocalState( const LocalVariables& lv );
#endif

  public:
    struct Data {
        typedef std::vector<VARIABLE_FLAG> FlagContainer; 
        typedef std::vector<double>        DataContainer;
        FlagContainer flags;
        DataContainer data;
#ifdef NDEBUG
        Data() : flags(0U), data(0U) {}
        Data( const Data& d ) : flags( d.flags ), data( d.data ) {}
        Data( Data&& d ) : flags( std::move(d.flags) ), data( std::move(d.data) ) {}
        Data& operator=( const Data& d )
          { if ( &d != this ) 
             { flags = d.flags; data = d.data; }
          return *this; }
        Data& operator=( Data&& ) = default;
#else // debugging: a lot more information is kept in storage
        uint32_t  scalars,            ///< scalar variables stored at the site this policy is associated with
                vectors,            ///< vector variables at this site
                tensors,            ///< tensor variables at this site
                arrays,             ///< array variables at this site
                flaggedArrays;      ///< flagged array variables at this site
        
        uint32_t  arrayLength,        ///< length of array variables associated with this site @todo only one size?
                flaggedArrayLength; ///< length of flagged array variables @todo only one size?

        Data() :
            flags              (0U),
            data               (0U),
            scalars            (0U),
            vectors            (0U),
            tensors            (0U),
            arrays             (0U),
            flaggedArrays      (0U),
            arrayLength        (0U),
            flaggedArrayLength (0U)
        {}

        Data( const Data& d ) :
            flags               ( d.flags ),
            data                ( d.data ),
            scalars             ( d.scalars ),
            vectors             ( d.vectors ),
            tensors             ( d.tensors ),
            arrays              ( d.arrays ),
            flaggedArrays       ( d.flaggedArrays ),
            arrayLength         ( d.arrayLength ),
            flaggedArrayLength  ( d.flaggedArrayLength )
        {}

        Data& operator=( const Data& d )
          {
            if ( &d != this ) {
                flags               = d.flags;
                data                = d.data;
                scalars             = d.scalars;
                vectors             = d.vectors;
                tensors             = d.tensors;
                arrays              = d.arrays;
                flaggedArrays       = d.flaggedArrays;
                arrayLength         = d.arrayLength;
                flaggedArrayLength  = d.flaggedArrayLength;
             }
            return *this;
         }
#endif /* debugging of Data */
      };

    /// copy or move supplied data into storage
    void LVS( const Data& data ) { data_ = data; }
    void LVS( Data&& data ) { data_ = std::move(data); }
    
    /// return copy of the data container for misc operations
    const Data LVS() const { return data_; }
  
  private:
    Data data_; ///< data and flag containers
};


// -----------------------------------------------------------------------------------------------------------------------
// METHODS that make use of static dispatching via the template<VARIABLE_TYPE ty,PLACEMENT pl> struct INDEX : public Index
// -----------------------------------------------------------------------------------------------------------------------
// SKM 20/6/2020
// nodes, elements, faces, interfaces
  
/**
double        Read    ( const csmp::INDEX<SCALAR,NODE>& ) const;
*/
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline double LocalVariableStorage<dim,STOREE>::Read( const csmp::INDEX<SCALAR,place>& idx ) const  
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
#ifndef NDEBUG
 assert( idx.index < data_.scalars );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.dataOffset < data_.data.size() );
#endif
    return data_.data[idx.dataOffset];
 }
 
/**

void            Read    ( const csmp::INDEX<SCALAR,NODE>&, ScalarVariable& )         const;

*/
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( const csmp::INDEX<SCALAR,place>& idx, ScalarVariable& sc ) const  
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( const csmp::INDEX<SCALAR,place>& idx, const ScalarVariable& sc )  
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
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
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<SCALAR,place>& idx ) const 
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
#endif
    return data_.flags[idx.flagOffset];
 }

/// ArrayVariable 
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<ARRAY,place>& idx ) const 
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
#endif
    return data_.flags[idx.flagOffset];
 }




/// Vector, Tensor, FlaggedArray variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<VECTOR,place>& idx, uint32_t i ) const
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
  assert( i < dim );
  assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    return data_.flags[ idx.flagOffset + i ];
 }

// tensor only has flags for its diagnao elements
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<TENSOR,place>& idx, uint32_t i ) const
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
  assert( i < dim );
  assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    return data_.flags[ idx.flagOffset + i ];
 }


template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<FLAGGEDARRAY,place>& idx, uint32_t i ) const
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
  const uint32_t array_length(idx.flagDepth);
  assert( i < array_length );
  assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    return data_.flags[ idx.flagOffset + i ];
 }
 

 
 
/// Scalar & Array variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<SCALAR,place>& idx, VARIABLE_FLAG flag ) 
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
#endif
    data_.flags[idx.flagOffset] = flag;
 }


template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<ARRAY,place>& idx, VARIABLE_FLAG flag ) 
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
 assert( idx.flagOffset < data_.flags.size() );
#endif
    data_.flags[idx.flagOffset] = flag;
 }


 
/// Vector variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<VECTOR,place>& idx, uint32_t i, VARIABLE_FLAG flag )
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
 assert( i < dim );
 assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    data_.flags[ idx.flagOffset + i ] = flag;
 }
 

/// Tensor variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<TENSOR,place>& idx, uint32_t i, VARIABLE_FLAG flag )
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
 assert( i < dim );
 assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    data_.flags[ idx.flagOffset + i ] = flag;
 }


/// FlaggedArray variable flag
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( const csmp::INDEX<FLAGGEDARRAY,place>& idx, uint32_t i, VARIABLE_FLAG flag )
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifdef VARIABLE_STORAGE_DEBUG
  const uint32_t array_length(idx.flagDepth);
  assert( i < array_length );
  assert( (idx.flagOffset+i) < data_.flags.size() );
#endif
    data_.flags[ idx.flagOffset + i ] = flag;
 }






/// Vector variable 
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( const csmp::INDEX<VECTOR,place>& idx, const VectorVariable<dim>& vc )  
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
 assert( idx.index < data_.vectors );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim-1) < data_.data.size() );
#endif
    for ( uint32_t i(0); i<dim; ++i ) {
         data_.flags[ idx.flagOffset+i ] = vc.Flag(i);
         data_.data[ idx.dataOffset+i ]  = vc[i];
      }
 }


/// Vector variable 
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( const csmp::INDEX<VECTOR,place>& idx, VectorVariable<dim>& vc ) const  
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
 assert( idx.index < data_.vectors );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim-1) < data_.data.size() );
#endif
    for ( uint32_t i(0); i<dim; ++i ) {
         vc.Flag(i) = data_.flags[ idx.flagOffset+i ];
         vc(i)      = data_.data[ idx.dataOffset+i ];
      }
 }
 

/// Tensor variable
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( const csmp::INDEX<TENSOR,place>& idx, const TensorVariable<dim>& ts )
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( idx.index < data_.tensors );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.flagOffset+dim-1) < data_.flags.size() );
  assert( (idx.dataOffset+dim*dim-1) < data_.data.size() );
#endif
    const uint32_t dataOffset(idx.dataOffset);
    const uint32_t flagOffset(idx.flagOffset);
    for ( auto i{0U}; i<dim; i++ )
      {
        data_.flags[ flagOffset+i ] = ts.Flag(i);
        for ( uint32_t j{0U}; j<dim; j++ )
          data_.data[ dataOffset+i*dim+j ] = ts(i,j);
      }
 }


/// Tensor variable
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( const csmp::INDEX<TENSOR,place>& idx, TensorVariable<dim>& ts ) const  
 {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
 assert( idx.index < data_.tensors );
#endif  
#ifdef VARIABLE_STORAGE_DEBUG
 assert( (idx.flagOffset+dim-1) < data_.flags.size() );
 assert( (idx.dataOffset+dim*dim-1) < data_.data.size() );
#endif
   const uint32_t dataOffset(idx.dataOffset);
   const uint32_t flagOffset(idx.flagOffset);
   for ( auto i{0U}; i<dim; i++ )
     {
       ts.Flag(i) = data_.flags[ flagOffset+i ] ;
       for ( uint32_t j{0U}; j<dim; j++ )
         ts(i,j) = data_.data[ dataOffset+i*dim+j ];
     }
 }


/// Array variable
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( const csmp::INDEX<ARRAY,place>& idx, const ArrayVariable& av )
  {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
  assert( av.Size() == idx.dataDepth );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset) < data_.flags.size() );
#endif
    const uint32_t data_offset( idx.dataOffset );
    const uint32_t flags_offset( idx.flagOffset );
    const uint32_t arraySize( idx.dataDepth );
    for( uint32_t i(0); i < arraySize; ++i )
      data_.data[ data_offset   + i ] = av[i];
    data_.flags[ flags_offset] = av.Flag();
  }


/// Array variable
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( const csmp::INDEX<ARRAY,place>& idx, ArrayVariable& av ) const
  {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    av.Resize( idx.dataDepth );
#ifndef NDEBUG
  assert( av.Size() == idx.dataDepth );
#endif
#ifdef VARIABLE_STORAGE_DEBUG
  assert( (idx.dataOffset+idx.dataDepth-1) < data_.data.size() );
  assert( (idx.flagOffset) < data_.flags.size() );
#endif
    const uint32_t data_offset( idx.dataOffset );
    const uint32_t flags_offset( idx.flagOffset );
    const uint32_t arraySize( idx.dataDepth );
    for( uint32_t i(0); i < arraySize; ++i )
      av(i)   = data_.data[ data_offset +  i ];
    av.Flag() = data_.flags[ flags_offset];
}

/// FlaggedArray variable
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( const csmp::INDEX<FLAGGEDARRAY,place>& idx, const FlaggedArrayVariable& av )
  {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( const csmp::INDEX<FLAGGEDARRAY,place>& idx, FlaggedArrayVariable& av ) const
  {
    static_assert( TypeMatchesVariablePlacement<STOREE,place>::value, "LocalVariableStorage: STOREE type does not match PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
  av.Resize( idx.dataDepth );
#ifndef NDEBUG
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
      av(i)     = data_.data[ data_offset +  i ];
      av.Flag(i)= data_.flags[ flags_offset + i ];
    }
}
  

// --------------------------------------------------------------------------------------------------------  
// integration point variables (will fail at COMPILE TIME when used for storees without integration points)
// --------------------------------------------------------------------------------------------------------  

/// Scalar variable value at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline double LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::INDEX<SCALAR,place>& idx ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( offset < data_.data.size() );
#endif
    return data_.data[offset];
  }


/// Scalar variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::INDEX<SCALAR,place>& idx, ScalarVariable& sc ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif
    sc.Flag() = data_.flags[flagOffset];
    sc        = data_.data[offset];
  }


/// Scalar variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::INDEX<SCALAR,place>& idx, const ScalarVariable& sc )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( offset < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif
    data_.flags[flagOffset] = sc.Flag();
    data_.data[offset] = sc();
  }


/// Scalar & Array variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<SCALAR,place>& idx ) const
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( flagOffset < data_.flags.size() );
#endif
    return data_.flags[flagOffset];
  }


/// Scalar & Array variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<ARRAY,place>& idx ) const
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( flagOffset < data_.flags.size() );
#endif
    return data_.flags[flagOffset];
  }



/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<VECTOR,place>& idx, uint32_t i ) const
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
    assert( ip < IPS_SI );
    assert( i < dim );
    assert( flagOffset+i < data_.flags.size() );
#endif

    return data_.flags[ flagOffset+i ];
  }
 
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<TENSOR,place>& idx, uint32_t i ) const
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
    assert( ip < IPS_SI );
    assert( i < dim );
    assert( flagOffset+i < data_.flags.size() );
#endif

    return data_.flags[ flagOffset+i ];
  }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>& idx, uint32_t i ) const
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
    assert( ip < IPS_SI );
    const uint32_t array_length(idx.flagDepth);
    assert( i < array_length );
    assert( flagOffset+i < data_.flags.size() );
#endif

    return data_.flags[ flagOffset+i ];
  }
   
  


/// Scalar & Array variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<SCALAR,place>& idx, VARIABLE_FLAG flag )
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );

#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( flagOffset < data_.flags.size() );
#endif

    data_.flags[flagOffset] = flag;
  }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<ARRAY,place>& idx, VARIABLE_FLAG flag )
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( flagOffset < data_.flags.size() );
#endif

    data_.flags[flagOffset] = flag;
  }


/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<VECTOR,place>& idx, uint32_t i, VARIABLE_FLAG flag )
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( i < dim );
  assert( flagOffset+i < data_.flags.size() );
#endif
    data_.flags[ flagOffset+i ] = flag;
  }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<TENSOR,place>& idx, uint32_t i, VARIABLE_FLAG flag )
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( i < dim );
  assert( flagOffset+i < data_.flags.size() );
#endif
    data_.flags[ flagOffset+i ] = flag;
  }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>& idx, uint32_t i, VARIABLE_FLAG flag )
  {
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
  const uint32_t array_length(idx.flagDepth);
  assert( i < array_length );
  assert( flagOffset+i < data_.flags.size() );
#endif
    data_.flags[ flagOffset+i ] = flag;
  }



/// Vector variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::INDEX<VECTOR,place>& idx, const VectorVariable<dim>& vc )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::INDEX<VECTOR,place>& idx, VectorVariable<dim>& vc ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::INDEX<TENSOR,place>& idx, const TensorVariable<dim>& ts )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::INDEX<TENSOR,place>& idx, TensorVariable<dim>& ts ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::INDEX<ARRAY,place>& idx, const ArrayVariable& av )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);
    const uint32_t arraySize( idx.dataDepth );

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( av.Size() == idx.dataDepth );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset< data_.flags.size() );
#endif
    for( uint32_t i(0); i < arraySize; ++i )
      data_.data[ offset + i ] = av[i];
    data_.flags[flagOffset] = av.Flag();
  }


/// Array variables
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::INDEX<ARRAY,place>& idx, ArrayVariable& av ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);
    const uint32_t arraySize( idx.dataDepth );

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( av.Size() == idx.dataDepth );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset < data_.flags.size() );
#endif
    for( uint32_t i(0); i < arraySize; ++i )
      av(i) = data_.data[ offset+i ];

    av.Flag( data_.flags[flagOffset] );
  }


/// FlaggedArray variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>& idx, const FlaggedArrayVariable& av )
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);
    const uint32_t arraySize( idx.dataDepth );

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( av.Size() == idx.dataDepth );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>& idx, FlaggedArrayVariable& av ) const
  {
    const uint32_t offset(DATA_OFFSET_IP);
    const uint32_t flagOffset(FLAG_OFFSET_IP);
    const uint32_t arraySize( idx.dataDepth );

    static_assert( place == ELEMENT_INTEGRATION_POINT || place == FACE_INTEGRATION_POINT || place == INTER_FACE_INTEGRATION_POINT, "LocalVariableStorage: STOREE type does not match IP PLACEMENT enumeration" );
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
#ifndef NDEBUG
  assert( ip < IPS_SI );
  assert( av.Size() == idx.dataDepth );
  assert( offset+arraySize-1 < data_.data.size() );
  assert( flagOffset+arraySize-1 < data_.flags.size() );
#endif
    for( uint32_t i(0); i < arraySize; ++i )
      {
        av(i) = data_.data[ offset+i ];
        av.Flag( i, data_.flags[flagOffset+i] );
      }
  }





// -------------------------------------------------------------------------------------------------------------------------  
// finite volume integration point variables (will fail at COMPILE TIME when used for storees without fv integration points)
// -------------------------------------------------------------------------------------------------------------------------  

/** 
    Scalar variable value at sector or facet integration points only.
*/
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline double LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>& idx ) const
  {
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    // offset to first instance of idx variable in the data vector
    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( offsetData.first < data_.data.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalDataDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalDataDepth; // for facet integration points only
    ///                                                                                    ^^^^^^^^
    return data_.data[ offsetData.first + sector_ip_offset ];
  }



/// Scalar variable at facet or sector integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>& idx, ScalarVariable& sc ) const
 {
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx);
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>& idx, const ScalarVariable& sc )
  {
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
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
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>& idx ) const
  {
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( offsetData.second < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;
    
    return data_.flags[offsetData.second + sector_ip_offset];
  }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<ARRAY,place>& idx ) const
  {
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
  assert( offsetData.second < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;
    
    return data_.flags[offsetData.second + sector_ip_offset];
  }



/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<VECTOR,place>& idx, uint32_t i ) const
  {
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( i < dim );
    assert( flagOffset+i < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;
    
    return data_.flags[ flagOffset + sector_ip_offset + i ];
 }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<TENSOR,place>& idx, uint32_t i ) const
  {
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( i < dim );
    assert( flagOffset+i < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;
    
    return data_.flags[ flagOffset + sector_ip_offset + i ];
 }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline VARIABLE_FLAG LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>& idx, uint32_t i ) const
  {
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    const uint32_t arraySize( idx.dataDepth );
    assert( i < arraySize );
    assert( flagOffset+i < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;
    
    return data_.flags[ flagOffset + sector_ip_offset + i ];
  }




/// Scalar & Array variable flag at facet or sector integration point (note that the Array has only a single flag)
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<SCALAR,place>& idx, VARIABLE_FLAG flag )
  {
    static_assert( TypeMatchesVariableType<ScalarVariable,SCALAR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( flagOffset < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    data_.flags[flagOffset + sector_ip_offset] = flag;
 }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<ARRAY,place>& idx, VARIABLE_FLAG flag )
  {
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t flagOffset(offsetData.second);
    
#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
  assert( flagOffset < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    data_.flags[flagOffset + sector_ip_offset] = flag;
 }



/// Vector, Tensor, FlaggedArray variable flag at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<VECTOR,place>& idx, uint32_t i, VARIABLE_FLAG flag )
  {
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( i < dim );
    assert( flagOffset+i < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    data_.flags[ flagOffset + sector_ip_offset + i ] = flag;
  }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<TENSOR,place>& idx, uint32_t i, VARIABLE_FLAG flag )
  {
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( i < dim );
    assert( flagOffset+i < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    data_.flags[ flagOffset + sector_ip_offset + i ] = flag;
  }

template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Status( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>& idx, uint32_t i, VARIABLE_FLAG flag )
  {
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    const uint32_t arraySize( idx.dataDepth );
    assert( i < arraySize );
    assert( flagOffset+i < data_.flags.size() );
#endif
    const uint32_t sector_ip_offset = (idx.place == SECTOR_INTEGRATION_POINT) ?
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvSector.totalFlagDepth :
                                    (sector_or_facet + ip) * idx.integrationPointVariables.ipvFacet.totalFlagDepth;

    data_.flags[ flagOffset + sector_ip_offset + i ] = flag;
  }



/// Vector variable at integration point
template<uint32_t dim, template<uint32_t> class STOREE>
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<VECTOR,place>& idx, const VectorVariable<dim>& vc )
  {
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<VECTOR,place>& idx, VectorVariable<dim>& vc ) const
  {
    static_assert( TemplateTypeMatchesVariableType<VectorVariable,VECTOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<TENSOR,place>& idx, const TensorVariable<dim>& ts )
  {
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<TENSOR,place>& idx, TensorVariable<dim>& ts ) const
  {
    static_assert( TemplateTypeMatchesVariableType<TensorVariable,TENSOR>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<ARRAY,place>& idx, const ArrayVariable& av )
  {
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);
    const uint32_t arraySize( idx.dataDepth );

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( av.Size() == idx.dataDepth );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<ARRAY,place>& idx, ArrayVariable& av ) const
  {
    static_assert( TypeMatchesVariableType<ArrayVariable,ARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);
    const uint32_t arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( av.Size() == idx.dataDepth );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Store( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>& idx, const FlaggedArrayVariable& av )
  {
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);
    const uint32_t arraySize( idx.dataDepth );

#ifndef NDEBUG  
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( av.Size() == idx.dataDepth );
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
template<PLACEMENT place> 
inline void LocalVariableStorage<dim,STOREE>::Read( uint32_t sector_or_facet, uint32_t ip, const csmp::INDEX<FLAGGEDARRAY,place>& idx, FlaggedArrayVariable& av ) const
  {
    static_assert( TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY>::value, "LocalVariableStorage: variable type does not match VARIABLE_TYPE enumeration" );
    static_assert( place == SECTOR_INTEGRATION_POINT || place == FACET_INTEGRATION_POINT ||
                   place == FACE_SECTOR_INTEGRATION_POINT || place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
                   place == FACE_FACET_INTEGRATION_POINT || idx.place == INTER_FACE_FACET_INTEGRATION_POINT,
                  "LocalVariableStorage: STOREE type does not match VARIABLE placement" );

    const STOREE<dim>* const storeePtr = static_cast<const STOREE<dim>*>(this);
    const std::pair<uint32_t, uint32_t> offsetData = lvsCompileTimeDispatch::containerOffset( storeePtr, idx );
    const uint32_t offset(offsetData.first);
    const uint32_t flagOffset(offsetData.second);
    const uint32_t arraySize( idx.dataDepth );
    av.Resize( idx.dataDepth );

#ifndef NDEBUG
    if ( place == SECTOR_INTEGRATION_POINT ) assert( ip < IPS_SE );
    if ( place == FACET_INTEGRATION_POINT ) assert( ip < IPS_FA );
    assert( sector_or_facet < storeePtr->Facets() || sector_or_facet < storeePtr->Sectors() );
    assert( ip < storeePtr->IntegrationPointsPerSector() || ip < storeePtr->IntegrationPointsPerFacet() );
    assert( av.Size() == idx.dataDepth );
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






} // end csmp

#endif
