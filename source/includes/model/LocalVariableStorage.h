#ifndef LOCAL_VARIABLE_STORAGE_H
#define LOCAL_VARIABLE_STORAGE_H

#include "Index.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"


#include <vector>

namespace csmp {

//template<size_t> class Element;
//template<size_t> class Face;
//template<size_t> class InterFace;

// CRT calls for elements/faces/interfaces
#define FE_SES static_cast<const STOREE*>(this)->Sectors()                                  ///< sectors of simplex storee
#define FE_FAS static_cast<const STOREE*>(this)->Facets()                                   ///< facets of simplex storee
#define FE_FVIPS_PER_SECTOR static_cast<const STOREE*>(this)->IntegrationPointsPerSector()  ///< integration points per sector
#define FE_FVIPS_PER_FACET static_cast<const STOREE*>(this)->IntegrationPointsPerFacet()    ///< integration points per facet
#define IPS_SI static_cast<const STOREE*>(this)->IntegrationPoints()                        ///< simplex integration points
#define IPS_SE FE_SES*FE_FVIPS_PER_SECTOR                                                   ///< sector integration points
#define IPS_FA FE_FAS*FE_FVIPS_PER_FACET                                                    ///< facet integration points


/** CSMP local/physical variable storage

@author P. Lang
@author S.K. Matthaei
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
template<size_t dim,class STOREE>
class LocalVariableStorage {
  public:
    // ctors, dtor and assignment
    LocalVariableStorage();
    ~LocalVariableStorage();
    LocalVariableStorage( const LocalVariables& lv );
    LocalVariableStorage( const LocalVariables& lv, const IntegrationPointVariables& iv );
    LocalVariableStorage( const LocalVariableStorage& );
    LocalVariableStorage( LocalVariableStorage&& );
    LocalVariableStorage& operator=( const LocalVariableStorage& );

    // size ops
    void            ResizePropertyStorage   ( const LocalVariables& lv );
    void            ResizePropertyStorage   ( const LocalVariables& lv, const IntegrationPointVariables& iv );
    void            ResizePropertyStorage   ( size_t dataComponents, size_t flagComponents );
    void            AddProperty             ( const csmp::Index& );
    void            DeleteProperty          ( const csmp::Index& );

    // local variables access
    bool            IsWithinRange( const csmp::Index&, double64, double64 ) const;
    double64        Read    ( const csmp::Index& )                          const;
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
    VARIABLE_FLAG   Status  ( const csmp::Index&, size_t ) const;           // vectors & tensors & flagged arrays
    void            Status  ( const csmp::Index&, VARIABLE_FLAG );          // scalars & arrays
    void            Status  ( const csmp::Index&, size_t, VARIABLE_FLAG );  // vectors & tensors & flagged arrays

    // integration point variables (will fail at COMPILE TIME when used for storees without integration points)
    bool            IsWithinRange( size_t ip, const csmp::Index&, double64, double64 )  const;
    double64        Read    ( size_t ip, const csmp::Index& )                           const;
    void            Read    ( size_t ip, const csmp::Index&, ScalarVariable& )          const;
    void            Read    ( size_t ip, const csmp::Index&, VectorVariable<dim>& )     const;
    void            Read    ( size_t ip, const csmp::Index&, TensorVariable<dim>& )     const;
    void            Read    ( size_t ip, const csmp::Index&, ArrayVariable& )           const;
    void            Read    ( size_t ip, const csmp::Index&, FlaggedArrayVariable& )    const;
    void            Store   ( size_t ip, const csmp::Index&, const ScalarVariable& );
    void            Store   ( size_t ip, const csmp::Index&, const VectorVariable<dim>& );
    void            Store   ( size_t ip, const csmp::Index&, const TensorVariable<dim>& );
    void            Store   ( size_t ip, const csmp::Index&, const ArrayVariable& );
    void            Store   ( size_t ip, const csmp::Index&, const FlaggedArrayVariable& );
    VARIABLE_FLAG   Status  ( size_t ip, const csmp::Index& ) const;                   // scalars & arrays
    VARIABLE_FLAG   Status  ( size_t ip, const csmp::Index&, size_t ) const;           // vectors & tensors  & flagged arrays
    void            Status  ( size_t ip, const csmp::Index&, VARIABLE_FLAG );          // scalars & arrays
    void            Status  ( size_t ip, const csmp::Index&, size_t, VARIABLE_FLAG );  // vectors & tensors  & flagged arrays

    // finite volume integration point variables (will fail at COMPILE TIME when used for storees without fv integration points)
    bool            IsWithinRange( size_t sector_or_facet, size_t ip, const csmp::Index&, double64, double64 )  const;
    double64        Read    ( size_t sector_or_facet, size_t ip, const csmp::Index& )                           const;
    void            Read    ( size_t sector_or_facet, size_t ip, const csmp::Index&, ScalarVariable& )          const;
    void            Read    ( size_t sector_or_facet, size_t ip, const csmp::Index&, VectorVariable<dim>& )     const;
    void            Read    ( size_t sector_or_facet, size_t ip, const csmp::Index&, TensorVariable<dim>& )     const;
    void            Read    ( size_t sector_or_facet, size_t ip, const csmp::Index&, ArrayVariable& )           const;
    void            Read    ( size_t sector_or_facet, size_t ip, const csmp::Index&, FlaggedArrayVariable& )    const;
    void            Store   ( size_t sector_or_facet, size_t ip, const csmp::Index&, const ScalarVariable&      );
    void            Store   ( size_t sector_or_facet, size_t ip, const csmp::Index&, const VectorVariable<dim>& );
    void            Store   ( size_t sector_or_facet, size_t ip, const csmp::Index&, const TensorVariable<dim>& );
    void            Store   ( size_t sector_or_facet, size_t ip, const csmp::Index&, const ArrayVariable&       );
    void            Store   ( size_t sector_or_facet, size_t ip, const csmp::Index&, const FlaggedArrayVariable&   );
    VARIABLE_FLAG   Status  ( size_t sector_or_facet, size_t ip, const csmp::Index& ) const;                   // scalars & arrays
    VARIABLE_FLAG   Status  ( size_t sector_or_facet, size_t ip, const csmp::Index&, size_t ) const;           // vectors & tensors  & flagged arrays
    void            Status  ( size_t sector_or_facet, size_t ip, const csmp::Index&, VARIABLE_FLAG );          // scalars & arrays
    void            Status  ( size_t sector_or_facet, size_t ip, const csmp::Index&, size_t, VARIABLE_FLAG );  // vectors & tensors  & flagged arrays

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
    struct Data
      {
        typedef std::vector<VARIABLE_FLAG> FlagContainer; 
        typedef std::vector<double64>      DataContainer;

        FlagContainer flags;
        DataContainer data;
#ifdef NDEBUG
        Data() : flags(0U), data(0U) {}
        Data( const Data& d ) : flags( d.flags ), data( d.data ) {}
        Data& operator=( const Data& d )
          { if ( &d != this ) 
             { flags = d.flags; data = d.data; }
          return *this; }
#else
        Data() :
            flags              (0U),
            data               (0U),
            scalars            (0U),
            vectors            (0U),
            tensors            (0U),
            arrays             (0U),
            arrayLength        (0U),
            flaggedArrays      (0U),
            flaggedArrayLength (0U)
        {}

        Data( const Data& d ) :
            flags               ( d.flags ),
            data                ( d.data ),
            scalars             ( d.scalars ),
            vectors             ( d.vectors ),
            tensors             ( d.tensors ),
            arrays              ( d.arrays ),
            arrayLength         ( d.arrayLength ),
            flaggedArrays       ( d.flaggedArrays ),
            flaggedArrayLength  ( d.flaggedArrayLength )
        {}

        // SKM ADDITION
        Data( Data&& d ) :
            flags{ d.flags },
            data{ d.data },
            scalars{ d.scalars },
            vectors{ d.vectors },
            tensors{ d.tensors },
            arrays{ d.arrays },
            arrayLength{ d.arrayLength },
            flaggedArrays{ d.flaggedArrays },
            flaggedArrayLength{ d.flaggedArrayLength }
        {}

        Data& operator=( const Data& d )
          {
            if ( &d != this ) {
                scalars             = d.scalars;
                vectors             = d.vectors;
                tensors             = d.tensors;
                arrays              = d.arrays;
                arrayLength         = d.arrayLength;
                flaggedArrays       = d.flaggedArrays;
                flaggedArrayLength  = d.flaggedArrayLength;
                flags               = d.flags;
                data                = d.data;
             }
            return *this;
         }

        // local variable state
        size_t  scalars,
                vectors,
                tensors,
                arrays,
                flaggedArrays;
        size_t  arrayLength,
                flaggedArrayLength;
#endif
      };

    /// copys supplied data into storage
    void LVS( const Data& data ) { data_ = data; }
    /// returns a copy of the data container for operations
    const Data LVS() const { return data_; }
  
  private:

    Data data_; ///< data and flag containers
};

} // end csmp

#endif
