#ifndef CSMP_INDEX_H
#define CSMP_INDEX_H

#include <iostream>
#include <climits>
#include <cstdio>
#include <set>

#include "CSMP_global_enumerations.h"
#include "CSMP_number_types.h"
#include "IndexTracker.h"
#include "LocalVariables.h"
#include "IntegrationPointVariables.h"

// additional index output statements and checks in variable storage
// #define VARIABLE_STORAGE_DEBUG

namespace csmp {

  /** @brief Accessor for discretized physical variables

  Link between PropertyDatabase and LocalVariableStorage.

  @todo SKM: (3-F) Provide scheme for compile-time variables
  @todo (2-D) Should become class
  */
  struct Index {
    Index(); 
    Index( VARIABLE_TYPE, PLACEMENT, size_t idx ); 
    Index( VARIABLE_TYPE, PLACEMENT, size_t idx, size_t dataDepth, size_t flagDepth, size_t dataOffset, size_t flagOffset, 
      const LocalVariables&, const IntegrationPointVariables& = IntegrationPointVariables(), size_t offsetFactorSimplex = 0, 
      size_t offsetFactorSector = 0, size_t ipFactorSimplex = 0, size_t ipFactorSector = 0, size_t ipFactorFacet = 0 ); 
    Index( const csmp::Index& );
    ~Index();

    Index&  operator=( const csmp::Index& );
    bool    operator==( const csmp::Index& ) const;
    bool    operator!=( const csmp::Index& ) const; 
    bool    operator<( const csmp::Index& ) const; 
    bool    IsDefined() const;

    void Attach( IndexTracker* indexTracker );
    void Detach();
    void UpdateData( const csmp::Index& idx );

    void Out(std::ostream& os) const;
    bool Out( FILE* fp ) const;
    bool In( FILE* fp );


    VARIABLE_TYPE               type;
    PLACEMENT                   place;
    size_t                      index;                      ///< For Scalars, Vectors, Tensors, Arrays, FlaggedArrays: the how many'th variable of its kind at specified placement
    size_t                      dataDepth;                  ///< Scalar:1 , Vector: dim, Tensor: dim*dim, Array:Size, FlaggedArray:Size
    size_t                      flagDepth;                  ///< Scalar:1 , Vector: dim, Tensor: dim, Array:1, FlaggedArray:Size
    size_t                      dataOffset;                 ///< Index in data container where data start
    size_t                      flagOffset;                 ///< Index in flag container where data start
    size_t                      offsetFactorSimplex;        ///< Factors used in integration point variable index arithmetic
    size_t                      offsetFactorSector;         ///< Factors used in integration point variable index arithmetic
    size_t                      ipFactorSimplex;            ///< Factors used in integration point variable index arithmetic
    size_t                      ipFactorSector;             ///< Factors used in integration point variable index arithmetic
    size_t                      ipFactorFacet;              ///< Factors used in integration point variable index arithmetic
    LocalVariables              localVariables;             ///< Description of state of physical variables at given placement (variables count etc..)
    IntegrationPointVariables   integrationPointVariables;  ///< Description of state of physical variables at given placement (variables count etc..)
    IndexTracker*               indexTracker;               ///< IndexTracker used for runtime updates of offsets (addition/removal of variables)
    };   


  /// Templatized version with allows for compile-time type selection
  template<VARIABLE_TYPE ty,PLACEMENT pl> 
  struct INDEX : public Index {
    explicit INDEX( size_t i ) : Index(ty,pl,i) {}
    };


  // type parsing (although slow)
  template<typename Var>  VARIABLE_TYPE variableType( const Var& );


  // console input and output
  std::ostream& operator<<( std::ostream&, const csmp::Index& );
  std::istream& operator>>( std::istream&, csmp::Index& );

  } // csmp




namespace csmp {

  /**
  @struct Index Index.h "main_library/Index.h"

  @author S.K. Matthaei
  @author P. Lang
  @author S. Geiger
  @author Stephen G. Roberts

  @section motivation Motivation

  Class encapsulates the three values that are required to retrieve a
  CSP variable from a model based on an unstructured mesh: The Type()
  distinguishes scalars, vectors and tensors from one-another since these
  have different storage requirements.

  The Place() specifies whether a variable is discretized on the nodes,
  elements, faces or on constraint points. This has implications as to
  how the variable is interpolated over the model domain.

  The Index() is used internally by CSMP for housekeeping. It corresponds
  to the number that the variable has in an alphabetically ordered
  map of variables of the same type and placement, in the same model.


  @section design Design Intent

  Create a universal accessor for discretized physical variables in a
  model.


  @section applicability Applicability

  Typically, a csmp::Index or Index object is initialized only once
  before it is used as an accessor to variables in a model. The
  initialization is done by the PropertyDatabase if so requested.

  Index is applicable whereever the parallel structure csmp::Index is
  applicable. The reason, why a class is present in addition to the
  structure is that this allows to assign default values when indices
  are constructed.


  @section collaborations Collaborations

  Index is known to the PropertyDatabase and all variable accessor
  and mutator functions. It is used by PropertyHandle and Operand objects,
  as well as Interrelations, Visitors and Algorithms to housekeep
  currently used variables.


  @section implementation Implementation

  Both as struct csmp::Index and this class.


  @section examples Application Examples

  To get an Index for a variable of interest in order to perform
  a computation, do the following:

  @code
  Index variable_key = prop_database.StorageKey("variable");
  @endcode

  */

  /** Default constructor

  Assigns an Index of scalar type on the model.
  The index number itself will be set to the maximum
  and hence interpreted as not defined.
  Inline implementation.
  */

  inline Index::Index() 
    : type(SCALAR), 
      place(UNDEFINED), 
      index(ULONG_MAX),
      dataDepth(0),
      flagDepth(0),
      dataOffset(0),
      flagOffset(0),
      offsetFactorSimplex(0),
      offsetFactorSector(0), 
      ipFactorSimplex(0),
      ipFactorSector(0),
      ipFactorFacet(0),
      localVariables(),
      integrationPointVariables(),
      indexTracker(NULL)
    {}


  /** Custom constructor

  Neither sets local nor integration point variables and has no option for offsets/depth ( i.e. = 0 )

  Assigns an Index with defined type and placement.
  Inline implementation.
  @param ty the variable type
  @param pl the placement
  @param idx the index

  */
  inline Index::Index( VARIABLE_TYPE ty, PLACEMENT pl, size_t idx ) 
    : type(ty), place(pl),
      index(idx), dataDepth(0), flagDepth(0), dataOffset(0), 
      flagOffset(0), offsetFactorSimplex(0), offsetFactorSector(0), 
      ipFactorSimplex(0), ipFactorSector(0), ipFactorFacet(0),
      localVariables(), integrationPointVariables(), indexTracker(NULL)
    {}


  /** Custom constructor

  This should be goto ctor, even though members may be initialized individually through struct
  Assigns an Index with defined type and placement.
  Inline implementation.

  @param ty the variable type
  @param pl the placement
  @param idx the index
  @param datadepth component count, size of variable (scalar 1, vector 2, tensor 3, array #components)
  @param dataoffset Offset in of data in variable storage container (implicit dependence on LocalVariableStorage)
  @param flagoffset Offset in of flag data in variable storage container (implicit dependence on LocalVariableStorage)
  @param lvs LocalVariables of placement pl
  @param ivs IntegrationPointVariables of elements
  */
  inline Index::Index( VARIABLE_TYPE ty, PLACEMENT pl, size_t idx, 
                        size_t datadepth, size_t flagdepth, size_t dataoffset, size_t flagoffset,
                        const LocalVariables& lvs, const IntegrationPointVariables& ivs, 
                        size_t offsetfactorsimplex, size_t offsetfactorsector, size_t ipfactorsimplex, 
                        size_t ipfactorsector, size_t ipfactorfacet ) 
    : type(ty), place(pl),
      index(idx), dataDepth(datadepth), flagDepth(flagdepth), dataOffset(dataoffset), 
      flagOffset(flagoffset), offsetFactorSimplex(offsetfactorsimplex), offsetFactorSector(offsetfactorsector), 
      ipFactorSimplex(ipfactorsimplex), ipFactorSector(ipfactorsector), ipFactorFacet(ipfactorfacet),
      localVariables(lvs), integrationPointVariables(ivs), indexTracker(NULL)
    {}

  /** Copy constructor

  Creates an index based on the provided parameter
  @param idx an existing csmp::Index object
  */
  inline Index::Index( const csmp::Index& idx ) 
    : type(idx.type), place(idx.place), index(idx.index),
      dataDepth(idx.dataDepth), flagDepth(idx.flagDepth), dataOffset(idx.dataOffset), flagOffset(idx.flagOffset),
      offsetFactorSimplex(idx.offsetFactorSimplex), offsetFactorSector(idx.offsetFactorSector), 
      ipFactorSimplex(idx.ipFactorSimplex), ipFactorSector(idx.ipFactorSector), ipFactorFacet(idx.ipFactorFacet),
      localVariables(idx.localVariables), integrationPointVariables(idx.integrationPointVariables),
      indexTracker(NULL)
    {
    // attach itself to IndexTracker and vice vers
    if(idx.indexTracker)
      {
        Detach();
        idx.indexTracker->Attach( this, &idx ); 
        Attach( idx.indexTracker );
      }
    }


  /// Assignment operator
  inline Index&  Index::operator=( const csmp::Index& idx ) {
    if ( this != &idx ) {
      type                      = idx.type;
      place                     = idx.place;
      index                     = idx.index;
      dataDepth                 = idx.dataDepth;
      flagDepth                 = idx.flagDepth;
      dataOffset                = idx.dataOffset;
      flagOffset                = idx.flagOffset;
      offsetFactorSimplex       = idx.offsetFactorSimplex;
      offsetFactorSector        = idx.offsetFactorSector;
      ipFactorSimplex           = idx.ipFactorSimplex;
      ipFactorSector            = idx.ipFactorSector;
      ipFactorFacet             = idx.ipFactorFacet;
      localVariables            = idx.localVariables;
      integrationPointVariables = idx.integrationPointVariables;
      /// attach itself to IndexTracker and vice versa.
      if( idx.indexTracker )
        { 
          Detach();
          idx.indexTracker->Attach( this, &idx );
          Attach( idx.indexTracker );
        }
      }        

    return *this;
    }

  /**
   @fn  inline void Index::UpdateData( const csmp::Index& idx )
  
   @brief Updates to the data described by idx without updating indexTracker
  
   @author  P. Lang
   @date  9/26/2012
  
   @param idx The index to be updated from
   */
  inline void Index::UpdateData( const csmp::Index& idx )
    {
    if ( this != &idx ) {
      type                      = idx.type;
      place                     = idx.place;
      index                     = idx.index;
      dataDepth                 = idx.dataDepth;
      flagDepth                 = idx.flagDepth;
      dataOffset                = idx.dataOffset;
      flagOffset                = idx.flagOffset;
      offsetFactorSimplex       = idx.offsetFactorSimplex;
      offsetFactorSector        = idx.offsetFactorSector;
      ipFactorSimplex           = idx.ipFactorSimplex;
      ipFactorSector            = idx.ipFactorSector;
      ipFactorFacet             = idx.ipFactorFacet;
      localVariables            = idx.localVariables;
      integrationPointVariables = idx.integrationPointVariables;
      }
    }




  /** Boolean equals operator

  @param i an csmp::Index object
  @return boolean result of the comparison
  @return true if equal
  */
  inline bool Index::operator==( const csmp::Index& i ) const 
    {
    return (  type==i.type && place==i.place && index==i.index && dataDepth==i.dataDepth && flagDepth==i.flagDepth && dataOffset==i.dataOffset 
      && flagOffset==i.flagOffset && offsetFactorSimplex==i.offsetFactorSimplex && offsetFactorSector==i.offsetFactorSector );
    }



  /** Boolean equals not operator

  @return boolean result of the comparison
  @return true if not equal
  */
  inline bool Index::operator!=( const csmp::Index& i ) const 
    {
    return ( type!=i.type || place!=i.place || index!=i.index || dataDepth!=i.dataDepth || flagDepth!=i.flagDepth || dataOffset!=i.dataOffset 
      || flagOffset!=i.flagOffset || offsetFactorSimplex!=i.offsetFactorSimplex || offsetFactorSector!=i.offsetFactorSector );
    }


  /** Boolean Function that checks if object has been defined

  @return true if an index has been assigned
  */
  inline bool Index::IsDefined() const {
    return (index != ULONG_MAX);
    }


  /// Registers an IndexTracker (does not detach from current!)
  inline void Index::Attach( IndexTracker* newTracker )
    {
    indexTracker = newTracker;
    }


  /// DTor: Detaches intself from IndexTracker
  inline Index::~Index()
    {
    if(indexTracker)
      indexTracker->Detach(this);
    }


  /// Detaches itself from IndexTracker (i.e. called upon destruction)
  inline void Index::Detach()
    {
    if(indexTracker)
      indexTracker->Detach(this);
    indexTracker = NULL;
    }


  /**
  @namespace csmp
  @brief Complex System Modeling Platform

  The Complex System Modelling Platform (CSMP++) is an object-oriented
  finite-element based application programmer interface (API),
  designed for the simulation of complex geological processes and their interactions.
  CSMP++ is available on a license for commercial use,
  and is available for free to the academic user.
  @author S.K. Matthaei
  */

  } // end csmp

#endif
