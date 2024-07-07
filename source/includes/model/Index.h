#ifndef CSMP_INDEX_H
#define CSMP_INDEX_H

#include <iostream>
#include <fstream>
#include <climits>
#include <cstdio>
#include <set>

#include "CSMP_global_enumerations.h"
#include "LocalVariables.h"
#include "IntegrationPointVariables.h"

// additional index output statements and checks in variable storage
// #define VARIABLE_STORAGE_DEBUG

namespace csmp {

class IndexTracker;

/** @brief Accessor for discretized physical variables

  Link between PropertyDatabase and LocalVariableStorage.
 
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

  @todo SKM: (3-F) complete scheme for compile-time variables
 
*/
struct Index {

    using int_type = LocalVariables::int_type;

    Index();
    Index( VARIABLE_TYPE, PLACEMENT, int_type idx );
  
    Index( VARIABLE_TYPE, PLACEMENT, int_type idx, int_type dataDepth, int_type flagDepth, int_type dataOffset, int_type flagOffset,
           const LocalVariables&, const IntegrationPointVariables& = IntegrationPointVariables(),
           int_type offsetFactorObject = 0, int_type offsetFactorSector = 0,
           int_type ipFactorObject = 0, int_type ipFactorSector = 0, int_type ipFactorFacet = 0 );
  
    Index( const csmp::Index& );
    Index( csmp::Index&& );
    ~Index();

    Index&  operator=( const csmp::Index& );
    Index&  operator=( csmp::Index&& );
    bool    operator==( const csmp::Index& ) const;
    bool    operator!=( const csmp::Index& ) const; 
    bool    operator<( const csmp::Index& ) const; 
    bool    IsDefined() const;

    void Attach( IndexTracker* indexTracker );
    void Detach();
    void UpdateData( const csmp::Index& idx );

    void Out() const;
    bool Out( std::fstream& fp ) const;
    bool In( std::fstream& fp );

    // TODO: make all these uint8_t (same size as the flags?
    VARIABLE_TYPE               type;
    PLACEMENT                   place;
    int_type                    index = UNSPECIFIED;        ///< For Scalars, Vectors, Tensors, Arrays, FlaggedArrays: the how many'th variable of its kind at specified placement
    int_type                    dataDepth;                  ///< Scalar:1 , Vector: dim, Tensor: dim*dim, Array:Size, FlaggedArray:Size
    int_type                    flagDepth;                  ///< Scalar:1 , Vector: dim, Tensor: dim, Array:1, FlaggedArray:Size
    int_type                    dataOffset;                 ///< Index in data container where data start
    int_type                    flagOffset;                 ///< Index in flag container where data start
    int_type                    offsetFactorCell;           ///< Factors used in integration point variable index arithmetic
    int_type                    offsetFactorSector;         ///< Factors used in integration point variable index arithmetic
    int_type                    ipFactorCell;               ///< Factors used in integration point variable index arithmetic
    int_type                    ipFactorSector;             ///< Factors used in integration point variable index arithmetic
    int_type                    ipFactorFacet;              ///< Factors used in integration point variable index arithmetic
    LocalVariables              localVariables;             ///< Description of state of physical variables at given placement (variables count etc..)
    IntegrationPointVariables   integrationPointVariables;  ///< Description of state of physical variables at given placement (variables count etc..)
    IndexTracker*               indexTracker;               ///< IndexTracker used for runtime updates of offsets (addition/removal of variables)
 };


  /// Templatized version allowing for compile-time type selection
template<VARIABLE_TYPE ty,PLACEMENT pl>
struct INDEX : public Index {
   static constexpr VARIABLE_TYPE VariableType = ty;
   static constexpr PLACEMENT VariablePlacement = pl;

   /// compile time construction of index to be used in factory implementations
   INDEX( int_type idx, int_type dataDepth, int_type flagDepth, int_type dataOffset, int_type flagOffset,
          const LocalVariables& lvars, const IntegrationPointVariables& ivars = IntegrationPointVariables(),
          int_type offsetFactorObject = 0, int_type offsetFactorSector = 0,
          int_type ipFactorObject = 0, int_type ipFactorSector = 0, int_type ipFactorFacet = 0 )
    : Index( ty, pl, idx, dataDepth, flagDepth, dataOffset, flagOffset,
             lvars, ivars,
             offsetFactorObject, offsetFactorSector,
             ipFactorObject, ipFactorSector, ipFactorFacet ) {}
              
   INDEX() : Index() {}
   explicit INDEX( int_type i ) : Index(ty,pl,i) {}
   explicit INDEX( csmp::Index& idx ) : Index(idx) {}
   explicit INDEX( csmp::Index&& idx ) : Index(idx) {}
};


// type parsing (although slow)
template<typename Var>  VARIABLE_TYPE variableType( const Var& );


// console input and output
std::ostream& operator<<( std::ostream&, const csmp::Index& );
std::istream& operator>>( std::istream&, csmp::Index& );

} // csmp


namespace csmp {

/**

  @namespace csmp
  @brief Complex System Modeling Platform

  The Complex System Modelling Platform (CSMP++) is an object-oriented
  finite-element based application programmer interface (API),
  designed for the simulation of complex geological processes and their interactions.
  CSMP++ is available on a license for commercial use,
  and is available for free to the academic user.
 
  @author S.K. Matthai
  @date 1994.
 
*/

} // end csmp

#endif
