// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_GLOBAL_ENUMERATIONS_H
#define CSMP_GLOBAL_ENUMERATIONS_H

#include <set>
#include <string>
#include <limits>
#include <cstdint>

namespace csmp {

/*
========================================
CSMP global constants and enumerations
========================================
*/

/// initialisation of const to maximum value that uint32_t  can take
constexpr size_t   NULL_IDX(std::numeric_limits<size_t>::max());
constexpr uint32_t NULL_IDX32U(std::numeric_limits<uint32_t>::max());
constexpr int      UNSPECIFIED(-1);

/**
@defgroup CSMP_global_enumerations
Enums available globally in the csmp Namespace
 */


/**
@addtogroup CSMP_global_enumerations
@{
*/

/// some fixed string sizes (max filename length on NTFS partition / arbitrary #<32k)
constexpr int NAME_STRING(256), INFO_STRING(1024);

/// for templates and specific initialisations (set compiler to treat enums as single byte types!)
enum ONE_BYTE_NUMBER : std::int8_t { ZERO, ONE, TWO, THREE,
                                     FOUR, FIVE, SIX, SEVEN, EIGHT, NINE, TEN,
                                     ELEVEN, TWELVE, THIRTEEN, FOURTEEN, FIVETEEN, SIXTEEN, SEVENTEEN,
                                     EIGHTEEN, NINETEEN, TWENTY, TWENTY_ONE, TWENTY_TWO, TWENTY_THREE,
                                     TWENTY_FOUR, TWENTY_FIVE, TWENTY_SIX, TWENTY_SEVEN, TWENTY_EIGHT,
                                     TWENTY_NINE, THIRTY, THIRTY_ONE, THIRTY_TWO, THIRTY_THREE, THIRTY_FOUR,
                                     THIRTY_FIVE, THIRTY_SIX, THIRTY_SEVEN, THIRTY_EIGHT, THIRTY_NINE, FOURTY,
                                     FOURTY_ONE, FOURTY_TWO, FOURTY_THREE, FOURTY_FOUR, FOURTY_FIVE, FOURTY_SIX,
                                     FOURTY_SEVEN, FOURTY_EIGHT, FOURTY_NINE, FIFTY, NOT_INITIALIZED=-1 };
/**
@}
*/

/**
@addtogroup CSMP_global_enumerations
@{ */

/// Discretized csmp variables
enum VARIABLE_TYPE : std::int8_t { SCALAR=1, VECTOR=2, TENSOR=3, ARRAY=4, FLAGGEDARRAY=5 };

//enum VARIABLE_FLAG : int_fast8_t
/// variable flag indicating treatment in computations (forced to be one-byte size because it is stored everywhere)
enum VARIABLE_FLAG : std::int8_t { PLAIN,          /**< modifyable, dependent or indep. var. */
                                   ANY,            /**< unspecified discriminator.  This is also used to identify if a variable has been assigned a value or not*/
                                   INIT_GUESS,     /**< convergence oriented not phys. meaningful, not checked  */
                                   INIT_COND,      /**< physically meaningful initial condition, checked */
                                   FIELD_DATA,     /**< (Geological) field data for comparison */
                                   PERIODIC,       /**< linked via ID to opposite side of model */
                                   DIRICH,         /**< Dirichlet boundary condition */
                                   NEUMANN,        /**< Neumann boundary condition */
                                   ROBIN,          /**< linear combination function and derivative values on the boundary of the domain */
                                   CONSTANT_FLUX   /**< Neumann is prescribed gradient, this one is for prescribed flux */
                                 };

/// Top down hierarchy of variable placements on mesh tree; finite volume quadrature pointe embedded
enum PLACEMENT : std::int8_t { UNDEFINED, // default
                               MODEL,     // increasing hight of tree
                               REGION,
                               BOUNDARY,
                               SPLIT_BOUNDARY,
                               ELEMENT,
                               ELEMENT_INTEGRATION_POINT,
                               SECTOR_INTEGRATION_POINT,
                               FACET_INTEGRATION_POINT,
                               FACE,
                               FACE_INTEGRATION_POINT,
                               FACE_SECTOR_INTEGRATION_POINT,
                               FACE_FACET_INTEGRATION_POINT,
                               EDGE, ///< only in 3D
                               INTER_FACE,
                               INTER_FACE_INTEGRATION_POINT,
                               INTER_FACE_SECTOR_INTEGRATION_POINT,
                               INTER_FACE_FACET_INTEGRATION_POINT,
                               NODE
                            };

/// to classify 
enum CELL_SHAPE : std::int8_t { POINT=0, LINE=1, SURFACE=2, VOLUME=3, HYPER_DIMENSIONAL=4 };

/// Side of lower-dimensional face or interface between two higher-dimensional elements ( INSIDE or OUTSIDE ) and a potential lower-dimensional parent element ( MIDDLE )
enum INTERFACE_SIDE : std::int32_t
{
    INSIDE  = -1,
    OUTSIDE =  1,
    MIDDLE  =  0
};



/**
    Classification of ManifoldType is about:
      •	how many independent regions meet at the Manifold
      •	and whether that meeting is symmetric or partial
      •	this classification applies to all nodes within Manifold simultaneously (including intervening ones)
      
    Since NodeManifolds also support node-centered finite volumes, even a STANDALONE NodeManifold is  physically meaningful.
    Even if:
      •	there is no InterFace
      •	there is no lower-dimensional FE entity
      •	there is only one branch

This is crucial.

  @note SplitBoundary objects exist only on the inside of models
  @note NodeManifolds exist only at SplitBoundary objects
  @attention the perimeter nodes of a SplitBoundary are not necessarily Manifolds; only if they contain lower-dim Region objects or SplitBoundaries cross
  @note The ManifoldType classification is for topologically collocated Nodes only!
  @note since SplitBoundaries are surfaces in 3D, surface is their highest dimension

*/
enum class ManifoldType : int8_t {
                                    STAND_ALONE,                ///<  collocated node set whose parent elements belong to exactly one subdomain each
                                    SPLIT_BOUNDARY,             ///<  2-node manifold along a SplitBoundary (most common)
                                    SPLIT_BOUNDARY_WITH_INTERNAL_MESH,
                                    SPLIT_BOUNDARY_CROSSING,    ///<  4-node manifold intersection of split boundaries
                                    MULTI_SB_CROSSING,          ///<  6-node cross of 3 SBs in
                                    SPLIT_BOUNDARY_END,         ///<  termination against Boundary or another SplitBoundary ( T-intersection)
                                };

/// converts classifiers to strings so that they can be printed
std::string parse( ManifoldType );



/**
      Geometric classification of nodes / points, BREP stands for boundary representation.
       @author SKM
       @date 9/07/2022
       
      In detail

      Value     Name   Description    Remeshing Constraint

      0 MESH_VERTEX Standard vertex inside a volume/surface. Fully mobile; can be deleted or moved to optimize quality.

      1 INTERIOR_POINT Where multiple lines cross or surfaces meet at a single point. Fixed. Cannot be removed or moved; defines the B-Rep skeleton.

      2 PERIMETER_POINT A point at the end of a line (0D boundary of a 1D line). Fixed. Defines the termination of a geometric feature.

      3 EXTERIOR_POINT  A point on the outside surface of the model. Usually redundant if Surface/Line flags exist, but marks "Hull" corners.

      4 INTERIOR_LINE  A vertex on a line located entirely inside a volume. Can slide along the line; cannot be moved off the line.

      5 PERIMETER_LINE A surface edge that exists inside the model (e.g., a hole's rim). Can slide along the perimeter; preserves the "cutout" shape.

      6 EXTERIOR_LINE An edge belonging to the outer hull of the model. Can slide along the edge; preserves the visual boundary.

      7 INTERIOR_SURFACE A vertex on a surface separating two volumes (Subdomain). Can move within the surface plane; preserves material interface.

      8 PERIMETER_SURFACE A surface forming an internal "island" or hull. Can move within the surface; preserves internal cavity shapes.

      9 EXTERIOR_SURFACE  A vertex on the outermost shell of the model. Can move within the surface; preserves the total model volume.

*/
enum TOPOTYPE : std::int8_t {
                                MESH_VERTEX,       ///< a point within the model volume
                                INTERIOR_POINT,    ///< a point where lines or surfaces touch or cross
                                PERIMETER_POINT,   ///< point at the end of a line inside a 2D model
                                EXTERIOR_POINT,    ///< on an outside surface of the model
                                INTERIOR_LINE,     ///< a line on the interior of the model
                                PERIMETER_LINE,    ///< a surface edge inside of the model
                                EXTERIOR_LINE,     ///< an edge of the model
                                INTERIOR_SURFACE,  ///< a surface within a 3D model
                                PERIMETER_SURFACE, ///< a surface forming the hull of an object inside of a 3D model
                                EXTERIOR_SURFACE   ///< a surface delimiting a 3D model
                            };

/// for testing  int8_t values prior to type-casting to TOPOTYPE
constexpr bool isValidTopotype(std::int8_t v) {
    switch (static_cast<TOPOTYPE>(v)) {
        case TOPOTYPE::MESH_VERTEX:
        case TOPOTYPE::INTERIOR_POINT:
        case TOPOTYPE::PERIMETER_POINT:
        case TOPOTYPE::EXTERIOR_POINT:
        case TOPOTYPE::INTERIOR_LINE:
        case TOPOTYPE::PERIMETER_LINE:
        case TOPOTYPE::EXTERIOR_LINE:
        case TOPOTYPE::INTERIOR_SURFACE:
        case TOPOTYPE::PERIMETER_SURFACE:
        case TOPOTYPE::EXTERIOR_SURFACE:
            return true;
        default:
            return false;
    }
}


/// converts classifiers to strings so that they can be printed
std::string parseTopology( TOPOTYPE );

/**
      Detects any Node that could be on the perimeter of a ModelSubDomain.
*/
template<uint32_t dim>
inline constexpr bool isPerimeterNode( TOPOTYPE gflag ) noexcept {
   if constexpr( dim == 3 ) {
        switch( gflag ) {
          case PERIMETER_SURFACE:
          case EXTERIOR_SURFACE:
//          case INTERIOR_SURFACE: false because this should not prompt the generation of a manifold
          case PERIMETER_LINE:
          case PERIMETER_POINT:
          case EXTERIOR_LINE:
          case EXTERIOR_POINT:
          case INTERIOR_POINT:
            return true;
          default:
            return false;
        }
     }
   else if constexpr( dim == 2 ) {
        switch( gflag ) {
          case PERIMETER_LINE:
          case PERIMETER_POINT:
          case EXTERIOR_LINE:
          case EXTERIOR_POINT:
          case INTERIOR_LINE:
          case INTERIOR_POINT:
            return true;
          default:
            return false;
        }
     }
   else {
         switch( gflag ) {
          case PERIMETER_POINT:
          case EXTERIOR_POINT:
          case INTERIOR_POINT:
            return true;
          default:
            return false;
        }
   }
} // end isPerimeterNode


template<uint32_t dim>
inline constexpr bool isPerimeterNodeInsideModel( TOPOTYPE gflag ) noexcept {
   if constexpr( dim == 3 ) {
        switch( gflag ) {
          case PERIMETER_SURFACE:
          case PERIMETER_LINE:
          case PERIMETER_POINT:
          case INTERIOR_SURFACE:
          case INTERIOR_LINE:
          case INTERIOR_POINT:
            return true;
          default:
            return false;
        }
     }
   else if constexpr( dim == 2 ) {
        switch( gflag ) {
          case PERIMETER_LINE:
          case PERIMETER_POINT:
          case INTERIOR_LINE:
          case INTERIOR_POINT:
            return true;
          default:
            return false;
        }
     }
   else {
         switch( gflag ) {
          case PERIMETER_POINT:
          case INTERIOR_POINT:
            return true;
          default:
            return false;
        }
   }
} // end isPerimeterNode





/// variable flag indicating treatment in computations
enum NORM_INDEX : std::int8_t
{
    POINT_DIFFERENCE           = -1,
    ABSOLUTE_POINT_DIFFERENCE  = 0,
    L1_NORM              = 1,
    MANHATTAN_NORM       = L1_NORM,
    TAXICAB_NORM         = L1_NORM,
    L2_NORM              = 2,
    EUCLIDIAN_NORM       = L2_NORM,
    LP_NORM              = 3,
    MAX_NORM             = 4,
    INF_NORM             = MAX_NORM,
    LINF_NORM            = MAX_NORM,
    L21_NORM             = 5,
    L22_NORM             = 6,
    FROBENIUS_NORM       = L22_NORM,
    HILBERT_SCHMIDT_NORM = L22_NORM,
    LPQ_NORM             = 7,
    INDUCED_L1_NORM      = 8,
    INDUCED_L2_NORM      = 9,
    INDUCED_MAX_NORM     = 11,
    INDUCED_LINF_NORM    = INDUCED_MAX_NORM
};

/** @} */


/**
@addtogroup CSMP_global_enumerations
@{ */

/// Fills set with all current PLACEMENT choices
void variablePlacementSet( std::set<PLACEMENT>& plSet );

/// Fills set with all finite-volume associated PLACEMENT choices
void fvVariablePlacementSet( std::set<PLACEMENT>& plSet );

/// Fills set with all VARIABLE_TYPE choices in current model
void variableTypeSet( std::set<VARIABLE_TYPE>& );

// conversion utilities

/// variable type
VARIABLE_TYPE  parseType( int type );
std::string    parseType( VARIABLE_TYPE );
VARIABLE_TYPE  parseType( const char* type );

class ScalarVariable;
template<uint32_t> class VectorVariable;
template<uint32_t> class TensorVariable;
class ArrayVariable;
class FlaggedArrayVariable;

// Helper to identify VectorVariable and TensorVariable (Template template parameters)
template <typename T> struct IsVector : std::false_type {};
template <uint32_t dim> struct IsVector<VectorVariable<dim>> : std::true_type {};

template <typename T> struct IsTensor : std::false_type {};
template <uint32_t dim> struct IsTensor<TensorVariable<dim>> : std::true_type {};

/**
    Csmp types are resolved at compile time via if-constexpr.
 */
template<typename Var>
constexpr VARIABLE_TYPE variableType( const Var& ) {
    using T = std::decay_t<Var>;

    if constexpr (std::is_same_v<T, ScalarVariable>) {
        return SCALAR;
    }
    else if constexpr (IsVector<T>::value) {
        return VECTOR;
    }
    else if constexpr (IsTensor<T>::value) {
        return TENSOR;
    }
    else if constexpr (std::is_same_v<T, ArrayVariable>) {
        return ARRAY;
    }
    else if constexpr (std::is_same_v<T, FlaggedArrayVariable>) {
        return FLAGGEDARRAY;
    }
    // double, integer, string etc.
    else {
        // tighter tests
        static_assert( std::floating_point<T> == true || std::integral<T>, "variable seems to be on of the inbuilt types like double or int");
        // This will trigger at compile-time if you pass a type not listed above.
        static_assert(sizeof(T) == 0, "variableType(const Var&): passed unsupported variable type argument");
        return SCALAR; // Fallback to satisfy return type (will never be reached)
    }
}


/// variable placement
PLACEMENT      intToPLACEMENT( int i ) noexcept;
PLACEMENT      parsePlacement( const char* placement ) noexcept;
std::string    parsePlacement( PLACEMENT ) noexcept;
bool           faceVariable( PLACEMENT ) noexcept;
bool           interFaceVariable( PLACEMENT ) noexcept;

/// determine from the element type whether the placement of the variable is Region, Boundary or SplitBoundary
template<uint32_t dim, template<uint32_t> class PLACE>
PLACEMENT      parsePlacement() noexcept;

bool           isPlacedOnIntegrationPoint( PLACEMENT ) noexcept;
/// variable flag
VARIABLE_FLAG  intToVARIABLE_FLAG( int i ) noexcept;
VARIABLE_FLAG  parseCondition( std::string_view );
VARIABLE_FLAG  parseStatus( const char* status ) noexcept;
std::string    parseStatus( VARIABLE_FLAG ) noexcept;

/// returns the side of the Face or InterFace element
std::string    parseSide( INTERFACE_SIDE ) noexcept;

/** @} */


} // end namespace csmp

#endif
