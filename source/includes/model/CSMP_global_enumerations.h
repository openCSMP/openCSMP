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

/// initialisation of const to maximum value that size_t can take
const size_t NULL_IDX(std::numeric_limits<size_t>::max());
const short  UNSPECIFIED(-1);

/**
@defgroup CSMP_global_enumerations
Enums available globally in the csmp Namespace
 */


/**
@addtogroup CSMP_global_enumerations
@{
*/

/// some fixed string sizes (max filename length on NTFS partition / arbitrary #<32k)
enum CSMP_STRING_LENGTHS { NAME_STRING=256, INFO_STRING=1024 };

/// for templates and specific initialisations (set compiler to treat enums as single byte types!)
enum ONE_BYTE_NUMBER { ZERO, ONE, TWO, THREE,
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
enum VARIABLE_TYPE { SCALAR=1, VECTOR=2, TENSOR=3, ARRAY=4, FLAGGEDARRAY=5 };

//enum VARIABLE_FLAG : int_fast8_t
/// variable flag indicating treatment in computations (forced to be one-byte size because it is stored everywhere)
enum VARIABLE_FLAG { PLAIN,          /**< modifyable, dependent or indep. var. */
                     ANY,            /**< unspecified discriminator.  This is also used to identify if a variable has been assigned a value or not*/
                     INIT_GUESS,     /**< convergence oriented not phys. meaningful, not checked  */
                     INIT_COND,      /**< physically meaningful initial condition, checked */
                     FIELD_DATA,     /**< (Geological) field data for comparison */
                     PERIODIC,       /**< linked via ID to opposite side of model */
                     DIRICH,         /**< Dirichlet boundary condition */
                     NEUMANN,        /**< Neumann boundary condition */
                     ROBIN,          /**< linear combination function and derivative values on the boundary of the domain */
                     CONSTANT_FLUX}; /**< Neumann is prescribed gradient, this one is for prescribed flux */

/// Top down hierarchy of variable placements on mesh tree; finite volume quadrature pointe embedded
enum PLACEMENT { UNDEFINED, // default
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
                 INTER_FACE,
                 INTER_FACE_INTEGRATION_POINT,
                 INTER_FACE_SECTOR_INTEGRATION_POINT,
                 INTER_FACE_FACET_INTEGRATION_POINT,
                 NODE,
                 SPLIT_NODE // TODO: check: is this still required in the Melbourne build
              };


/// Side of lower-dimensional face or interface between two higher-dimensional elements ( INSIDE or OUTSIDE ) and a potential lower-dimensional parent element ( MIDDLE )
enum INTERFACE_SIDE
{
    INSIDE  = -1,
    OUTSIDE =  1,
    MIDDLE  =  0
};


/// variable flag indicating treatment in computations
enum NORM_INDEX
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
template<typename csmp_type>
VARIABLE_TYPE  variableType();

/// variable placement
PLACEMENT      intToPLACEMENT( int i );
PLACEMENT      parsePlacement( const char* placement );
std::string    parsePlacement( PLACEMENT );
bool           faceVariable( PLACEMENT );
bool           interFaceVariable( PLACEMENT );

/// determine from the element type whether the placement of the variable is Region, Boundary or SplitBoundary
template<size_t dim, template<size_t> class PLACE>
PLACEMENT      parsePlacement();

bool           isPlacedOnIntegrationPoint( PLACEMENT );
/// variable flag
VARIABLE_FLAG  intToVARIABLE_FLAG( int i );
VARIABLE_FLAG  parseCondition( std::string& s );
VARIABLE_FLAG  parseStatus( const char* status );
std::string    parseStatus( VARIABLE_FLAG );

/// returns the side of the Face or InterFace element
std::string    parseSide( INTERFACE_SIDE );

/** @} */


} // end namespace csmp

#endif
