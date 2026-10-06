#include "CSMP_global_enumerations.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include "Region.h"
#include "Boundary.h"
#include "Edge.h"
#include "SplitBoundary.h"
#include "Exception.h"

using namespace std;

namespace csmp {

  void variablePlacementSet( set<PLACEMENT>& plSet )
  {
    plSet.insert(NODE);
    plSet.insert(FACE);
    plSet.insert(FACE_INTEGRATION_POINT);
    plSet.insert(FACE_SECTOR_INTEGRATION_POINT);
    plSet.insert(FACE_FACET_INTEGRATION_POINT);
    plSet.insert(INTER_FACE);
    plSet.insert(INTER_FACE_INTEGRATION_POINT);
    plSet.insert(INTER_FACE_SECTOR_INTEGRATION_POINT);
    plSet.insert(INTER_FACE_FACET_INTEGRATION_POINT);
    plSet.insert(ELEMENT);
    plSet.insert(ELEMENT_INTEGRATION_POINT);
    plSet.insert(SECTOR_INTEGRATION_POINT);
    plSet.insert(FACET_INTEGRATION_POINT);
    plSet.insert(REGION);
    plSet.insert(BOUNDARY);
    plSet.insert(SPLIT_BOUNDARY);
    plSet.insert(MODEL);
  }

  void fvVariablePlacementSet( set<PLACEMENT>& plSet )
  {
    plSet.insert(FACE_SECTOR_INTEGRATION_POINT);
    plSet.insert(FACE_FACET_INTEGRATION_POINT);
    plSet.insert(INTER_FACE_SECTOR_INTEGRATION_POINT);
    plSet.insert(INTER_FACE_FACET_INTEGRATION_POINT);
    plSet.insert(SECTOR_INTEGRATION_POINT);
    plSet.insert(FACET_INTEGRATION_POINT);
  }



  void variableTypeSet( set<VARIABLE_TYPE>& vtSet )
  {
    vtSet.insert(SCALAR);
    vtSet.insert(VECTOR);
    vtSet.insert(TENSOR);
    vtSet.insert(ARRAY);
    vtSet.insert(FLAGGEDARRAY);
  }

// VARIABLE TYPE

VARIABLE_TYPE  parseType( int type )
 {
    if ( type == SCALAR )        return SCALAR;
    if ( type == VECTOR )        return VECTOR;
    if ( type == TENSOR )        return TENSOR;
    if ( type == ARRAY )         return ARRAY;
    if ( type == FLAGGEDARRAY )  return FLAGGEDARRAY;

    cout <<"\nparseType(int): unable to parse: '";
    cout << type <<"' returning -1"<< endl;
    return static_cast<VARIABLE_TYPE>(-1);
 }

VARIABLE_TYPE  parseType( const char* type )
 {
    string  stype(type);
    if ( stype == "SCALAR" )        return SCALAR;
    if ( stype == "VECTOR" )        return VECTOR;
    if ( stype == "TENSOR" )        return TENSOR;
    if ( stype == "ARRAY" )         return ARRAY;
    if ( stype == "FLAGGEDARRAY" )  return FLAGGEDARRAY;

    if ( stype == "scalar" )        return SCALAR;
    if ( stype == "vector" )        return VECTOR;
    if ( stype == "tensor" )        return TENSOR;
    if ( stype == "array" )         return ARRAY;
    if ( stype == "flaggedarray" )  return FLAGGEDARRAY;

    cout <<"\nparseType(const char*): unable to parse: '";
    cout << type <<"' returning -1"<< endl;
    return static_cast<VARIABLE_TYPE>(-1);
 }

std::string  parseType( VARIABLE_TYPE vtype )
 {
    if ( vtype == SCALAR )          return string("SCALAR");
    if ( vtype == VECTOR )          return string("VECTOR");
    if ( vtype == TENSOR )          return string("TENSOR");
    if ( vtype == ARRAY )           return string("ARRAY");
    if ( vtype == FLAGGEDARRAY )    return string("FLAGGEDARRAY");

    cout <<"\nparseType(VARIABLE_TYPE): unable to parse VARIABLE_TYPE: "<< vtype << endl;
    return string("NOT_SPECIFIED");
 }

template<typename csmp_type>
VARIABLE_TYPE variableType()
{
    return SCALAR;
}

template<>
VARIABLE_TYPE variableType<csmp::ScalarVariable>()
{
    return SCALAR;
}

template<>
VARIABLE_TYPE variableType<csmp::VectorVariable<1U> >()
{
    return VECTOR;
}

template<>
VARIABLE_TYPE variableType<csmp::VectorVariable<2U> >()
{
    return VECTOR;
}

template<>
VARIABLE_TYPE variableType<csmp::VectorVariable<3U> >()
{
    return VECTOR;
}

template<>
VARIABLE_TYPE variableType<csmp::TensorVariable<1U> >()
{
    return TENSOR;
}

template<>
VARIABLE_TYPE variableType<csmp::TensorVariable<2U> >()
{
    return TENSOR;
}

template<>
VARIABLE_TYPE variableType<csmp::TensorVariable<3U> >()
{
    return TENSOR;
}

template<>
VARIABLE_TYPE variableType<csmp::ArrayVariable>()
{
    return ARRAY;
}

template<>
VARIABLE_TYPE variableType<csmp::FlaggedArrayVariable>()
{
    return FLAGGEDARRAY;
}

// VARIABLE PLACEMENT

template<>  PLACEMENT parsePlacement<1U,Node>() noexcept { return NODE; }
template<>  PLACEMENT parsePlacement<2U,Node>() noexcept { return NODE; }
template<>  PLACEMENT parsePlacement<3U,Node>() noexcept { return NODE; }

template<>  PLACEMENT parsePlacement<1U,Element>() noexcept { return ELEMENT; }
template<>  PLACEMENT parsePlacement<2U,Element>() noexcept { return ELEMENT; }
template<>  PLACEMENT parsePlacement<3U,Element>() noexcept { return ELEMENT; }

template<>  PLACEMENT parsePlacement<1U,Face>() noexcept { return FACE; }
template<>  PLACEMENT parsePlacement<2U,Face>() noexcept { return FACE; }
template<>  PLACEMENT parsePlacement<3U,Face>() noexcept { return FACE; }

template<>  PLACEMENT parsePlacement<1U,Edge>() noexcept { return FACE; }
template<>  PLACEMENT parsePlacement<2U,Edge>() noexcept { return FACE; }
template<>  PLACEMENT parsePlacement<3U,Edge>() noexcept { return FACE; }

template<>  PLACEMENT parsePlacement<1U,InterFace>() noexcept { return INTER_FACE; }
template<>  PLACEMENT parsePlacement<2U,InterFace>() noexcept { return INTER_FACE; }
template<>  PLACEMENT parsePlacement<3U,InterFace>() noexcept { return INTER_FACE; }

template<>  PLACEMENT parsePlacement<1U,Region>() noexcept { return REGION; }
template<>  PLACEMENT parsePlacement<2U,Region>() noexcept { return REGION; }
template<>  PLACEMENT parsePlacement<3U,Region>() noexcept { return REGION; }

template<>  PLACEMENT parsePlacement<1U,Boundary>() noexcept { return BOUNDARY; }
template<>  PLACEMENT parsePlacement<2U,Boundary>() noexcept { return BOUNDARY; }
template<>  PLACEMENT parsePlacement<3U,Boundary>() noexcept { return BOUNDARY; }

template<>  PLACEMENT parsePlacement<1U,SplitBoundary>() noexcept { return SPLIT_BOUNDARY; }
template<>  PLACEMENT parsePlacement<2U,SplitBoundary>() noexcept { return SPLIT_BOUNDARY; }
template<>  PLACEMENT parsePlacement<3U,SplitBoundary>() noexcept { return SPLIT_BOUNDARY; }

template<>  PLACEMENT parsePlacement<1U,Model>() noexcept { return MODEL; }
template<>  PLACEMENT parsePlacement<2U,Model>() noexcept { return MODEL; }
template<>  PLACEMENT parsePlacement<3U,Model>() noexcept { return MODEL; }



PLACEMENT parsePlacement( const char* placement ) noexcept
{
    if ( !placement ) return MODEL;

    // normalise: lowercase + spaces → underscores
    // collapses all variants into a single canonical form:
    //   "Element Integration Point" → "element_integration_point"
    //   "ELEMENT_INTEGRATION_POINT" → "element_integration_point"
    //   "element integration point" → "element_integration_point"
    //   "InterFace"                 → "interface"
    std::string n( placement );
    std::transform( n.begin(), n.end(), n.begin(),
                    []( unsigned char c ) {
                        return c == ' ' ? '_' : static_cast<char>( std::tolower(c) );
                    } );

    using sv = std::string_view;

    constexpr std::array<std::pair<sv, PLACEMENT>, 38> lookup {{
        { sv{"node"},                                NODE                                },
        { sv{"element"},                             ELEMENT                             },
        { sv{"elmt"},                                ELEMENT                             },
        { sv{"ipoint"},                              ELEMENT_INTEGRATION_POINT           },
        { sv{"cpoint"},                              ELEMENT_INTEGRATION_POINT           },
        { sv{"integration_point"},                   ELEMENT_INTEGRATION_POINT           },
        { sv{"integrationpoint"},                    ELEMENT_INTEGRATION_POINT           },
        { sv{"element_integration_point"},           ELEMENT_INTEGRATION_POINT           },
        { sv{"simplex_integration_point"},           ELEMENT_INTEGRATION_POINT           },
        { sv{"simplexintegrationpoint"},             ELEMENT_INTEGRATION_POINT           },
        { sv{"sector_integration_point"},            SECTOR_INTEGRATION_POINT            },
        { sv{"sectorintegrationpoint"},              SECTOR_INTEGRATION_POINT            },
        { sv{"facet_integration_point"},             FACET_INTEGRATION_POINT             },
        { sv{"facetintegrationpoint"},               FACET_INTEGRATION_POINT             },
        { sv{"face"},                                FACE                                },
        { sv{"face_integration_point"},              FACE_INTEGRATION_POINT              },
        { sv{"faceintegrationpoint"},                FACE_INTEGRATION_POINT              },
        { sv{"face_sector_integration_point"},       FACE_SECTOR_INTEGRATION_POINT       },
        { sv{"facesectorintegrationpoint"},          FACE_SECTOR_INTEGRATION_POINT       },
        { sv{"face_facet_integration_point"},        FACE_FACET_INTEGRATION_POINT        },
        { sv{"facefacetintegrationpoint"},           FACE_FACET_INTEGRATION_POINT        },
        { sv{"interface"},                           INTER_FACE                          },
        { sv{"inter_face"},                          INTER_FACE                          },
        { sv{"iface"},                               INTER_FACE                          },
        { sv{"inter_face_integration_point"},        INTER_FACE_INTEGRATION_POINT        },
        { sv{"interface_integration_point"},         INTER_FACE_INTEGRATION_POINT        },
        { sv{"interfaceintegrationpoint"},           INTER_FACE_INTEGRATION_POINT        },
        { sv{"inter_face_sector_integration_point"}, INTER_FACE_SECTOR_INTEGRATION_POINT },
        { sv{"interface_sector_integration_point"},  INTER_FACE_SECTOR_INTEGRATION_POINT },
        { sv{"interfacesectorintegrationpoint"},     INTER_FACE_SECTOR_INTEGRATION_POINT },
        { sv{"inter_face_facet_integration_point"},  INTER_FACE_FACET_INTEGRATION_POINT  },
        { sv{"interface_facet_integration_point"},   INTER_FACE_FACET_INTEGRATION_POINT  },
        { sv{"interfacefacetintegrationpoint"},      INTER_FACE_FACET_INTEGRATION_POINT  },
        { sv{"region"},                              REGION                              },
        { sv{"boundary"},                            BOUNDARY                            },
        { sv{"split_boundary"},                      SPLIT_BOUNDARY                      },
        { sv{"splitboundary"},                       SPLIT_BOUNDARY                      },
        { sv{"model"},                               MODEL                               },
    }};

    const sv ns( n );
    for ( const auto& p : lookup )
        if ( p.first == ns ) return p.second;

    std::cerr << "\nparsePlacement: unable to parse: '"
              << placement << "' — returning MODEL\n";
    return MODEL;
}




/**
    @brief Converts a PLACEMENT enum value to its string representation.

    @param placement  The PLACEMENT value to convert.
    @return           String name of the enum value, or "MODEL" if not recognised
                      (with a diagnostic message to stderr).
*/
std::string parsePlacement( PLACEMENT placement ) noexcept
{
    switch ( placement ) {
        case NODE:                                return "NODE";
        case ELEMENT:                             return "ELEMENT";
        case ELEMENT_INTEGRATION_POINT:           return "ELEMENT_INTEGRATION_POINT";
        case SECTOR_INTEGRATION_POINT:            return "SECTOR_INTEGRATION_POINT";
        case FACET_INTEGRATION_POINT:             return "FACET_INTEGRATION_POINT";
        case FACE:                                return "FACE";
        case FACE_INTEGRATION_POINT:              return "FACE_INTEGRATION_POINT";
        case FACE_SECTOR_INTEGRATION_POINT:       return "FACE_SECTOR_INTEGRATION_POINT";
        case FACE_FACET_INTEGRATION_POINT:        return "FACE_FACET_INTEGRATION_POINT";
        case INTER_FACE:                          return "INTER_FACE";
        case INTER_FACE_INTEGRATION_POINT:        return "INTER_FACE_INTEGRATION_POINT";
        case INTER_FACE_SECTOR_INTEGRATION_POINT: return "INTER_FACE_SECTOR_INTEGRATION_POINT";
        case INTER_FACE_FACET_INTEGRATION_POINT:  return "INTER_FACE_FACET_INTEGRATION_POINT";
        case REGION:                              return "REGION";
        case BOUNDARY:                            return "BOUNDARY";
        case SPLIT_BOUNDARY:                      return "SPLIT_BOUNDARY";
        case MODEL:                               return "MODEL";
        default:
            std::cerr << "\nparsePlacement(PLACEMENT): unable to parse: "
                      << static_cast<int>( placement ) << '\n';
            return "MODEL";
    }
}

 
 
 
 
bool faceVariable(PLACEMENT place) noexcept {
    static constexpr std::array<PLACEMENT, 4> faceTypes{
        FACE,
        FACE_INTEGRATION_POINT,
        FACE_FACET_INTEGRATION_POINT,
        FACE_SECTOR_INTEGRATION_POINT
    };
    
    return std::ranges::contains(faceTypes, place);
}


bool interFaceVariable(PLACEMENT place) noexcept {
    static constexpr std::array<PLACEMENT, 4> interFaceTypes{
        INTER_FACE,
        INTER_FACE_INTEGRATION_POINT,
        INTER_FACE_FACET_INTEGRATION_POINT,
        INTER_FACE_SECTOR_INTEGRATION_POINT
    };
    
    return std::ranges::contains(interFaceTypes, place);
}


bool isPlacedOnIntegrationPoint(PLACEMENT p) noexcept {
    static constexpr std::array<PLACEMENT, 9> integrationPoints{
        ELEMENT_INTEGRATION_POINT,
        SECTOR_INTEGRATION_POINT,
        FACET_INTEGRATION_POINT,
        FACE_INTEGRATION_POINT,
        FACE_FACET_INTEGRATION_POINT,
        FACE_SECTOR_INTEGRATION_POINT,
        INTER_FACE_INTEGRATION_POINT,
        INTER_FACE_SECTOR_INTEGRATION_POINT,
        INTER_FACE_FACET_INTEGRATION_POINT
    };
    
    return std::ranges::contains(integrationPoints, p);
}
 
 

/**
    @brief Converts an integer to the corresponding PLACEMENT enum value.

    The integer encoding matches the internal CSMP binary storage format.
    Note that value 4 is not assigned — it is a legacy gap in the encoding.

    @param i  Integer to convert.
    @return   Corresponding PLACEMENT value, or MODEL if not recognised
              (with a diagnostic message to stderr).
*/
PLACEMENT intToPLACEMENT( int i ) noexcept
{
    switch ( i ) {
        case  0: return NODE;
        case  1: return ELEMENT_INTEGRATION_POINT;
        case  2: return SECTOR_INTEGRATION_POINT;
        case  3: return FACET_INTEGRATION_POINT;
        // case 4: not assigned — legacy gap in encoding
        case  5: return FACE;
        case  6: return INTER_FACE;
        case  7: return ELEMENT;
        case  8: return REGION;
        case  9: return BOUNDARY;
        case 10: return SPLIT_BOUNDARY;
        case 11: return MODEL;
        case 12: return FACE_INTEGRATION_POINT;
        case 13: return FACE_SECTOR_INTEGRATION_POINT;
        case 14: return FACE_FACET_INTEGRATION_POINT;
        case 15: return INTER_FACE_INTEGRATION_POINT;
        case 16: return INTER_FACE_SECTOR_INTEGRATION_POINT;
        case 17: return INTER_FACE_FACET_INTEGRATION_POINT;
        default:
            std::cerr << "\nintToPlacement: unable to parse integer: "
                      << i << '\n';
            return MODEL;
    }
}




// VARIABLE FLAG

/**
    @brief Converts an integer to the corresponding VARIABLE_FLAG enum value.

    The integer encoding matches the internal CSMP binary storage format.

    @param i  Integer to convert.
    @return   Corresponding VARIABLE_FLAG value, or ANY if not recognised
              (with a diagnostic message to stderr).
*/
VARIABLE_FLAG intToVARIABLE_FLAG( int i ) noexcept
{
    switch ( i ) {
        case 0: return PLAIN;
        case 1: return ANY;
        case 2: return INIT_GUESS;
        case 3: return INIT_COND;
        case 4: return FIELD_DATA;
        case 5: return PERIODIC;
        case 6: return DIRICH;
        case 7: return NEUMANN;
        case 8: return ROBIN;
        case 9: return CONSTANT_FLUX;
        default:
            std::cerr << "\nintToVARIABLE_FLAG: unable to parse integer: "
                      << i << '\n';
            return ANY;
    }
}



VARIABLE_FLAG parseStatus( const char* status ) noexcept
{
    using sv = std::string_view;
    constexpr std::array lookup {
        std::pair{sv{"ANY"},               ANY        },
        std::pair{sv{"any"},               ANY        },
        std::pair{sv{"PLAIN"},             PLAIN      },
        std::pair{sv{"Plain"},             PLAIN      },
        std::pair{sv{"plain"},             PLAIN      },
        std::pair{sv{"DIRICH"},            DIRICH     },
        std::pair{sv{"DIRICHLET"},         DIRICH     },
        std::pair{sv{"Dirichlet"},         DIRICH     },
        std::pair{sv{"NEUMANN"},           NEUMANN    },
        std::pair{sv{"Neumann"},           NEUMANN    },
        std::pair{sv{"PERIODIC"},          PERIODIC   },
        std::pair{sv{"periodic"},          PERIODIC   },
        std::pair{sv{"ROBIN"},             ROBIN      },
        std::pair{sv{"Robin"},             ROBIN      },
        std::pair{sv{"FIELD_DATA"},        FIELD_DATA },
        std::pair{sv{"field_data"},        FIELD_DATA },
        std::pair{sv{"INIT_GUESS"},        INIT_GUESS },
        std::pair{sv{"initial_guess"},     INIT_GUESS },
        std::pair{sv{"INIT_COND"},         INIT_COND  },
        std::pair{sv{"initial_condition"}, INIT_COND  },
    };

    if ( !status ) return ANY;
    const sv s( status );
    for ( const auto& p : lookup )
        if ( p.first == s ) return p.second;

    std::cerr << "\nparseStatus: unable to parse: '" << status << "'\n";
    return ANY;
}


std::string parseStatus(VARIABLE_FLAG status) noexcept {
    static const std::unordered_map<VARIABLE_FLAG, std::string> reverseStatusMap{
        {PLAIN, "PLAIN"},
        {DIRICH, "DIRICH"},
        {NEUMANN, "NEUMANN"},
        {INIT_GUESS, "INIT_GUESS"},
        {INIT_COND, "INIT_COND"},
        {ROBIN, "ROBIN"},
        {FIELD_DATA, "FIELD_DATA"},
        {PERIODIC, "PERIODIC"},
        {ANY, "ANY"},
    };

    if (auto it = reverseStatusMap.find(status); it != reverseStatusMap.end()) {
        return it->second;
    }

    std::cerr << "\nparseStatus(VARIABLE_FLAG): unable to parse: " 
              << static_cast<int>(status) <<"'\n";
 
    return "ANY";
 }


/**
    @brief Interprets a string as a CSMP variable boundary condition flag.

    Input is normalised to uppercase before matching so the comparison
    is case-insensitive. Both underscore and space-separated variants
    are accepted (e.g. "INITIAL_CONDITION" and "INITIAL CONDITION").

    @param s  String to interpret (e.g. "DIRICH", "Neumann", "plain").
    @return   Corresponding VARIABLE_FLAG value.
    @throws   csmp::Exception (FATAL_ERROR) if the string is not recognised.
*/
VARIABLE_FLAG parseCondition( std::string_view s )
{
    using sv = std::string_view;

    // compile-time lookup table — no runtime initialisation, no hashing, no UB
    constexpr std::array<std::pair<std::string_view, VARIABLE_FLAG>, 20> lookup {{
        std::pair{sv{"ANY"},               ANY           },
        std::pair{sv{"PLAIN"},             PLAIN         },
        std::pair{sv{"DIRICH"},            DIRICH        },
        std::pair{sv{"DIRICHLET"},         DIRICH        },
        std::pair{sv{"NEUMANN"},           NEUMANN       },
        std::pair{sv{"NATURAL"},           NEUMANN       },
        std::pair{sv{"ROBIN"},             ROBIN         },
        std::pair{sv{"PERIODIC"},          PERIODIC      },
        std::pair{sv{"INIT_GUESS"},        INIT_GUESS    },
        std::pair{sv{"INITIAL_GUESS"},     INIT_GUESS    },
        std::pair{sv{"init_guess"},        INIT_GUESS    },
        std::pair{sv{"INIT_COND"},         INIT_COND     },
        std::pair{sv{"INITIAL_CONDITION"}, INIT_COND     },
        std::pair{sv{"init_cond"},         INIT_COND     },
        std::pair{sv{"FIELD_DATA"},        FIELD_DATA    },
        std::pair{sv{"field_data"},        FIELD_DATA    },
        std::pair{sv{"CONST_FLUX"},        CONSTANT_FLUX },
        std::pair{sv{"CONSTANT_FLUX"},     CONSTANT_FLUX },
        std::pair{sv{"const_flux"},        CONSTANT_FLUX },
    }};

    // normalise input to uppercase
    std::string n( s );
    std::ranges::transform( n, n.begin(),
                            []( unsigned char c ) { return std::toupper(c); } );

    const sv normalised( n );
    for ( const auto& p : lookup )
        if ( p.first == normalised ) return p.second;

    throw csmp::Exception( FATAL_ERROR,
                           "parseCondition:",
                           n.c_str(),
                           ": boundary condition type specifier was not recognised" );
}




/**
    @brief Converts an INTERFACE_SIDE enum value to its string representation.

    @param side  The INTERFACE_SIDE value to convert.
    @return      String name of the enum value, or "UNDEFINED" if not recognised.
*/
std::string parseSide( INTERFACE_SIDE side ) noexcept
{
    switch ( side ) {
        case INSIDE:  return "INSIDE";
        case OUTSIDE: return "OUTSIDE";
        case MIDDLE:  return "MIDDLE";
        default:      return "UNDEFINED";
    }
}




/**
    STAND_ALONE,                ///< multiplicated point
    SPLIT_BOUNDARY,             ///<  2-node manifold along a SplitBoundary (most common)
    SPLIT_BOUNDARY_CROSSING,    ///<  4-node manifold intersection of split boundaries
    MULTI_SB_CROSSING,          ///<  6-node cross of 3 SBs in
    SPLIT_BOUNDARY_TERMINATION, ///<  T-intersection of SBs are termination of SB against Boundary
    SPLIT_BOUNDARY_END,         ///<  termination against model boundary
    
    @note if the split boundary cannot be determined, method returns SPLIT_BOUNDARY
*/
string parse( ManifoldType topology )
 {
   switch( topology ) {
        case ManifoldType::STAND_ALONE: return "STAND_ALONE";
        case ManifoldType::SPLIT_BOUNDARY: return "SPLIT_BOUNDARY";
        case ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH: return "SPLIT_BOUNDARY_WITH_INTERNAL_MESH";
        case ManifoldType::SPLIT_BOUNDARY_CROSSING: return "SPLIT_BOUNDARY_CROSSING";
        case ManifoldType::MULTI_SB_CROSSING: return "MULTI_SB_CROSSING";
//        case ManifoldType::SPLIT_BOUNDARY_TERMINATION: return "SPLIT_BOUNDARY_TERMINATION";
        case ManifoldType::SPLIT_BOUNDARY_END:  return "SPLIT_BOUNDARY_END";
        //default: return "NOT_CLASSIFIED";
      }
    return "SPLIT_BOUNDARY";
 }


/**
       Topologic qualifiers include:
       
      MESH_VERTEX,        - a point within the model volume
      INTERSECTION_POINT  - a point where lines cross or multiple surfaces intersect
      PERIMETER_POINT,    - point at the end of a line inside a 2D model
      EXTERIOR_POINT,     - on an outside surface of the model
      INTERIOR_LINE,      - a line on the interior of the model
      PERIMETER_LINE,     - a surface edge inside of the model
      EXTERIOR_LINE,      - an edge of the model
      INTERSECTION_LINE,  - belonging to multiple surfaces
      INTERIOR_SURFACE,   - a surface withing the model
      PERIMETER_SURFACE,  - a surface forming the hull of an object inside of the model
      EXTERIOR_SURFACE    - a surface delimiting the model
      
*/
std::string parseTopology( TOPOTYPE topology )
 {
   switch( topology ) {
        case MESH_VERTEX: return "MESH_VERTEX";
        case INTERIOR_POINT: return "INTERIOR_POINT";
        case PERIMETER_POINT: return "PERIMETER_POINT";
        case EXTERIOR_POINT: return "EXTERIOR_POINT";
        case INTERIOR_LINE: return "INTERIOR_LINE";
        case PERIMETER_LINE: return "PERIMETER_LINE";
        case EXTERIOR_LINE:  return "EXTERIOR_LINE";
        case INTERIOR_SURFACE: return "INTERIOR_SURFACE";
        case PERIMETER_SURFACE: return "PERIMETER_SURFACE";
        case EXTERIOR_SURFACE: return "EXTERIOR_SURFACE";
        default:
//          cerr <<"\n"<<"parseTopology: error: topology not parsed, returning 'NOT CLASSIFIED'"<< endl;
          return "NOT_CLASSIFIED";
      }
 }




} // end csmp
