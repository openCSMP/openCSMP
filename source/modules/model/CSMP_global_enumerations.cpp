// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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

// Helper: normalise string (lowercase)
static std::string normalise(std::string_view str) noexcept {
    std::string result(str);
    std::ranges::transform(result, result.begin(), 
                          [](unsigned char c) { return std::tolower(c); });
    return result;
}
// Static map for placement parsing
PLACEMENT parsePlacement(const char* placement) noexcept {
    static const std::unordered_map<std::string, PLACEMENT> placementMap{
        // Uppercase variants
        {"NODE", NODE},
        {"ELEMENT_INTEGRATION_POINT", ELEMENT_INTEGRATION_POINT},
        {"INTEGRATION_POINT", ELEMENT_INTEGRATION_POINT},
        {"SECTOR_INTEGRATION_POINT", SECTOR_INTEGRATION_POINT},
        {"FACET_INTEGRATION_POINT", FACET_INTEGRATION_POINT},
        {"FACE_INTEGRATION_POINT", FACE_INTEGRATION_POINT},
        {"FACE_SECTOR_INTEGRATION_POINT", FACE_SECTOR_INTEGRATION_POINT},
        {"FACE_FACET_INTEGRATION_POINT", FACE_FACET_INTEGRATION_POINT},
        {"INTER_FACE_INTEGRATION_POINT", INTER_FACE_INTEGRATION_POINT},
        {"INTER_FACE_SECTOR_INTEGRATION_POINT", INTER_FACE_SECTOR_INTEGRATION_POINT},
        {"INTER_FACE_FACET_INTEGRATION_POINT", INTER_FACE_FACET_INTEGRATION_POINT},
        {"FACE", FACE},
        {"INTER_FACE", INTER_FACE},
        {"INTERFACE", INTER_FACE},
        {"ELEMENT", ELEMENT},
        {"REGION", REGION},
        {"BOUNDARY", BOUNDARY},
        {"SPLIT_BOUNDARY", SPLIT_BOUNDARY},
        {"MODEL", MODEL},
        
        // Lowercase variants
        {"node", NODE},
        {"cpoint", ELEMENT_INTEGRATION_POINT},
        {"ipoint", ELEMENT_INTEGRATION_POINT},
        {"integration point", ELEMENT_INTEGRATION_POINT},
        {"integration_point", ELEMENT_INTEGRATION_POINT},
        {"simplex integration point", ELEMENT_INTEGRATION_POINT},
        {"simplex_integration_point", ELEMENT_INTEGRATION_POINT},
        {"element integration point", ELEMENT_INTEGRATION_POINT},
        {"element_integration_point", ELEMENT_INTEGRATION_POINT},
        {"sector integration point", SECTOR_INTEGRATION_POINT},
        {"sector_integration_point", SECTOR_INTEGRATION_POINT},
        {"facet integration point", FACET_INTEGRATION_POINT},
        {"facet_integration_point", FACET_INTEGRATION_POINT},
        {"face integration point", FACE_INTEGRATION_POINT},
        {"face_integration_point", FACE_INTEGRATION_POINT},
        {"face sector integration point", FACE_SECTOR_INTEGRATION_POINT},
        {"face_sector_integration_point", FACE_SECTOR_INTEGRATION_POINT},
        {"face facet integration point", FACE_FACET_INTEGRATION_POINT},
        {"face_facet_integration_point", FACE_FACET_INTEGRATION_POINT},
        {"interface integration point", INTER_FACE_INTEGRATION_POINT},
        {"inter face integration point", INTER_FACE_INTEGRATION_POINT},
        {"inter_face_integration_point", INTER_FACE_INTEGRATION_POINT},
        {"interface sector integration point", INTER_FACE_SECTOR_INTEGRATION_POINT},
        {"inter face sector integration point", INTER_FACE_SECTOR_INTEGRATION_POINT},
        {"inter_face_sector_integration_point", INTER_FACE_SECTOR_INTEGRATION_POINT},
        {"inter face facet integration point", INTER_FACE_FACET_INTEGRATION_POINT},
        {"interface facet integration point", INTER_FACE_FACET_INTEGRATION_POINT},
        {"inter_face_facet_integration_point", INTER_FACE_FACET_INTEGRATION_POINT},
        {"face", FACE},
        {"interface", INTER_FACE},
        {"inter face", INTER_FACE},
        {"inter_face", INTER_FACE},
        {"iface", INTER_FACE},
        {"element", ELEMENT},
        {"elmt", ELEMENT},
        {"region", REGION},
        {"boundary", BOUNDARY},
        {"splitboundary", SPLIT_BOUNDARY},
        {"model", MODEL},
        
        // PascalCase variants
        {"Node", NODE},
        {"IntegrationPoint", ELEMENT_INTEGRATION_POINT},
        {"SimplexIntegrationPoint", ELEMENT_INTEGRATION_POINT},
        {"FacetIntegrationPoint", FACET_INTEGRATION_POINT},
        {"SectorIntegrationPoint", SECTOR_INTEGRATION_POINT},
        {"FaceIntegrationPoint", FACE_INTEGRATION_POINT},
        {"FaceFacetIntegrationPoint", FACE_FACET_INTEGRATION_POINT},
        {"FaceSectorIntegrationPoint", FACE_SECTOR_INTEGRATION_POINT},
        {"InterFaceIntegrationPoint", INTER_FACE_INTEGRATION_POINT},
        {"InterFaceFacetIntegrationPoint", INTER_FACE_FACET_INTEGRATION_POINT},
        {"InterFaceSectorIntegrationPoint", INTER_FACE_SECTOR_INTEGRATION_POINT},
        {"Face", FACE},
        {"InterFace", INTER_FACE},
        {"Element", ELEMENT},
        {"Region", REGION},
        {"Boundary", BOUNDARY},
        {"SplitBoundary", SPLIT_BOUNDARY},
        {"Model", MODEL},
    };
    if (auto it = placementMap.find(placement); it != placementMap.end()) {
        return it->second;
    }
    std::cerr << "\nparsePlacement(const char*): unable to parse: '" 
              << placement << "' returning MODEL\n";
    return MODEL;
}


std::string parsePlacement(PLACEMENT placement) noexcept {
    static const std::unordered_map<PLACEMENT, std::string> reverseMap{
        {NODE, "NODE"},
        {ELEMENT_INTEGRATION_POINT, "ELEMENT_INTEGRATION_POINT"},
        {FACET_INTEGRATION_POINT, "FACET_INTEGRATION_POINT"},
        {SECTOR_INTEGRATION_POINT, "SECTOR_INTEGRATION_POINT"},
        {FACE_INTEGRATION_POINT, "FACE_INTEGRATION_POINT"},
        {FACE_FACET_INTEGRATION_POINT, "FACE_FACET_INTEGRATION_POINT"},
        {FACE_SECTOR_INTEGRATION_POINT, "FACE_SECTOR_INTEGRATION_POINT"},
        {INTER_FACE_INTEGRATION_POINT, "INTER_FACE_INTEGRATION_POINT"},
        {INTER_FACE_FACET_INTEGRATION_POINT, "INTER_FACE_FACET_INTEGRATION_POINT"},
        {INTER_FACE_SECTOR_INTEGRATION_POINT, "INTER_FACE_SECTOR_INTEGRATION_POINT"},
        {FACE, "FACE"},
        {INTER_FACE, "INTER_FACE"},
        {ELEMENT, "ELEMENT"},
        {REGION, "REGION"},
        {BOUNDARY, "BOUNDARY"},
        {SPLIT_BOUNDARY, "SPLIT_BOUNDARY"},
        {MODEL, "MODEL"},
    };

    if (auto it = reverseMap.find(placement); it != reverseMap.end()) {
        return it->second;
    }

    std::cerr << "\nparsePlacement(PLACEMENT): unable to parse: " 
              << static_cast<int>(placement) << '\n';
    return "MODEL";
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
 
 

PLACEMENT intToPLACEMENT(int i) noexcept {
    static const std::unordered_map<int, PLACEMENT> intMap{
        {0, NODE},
        {1, ELEMENT_INTEGRATION_POINT},
        {2, SECTOR_INTEGRATION_POINT},
        {3, FACET_INTEGRATION_POINT},
        {5, FACE},
        {6, INTER_FACE},
        {7, ELEMENT},
        {8, REGION},
        {9, BOUNDARY},
        {10, SPLIT_BOUNDARY},
        {11, MODEL},
        {12, FACE_INTEGRATION_POINT},
        {13, FACE_SECTOR_INTEGRATION_POINT},
        {14, FACE_FACET_INTEGRATION_POINT},
        {15, INTER_FACE_INTEGRATION_POINT},
        {16, INTER_FACE_SECTOR_INTEGRATION_POINT},
        {17, INTER_FACE_FACET_INTEGRATION_POINT},
    };

    if (auto it = intMap.find(i); it != intMap.end()) {
        return it->second;
    }

    std::cerr << "\nintToPlacement(int): unable to parse integer: " << i << '\n';
    return MODEL;
}



// VARIABLE FLAG

VARIABLE_FLAG intToVARIABLE_FLAG(int i) noexcept {
    static const std::unordered_map<int, VARIABLE_FLAG> flagMap{
        {0, PLAIN},
        {1, ANY},
        {2, INIT_GUESS},
        {3, INIT_COND},
        {4, FIELD_DATA},
        {5, PERIODIC},
        {6, DIRICH},
        {7, NEUMANN},
        {8, ROBIN},
        {9, CONSTANT_FLUX},
    };

    if (auto it = flagMap.find(i); it != flagMap.end()) {
        return it->second;
    }

    std::cerr << "\nintToVARIABLE_FLAG(int): unable to parse integer: " << i << '\n';
    return ANY;
}


VARIABLE_FLAG parseStatus(const char* status) noexcept {
    static const std::unordered_map<std::string, VARIABLE_FLAG> statusMap{
        {"ANY", ANY},
        {"any", ANY},
        {"PLAIN", PLAIN},
        {"Plain", PLAIN},
        {"plain", PLAIN},
        {"DIRICH", DIRICH},
        {"DIRICHLET", DIRICH},
        {"Dirichlet", DIRICH},
        {"NEUMANN", NEUMANN},
        {"Neumann", NEUMANN},
        {"PERIODIC", PERIODIC},
        {"periodic", PERIODIC},
        {"ROBIN", ROBIN},
        {"Robin", ROBIN},
        {"FIELD_DATA", FIELD_DATA},
        {"field data", FIELD_DATA},
        {"INIT_GUESS", INIT_GUESS},
        {"initial guess", INIT_GUESS},
        {"INIT_COND", INIT_COND},
        {"initial condition", INIT_COND},
    };

    if (auto it = statusMap.find(status); it != statusMap.end()) {
        return it->second;
    }

    std::cerr << "\nparseStatus(const char*): unable to parse: '" << status << "'\n";
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

Interprets the supplied string as an essential (boundary) condition.
Possible values are given in 'CSMP_definitions.h'
and MULTIPLE. The last type can be used to specify Cauchy or
Robbin conditions.

@section arguments Input Arguments

The string that shall be interpreted as boundary condition.
*/
VARIABLE_FLAG parseCondition(std::string_view s) {
    static const std::unordered_map<std::string, VARIABLE_FLAG> conditionMap{
        {"PLAIN", PLAIN},
        {"ANY", ANY},
        {"INIT_GUESS", INIT_GUESS},
        {"INITIAL GUESS", INIT_GUESS},
        {"INIT_COND", INIT_COND},
        {"INITIAL CONDITION", INIT_COND},
        {"FIELD_DATA", FIELD_DATA},
        {"FIELD DATA", FIELD_DATA},
        {"PERIODIC", PERIODIC},
        {"DIRICH", DIRICH},
        {"DIRICHLET", DIRICH},
        {"NEUMANN", NEUMANN},
        {"NATURAL", NEUMANN},
        {"ROBIN", ROBIN},
        {"CONST_FLUX", CONSTANT_FLUX},
        {"CONSTANT_FLUX", CONSTANT_FLUX},
        {"CONSTANT FLUX", CONSTANT_FLUX},
    };

    // Normalise input to uppercase
    std::string normalized(s);
    std::ranges::transform(normalized, normalized.begin(),
                          [](unsigned char c) { return std::toupper(c); });

    if (auto it = conditionMap.find(normalized); it != conditionMap.end()) {
        return it->second;
    }

    throw csmp::Exception(FATAL_ERROR,
                         "parseCondition:",
                         normalized.c_str(),
                         ": boundary condition type specifier was not recognised");
}



/**
    To print enum INTERFACE_SIDE.
*/
std::string parseSide(INTERFACE_SIDE side) noexcept {
    static const std::unordered_map<INTERFACE_SIDE, std::string> sideMap{
        {INSIDE, "INSIDE"},
        {OUTSIDE, "OUTSIDE"},
        {MIDDLE, "MIDDLE"},
    };

    if (auto it = sideMap.find(side); it != sideMap.end()) {
        return it->second;
    }

    return "UNDEFINED";
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
