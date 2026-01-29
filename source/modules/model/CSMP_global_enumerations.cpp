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


template<uint32_t dim, template<uint32_t> class PLACE>
PLACEMENT parsePlacement() {
    return UNDEFINED;
 }

template<>  PLACEMENT parsePlacement<1U,Node>() { return NODE; }
template<>  PLACEMENT parsePlacement<2U,Node>() { return NODE; }
template<>  PLACEMENT parsePlacement<3U,Node>() { return NODE; }

template<>  PLACEMENT parsePlacement<1U,Element>() { return ELEMENT; }
template<>  PLACEMENT parsePlacement<2U,Element>() { return ELEMENT; }
template<>  PLACEMENT parsePlacement<3U,Element>() { return ELEMENT; }

template<>  PLACEMENT parsePlacement<1U,Face>() { return FACE; }
template<>  PLACEMENT parsePlacement<2U,Face>() { return FACE; }
template<>  PLACEMENT parsePlacement<3U,Face>() { return FACE; }

template<>  PLACEMENT parsePlacement<1U,Edge>() { return FACE; }
template<>  PLACEMENT parsePlacement<2U,Edge>() { return FACE; }
template<>  PLACEMENT parsePlacement<3U,Edge>() { return FACE; }

template<>  PLACEMENT parsePlacement<1U,InterFace>() { return INTER_FACE; }
template<>  PLACEMENT parsePlacement<2U,InterFace>() { return INTER_FACE; }
template<>  PLACEMENT parsePlacement<3U,InterFace>() { return INTER_FACE; }

template<>  PLACEMENT parsePlacement<1U,Region>() { return REGION; }
template<>  PLACEMENT parsePlacement<2U,Region>() { return REGION; }
template<>  PLACEMENT parsePlacement<3U,Region>() { return REGION; }

template<>  PLACEMENT parsePlacement<1U,Boundary>() { return BOUNDARY; }
template<>  PLACEMENT parsePlacement<2U,Boundary>() { return BOUNDARY; }
template<>  PLACEMENT parsePlacement<3U,Boundary>() { return BOUNDARY; }

template<>  PLACEMENT parsePlacement<1U,SplitBoundary>() { return SPLIT_BOUNDARY; }
template<>  PLACEMENT parsePlacement<2U,SplitBoundary>() { return SPLIT_BOUNDARY; }
template<>  PLACEMENT parsePlacement<3U,SplitBoundary>() { return SPLIT_BOUNDARY; }

template<>  PLACEMENT parsePlacement<1U,Model>() { return MODEL; }
template<>  PLACEMENT parsePlacement<2U,Model>() { return MODEL; }
template<>  PLACEMENT parsePlacement<3U,Model>() { return MODEL; }

PLACEMENT  parsePlacement( const char* placement )
 {
    string  splace(placement);
    if ( splace == "NODE"                                   )   return NODE;
    if ( splace == "ELEMENT_INTEGRATION_POINT"              )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "INTEGRATION_POINT"                      )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "ELEMENT_INTEGRATION_POINT"              )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "SECTOR_INTEGRATION_POINT"               )   return SECTOR_INTEGRATION_POINT;
    if ( splace == "FACET_INTEGRATION_POINT"                )   return FACET_INTEGRATION_POINT;
    if ( splace == "FACE_INTEGRATION_POINT"                 )   return FACE_INTEGRATION_POINT;
    if ( splace == "FACE_SECTOR_INTEGRATION_POINT"          )   return FACE_SECTOR_INTEGRATION_POINT;
    if ( splace == "FACE_FACET_INTEGRATION_POINT"           )   return FACE_FACET_INTEGRATION_POINT;
    if ( splace == "INTER_FACE_INTEGRATION_POINT"           )   return INTER_FACE_INTEGRATION_POINT;
    if ( splace == "INTER_FACE_SECTOR_INTEGRATION_POINT"    )   return INTER_FACE_SECTOR_INTEGRATION_POINT;
    if ( splace == "INTER_FACE_FACET_INTEGRATION_POINT"     )   return INTER_FACE_FACET_INTEGRATION_POINT;
    if ( splace == "FACE"                                   )   return FACE;
    if ( splace == "INTER_FACE"                             )   return INTER_FACE;
    if ( splace == "INTERFACE"                              )   return INTER_FACE;
    if ( splace == "ELEMENT"                                )   return ELEMENT;
    if ( splace == "REGION"                                 )   return REGION;
    if ( splace == "BOUNDARY"                               )   return BOUNDARY;
    if ( splace == "SPLIT_BOUNDARY"                         )   return SPLIT_BOUNDARY;
    if ( splace == "MODEL"                                  )   return MODEL;

    if ( splace == "node"                                   )   return NODE;
    if ( splace == "cpoint"                                 )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "ipoint"                                 )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "integration point"                      )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "integration_point"                      )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "simplex integration point"              )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "simplex_integration_point"              )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "element integration point"              )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "element_integration_point"              )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "sector integration point"               )   return SECTOR_INTEGRATION_POINT;
    if ( splace == "sector_integration_point"               )   return SECTOR_INTEGRATION_POINT;
    if ( splace == "facet integration point"                )   return FACET_INTEGRATION_POINT;
    if ( splace == "facet_integration_point"                )   return FACET_INTEGRATION_POINT;
    if ( splace == "face integration point"                 )   return FACE_INTEGRATION_POINT;
    if ( splace == "face_integration_point"                 )   return FACE_INTEGRATION_POINT;
    if ( splace == "face sector integration point"          )   return FACE_SECTOR_INTEGRATION_POINT;
    if ( splace == "face_sector_integration_point"          )   return FACE_SECTOR_INTEGRATION_POINT;
    if ( splace == "face facet integration point"           )   return FACE_FACET_INTEGRATION_POINT;
    if ( splace == "face_facet_integration_point"           )   return FACE_FACET_INTEGRATION_POINT;
    if ( splace == "facet integration point"                )   return FACET_INTEGRATION_POINT;
    if ( splace == "facet_integration_point"                )   return FACET_INTEGRATION_POINT;
    if ( splace == "interface integration point"            )   return INTER_FACE_INTEGRATION_POINT;
    if ( splace == "inter face integration point"           )   return INTER_FACE_INTEGRATION_POINT;
    if ( splace == "inter_face_integration_point"           )   return INTER_FACE_INTEGRATION_POINT;
    if ( splace == "interface sector integration point"     )   return INTER_FACE_SECTOR_INTEGRATION_POINT;
    if ( splace == "inter face sector integration point"    )   return INTER_FACE_SECTOR_INTEGRATION_POINT;
    if ( splace == "inter_face_sector_integration_point"    )   return INTER_FACE_SECTOR_INTEGRATION_POINT;
    if ( splace == "inter face facet integration point"     )   return INTER_FACE_FACET_INTEGRATION_POINT;
    if ( splace == "interface facet integration point"      )   return INTER_FACE_FACET_INTEGRATION_POINT;
    if ( splace == "inter_face_facet_integration_point"     )   return INTER_FACE_FACET_INTEGRATION_POINT;
    if ( splace == "face"                                   )   return FACE;
    if ( splace == "interface"                              )   return INTER_FACE;
    if ( splace == "inter face"                             )   return INTER_FACE;
    if ( splace == "inter_face"                             )   return INTER_FACE;
    if ( splace == "iface"                                  )   return INTER_FACE;
    if ( splace == "element"                                )   return ELEMENT;
    if ( splace == "elmt"                                   )   return ELEMENT;
    if ( splace == "region"                                 )   return REGION;
    if ( splace == "boundary"                               )   return BOUNDARY;
    if ( splace == "splitboundary"                          )   return SPLIT_BOUNDARY;
    if ( splace == "model"                                  )   return MODEL;

    if ( splace == "Node"                                   )   return NODE;
    if ( splace == "IntegrationPoint"                       )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "SimplexIntegrationPoint"                )   return ELEMENT_INTEGRATION_POINT;
    if ( splace == "FacetIntegrationPoint"                  )   return FACET_INTEGRATION_POINT;
    if ( splace == "SectorIntegrationPoint"                 )   return SECTOR_INTEGRATION_POINT;
    if ( splace == "FaceIntegrationPoint"                   )   return FACE_INTEGRATION_POINT;
    if ( splace == "FaceFacetIntegrationPoint"              )   return FACE_FACET_INTEGRATION_POINT;
    if ( splace == "FaceSectorIntegrationPoint"             )   return FACE_SECTOR_INTEGRATION_POINT;
    if ( splace == "InterFaceIntegrationPoint"              )   return INTER_FACE_INTEGRATION_POINT;
    if ( splace == "InterFaceFacetIntegrationPoint"         )   return INTER_FACE_FACET_INTEGRATION_POINT;
    if ( splace == "InterFaceSectorIntegrationPoint"        )   return INTER_FACE_SECTOR_INTEGRATION_POINT;
    if ( splace == "Face"                                   )   return FACE;
    if ( splace == "InterFace"                              )   return INTER_FACE;
    if ( splace == "Element"                                )   return ELEMENT;
    if ( splace == "Region"                                 )   return REGION;
    if ( splace == "Boundary"                               )   return BOUNDARY;
    if ( splace == "SplitBoundary"                          )   return SPLIT_BOUNDARY;
    if ( splace == "Model"                                  )   return MODEL;

    std::cout <<"\nparsePlacement(const char*): unable to parse: '";
    std::cout << placement <<"' returning MODEL"<< endl;
    return MODEL;
 }

std::string  parsePlacement( PLACEMENT splace )
 {
    if ( splace == NODE                               ) return string("NODE");
    if ( splace == ELEMENT_INTEGRATION_POINT          ) return string("ELEMENT_INTEGRATION_POINT");
    if ( splace == FACET_INTEGRATION_POINT            ) return string("FACET_INTEGRATION_POINT");
    if ( splace == SECTOR_INTEGRATION_POINT           ) return string("SECTOR_INTEGRATION_POINT");
    if ( splace == FACE_INTEGRATION_POINT             ) return string("FACE_INTEGRATION_POINT");
    if ( splace == FACE_FACET_INTEGRATION_POINT       ) return string("FACE_FACET_INTEGRATION_POINT");
    if ( splace == FACE_SECTOR_INTEGRATION_POINT      ) return string("FACE_SECTOR_INTEGRATION_POINT");
    if ( splace == INTER_FACE_INTEGRATION_POINT       ) return string("INTER_FACE_INTEGRATION_POINT");
    if ( splace == INTER_FACE_FACET_INTEGRATION_POINT ) return string("INTER_FACE_FACET_INTEGRATION_POINT");
    if ( splace == INTER_FACE_SECTOR_INTEGRATION_POINT) return string("INTER_FACE_SECTOR_INTEGRATION_POINT");
    if ( splace == FACE                               ) return string("FACE");
    if ( splace == INTER_FACE                         ) return string("INTER_FACE");
    if ( splace == ELEMENT                            ) return string("ELEMENT");
    if ( splace == REGION                             ) return string("REGION");
    if ( splace == BOUNDARY                           ) return string("BOUNDARY");
    if ( splace == SPLIT_BOUNDARY                     ) return string("SPLIT_BOUNDARY");
    if ( splace == MODEL                              ) return string("MODEL");

    cout <<"\nparsePlacement(PLACEMENT): unable to parse const char*: "<< splace << endl;
    return string("MODEL");
 }
 
 
 
 
bool  faceVariable( PLACEMENT place )
 {
    if ( place == FACE ) return true;
    if ( place == FACE_INTEGRATION_POINT ) return true;
    if ( place == FACE_FACET_INTEGRATION_POINT ) return true;
    if ( place == FACE_SECTOR_INTEGRATION_POINT ) return true;
    return false;
 }
 
 
bool  interFaceVariable( PLACEMENT place )
 {
    if ( place == INTER_FACE ) return true;
    if ( place == INTER_FACE_INTEGRATION_POINT ) return true;
    if ( place == INTER_FACE_FACET_INTEGRATION_POINT ) return true;
    if ( place == INTER_FACE_SECTOR_INTEGRATION_POINT ) return true;
    return false;
 }
 
 

PLACEMENT intToPLACEMENT( int i )
  {
    if ( i == 0  )  return NODE;
    if ( i == 1  )  return ELEMENT_INTEGRATION_POINT;
    if ( i == 2  )  return SECTOR_INTEGRATION_POINT;
    if ( i == 3  )  return FACET_INTEGRATION_POINT;
    if ( i == 5  )  return FACE;
    if ( i == 6  )  return INTER_FACE;
    if ( i == 7  )  return ELEMENT;
    if ( i == 8  )  return REGION;
    if ( i == 9  )  return BOUNDARY;
    if ( i == 10 )  return SPLIT_BOUNDARY;
    if ( i == 11 )  return MODEL;
    if ( i == 12 )  return FACE_INTEGRATION_POINT;
    if ( i == 13 )  return FACE_SECTOR_INTEGRATION_POINT;
    if ( i == 14 )  return FACE_FACET_INTEGRATION_POINT;
    if ( i == 15 )  return INTER_FACE_INTEGRATION_POINT;
    if ( i == 16 )  return INTER_FACE_SECTOR_INTEGRATION_POINT;
    if ( i == 17 )  return INTER_FACE_FACET_INTEGRATION_POINT;

    cout <<"\nintToPlacement(int): unable to parse integer: "<< i << endl;
    return MODEL;
  }


bool isPlacedOnIntegrationPoint( PLACEMENT p )
 {
     if ( p == ELEMENT_INTEGRATION_POINT ) return true;
     if ( p == SECTOR_INTEGRATION_POINT ) return true;
     if ( p == FACET_INTEGRATION_POINT ) return true;

     if ( p == FACE_INTEGRATION_POINT ) return true;
     if ( p == FACE_FACET_INTEGRATION_POINT ) return true;
     if ( p == FACE_SECTOR_INTEGRATION_POINT ) return true;

     if ( p == INTER_FACE_INTEGRATION_POINT ) return true;
     if ( p == INTER_FACE_SECTOR_INTEGRATION_POINT ) return true;
     if ( p == INTER_FACE_FACET_INTEGRATION_POINT ) return true;
   
     return false;
 }



// VARIABLE FLAG

/**
    Options numbered in this sequence, starting with zero are:
    PLAIN,   ANY, INIT_GUESS,  INIT_COND,  FIELD_DATA,  PERIODIC,   DIRICH,   NEUMANN,   ROBIN,  CONSTANT_FLUX
*/
VARIABLE_FLAG intToVARIABLE_FLAG( int i )
  {
    if      ( i == 0 ) return PLAIN;
    else if ( i == 1 ) return ANY;
    else if ( i == 2 ) return INIT_GUESS;
    else if ( i == 3 ) return INIT_COND;
    else if ( i == 4 ) return FIELD_DATA;
    else if ( i == 5 ) return PERIODIC;
    else if ( i == 6 ) return DIRICH;
    else if ( i == 7 ) return NEUMANN;
    else if ( i == 8 ) return ROBIN;
    else if ( i == 9 ) return CONSTANT_FLUX;
    else
    std::cerr <<"\nintToVARIABLE_FLAG(int): unable to parse integer: "<< i << std::endl;
    return ANY;
  }

/**
    Converts the text string into any of the possible csmp::VARIABLE_FLAG
    values.

    @return if the flag cannot be parsed, the value ANY is returned.
*/
VARIABLE_FLAG  parseStatus( const char* status )
 {
    std::string  sstatus(status);
    // to increase performance, the flag values are listed in an order of decreasing likelihood
    if ( sstatus == "ANY" )        return ANY;
    if ( sstatus == "any" )        return ANY;
    if ( sstatus == "PLAIN" )      return PLAIN;
    if ( sstatus == "Plain" )      return PLAIN;
    if ( sstatus == "plain" )      return PLAIN;
    if ( sstatus == "DIRICH" )     return DIRICH;
    if ( sstatus == "DIRICHLET" )  return DIRICH;
    if ( sstatus == "Dirichlet" )  return DIRICH;
    if ( sstatus == "NEUMANN" )    return NEUMANN;
    if ( sstatus == "Neumann" )    return NEUMANN;
    if ( sstatus == "PERIODIC" )   return PERIODIC;
    if ( sstatus == "periodic" )   return PERIODIC;
    if ( sstatus == "ROBIN" )      return ROBIN;
    if ( sstatus == "Robin" )      return ROBIN;
    if ( sstatus == "FIELD_DATA" ) return FIELD_DATA;
    if ( sstatus == "field data" ) return FIELD_DATA;
    if ( sstatus == "INIT_GUESS" )    return INIT_GUESS;
    if ( sstatus == "initial guess" ) return INIT_GUESS;
    if ( sstatus == "INIT_COND" )         return INIT_COND;
    if ( sstatus == "initial condition" ) return INIT_COND;

    std::cerr <<"\nparseStatus(const char*): unable to parse: '"<< status <<"'"<< std::endl;
    return ANY;
 }


/**
    Convert VARIABLE_FLAG from an enumeration to a string so that it can be printed to stdout.
*/
std::string  parseStatus( VARIABLE_FLAG status )
 {
    if ( status == PLAIN ) return std::string("PLAIN");
    if ( status == DIRICH ) return std::string("DIRICH");
    if ( status == NEUMANN ) return std::string("NEUMANN");
    if ( status == INIT_GUESS ) return std::string("INIT_GUESS");
    if ( status == INIT_COND ) return std::string("INIT_COND");
    if ( status == ROBIN ) return std::string("ROBIN");
    if ( status == FIELD_DATA ) return std::string("FIELD_DATA");
    if ( status == PERIODIC ) return std::string("PERIODIC");
    if ( status == ANY ) return std::string("ANY");

    std::cerr <<"\nparseStatus(VARIABLE_FLAG): unable to parse const char*"<< status << std::endl;
    return std::string("ANY");
 }

/**

Interprets the supplied string as an essential (boundary) condition.
Possible values are given in 'CSMP_definitions.h'
and MULTIPLE. The last type can be used to specify Cauchy or
Robbin conditions.

@section arguments Input Arguments

The string that shall be interpreted as boundary condition.
*/
VARIABLE_FLAG  parseCondition( std::string& s )
{
    std::string S(s);
    std::transform(S.begin(),S.end(),S.begin(),::toupper);

    if ( S == "PLAIN"                                                          ) return PLAIN;
    if ( S == "ANY"                                                            ) return ANY;
    if ( S == "INIT_GUESS" || S == "INITIAL GUESS"                             ) return INIT_GUESS;
    if ( S == "INIT_COND"  || S == "INITIAL CONDITION"                         ) return INIT_COND;
    if ( S == "FIELD_DATA" || S == "FIELD DATA"                                ) return FIELD_DATA;
    if ( S == "PERIODIC"                                                       ) return PERIODIC;
    if ( S == "DIRICH"     || S == "DIRICHLET"                                 ) return DIRICH;
    if ( S == "NEUMANN"    || S == "NATURAL"                                   ) return NEUMANN;
    if ( S == "ROBIN"                                                          ) return ROBIN;
    if ( S == "CONST_FLUX" || S == "CONSTANT_FLUX"     || S == "CONSTANT FLUX" ) return CONSTANT_FLUX;

    throw csmp::Exception( FATAL_ERROR,
                           "parseCondition:",
                           S.c_str(),
                           ": boundary condition type specifier was not recognized");
    return PLAIN;
}



/**
    To print enum INTERFACE_SIDE.
*/
string parseSide( INTERFACE_SIDE side )
 {
    if ( side == INSIDE )  return "INSIDE";
    if ( side == OUTSIDE ) return "OUTSIDE";
    if ( side == MIDDLE )  return "MIDDLE";
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
        case ManifoldType::SPLIT_BOUNDARY_TERMINATION: return "SPLIT_BOUNDARY_TERMINATION";
        case ManifoldType::SPLIT_BOUNDARY_END:  return "SPLIT_BOUNDARY_END";
        //default: return "NOT_CLASSIFIED";
      }
    return "SPLIT_BOUNDARY";
 }


/**
       Topologic qualifiers include:
       
      MESH_VERTEX,              - a point within the model volume
      INTERSECTION_POINT  - a point where lines cross or multiple surfaces intersect
      PERIMETER_POINT,       - point at the end of a line inside a 2D model
      EXTERIOR_POINT,         - on an outside surface of the model
      INTERIOR_LINE,            - a line on the interior of the model
      PERIMETER_LINE,        -  a surface edge inside of the model
      EXTERIOR_LINE,          - an edge of the model
      INTERSECTION_LINE, - belonging to multiple surfaces
      INTERIOR_SURFACE,    - a surface withing the model
      PERIMETER_SURFACE,    - a surface forming the hull of an object inside of the model
      EXTERIOR_SURFACE       - a surface delimiting the model
      
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
        //default: return "NOT_CLASSIFIED";
      }
    return "STAND_ALONE";
 }




} // end csmp
