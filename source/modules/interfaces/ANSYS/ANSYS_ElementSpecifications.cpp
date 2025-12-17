#include "ANSYS_ElementSpecifications.h"

using namespace std;

namespace csmp {

/** Returns the ANSYS integer flag corresponding to the ANSYS element name
string.  
*/
int8_t  ANSYS_ElementSpecifications::ANSYS_Type( const string& icem_FEtype )
 {
    // bar elements
    if ( icem_FEtype == "BAR_2" )    return 2;
    if ( icem_FEtype == "BAR_3" )    return 3;
    // trapezoidal elements (quads)
    if ( icem_FEtype == "QUAD_4" )   return 14;
    if ( icem_FEtype == "QUAD_4_X" ) return 15;
    if ( icem_FEtype == "QUAD_8" )   return 16;
    if ( icem_FEtype == "QUAD_8_X" ) return 17;
    if ( icem_FEtype == "QUAD_9" )   return 19;
    // triangles
    if ( icem_FEtype == "TRI_3" )    return 8;
    if ( icem_FEtype == "TRI_3_X" )  return 9;
    if ( icem_FEtype == "TRI_6" )    return 10;
    if ( icem_FEtype == "TRI_6_X" )  return 11;
    // tetrahedral elements
    if ( icem_FEtype == "TETRA_4" )  return 4;
    if ( icem_FEtype == "TETRA_10" ) return 5;
    // hexahedral elements
    if ( icem_FEtype == "HEXA_8" )   return 6;
    if ( icem_FEtype == "HEXA_20" )  return 7;
    if ( icem_FEtype == "HEXA_27" )  return 20;
    // pyramids
    if ( icem_FEtype == "PYRA_5" )   return 18;
    if ( icem_FEtype == "PYRA_13" )  return 24;
    if ( icem_FEtype == "PYRA_14" )  return 22;
    // pentahedra
    if ( icem_FEtype == "PENTA_6" )  return 12;
    if ( icem_FEtype == "PENTA_15" ) return 13;
    if ( icem_FEtype == "PENTA_18" ) return 21;
    if ( icem_FEtype == "POLYGON" )  return 23;

    cout <<"\nANSYS_ElementSpecifications::ANSYS_Type: ";
    cout <<" Element type appears different from ANSYS types ?";
    cout <<"\nType: "<< icem_FEtype << endl;
         
    return 0;
    
 } // end ANSYS_Type




    
    



/** Returns the number of nodes of the ANSYS element type as identified by
the ANSYS element number.  

ANSYS defines the following element identification strings and flags:  

@code
element type     identifier
#define BAR_2		2	#define TETRA_4		4
#define BAR_3		3	#define TETRA_10	5

#define HEXA_8		6	#define QUAD_4		14	#define QUAD_8_X	17
#define HEXA_20		7	#define QUAD_4_X	15	#define QUAD_9		19
#define HEXA_27		20	#define QUAD_8		16

#define TRI_3		8
#define TRI_3_X		9
#define TRI_6		10
#define TRI_6_X		11

#define PENTA_6		12	#define PYRA_5		18
#define PENTA_15	13	#define PYRA_14		22
#define PENTA_18	21	#define PYRA_13		24

#define POLYGON		23
@endcode
*/
size_t ANSYS_ElementSpecifications::NodesPerElementOfType( int8_t etype )
 {
    if ( etype ==  2 ) return  2U; // bar element
    if ( etype ==  3 ) return  3U; // 
    if ( etype ==  4 ) return  4U; // tetrahedron
    if ( etype ==  5 ) return 10U; // 
    if ( etype ==  6 ) return  8U; // hexahedron
    if ( etype ==  7 ) return 20U; // 
    if ( etype == 20 ) return 27U; // 
    if ( etype == 14 ) return  4U; // quadrilateral
    if ( etype == 15 ) return  5U; // 
    if ( etype == 16 ) return  8U; // 
    if ( etype == 17 ) return  9U; // 
    if ( etype == 19 ) return  9U; // 
    if ( etype ==  8 ) return  3U; // linear triangle
    if ( etype ==  9 ) return  4U; // 
    if ( etype == 10 ) return  6U; // 
    if ( etype == 11 ) return  7U; // 
    if ( etype == 12 ) return  6U; // pentahedron = prism element
    if ( etype == 13 ) return 15U; // 
    if ( etype == 21 ) return 18U; // 
    if ( etype == 18 ) return  5U; // pyramid
    if ( etype == 22 ) return 14U; // 
    if ( etype == 24 ) return 13U; // 
    if ( etype == 23 ) return  0U; // polygon
    
    return 0U; // unknown element type
 
 } // end  
    



/** Returns the number of neighbors of the input ANSYS finite element type.
*/
size_t ANSYS_ElementSpecifications::NeighborsPerElementOfType( int8_t etype )
 {
    if ( etype ==  2 ) return  2U; // bar element
    if ( etype ==  3 ) return  2U; // 
    if ( etype ==  4 ) return  4U; // tetrahedron
    if ( etype ==  5 ) return  4U; // 
    if ( etype ==  6 ) return  6U; // hexahedron
    if ( etype ==  7 ) return  6U; // 
    if ( etype == 20 ) return  6U; // 
    if ( etype == 14 ) return  4U; // quadrilateral
    if ( etype == 15 ) return  4U; // 
    if ( etype == 16 ) return  4U; // 
    if ( etype == 17 ) return  4U; // 
    if ( etype == 19 ) return  4U; // 
    if ( etype ==  8 ) return  3U; // linear triangle
    if ( etype ==  9 ) return  3U; // 
    if ( etype == 10 ) return  3U; // 
    if ( etype == 11 ) return  3U; // 
    if ( etype == 12 ) return  5U; // pentahedron = prism element
    if ( etype == 13 ) return  5U; // 
    if ( etype == 21 ) return  5U; // 
    if ( etype == 18 ) return  5U; // pyramid
    if ( etype == 22 ) return  5U; // 
    if ( etype == 24 ) return  5U; // 
    if ( etype == 23 ) return 10U; // polygon
   
    cout <<"\nANSYS_ElementSpecifications::NeighborsPerElementOfType: ";
    cout <<"unable to parse element type, returning 0"<< endl; 
    return 0U; // unknown element type
 
 } // end  
    



size_t ANSYS_ElementSpecifications::NodesPerFaceForElementOfType( int8_t etype,
                                                                    size_t face )
 {
    // bar element
    if ( etype ==  2 or etype == 3 ) { 
         assert( face < 2U );
         return 1U;
      }
    // tetrahedron
    if ( etype ==  3  or etype == 4 or etype == 5 ) { 
         assert( face < 4U );
         return 3U;
      }
    // hexahedron
    if ( etype ==  6 or etype == 7 or etype == 20 ) {
         assert( face < 6U );
         return 4U;
      }
    // quadrilateral
    if ( etype == 14 or etype == 15 or etype == 16 or etype == 17 or etype == 19 ) {
         assert( face < 4U );
         return 2U;
      }
    // linear triangle
    if ( etype == 8 or etype == 9 or etype == 10 or etype == 11 ) {
         assert( face < 4U );
         return 2U;
      }    
    // pentahedron = prism element
    if ( etype == 12 or etype == 13 or etype == 21 ) {
         assert( face < 5U );
         if ( face<=2U ) return 4U;
         return 3U;
      }

    // polygon
    if ( etype == 23 ) return ULONG_MAX; // polygon
    
    cout <<"\nANSYS_ElementSpecifications::NodesPerFaceForElementOfType: ";
    cout <<"unable to parse element type, returning ULONG_MAX"<< endl; 
    return ULONG_MAX; // unknown number of faces
    
 } // end NodesPerFaceForElementOfType
 


 

  
size_t ANSYS_ElementSpecifications::FaceNodeForElementOfType( int8_t etype,
                                                              size_t face, 
                                                              size_t face_node )
 {
    // bar element
    if ( etype ==  2 or etype == 3 ) { 
         assert( face < 2U );
         assert( face_node == 0U );
         // face 1 is adjacent to node 0 and 1 to node 1
         if ( face==0U ) return 0U;
         if ( face==1U ) return 1U;
      }
    // tetrahedron
    if ( etype ==  3  or etype == 4 or etype == 5 ) { 
         assert( face < 4U );
         assert( face_node < 3U );
         if ( face==0U ) {
              if ( face_node==0U ) return 1U;
              if ( face_node==1U ) return 2U;
              if ( face_node==2U ) return 3U;
           }
         if ( face==1U ) {
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 2U;
              if ( face_node==2U ) return 3U;
           }
         if ( face==2U ) {
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 1U;
              if ( face_node==2U ) return 3U;
           }
         if ( face==3U ) { 
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 1U;
              if ( face_node==2U ) return 2U;
           }
      }
    // hexahedron
    if ( etype ==  6 or etype == 7 or etype == 20 ) {
         assert( face < 6U );
         assert( face_node < 4U );
         if ( face==0U ) {
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 3U;
              if ( face_node==2U ) return 2U;
              if ( face_node==3U ) return 1U;
           }
         if ( face==1U ) {
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 1U;
              if ( face_node==2U ) return 5U;
              if ( face_node==3U ) return 4U;
           }
         if ( face==2U ) {
              if ( face_node==0U ) return 1U;
              if ( face_node==1U ) return 2U;
              if ( face_node==2U ) return 6U;
              if ( face_node==3U ) return 5U;
           }
         if ( face==3U ) { 
              if ( face_node==0U ) return 2U;
              if ( face_node==1U ) return 3U;
              if ( face_node==2U ) return 7U;
              if ( face_node==3U ) return 6U;
           }
         if ( face==4U ) {
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 4U;
              if ( face_node==2U ) return 7U;
              if ( face_node==3U ) return 3U;
           }
         if ( face==5U ) { 
              if ( face_node==0U ) return 4U;
              if ( face_node==1U ) return 5U;
              if ( face_node==2U ) return 6U;
              if ( face_node==3U ) return 7U;
           }
      }
    
    // quadrilateral
    if ( etype == 14 or etype == 15 or etype == 16 or etype == 17 or etype == 19 ) {
         assert( face < 4U );
         assert( face_node < 2U );
         if ( face==0U ) {
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 1U;
           }
         if ( face==1U ) {
              if ( face_node==0U ) return 1U;
              if ( face_node==1U ) return 2U;
           }
         if ( face==2U ) {
              if ( face_node==0U ) return 2U;
              if ( face_node==1U ) return 3U;
           }
         if ( face==3U ) { 
              if ( face_node==0U ) return 3U;
              if ( face_node==1U ) return 0U;
           }
      }
      
    // linear triangle
    if ( etype == 8 or etype == 9 or etype == 10 or etype == 11 ) {
         assert( face < 3U );
         assert( face_node < 2U );
         if ( face==0U ) {
              if ( face_node==0U ) return 1U;
              if ( face_node==1U ) return 2U;
           }
         if ( face==1U ) {
              if ( face_node==0U ) return 2U;
              if ( face_node==1U ) return 0U;
           }
         if ( face==2U ) {
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 1U;
           }
      }    
    // pentahedron = prism element
    if ( etype == 12 or etype == 13 or etype == 21 ) {
         assert( face < 5U );
         if ( face==0U ) { 
              assert( face_node < 3U );
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 2U;
              if ( face_node==2U ) return 1U;
           }
         if ( face==1U ) {
              assert( face_node < 4U );
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 1U;
              if ( face_node==2U ) return 4U;
              if ( face_node==3U ) return 3U;
           }
         if ( face==2U ) {
              assert( face_node < 4U );
              if ( face_node==0U ) return 1U;
              if ( face_node==1U ) return 2U;
              if ( face_node==2U ) return 5U;
              if ( face_node==3U ) return 4U;
           }
         if ( face==3U ) {
              assert( face_node < 4U );
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 3U;
              if ( face_node==2U ) return 5U;
              if ( face_node==3U ) return 2U;
           }
         if ( face==4U ) {
              assert( face_node < 3U );
              if ( face_node==0U ) return 3U;
              if ( face_node==1U ) return 4U;
              if ( face_node==2U ) return 5U;
           }
       }

    // pyramid
    if ( etype == 18 or etype == 22 or etype == 24 ) {
         assert( face < 5U );
         if ( face==0U ) { 
              assert( face_node < 3U );
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 1U;
              if ( face_node==2U ) return 4U;
           }
         if ( face==1U ) {
              assert( face_node < 4U );
              if ( face_node==0U ) return 1U;
              if ( face_node==1U ) return 2U;
              if ( face_node==2U ) return 4U;
           }
         if ( face==2U ) {
              assert( face_node < 4U );
              if ( face_node==0U ) return 2U;
              if ( face_node==1U ) return 3U;
              if ( face_node==2U ) return 4U;
           }
         if ( face==3U ) {
              assert( face_node < 4U );
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 4U;
              if ( face_node==2U ) return 3U;
           }
         if ( face==4U ) {
              assert( face_node < 3U );
              if ( face_node==0U ) return 0U;
              if ( face_node==1U ) return 3U;
              if ( face_node==2U ) return 2U;
              if ( face_node==3U ) return 1U;
           }
      }

    return ULONG_MAX;

 } // end FaceNodeForElementOfType
  
 


size_t ANSYS_ElementSpecifications::FacesPerElementOfType( int8_t ANSYS_finite_element_type )
 {
    // since these are equivalent numbers
    return NeighborsPerElementOfType( ANSYS_finite_element_type );
    
 } // end FacesPerElementOfType




/** Returns true if the ANSYS input element type is a line. Possible ANSYS
types are:  
 
@code
element type     identifier
#define BAR_2		2
#define BAR_3		3
@endcode*/
bool  ANSYS_ElementSpecifications::LineElement( int8_t etype )
 {
    if ( etype == 2 || etype == 3 ) return true;
    
    return false;
 }
 
 
 
/** Returns true if the ANSYS input element type is a surface element.
Possible surface element types are:  
 
@code
element type     identifier
#define QUAD_4		14	#define QUAD_8_X	17
#define QUAD_4_X	15	#define QUAD_9		19
#define QUAD_8		16

#define TRI_3		  8
#define TRI_3_X		9
#define TRI_6		 10
#define TRI_6_X	 11
@endcode
*/
bool  ANSYS_ElementSpecifications::SurfaceElement( int8_t etype )
 {
    if ( etype ==  8 || etype ==  9 || etype == 10 || 
         etype == 11 || etype == 14 || etype == 15 || 
         etype == 16 || etype == 17 || etype == 19 ) return true;
    
    return false;
 }
 
 
 
/** Returns true if the ANSYS input element type is a volumetric element.
Possible ANSYS volumetric element types are:  

@code
element type     identifier
#define TETRA_4		4	#define TETRA_10	5

#define HEXA_8		6	
#define HEXA_20		7	
#define HEXA_27		20	

#define PENTA_6		12	#define PYRA_5		18
#define PENTA_15	13	#define PYRA_14		22
#define PENTA_18	21	#define PYRA_13		24
@endcode */
bool  ANSYS_ElementSpecifications::VolumeElement( int8_t etype )
 {
    if ( etype ==  4 || etype ==  5 || etype ==  6 || 
         etype ==  7 || etype == 12 || etype == 13 || 
         etype == 18 || etype == 20 || etype == 21 ||
         etype == 22 || etype == 24 ) return true;
    
    return false;
 }
 

 
/**
 
The element type used defines the spatial dimension of the model. It is
indicative only of the minimum dimension since lines and surfaces can 
also exist in 3D space, which the name makes explicit.   

@section arguments Input Arguments 

Method expects a string argument that contains the ANSYS element type
to deduce the spatial dimension from.  

@return size_t The minimum spatial dimension is returned as an unsigned integer ranging
from 1(=1D) to 3(=3D).  

*/
size_t  ANSYS_ElementSpecifications::MinimumSpatialDimension( const string& etype )
 {
    if ( LineElement(etype) )    return 1U;
    if ( SurfaceElement(etype) ) return 2U;
    
    return 3U;
 }
 
 
 
 
/** Returns true if the ANSYS element name corresponds to a line element.
Possible line element types are:  

@code
element type     identifier
BAR_2		2
BAR_3		3
@endcode 

@return bool whether the element type is a line element.

*/
bool  ANSYS_ElementSpecifications::LineElement( const string& etype )
 {
    if ( etype == "BAR_2" || etype == "BAR_3" ) return true;
    
    return false;
 }
 
 
 
/**
 
Returns true is the ANSYS input element name is a surface element. Possible
names of ANSYS surface elements are:  

@code
element type     identifier
TRI_3		8
TRI_3_X		9
TRI_6		10
TRI_6_X		11
QUAD_4		14					
QUAD_4_X	15					
QUAD_8		16
QUAD_8_X	17
QUAD_9		19
@endcode */
bool  ANSYS_ElementSpecifications::SurfaceElement( const string& etype )
 {
    if ( etype == "TRI_3"   || etype == "TRI_3_X"   || etype == "TRI_6"    || 
         etype == "TRI_6_X" || etype == "QUAD_4"    || etype == "QUAD_4_X" || 
         etype == "QUAD_8"  || etype == "QUAD_8_X"  || etype == "QUAD_9" ) 
      return true;
    
    return false;
 }



 
/**
 
Return true if the ANSYS input element name corresponds to a volumetric 
element. Possible ANSYS volumetric element types are: 
 
ANSYS defines the following element identification strings and flags

@code
element type     identifier
#define TETRA_4		4	#define TETRA_10	5

#define HEXA_8		6	
#define HEXA_20		7	
#define HEXA_27		20	

#define PENTA_6		12	#define PYRA_5		18
#define PENTA_15	13	#define PYRA_14		22
#define PENTA_18	21	#define PYRA_13		24

@endcode */
bool  ANSYS_ElementSpecifications::VolumeElement( const string& etype )
 {
    if ( etype == "TETRA_4"  || etype == "TETRA_10" || etype == "HEXA_8"  || 
         etype == "HEXA_20"  || etype == "HEXA_27"  || etype == "PENTA_6" || 
         etype == "PENTA_15" || etype == "PENTA_18" || etype == "PYRA_5"  || 
         etype == "PYRA_14"  || etype == "PYRA_13" ) return true;
    
    return false;
 }    
    



/**
 
Method deduces order of element interpolation functions from the 
input finite-element type. 

@section arguments Input Arguments 

The ANSYS type name of the finite element that shall be queried.  

@return size_t giving the order of interpolation (1=linear, 2=quadratic etc.)

The interpolation order is returned as an unsigned integer value:
1 = linear, 2 = quadratic, 3 = cubic.  
*/
size_t  ANSYS_ElementSpecifications::InterpolationOrder( const std::string& etype )
 {
     if ( LinearElement( ANSYS_Type(etype) ) )    return 1U;
     if ( QuadraticElement( ANSYS_Type(etype) ) ) return 2U;
     if ( CubicElement( ANSYS_Type(etype) ) )     return 3U;

     return 0U; // order not identified, probably constant as in a polygon
 }


size_t  ANSYS_ElementSpecifications::InterpolationOrder( int8_t etype )
  {
    if ( LinearElement( etype ) )    return 1U;
    if ( QuadraticElement( etype ) ) return 2U;
    if ( CubicElement( etype ) )     return 3U;

    return 0U; // order not identified, probably constant as in a polygon
  }


/**
 
Returns whether the ANSYS element type specifier denotes a finite
element type that uses linear basis functions. ANSYS element types with
such characteristics are:  

@code
linear
==============
BAR_2		 2
TETRA_4		 4
HEXA_8		 6
TRI_3		 8
TRI_3_X		 9
PENTA_6		12
QUAD_4		14
QUAD_4_X	15
PYRA_5		18
@endcode */
bool  ANSYS_ElementSpecifications::LinearElement( int8_t etype )
 {
    if ( etype ==  2 || etype ==  4 || 
         etype ==  6 || etype ==  8 || etype ==  9 || 
         etype == 12 || etype == 14 || etype == 15 ||
         etype == 18 ) return true;
 
    return false;
 }
 
 
 
 
/**
 
Returns whether the ANSYS element type specifier denotes a finite
element type that uses quadratic basis functions. ANSYS element types with
such characteristics are:  

@code
quadratic
==============
TETRA_10	 5
BAR_3 3
HEXA_20		 7
TRI_6		10
TRI_6_X		11
PENTA_15	13
QUAD_8		16
QUAD_8_X	17
QUAD_9		19
PENTA_18	21
PYRA_14		22
PYRA_13		24
@endcode*/
bool  ANSYS_ElementSpecifications::QuadraticElement( int8_t etype )
 {
    if ( etype ==  5 || etype ==  7 || etype == 10 || 
         etype == 11 || etype == 13 || etype == 16 || 
         etype == 17 || etype == 19 || etype == 21 ||
         etype == 22 || etype == 24 || etype == 3 ) return true;

    return false;
 }
 

 
/**
 
Returns whether the ANSYS element type specifier denotes a finite
element type that uses cubic basis functions. ANSYS element types with
such characteristics are:  

@code
cubic
==============
HEXA_32		32
@endcode*/
bool  ANSYS_ElementSpecifications::CubicElement( int8_t etype )
 {
    if ( etype == 32 ) return true;
    return false;
 }



/**
 
Returns a list of line type elements that includes both the ANSYS types
and the corresponding CSP types (see FiniteElement.h). The output
strings are:  

@code
ANSYS element type
BAR_2  BAR_3  BAR_4
CSMP types
LINEAR_BAR  
ISOPARAMETRIC_LINEAR_BAR
ISOPARAMETRIC_QUADRATIC_BAR
ISOPARAMETRIC_CUBIC_BAR
@endcode */
void  ANSYS_ElementSpecifications::LineElements( list<string>& line_elements )
 {
    if ( !line_elements.empty() ) 
      line_elements.erase( line_elements.begin(), line_elements.end() );
      
    line_elements.push_back("BAR_2");
    line_elements.push_back("BAR_3");
    line_elements.push_back("BAR_4");
    // CSMP_TYPES
    line_elements.push_back("LINEAR_BAR");
    line_elements.push_back("ISOPARAMETRIC_LINEAR_BAR");
    line_elements.push_back("ISOPARAMETRIC_QUADRATIC_BAR");
    line_elements.push_back("ISOPARAMETRIC_CUBIC_BAR");
 }
 
 
 
 
/**

Returns a list of the ANSYS surface element type names and the
corresponding CSP types. These are: 
 
@code
ANSYS element type     identifier
TRI_3		8
TRI_3_X		9
TRI_6		10
TRI_6_X		11
QUAD_4		14					
QUAD_4_X	15					
QUAD_8		16
QUAD_8_X	17
QUAD_9		19
CSP types
LINEAR_TRIANGLE
LINEAR_TRIANGLE3D
ISOPARAMETRIC_LINEAR_TRIANGLE
BARYCENTRIC_LINEAR_TRIANGLE
QUADRATIC_TRIANGLE
BARYCENTRIC_QUADRATIC_TRIANGLE
ISOPARAMETRIC_QUADRATIC_TRIANGLE
ISOPARAMETRIC_CUBIC_TRIANGLE
ISOPARAMETRIC_LINEAR_QUADRILATERAL
ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL
ISOPARAMETRIC_QUADRATIC_QUADRILATERAL
ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL
ISOPARAMETRIC_CUBIC_QUADRILATERAL
@endcode*/
void  ANSYS_ElementSpecifications::SurfaceElements( list<string>& surf_elements )
 {
    if ( !surf_elements.empty() ) 
      surf_elements.erase( surf_elements.begin(), surf_elements.end() );
      
    surf_elements.push_back("TRI_3");
    surf_elements.push_back("TRI_3_X");
    surf_elements.push_back("TRI_6");
    surf_elements.push_back("TRI_6_X");
    surf_elements.push_back("QUAD_4");
    surf_elements.push_back("QUAD_4_X");
    surf_elements.push_back("QUAD_8");
    surf_elements.push_back("QUAD_8_X");
    surf_elements.push_back("QUAD_9");
    // CSMP Types
	surf_elements.push_back("LINEAR_RECTANGLE");
    surf_elements.push_back("LINEAR_TRIANGLE");
    surf_elements.push_back("LINEAR_TRIANGLE3D");
    surf_elements.push_back("ISOPARAMETRIC_LINEAR_TRIANGLE");
    surf_elements.push_back("BARYCENTRIC_LINEAR_TRIANGLE");
    surf_elements.push_back("QUADRATIC_TRIANGLE");
    surf_elements.push_back("BARYCENTRIC_QUADRATIC_TRIANGLE");
    surf_elements.push_back("ISOPARAMETRIC_QUADRATIC_TRIANGLE");
    surf_elements.push_back("ISOPARAMETRIC_CUBIC_TRIANGLE");
    surf_elements.push_back("ISOPARAMETRIC_LINEAR_QUADRILATERAL");
    surf_elements.push_back("ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL");
    surf_elements.push_back("ISOPARAMETRIC_QUADRATIC_QUADRILATERAL");
    surf_elements.push_back("ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL");
    surf_elements.push_back("ISOPARAMETRIC_CUBIC_QUADRILATERAL");
 }


 

/**
 
Returns a list of the ANSYS element type names that denote volumetric 
elements and the corresponding CSP types. These are:  

@code
ANSYS element type     identifier
#define TETRA_4		4	#define TETRA_10	5
#define HEXA_8		6	
#define HEXA_20		7	
#define HEXA_27		20	
#define PENTA_6		12	#define PYRA_5		18
#define PENTA_15	13	#define PYRA_14		22
#define PENTA_18	21	#define PYRA_13		24
CSMP element type
LINEAR_TETRAHEDRON
ISOPARAMETRIC_LINEAR_TETRAHEDRON
ISOPARAMETRIC_LINEAR_HEXAHEDRON
ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20
ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27
ISOPARAMETRIC_LINEAR_PYRAMID
ISOPARAMETRIC_QUADRATIC_PYRAMID13
ISOPARAMETRIC_QUADRATIC_PYRAMID14
ISOPARAMETRIC_LINEAR_PRISM
ISOPARAMETRIC_LINEAR_PRISM
ISOPARAMETRIC_QUADRATIC_PRISM15
ISOPARAMETRIC_QUADRATIC_PRISM18
@endcode*/
void  ANSYS_ElementSpecifications::VolumeElements( list<string>& vol_elements )
 {
    if ( !vol_elements.empty() ) 
      vol_elements.erase( vol_elements.begin(), vol_elements.end() );
    
    // ANSYS types
    vol_elements.push_back("TETRA_4");
    vol_elements.push_back("TETRA_10");
    vol_elements.push_back("HEXA_8");
    vol_elements.push_back("HEXA_20");
    vol_elements.push_back("HEXA_27");
    vol_elements.push_back("PENTA_6");
    vol_elements.push_back("PENTA_15");
    vol_elements.push_back("PENTA_18");
    vol_elements.push_back("PYRA_5");
    vol_elements.push_back("PYRA_14");
    vol_elements.push_back("PYRA_13");
    // CSMP_TYPES
    vol_elements.push_back("LINEAR_TETRAHEDRON");
    vol_elements.push_back("ISOPARAMETRIC_LINEAR_TETRAHEDRON");
    vol_elements.push_back("ISOPARAMETRIC_LINEAR_HEXAHEDRON");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27");
    vol_elements.push_back("ISOPARAMETRIC_LINEAR_PYRAMID");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_PYRAMID13");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_PYRAMID14");
    vol_elements.push_back("ISOPARAMETRIC_LINEAR_PRISM");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_PRISM15");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_PRISM18");
 }




/**
 
Returns the CSP finite element type name that corresponds to the ANSYS type
as identified by an integer identifier. The user needs to specify 
whether an isoparametric or globally interpolated element shall be 
returned because ANSYS elements can be converted into both.  

Note however that ANSYS fits quadratic or cubic elements to the bounding
curves. Thus, global interpolation will give incorrect results as
it works only for straight-sided elements.  
*/
std::string  ANSYS_ElementSpecifications::CSMP_TypeNameFrom_ANSYS_Type( int8_t etype, bool isoparametric, uint32_t dim )
{
    return parseFiniteElementType( CSMP_TypeFrom_ANSYS_Type( etype, isoparametric, dim ) );
}

std::string  ANSYS_ElementSpecifications::CSMP_TypeNameFrom_ANSYS_TypeName( const std::string& etype, bool isoparametric, uint32_t dim )
{
    return parseFiniteElementType( CSMP_TypeFrom_ANSYS_TypeName( etype, isoparametric, dim ) );
}

std::string  parse_ANSYS_BoundaryFlag( int8_t i )
 {
    if ( i == IRREGULAR_OUTSIDE ) return string("IRREGULAR");
    if ( i == LEFT_OUTSIDE )    return string("LEFT");
    if ( i == RIGHT_OUTSIDE )   return string("RIGHT");
    if ( i == BOTTOM_OUTSIDE )  return string("BOTTOM");
    if ( i == TOP_OUTSIDE )     return string("TOP");
    if ( i == FRONT_OUTSIDE )   return string("FRONT");
    if ( i == BACK_OUTSIDE )    return string("BACK");
    if ( i == CNR_MIN )         return string("CNR1");
    if ( i == CNR_X )           return string("CNR2");
    if ( i == CNR_XY )          return string("CNR3");
    if ( i == CNR_Y )           return string("CNR4");
    if ( i == CNR_Z )           return string("CNR5");
    if ( i == CNR_XZ )          return string("CNR6");
    if ( i == CNR_MAX )         return string("CNR7");
    if ( i == CNR_YZ )          return string("CNR8");
    if ( i == BACK_BOTTOM )     return string("EDGE1");
    if ( i == BACK_RIGHT )      return string("EDGE2");
    if ( i == BACK_TOP )        return string("EDGE3");
    if ( i == BACK_LEFT )       return string("EDGE4");
    if ( i == BOTTOM_LEFT )     return string("EDGE5");
    if ( i == BOTTOM_RIGHT )    return string("EDGE6");
    if ( i == TOP_RIGHT )       return string("EDGE7");
    if ( i == TOP_LEFT )        return string("EDGE8");
    if ( i == FRONT_BOTTOM )    return string("EDGE9");
    if ( i == FRONT_RIGHT )     return string("EDGE10");
    if ( i == FRONT_TOP )       return string("EDGE11");
    if ( i == FRONT_LEFT )      return string("EDGE12");
    if ( i == REGION_BOUNDARY ) return string("INTERNAL");

    cout <<"\nparse_ANSYS_BoundaryFlag(int): unable to parse CSMP_BOUNDARY_FLAG: "<< i << endl;
    return string("parse_ANSYS_BoundaryFlag(CSMP_BOUNDARY_FLAG): not resolved.");
 }



/**
 
Returns the CSMP finite element type enum identifier that corresponds to the
ANSYS type as identified by an integer identifier. The user needs to specify 
whether an isoparametric or globally interpolated element shall be 
returned because ANSYS elements can be converted into both.  

Note however that ANSYS fits quadratic or cubic elements to the bounding
curves. Thus, global interpolation will give incorrect results as
it works only for straight-sided elements.  
 */
CSMP_FEM_TYPE  ANSYS_ElementSpecifications::CSMP_TypeFrom_ANSYS_Type( int8_t etype, bool isoparametric, uint32_t dim )
 {
    if ( isoparametric ) {
	    // bar elements
	    if ( etype ==  2 ) return ISOPARAMETRIC_LINEAR_BAR;                           // BAR_2;
	    if ( etype ==  3 ) return ISOPARAMETRIC_QUADRATIC_BAR;                        // BAR_3;
	    // trapezoidal elements (quads)
	    if ( etype == 14 ) return ISOPARAMETRIC_LINEAR_QUADRILATERAL;                 // QUAD_4;
	    if ( etype == 15 ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL;     // QUAD_4_X;
	    if ( etype == 16 ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL;              // QUAD_8;
	    if ( etype == 17 ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL;  // QUAD_8_X;
	    if ( etype == 19 ) return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9;             // QUAD_9;
	    // triangles
	    if ( etype ==  8 ) return ISOPARAMETRIC_LINEAR_TRIANGLE;                      // TRI_3;
	    if ( etype ==  9 ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE;        // TRI_3_X;
	    if ( etype == 10 ) return ISOPARAMETRIC_QUADRATIC_TRIANGLE;                   // TRI_6;
	    if ( etype == 11 ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE;       // TRI_6_X;
	    // tetrahedral elements
	    if ( etype ==  4 ) return ISOPARAMETRIC_LINEAR_TETRAHEDRON;                   // TETRA_4;
	    if ( etype ==  5 ) return ISOPARAMETRIC_QUADRATIC_TETRAHEDRON;                // TETRA_10;
	    // hexahedral elements
	    if ( etype ==  6 ) return ISOPARAMETRIC_LINEAR_HEXAHEDRON;       // HEXA_8;
	    if ( etype ==  7 ) return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20;  // HEXA_20;
	    if ( etype == 20 ) return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27;  // HEXA_27;
	    // pyramids
	    if ( etype == 18 ) return ISOPARAMETRIC_LINEAR_PYRAMID;          // PYRA_5;
	    if ( etype == 24 ) return ISOPARAMETRIC_QUADRATIC_PYRAMID13;     // PYRA_13;
	    if ( etype == 22 ) return ISOPARAMETRIC_QUADRATIC_PYRAMID14;     // PYRA_14;
	    // pentahedra
	    if ( etype == 12 ) return ISOPARAMETRIC_LINEAR_PRISM;            // PENTA_6;
	    if ( etype == 13 ) return ISOPARAMETRIC_QUADRATIC_PRISM15;       // PENTA_15;
	    if ( etype == 21 ) return ISOPARAMETRIC_QUADRATIC_PRISM18;       // PENTA_18;
      }
      
    // straightsided elements using the global coordinate system and exact integration   
    // bar elements
    if ( etype ==  2 ) return LINEAR_BAR;                                         // BAR_2;
    if ( etype ==  3 ) return QUADRATIC_BAR;                        // BAR_3;
    // trapezoidal elements (quads)
    if ( etype == 14 ) return LINEAR_RECTANGLE;                 // QUAD_4;
//    if ( etype == 15 ) return BARYCENTRIC_LINEAR_QUADRILATERAL;     // QUAD_4_X;
//    if ( etype == 16 ) return QUADRATIC_QUADRILATERAL;              // QUAD_8;
//    if ( etype == 17 ) return BARYCENTRIC_QUADRATIC_QUADRILATERAL;  // QUAD_8_X;
//    if ( etype == 19 ) return QUADRATIC_QUADRILATERAL9;             // QUAD_9;
    // triangles
    if ( etype ==  8 and dim == 2U ) return LINEAR_TRIANGLE;          // TRI_3;
    if ( etype ==  8 and dim == 3U ) return LINEAR_TRIANGLE3D;        // TRI_3;
    if ( etype ==  9 ) return BARYCENTRIC_LINEAR_TRIANGLE;        // TRI_3_X;
    if ( etype == 10 ) return QUADRATIC_TRIANGLE;                   // TRI_6;
    if ( etype == 11 ) return BARYCENTRIC_QUADRATIC_TRIANGLE;       // TRI_6_X;
    // tetrahedral elements
    if ( etype ==  4 ) return LINEAR_TETRAHEDRON;                   // TETRA_4;
    if ( etype ==  5 ) return QUADRATIC_TETRAHEDRON;                // TETRA_10;
    // hexahedral elements
    if ( etype ==  6 ) return LINEAR_CUBOID;       // HEXA_8;
//    if ( etype ==  7 ) return QUADRATIC_HEXAHEDRON20;  // HEXA_20;
//    if ( etype == 20 ) return QUADRATIC_HEXAHEDRON27;  // HEXA_27;
    // pyramids
//    if ( etype == 18 ) return LINEAR_PYRAMID;          // PYRA_5;
//    if ( etype == 24 ) return QUADRATIC_PYRAMID13;     // PYRA_13;
//    if ( etype == 22 ) return QUADRATIC_PYRAMID14;     // PYRA_14;
    // pentahedra
//    if ( etype == 12 ) return LINEAR_PRISM;            // PENTA_6;
//    if ( etype == 13 ) return QUADRATIC_PRISM15;       // PENTA_15;
//    if ( etype == 21 ) return QUADRATIC_PRISM18;       // PENTA_18;

    cout <<"\nANSYS_ElementSpecifications::CSP_TypeFromANSYSType";
    if ( isoparametric ) cout <<" (isoparametric): ";
    else                 cout <<": ";
    cout <<" Element type specifier cannot be parsed ";
    cout <<"\nType specifier: "<< etype << std::endl;
    return UNKNOWN;
    
 } // end CSMP_TypeFrom_ANSYS_Type (int, bool)




/**
 
Returns the CSMP finite element type enum identifier that corresponds to the
ANSYS type as identified by a character string. The user needs to specify 
whether an isoparametric or globally interpolated element shall be 
returned because ANSYS elements can be converted into both.  

Note however that ANSYS fits quadratic or cubic elements to the bounding
curves. Thus, global interpolation will give incorrect results as
it works only for straight-sided elements.  
 */
CSMP_FEM_TYPE  ANSYS_ElementSpecifications::CSMP_TypeFrom_ANSYS_TypeName( const std::string& ANSYS_element_type, bool isoparametric, uint32_t dim )
 {
    string etype(ANSYS_element_type);
 
    if ( isoparametric ) {
        // bar elements
        if ( etype == "BAR_2" )    return ISOPARAMETRIC_LINEAR_BAR;
        if ( etype == "BAR_3" )    return ISOPARAMETRIC_QUADRATIC_BAR;
        // trapezoidal elements (quads)
        if ( etype == "QUAD_4" )   return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
        if ( etype == "QUAD_4_X" ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL;
        if ( etype == "QUAD_8" )   return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL;
        if ( etype == "QUAD_8_X" ) return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL;
        if ( etype == "QUAD_9" )   return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9;
        // triangles
        if ( etype == "TRI_3" )    return ISOPARAMETRIC_LINEAR_TRIANGLE;
        if ( etype == "TRI_3_X" )  return ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE;
        if ( etype == "TRI_6" )    return ISOPARAMETRIC_QUADRATIC_TRIANGLE;
        if ( etype == "TRI_6_X" )  return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE;
        // tetrahedral elements
        if ( etype == "TETRA_4" )  return ISOPARAMETRIC_LINEAR_TETRAHEDRON;
        if ( etype == "TETRA_10" ) return ISOPARAMETRIC_QUADRATIC_TETRAHEDRON;
        // hexahedral elements
        if ( etype == "HEXA_8" )   return ISOPARAMETRIC_LINEAR_HEXAHEDRON;
        if ( etype == "HEXA_20" )  return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20;
        if ( etype == "HEXA_27" )  return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27;
        // pyramids
        if ( etype == "PYRA_5" )   return ISOPARAMETRIC_LINEAR_PYRAMID;
        if ( etype == "PYRA_13" )  return ISOPARAMETRIC_QUADRATIC_PYRAMID13;
        if ( etype == "PYRA_14" )  return ISOPARAMETRIC_QUADRATIC_PYRAMID14;
        // pentahedra
        if ( etype == "PENTA_6" )  return ISOPARAMETRIC_LINEAR_PRISM;
        if ( etype == "PENTA_15" ) return ISOPARAMETRIC_QUADRATIC_PRISM15;
        if ( etype == "PENTA_18" ) return ISOPARAMETRIC_QUADRATIC_PRISM18;
      }

    // STRAIGHT-SIDED ELEMENTS WITH A GLOBAL COORDINATE SYSTEM

    // bar elements
    if ( etype == "BAR_2" )    return LINEAR_BAR;
    if ( etype == "BAR_3" )    return QUADRATIC_BAR;
    // trapezoidal elements (quads)
    if ( etype == "QUAD_4" )   return LINEAR_RECTANGLE;
//    if ( etype == "QUAD_4_X" ) return BARYCENTRIC_LINEAR_QUADRILATERAL;
//    if ( etype == "QUAD_8" )   return QUADRATIC_QUADRILATERAL;
//    if ( etype == "QUAD_8_X" ) return BARYCENTRIC_QUADRATIC_QUADRILATERAL;
//    if ( etype == "QUAD_9" )   return QUADRATIC_QUADRILATERAL9;
    // triangles
    if ( etype == "TRI_3" and dim == 3U ) return LINEAR_TRIANGLE3D;
    if ( etype == "TRI_3" and dim == 2U ) return LINEAR_TRIANGLE;
    if ( etype == "TRI_3_X" )  return BARYCENTRIC_LINEAR_TRIANGLE;
    if ( etype == "TRI_6" )    return QUADRATIC_TRIANGLE;
    if ( etype == "TRI_6_X" )  return BARYCENTRIC_QUADRATIC_TRIANGLE;
    // tetrahedral elements
    if ( etype == "TETRA_4" )  return LINEAR_TETRAHEDRON;
    if ( etype == "TETRA_10" ) return QUADRATIC_TETRAHEDRON;
    // hexahedral elements
    if ( etype == "HEXA_8" )   return LINEAR_CUBOID;
//    if ( etype == "HEXA_20" )  return QUADRATIC_HEXAHEDRON20;
//    if ( etype == "HEXA_27" )  return QUADRATIC_HEXAHEDRON27;
    // pyramids
//    if ( etype == "PYRA_5" )   return LINEAR_PYRAMID;
//    if ( etype == "PYRA_13" )  return QUADRATIC_PYRAMID13;
//    if ( etype == "PYRA_14" )  return QUADRATIC_PYRAMID14;
    // pentahedra
//    if ( etype == "PENTA_6" )  return LINEAR_PRISM;
//    if ( etype == "PENTA_15" ) return QUADRATIC_PRISM15;
//    if ( etype == "PENTA_18" ) return QUADRATIC_PRISM18;

    cout <<"\nANSYS_ElementSpecifications::CSMP_TypeFrom_ANSYS_Type: ";
    cout <<" Element type is unknown to CSMP ";
    cout <<"\nType string: "<< etype << endl;
    
    return UNKNOWN;
    
 } // end CSMP_TypeFromANSYSType (char* s, bool)


} // end namespace csmp
