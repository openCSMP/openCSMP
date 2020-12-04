#ifndef SKUA_FINITE_ELEMENT_MESH_INTERFACE_H
#define SKUA_FINITE_ELEMENT_MESH_INTERFACE_H

#include "CSMP_definitions.h"
#include "CSMP_ElementSpecifications.h"

namespace csmp {

class ModelTopology;
class SKUA_ElementSpecifications;
template<size_t> class VSet;
template<size_t> class Model;

/**
 
@brief File based interface between SKUA and CSMP.

To provide a direct file-based interface between the SKUA's  FiniteElementMesher and CSMP,
allowing to read parts of the input mesh selectively, assign material properties, boundary conditions
and to rigorously test the incorporated mesh for its compliance with
CSMP requirements.
 
@author S.K. Matthai
@author Thomas Jerome
@date 2020
 
@section participants Participants

Part of  the file format of this interface is modelled on the SKUA interface using the classes VSet, FEM_Primitive and 
ModelTopology to store the information retrieved from file. 
 
 
@section implementation Implementation

Two files: one human readable text (.asc = ASCII) file and a binary (.dat) file
that can be transferred across platforms (Windows OSX, linux etc/).
Both files have the same name (model name) and are distinguished only by their extension.

The regions and material file ('*.asc'):

- model name (line break)

- Comment line (header): 1. filename, 2. date of generaton

- Number of labelled objects

- Table with objects, element types, material ID values (rock types), number of elements (cells)

@code
  Object   element-type   material ID   number-of-elements
  LIV          TETRA_4      3                  760
  BOXS         TRI_3        0                  2030
  CURVES       BAR_2        0                  20
@endcode

The material ID facilitates the association of 'rock types' with model regions.
  
- List of elements making up each region (numbered 0...n-1),  preceded by region name and 
number of elements of region.

@code
  UNIT_1 6
  12 13 14 15 16 17
@endcode

The binary file ('*.dat') contains: 

1. NODE COORDINATE DATA (double64):

PX, PY, PZ records: preceded by a single unsigned long indicating the size of
these records which (individually) contain (with or without line breaks):

px: -> x-coordinates of the nodes  
py: -> y-coordinates of the nodes
pz: -> z-coordinates of the nodes

(all entries in these records have the format double64 and
the node-IDs are implicit from 0 to nodes-1):

Such records could be written as shown in the C++ code example
below:

@code
double64*        ptr; // data pointer
...
unsigned long  n0(0), n1(1);

if ( (records=px.size()) > 0 && (ptr=const_cast<double64*>(px.Data())) != NULL ) 
  {
     fwrite( (void*) &records, sizeof(unsigned long), 1, fp );
     fwrite( (void*) ptr, sizeof(double64), records, fp );
  }
else fwrite( (void*) &n0, sizeof(unsigned long), 1, fp );
@endcode

2. NODE FLAGS (int)

PBFLAGS: Same record length as the record set above but in int32 format.

- if a node lies on a model boundary, a negative integer value
  is assigned to identify that boundary uniquely in terms of the 
  named surface it belongs to. If the node is not located on a 
  model boundary, its int value is set to zero. If lines or curves were 
  specified on the model boundaries, edge nodes along these lines are 
  flagged by unique negative integers by analogy to the identifiers for
  the edges in the SKUA topology.
  
  A boundary and edge numbering scheme for box-shaped models in CSP has
  the following values (adhering to a righthand-rule coordinate system
  where x points to the right, y points upward and z points to the front:

@code
#define BACK_BOTTOM       -16 // BACK and BOTTOM
#define BACK_RIGHT        -17 // BACK and RIGHT
#define BACK_TOP          -18 // BACK and TOP
#define BACK_LEFT         -19 // BACK and LEFT
#define BOTTOM_RIGHT      -20 // BOTTOM and RIGHT
#define TOP_RIGHT         -21 // TOP and RIGHT
#define TOP_LEFT          -22 // TOP and LEFT
#define BOTTOM_LEFT       -23 // BOTTOM and LEFT
#define FRONT_BOTTOM      -24 // FRONT and BOTTOM
#define FRONT_RIGHT       -25 // FRONT and RIGHT
#define FRONT_TOP         -26 // FRONT and TOP
#define FRONT_LEFT        -27 // FRONT and LEFT
@endcode
   
    
3. BOUNDARY CONDITIONS APPLIED TO NODES (double64)

PBVALS: Analogous to the record sets above but in double64 format:

- Corresponding to the non-zero flags from above, where a node is 
  located at a model boundary, a boundary value as previously assigned 
  to the SKUA topological model can be stored in this data record. 
  If no value was assigned, this record will just hold zeroes.


4. TYPE OF FINITE ELEMENTS (enum -> unsigned int)

PELMT (0...e-1): Record is preceded by an unsigned long indicating 
record size (=number of finite elements in the mesh). 'pelmt'
consists of a single record of (all as unsigned longs):

- element-type identifiers (unsigned int) for each 
  finite element in the mesh. How they correspond to element
  type strings is specified further below. 
  
  Dieser Eintrag vermeidet Doppeldeutigkeiten, die dadurch entstehen
  koennten, dass es verschiedene Elementtypen mit der gleichen Anzahl
  von Knotenpunkten gibt.   


5. MESH CONNECTIVITY (unsigned long) 

PLIST: Record preceded by an unsigned long indicating record size,
a single record follows containing (all as unsigned longs):

- node ids of each element in
  counter-clockwise order (according to node numbering on reference elements). 
  The enlisted node numbers range from 0...nodes-1 as above.

The elements are numbered implicitly and continuously by their position 
in the plist (there must be no jumps in the element ids).

    
PFVERTS: Record preceded by unsigned long indicating record size, then a
single record (all long, where negative numbers denote model boundaries):
    
- neighbour elements
  of each element in an order defined for each existing element 
  type (see further below, for tetrahedra and triangles, neighbors are
  those elements touching the element's faces opposite of its corner nodes). 
  Neighbor elements are denoted 0...element-1. Where there are no neighbour 
  elements adjacent to the element's faces because the element lies on a model 
  boundary, a negative long is used to identify this boundary surface.  
    
    
6. MATERIAL IDENTIFIER FOR EACH ELEMENT = ROCK TYPE
    
PMTRL: preceded by number of elements, associates each element with the 
material flag from the '*.asc' file such that properties can be assigned.  

Comments: leading hash key (# everything thereafter is ignored) 


7. ELEMENT AND NODE PROPERTIES

Following the material ID record of integers, the property data block starts.
Since the properties have different placements (element, node etc.) and types (scalar, vector, tensor, array, flagged-array),
this information must be specified alongside with the property name.
Furthermore, by contrast with VectorVariable (dim) and TensorVariable (dim x dim), the size of which is defined by the model dimension, dim,
ArrayVariable and FlaggedArrayVariable have a length (number of entries) that is defined at runtime
and therefore has to be specified.
While each ArrayVariable has only one flag, each FlaggedArrayVariable has one flag per value.
For this reason, flag strides and data strides are provided for all variables.

In summary, a complete record for a CSMP variable consists of:

name(string)    placement(enum PLACEMENT)   type (VARIABLE_TYPE)    (array length (int32))
number of data entries (int32)
flag values (int32)
property values (double64)

Since the names can contain whitespace, they are placed in quotation marks.

@section examples Application Examples
 
Human readable geometry descrition file (*.asc): 

@code
example1
'icem1.asc' generated by SKM on 8/6/2001
7 # number of regions that make up the simulation model
# region_name   element_type   material_ID   number-of-elements-in-region
soil			     TRI_3					1					6
rock			     TRI_3					2					1
air				     TRI_3					3					5
left_boundary	 BAR_2					1					2
right_boundary BAR_2					3					2
top_boundary	 BAR_2					3					2
bottom_boundary	BAR_2					1					2
# now the elements which make up each object are listed in sequence
soil 6
2 4 7 8 9 10
rock 1
5
air 5
0 1 3 6 11
left_boundary 2
12 13
right_boundary 2
18 19
top_boundary 2
14 15
bottom_boundary 2
16 17
@endcode

Mesh connectivity file (*.dat) in text ASCII version: 

@code
11 # px, py, pz
0. 5. 10. 0. 5. 10. 2.5 7.5 0. 5. 10.
10. 10. 10. 6. 7. 5. 4. 4. 0. 0. 0.
0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
# pbflags
-13 -5 -14 -2 0 -3 0 0 -8 -4 -10
# pbvals
0. 5. 10. 0. 0. 10. 0. 0. 0. 5. 10.
20 # pelements
3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2
42 # plist
0 4 1  1 4 2  0 3 4  3 6 4  4 6 7  4 7 5  3 8 6
6 8 9  6 9 7  7 9 10  7 10 5  0 3  3 8  0 1  1 2
8 9  9 10  2 5  5 10
42 # pfverts
1 -5 2  3 -5 0  4 0 -2  5 2 7  9 6 4  11 3 5  -3 1 6
8 4 -2  -4 9 7  10 5 8  11 9 -4  -3 6 10
-2 2  -2 7  0 -5  1 -5  -4 8  -4 10  3 -3  11 -3
20 # pmaterial
3 3 1 3 1 2 3 1 1 1 1 3 1 1 3 3 1 1 3 3
1 # element and node properties
"permeability"  "element" "scalar"
1 1
20
1 (ANY) 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
20
1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 
1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12 1.0e-12  
@endcode

*/
class SKUA_FiniteElementMeshInterface {
  public:
      explicit SKUA_FiniteElementMeshInterface( bool create_isoparametric_element_mesh=true );
      ~SKUA_FiniteElementMeshInterface();
      
      /// Reading SKUA mesh from .asc and .dat files; deduces irregular if more than the top surface is warped (IRREGULAR)
      template<size_t dim>
      void Read_SKUA_Mesh( const std::string& mesh_file_set,
                            VSet<dim>&,
                            ModelTopology&,
                            bool binary_input_file );
  protected:
    /// ascii format
    template<size_t dim>
    void ReadMeshASCII( const std::string& meshfile,
                        VSet<dim>&,
                        ModelTopology& );

    /// binary format - use  binary interfaces for efficient reading of large files
    template<size_t dim>
    void ReadMeshBinary( const std::string& meshfile,
                         VSet<dim>&,
                         ModelTopology& );
  private:
    /// reading comments
    bool IsCommentLine( char* str ) const;
    bool SkipPotentialComment( std::ifstream& ) const;
    void AdvancePastCommentLine( std::ifstream& ) const;
    
    /// ASCII file ( regions and corresponding element id's )
    bool ReadTitleASCII( std::ifstream& ifs, std::string& title );
    bool ReadRegionsAndElementTypesASCII( std::ifstream& ifs );

    /// ASCI data file
    template<size_t dim>
    bool ReadNodeCoordinatesASCII( std::ifstream& ifs, VSet<dim>& );
    template<size_t dim>
    bool ReadBoundaryFlagsAndConditionsASCII( std::ifstream& ifs, VSet<dim>& );
    template<size_t dim>
    bool ReadPelementASCII( std::ifstream&, VSet<dim>& );
    template<size_t dim>
    bool ReadPlistASCII( std::ifstream&, VSet<dim>& );
    template<size_t dim>
    bool ReadPfvertsASCII( std::ifstream&, VSet<dim>& );
    template<size_t dim>
    bool ReadPmaterialASCII( std::ifstream&, VSet<dim>& );
    template<size_t dim>
    bool ReadPropertyRecordsASCII( std::ifstream&, VSet<dim>& );

    /// Binary data file
    template<size_t dim>
    bool ReadNodeCoordinatesBinary( FILE*, VSet<dim>& );
    template<size_t dim>
    bool ReadBoundaryFlagsAndConditionsBinary( FILE*, VSet<dim>& );
    template<size_t dim>
    bool ReadPelementBinary( FILE*, VSet<dim>& );
    template<size_t dim>
    bool ReadPlistBinary( FILE*, VSet<dim>& );
    template<size_t dim>
    bool ReadPfvertsBinary( FILE*, VSet<dim>& );
    template<size_t dim>
    bool ReadPmaterialBinary( FILE*, VSet<dim>& );
//    TODO: template<size_t dim>
//    bool ReadPropertyRecordBinary( FILE*, VSet<dim>& );

    /// removes any object and data records potentially held in the interface
    void Clear();

  private:
    // storage for geometric objects and their element types
    //  objectname  elmttype  identifier (volumes, surfaces, curves in file) 
    std::string                                      file_header_;
    std::string                                      boundary_tag_;
    typedef CSMP_ElementSpecifications               csmp_elmt_specs;
    //          region name   type of elements
    std::multimap<std::string,std::string>           object_specs_;
    //          region name   element ids
    std::multimap<std::string,std::vector<size_t> >  object_elements_;
    bool  isoparametric_;
};


/// returns the corresponding CSMP element type, taking into account whether an isoparametric FEM formulation is used
CSMP_FEM_TYPE convertSKUA_ElementType( int32 etype, bool isoparametric );



} // end namespace csmp

#endif
