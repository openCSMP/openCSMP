#ifndef ANSYS_INTERFACE_H
#define ANSYS_INTERFACE_H

#include "CSMP_definitions.h"
#include "ANSYS_ElementSpecifications.h"
#include "CSMP_ElementSpecifications.h"

namespace csmp {

class ModelTopology;
class ANSYS_ElementSpecifications;
template<size_t> class VSet;
template<size_t> class Model;

/**
 
@brief File based interface between ANSYS and CSMP.

To provide a file-based interface between the TETRA meshing tool from
ANSYS Technologies Inc. and CSMP which also allows to extract parts of the
input mesh selectively, assign material properties, boundary conditions
and to rigorously test the incorporated mesh for its compliance with
CSMP requirements.
 
@author S.K. Matthai
@author S. Geiger
@author S.G. Roberts
@date 2001
 
@section design Design Intent
 
File reading has been broken up into modules such that with possible 
changes only small bits of code must be recreated. Also files have been
separated into text and binary input, such that the user-readable part
is readily accessible in textformat while the bulk of the data reside
in an efficiently readable binary file.
 

@section applicability Applicability

To files which have been written with a specific output module created
specifically for CSP by ANSYS technologies (Ulrike Wolf, Berlin).
 
 
@section participants Participants

The ANSYS interface uses the classes VSet, FEM_Primitive and 
ModelTopology to store the information in which has been retrieved 
from file. 
 
 
@section implementation Implementation

To keep things confidential, this is all in German ! 
 
2 Dateien, eine Text, die andere binaer (vorerst nur im Textformat zum Testen). 
Beide mit dem gleichen Namen aber die erstere mit der Endung '*.asc' und die
zweite mit der Endung '*.dat'. 

Die Textdatei ('*.asc'):

- der Modellname (eine Zeile)

- Kommentarzeile (ehemaliger header): 1. Name der Datei, 2. wann sie generiert
  wurde (alles in einer Zeile)

- die Anzahl der mit Namen bezeichneten Objekte

- eine Tabelle aller Objekte und ihre Charakteristika wie folgt:

@code
  Objektname   Elementtyp   Materialkennzahl   Anzahl-Elemente
  LIV          TETRA_4      3                  760
  BOXS         TRI_3        0                  2030
  CURVES       BAR_2        0                  20
@endcode

Die Materialkennzahl wird benoetigt um Objekte mit verschiedenen Teilvolumina zu 
ermoeglichen.
  
- Listen der Elemente (numeriert 0...n-1), die die verschiedenen Objekte bilden,
  angefuehrt vom Objektnamen und Elementzahl z.B.:

@code
  LIV_SMALL 6
  12 13 14 15 16 17
@endcode

Diese Teilinformationen der Gesammtausgabe werden immer so wenig Platz einnehmen, 
dass man sie bequem in einem Textfile abspeichern kann. 
Ausserdem gewinne ich auf diese Weise schnell einen Ueberblick darueber
was in dem output file noch fehlt, bzw. was ich vergessen habe in ANSYSTETRA zu benennen.
Auch die Grenzflaechen quaderfoermiger Modelle koennen hier entsprechend der 
CSP Boundary flags als TOP, BOTTOM, LEFT, RIGHT, FRONT, BACK benannt werden.

Alle anderen Informationen lassen sich dann ganz einfach
als 'arrays' in die binaere Datei abspeichern. Diese arrays muessen aber mit der
Angabe der Anzahl ihrer Elemente beginnen damit ich sie effizient lesen kann. Dies
gillt auch fuer die vorlauefige Textversion dieses Files.


Die binaere Datei ('*.dat') enthaelt: 


1. Basic node coordinate data (double64) (wie bereits gemacht):
-------------------------------------------------------------

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

2. Node flags (int) (wie bereits gemacht aber die Kanten muessen noch 
   benannt werden)
   ------------------------------------------------------------------

PBFLAGS: Same record length as the record set above but in int32 format.

- if a node lies on a model boundary, a negative integer value
  is assigned to identify that boundary uniquely in terms of the 
  named surface it belongs to. If the node is not located on a 
  model boundary, its int value is set to zero. If lines or curves were 
  specified on the model boundaries, edge nodes along these lines are 
  flagged by unique negative integers by analogy to the identifiers for
  the edges in the ANSYS topology.
  
  A boundary and edge numbering scheme for box-shaped models in CSP has
  the following values (adhering to a righthand-rule coordinate system
  where x points to the right, y points upward and z points to the front:

hier die Kantenbezeichnungen:

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
    
3. Nodal boundary conditions (double64) (wie bereits gemacht)
-----------------------------------------------------------

PBVALS: Analogous to the record sets above but in double64 format:

- Corresponding to the non-zero flags from above, where a node is 
  located at a model boundary, a boundary value as previously assigned 
  to the ANSYS topological model can be stored in this data record. 
  If no value was assigned, this record will just hold zeroes.


4. Finite element types used for each element (unsigned int)
------------------------------------------------------------

PELMT (0...e-1): Record is preceded by an unsigned long indicating 
record size (=number of finite elements in the mesh). 'pelmt'
consists of a single record of (all as unsigned longs):

- element-type identifiers (unsigned int) for each 
  finite element in the mesh. How they correspond to element
  type strings is specified further below. 
  
  Dieser Eintrag vermeidet Doppeldeutigkeiten, die dadurch entstehen
  koennten, dass es verschiedene Elementtypen mit der gleichen Anzahl
  von Knotenpunkten gibt.   


5. Mesh connectivity data (unsigned long) (wie bereits gemacht)
---------------------------------------------------------------
(Hier muss aber jeder einzelne record mit seiner groesse beginnen)

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
    
    
6. Materialzugehoerigkeiten fuer alle Elemente 
----------------------------------------------
    
PMTRL: preceded by number of elements, associates each element with the 
material flag from the '*.asc' file such that properties can be assigned.  

Kommentare: Beginnen immer mit dem Gartenzaun (# sign) nachdem alles weitere
ignoriert wird.
 
@section examples Application Examples
 
Human readable geometry descrition file (*.asc): 

@code
icem_example1
'icem1.asc' generated by SKM on 8/6/2001
7 # Anzahl Families
# Objektname   Elementtyp   Materialkennzahl   Anzahl-Elemente
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

Mesh connectivity file (*.dat) in text format: 

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
@endcode

*/
class ANSYS_Interface {

public:

    explicit ANSYS_Interface( bool create_isoparametric_element_mesh=true);

    ~ANSYS_Interface();
    
    /// Reading function

    template<size_t dim>
    void Read_ANSYS_Mesh( const std::string& mesh_file_set,
                          VSet<dim>&,
                          ModelTopology&,
                          bool binary_input_file,
                          bool irregular_mesh );

    /// turn this option on or off
    void AssignMaterialPropertiesInteractively( bool ass=true ); 

  protected:

    /// ascii format
    template<size_t dim>
    void ReadMeshASCII( const std::string& meshfile,
                        VSet<dim>&,
                        ModelTopology&,
                        bool irregular_mesh );

    /// binary format
    /// use the binary interfaces for efficient reading of large files
    template<size_t dim>
    void ReadMeshBinary( const std::string& meshfile,
                         VSet<dim>&,
                         ModelTopology&,
                         bool irregular_mesh );

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

    void Clear();

    bool  interactive_property_assignment_;
    bool  isoparametric_;

    // storage for geometric objects and their element types
    //  objectname  elmttype  identifier (volumes, surfaces, curves in file) 
    std::string                                       file_header_;
    std::string                                       boundary_tag_;
    typedef ANSYS_ElementSpecifications               elmt_specs;       ///< not local variable but shorthand
    typedef CSMP_ElementSpecifications                csmp_elmt_specs;
    //          region name   type of elements
    std::multimap<std::string,std::string>            object_specs_;
    //          region name   element ids
    std::multimap<std::string,std::vector<size_t> >   object_elements_;
};

/// ANSYS - to - CSMP FE's
template<size_t dim>
void Convert_ANSYS_To_CSMP_FiniteElementTypes( VSet<dim>&,bool);
void Convert_ANSYS_To_CSMP_FiniteElementTypes( std::multimap<std::string,std::string>&,bool,size_t);

class ANSYS_ModelSettings
{
public:

    ANSYS_ModelSettings( const std::string& mesh_file_prefix );
    ~ANSYS_ModelSettings();

    ANSYS_ModelSettings( const ANSYS_ModelSettings& );
    ANSYS_ModelSettings& operator=( const ANSYS_ModelSettings& );

    void MeshSetup( const std::string& regions_file_prefix,
                    bool irregular_mesh    = true,    /* true = non-box shaped model, false = box shaped model */
                    bool binary_file       = true,    /* true = binary, false = ascii */
                    bool create_boundaries = true );  /* true = create boundaries around model, false = do not create boundaries */

    void MeshSetup( const std::set<std::string>& regions,
                    bool irregular_mesh    = true,    /* true = non-box shaped model, false = box shaped model */
                    bool binary_file       = true,    /* true = binary, false = ascii */
                    bool create_boundaries = true );  /* true = create boundaries around model, false = do not create boundaries */

public:

    std::string mesh_file_prefix_;    // name of model or main mesh file
    std::set<std::string> regions_;   // name of regions to be used
    bool irregular_mesh_;     /* true = non-box shaped model, false = box shaped model */
    bool binary_file_;        /* true = binary, false = ascii */
    bool create_boundaries_;  /* true = create boundaries around model, false = do not create boundaries */
};

} // end namespace csmp

#endif
