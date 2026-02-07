#ifndef CSMP_MODEL_H
#define CSMP_MODEL_H

#include "CSMP_definitions.h"
#include "Box.h"
#include "VSet.h"
#include "PropertyDatabase.h"
#include "FiniteElementManager.h"
#include "MeshManager.h"
#include "RegionInterface.h"
#include "BoundaryInterface.h"
#include "SplitBoundaryInterface.h"
#include "LocalVariableStorage.h"

namespace csmp {

class PropertyConstraints;
class ModelTopology;
class Standard_IO_Handler;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class Element;
template<uint32_t> class VSet;
template<uint32_t> class Interrelation;
template<uint32_t> class Visitor;
template<uint32_t> class FiniteVolumeStencilManager;
template<typename> class FEM_Data;

template<uint32_t,template<uint32_t> class> class PDE_Integrator;
template<uint32_t,template<uint32_t> class> class PDE_Integrator_UoM;


/**
@brief Model the playground for the physics of interest

@author S.K. Matthai
@author Stephen G. Roberts
@date 1999

@section motivation Motivation

Finite-element models typically are characterized by arrays
of data that need to be accessed via complex indexing operations and whose
interconnectedness is not explicit. Also it is relatively error
prone to change the size of these arrays at runtime and to update the
connectivity of represented finite elements. Tohe design of the Model
class was therefore motivated by the need to hide this complexity from the
user such that he/she can focus their undivided attention on the geological
complexity of the problem at hand.


@section design Design Intent

The design intent was to have an object representation
of the mass and space occupancy of a real world system of interest.
One should be able to work with and observe this system through the interfaces
of the Model object.


@section applicability Applicability

A Model representation can be built for any geometrical model which
represents a finite-element discretization of a real-world system.
Computations on a Model are carried out through its interface Apply().
Through this interface you can carry out finite-element computations
specified by Algorithms or own types derived from the classes Algorithm,
and 'Algorithm'.  You can also calculate interrelations among variables
(including the dependent variables that are computed at each timestep).
You do this with your own classes that you derive from the Interrelation
base class.  An example for such an Interrelation-derived calculation would
be to calculate a permeability that is dependent on fluid pressure.  Any
derived algorithm or interrelation can also be restricted to Regions of
Elements inside the Model.


@section structure Structure

The Model is an agglomerate of objects representing the finite-element mesh,
implemented as a connectivity scheme and a storage scheme for physical
variables that are assigned to nodes, constraint points, or elements
themselves.  Regions of elements that make up geological entities in the
model are referred to as Region objects.

The Model contains a class hierarchy of Region, Element, IntegrationPoint,
and Node objects. Their connections reflect the connectivity of the mesh and
they are handled internally by the MeshManager. The complementary
MemoryManager object manages the storage for the distributed physical
variables of a computation. The PropertyDatabase object keeps track of the
existing properties, their storage locations and specifications.


@section participants Participants

The Model functionality is instrumentalized through instances of Node,
IntegrationPoint, and Element classes managed by the MeshManager;
basic CSMP variables managed by the MemoryManager;
a PropertyDatabase, and a map of Regions.


@section collaborations Collaborations

In a typical CSMP simulation, the Model collaborates with Algorithm,
Interrelation, and Visitor objects via its Apply() and Accept() interfaces,
respectively. For data transfer, the Model can also exchange properties
with a FiniteDifferenceGrid object using the Element interpolation functions.


@section consequences Consequences

A Model encapsulates the discretization of a geological object with
associated properties. It thereby creates an interface to this object which
allows to carry out computations on the object, modify property values,
input and output these to other tools, and to address and manipulate sub-
regions of the object. Sub-regions can be identified on the
basis of property values and they can be associated with a name. Once this
is done, most of the Model interface can also be used for
individual Regions.


@section implementation Implementation

The public interfaces of the Model invoke a data-access process for
the discretized real-world model. The key steps in this process are:

(1) The Property Database is queried for the
variable specifications of a variable defined as an input string (for
instance "permeability"). The database returns a csmp::Index that uniquely
identifies the variable for efficient access in the computations.

(2) The placement of the variable (Node, IntegrationPoint or Element)
determines the depth of the search for the variable in the hierarchy of the
mesh in the Model object.

(3) Once a variable location object (Node, IntegrationPoint, or Element)
is found, its ID is used to retrieve the value or status of the variable from the
MemoryManager via methods like Read() or Store().  These methods are
overloaded to retrieve Scalar-, Vector-, and TensorVariable<dim> instances.

(4) Within the MemoryManager object the variables live in
STL vectors that were instantiated for the specific variable types.  The length
of such vectors depends on the placement of the variables.  Originally the
MemoryManager object is build for the variables that were specified in
the variable database from your input file.  When you create a new variable at
runtime (using PropertyHandle objects), a new vector is build and inserted
for this variable. The addresses of other variables remain valid in
this process. Thus, you can efficiently create new variables at runtime.


@section examples Application Examples

In the following example, the element type LinearTriangle is set as
default finite-element. A Triangulator mesh generator object is used by
the Model constructor to generate a mesh of triangular finite-elements.
This mesh is build from regular-gridded input data for which the user is
prompted and the mesh in used further to buid a Model object named
'model'. In this process, the default variable text file 'CSP_variables.txt'
is read to initialize the variable database. The porosity of the new model is
set to a uniform value of 3% and DIRICH(let) boundary conditions are applied
at the model top (cross-sectional model). Now, instances of subclasses of
Algorithms, Interrelations and Visitors are applied by 'passing' them
to the Model. This invokes global (mesh-wide) and local (finite-element
restricted) computations. Finally the variable 'temperature' is output
to the HDF file 'computed_temperature'.


@code
Triangulator mesher;

Model   model( "CM-simulation", mesher );

model.InputUniformValueWhere ( PLAIN, "porosity", 0.03 );

model.AssignBoundaryValues( TOP, "temperature", DIRICH, top_T, top_T );

model.Apply( interrelation_subclass );

model.Apply( csp_algorithm_subclass );

model.OutputDataToHDF ( "computed_temperature", "temperature" );
@endcode


@todo (3) Put Apply(PDE_Int) back into domain classes (from Model to Boundaries etc...)
@todo (3) Replace references to groupMap_ and uniqueGroupMap_ in Model.cpp by corresponding Interface functionality
@todo (3) Test binary IO of SplitBoundaries
*/
template<uint32_t dim>
class Model : public RegionInterface<dim, Model>,
              public BoundaryInterface<dim, Model>,
              public SplitBoundaryInterface<dim, Model>,
              public LocalVariableStorage<dim, Model> 
{

public:
  // class Model is not copy constructable

  /// constructs model with subdomains (Region, Boundary, SplitBoundary), variables file name with extension "*-variables.txt" where * is the name of the model
  Model( ModelTopology& mesh_topology, VSet<dim>& mesh, const char* var_file, bool use_regions_file_if_any_to_select_domains_to_keep );

  /// Reconstructor:  reads model from set of CSMP's native binary files
  explicit Model( const std::string& binaryFiles );

  /// Reconstructor:  reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
  Model( const std::string& binaryFileName, const std::string& variable_txt_file );

  /// Reconstructor: reads model from set of CSMP's native binary file, but only reading the specified subset of variables
  Model( const std::string& binaryFileName, const std::set<std::string>& subset_variables );

  /// using the supplied polygonal data constructs unnamed single-domain model without regions or boundaries
  Model( VSet<dim>& polygonal_dataset, const char* var_file );

  /// constructs purely topological unnamed single-domain model without regions, boundaries nor variable storage
  explicit Model( VSet<dim>& );

  /// destructor that needs to be overloaded when a subclass is derived from model
  virtual ~Model();

  /// returns const (read-only) reference to database object where all variable access and type specifications are stored
  const PropertyDatabase<dim>&  Database() const;

  /// returns reference to database object where all variable access and type specifications are stored
  PropertyDatabase<dim>&  Database();

  /// read-only access the mesh container (nodes, elements, faces, interfaces and associated variable storage)
  const MeshManager<dim>&  Mesh() const;

  /// access the mesh container (nodes, elements, faces, interfaces and associated variable storage)
  MeshManager<dim>&  Mesh();

  /// read-only access to the finite element types that are needed to support the current mesh
  const  FiniteElementManager&  FE_Manager() const;
  
  /// prompts the finite volume manager to connect the elements and faces with finite volume stencils 
  void InstantiateFiniteVolumes();

  PLACEMENT Placement() const { return MODEL; }


  // -----------------------------------------------
  // Binary input/output
  // -----------------------------------------------

  /// writes entire model with associated properties to disk; non-constant because this involves region creation; not const because mesh is updated
  void OutputToBinaryFile( const char* ); // not const, because cells might get reshuffled

  /// reads model written by OutputToBinaryFile() including all associated properties; if subset of variables is not empty only these will be read
  void InputFromBinaryFile( const char* model_name, const std::set<std::string>& subset_variables );
  
  /// computes deques of numbered Node, Element, Face and InterFace objects, and outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
  void OutputMeshTo( VSet<dim>&, bool get_indices_from_stored_variables=false );

  /// writes discretised variable to generic variable container
  template<class Var>
  void OutputVariableTo( const char* var, FEM_Data<Var>& ) const;

  /// initialised the input property using the data supplied via the VSet
  template<class T>
  void InputVariableFrom( const char* input_prop, const FEM_Data<T>& );
  
  /// inputs all discretised variables stored in the supplied VSet into the model
  void InputVariablesFrom( const VSet<dim>& );

  // ------------------------------------------------------------------------
  // Property interface
  // ------------------------------------------------------------------------

  /// inserts (if new) variable into the database and creates storage for it on the entities where it shall be discretized
  csmp::Index  CreateProperty( const char* new_prop, const char* notation, const char* unit,
                               VARIABLE_TYPE type = SCALAR, PLACEMENT place = NODE,
                               uint32_t vsize = 1, double vmin = -1.0e+30, double vmax = 1.0e+30,
                               std::string usage = "???" );

  /// deletes property from the database and the distributed containers all across the model
  void DeleteProperty( const char* property );

  /// renumbers everything, using stored  "node number" etc. if available, else MeshManager is prompted to created unique numbering
  void IndexByPropertyValues();

  /// sets the values of the distributed variable all across the model; to enter scalar value use makeScalar(flag,value) helper function
  template<class T>
  void InputPropertyValue( const char* input_prop, const T& value );

  /// sets the value of the property on those sites (nodes, elements) whose AtBoundary() function matches the BOX_BOUNDARY value
  template<class T>
  void InputBoundaryValue( BOX_BOUNDARY, const char* input_prop, const T& value );

  /// uses the vector of variable flags (ANY, DIRICH etc.) to change the status of the target variable at the given BOX_BOUNDARY
  void InputBoundaryFlags( BOX_BOUNDARY, const char* property, const std::vector<VARIABLE_FLAG>& flags );

  /// the supplied 4 scalar node variable values (2 in 2D) are linearly interpolated across the square (line) boundary of the box-shaped model
  void InterpolateBoundaryValues( BOX_BOUNDARY, const char* input_prop, const std::vector<ScalarVariable>& bvalues );

  /// the supplied 4 vector node variable values (2 in 2D) are linearly interpolated across the square (line) boundary of the box-shaped model
  void InterpolateBoundaryValues( BOX_BOUNDARY, const char* input_prop, const std::vector<VectorVariable<dim> >& bvalues );

  /// replaces the values of the target propery (to) with the values of property (from); both variables must have the same type and placement
  void CopyReplace( const char* from, const char* to );

  /// computes the (constant valued) first spatial derivative of the node property and assigns it to the element property / barycenter of the element
  bool CopyGradientOfProperty_A_To_B( const char* node_prop, const char* element_prop );

  /// linearly interpolates the value of the integration point property to the barycentre of element; result is stored as element property
  void InterpolateIntegrationPointToCellProperty( const char* ipoint_prop, const char* eprop );

  /// interpolates node property values to the barycentre of element and stores results in element property
  void InterpolateNodeToCellProperty( const char* nprop, const char* eprop, bool verbose = true );

  /// interpolates distributed node property values the quadrature points of numerically integrated finite elements
  void InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* ipoint_prop );

  /// linearly extrapolates the values of the integration point variable to the element nodes where an averaging with the neighbor elements is performed
  void ExtrapolateIntegrationPointToNodeProperty( const char* ipoint_prop, const char* eprop );

  /// piecewise constant element property values are extrapolated to nodes using a choice of averaging schemes (1/distance vs. element-volume weighted)
  void ExtrapolateCellToNodeProperty( const char* eprop, const char* nprop, bool by_distance = true );

  /// changes the flags of the target variable all across the model
  void ChangePropertyStatus( const char* input_prop, VARIABLE_FLAG new_status );

  /// where the values of the target property are within the given range the flag of the target variables are changed to the new status
  void ChangePropertyStatusWhere( const char* var, double min, double max, VARIABLE_FLAG new_status );

  /// returns the opposite corners of the bounding box that encloses the model
  void MinMaxCoordinates( Point<dim>& xyz_min, Point<dim>& xyz_max ) const;

  /// returns the value range of the target property within the entire model
  void MinMaxOf( const char* prop, double& min, double& max ) const;
  
  /// checks the value range of a variable against the range specified in its Parameter record in the PropertyDatabse
  bool IsWithinRange( const char* prop, const char* model_subdomain="Model" ) const;

  /// permits to transfer node coordinate components to the target scalar node variable; char options are 'x', 'y', 'z'
  void AssignNodeCoordinatesTo( const char* scalar_variable, char coord ); // x, y, z

  /// transfers  node coordinates to  target VectorVariable<dim>
  void AssignNodeCoordinatesTo( const char* vector_variable );

  /// shifts the node coordinates by an amount that is determined by the value of the target VectorVariable<dim>
  void MoveNodeCoordinatesBy( const char* vector_variable );

  /// 'characteristic' options are: volume, inner radius, and aspect ratio
  void AssignCellCharacteristicsTo( const char* characteristic, const char* var );

  // ------------------------------------------------------------------------
  // Interrelations, Visitors and Algorithms
  // (to apply these to specific regions or boundaries access these directly)
  // ------------------------------------------------------------------------

  /// prompts the region-by-region, element-by-element or node-by-node calculation of the target variable via the interrelation
  void Apply( Interrelation<dim>& relation, const char* region = "Model" );

  /// support of the visitor pattern giving visitors access to the model
  void Accept( csmp::Visitor<dim>& );

  /// application of integration scheme to model, subregions thereof or boundary or split-boundary objects
  void Apply( PDE_Integrator<dim,Element>&, bool debug = false );
  void Apply( PDE_Integrator<dim,Face>&, bool debug = false );

  /// application of integration scheme to a particular region, boundary of split-boundary identified by name
  void Apply( PDE_Integrator<dim,Element>&, const char* region_name, bool debug = false );
  void Apply( PDE_Integrator<dim,Face>&, const std::string& boundary_name, bool debug = false );

  // ----------------------------------------
  // Screen output
  // ----------------------------------------

  /// the name of the computational model
  const char* Name() const;
  void        Name( const char* );

  /// console output
  void OutputVariableToScreen( const char* prop ) const;
  void Out() const;

  void Verbose( bool verbose );
  bool Verbose() const;

protected:

  /// default constructor that is custom-made
  Model();
  
  /// to construct model  as a base class to a derived model built from an external dataset using a variable file in ASCII format
  explicit Model( const char* complete_variables_file_name );

  /// builds with model without specifically named subdomains other than "Model".  If model is box-shaped corresponding boundary flags will exist.
  void Initialize( VSet<dim>& );

  /// builds model treating any topologic entities from regions file as unique regions; later converts correspondingly labelled ones into boundaries and split boundaries
  void Initialize( const char* regions_file_prefix, ///< use to select regions in ANSYS or other mesh
                   ModelTopology&, VSet<dim>& );    ///< must not contain any Face or InterFace objects which will be built by this method

  /// builds model  with regions, boundaries, and splitboundaries as identified by "BOX_BOUNDARY", "boundary" or "splitboundary" strings in the domain names
  void Initialize( ModelTopology&, ///<  must store Region, Boundary, and SplitBoundary objects 
                   VSet<dim>& );   ///< must store corresponding Element, Face and InterFace objects matching the idx integers stored in topology

  void InitializeLocalVariableStorage();
  bool UpdateSubdomainPropertyStorage();

  static std::string BinaryVsetFileName( const char* base_file_name );
  static std::string BinaryRegionsFileName( const char* base_file_name );
  static std::string BinaryBoundariesFileName( const char* base_file_name );
  static std::string BinarySplitBoundariesFileName( const char* base_file_name );
  static std::string BinaryVariablesFileName( const char* base_file_name );

private:

  /// prevent accidential copy construction of large Model object
  Model( const Model& ) = delete;
  Model& operator=( const Model& ) = delete;

  /// checks wether elements have their correct neighbors and are in the expected model domains; @return number of major errors encoutered.
  int32_t CheckElementConnectivity();

  std::string            model_name_;       ///< name of simulation model
  PropertyDatabase<dim>  database_;         ///< where variable specifications are stored
  MeshManager<dim>       mesh_manager_;     ///< stores mesh: all Node, Element, Face, InterFace objects
  bool                   verbose_ = false;  ///< for detailed screen output todo: replace with global verbose singleton
};

// SUPPORTING FUNCTIONS

/// attempts to return the spatial dimension of the model stored in the file (1-3D)
uint32_t spatialDimensionOfModel( const char* csmp_binary );

/// returns the extent of the model in the x,y,z dimensions and reports this back as a string
std::string  boundingBox( const Model<3U>& sg, double& dim_x, double& dim_y, double& dim_z );

/// returns intermediate (true) or maximum (false) model dimensions
template<uint32_t  dim>
double  printModelDimensions( const Model<dim>&, bool intermed_or_max = false );

/// calculates the center of gravity of the model
template<uint32_t  dim>
Point<dim>  centerOfGravity( const Model<dim>& );

/// prints range to screen; returns either min(arg=false) or maximum variable value (default)
template<uint32_t  dim>
double  printRangeOfVariable( const Model<dim>&,
                              const char* var, bool print_maximum = true );

/// prints range of target variable in model to screen and logs it to IO handler
template<uint32_t  dim>
double  printRangeOfVariable( const Model<dim>&,
                              Standard_IO_Handler& io, const char* var,
                              bool max_instead_of_min = true );

/// prints range of target variable within specific model subdomain (Region, Boundary or SplitBoundary)
template<uint32_t  dim>
double  printRangeOfVariable( const Model<dim>&,
                              const char* subdomain, const char* var, bool max_or_min = true );

/// prints range of target variable within specific model subdomain and logs it to IO handler
template<uint32_t  dim>
double  printRangeOfVariable( const Model<dim>&,
                              Standard_IO_Handler&,
                              const char* group, const char* var,
                              bool max_instead_of_min = true );

/// prints min/max values stored in supplied vector
void printRangeOf( const std::vector<std::pair<double, double> >& );

/// prints min/max values stored in supplied vector of vectors
void printRangeOfVectorOfVectors( const std::vector<std::vector<double> >& );

/// converts the VARIABLE_FLAG flag(s) of a variable into integer values stored in its number part
template<uint32_t dim>
void flagToNumber( Model<dim>&, const char* variable );

/// converts the VARIABLE_FLAG flag(s) of first variable into integer values stored in the second variable
template<uint32_t dim>
void flagToNumber( Model<dim>&, const char* flag_variable, const char* number_variable );

/// converts BOX_BOUNDARY flag to variable values stored in specified variable
template<uint32_t dim>
void boundaryFlagToNumber( Model<dim>& model, const char* box_boundary_var );

/// converts the model topology flags of the nodes (TopoType) into node variable values that can be visualised
template<uint32_t dim>
void topoTypeToNumber( Model<dim>&, const char* TopoType_variable );

/// using random number generator, adds percentage of Gaussian noise to variable values
template<uint32_t dim>
void randomPerturb( Model<dim>&, const char* prop, double by_percent_of_max_value );

/// in target region, element variable is extrapolated to node and back as many times as indicated by n_smoothing_cycles
template<uint32_t dim>
void smoothElementVariable( Model<dim>&, const char*, const char* element_var, const char* temp_node_var, uint32_t n_smoothing_cycles );

/// replaces no-data values of target variable with nearest-neighbor values until there are none left, by default NAN's are no-data values
template<uint32_t dim>
void nearestNeighborFill( Model<dim>&, const char* target_region, const char* variable, double no_data_value );

/// imposes either an upper- or a lower limit on the variable in the region of interest
template<uint32_t dim>
void imposeLimitOn( Model<dim>& model, const char* region, const char* variable, bool upper_limit, double limit_value );

/// mapping node coordinates to a node variable
template<uint32_t dim>
void assignNodeCoordinatesTo( Model<dim>& sg, const char coordinate, const char* nodal_variable );

/// for Triangulator meshes: for meshes created from pixels, finds outlier triangles (3 nodes at region boundaries) and flips their property values
void stripDomainEdgesFor( Model<2U>&, const char* el_prop );


// MESH MANAGEMENT UTILITIES

/// initialise TOPOTYPE flags and returns # nodes on faces, edges, and vertices (primitive types of a B-REP)
template<uint32_t dim>
tuple<size_t,size_t,size_t> initialise_BREP_TopologyFlags( Model<dim>& );

/// compares the mesh connectivity in the model with that of the input vset; returns true if both have the same
template<uint32_t dim>
bool compareConnectivity( const Model<dim>&, const VSet<dim>& );



} // end namespace csmp

#endif
