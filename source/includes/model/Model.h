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
#include <string>

namespace csmp {

class PropertyConstraints;
class ModelTopology;
template<size_t> class Face;
template<size_t> class InterFace;
template<size_t> class Element;
template<size_t> class VSet;
template<size_t> class Interrelation;
template<size_t> class Visitor;
template<size_t> class FiniteVolumeStencilManager;

template<size_t,template<size_t> class> class PDE_Integrator;


template<size_t dim>
class Model : public RegionInterface<dim,Model>,
              public BoundaryInterface<dim,Model>,
              public SplitBoundaryInterface<dim,Model>,
              public LocalVariableStorage<dim,Model<dim> > ///< @todo FIX TEMPLATE-TEMPLATE parameter
{

public:

    /// using the supplied polygonal data constructs unnamed single-domain model without regions or boundaries
    Model( VSet<dim>&, const char* var_file, bool isoparametric=false, bool binaryVariablesFile = false );
  
    /// using the supplied polygonal data constructs unnamed single-domain model without regions, boundaries nor variable storage
    Model( VSet<dim>&, bool isoparametric=false );
  
    /// constructs un-named multi-domain model
    Model( ModelTopology&, VSet<dim>&, const char* var_file, bool binaryVariablesFile = false );
    Model( ModelTopology&, VSet<dim>& );

    /// to read model from set of CSMP native binary files
    explicit Model( const std::string& binaryFiles );

    /// to read model from set of CSMP native binary files
    Model( const std::string& varFile, bool binary );

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

    /// access to the finite element types that are needed to support the current mesh
    FiniteElementManager&  FE_Manager();
  
    /// read only access to low-level finite volume functionality
    const  FiniteVolumeStencilManager<dim>*  FV_Manager() const;

    /// access to low-level finite volume functionality
    FiniteVolumeStencilManager<dim>*  FV_Manager();

    PLACEMENT Placement() const { return MODEL; }

    /// connects the finite-volume stencil pointers of the elements to the stencils after initialising them
    void InstantiateFiniteVolumes();

    // -----------------------------------------------
    // Binary input/output
    // -----------------------------------------------

    // NEW output and input interfaces
    /// writes entire model with associated properties to disk; non-constant because this involves region creation
    void OutputToBinaryFile( const char* ) const;
  
    // NEW
    /// reads model written by OutputToDisk() including all associated properties
    void InputFromBinaryFile( const char* );

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
    csmp::Index  CreateProperty(const char* new_prop, const char* unit,
                                VARIABLE_TYPE type=SCALAR, PLACEMENT place=NODE,
                                size_t vsize=1 , double64 vmin = -1.0e+30, double64 vmax = 1.0e+30,
                                std::string usage = "???");

    /// deletes property from the database and the distributed containers all across the model
    void DeleteProperty( const char* property );

    /// renumbers everything, starting face and interface numbers after element index max; TODO: deprecate
    size_t UpdateIndices() const;

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
    void InterpolateIntegrationPointToElementProperty( const char* ipoint_prop, const char* eprop );

    /// interpolates node property values to the barycentre of element and stores results in element property
    void InterpolateNodeToElementProperty( const char* nprop, const char* eprop, bool verbose = true );

    /// interpolates distributed node property values the quadrature points of numerically integrated finite elements
    void InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* ipoint_prop );
  
    /// linearly extrapolates the values of the integration point variable to the element nodes where an averaging with the neighbor elements is performed
    void ExtrapolateIntegrationPointToNodeProperty( const char* ipoint_prop, const char* eprop );
  
    /// piecewise constant element property values are extrapolated to nodes using a choice of averaging schemes (1/distance vs. element-volume weighted)
    void ExtrapolateElementToNodeProperty( const char* eprop, const char* nprop, bool by_distance = true );

    /// changes the flags of the target variable all across the model
    void ChangePropertyStatus( const char* input_prop, VARIABLE_FLAG new_status );

    /// where the values of the target property are within the given range the flag of the target variables are changed to the new status
    void ChangePropertyStatusWhere( const char* var, double64 min, double64 max, VARIABLE_FLAG new_status );

    /// returns the opposite corners of the bounding box that encloses the model
    void MinMaxCoordinates( Point<dim>& xyz_min, Point<dim>& xyz_max ) const;

    /// returns the value range of the target property within the entire model
    void MinMaxOf( const char* prop, double64& min, double64& max ) const;

    /// permits to transfer node coordinate components to the target scalar node variable; char options are 'x', 'y', 'z'
    void AssignNodeCoordinatesTo( const char* scalar_variable, char coord ); // x, y, z

    /// transfers the node coordinates to the target VectorVariable<dim>
    void AssignNodeCoordinatesTo( const char* vector_variable );

    /// shifts the node coordinates by an amount that is determined by the value of the target VectorVariable<dim>
    void MoveNodeCoordinatesBy( const char* vector_variable );

    /// 'characteristic' options are: volume, inner radius, and aspect ratio
    void AssignElementCharacteristicsTo( const char* characteristic, const char* var );

    // ------------------------------------------------------------------------
    // Interrelations, Visitors and Algorithms
    // (to apply these to specific regions or boundaries access these directly)
    // ------------------------------------------------------------------------
  
    /// prompts the region-by-region, element-by-element or node-by-node calculation of the target variable via the interrelation
    void Apply( Interrelation<dim>& relation, const char* region="Model" );

    /// support of the visitor pattern giving visitors access to the model
    void Accept( csmp::Visitor<dim>& );

    /// application of integration scheme to model, subregions thereof or boundary or split-boundary objects
    void Apply( PDE_Integrator<dim,csmp::Region>&, bool debug=false );
    void Apply( PDE_Integrator<dim,csmp::Boundary>&, bool debug=false );
    void Apply( PDE_Integrator<dim,csmp::SplitBoundary>&, bool debug=false );

    /// application of integration scheme to a particular region, boundary of split-boundary identified by name
    void Apply( PDE_Integrator<dim,csmp::Region>&, const char* region_name, bool debug=false );
    void Apply( PDE_Integrator<dim,csmp::Boundary>&, const std::string& boundary_name, bool debug=false );
    void Apply( PDE_Integrator<dim,csmp::SplitBoundary>&, const std::string& splitboundary_name, bool debug=false );

    // ----------------------------------------
    // Screen output
    // ----------------------------------------

    /// the name of the computational model
    const char* Name() const;
    void        Name( const char* );

    /// console output
    void OutputVariableToScreen( const char* prop ) const;
    void Out() const;

    void Verbose(bool verbose);
    bool Verbose();

  protected:

    Model();

    /// build model from scratch
    void Initialize( bool isoparametric_elements,
                     VSet<dim>& vset,
                     bool create_boundaries,
                     bool non_box_shaped_model);

    /// builds model from scratch including region information from file (this method is used by ANSYS_Model3D) 
    void Initialize( const char* regions_file_prefix,
                     ModelTopology& mesh_topology,
                     VSet<dim>& vset,
                     bool create_boundaries,
                     bool non_box_shaped_model );

    /// builds model from scratch without any region information; the only (unique) region will be 'Model'
    void Initialize( ModelTopology& mesh_topology,
                     VSet<dim>& vset,
                     bool create_boundaries,
                     bool non_box_shaped_model);

    void InitializeLocalVariableStorage();
    bool UpdateSubdomainPropertyStorage();

    static std::string BinaryVsetFileName( const char* base_file_name );
    static std::string BinaryRegionsFileName( const char* base_file_name );
    static std::string BinaryBoundariesFileName( const char* base_file_name );
    static std::string BinarySplitBoundariesFileName( const char* base_file_name );
    static std::string BinaryVariablesFileName( const char* base_file_name );

  private:

    /// prevent accidential copy construction of large object
    Model( const Model& );
    Model& operator=( const Model& );
  
    std::string                       model_name_;       ///< name of simulation model
    PropertyDatabase<dim>             database_;         ///< where variable specifications are stored
    FiniteElementManager              fem_manager_;      ///< current FiniteElement objects in model
    MeshManager<dim>                  mesh_manager_;     ///< stores mesh: all Node, Element, Face, InterFace objects
    FiniteVolumeStencilManager<dim>*  fvStencilManager_; ///< current finite volume specifications
    bool                              verbose_;          ///< for detailed screen output todo: replace with global verbose singleton
};

/// returns the extent of the model in the x,y,z dimensions and reports this back as a string
std::string  boundingBox( const Model<3U>& sg, double64& dim_x, double64& dim_y, double64& dim_z );

/**

@class Model Model "main_library/Model.h"

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



/**

Database() comes in a constant and a volatile version,
giving you access to the PropertyDatabase object inside the Model.

@return Either a constant or a volatile reference to the PropertyDatabase object.
 */
template<size_t dim>
inline PropertyDatabase<dim>&  Model<dim>::Database() { return database_; }

template<size_t dim>
inline const PropertyDatabase<dim>& Model<dim>::Database() const { return database_; }

template<size_t dim>
inline const FiniteElementManager&  Model<dim>::FE_Manager() const
 { return fem_manager_; }

template<size_t dim>
inline FiniteElementManager&  Model<dim>::FE_Manager()
 { return fem_manager_; }


template<size_t dim>
inline const FiniteVolumeStencilManager<dim>* Model<dim>::FV_Manager() const
{ return fvStencilManager_; }

template<size_t dim>
inline FiniteVolumeStencilManager<dim>*  Model<dim>::FV_Manager()
{ return fvStencilManager_; }


/**

Mesh() comes in a constant and in a volatile version, giving you
access to the MeshManager object inside of the Model object.

@return Either a constant or a volatile reference to the MeshManager.
*/
template<size_t dim>
inline MeshManager<dim>&  Model<dim>::Mesh() { return mesh_manager_; }

template<size_t dim>
inline const MeshManager<dim>&  Model<dim>::Mesh() const { return mesh_manager_; }






} // end namespace csmp

#endif


