#include "Model.h"
#include "VSet.h"
#include "VSetConverter.h"
#include "ModelTopology.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "Interrelation.h"
#include "Visitor.h"
#include "PDE_Integrator.h"
#include "PropertyConstraints.h"
#include "FEM_Data.h"
#include "Box.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"
#include "variableOperations.h"
#include "MeshManagementUtilities.h"
#include "binaryReadWrite.h"
#include "ModelTime.h"
#include "FiniteVolumeStencilManager.h"
#include "UnionFind.h"
#include "plf_colony.h"
#include "VTK_Interface.h"

#include <chrono>
#include <filesystem>

// #define CSMP_MODEL_DEBUG

using namespace std;

namespace csmp {

// PROTECTED CONSTRUCTORS

/**
    Constructs a completely empty model.
    It is used by one of the ANSYS Model interfaces.
 
    @note model will contain a single region called 'model' and it will be considered unique.
    
    @todo SKM deprecate; make the default constructor private to comply with inheritance principles
*/
template<uint32_t dim>
Model<dim>::Model()
  : model_name_("to be named"),
    database_(),
    verbose_(true)
  {
  }


/**
Constructor used for base class construction in subclasses of the model.
Example: ANSYS_Model3D where model construction input data are read from
text file, initializing VSet and ModelTopology classes that serve as an input
for model construction.

@note this is not a stand-alone constructor and leaves the model in an incomplete state.
*/
template<uint32_t dim>
Model<dim>::Model( const char* varTextFile )
  : model_name_( "to be named" ),
    database_( varTextFile ),
    verbose_( true )
{
  InitializeLocalVariableStorage();
}


// PUBLIC CONSTRUCTORS



/**
    Re-constructor: builds model from binary file set, using a binary variables file that has the same name as the Model.
    
    @author SKM
    @date 11/9/2019
    
*/
template<uint32_t dim>
Model<dim>::Model( const string& binaryFileName )
  : model_name_( binaryFileName ),
    database_( BinaryVariablesFileName( binaryFileName.c_str() ).c_str(), set<string>() )
{
  InitializeLocalVariableStorage();
  set<string> empty_set; // prompts model to read all the variables contained in the binary
  InputFromBinaryFile( binaryFileName.c_str(), empty_set );
}



/**
     Reconstructor:  reads model from CSMP's native binary files, but ignoring the binary variable file that comes with the file set.
     Only the variable whose values are actually stored in the VSet are created.
     Extra (additional) storage is created for the variables in the supplied variable file if these are not already contained in
     the VSet.
     
     @param binaryFileName prefix (model name) used in all the binary files relating to the model.
     @param variable_txt_file '-variables.txt' with the specifications of the variables that are desired in addition
     to those that are contained in the VSet.
     
     @author SKM
     @date 26/9/2022
 */
template<uint32_t dim>
Model<dim>::Model( const string& binaryFileName, const string& variable_txt_file )
  : model_name_( binaryFileName ),
    database_( BinaryVariablesFileName( variable_txt_file.c_str() ).c_str(), set<string>() )
{
  InitializeLocalVariableStorage();
  set<string> empty_set; // prompts model to read all the variables contained in the binary
  InputFromBinaryFile( binaryFileName.c_str(), empty_set );
}


/**
    Re-constructor: builds model from binary file set, using only the variables specified in the subset, but with the definitions
    that these variables have in the binary.
    
    @param binaryFileName  '_variables.dat' file with the same name as the model that is stored in the binary file set
    @param subset_variables subset of variables that must be defined by name in the  '_variables.dat' file supplied with the model

@author Junchul Kim
@author SKM refactored 2022
@date 2019

*/
template<uint32_t dim>
Model<dim>::Model( const string& binaryFileName, const set<string>& subset_variables )
  : model_name_( binaryFileName ),
    database_( BinaryVariablesFileName( binaryFileName.c_str() ).c_str(), subset_variables )
{
  InitializeLocalVariableStorage();
  InputFromBinaryFile( binaryFileName.c_str(), subset_variables );
}




/**
Builds 2 and 3D, finite-element meshes according to the
element specifications given in the supplied VSet. The vset will contain
elements of different types, for instance, tetrahedra and prism elements. If,
these are not already stored in the VSet, different material properties
can be assigned to regions inside the the calculation.

@param vset A VSet is supplied as constructor argument and should contain a
finite-element mesh and associated properties with names which must correspond
to the properties specified in the Property input file
(default: CSMP_variables.txt).  Note that VSet is not const because it may be shrunk in construction process).

@section application Application

Constructor is used when an ANSYS Model is built from topology and VData.

@attention No boundaries or regions can be created here because there is no way to store corresponding
 information in the VSet. Use a constructor with Model topology to achieve this.
 
*/
template<uint32_t dim>
Model<dim>::Model( VSet<dim>& vset, const char* var_file )
  : model_name_( "to be named" ),
    database_( var_file ),
    verbose_( true )
{
  Initialize( vset );

} // end VSet constructor



/**
   Build model without variable storage.

   @attention No boundaries or regions are created either because there is no way to store corresponding
     information in the VSet. Use a constructor with Model topology to achieve this.

*/
template<uint32_t dim>
Model<dim>::Model( VSet<dim>& vset )
  : model_name_( "to be named" ),
    database_(),
    verbose_(true)
{
  Initialize( vset );

} // end VSet constructor




/**
Constructor is similar to previous one, but only a subset of the model
geometry stored in the vset will be used according to the specifications
given by the data of the model topology object.

@param vset A VSet is supplied as constructor argument and should contain a
finite-element mesh and associated properties with names which must correspond
to the properties specified in the Property input file
(default: CSMP_variables.txt).

@param mesh_topology a ModelTopology object which should be initialized
as is done for instance by methods of the ANSYS_Interface. The
model topology object has a rich interface which allows the user to
select sub portions of the model stored in the vset.

@section implementation Implementation

The constructor first initializes the MeshManager and then the
MemoryManager with the data from the VSet. Then property data
stored in the VSet as FEM_Data are used to initialize some
of the properties for which memory was allocated for by the
MemoryManager.

@section application Application

Constructor will build a Model for 2D meshes from 'Triangle' as well
as 3D meshes from GoCad or NASTRAN input data.

@section messages Messages

If no properties are supplied, the constructor will issue a Warning message.
In this case, property values have to be assigned before a computation is
carried out since all basic variables are initialized to NAN (not a number)
by default.

@attention per default this constructor will not create any boundaries

@attention the regions file is used only if there is one

*/
template<uint32_t dim>
Model<dim>::Model( ModelTopology& mesh_topology, VSet<dim>& vset,
                   const char* var_file,
                   bool treat_domains_as_regions_and_use_regions_file )
  : model_name_( mesh_topology.ModelName() ),
    database_( var_file ),
    verbose_( true )
{
   if ( treat_domains_as_regions_and_use_regions_file ) {
        const string regions_file_prefix(mesh_topology.ModelName());
        // selects domains via regions file and converts lower-dimensional regions into Boundary objects
        // (use if there are no Face or InterFace objects stored in VSet)
        Initialize( regions_file_prefix.c_str(), mesh_topology, vset );
     }
   // identifies Region, Boundary, and SplitBoundary objects by their names, expecting that corresponding
   // Element, Face and or Interface objects exist in VSet
   else Initialize( mesh_topology, vset );

} // end VSet/ModelTopology constructor








/**  
Initialise( topology, vset, create boundary ...)

Performs the following steps:

0. Using the regions file, selects wanted mesh regions topology and corresponding elements from vset

1. Reduces mesh size and element types to those left over in the topology object

2. initializing the finite-element manager

3. builds finite element mesh and property storage -> done by MeshManager

4. forms unique Region called Model (is in unique regions if there are no other unique regions) else in non-unique Regions

5. Tests with a flood-fill whether the model is contiguous

6. Assigns node and element material properties and variable values to model -> done by InputVariablesFrom(vset)

7. Associates supplied subregions with regions (model subdomains) -> done by FormRegionsFrom(topology)

8. Forms Boundaries -> done by EstablishBoundaries() or EstablishBoxBoundaries() if the model is box-shaped

9. Adds required property storage for regions and boundaries (however their properties are not initialised here)

@attention MOST COMMONLY USED MODEL CONSTRUCTION METHOD FOR  INPUT FROM ANSYS via the CSP interface
 
@attention a fully valid VSet is expected by this method.

@note should only be used for models created externally.
*/
template<uint32_t dim>
void Model<dim>::Initialize( const char* regions_file_prefix, ///< normally this would be equivalent to the model name
                             ModelTopology& mesh_topology,    ///< stores tbe model name as well as the regions in the form of element idx
                             VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( vset.Faces() > 0 || vset.Interfaces() > 0 )
       csmp_error.Note( FATAL_ERROR, "Model<dim>::Initialize(regionfile,ModelTopology,VSet):",
                         "method works only for vsets without Face or InterFace objects as it creates them by itself from lower-dimensiona regions.");
      
     if ( !vset.WithNeighbourConnectivity() )
       csmp_error.Note( FATAL_ERROR, "Model<dim>::Initialize(regionfile,ModelTopology,VSet):", "'pfverts' array is missing.");

    // 1. eliminating unwanted mesh regions from topology and vset, rebuilding boundary flags, check   element numbering etc.
    string prefix( regions_file_prefix );
    if ( filesystem::exists( prefix + "-regions.txt" ) )
      {
        mesh_topology.ReduceToDomains( regions_file_prefix );
        //    reducing the element data to the desired elements specified in the topology object
        //    if the element numbers in the two are different.
        if ( mesh_topology.Cells() != vset.Elements() + vset.Faces() + vset.Interfaces() ) {
            map<size_t,size_t>  old_and_new_elmtids;
            mesh_topology.CreateNewCellNumbers( old_and_new_elmtids );
            vset.ReduceTo( old_and_new_elmtids );
            old_and_new_elmtids.clear();
          }
      }
      
     // the VSet must be correct calling initialise
     if ( (vset.PfvertsBegin() == vset.PfvertsEnd()) )
         csmp_error.Note( FATAL_ERROR, "Model<dim>::Initialize(regionfile,ModelTopology,VSet):", "'pfverts' array is missing.");

    // 2. building the finite element mesh and property storage
    const bool with_FV_variables = (finiteVolumeVariables(Database()) > 0 ) ? true : false;
    mesh_manager_.Initialize( Database(), vset, with_FV_variables );
    const bool contiguous_model( mesh_manager_.IsContiguous() );
    if ( !contiguous_model )
      csmp_error.Note( INFO, "Model<dim>::Initialize(regionfile,ModelTopology,VSet):",
                         "model contains disconnected mesh patches - will attempt to connect them with SplitBoundary objects." );

    // 3. assigning properties to mesh; this does not depend on regions, but region formation may depend on variable values
    InputVariablesFrom( vset );

    // 4. forming default computational domain called "Model" and regions
    const bool place_into_unique_regions{ mesh_topology.ModelDomains() == 0 };
    const size_t elmts = this->FormModelRegion( place_into_unique_regions );
    if ( elmts == 0U )
      csmp_error.Note( FATAL_ERROR, "Model<dim>::Initialize(regionfile,ModelTopology,VSet):", "Region 'Model' has zero elements.");

    // 5. associating supplied subregions with regions (model subdomains)
    const bool ignore_domain_identification_by_name{true};
    this->FormRegionsFrom( mesh_topology, ignore_domain_identification_by_name );
	
    // 6. forming Boundaries
    //    if the model is box-shaped (albeit perhaps with an irregular top surface)
    if (  mesh_topology.BoxShapedModel() ) {
         this->EstablishBoxBoundaries();
         // (re)creating the box-boundary flags (needs respective Boundary objects: see Box.h")
         cout << "\nModel<dim>::Initialize: Since this is a box-shaped model, also, corresponding AT_BOUNDARY flags were created...\n";
      }
    // irregularly shaped models
    else {
          if ( contiguous_model )
            this->EstablishBoundariesFromRegions();
          else
            csmp_error.Note( ERROR, "Model::Intialise(regionfile,ModelTopology,VSet)", "discontiguous model not handled yet");

         // reporting
         this->RegionsOut();
         // this->BoundariesOut(); - was already reported when these were created
      }

    // 7. forming SplitBoundaries if a discontiguous model was detected
    if ( !contiguous_model ) {
         this->DetectAndCreateSplitBoundaries();
         // reporting which boundaries were created
         this->SplitBoundariesOut();
      }

    // 8. adding property storage to the Model
    InitializeLocalVariableStorage();  // for the model
    UpdateSubdomainPropertyStorage();  // for its regions, boundaries and splitboundaries

#ifdef CSMP_MODEL_DEBUG
integrityCheck<dim,Element>( mesh_manager_.CellsBegin(), mesh_manager_.CellsEnd() );
if ( mesh_manager_.Faces() > 0 )
  integrityCheck<dim,Face>( mesh_manager_.FacesBegin(), mesh_manager_.FacesEnd() );
if ( mesh_manager_.InterFaces() > 0 )
  integrityCheck<dim,InterFace>( mesh_manager_.InterFacesBegin(), mesh_manager_.InterFacesEnd() );
#endif

    cout << "\n==================================================================================================";
    cout << "\nModel '"<< this->Name() <<"' has been established successfully ";
    if ( this->Mesh().Elements() > 0 ) {
         size_t volume_elmts{0U}, surface_elmts{0U}, line_elmts{0U};
         cout <<"(total cells "<< currentCellTypes( this->Mesh(), ELEMENT, volume_elmts, surface_elmts, line_elmts );
         cout <<", nodes "<< this->Mesh().Nodes() <<")";
         cout <<"\n\t\t\t("<< Mesh().Elements() <<" elements: volumes "<< volume_elmts <<", surfaces "<< surface_elmts <<", lines "<< line_elmts <<")";
      }
    if ( this->Mesh().Faces() > 0 ) {
         size_t volume_faces{0U}, surface_faces{0U}, line_faces{0U};
         currentCellTypes( this->Mesh(), FACE, volume_faces, surface_faces, line_faces );
         assert( volume_faces == 0U );
         cout <<"\n\t\t\t("<< Mesh().Faces() <<" faces: surfaces "<< surface_faces <<", lines "<< line_faces <<")";
      }
    if ( this->Mesh().InterFaces() > 0 ) {
         size_t volume_ifaces{0U}, surface_ifaces{0U}, line_ifaces{0U};
         currentCellTypes( this->Mesh(), INTER_FACE, volume_ifaces, surface_ifaces, line_ifaces );
         assert( volume_ifaces == 0U );
         cout <<"\n\t\t\t("<< Mesh().InterFaces() <<" interfaces: surfaces "<< surface_ifaces <<", lines "<< line_ifaces <<")";
      }
    cout << "\n==================================================================================================";
    cout << endl;

} // end Initialize (regionfile,ModelTopology,VSet)





/**
    NEW (2022)! - builds model assuming that all information about regions, boundaries or split boundaries is stored in from ModelTopology.
    
    @note If the Model topology object enlists a Boundary object, the corresponding stored idx values are interpreted as Face ids,
    noting that the all entries in the VSet are numbered consecutively and continuously starting with Element followed by Face and InterFace objects.
*/
template<uint32_t dim>
void Model<dim>::Initialize( ModelTopology& mesh_topology, VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( !vset.WithNeighbourConnectivity() )
      csmp_error.Note( FATAL_ERROR, "Model<dim>::Initialize(ModelTopology,VSet):", "'pfverts' array is missing.");

     // the VSet must be correct calling initialise
     if ( (vset.PfvertsBegin() == vset.PfvertsEnd()) )
         csmp_error.Note( FATAL_ERROR, "Model<dim>::Initialize(ModelTopology,VSet):", "'pfverts' array is missing.");

    // 1. building the finite element mesh and property storage
    const bool with_FV_variables = (finiteVolumeVariables(Database()) > 0 ) ? true : false;
    mesh_manager_.Initialize( Database(), vset, with_FV_variables );
    const bool contiguous_model( mesh_manager_.IsContiguous() );
    if ( !contiguous_model )
      csmp_error.Note( INFO, "Model<dim>::Initialize(ModelTopology,VSet):",
                        "model contains disconnected mesh patches. They will be connected with SplitBoundary objects." );

    // 2. assigning properties to mesh; this does not depend on regions, but region formation may depend on variable values
    InputVariablesFrom( vset );

    // 3. forming default computational domain called "Model" and regions
    const bool place_into_unique_regions{ mesh_topology.ModelDomains() == 0 };
    const size_t elmts = this->FormModelRegion( place_into_unique_regions );
    if ( elmts == 0U )
      csmp_error.Note( FATAL_ERROR, "Model<dim>::Initialize(ModelTopology,VSet):", "Region 'Model' has zero elements.");

    // 4. associating supplied subregions with regions (model subdomains)
    this->FormRegionsFrom( mesh_topology );
    this->FormBoundariesFrom( mesh_topology );
    this->FormSplitBoundariesFrom( mesh_topology );
	
    this->RegionsOut();
    this->BoundariesOut();
    this->SplitBoundariesOut();

    // 5. adding property storage to the Model
    InitializeLocalVariableStorage();  // for the model
    UpdateSubdomainPropertyStorage();  // for its regions, boundaries and splitboundaries

#ifdef CSMP_MODEL_DEBUG
integrityCheck<dim,Element>( mesh_manager_.CellsBegin(), mesh_manager_.CellsEnd() );
if ( mesh_manager_.Faces() > 0 )
  integrityCheck<dim,Face>( mesh_manager_.FacesBegin(), mesh_manager_.FacesEnd() );
if ( mesh_manager_.InterFaces() > 0 )
  integrityCheck<dim,InterFace>( mesh_manager_.InterFacesBegin(), mesh_manager_.InterFacesEnd() );
#endif

    cout << "\n==================================================================================================";
    cout << "\nModel '"<< this->Name() <<"' has been established successfully ";
    if ( this->Mesh().Elements() > 0 ) {
         size_t volume_elmts{0U}, surface_elmts{0U}, line_elmts{0U};
         cout <<"(total cells "<< currentCellTypes( this->Mesh(), ELEMENT, volume_elmts, surface_elmts, line_elmts );
         cout <<", nodes "<< this->Mesh().Nodes() <<")";
         cout <<"\n\t\t\t("<< Mesh().Elements() <<" elements: volumes "<< volume_elmts <<", surfaces "<< surface_elmts <<", lines "<< line_elmts <<")";
      }
    if ( this->Mesh().Faces() > 0 ) {
         size_t volume_faces{0U}, surface_faces{0U}, line_faces{0U};
         currentCellTypes( this->Mesh(), FACE, volume_faces, surface_faces, line_faces );
         assert( volume_faces == 0U );
         cout <<"\n\t\t\t("<< Mesh().Faces() <<" faces: surfaces "<< surface_faces <<", lines "<< line_faces <<")";
      }
    if ( this->Mesh().InterFaces() > 0 ) {
         size_t volume_ifaces{0U}, surface_ifaces{0U}, line_ifaces{0U};
         currentCellTypes( this->Mesh(), INTER_FACE, volume_ifaces, surface_ifaces, line_ifaces );
         assert( volume_ifaces == 0U );
         cout <<"\n\t\t\t("<< Mesh().InterFaces() <<" interfaces: surfaces "<< surface_ifaces <<", lines "<< line_ifaces <<")";
      }
    cout << "\n==================================================================================================";
    cout << endl;

  
} // end Initialize (VSet / ModelTopology)


// Testing the Regions that will become boundaries OK
//this->Region("BACK").NodeAttributesToCSV();
//this->Region("RIGHT").NodeAttributesToCSV();
//this->Region("TOP").NodeAttributesToCSV();
//this->Region("LEFT").NodeAttributesToCSV();
//this->Region("BOTTOM").NodeAttributesToCSV();
//this->Region("FRONT").NodeAttributesToCSV();

// DEBUGGING
// -----------------------------------------------------------------------------------------
//VTK_Interface<dim>  vtk_output;
//if ( this->Database().IsDefined("permeability") )
//  vtk_output.OutputDataToVTK( *this, this->Name(), string("permeability"), 1, true );
//else
//  vtk_output.OutputDataToVTK( *this, this->Name(), string("element variable 1"), 1, true );
// -----------------------------------------------------------------------------------------






/**
    Initialises model from VSet. Very similar to Initialise(VSet,ModelTopology), but without
    the creation of regions other than 'Model'.

 @attention a fully (boundary) flagged valid VSet is expected by this method.

*/
template<uint32_t dim>
void Model<dim>::Initialize( VSet<dim>& vset )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( vset.Faces() > 0 ||
       vset.Interfaces() > 0 )
    csmp_error.Note( FATAL_ERROR, "Model<dim>::Initialize(VSet):", "Face and InterFace objects not handled by this method.");

  // 1. checking for neighbor connectivity and boundary flags
  // the VSet must be correct calling initialise
  if ( (vset.PfvertsBegin() == vset.PfvertsEnd()) ) {
       csmp_error.Note( WARNING, "Model<dim>::Initialize(VSet):", "'pfverts' array is missing; establishing it now.");
       if constexpr ( dim == 2 )
         vset.EstablishElementConnectivity2D();
       if constexpr ( dim == 3 )
         vset.EstablishElementConnectivity3D();
    }
  // 2. checking whether BOX_BOUNDARY flags are there which are essential for a model without boundary domains
  if ( vset.BFlags() <= 2 ) {
       // computing a tolerance for the identification of BOX boundaries from the model coordinates
       double tolerance = fabs( vset.X_Range().second );
       tolerance = max( tolerance, fabs( vset.Y_Range().second ) );
       tolerance = max( tolerance, fabs( vset.Z_Range().second ) );
       tolerance *= 1.0e-5;
       csmp_error.Note( WARNING, "Model<dim>::Initialize(VSet):", "'BOX_BOUNDARY' flags are incomplete.");
       VSetConverter<dim>().EstablishBoundaryFlagsForBoxModel( vset, tolerance );
    }

  // 3. building the finite element mesh and property storage
  const bool with_FV_variables = (finiteVolumeVariables(Database()) > 0 ) ? true : false;
  mesh_manager_.Initialize( Database(), vset, with_FV_variables );
  
  // 4. forming default computational domain called "Model" or contiguous multiple domains called "Model_#n"
  const bool unique(true);
  const size_t elmts = this->FormModelRegion( unique );
  if ( elmts == 0U )
    csmp_error.Note( FATAL_ERROR, "Model<dim>::Initialize(VSet):", "Region 'Model' has zero elements.");
  
  cout << "\nModel<dim>::Initialize (VSet): mesh has been built successfully..." << endl;

  // 5. assigning properties to mesh
  InputVariablesFrom( vset );

  cout << "\nModel<dim>::Initialize(VSet): ";
  cout << "Mesh has been built successfully..." << endl;

  // 6. Adding potentially required property storage
  InitializeLocalVariableStorage();
  UpdateSubdomainPropertyStorage();

#ifdef CSMP_MODEL_DEBUG
integrityCheck<dim,Element>( mesh_manager_.CellsBegin(), mesh_manager_.CellsEnd() );
if ( mesh_manager_.Faces() > 0 )
  integrityCheck<dim,Face>( mesh_manager_.FacesBegin(), mesh_manager_.FacesEnd() );
if ( mesh_manager_.InterFaces() > 0 )
  integrityCheck<dim,InterFace>( mesh_manager_.InterFacesBegin(), mesh_manager_.InterFacesEnd() );
#endif

  cout << "\n===========================================================";
  cout << "\nModel: '"<< Name() <<"' has been constructed successfully!";
  cout << "\n===========================================================";
  cout << endl;
  
} // end Initialize(VSet)




/**
     Instantiates FiniteVolumeStencilManager and initialises FiniteVolumeStencils
     if correpoding pointers are NULL.
     
     delegates this step to MeshManager::InitializeFiniteVolumeStencils()
 
*/
template<uint32_t dim>
void Model<dim>::InstantiateFiniteVolumes()
{
   const bool assign_stencils_to_elements{true};
   Mesh().InitializeFiniteVolumeStencils( Database(), assign_stencils_to_elements );
}





template<uint32_t dim>
const char* Model<dim>::Name() const
{
  return model_name_.c_str();
}


template<uint32_t dim>
void Model<dim>::Name( const char* new_name )
{
  model_name_ = new_name;
}





/// removes dynamically allocated finite volume stencil manager
template<uint32_t dim>
Model<dim>::~Model()
{
}



template<uint32_t dim>
void Model<dim>::Verbose( bool verbose )
{
  verbose_ = verbose;
}



template<uint32_t dim>
bool Model<dim>::Verbose() const
{
  return verbose_;
}



template<uint32_t dim>
string Model<dim>::BinaryVsetFileName( const char* base_file_name )
{
  string fullFileName( base_file_name );
  fullFileName.append( ".vset" );
  return fullFileName;
}

template<uint32_t dim>
string Model<dim>::BinaryRegionsFileName( const char* base_file_name )
{
  string fullFileName( base_file_name );
  fullFileName.append( "_regions.dat" );
  return fullFileName;
}

template<uint32_t dim>
string Model<dim>::BinaryBoundariesFileName( const char* base_file_name )
{
  string fullFileName( base_file_name );
  fullFileName.append( "_boundaries.dat" );
  return fullFileName;
}

template<uint32_t dim>
string Model<dim>::BinarySplitBoundariesFileName( const char* base_file_name )
{
  string fullFileName( base_file_name );
  fullFileName.append( "_splitboundaries.dat" );
  return fullFileName;
}

template<uint32_t dim>
string Model<dim>::BinaryVariablesFileName( const char* base_file_name )
{
  string fullFileName( base_file_name );
  fullFileName.append( "_variables.dat" );
  return fullFileName;
}




/// attempts to return the spatial dimension of the model stored in the file (1-3D)
uint32_t spatialDimensionOfModel( const char* csmp_binary )
 {
    throw csmp::Exception( ERROR, "spatialDimensionOfModel(binary file)", "method not implemented yet.");
    return 3U;
    
 } // end dimensionModelInBinaryFile



  /// computes deques of numbered Node, Element, Face and InterFace objects, and outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
template<uint32_t dim>
void  Model<dim>::OutputMeshTo( VSet<dim>& vset, bool get_indices_from_stored_variables )
 {
    if ( get_indices_from_stored_variables ) this->IndexByPropertyValues();
    mesh_manager_.OutputMeshTo( vset,  get_indices_from_stored_variables );
 }




/**
OutputVariableTo() is an overloaded method which will output any kind
of property to the supplied FEM_Data object.

@param out_var the name of the output variable.

The second method argument is a reference to a FEM_Data object which
will store the output property values. All previous property values in
this object will be deleted.

@section application Application

The method is used, for instance, to add property information to a
VSet which is output from the Model. The output FEM_Data
object can however also be used for variable visualization or transfer
to another Model object.

@section messages Messages

When OutputVariableTo() is called one must ascertain that the property
type matches the supplied FEM_Data template, else a fatal error will
be reported.
*/
template<uint32_t dim>
template<class Var>
void  Model<dim>::OutputVariableTo( const char* out_var, FEM_Data<Var>& data ) const
{
  this->Region( "Model" ).OutputVariableTo( out_var, data );

} // end OutputVariableTo

template void Model<1U>::OutputVariableTo( const char*, FEM_Data<ScalarVariable>& ) const;
template void Model<1U>::OutputVariableTo( const char*, FEM_Data<VectorVariable<1U> >& ) const;
template void Model<1U>::OutputVariableTo( const char*, FEM_Data<TensorVariable<1U> >& ) const;
template void Model<1U>::OutputVariableTo( const char*, FEM_Data<ArrayVariable>& ) const;
template void Model<1U>::OutputVariableTo( const char*, FEM_Data<FlaggedArrayVariable>& ) const;

template void Model<2U>::OutputVariableTo( const char*, FEM_Data<ScalarVariable>& ) const;
template void Model<2U>::OutputVariableTo( const char*, FEM_Data<VectorVariable<2U> >& ) const;
template void Model<2U>::OutputVariableTo( const char*, FEM_Data<TensorVariable<2U> >& ) const;
template void Model<2U>::OutputVariableTo( const char*, FEM_Data<ArrayVariable>& ) const;
template void Model<2U>::OutputVariableTo( const char*, FEM_Data<FlaggedArrayVariable>& ) const;

template void Model<3U>::OutputVariableTo( const char*, FEM_Data<ScalarVariable>& ) const;
template void Model<3U>::OutputVariableTo( const char*, FEM_Data<VectorVariable<3U> >& ) const;
template void Model<3U>::OutputVariableTo( const char*, FEM_Data<TensorVariable<3U> >& ) const;
template void Model<3U>::OutputVariableTo( const char*, FEM_Data<ArrayVariable>& ) const;
template void Model<3U>::OutputVariableTo( const char*, FEM_Data<FlaggedArrayVariable>& ) const;







/**
Reads the values of scalar, vector, tensor, array and flagged array variables stored
in the VSet and assigns them to the current model.

All the placements node, element, element integration point, finite volue facet and
sector integration point as well as region are read from file.
However, the region placement works only if there is only a single region in the
model.
*/
template<uint32_t dim>
void Model<dim>::InputVariablesFrom( const VSet<dim>& vset )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !vset.DataEmpty() ) mesh_manager_.InputStoredVariablesFrom( Database(), vset );
  else ErrorHandler::Instance().Note( INFO, "Model<dim>::InputVariablesFrom:", "No properties found in VSet." );

  // properties and values stored on the model itself
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    // apart from the name string key in the map, PropertyData contains the most important variable specifications
    if ( (*pit).second.Placement() != MODEL ) continue;
    // some checks
    assert( Database().IsDefined( (*pit).first.c_str() ) );
    const csmp::Index key( Database().StorageKey( (*pit).first.c_str() ) );
    assert( key.place == MODEL );
    assert( (*pit).second.Size() / key.dataDepth == 1U );

    switch ( key.type )
      {
        case SCALAR: {
              ScalarVariable value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        case VECTOR: {
              VectorVariable<dim> value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        case TENSOR: {
              TensorVariable<dim> value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        case ARRAY: {
              ArrayVariable value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        case FLAGGEDARRAY: {
              FlaggedArrayVariable value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        default:
          csmp_error.Note( ERROR, "Model<dim>::InputVariablesFrom:",
                             (*pit).first, "type of Model variable not recognized." );
      }
  }
  
} // end InputVariablesFrom




/**
InputVariableFrom() lets you input variable data stored in a FEM_Data
template class object to a Model variable. class T here is a place
holder for the data type which may be a double or any basic CSP variable.

@param input_prop the name of the variable to which the data
are to be assigned.
@param vdata is a reference to the FEM_Data data object.

@section implementation Implementation

How much do I await the day when member templates defined outside of class
declarations will be mastered by the gang of C++ compilers. Since this is
not currently the case you have to bear with overloaded method
cluttered class interfaces.

@section application Application

FEM_Data class objects are convinient containers for CSP property data
and are used also by peripheral programs such as the OpenGL viewing tool.
They can also be used to compare values of different properties of the
same type and placement using STL algorithms.

@section messages Messages

InputVariableFrom() will report an error and return without completing its
task, if the variable type or placement in the Property database does not
match the specifications of the FEM_Data dataset.

@attention all elements and nodes need to be numbered uniquely and in the same
fashion as they were output to file for this method to work correctly.
*/
template<uint32_t dim>
template<class Var>
void  Model<dim>::InputVariableFrom( const char* input_prop, const FEM_Data<Var>& vdata )
{
  const string master_region( "Model" );
  this->Region( master_region.c_str() ).InputVariableFrom( input_prop, vdata );

} // end InputVariableFrom

template void Model<1U>::InputVariableFrom( const char*, const FEM_Data<ScalarVariable>& );
template void Model<1U>::InputVariableFrom( const char*, const FEM_Data<VectorVariable<1U> >& );
template void Model<1U>::InputVariableFrom( const char*, const FEM_Data<TensorVariable<1U> >& );
template void Model<1U>::InputVariableFrom( const char*, const FEM_Data<ArrayVariable>& );
template void Model<1U>::InputVariableFrom( const char*, const FEM_Data<FlaggedArrayVariable>& );

template void Model<2U>::InputVariableFrom( const char*, const FEM_Data<ScalarVariable>& );
template void Model<2U>::InputVariableFrom( const char*, const FEM_Data<VectorVariable<2U> >& );
template void Model<2U>::InputVariableFrom( const char*, const FEM_Data<TensorVariable<2U> >& );
template void Model<2U>::InputVariableFrom( const char*, const FEM_Data<ArrayVariable>& );
template void Model<2U>::InputVariableFrom( const char*, const FEM_Data<FlaggedArrayVariable>& );

template void Model<3U>::InputVariableFrom( const char*, const FEM_Data<ScalarVariable>& );
template void Model<3U>::InputVariableFrom( const char*, const FEM_Data<VectorVariable<3U> >& );
template void Model<3U>::InputVariableFrom( const char*, const FEM_Data<TensorVariable<3U> >& );
template void Model<3U>::InputVariableFrom( const char*, const FEM_Data<ArrayVariable>& );
template void Model<3U>::InputVariableFrom( const char*, const FEM_Data<FlaggedArrayVariable>& );






/**
     Sets the Idx values of the nodes, elements, faces, and interfaces of the model according to the stored values, if any.
     @code
     "node number", "element number", "face number", "interface number"
     @endcode
     
     @attention if elements were created or deleted, then there may now be gaps in the numbering.
*/
template<uint32_t dim>
void Model<dim>::IndexByPropertyValues()
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    MeshManager<dim>& mesh(this->Mesh());
    const bool in_a_single_sequence{true};
    
    if ( Database().IsDefined( "node number" ) ) {
         const csmp::Index key = this->Database().StorageKey("node number");
         const size_t n_nodes{mesh.Nodes()};
         const typename plf::colony<csmp::Node<dim>>::iterator nodes_end(mesh.NodesEnd());
         for ( auto nit=mesh.NodesBegin(); nit!=nodes_end; ++nit ) {
              const double number{ (*nit).Read(key) };
              if ( !isnan(number) ) {
                   const size_t node_number = static_cast<size_t>(number);
                   if ( node_number >= n_nodes )
                     csmp_error.Note( WARNING, "Model::IndexByPropertyValues:",
                                               "'node number' exceeds range of available nodes:", to_string(node_number) );
                   (*nit).Idx( node_number );
                }
              else {
                     (*nit).Idx( numeric_limits<size_t>::max() );
                     csmp_error.Note( WARNING, "Model::IndexByPropertyValues:",
                                     "double turned into 'node number' was a NaN; setting value to MAX");
                }
           }
      }
    else {
         csmp_error.Note( ERROR, "Model::IndexByPropertyValues:",
                                  "'node number' is not defined; default unique contiguous numbering will be used" );
         mesh.AssignUniqueNumbers( in_a_single_sequence );
         return;
      }
    
    if ( Database().IsDefined( "element number" ) ) {
         const csmp::Index key = this->Database().StorageKey("element number");
         const size_t n_elmts{mesh.Elements()};
         const typename plf::colony<csmp::Element<dim>>::iterator elmts_end(mesh.ElementsEnd());
         for ( auto it=mesh.ElementsBegin(); it!=elmts_end; ++it ) {
              const double number{ (*it).Read(key) };
              if ( !isnan(number) ) {
                   const size_t elmt_number = static_cast<size_t>(number);
                   if ( elmt_number >= n_elmts )
                     csmp_error.Note( WARNING, "Model::IndexByPropertyValues:",
                                       "'element number' exceeds range of available elements:", to_string(elmt_number) );
                   (*it).Idx( elmt_number );
                }
              else {
                     (*it).Idx( numeric_limits<size_t>::max() );
                     csmp_error.Note( WARNING, "Model::IndexByPropertyValues:",
                                     "double turned into 'element number' was a NaN; setting value to MAX");
                }
           }
      }
    else {
         csmp_error.Note( ERROR, "Model::IndexByPropertyValues:",
                                  "'element number' is not defined; default unique contiguous numbering will be used" );
         mesh.AssignUniqueNumbers( in_a_single_sequence );
         return;
      }
    
    // faces
    if ( Mesh().Faces() > 0 ) {
        if ( Database().IsDefined( "face number" ) ) {
             const csmp::Index key = this->Database().StorageKey("face number");
             const size_t n_faces{mesh.Faces()};
             const typename plf::colony<csmp::Face<dim>>::iterator faces_end(mesh.FacesEnd());
             for ( auto it=mesh.FacesBegin(); it!=faces_end; ++it ) {
                  const double fp_val = (*it).Read(key);
                  if ( isnan( fp_val) )
                    csmp_error.Note( ERROR, "Model::IndexByPropertyValues:",
                                  "'face number' value is NaN, replacing with 'uint32_' max." );
                  const size_t face_number = ( isnan(fp_val) ) ? numeric_limits<uint32_t>::max() : static_cast<uint32_t>(fp_val);
                  if ( face_number >= n_faces + mesh.Elements() )
                    csmp_error.Note( WARNING, "Model::IndexByPropertyValues:",
                                      "'face number' exceeds range of available faces:",
                                             to_string(face_number) );
                  (*it).Idx( face_number );
               }
          }
        else {
             csmp_error.Note( ERROR, "Model::IndexByPropertyValues:",
                                     "'face number' is not defined; default unique contiguous numbering will be used" );
             mesh.AssignUniqueNumbers( in_a_single_sequence );
             return;
          }
      }

    // interfaces
    if ( Mesh().InterFaces() > 0 ) {
        if ( Database().IsDefined( "interface number" ) ) {
             const csmp::Index key = this->Database().StorageKey("interface number");
             const size_t n_ifaces{mesh.InterFaces()};
             const size_t n_all_cells{ n_ifaces + mesh.Faces() + mesh.Elements() };
             const typename plf::colony<csmp::InterFace<dim>>::iterator ifaces_end(mesh.InterFacesEnd());
             for ( auto it=mesh.InterFacesBegin(); it!=ifaces_end; ++it ) {
                  const size_t iface_number = static_cast<uint32_t>((*it).Read(key));
                  if ( iface_number >= n_all_cells )
                    csmp_error.Note( WARNING, "Model::IndexByPropertyValues:",
                                      "'interface number' exceeds range of available faces:",
                                            to_string(iface_number) );
                  (*it).Idx( iface_number );
               }
          }
        else {
             csmp_error.Note( ERROR, "Model::IndexByPropertyValues:",
                                      "'interface number' is not defined; default unique contiguous numbering will be used" );
             mesh.AssignUniqueNumbers( in_a_single_sequence );
          }
      }

 } // end IndexByPropertyValues





/**
Add a property at runtime
@warning Any runtime change in variables invalidates existing Index objects! They need to be refreshed.
@todo (1-F) Support for BOUNDARY... / SPLIT_BOUNDARY (create overload taking index)
*/
template<uint32_t dim>
csmp::Index  Model<dim>::CreateProperty( const char* new_prop,
                                         const char* notation,
                                         const char* unit,
                                         VARIABLE_TYPE vtype,
                                         PLACEMENT vplace,
                                         uint32_t vsize,
                                         double vmin,
                                         double vmax,
                                         string usage )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( database_.IsDefined( new_prop ) ) {
    csmp_error.Note( WARNING, "Model<dim>::CreateProperty:", new_prop, "property already exists." );
    return database_.StorageKey( new_prop );
  }

  const uint32_t  prop_index = database_.VariableCount( vplace, vtype );
  csmp::Index     prop_key   = database_.AddProperty( new_prop, notation, unit, prop_index, vtype, vplace, vsize, vmin, vmax, usage );

  if ( vplace == NODE ) {
    csmp::Region<dim>&  gref( this->Region( "Model" ) );
    for ( auto eit = gref.NodesBegin(); eit != gref.NodesEnd(); eit++ )
      (*eit)->AddProperty( prop_key );
  }
  else if ( vplace == ELEMENT || vplace == ELEMENT_INTEGRATION_POINT || vplace == SECTOR_INTEGRATION_POINT || vplace == FACET_INTEGRATION_POINT ) {
    csmp::Region<dim>&  gref( this->Region( "Model" ) );
    for ( auto eit = gref.CellsBegin(); eit != gref.CellsEnd(); eit++ )
      (*eit)->AddProperty( prop_key );
  }
  else if ( vplace == FACE || vplace == FACE_INTEGRATION_POINT || vplace == FACE_SECTOR_INTEGRATION_POINT || vplace == FACE_FACET_INTEGRATION_POINT ) {
    for ( auto git = this->BoundariesBegin(); git != this->BoundariesEnd(); ++git )
      for ( auto eit = (*git).second.CellsBegin(); eit != (*git).second.CellsEnd(); eit++ )
        (*eit)->AddProperty( prop_key );
  }
  else if ( vplace == INTER_FACE || vplace == INTER_FACE_INTEGRATION_POINT || vplace == INTER_FACE_SECTOR_INTEGRATION_POINT || vplace == INTER_FACE_FACET_INTEGRATION_POINT ) {
    for ( auto git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); ++git )
      for ( auto eit = (*git).second.CellsBegin(); eit != (*git).second.CellsEnd(); eit++ )
        (*eit)->AddProperty( prop_key );
  }
  else if ( vplace == MODEL ) {
     LocalVariableStorage<dim,Model>::AddProperty( prop_key );
  }
  else if ( vplace == REGION ) {
    // each region has its unique property value even if they are overlapping!
    for ( auto git = this->UniqueRegionsBegin(); git != this->UniqueRegionsEnd(); git++ )
      (*git).second.AddProperty( prop_key );
    for ( typename map<string, csmp::Region<dim> >::iterator
          git = this->RegionsBegin(); git != this->RegionsEnd(); git++ )
      (*git).second.AddProperty( prop_key );
  }
  else if ( vplace == BOUNDARY ) {
    for ( typename map<string, csmp::Boundary<dim> >::iterator
          git = this->BoundariesBegin(); git != this->BoundariesEnd(); git++ )
      (*git).second.AddProperty( prop_key );
    for ( typename map<string, csmp::SplitBoundary<dim> >::iterator
          git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); git++ )
      (*git).second.AddProperty( prop_key );
  }
  else
    throw csmp::Exception( ERROR, "Model<dim>::CreateProperty:",
                           new_prop, "placement not supported yet, nothing was done." );

  if ( Verbose() )
    cout << "\nModel<dim>::CreateProperty: successfully created the new property '" << new_prop << "'\n";
  return prop_key;

} // end CreateProperty





/**
Removes a property from the model.
@warning All indices have to be updated.
@todo (1-F) Support for new placements
@warning Existing Index objects may be invalidated after invoking this function
*/
template<uint32_t dim>
void  Model<dim>::DeleteProperty( const char* property )
{

  if ( !database_.IsDefined( property ) ) {
    throw csmp::Exception( INFO,
                           "Model<dim>::CreateProperty",
                           property,
                           "property does not exist" );
    return;
  }

  csmp::Index  prop_key = database_.StorageKey( property );
  csmp::Region<dim>&  gref( this->Region( "Model" ) );

  if ( prop_key.place == NODE ) {
    for ( auto eit = gref.NodesBegin(); eit != gref.NodesEnd(); eit++ )
      (*eit)->DeleteProperty( prop_key );
  }
  else if ( prop_key.place == ELEMENT || prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACET_INTEGRATION_POINT ) {
    for ( auto eit = gref.CellsBegin(); eit != gref.CellsEnd(); eit++ )
      (*eit)->DeleteProperty( prop_key );
  }
  else if ( prop_key.place == FACE || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT ) {
    for ( auto git = this->BoundariesBegin(); git != this->BoundariesEnd(); ++git )
      for ( auto eit = (*git).second.CellsBegin(); eit != (*git).second.CellsEnd(); eit++ )
        (*eit)->DeleteProperty( prop_key );
  }
  else if ( prop_key.place == INTER_FACE || prop_key.place == INTER_FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT || prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
    for ( auto git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); ++git )
      for ( auto eit = (*git).second.CellsBegin(); eit != (*git).second.CellsEnd(); eit++ )
        (*eit)->DeleteProperty( prop_key );
  }
  else if ( prop_key.place == MODEL ) {
     LocalVariableStorage<dim,Model>::DeleteProperty( prop_key );
  }
  else if ( prop_key.place == REGION ) {
    // each region has its unique property value even if they are overlapping!
    for ( auto git = this->UniqueRegionsBegin(); git != this->UniqueRegionsEnd(); git++ )
      (*git).second.DeleteProperty( prop_key );
    for ( auto git = this->RegionsBegin(); git != this->RegionsEnd(); git++ )
      (*git).second.DeleteProperty( prop_key );
  }
  else if ( prop_key.place == BOUNDARY ) {
    for ( auto git = this->BoundariesBegin(); git != this->BoundariesEnd(); git++ )
      (*git).second.DeleteProperty( prop_key );
    for ( auto git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); git++ )
      (*git).second.DeleteProperty( prop_key );
  }
  else
    throw csmp::Exception( ERROR, "Model<dim>::DeleteProperty", "Placement not supported yet" );

  database_.DeleteProperty( property );

} // end DeleteProperty




/**
In a finite-element mesh, each node is shared by a variable number of
elements which also may differ in their area. This method extrapolates
element properties to node properties by calculating the element property
average of all elements which share the node, weighted by element area
for each node of the mesh.

'by_distance' means that weighting is not by the area/volume of the neighbouring
element but by the distance of its barycenter to the node to which the
extrapolation is done.


@param eprop The name of the element property which is extrapolated,
@param nprop the name of the node property to which the extrapolation
is written to.

@section implementation Implementation

The element property values are weighted by the area of the element and then
added for each node to a property array vector, the number of additions to each
node is counted.

@section application Application

To obtain higher order (more continuous) velocity fields, the
TransportVisitor object uses the element property 'velocity' and the
'nodal velocity' at the same time. Other applications include interrelations
where the output property is a node variable, but which require element
properties in the calculations.

While the extrapolation approach used in this method is topologically
sound, the method may still produce somewhat patchy looking results in meshes
where the element size varies greatly and a node lies on the boundary of
2 regions with strongly differing material properties.

@section messages Messages

The method checks the placement of the specified variables and reports an
error if it does not match the specifications from above. In this case
it will return without carrying out the extrapolation.
*/
template<uint32_t dim>
void  Model<dim>::ExtrapolateCellToNodeProperty( const char* eprop, const char* nprop, bool by_distance )
{
  this->Region( "Model" ).ExtrapolateCellToNodeProperty( eprop, nprop, by_distance );

} // end ExtrapolateCellToNodeProperty





/**
The status of the distributed physical variable as identified by the string
char*, is changed from the first VARIABLE_FLAG argument to the second. This
is done only where the distributed physical variable has a status equivalent
to that given by the first VARIABLE_FLAG argument of this method.

@param input_prop Name of property the status of which shall be changed at where it is
flagged equivalent to argument two.

@section application Application

To drop Dirichlet boundary conditions from a previous timestep etc.

@section messages Messages

If no assignments can be made by ChangeConditionFromTo(), because no
object properties with the desired flags are found, an error is reported.
*/
template<uint32_t dim>
void Model<dim>::ChangePropertyStatus( const char* input_prop, VARIABLE_FLAG new_status )
{
  this->Region( "Model" ).ChangePropertyStatus( input_prop, new_status, COMPLETE );

} // end ChangePropertyStatusTo




/**
Where the specified property has a value within the range given by min-max,
its status is changed.

The status of the distributed physical variable as identified by the
string char*, is changed to the status defined by the VARIABLE_FLAG argument
of this method, if the value of the distributed variable lies within the
range given as open interval by the fTs min and max. This method
only works for scalar and vector variables. The status of vector variables
is changed only if all elements of the VectorVariable<dim> are within the
defined value range.

@section application Application

Data in the specified range may be preserved if set to a DIRICH value or
one may choose to form a group on the basis of their status or by
intersecting a Region with such a group.

@section messages Messages

Thus far this method cannot change the status of TensorVariable<dim>'s and
will therefore issue a warning if this is attempted.

Also if no status changes were made, the method will report this.
*/
template<uint32_t dim>
void Model<dim>::ChangePropertyStatusWhere( const char* var,
                                            double vmin, double vmax,
                                            VARIABLE_FLAG new_status )
{
  this->Region( "Model" ).ChangePropertyStatusWhere( var, new_status, vmin, vmax );

} // end ChangePropertyStatusWhere



/**
The User's guide describes the Visitor design pattern which was
adopted from Gamma et al. 1994 (Design Patterns, Addison Wesley, p. 331).
Depending on its internally specified ApplicationLevel, Accept() either applies
the Visitor to the entire Model or passes it down to members of the Model
hierarchy (Region > > Boundary > Element > Face > InterFace > Node.
In the last case, the visitor is passed to all nodes of the
entire Model.

@param v Any visitor derived from the Visitor base class.

@section implementation Implementation

The method consists a of switch( visitor.ApplicationLevel() ) statement with
iterator loops through the different object collections.
These deliver the visitor to its ApplicationTarget().

@section application Application

Standard method to give visitors access to the Model.
@attention visitors will only be applied to the unique regions in the model
in order to avoid duplicate application.
If a non-unique region is the target, the visitor should be passed directly to it.

@section messages Messages

Accept() will report an error if the visitors destination is a specific region, which
is not defined.

@todo SKM (1) Reorganize this so that no confusion remains where the visitor goes.
*/
template<uint32_t dim>
void Model<dim>::Accept( csmp::Visitor<dim>& v )
{
  switch ( v.ApplicationLevel() ) {
    // means direct access to elements and nodes without being able to fetch
    // data from regions or boundaries
    case MODEL: {
      csmp::Region<dim>&  mref( this->Region( "Model" ) );
      v.Visit( this );
      switch ( v.ApplicationTarget() ) {
        case ELEMENT: {
          for ( auto el_it = mref.CellsBegin(); el_it != mref.CellsEnd(); el_it++ )
            (*el_it)->Accept( v );
        }
                      return;
        case NODE: {
          for ( auto nd_it = mref.NodesBegin(); nd_it != mref.NodesEnd(); nd_it++ )
            (*nd_it)->Accept( v );
        }
                   return;
        case FACE: { // all faces contained in the model live inside of the boundaries
          for ( auto bit = this->BoundariesBegin(); bit != this->BoundariesEnd(); bit++ )
            for ( auto it = (*bit).second.CellsBegin(); it != (*bit).second.CellsEnd(); it++ )
              (*it)->Accept( v );
        }
                   return;
        case INTER_FACE:
          for ( auto bit = this->SplitBoundariesBegin(); bit != this->SplitBoundariesEnd(); bit++ )
            for ( auto it = (*bit).second.CellsBegin(); it != (*bit).second.CellsEnd(); it++ )
              (*it)->Accept( v );
          return;
          // access is granted to all boundaries of each type
        case SPLIT_BOUNDARY:
          for ( auto it( this->SplitBoundariesBegin() ); it != this->SplitBoundariesEnd(); ++it )
            (*it).second.Accept( v );
          break;
        case BOUNDARY:
          for ( auto it( this->BoundariesBegin() ); it != this->BoundariesEnd(); ++it )
            (*it).second.Accept( v );
          break;
        default:
          throw csmp::Exception( WARNING, "Model<dim>::Accept:",
                                 "application target of visitor unresolved; nothing was done." );
      }
    }
                return;
    case REGION: { // processing all unique regions except for the "Model"
      typename map<string, csmp::Region<dim> >::iterator  rit = this->UniqueRegionsBegin();
      // if there is only the region "Model"
      if ( this->Regions() == 1U and (*rit).first == "Model" ) {
        (*rit).second.Accept( v );
        return;
      }
      else {
        // visitor is applied to all unique regions except for "Model"
        while ( rit != this->UniqueRegionsEnd() ) {
          if ( (*rit).first != "Model" ) (*rit).second.Accept( v );
          rit++;
        }
      }
    }
                 return;
    case SPLIT_BOUNDARY:
      for ( typename Model<dim>::splitBoundaryIterator it( this->SplitBoundariesBegin() );
      it != this->SplitBoundariesEnd();
        ++it )
        (*it).second.Accept( v );
      break;
    case BOUNDARY:
      if ( v.ApplicationTarget() == FACE || v.ApplicationTarget() == BOUNDARY )
        for ( auto it = this->BoundariesBegin(); it != this->BoundariesEnd(); ++it )
          (*it).second.Accept( v );
      // TODO: check wether this is redundant
      else if ( v.ApplicationTarget() == INTER_FACE || v.ApplicationTarget() == SPLIT_BOUNDARY )
        for ( auto it = this->SplitBoundariesBegin(); it != this->SplitBoundariesEnd(); ++it )
          (*it).second.Accept( v );
      return;
    default:
      throw csmp::Exception( WARNING, "Model<dim>::Accept",
                             "application level of visitor unresolved, nothing done..." );
  }
} // end Accept




/**

Allows assignment of uniform property values  for all CSMP variable types.
Any previous values are overwritten irrespective of their flags and will receive the flag of the input variable.
It follows that  potential boundary conditions need to be reassigned after this method has been called.

Note that the method checks the input variables against the ranges
specified in the variables database (see '*-variables.txt' file). If the values
fall out of these ranges, the method will reset values to the
nearest range bound.

@param input_prop The name of the input property (which is checked against the property
database),
@param value and the initialized CSP basic variable of scalar, vector, or tensor type.

@section messages Messages

If the variable name is not defined or the variable value falls outside
of the range that was specified in the current '-variables.txt' file
the method will report errors.

Similarly when the property placement (i.e. element, node etc.) is
not recognized, an error will be reported.

*/
template<uint32_t dim>
template<class T>
void Model<dim>::InputPropertyValue( const char* input_prop, const T& value )
{
  const csmp::Index prop_key = this->database_.StorageKey( input_prop );
  // properties stored on the model
  if ( prop_key.place == MODEL ) {
       this->Store( prop_key, value );
       return;
    }
    
  // all regions
  if ( prop_key.place == REGION ) {
      for ( auto it = this->UniqueRegionsBegin(); it != this->UniqueRegionsEnd(); ++it )
        (*it).second.Store( prop_key, value );
      for ( auto it = this->RegionsBegin(); it != this->RegionsEnd(); ++it )
        (*it).second.Store( prop_key, value );
      return;
    }
    
  // all boundaries  
  if ( prop_key.place == BOUNDARY ) {
       for ( auto it = this->BoundariesBegin(); it != this->BoundariesEnd(); ++it )
         (*it).second.Store( prop_key, value );
      return;
   }
  
  // split boundaries 
  if ( prop_key.place == SPLIT_BOUNDARY ) {
       for ( auto it = this->SplitBoundariesBegin(); it != this->SplitBoundariesEnd(); ++it )
         (*it).second.Store( prop_key, value );
      return;
   }
   
  // node and element properties are handled by direct access to model domain 
  if ( prop_key.place == ELEMENT || prop_key.place == ELEMENT_INTEGRATION_POINT ||
       prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == SECTOR_INTEGRATION_POINT ||
       prop_key.place == NODE ) {
      this->Region( "Model" ).InputPropertyValue( input_prop, value, COMPLETE );
      return;
   }
  
  // properties discretised only on boundaries
  if ( prop_key.place == FACE || prop_key.place == FACE_INTEGRATION_POINT ||
       prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT )
    {
       for ( auto it = this->BoundariesBegin(); it != this->BoundariesEnd(); ++it )
         (*it).second.InputPropertyValue( input_prop, value, COMPLETE );
       return;
    }
  
  // properties discretised on split boundaries  
  if ( prop_key.place == INTER_FACE || prop_key.place == INTER_FACE_INTEGRATION_POINT ||
       prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
    {
       for ( auto it = this->SplitBoundariesBegin(); it != this->SplitBoundariesEnd(); ++it )
         (*it).second.InputPropertyValue( input_prop, value, COMPLETE );
       return;
    }
  
  // if something fell through the cracks
  throw csmp::Exception( ERROR, "Model<dim>::InputPropertyValue", input_prop,
                        "placement was not recognised. No assignments were made");
  
} // end InputPropertyValue



  // instantiations of extra member function templates

  // 3D
  // scalar variables
template void Model<3U>::InputPropertyValue<ScalarVariable>( const char*, const ScalarVariable& );

// vector variables
template void Model<3U>::InputPropertyValue<VectorVariable<3U> >( const char*, const VectorVariable<3U>& );

// tensor variables
template void Model<3U>::InputPropertyValue<TensorVariable<3U> >( const char*, const TensorVariable<3U>& );

// array variable
template void Model<3U>::InputPropertyValue<ArrayVariable>( const char*, const ArrayVariable& );

// flagged array variable
template void Model<3U>::InputPropertyValue<FlaggedArrayVariable>( const char*, const FlaggedArrayVariable& );

// 2D
// scalar variables
template void Model<2U>::InputPropertyValue<ScalarVariable>( const char*, const ScalarVariable& );

// vector variables
template void Model<2U>::InputPropertyValue<VectorVariable<2U> >( const char*, const VectorVariable<2U>& );

// tensor variables
template void Model<2U>::InputPropertyValue<TensorVariable<2U> >( const char*, const TensorVariable<2U>& );

// array variable
template void Model<2U>::InputPropertyValue<ArrayVariable>( const char*, const ArrayVariable& );

// flagged array variable
template void Model<2U>::InputPropertyValue<FlaggedArrayVariable>( const char*, const FlaggedArrayVariable& );

// 1D
// scalar variables
template void Model<1U>::InputPropertyValue<ScalarVariable>( const char*, const ScalarVariable& );

// vector variables
template void Model<1U>::InputPropertyValue<VectorVariable<1U> >( const char*, const VectorVariable<1U>& );

// tensor variables
template void Model<1U>::InputPropertyValue<TensorVariable<1U> >( const char*, const TensorVariable<1U>& );

// array variable
template void Model<1U>::InputPropertyValue<ArrayVariable>( const char*, const ArrayVariable& );

// flagged array variable
template void Model<1U>::InputPropertyValue<FlaggedArrayVariable>( const char*, const FlaggedArrayVariable& );



/**
Given a box-shaped model, ie. a model whose boundaries are rectangular,
this method assigns uniform variable values to each of these specific boundaries
including edges and model corner points.

Assignments to individual edges or corners are also supported.

@attention elements are considered boundary elements if at least one
of their faces is on this boundary.

@attention each side boundary shares edges and corner points with the adjacent one.
Therefore methods like 'belongsToSide" must be used to capture all the nodes that belong
to an individual boundary.

@test updated by SKM 29/5/14.
*/
template<uint32_t dim>
template<class T>
void Model<dim>::InputBoundaryValue( BOX_BOUNDARY boundary, const char* input_prop, const T& value )
{
  const csmp::Index  prop_key = database_.StorageKey( input_prop );

  if ( prop_key.type == TENSOR )
    throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryValue: input variable:", input_prop,
                           "method has not been implemented for TensorVariables." );

  if ( prop_key.place != NODE )
    throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryValue: input variable:", input_prop,
                           "method can only be applied to node variables." );

  csmp::Region<dim>&  model_domain( this->Region( "Model" ) );

  if ( isSide( boundary ) ) {
    for ( typename vector<Node<dim>*>::const_iterator
          nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ )
      if ( belongsToSide( boundary, (*nit)->AtBoundary() ) ) (*nit)->Store( prop_key, value );
    return;
  }

  if ( isEdge( boundary ) ) {
    for ( typename vector<Node<dim>*>::const_iterator
          nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ )
      if ( belongsToEdge( boundary, (*nit)->AtBoundary() ) ) (*nit)->Store( prop_key, value );
    return;
  }

  if ( isCorner( boundary ) ) {
    for ( typename vector<Node<dim>*>::const_iterator
          nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ )
      if ( (*nit)->AtBoundary() == boundary ) (*nit)->Store( prop_key, value );
    return;
  }

  throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryValue: ",
                         "boundary flag could not be parsed." );

} // end InputBoundaryValue


template void Model<1U>::InputBoundaryValue<ScalarVariable >( BOX_BOUNDARY, const char*, const ScalarVariable& );
template void Model<2U>::InputBoundaryValue<ScalarVariable >( BOX_BOUNDARY, const char*, const ScalarVariable& );
template void Model<3U>::InputBoundaryValue<ScalarVariable >( BOX_BOUNDARY, const char*, const ScalarVariable& );

template void Model<1U>::InputBoundaryValue<VectorVariable<1U> >( BOX_BOUNDARY, const char*, const VectorVariable<1U>& );
template void Model<2U>::InputBoundaryValue<VectorVariable<2U> >( BOX_BOUNDARY, const char*, const VectorVariable<2U>& );
template void Model<3U>::InputBoundaryValue<VectorVariable<3U> >( BOX_BOUNDARY, const char*, const VectorVariable<3U>& );







template<uint32_t dim>
void Model<dim>::InputBoundaryFlags( BOX_BOUNDARY boundary, const char* property, const vector<VARIABLE_FLAG>& flags )
{
  const csmp::Index  prop_key = database_.StorageKey( property );

  if ( flags.empty() )
    throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags: ",
                           "boundary flag vector is empty, nothing was done" );

  if ( prop_key.type == TENSOR )
    throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags: ",
                           "method has not been implemented for TensorVariables" );

  if ( prop_key.place != NODE )
    throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags: ",
                           "method can only be applied to node variables" );

  csmp::Region<dim>&  model_domain( this->Region( "Model" ) );

  if ( prop_key.type == SCALAR ) {
    if ( isSide( boundary ) ) {
      for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ )
        if ( belongsToSide( boundary, (*nit)->AtBoundary() ) ) (*nit)->Status( prop_key, flags[0] );
      return;
    }

    if ( isEdge( boundary ) ) {
      for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ )
        if ( belongsToEdge( boundary, (*nit)->AtBoundary() ) ) (*nit)->Status( prop_key, flags[0] );
      return;
    }

    if ( isCorner( boundary ) ) {
      for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ )
        if ( (*nit)->AtBoundary() == boundary ) (*nit)->Status( prop_key, flags[0] );
      return;
    }
  }
  else { // vector variables
    if ( flags.size() != dim )
      throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags (vector variable): ",
                             "boundary flag vector has the wrong number of entries" );
    if ( isSide( boundary ) ) {
      for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ )
        if ( belongsToSide( boundary, (*nit)->AtBoundary() ) )
          for ( auto i{0U}; i<dim; i++ ) (*nit)->Status( prop_key, i, flags[i] );
      return;
    }

    if ( isEdge( boundary ) ) {
      for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ )
        if ( belongsToEdge( boundary, (*nit)->AtBoundary() ) )
          for ( auto i{0U}; i<dim; i++ ) (*nit)->Status( prop_key, i, flags[i] );
      return;
    }

    if ( isCorner( boundary ) ) {
      for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ )
        if ( (*nit)->AtBoundary() == boundary )
          for ( auto i{0U}; i<dim; i++ ) (*nit)->Status( prop_key, i, flags[i] );
      return;
    }
  }
  throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags: ",
                         "boundary flag could not be handled." );
} // end InputBoundaryFlags





/**
Uses linear interpolation to interpolate boundary values across model boundary.

Last argument is a vector because there may either be two or four boundary endpoints
depending on whether a 2 or three-dimensional model is used.

@attention This method is about to be deprecated.
*/
template<uint32_t dim>
void Model<dim>::InterpolateBoundaryValues( BOX_BOUNDARY side, const char* input_prop,
                                            const vector<ScalarVariable>& bvalues )
{
  const csmp::Index prop_key = database_.StorageKey( input_prop );

  if ( prop_key.type != SCALAR )
    throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                           "this method can only be applied to scalar variables" );

  if ( prop_key.place != NODE )
    throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                           "method can only be applied to node variables" );

  // checking that the boundary flags are all consistent
  assert( !bvalues.empty() );
  ScalarVariable  res( bvalues[0] ); // getting the variable flag
  for ( auto vit = bvalues.begin(); vit != bvalues.end(); vit++ )
    if ( (*vit).Flag() != res.Flag() )
      throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                             "inconsistent flagging of scalar boundary values" );

  csmp::Region<dim>&  model_domain( this->Region( "Model" ) );

  assert( bvalues.size() >= 2U );
  double  v1 = bvalues[0]();
  double  v2 = bvalues[1]();

  csmp::Point<dim>  xyz_min, xyz_max;
  MinMaxCoordinates( xyz_min, xyz_max );

  // 2D case
  if constexpr ( dim == 2U ) {
    assert( isSide( side ) );
    boundaryMinMaxCoordinates( side, xyz_min, xyz_max );
    for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
      res() = linearInterpolate( make_pair( xyz_min, v1 ), make_pair( xyz_max, v2 ), (*nit)->Coordinate() );
      (*nit)->Store( prop_key, res );
    }
    return;

  } // end 2D

  assert( !isCorner( side ) );

  if ( isSide( side ) ) {
    assert( bvalues.size() == 4U );
    double  v3 = bvalues[2]();
    double  v4 = bvalues[3]();

    switch ( side )
    {
      case LEFT:   // YZ PLANE
        boundaryMinMaxCoordinates( LEFT, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res() = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
          (*nit)->Store( prop_key, res );
        }
        break;

      case RIGHT:  // YZ PLANE
        boundaryMinMaxCoordinates( RIGHT, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res() = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
          (*nit)->Store( prop_key, res );
        }
        break;

      case BACK:   // XY PLANE
        boundaryMinMaxCoordinates( BACK, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res() = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
          (*nit)->Store( prop_key, res );
        }
        break;

      case FRONT:  // XY PLANE
        boundaryMinMaxCoordinates( FRONT, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res() = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
          (*nit)->Store( prop_key, res );
        }
        break;

      case TOP:    // XZ PLANE
        boundaryMinMaxCoordinates( TOP, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res() = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
          (*nit)->Store( prop_key, res );
        }
        break;

      case BOTTOM: // XZ PLANE
        boundaryMinMaxCoordinates( BOTTOM, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res() = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
          (*nit)->Store( prop_key, res );
        }
        break;

      case IRREGULAR:
        // natural neighbor interpolation might solve this case
        throw csmp::Exception( WARNING, "Model::AssignBoundaryValues",
                               "IRREGULAR boundaries are not handled yet. Nothing is done..." );
        break;

      default:
        throw csmp::Exception( ERROR, "Model::AssignBoundaryValues",
                               "type of boundary could not be identified. Nothing is done..." );
    } // end side
  }
  else {  // if edge
    boundaryMinMaxCoordinates( side, xyz_min, xyz_max );
    for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
      res() = linearInterpolate( make_pair( xyz_min, v1 ), make_pair( xyz_max, v2 ), (*nit)->Coordinate() );
      (*nit)->Store( prop_key, res );
    }
  }

} // end InterpolateBoundaryValues (scalar version)






  ///  vector variable version of previous method
template<uint32_t dim>
void Model<dim>::InterpolateBoundaryValues( BOX_BOUNDARY side, const char* input_prop,
                                            const vector<VectorVariable<dim> >& bvalues )
{
  const csmp::Index  prop_key = database_.StorageKey( input_prop );

  if ( prop_key.type == TENSOR )
    throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                           "method has not been implemented for tensor variables" );

  if ( prop_key.place != NODE )
    throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                           "method can only be applied to node variables" );

  // checking that the boundary flags are all consistent
  assert( !bvalues.empty() );
  VectorVariable<dim>  res( bvalues[0] ); // getting the flags
  for ( typename vector<VectorVariable<dim> >::const_iterator
        vit = bvalues.begin(); vit != bvalues.end(); vit++ )
    for ( auto i{0U}; i<dim; i++ )
      if ( (*vit).Flag( i ) != res.Flag( i ) )
        throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                               "inconsistent flagging of boundary vector variables" );

  csmp::Region<dim>&  model_domain( this->Region( "Model" ) );

  assert( bvalues.size() >= 2U );
  VectorVariable<dim>  v1 = bvalues[0];
  VectorVariable<dim>  v2 = bvalues[1];

  csmp::Point<dim>  xyz_min, xyz_max;
  MinMaxCoordinates( xyz_min, xyz_max );

  // 2D case
  if ( dim == 2U ) {
    assert( isSide( side ) );
    boundaryMinMaxCoordinates( side, xyz_min, xyz_max );
    for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
      res( 0 ) = linearInterpolate( make_pair( xyz_min, v1[0] ), make_pair( xyz_max, v2[0] ), (*nit)->Coordinate() );
      res( 1 ) = linearInterpolate( make_pair( xyz_min, v1[1] ), make_pair( xyz_max, v2[1] ), (*nit)->Coordinate() );
      (*nit)->Store( prop_key, res );
    }
    return;

  } // end 2D

  assert( !isCorner( side ) );

  if ( isSide( side ) ) {
    assert( bvalues.size() == 4U );
    VectorVariable<dim>  v3 = bvalues[2];
    VectorVariable<dim>  v4 = bvalues[3];

    switch ( side )
    {
      case LEFT:   // YZ PLANE
        boundaryMinMaxCoordinates( LEFT, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res( 0 ) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
          res( 1 ) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
          res( 2 ) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
          (*nit)->Store( prop_key, res );
        }
        break;

      case RIGHT:  // YZ PLANE
        boundaryMinMaxCoordinates( RIGHT, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res( 0 ) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
          res( 1 ) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
          res( 2 ) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
          (*nit)->Store( prop_key, res );
        }
        break;

      case BACK:   // XY PLANE
        boundaryMinMaxCoordinates( BACK, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res( 0 ) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
          res( 1 ) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
          res( 2 ) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
          (*nit)->Store( prop_key, res );
        }
        break;

      case FRONT:  // XY PLANE
        boundaryMinMaxCoordinates( FRONT, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res( 0 ) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
          res( 1 ) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
          res( 2 ) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
          (*nit)->Store( prop_key, res );
        }
        break;

      case TOP:    // XZ PLANE
        boundaryMinMaxCoordinates( TOP, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res( 0 ) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
          res( 1 ) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
          res( 2 ) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
          (*nit)->Store( prop_key, res );
        }
        break;

      case BOTTOM: // XZ PLANE
        boundaryMinMaxCoordinates( BOTTOM, xyz_min, xyz_max );
        for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
          res( 0 ) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
          res( 1 ) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
          res( 2 ) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
          (*nit)->Store( prop_key, res );
        }
        break;

      case IRREGULAR:
        // natural neighbor interpolation might solve this case
        throw csmp::Exception( WARNING, "Model::AssignBoundaryValues",
                               "IRREGULAR boundaries are not handled yet. Nothing is done..." );
        break;

      default:
        throw csmp::Exception( ERROR, "Model::AssignBoundaryValues",
                               "type of boundary could not be identified. Nothing is done..." );
    } // end side
  }
  else {  // if edge
    boundaryMinMaxCoordinates( side, xyz_min, xyz_max );
    for ( auto nit = model_domain.NodesBegin(); nit != model_domain.NodesEnd(); nit++ ) {
      res( 0 ) = linearInterpolate( make_pair( xyz_min, v1[0] ), make_pair( xyz_max, v2[0] ), (*nit)->Coordinate() );
      res( 1 ) = linearInterpolate( make_pair( xyz_min, v1[1] ), make_pair( xyz_max, v2[1] ), (*nit)->Coordinate() );
      res( 2 ) = linearInterpolate( make_pair( xyz_min, v1[2] ), make_pair( xyz_max, v2[2] ), (*nit)->Coordinate() );
      (*nit)->Store( prop_key, res );
    }
  }

} // end InterpolateBoundaryValues (vector version)



/**
CopyReplace() replaces the values of the variable 'to' with the values
of the variable 'from'. Note, that this operation can only be defined for
variables which have the same type and the same placement. Note also, that
variable flags such as DIRICH are also copied to the target variable.

@param from The variable names of the physical variable (property) which is copied
@param to the variable name that will be replaced.

@section application Application

Use CopyReplace(), for instance to make backups of variable values before
applying time-dependent calculations or to calculate time derivatives
(via self-defined interrelations) from the second before last, the last and
the present variable value.

@section messages Messages

The method will return without modifying any variable values if the type
or the placement of the input variables are not the same. In this case
an error will be reported.
*/
template<uint32_t dim>
void Model<dim>::CopyReplace( const char* from, const char* to )
{
  this->Region( "Model" ).CopyReplace( from, to );
} // end CopyReplace




/**
Provided that property 'a' is a scalar node variable and property 'b' is a
vector variable placed on the element, CopyGradientOfProperty_A_To_B() will
calculate the gradient of property 'a' for each element and assign the
result to the element variable 'b'.

@param prop_a specifies the scalar node variable of which the gradient will be calculated
second string identifies the vector variable placed on the element which
will store the calculated gradient of property 'a'.

@section application Application

The definition of many geophysical or geochemical interrelations requires
a knowledge of property gradients which can be calculated with this method,
using the element interpolation functions of type of finite element which
was used to build the mesh. The method can also be used for post-processing
following the application of algorithms.

@section messages Messages

CopyGradientOfProperty_A_To_B() will return without completing any
calculations, if the input variables do not comply with the specifications
outlined above.
*/
template<uint32_t dim>
bool  Model<dim>::CopyGradientOfProperty_A_To_B( const char* prop_a, const char* prop_b )
{
  return this->Region( "Model" ).CopyGradientOfProperty_A_To_B( prop_a, prop_b );
}





/**
InterpolateNodePropertyToCellProperty() interpolates node property to the
'barycenter' of the triangle. The resulting
value is different from the result obtained by applying the interrelation
subclass NodeToCellProperty. The Calculate() method of the latter assigns
the average value of the 3 element nodes to the element property 'eprop'.

@section arguments Input Arguments

The two string arguments define the name of the node property which is
interpolated and the name of the element property to which the
interpolated value is written.

@section application Application

This method was developed for two-dimensional convection calculations
where the fluid density as element property must be exactly the same for in
the two adjacent triangles which make up a square in regular-gridded meshes.

@section messages Messages

If the variable placement or type of the specified properties fails to
match the specifications outlined above,
InterpolateNodePropertyToCellProperty() will report an error and
return without completing its task.
*/
template<uint32_t dim>
void  Model<dim>::InterpolateNodeToCellProperty( const char* nprop, const char* eprop, bool verbose )
{
  this->Region( "Model" ).InterpolateNodeToCellProperty( nprop, eprop );

  if ( verbose ) {
    cout << "\nModel<" << dim << ">::InterpolateNodeToCellProperty: ";
    cout << "'" << nprop << "' has been successfully interpolated to '" << eprop << "'." << endl;
  }

} // end InterpolateNodeToCellProperty





template<uint32_t dim>
void  Model<dim>::InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* cpprop )
{
  this->Region( "Model" ).InterpolateNodeToIntegrationPointProperty( nprop, cpprop );

  cout << "\nModel<" << dim << ">::InterpolateNodeToIntegrationPointProperty: ";
  cout << "'" << nprop << "' has been successfully interpolated to '" << cpprop << "'." << endl;

} // end InterpolateNodeToIntegrationPointProperty



/**
This method interpolates the specified constraint point property
to the target element property.

@param cprop The name of the constraint point property which shall be interpolated to the
@param eprop element property.

The results are returned to the Model.

@section implementation Implementation

Currently the method averages the IntegrationPoint variable values
to find the element property value.

@section application Application

To visualise constraint point properties as CELL_CENTERED variables in
VTK they have to have a unique value in the element = CELL.

Alternatively, one could integrate the property over the element and
divide the result value by the element area. However, this would work
only for scalar properties.

@section messages Messages

Consistency checks are performed on the placement and type of the
input variables.
*/
template<uint32_t dim>
void  Model<dim>::InterpolateIntegrationPointToCellProperty( const char* cprop, const char* eprop )
{
  this->Region( "Model" ).InterpolateIntegrationPointToCellProperty( cprop, eprop );

  cout << "\nModel<" << dim << ">::InterpolateIntegrationPointToCellProperty: ";
  cout << "'" << cprop << "' has been successfully interpolated to '" << eprop << "'." << endl;

} // end InterpolateIntegrationPointToCellProperty



/**
Using the element interpolation functions, this method extrapolates the
desired constraint point property to the nodes. Constributions of adjacent
elements are averaged but no weighting by element size or proximity of
barycentre to the node is applied.

@param cprop The names of the targeted constraint point
@param nprop and node variables.

The result of the extrapolation is returned into the Model property
storage.

@section implementation Implementation

The method depends on a corresponding function of the FiniteElement
which is not implemented for all finite element types but for all
isoparametric elements. In quadratic isoparametric elements the method
uses a linear extrapolation which is justified for stresses or strains,
for instance, because they define a trilinear field on the
IntegrationPoints.


@section messages Messages

A consistency check on variable type and placement is performed.
*/
template<uint32_t dim>
void  Model<dim>::ExtrapolateIntegrationPointToNodeProperty( const char* cprop, const char* nprop )
{
  this->Region( "Model" ).ExtrapolateIntegrationPointToNodeProperty( cprop, nprop );

  cout << "\nModel<" << dim << ">::ExtrapolateIntegrationPointToNodeProperty: ";
  cout << "'" << cprop << "' has been successfully interpolated to '" << nprop << "'." << endl;

} // end ExtrapolateIntegrationPointToNodeProperty








/**
MinMaxCoordinates() returns the coordinate extrema of a supplied list
of either Node, IntegrationPoint, or Element objects into two supplied
VectorVariable<dim>s. This is equivalent to finding a bounding box which
contains the enlisted objects given that they all lie in the same plane.
If this is not the case, the results refer to a rectangular bounding box.

The two VectorVariable<dim> arguments will store the mininum and maximum
coordinates of the bounding rectangle (2D) or bounding box (3D) which
is defined by the object coordinates.

@section implementation Implementation

Depending on the specified placement, the method loops over the Node,
IntegrationPoint, or Element objects to find their coordinate extrema.
In the case of Element objects, the method uses the BaryCentre (centre
of gravity) of the elements as coordinate point.

@section application Application

Method is used by AssignBoundaryValues().

@section messages Messages

MinMaxCoordinates() will report an error if the placement was input
incorrectly.
*/
template<uint32_t dim>
void Model<dim>::MinMaxCoordinates( Point<dim>& xyz_min,
                                    Point<dim>& xyz_max ) const
{
  this->Region( "Model" ).MinMaxCoordinates( xyz_min, xyz_max );

} // end MinMaxCoordinates





/**
OutputMesh() will print to 'stdout' the mesh connectivity and boundary
flags of each Node, IntegrationPoint, and Element.

@section implementation Implementation

Method calls the Out() interface of the MeshManager.

@section application Application

For model testing, for instance, to test whether the boundaries of a model
have been flagged correctly. Use only for very small meshes since the text
output of large meshes will be prohibitively large.
*/
template<uint32_t dim>
void Model<dim>::Out() const
{
  cout << "\n\n\n\nModel<" << dim << ">::Out: ";
  database_.Out();
  mesh_manager_.Out();
  
  this->RegionsOut();
  this->BoundariesOut();
  this->SplitBoundariesOut();

} // end Out






/**
AssignCellCharacteristicsTo() allows to assign a number of Element
characteristics as identified by strings (second argument) to scalar physical
variables. These characteristics are:

"volume" (3D)<br>
"inner radius"<br>
"aspect ratio" (longest / shortest segment)

@param characteristic The element characteristic is assigned to
@param var a property of the user's choice and as identified by its name in the variable database.

@section implementation Implementation

The method was implemented mainly for two-dimensional calculations.
Therefore, height and width only work in two-dimensional calculations.
The characteristics are obtained from the FiniteElement which is
bridged to the current Element that is queried.

@section application Application

The characteristic "inner radius" is commonly used to find the
appropriate resolution for an advection or visualization grid on which a
variable is to be mapped on.

AssignCellCharacteristicsTo() can be used to:
- test the shape of elements in a mesh for their suitability for a
computation (aspect ratio).
- testing how skewed the elements got by deformation
- integrate element properties in specific calculations, for instance,
if a fracture is just one element wide, the minimum height of
these fracture elements may be equivalent to the fracture aperture.

@section messages Messages

Firstly, the method will report an error and return, if the target
property to which the element characteristic shall be assigned to is not
an element variable.

Further errors may be reported if the user tries to use the characteristic
area in a 3D computation, or volume in a 2D computation. Also, element
height and width are thus far only available in 2D.
*/
template<uint32_t dim>
void Model<dim>::AssignCellCharacteristicsTo( const char* characteristic, const char* var )
{
  this->Region( "Model" ).AssignCellCharacteristicsTo( characteristic, var );

} // end



/**
Displaces the node coordinates by displacements supplied as through the
VECTOR variable (dim = dimensions of the model). For each element, the
nodes of which were displaced, the private boolean variable shape_to_date
is set to false such that its volume is newly calculated
once it is requested.

@param vector_variable The name of the nodal vector variable which holds the node
coordinate displacement.

@section implementation Implementation

After the assignment of the nodal displacements, MoveNodeCoordinatesBy()
will set the Element state variable 'shape_to_data' to false. This will
prompt re-calculation of elemental volumes if these are queried in
subsequent computations.

@section application Application

If a mesh shall be deformed using the diplacements of a deformation
calculation stored in a vector variable, ChangeNodeCoordinatesTo() can
be used displace the node coordinates by these displacements.

@section messages Messages

Due to the total garbage results that may arise,
MoveNodeCoordinatesBy() will halt the simulation reporting a fatal
error, if the target property is node a node or vector type property.
*/
template<uint32_t dim>
void Model<dim>::MoveNodeCoordinatesBy( const char* vector_variable )
{
  this->Region( "Model" ).MoveNodeCoordinatesBy( vector_variable );
}


/**
Assigns node coordinates to the supplied node variable which may be either
of scalar or vector type. If either the X, Y, or Z coordinate are to be
assigned to a scalar variable an index from 0...dim-1 must be supplied
to indicate which coordinate axis shall be mapped.

@param vector_variable The name of the physical variable to which the node coordinate shall
be assigned to, and, if a scalar property shall be initialized with
coordinate values, an int index to indicate which coordinate axis is
desired (0=X, 1=Y, 2=Z).

@section implementation Implementation

Method uses the Coordinates() and X(), Y(), and Z() interfaces of the
Node class.

@section application Application

To arrive at coordinate values for a model where topography is
important etc.

@section messages Messages

The method will report a fatal error and terminate the computation, if
the target variable is not placed on the node or if it is not of the
correct type.
*/
template<uint32_t dim>
void Model<dim>::AssignNodeCoordinatesTo( const char* vector_variable )
{
  this->Region( "Model" ).AssignNodeCoordinatesTo( vector_variable );
}



template<uint32_t dim>
void Model<dim>::AssignNodeCoordinatesTo( const char* scalar_variable, char c )
{
  this->Region( "Model" ).AssignNodeCoordinatesTo( scalar_variable, c );
}



/**
An Algorithm object is passed either directly to the Model or to a
subdomain of it to which it applies. In this process, the Algorithm is
provided with iterators to the target collection of (default = Elements)
to accumulate data, assign boundary conditions and to carry out a
global matrix inversion using the 'Solver' or any other calculation tool.

@param problem An initialized Algorithm object, i.e. an algorithm to which partial
differential operators have been assigned (see MathOperator... objects).

@section implementation Implementation

The Apply() method calls on the interface of the Algorithm object to carry
out a sequence of operations which will build, initialize, and invert
a global solution matrix in order to find distributed values of the
dependent variable(s) of interest. In detail, this sequence will
vary, dependent on whether the calculation establishes transient or
steady-state variable values. The more complicated transient calculation
will call (see interface documentation of the Algorithm class):


@code
EstablishMatrixSetup()       to allocate the solution matrix etc.
Accumulate()                 to assemble the solution matrix.
AssignInitialConditions()    to put the present model state into the
righthand vector
LateAccumulate()             to add source or sink rates into the
righthandside
AssignEssentialConditions()  assign Neumann or Dirichlet conditions
Solve()                      invert the solution matrix with the Solver
OutputResults()              test and map the solution vector back into
the Model / MemoryManager storage
PostProcess()                calculate further variable values on the
basis of the calculation results
@endcode


@section application Application

Standard interface for 2 and 3D global computations.

@section messages Messages

The interface methods of the Algorithm class Accumulate() and Solve()
will report the progress of a computation in a fashion that depends on
the type of problem which is solved and on the solution technique. Please
refer to the documentation of the Algorithm interface and the Solver to
learn more about this specific output.
*/
template<uint32_t dim>
void Model<dim>::Apply( PDE_Integrator<dim, csmp::Region>& problem, bool debug )
{
  problem.IntegrateOver( this->Region( "Model" ), debug );

}


/// application to all boundaries
template<uint32_t dim>
void Model<dim>::Apply( PDE_Integrator<dim, csmp::Boundary>& problem, bool debug )
{
  for ( typename map<string, csmp::Boundary<dim> >::iterator it = this->BoundariesBegin(); it != this->BoundariesEnd(); ++it )
    problem.IntegrateOver( (*it).second, debug );
}



/**
specific regions
*/
template<uint32_t dim>
void Model<dim>::Apply( PDE_Integrator<dim, csmp::Region>& problem, const char* region_name, bool debug )
{
  problem.IntegrateOver( this->Region( region_name ), debug );

}


/// PDE solution applied to a specific boundary
template<uint32_t dim>
void Model<dim>::Apply( PDE_Integrator<dim, csmp::Boundary>& problem, const string& boundary_name, bool debug )
{
  problem.IntegrateOver( this->Boundary( boundary_name ), debug );

} // Apply (PDE_Integrator)








/**
Applies a calculation recipe which is specified in an Interrelation
subclass to every element/face/node etc. in the current Model.
If the Interrelation applies to a single Region (subregion of the Model)
it is passed on to the interface of the corresponding Region object.

@param relation A user-defined subclass of an Interrelation base class.

@section implementation Implementation

Executes the method Apply() on interrelation objects derived from the
Interrelation base class. Apply() will in turn call Calculate() in the
class derived from the Interrelation base class (subclass). Remember that
Interrelations are defined as calculations that do not require the
inversion of a global matrix.

@section application Application

Examples of interrelations include the calculation of strains from
displacements, the horizontal and vertical components of the
fluid-flow velocity from nodal pressures and/or elemental permeabilities
and porosities.

@section messages Messages

Please consult the documentation of the Interrelation class and the
derived subclasses for a description of the output which results when
these objects are applied to the Model or a Region.
*/
template<uint32_t dim>
void Model<dim>::Apply( Interrelation<dim>& relation, const char* region )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  if ( this->uniqueRegionMap_.empty() ) {
    csmp_error.Note( ERROR, "Model<dim>::Apply",
                       "cannot apply Interrelation because there is no unique group in model" );
    return;
  }

  // if "Model" is the target
  if ( strcmp( region, "Model" ) == 0 ) {
    relation.Apply( *this );
    return;
  }

  // if "Model" is the only unique group, the interrelation is applied to it
  // else it is passed to all unique regions
  if ( this->uniqueRegionMap_.find( region ) != this->uniqueRegionMap_.end() ) {
    this->Region( region ).Apply( relation );
    return;
  }

  for ( auto& it : this->uniqueRegionMap_ )
    it.second.Apply( relation );

} // end Apply (Interrelation)




/**
Prints the values and flags of the distributed variable of interest on
stdout.

@param prop The name of the physical variable of interest.

@section application Application

To check variable values in small test problems.

@section messages Messages

OutputVariableToScreen() will report the variable type, its placement,
and property index. Then it will print the host object IDs followed
by the variable flags and values.
*/
template<uint32_t dim>
void Model<dim>::OutputVariableToScreen( const char* prop ) const
{
  this->Region( "Model" ).OutputVariableToScreen( prop );

} // end OutputVariableToScreen












/**
Finds the minimum and maximum values of a physical variable which is
distributed over the current finite-element mesh. If the
property is a vector variable, the vector length, and if it is a
tensor variable, the tensors determinant will be reported.

@param prop The name of the physical variable

The results of the range search will be returned into the second and third
method arguments. If the variable has not been initialized, the value
NAN (not a number) is returned.

@section application Application

For example, MinMaxOf() is useful to find maximum property values which
do not depend on Dirichlet boundary conditions but on applied sources
and sinks.

@section messages Messages

MinMaxOf() will report an error if the property of interest is unknown
to the PropertyDatabase.
*/
template<uint32_t dim>
void Model<dim>::MinMaxOf( const char* prop, double& vmin, double& vmax ) const
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !database_.IsDefined( prop ) ) {
    csmp_error.Note( ERROR, "Model<dim>::MinMaxOf", prop, "is undefined; nothing could be done" );
    return;
  }
  csmp::Index  prop_key( database_.StorageKey( prop ) );

  // properties / variables placed on the model
  if ( prop_key.place == MODEL ) {
    if ( prop_key.type == SCALAR ) vmin = vmax = this->Read( prop_key );
    else if ( prop_key.type == VECTOR ) {
      VectorVariable<dim>  vc;
      this->Read( prop_key, vc );
      vmin = vmax = vc.Length();
    }
    else if ( prop_key.type == TENSOR ) {
      TensorVariable<dim>  ts;
      this->Read( prop_key, ts );
      // recovering and sorting to find minimum and maximum Eigen values
      minMaxEigenValues( ts, vmin, vmax );
    }
    else if ( prop_key.type == ARRAY ) {
      ArrayVariable  a( prop_key.dataDepth );
      this->Read( prop_key, a );
      a.MinMax( vmin, vmax );
    }
    else if ( prop_key.type == FLAGGEDARRAY ) {
      FlaggedArrayVariable  fa( prop_key.dataDepth );
      this->Read( prop_key, fa );
      fa.MinMax( vmin, vmax );
    }
    return;
  }


  // standard properties present everywhere in the model
  if ( prop_key.place == NODE || prop_key.place == ELEMENT || prop_key.place == ELEMENT_INTEGRATION_POINT ||
       prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == SECTOR_INTEGRATION_POINT ) {
    this->Region( "Model" ).MinMaxOf( prop_key, vmin, vmax );
    return;
  }


  // unique and non-unique regions 
  if ( prop_key.place == REGION ) {
    typename map<string, csmp::Region<dim> >::const_iterator  git( this->UniqueRegionsBegin() );
    (*git).second.MinMaxOf( prop_key, vmin, vmax );
    double  gmin( vmin ), gmax( vmax );
    while ( git != this->UniqueRegionsEnd() ) {
      (*git).second.MinMaxOf( prop_key, vmin, vmax );
      gmin = min( gmin, vmin );
      gmax = max( gmax, vmax );
      git++;
    }
    for ( auto ngit = this->RegionsBegin(); ngit != this->RegionsEnd(); ngit++ ) {
      (*ngit).second.MinMaxOf( prop_key, vmin, vmax );
      gmin = min( gmin, vmin );
      gmax = max( gmax, vmax );
    }
    vmin = gmin;
    vmax = gmax;
    return;
  }


  // boundaries
  if ( prop_key.place == BOUNDARY || prop_key.place == FACE || prop_key.place == FACE_INTEGRATION_POINT ||
       prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT ) {
    typename map<string, csmp::Boundary<dim> >::const_iterator  git( this->BoundariesBegin() );
    (*git).second.MinMaxOf( prop_key, vmin, vmax );
    double  gmin( vmin ), gmax( vmax );
    while ( git != this->BoundariesEnd() ) {
      (*git).second.MinMaxOf( prop_key, vmin, vmax );
      gmin = min( gmin, vmin );
      gmax = max( gmax, vmax );
      git++;
    }
    vmin = gmin;
    vmax = gmax;
    return;
  }


  // split boundaries
  if ( prop_key.place == SPLIT_BOUNDARY || prop_key.place == INTER_FACE || prop_key.place == INTER_FACE_INTEGRATION_POINT ||
       prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
    typename map<string, csmp::SplitBoundary<dim> >::const_iterator  git( this->SplitBoundariesBegin() );
    (*git).second.MinMaxOf( prop_key, vmin, vmax );
    double  gmin( vmin ), gmax( vmax );
    while ( git != this->SplitBoundariesEnd() ) {
      (*git).second.MinMaxOf( prop_key, vmin, vmax );
      gmin = min( gmin, vmin );
      gmax = max( gmax, vmax );
      git++;
    }
    vmin = gmin;
    vmax = gmax;
  }

} // end MinMaxOf




/**
      For vector and tensor variables, this method uses the length and the range of Eigenvalues as a measure.
      For array variables, the min and max values (L1 norm) are used.
*/
template<uint32_t dim>
bool Model<dim>::IsWithinRange( const char* var_name, const char* model_subdomain ) const
 {
    double  omin(DBL_MAX), omax(DBL_MIN), pmin, pmax;
    
    // is the subdomain a region, boundary, or split boundary?
    if ( this->ContainsRegion( string(model_subdomain) ) )
      this->Region(model_subdomain).MinMaxOf( var_name, omin, omax );
    else if ( this->ContainsBoundary( string(model_subdomain) ) )
      this->Boundary(model_subdomain).MinMaxOf( var_name, omin, omax );
    else if ( this->ContainsSplitBoundary( string(model_subdomain) ) )
      this->SplitBoundary(model_subdomain).MinMaxOf( var_name, omin, omax );    
    
    Database().RangeOf(  var_name, pmin, pmax );
    if ( omin >= pmin && omax <= pmax ) return true;

    return false;
    
 } // end IsWithinRange





  /// Instantiates MODEL LocalVariables
template<uint32_t dim>
void Model<dim>::InitializeLocalVariableStorage()
{
  this->ResizePropertyStorage( Database().LocalVariablesAt( MODEL ) );
}


/**
Property storage for the geometric primitives that form the mesh is created
before Regions, Boundaries, and / or SplitBoundaries are formed.

This method fulfills the job at a later stage, updating only if such
variables exist.
*/
template<uint32_t dim>
bool Model<dim>::UpdateSubdomainPropertyStorage()
{
  bool  increased_storage( false );

  // 1. Regions
  if ( this->Database().VariableCount( REGION ) > 0 ) {
    for ( typename map<string, csmp::Region<dim> >::iterator
          git = this->UniqueRegionsBegin(); git != this->UniqueRegionsEnd(); git++ )
      (*git).second.ResizePropertyStorage( this->Database().LocalVariablesAt( REGION ) );
    for ( typename map<string, csmp::Region<dim> >::iterator
          git = this->RegionsBegin(); git != this->RegionsEnd(); git++ )
      (*git).second.ResizePropertyStorage( this->Database().LocalVariablesAt( REGION ) );
    increased_storage = true;
  }

  // 2. Boundaries
  if ( this->Database().VariableCount( BOUNDARY ) > 0 ) {
    for ( typename map<string, csmp::Boundary<dim> >::iterator
          git = this->BoundariesBegin(); git != this->BoundariesEnd(); git++ )
      (*git).second.ResizePropertyStorage( this->Database().LocalVariablesAt( BOUNDARY ) );
  }

  // 3. SplitBoundaries
  if ( this->Database().VariableCount( SPLIT_BOUNDARY ) > 0 ) {
    for ( typename map<string, csmp::SplitBoundary<dim> >::iterator
          git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); git++ )
      (*git).second.ResizePropertyStorage( this->Database().LocalVariablesAt( SPLIT_BOUNDARY ) );
  }

  return increased_storage;

}  // end InitializeRegionProperties





/**
Database() comes in a constant and a volatile version,
giving you access to the PropertyDatabase object inside the Model.

@return Either a constant or a volatile reference to the PropertyDatabase object.
*/
template<uint32_t dim>
PropertyDatabase<dim>&  Model<dim>::Database() { return database_; }

template<uint32_t dim>
const PropertyDatabase<dim>& Model<dim>::Database() const { return database_; }

template<uint32_t dim>
const FiniteElementManager&  Model<dim>::FE_Manager() const
{ return mesh_manager_.FiniteElements(); }



/**
Mesh() comes in a constant and in a volatile version, giving you
access to the MeshManager object inside of the Model object.

@return Either a constant or a volatile reference to the MeshManager.
*/
template<uint32_t dim>
MeshManager<dim>&  Model<dim>::Mesh() { return mesh_manager_; }

template<uint32_t dim>
const MeshManager<dim>&  Model<dim>::Mesh() const { return mesh_manager_; }








// ================================================================================================================
//
//            NEW BINARY FILE INPUT / OUTPUT
//
// ================================================================================================================


/**
writes entire model with associated properties / variables to CSMP native set of binary files.
*/
template<uint32_t dim>
void Model<dim>::OutputToBinaryFile( const char* file_string )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 0. reporting
  double& model_time( ModelTime::Instance().modelTime );
  cout << "\nModel<" << dim << ">::OutputToBinaryFile: Saving model '" << Name();
  cout << "' at current time level, t = " << model_time << " secs." << endl;

  // 1. mesh output: creating a VSet including Face and InterFace objects
  //   (elements, faces, and interfaces are numbered in a single continous sequence)
  VSet<dim>  vset;
  mesh_manager_.OutputMeshTo( vset );

  // 2. property output into VSet including Face and InterFace data
  mesh_manager_.OutputStoredVariablesTo( Database(), vset );

  // adding properties stored on "Model"
  map<string, Index>  properties;
  Database().ListVariables( MODEL, properties );

  // for all model properties
  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      // setting the specifications for the property storage (no memory allocation yet)
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( (*pit).second.flagDepth );
      const size_t data_capacity( (*pit).second.dataDepth );
      // allocating memory to store the property flags and values
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
        {
          case SCALAR: {
                ScalarVariable value;
                this->Read( (*pit).second, value );
                pushBack( data, value );
              }
                break;
          case VECTOR: {
                VectorVariable<dim> value;
                this->Read( (*pit).second, value );
                pushBack( data, value );
              }
            break;
          case TENSOR: {
                TensorVariable<dim> value;
                this->Read( (*pit).second, value );
                pushBack( data, value );
              }
            break;
          case ARRAY: {
                ArrayVariable value;
                this->Read( (*pit).second, value );
                pushBack( data, value );
              }
            break;
          case FLAGGEDARRAY: {
                FlaggedArrayVariable value;
                this->Read( (*pit).second, value );
                pushBack( data, value );
              }
            break;
          default:
            csmp_error.Note( ERROR, "Model<dim>::OutputToBinaryFile:",
                               (*pit).first, "type of Model variable not recognized." );
        }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }

  // writing the VSet to binary file
  vset.OutputTo( BinaryVsetFileName( file_string ).c_str(), model_time );

  // 3. regions: unique and then the non-unique regions
  this->OutputRegionsToBinary( BinaryRegionsFileName( file_string ).c_str() );

  // 4. boundaries
  this->OutputBoundariesToBinary( BinaryBoundariesFileName( file_string ).c_str() );

  // 5. splitboundaries
  this->OutputSplitBoundariesToBinary( BinarySplitBoundariesFileName(file_string).c_str() );
  
  // node manifolds are deduced from the connectivty stored in the VSet and handled inside the MeshManager

  // 6. variable specifications through the database
  Database().BinaryOut( BinaryVariablesFileName( file_string ).c_str() );
  cout << "\nModel<" << dim << ">::OutputToBinaryFile: Output of model '" << Name();
  cout << "' to CSMP binaries completed successfully.\n\n";

} // end OutputToBinaryFile




/**
     Reads model written by OutputToBinaryFile() including all associated properties or a subset of variables
     
     @param model_name of the saved model that is to be retrieved from the binary file
     @param subset_variables selected variables that shall be read from the binary file; if empty, all variables are read
     
     @attention If the variable subset is empty and variables are encoutered in the VSet which are not yet contained in the current database,
     extra storage is created for these.
     
     @attention if any of the variables in the subset is not known to the database an error is reported.
*/
template<uint32_t dim>
void Model<dim>::InputFromBinaryFile( const char* model_name, const set<string>& subset_variables )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. The binary data are read into VSet
  double& model_time( ModelTime::Instance().modelTime );
  VSet<dim>  vset;
  cout << "\nModel<" << dim;
  cout << ">::InputFromBinaryFile: Reading '" << model_name;
  cout << "' from VSet... " << endl;
  vset.InputFrom( BinaryVsetFileName( model_name ).c_str(), model_time, subset_variables );

  // 2.1 If there are variables with meaningful values stored in the VSet which are not contained in the database
  //     new entries are created for these (this option applies only if the subset is empty).
  if ( subset_variables.empty() ) {
       int counter{0};
       for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
         if ( Database().IsDefined( (*pit).first.c_str() ) == false ) {
              double vmin, vmax;
              string notation{"abbrev"}; notation += counter++;
              (*pit).second.MinMaxOf( vmin, vmax );
              // if we are dealing with a range of real numbers
              if ( !isnan(vmin) && !isnan(vmax) && vmin < vmax )
                // prop name, prop unit, VARIABLE_TYPE, PLACEMENT, vsize, vmin, vmax, usage
                Database().AddProperty( (*pit).first.c_str(), notation.c_str(), "SI",
                                        (*pit).second.Type(), (*pit).second.Placement(), (*pit).second.Components(),
                                         vmin, vmax, "property read from VSet" );
           }
    }

  // 2.2 build finite element mesh and associated property storage, initialising 'mtrl' identifiers and boundary flags
  const bool with_FV_variables = (finiteVolumeVariables(Database()) > 0 ) ? true : false;
  mesh_manager_.Initialize( database_, vset, with_FV_variables );

  // 3. assigning properties to mesh (this reads in the properties output to file via Region::OutputTo(VSet) )
  if ( !vset.DataEmpty() ) mesh_manager_.InputStoredVariablesFrom( Database(), vset );
  else ErrorHandler::Instance().Note( INFO, "Model<dim>::InputFromBinaryFile:", "No properties found in VSet." );

  // 4. properties and values stored on the model itself
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    // apart from the name string key in the map, PropertyData contains the most important variable specifications
    if ( (*pit).second.Placement() != MODEL ) continue;
    // if we are only considering variables in the supplied subset
    if ( !subset_variables.empty() && subset_variables.find( (*pit).first ) == subset_variables.end() ) continue;
    // more checks
    assert( Database().IsDefined( (*pit).first.c_str() ) );
    const csmp::Index key( Database().StorageKey( (*pit).first.c_str() ) );
    assert( key.place == MODEL );
    assert( (*pit).second.Size() / key.dataDepth == 1U );

    switch ( key.type )
      {
        case SCALAR: {
              ScalarVariable value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        case VECTOR: {
              VectorVariable<dim> value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        case TENSOR: {
              TensorVariable<dim> value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        case ARRAY: {
              ArrayVariable value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        case FLAGGEDARRAY: {
              FlaggedArrayVariable value;
              read( (*pit).second, 0, value );
              this->Store( key, value );
            }
          break;
        default:
          csmp_error.Note( ERROR, "Model<dim>::InputFromBinaryFile:",
                             (*pit).first, "type of Model variable not recognized." );
      }
  }

  // 5. reconstruction of the regions
  this->InputRegionsFromBinary( BinaryRegionsFileName( model_name ).c_str(), subset_variables );

  // making sure that the computational region has been built
  if ( mesh_manager_.IsContiguous() )
    if ( !this->ContainsRegion( "Model" ) )
      throw csmp::Exception( ERROR, "Model<>::InputFromBinaryFile", "Root region 'Model' is not present." );

  // 6. reconstructing the boundaries, if any
  if ( vset.Faces() > 0 )
    this->InputBoundariesFromBinary( BinaryBoundariesFileName( model_name ).c_str(), subset_variables );
  
  // 7. reconstructing the splitboundaries, if any
  if ( vset.Interfaces() > 0 )
    this->InputSplitBoundariesFromBinary( BinarySplitBoundariesFileName(model_name).c_str(), subset_variables );

  cout << "\nModel<" << dim << ">::InputFromBinaryFile: input from binaries (file set: " << model_name << ") completed successfully.\n\n";

} // end InputFromBinaryFile







  // ================================================================================================================
  //
  //            GLOBAL FUNCTIONS INVOLVING MODEL
  //
  // ================================================================================================================




/**
Calculates the x, y, z extent of the model and returns these lengths into its arguments.
The return value is a string that contains the dimensions with explanations.
*/
string  boundingBox( const Model<3U>& sg, double& dim_x, double& dim_y, double& dim_z )
{
  Point<3U>  xyz_min, xyz_max;
  sg.MinMaxCoordinates( xyz_min, xyz_max );
  dim_x = xyz_max[0] - xyz_min[0];
  dim_y = xyz_max[1] - xyz_min[1];
  dim_z = xyz_max[2] - xyz_min[2];

  string  dimensions;
  char num[30];
  dimensions += "x-length (m): ";
  sprintf( num, "%lf", dim_x );
  dimensions += num;
  dimensions += ",  ";
  dimensions += "y-length: ";
  sprintf( num, "%lf", dim_y );
  dimensions += num;
  dimensions += ",  ";
  dimensions += "z-length: ";
  sprintf( num, "%lf", dim_z );
  dimensions += num;
  dimensions += " (box-shaped model). ";

  return dimensions;

} // end boundingBox



/**

Method obtains the range of values for the physical variable that it
is prompted for by the user. For vector and tensor variables, the
length and the minimum/maximum eigenvalues are returned (check whether
the latter actually happens).

With the optional boolean argument (default=true), the user can determin
the return value. As the default, the maximum obtained value is returned;
else the minimum.

@section arguments Input Arguments

The model to be examined, the name of the target physical variable, and
targeted return value (true->maximum, false->minimum of target variable).

The second version of this method also takes an I/O handler as argument
in order to log the calculated values to file etc.

@return Either the maximum (default) or the minimum value of the target variable.
The result is printed to the screen.

*/
template<uint32_t  dim>
double printRangeOfVariable( const Model<dim>& sg,
                               const char* var, bool max_or_min )
 {
     double pmin, pmax;
     sg.MinMaxOf( var, pmin, pmax );
     cout << scientific << setprecision(5) <<"\nRange of variable ["<< sg.Database().Unit(var) <<"]: '";
     cout << var <<"': "<< pmin <<" to "<< pmax << endl;
          
     if ( !max_or_min ) return pmin;
     return pmax;
 }



template<uint32_t  dim>
double printRangeOfVariable( const Model<dim>& sg,
                               Standard_IO_Handler& io,
                               const char* var, bool max_or_min )
 {
     double  pmin, pmax;
     sg.MinMaxOf( var, pmin, pmax );
     cout << scientific << setprecision(5) <<"\nRange of variable ["<< sg.Database().Unit(var) <<"]: '";
     cout << var <<"': "<< pmin <<" to "<< pmax << endl;
     
     // recording the measured variable value range at given timestep
     double& model_time( ModelTime::Instance().modelTime );
     char   info[100];
     sprintf( info, "%lf", model_time );
     string var_info(info);
     var_info += " secs, range of'";
     var_info += var;
     var_info += "' [";
     var_info += sg.Database().Unit(var);
     var_info += "]: ";
     sprintf( info, "%lf", pmin );
     var_info += info;
     var_info += " to ";
     sprintf( info, "%lf", pmax );
     var_info += info;
     
     io.RecordInformation( var_info );
     
     if ( !max_or_min ) return pmin;
     return pmax;
 }




/// as above, but for individual model regions
template<uint32_t  dim>
double printRangeOfVariable( const Model<dim>& sg,
                               const char* group, const char* var, bool max_or_min )
 {
     double  pmin, pmax;
     const PropertyDatabase<dim>& p_ref = sg.Database();
     const PLACEMENT place = sg.Database().Placement(var);
     
     if ( sg.ContainsRegion(group) && !faceVariable(place) && !interFaceVariable(place) )
       sg.Region( group ).MinMaxOf( var, pmin, pmax );
     else if ( sg.ContainsBoundary(group) ) sg.Boundary( group ).MinMaxOf( var, pmin, pmax );
     else if ( sg.ContainsSplitBoundary(group) ) sg.SplitBoundary( group ).MinMaxOf( var, pmin, pmax );
     else {
          cerr <<"\nprintRangeOfVariable: '"<< group <<"' does not exist."<< endl;
          return numeric_limits<double>::signaling_NaN();
       }
     cout << scientific << setprecision(5) <<"\nRange of variable ["<< p_ref.Unit(var) <<"]: '";
     cout << var <<"' in subdomain of model '"<< group <<"': "<< pmin <<" to "<< pmax << endl;
          
     if ( !max_or_min ) return pmin;
     return pmax;
 }




template<uint32_t  dim>
double printRangeOfVariable( const Model<dim>& sg,
                                Standard_IO_Handler& io,
                                const char* group,
                                const char* var, bool max_or_min )
 {
     double& model_time( ModelTime::Instance().modelTime );
     double         pmin, pmax;
     const PropertyDatabase<dim>& p_ref = sg.Database();
     const PLACEMENT place = sg.Database().Placement(var);

     if ( sg.ContainsRegion(group) && !faceVariable(place) && !interFaceVariable(place) )
       sg.Region( group ).MinMaxOf( var, pmin, pmax );
     else if ( sg.ContainsBoundary(group) ) sg.Boundary( group ).MinMaxOf( var, pmin, pmax );
     else if ( sg.ContainsSplitBoundary(group) ) sg.SplitBoundary( group ).MinMaxOf( var, pmin, pmax );
     else {
          cerr <<"\nprintRangeOfVariable: '"<< group <<"' does not exist."<< endl;
          return numeric_limits<double>::signaling_NaN();
       }
     cout << scientific << setprecision(5) <<"\nRange of variable ["<< p_ref.Unit(var) <<"]: '";
     cout << var <<"': "<< pmin <<" to "<< pmax <<" in subdomain of model '"<< group <<"'"<< endl;
     
     // recording the measured variable value range at given timestep
     char info[100];
     sprintf( info, "%lf", model_time );
     string var_info(info);
     var_info += info;
     var_info += ", region: ";
     var_info += group;
     var_info += ", secs, range of '";
     var_info += var;
     var_info += "' [";
     var_info += p_ref.Unit(var);
     var_info += "]: ";
     sprintf( info, "%lf", pmin );
     var_info += info;
     var_info += " to ";
     sprintf( info, "%lf", pmax );
     var_info += info;
     
     io.RecordInformation( var_info );
     
     if ( !max_or_min ) return pmin;
     return pmax;
 }






/**

Method prints the physical dimensions of the model and returns either
the maximum or the intermediate axis, depending on the value of its
second argument.

returns intermediate (true) or maximum (false) model dimensions.

@section arguments Input Arguments

The current model and a boolean flag. For 'true' the intermediate axis
is returned, if 'false' the long axis is returned.

@return The intermediate or long axis of the current model.
*/
template<uint32_t  dim>
double  printModelDimensions( const Model<dim>& sg, bool intermed_or_max )
 {
    Point<dim> xyz_min, xyz_max;
    sg.MinMaxCoordinates( xyz_min, xyz_max );
    cout <<"\nprintModelDimensions: Dimensions of model (meters): "<< endl;
    cout <<"xmin, xmax (horizontal right):    "<< xyz_min[0] <<" "<< xyz_max[0] << endl;
    if constexpr ( dim != 1U ) cout <<"ymin, ymax (vertical upward):     "<< xyz_min[1] <<" "<< xyz_max[1] << endl;
    if constexpr ( dim == 3U ) cout <<"zmin, zmax (horizontal to front): "<< xyz_min[2] <<" "<< xyz_max[2] << endl << endl;

    set<double,greater<double> >  axis;
    axis.insert( xyz_max[0] - xyz_min[0] );
    if constexpr ( dim != 1U ) axis.insert( xyz_max[1] - xyz_min[1] );
    if constexpr ( dim == 3U ) axis.insert( xyz_max[2] - xyz_min[2] );
    
    set<double,greater<double> >::const_iterator  it = axis.begin();
    
    if ( !intermed_or_max ) return *it;
    
    if ( axis.size() >= 2U ) it++;
    
    return *it;

 } // end printModelDimensions




/**
    calculates the center of gravity of the model by averaging
    the barycenter locations of all highest-dimensional elements.
*/
template<uint32_t  dim>
Point<dim>  centerOfGravity( const Model<dim>& model )
 {
    const Region<dim>& mref(model.Region("Model"));
    auto it(mref.CellsBegin());
    Point<dim>  center((*it)->BaryCenter());
    double    counter(0.);
    it++;
   
    while( it != mref.CellsEnd() ) {
         if ( dim == 3U ) {
               if ( (*it)->FE()->IsVolume() ) {
                    center += (*it)->BaryCenter();
                    counter += 1.;
                 }
            }
         else if ( dim == 2U ) {
               if ( (*it)->FE()->IsSurface() ) {
                    center += (*it)->BaryCenter();
                    counter += 1.;
                 }
            }
         else /* 1D */ {
                    center += (*it)->BaryCenter();
                    counter += 1.;
            }
         it++;
      }
    center /= counter;
    return center;

 } // end CenterOfGravity

template Point<1U>  centerOfGravity( const Model<1U>& );
template Point<2U>  centerOfGravity( const Model<2U>& );
template Point<3U>  centerOfGravity( const Model<3U>& );



// template instantiations
template
double  printModelDimensions( const Model<1U>& sg, bool intermed_or_max );

template
double  printRangeOfVariable( const Model<1U>& sg,
                                const char* var, bool max_or_min );
template
double  printRangeOfVariable( const Model<1U>& sg,
	                              Standard_IO_Handler& io, const char* var,
	                              bool max_instead_of_min );
template
double  printRangeOfVariable( const Model<1U>& sg,
                                const char* group, const char* var, bool max_or_min );
template
double  printRangeOfVariable( const Model<1U>& sg,
                                Standard_IO_Handler& io,
                                const char* group, const char* var,
                                bool max_instead_of_min );

template
double  printModelDimensions( const Model<2U>& sg, bool intermed_or_max );

template
double  printRangeOfVariable( const Model<2U>& sg,
                                const char* var, bool max_or_min );
template
double  printRangeOfVariable( const Model<2U>& sg,
                                Standard_IO_Handler& io, const char* var,
                                bool max_instead_of_min );
template
double  printRangeOfVariable( const Model<2U>& sg,
                                const char* group, const char* var, bool max_or_min );
template
double  printRangeOfVariable( const Model<2U>& sg,
                                Standard_IO_Handler& io,
                                const char* group, const char* var,
                                bool max_instead_of_min );

template
double  printModelDimensions( const Model<3U>& sg, bool intermed_or_max );

template
double  printRangeOfVariable( const Model<3U>& sg,
                                const char* var, bool max_or_min );
template
double  printRangeOfVariable( const Model<3U>& sg,
                                Standard_IO_Handler& io, const char* var,
                                bool max_instead_of_min );
template
double  printRangeOfVariable( const Model<3U>& sg,
                                const char* group, const char* var, bool max_or_min );
template
double  printRangeOfVariable( const Model<3U>& sg,
                                Standard_IO_Handler& io,
                                const char* group, const char* var,
                                bool max_instead_of_min );







/**
     Smoothes scalar element variable by extrapolating it to the nodes and back-interpolating it to barycenters
     Apart from the mname of element variable to be smoothed, the name of the temporary node variable needs to be specified.
     Uses the dummy variable 'dummy node' to store the interim result
*/
template<uint32_t dim>
void smoothElementVariable( Model<dim>& model, const char* region, const char* element_var, const char* temp_node_var, uint32_t n_smoothing_cycles )
 {
    csmp::Index eprop_key = model.Database().StorageKey(element_var);
    csmp::Index nprop_key = model.Database().StorageKey(temp_node_var);
    assert( nprop_key.place == NODE );

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    if ( eprop_key.place != ELEMENT ) {
         csmp_error.Note( ERROR, "smoothElementVariable", "Smoothed variable must be placed on the element" );
         return;
      }
    if ( nprop_key.place != NODE ) {
         csmp_error.Note( ERROR, "smoothElementVariable", "Temporary variable must be placed on the node" );
         return;
      }
    if ( nprop_key.type != eprop_key.type ) {
         csmp_error.Note( ERROR, "smoothElementVariable", "Smoothed and temporary variable must have the same type" );
         return;
      }
    if ( n_smoothing_cycles == 0 ) {
         csmp_error.Note( WARNING, "smoothElementVariable", "smoothing cycles=0; nothing was done" );
         return;
      }
   
    // smoothing
    Region<dim>& ref = model.Region(region);
   
    for ( auto i{0U}; i<n_smoothing_cycles; i++ ) {
         ref.ExtrapolateCellToNodeProperty( element_var, temp_node_var );
         ref.InterpolateNodeToCellProperty( temp_node_var, element_var );
      }

 } // end smoothElementVariable

template void smoothElementVariable( Model<1U>&, const char*, const char*, const char*, uint32_t );
template void smoothElementVariable( Model<2U>&, const char*, const char*, const char*, uint32_t );
template void smoothElementVariable( Model<3U>&, const char*, const char*, const char*, uint32_t );




/**

Randomly perturbs the values of a scalar target property by subtracting an
amount which varies between minus zero and the specified percentage of the
original maximum value of the target property.

@section arguments Input Arguments

RandomPerturb() requires the name of the property which shall be perturbed
and the percentage of the original maximum value of the property by which
the property shall be perturbed, in order to operate.

@section application Application

Processes which are critically dependent on initial conditions can profit
from a 'noisy' input signal when one tries to simulate natural behaviour.

@section messages Messages

RandomPerturb() only handles scalar variables and it will therefore report
an error and return without executing when one tries to perturb a vector or
tensor variable.

*/
template<uint32_t dim>
void randomPerturb( Model<dim>& sg, const char* prop, double by_percent_of_max_value )
 {
    csmp::Index prop_key = sg.Database().StorageKey(prop);
    Region<dim>&  sgroup(sg.Region("Model"));
    
    if ( prop_key.type != SCALAR )
      throw csmp::Exception( ERROR, "Model::RandomPerturb",
                                     "Can only perturb scalar values so far" );

    double dmin, dmax;
    ScalarVariable  sc;
    sgroup.MinMaxOf( prop, dmin, dmax );
    
    random_device rd;
    // seed value is designed specifically to make initialization
    // parameters of mt19937 (instance of mersenne_twister_engine<>)
    // different across executions of application
    mt19937::result_type seed = rd() ^ (
            (mt19937::result_type)
            chrono::duration_cast<chrono::seconds>(
                chrono::system_clock::now().time_since_epoch()
                ).count() +
            (mt19937::result_type)
            chrono::duration_cast<chrono::microseconds>(
                chrono::high_resolution_clock::now().time_since_epoch()
                ).count() );

    mt19937 gen(seed);

    // generating floating point values
    uniform_real_distribution<double> rngen(0, by_percent_of_max_value * dmax * 0.01);

    switch( prop_key.place )
      {
         case NODE:
              for ( typename vector<Node<dim>*>::const_iterator
                    nit=sgroup.NodesBegin(); nit!=sgroup.NodesEnd(); nit++ )
                {
                   (*nit)->Read( prop_key, sc );
                   sc -= rngen(gen);
                   (*nit)->Store( prop_key, sc );
                }
           break;
         case ELEMENT_INTEGRATION_POINT:
              for ( typename vector<Element<dim>*>::const_iterator
                    eit=sgroup.CellsBegin(); eit!=sgroup.CellsEnd(); eit++ )
                for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                {
                   (*eit)->Read( i, prop_key, sc );
                   sc -= rngen(gen);
                   (*eit)->Store( i, prop_key, sc );
                }
           break;
         case ELEMENT:
              for ( typename vector<Element<dim>*>::const_iterator
                    eit=sgroup.CellsBegin(); eit!=sgroup.CellsEnd(); eit++ )
                {
                   (*eit)->Read( prop_key, sc );
                   sc -= rngen(gen);
                   (*eit)->Store( prop_key, sc );
                }
           break;
         default:
           cout <<"\nrandomPerturb: property placement not handled."<< endl;
      }
      
 } // end RandomPerturb

template void randomPerturb( Model<1U>&, const char*, double );
template void randomPerturb( Model<2U>&, const char*, double );
template void randomPerturb( Model<3U>&, const char*, double );






/**
    convert the flag(s) of a variable into integer values stored in its number part
    
    @attention works only for scalars and basic property placements.
    
    @author SKM 21/5/2014
*/
template<uint32_t dim>
void flagToNumber( Model<dim>& model, const char* variable )
 {
    csmp::Region<dim>&  mref(model.Region("Model"));
    csmp::Index  prop_key = model.Database().StorageKey(variable);
   
    if ( prop_key.type != SCALAR )
      throw csmp::Exception( ERROR, "flagToNumber:", "method has not been implemented yet" );

    switch( prop_key.place )
      {
         case NODE:
              for ( typename vector<Node<dim>*>::const_iterator
                    nit=mref.NodesBegin(); nit!=mref.NodesEnd(); nit++ )
                {
                   // overwrites variable value with integer value of its flag enum
                   double value = static_cast<double>( (*nit)->Status(prop_key) );
                   (*nit)->Store( prop_key, makeScalar( (*nit)->Status(prop_key), value ) );
                }
           break;
         case ELEMENT_INTEGRATION_POINT:
              for ( typename vector<Element<dim>*>::const_iterator
                    eit=mref.CellsBegin(); eit!=mref.CellsEnd(); eit++ )
                for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ )
                {
                   double value = static_cast<double>( (*eit)->Status(prop_key) );
                   (*eit)->Store( prop_key, makeScalar( (*eit)->Status(prop_key), value ) );
                }
           break;
         case ELEMENT:
              for ( typename vector<Element<dim>*>::const_iterator
                    eit=mref.CellsBegin(); eit!=mref.CellsEnd(); eit++ )
                {
                   double value = static_cast<double>( (*eit)->Status(prop_key) );
                   (*eit)->Store( prop_key, makeScalar( (*eit)->Status(prop_key), value ) );
                }
           break;
         default:
           cout <<"\nflagToNumber: property placement not handled."<< endl;
      }
   
 } // end flagToNumber

template void flagToNumber( Model<1U>&, const char* );
template void flagToNumber( Model<2U>&, const char* );
template void flagToNumber( Model<3U>&, const char* );




/**
    convert the flag(s) of first variable into double values stored in the second variable
    
    @attention works only for node-property placement.
    
    @author SKM 8/12/2016
*/
template<uint32_t dim>
void flagToNumber( Model<dim>& model, const char* flag_variable, const char* value_variable )
 {
    csmp::Region<dim>&  mref(model.Region("Model"));
    csmp::Index  flag_key = model.Database().StorageKey(flag_variable);  // input
    csmp::Index  prop_key = model.Database().StorageKey(value_variable); // output
 
    if ( flag_key.type != prop_key.type )
      throw csmp::Exception( ERROR, "flagToNumber:", "flag and value variables must be of the same type." );

    if ( flag_key.place != prop_key.place )
      throw csmp::Exception( ERROR, "flagToNumber:", "flag and value variables must have the same placement." );

    if ( flag_key.type != SCALAR and flag_key.type != VECTOR )
      throw csmp::Exception( ERROR, "flagToNumber:", "method handles only scalar and vector variables." );
  
    switch( prop_key.place )
      {
         case NODE:
              if ( flag_key.type == SCALAR ) {
                  for ( typename vector<Node<dim>*>::const_iterator
                        nit=mref.NodesBegin(); nit!=mref.NodesEnd(); nit++ )
                    {
                       // retrieves status of the flag variable
                       double value = static_cast<double>( (*nit)->Status(flag_key) );
                       // overwrites value of value variable with integer value of its flag enum
                       (*nit)->Store( prop_key, makeScalar( (*nit)->Status(prop_key), value ) );
                    }
                }
              else if ( flag_key.type == VECTOR ) {
                  VectorVariable<dim> vc;
                  for ( typename vector<Node<dim>*>::const_iterator
                        nit=mref.NodesBegin(); nit!=mref.NodesEnd(); nit++ )
                    {
                       (*nit)->Read( flag_key, vc );
                       for ( auto i{0U}; i<dim; ++ i )
                         vc(i) = static_cast<double>( vc.Flag(i) );
                       (*nit)->Store( prop_key, vc );
                    }
                }
           break;
         default:
           throw csmp::Exception( ERROR, "flagToNumber:", "property placement not handled." );
      }
   
 } // end flagToNumber

template void flagToNumber( Model<1U>&, const char*, const char* );
template void flagToNumber( Model<2U>&, const char*, const char* );
template void flagToNumber( Model<3U>&, const char*, const char* );












/**

StripDomainEdgesFor() tests the spatial distribution of a scalar input element
property for outliers and removes these. Property outliers are elements
that constitute a property value boundary with >=two of their faces.
When such elements are detected, their property value is set to the
average of the surrounding elements. The property is changed only if
the outside property is either smaller or greater than the element
property.

@section arguments Input Arguments

The name of the property whose variations over the mesh shall be defined
by relatively smooth boundaries.

@section implementation Implementation

The method using the connections among elements to test whether the
element represents a property outlier.

@section application Application

The method is used for regular meshes which were created from pixel-type
input data. In this case the method allows to capitalize on the element
splits which were created by the Triangulator meshing tool along
property boundaries. The result are boundaries with 45o segments that
superseed the stepwise property boundaries of the original mesh.

Numerous calls to StripDomainEdgesFor() also allow to erode regions
defined by stepwise property variations.

@section messages Messages

The method will always warn the user that the model properties are
modified. If the target property is not an element property, the
simulation will be halted by a fatal error.

If the target property is not a scalar variable the method will return
without modifying the target property and it will report a warning.

 */
void stripDomainEdgesFor( Model<2U>& sg, const char* el_prop )
 {
     const csmp::Index  prop_key = sg.Database().StorageKey(el_prop);

     if ( prop_key.place != ELEMENT )
       throw csmp::Exception( FATAL_ERROR, "stripDomainEdgesFor<2U>::StripDomainEdgesFor",
                                    "The requested property is not an element variable");

     if ( prop_key.type != SCALAR ) {
          throw csmp::Exception( WARNING, "stripDomainEdgesFor<doubleoat,2U>::StripDomainEdgesFor",
                                   "only SCALAR variables are handled so far");
          return;
       }
     
    map<size_t,ScalarVariable > new_sc_data;
    ScalarVariable              sc;

    csmp::Region<2>&  model_domain(sg.Region("Model"));

    for ( auto n=0U; n<model_domain.Cells(); n++ )
       {
          //  for elements that are not located at model boundary
          if ( atBoundary( model_domain.E(n) ) == NOT )
            {
               // getting the scalar variable data
               model_domain.E(n)->Read( prop_key, sc );
              
               // checking whether element-property should be changed
               // because the element is located at a region boundary
               // ---------------------------------------------------
               // 1. counting the surrounding values that are different from el-value
               double     sc_sum(0U);
               unsigned int counter(0U);
               for ( auto i{0U}; i<model_domain.E(n)->Neighbors(); i++ ) {
                   assert( model_domain.E(n)->Neighbor(i) != nullptr );
                   if ( sc() > model_domain.E(n)->Neighbor(i)->Read( prop_key ) ) {
                        sc_sum += model_domain.E(n)->Neighbor(i)->Read( prop_key );
                        counter++;
                     }
                 }
               // if more than 2 neighbors have a different property value, this value
               // is assigned to the element
               // TODO: if were are not dealing with triangular elements, this number (2U) is not correct
               if ( counter >= 2U ) sc = sc_sum / static_cast<double>(counter);
          
               // storing the new values of only those elements that must be changed
               new_sc_data[ n ] = sc;
            }
       }
    
     // modyfying those elements that were found to be isolated
     // this implies that isolated squares are removed
     for ( map<size_t,ScalarVariable >::iterator
           sc_it=new_sc_data.begin(); sc_it!=new_sc_data.end(); sc_it++ )
       model_domain.E( (*sc_it).first )->Store( prop_key, (*sc_it).second );
       
     cout <<"\n\nstripDomainEdgesFor<2U>::StripDomainEdgesFor: "<< new_sc_data.size() <<" '"<< el_prop;
     cout <<"' domain-edge elements have been modified to create a smoother boundary."<< endl;
            
   } // end StripRoughDomainEdgesFor











/**
     Assigns chosen node coordinate (x or y or z) to the target node variable.
*/
template<uint32_t dim>
void assignNodeCoordinatesTo( Model<dim>& sg, const char coordinate, const char* node_var )
 {
      csmp::Index nvar_key = sg. Database().StorageKey(node_var);
      assert( nvar_key.place == NODE );
      
      Region<dim>&  sgref(sg.Region("Model"));

      const typename vector<Node<dim>* >::const_iterator  nodesEnd(sgref.NodesEnd());
      
      if ( coordinate == 'x' or coordinate == 'X' )
        for ( typename vector<Node<dim>* >::const_iterator nit=sgref.NodesBegin(); nit!=nodesEnd; ++nit )
          (*nit)->Store( nvar_key, makeScalar( (*nit)->Status(nvar_key), (*nit)->x() ) );
        
      if ( dim > 1 and (coordinate == 'y' or coordinate == 'Y') )
        for ( typename vector<Node<dim>* >::const_iterator nit=sgref.NodesBegin(); nit!=nodesEnd; ++nit )
          (*nit)->Store( nvar_key, makeScalar( (*nit)->Status(nvar_key), (*nit)->y() ) );

      if ( dim > 2 and (coordinate == 'z' or coordinate == 'Z') )
        for ( typename vector<Node<dim>* >::const_iterator nit=sgref.NodesBegin(); nit!=nodesEnd; ++nit )
          (*nit)->Store( nvar_key, makeScalar( (*nit)->Status(nvar_key ), (*nit)->z() ) );
 
 }  // end

template void assignNodeCoordinatesTo( Model<1U>&, const char, const char* );
template void assignNodeCoordinatesTo( Model<2U>&, const char, const char* );
template void assignNodeCoordinatesTo( Model<3U>&, const char, const char* );




/**

Function evaluates that the VSet connectivity is exactly the same
as the data in the current Model!

The Model is used as the reference case.

*/
template<uint32_t dim>
bool compareConnectivity( const Model<dim>& sg, const VSet<dim>& vset )
 {
    bool correct(true);
   
    const Region<dim>&  gref(sg.Region("Model"));
    if ( gref.Cells() != vset.Elements() ) cout <<"\ncompareConnectivity: element number mismatch."<< endl;
    if ( gref.Nodes() != vset.Vertices() ) cout <<"\ncompareConnectivity: node number mismatch."<< endl;
  
    // 1. plist
    for ( size_t i=0U; i<gref.Cells(); i++ )
      {
         for ( uint32_t j{0U}; j<gref.E(i)->Nodes(); j++ )
           if ( gref.E(i)->N(j)->Idx() != vset.Plist( gref.E(i)->Idx(), j ) ) {
                 cerr <<"\ncompareConnectivity: plist inconsistency: sg node id: "<< gref.E(i)->N(j)->Idx();
                 cerr <<" vs. vset nid: "<< vset.Plist( gref.E(i)->Idx(), j );
                 correct = false;
             }
      }
    
    // 2. pfverts
    for ( size_t i=0U; i<gref.Cells(); i++ )
      {
         for ( uint32_t j{0U}; j<gref.E(i)->Neighbors(); j++ )
           if ( gref.E(i)->Neighbor(j) and
                static_cast<int32_t>(gref.E(i)->Neighbor(j)->Idx()) != vset.Pfvert( gref.E(i)->Idx(), j ) ) {
                 cerr <<"\ncompareConnectivity: plist inconsistency: sg node id: "<< gref.E(i)->Neighbor(j)->Idx();
                 cerr <<" vs. vset nid: "<< vset.Pfvert( gref.E(i)->Idx(), j );
                 correct = false;
             }
      }
   
   return correct;
    
 } // end compare

template bool compareConnectivity<2U>( const Model<2U>&, const VSet<2U>& );
template bool compareConnectivity<3U>( const Model<3U>&, const VSet<3U>& );








/**
    Imposes a user-defined upper or lower limit on the value of the variable of interest.
    
    For a vector variable, its length gets scaled to the limit value.
    For a tensor variable nothing can be done yet, so an exception is thrown.
    
    @author SKM 7/9/2014
*/
template<uint32_t dim>
void imposeLimitOn( Model<dim>& model, const char* region, const char* variable, bool upper_limit, double limit_value )
 {
    Region<dim>&  rref(model.Region(region));
    csmp::Index   prop_key(model.Database().StorageKey(variable));
    double      min, max;
    model.Database().RangeOf( variable, min, max );
   
    if ( upper_limit && limit_value > max ) {
         cerr <<"\nIntended upper limit on variable '"<< variable <<"' exceeds that defined in database: ";
         cerr << limit_value <<" vs. "<< max << endl;
         throw csmp::Exception( ERROR, "imposeLimitOn:", "user-defined limit is out of bounds specified in variable database." );
      }
    if ( !upper_limit && limit_value < min ) {
         cerr <<"\nIntended lower limit on variable '"<< variable <<"' is lower than that defined in database: ";
         cerr << limit_value <<" vs. "<< min << endl;
         throw csmp::Exception( ERROR, "imposeLimitOn:", "user-defined limit is out of bounds specified in variable database." );
      }
   
    if ( prop_key.type == SCALAR ) {
        ScalarVariable  sc;
        if ( upper_limit )
          switch( prop_key.place )
            {
               case MODEL:
                  model.Store( prop_key, makeScalar( model.Status(prop_key), std::min(limit_value,rref.Read(prop_key)) ) );
                 break;
               case REGION:
                  rref.Store( prop_key, makeScalar( rref.Status(prop_key), std::min(limit_value,rref.Read(prop_key)) ) );
                 break;
               case ELEMENT:
                  for ( auto it=rref.CellsBegin();  it!=rref.CellsEnd(); ++it ) {
                      double val = (*it)->Read( prop_key );
                      (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key), std::min(limit_value,val) ) );
                   }
                 break;
               case ELEMENT_INTEGRATION_POINT:
                  for ( auto it=rref.CellsBegin();  it!=rref.CellsEnd(); ++it )
                    for ( auto i{0U}; i<(*it)->IntegrationPoints(); i++ ) {
                         double val = (*it)->Read( i, prop_key );
                         (*it)->Store( i, prop_key, makeScalar( (*it)->Status(i,prop_key), std::min(limit_value,val) ) );
                      }
                 break;
               case NODE:
                  for ( auto it=rref.NodesBegin(); it!=rref.NodesEnd(); ++it ) {
                      double val = (*it)->Read( prop_key );
                      (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key), std::min(limit_value,val) ) );
                   }
                 break;
               default:
                 throw csmp::Exception( ERROR, "imposeLimitOn:", "variable placement not recognized." );
            }
          else // if a lower limit shall be imposed
          switch( prop_key.place )
            {
               case MODEL:
                  model.Store( prop_key, makeScalar( model.Status(prop_key), std::max(limit_value,rref.Read(prop_key)) ) );
                 break;
               case REGION:
                  rref.Store( prop_key, makeScalar( rref.Status(prop_key), std::max(limit_value,rref.Read(prop_key)) ) );
                 break;
               case ELEMENT:
                  for ( auto it=rref.CellsBegin();  it!=rref.CellsEnd(); ++it ) {
                      double val = (*it)->Read( prop_key );
                      (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key), std::max(limit_value,val) ) );
                   }
                 break;
               case ELEMENT_INTEGRATION_POINT:
                  for ( auto it=rref.CellsBegin();  it!=rref.CellsEnd(); ++it )
                    for ( auto i{0U}; i<(*it)->IntegrationPoints(); i++ ) {
                         double val = (*it)->Read( i, prop_key );
                         (*it)->Store( i, prop_key, makeScalar( (*it)->Status(i,prop_key), std::max(limit_value,val) ) );
                      }
                 break;
               case NODE:
                  for ( auto it=rref.NodesBegin(); it!=rref.NodesEnd(); ++it ) {
                      double val = (*it)->Read( prop_key );
                      (*it)->Store( prop_key, makeScalar( (*it)->Status(prop_key), std::max(limit_value,val) ) );
                   }
                 break;
               default:
                 throw csmp::Exception( ERROR, "imposeLimitOn:", "variable placement not recognized." );
            }
      }
   
    // for a vector variable, its length gets scaled to the limit value
    else if ( prop_key.type == VECTOR ) {
        VectorVariable<dim>  vc;
        switch( prop_key.place )
          {
             case MODEL: {
                    model.Read( prop_key, vc );
                    const double vmagnitude = vc.Length();
                    assert( vmagnitude > 0. );
                    // if the vector is too long it gets scaled back
                    if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                    else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                    model.Store( prop_key, vc );
                 }
               break;
             case REGION: {
                    rref.Read( prop_key, vc );
                    const double vmagnitude = vc.Length();
                    assert( vmagnitude > 0. );
                    if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                    else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                    rref.Store( prop_key, vc );
                 }
               break;
             case ELEMENT:
                for ( auto it=rref.CellsBegin();  it!=rref.CellsEnd(); ++it ) {
                    (*it)->Read( prop_key, vc );
                    const double vmagnitude = vc.Length();
                    assert( vmagnitude > 0. );
                    if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                    else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                    (*it)->Store( prop_key, vc );
                 }
               break;
             case ELEMENT_INTEGRATION_POINT:
                for ( auto it=rref.CellsBegin();  it!=rref.CellsEnd(); ++it )
                  for ( auto i{0U}; i<(*it)->IntegrationPoints(); i++ ) {
                       (*it)->Read( i, prop_key, vc );
                        const double vmagnitude = vc.Length();
                        assert( vmagnitude > 0. );
                        if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                        else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                       (*it)->Store( i, prop_key, vc );
                    }
               break;
             case NODE:
                for ( auto it=rref.NodesBegin(); it!=rref.NodesEnd(); ++it ) {
                    (*it)->Read( prop_key, vc );
                    const double vmagnitude = vc.Length();
                    assert( vmagnitude > 0. );
                    if ( upper_limit and vmagnitude > max ) vc /= (vmagnitude / max);
                    else if ( vmagnitude < min ) vc *= (min / vmagnitude);
                    (*it)->Store( prop_key, vc );
                 }
               break;
             default:
               throw csmp::Exception( ERROR, "imposeLimitOn:", "variable placement not recognized." );
          }
      }
    else throw csmp::Exception( ERROR, "imposeLimitOn:", "variable type not recognized." );

 } // end imposeLimitOn

template void imposeLimitOn( Model<1U>&, const char*, const char*, bool, double );
template void imposeLimitOn( Model<2U>&, const char*, const char*, bool, double );
template void imposeLimitOn( Model<3U>&, const char*, const char*, bool, double );









  // explicit instantiations
template class RegionInterface<1U, Model>;
template class RegionInterface<2U, Model>;
template class RegionInterface<3U, Model>;

template class BoundaryInterface<1U, Model>;
template class BoundaryInterface<2U, Model>;
template class BoundaryInterface<3U, Model>;

template class SplitBoundaryInterface<1U, Model>;
template class SplitBoundaryInterface<2U, Model>;
template class SplitBoundaryInterface<3U, Model>;

template class Model<1U>;
template class Model<2U>;
template class Model<3U>;

} // end namespace csmp
