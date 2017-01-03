#include "Model.h"
#include "VSet.h"
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
#include "CSMP_highLevelUtilities.h"
#include "binaryReadWrite.h"
#include "ModelTime.h"
#include "FiniteVolumeStencilManager.h"

#include <cstring>

using namespace std;

namespace csmp {

// PROTECTED CONSTRUCTORS

/**
    Constructs a completely empty model.
    It is used by one of the ANSYS Model interfaces.
 
    @note model will contain a single region called 'model' and it will be considered unique.
    
    @todo SKM deprecate; make the default constructor private to comply with inheritance principles
*/
template<size_t dim>
Model<dim>::Model()
  : database_(),
    fvStencilManager_(NULL),
    model_name_("undefined"),
    verbose_(true)
  {
  }


/**
    Constructor used for base class construction in subclasses of the model.
    Example: ANSYS_Model3D where model construction input data are read from
    file, initializing VSet and ModelTopology classes that serve as an input
    for model construction.
    
    @note this is not a stand-alone constructor and leaves the model in an incomplete state.
*/
template<size_t dim>
Model<dim>::Model( const std::string& varFile, bool binary )
  : database_( varFile.c_str(), binary ),
    fvStencilManager_(NULL),
    model_name_("undefined"),
    verbose_(true)
{
   InitializeLocalVariableStorage();
}


// PUBLIC CONSTRUCTORS


/**
    Full input from CSMP-native binary file.
*/
template<size_t dim>
Model<dim>::Model( const std::string& binaryFileNames )
  : database_( BinaryVariablesFileName(binaryFileNames.c_str()).c_str(), true ),
    model_name_(binaryFileNames),
    fvStencilManager_(NULL)
{
  InitializeLocalVariableStorage();
//  InputFromDisk(binaryFileNames); // WORKS, but all regions are rebuilt from scratch
  InputFromBinaryFile(binaryFileNames.c_str());
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
  (default: CSMP_variables.txt).

  @section application Application

  Constructor is used when an ANSYS Model is built from topology and VData.
*/
template<size_t dim>
Model<dim>::Model( VSet<dim>& vset, const char* var_file, bool isoparametric_elements, bool binaryVariablesFile )
  : database_( var_file, binaryVariablesFile ),
    model_name_("undefined"),
    fvStencilManager_(NULL),
    verbose_(true)
 {
    Initialize( isoparametric_elements, vset,
                false /* do not create boundaries */,
                false /* regular boundaries */);

 } // end VSet constructor


template<size_t dim>
Model<dim>::Model( VSet<dim>& vset, bool isoparametric_elements )
  : database_(),
    model_name_("undefined"),
    fvStencilManager_(NULL)
 {
    Initialize( isoparametric_elements, vset,
                false /* do not create boundaries */,
                false /* regular boundaries */);

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
  
  @attention this constructor will not create any boundaries
*/
template<size_t dim>
Model<dim>::Model( ModelTopology& mesh_topology, VSet<dim>& vset, const char* var_file, bool binaryVariablesFile )
 : database_( var_file, binaryVariablesFile ),
   fvStencilManager_(NULL),
   model_name_("undefined"),
   verbose_(true)
 {
    Initialize( mesh_topology, vset,
                false /* do not create boundries */,
                false /* regular boundaries */ );

 } // end VSet/ModelTopology constructor



template<size_t dim>
Model<dim>::Model( ModelTopology& mesh_topology, VSet<dim>& vset )
 : database_( ),
   model_name_("undefined"),
   fvStencilManager_(NULL)
 {
    Initialize( mesh_topology, vset,
                false /* do not create boundaries */,
                false /* regular boundaries */);
   
 } // end VSet/ModelTopology constructor









/**  Initialise( topology, vset, create boundary ...)

   Performs the following steps:
   
   0. eliminates unwanted mesh regions from topology and corresponding elements from vset
   
   1. Reduces mesh specifications to actual desired element types as specified by the topology object
   
   2. initializing the finite-element manager
   
   3. builds finite element mesh and property storage -> done by MeshManager
   
   4. forms unique root Region called Model (is in unique regions if there are no other unique regions)
   
   5. Tests with a flood-fill whether the model is contiguous
 
   6. Assigns node and element material properties and variable values to model -> done by InputVariablesFrom(vset)
   
   7. Associates supplied subregions with regions (model subdomains) -> done by FormRegionsFrom(topology)
   
   8. Forming Boundaries -> done by EstablishBoundaries()
   
   9. Adding potentially required property storage for regions and boundaries (however these properties are not initialised here

   @attention MOST COMMONLY USED MODEL CONSTRUCTION  FROM EXTERNAL DATA METHOD - including ANSYS_Model3D

*/
template<size_t dim>
void Model<dim>::Initialize( const char* regions_file_prefix,
                             ModelTopology& mesh_topology,
                             VSet<dim>& vset,
                             bool create_boundaries,
                             bool non_box_shaped_model )
{
    // 1. eliminating the unwanted mesh regions from topology and vset
    mesh_topology.ReduceToRegions( regions_file_prefix );
  
    // 2. building the model
    Initialize( mesh_topology, vset, create_boundaries, non_box_shaped_model );
  
} // end Initialize (with regions from file)





/**
    custom constructor
*/
template<size_t dim>
void Model<dim>::Initialize( ModelTopology& mesh_topology,
                             VSet<dim>& vset,
                             bool create_boundaries,
                             bool non_box_shaped_model )
{
    // 1. Reducing the mesh data to the desired element types as specified
    //    by the topology object
    map<size_t,size_t>  old_and_new_elmtids;
    mesh_topology.CreateNewElementNumbers( old_and_new_elmtids );
    vset.ReduceTo( old_and_new_elmtids );
    old_and_new_elmtids.clear();

    // 2. initializing the finite-element manager true=isoparametric
    fem_manager_.InitializeElements( dim,
                                     mesh_topology.InterpolationOrder(),
                                     mesh_topology.IsoparametricElements() );

    // 3. building the finite element mesh and property storage
    mesh_manager_.Initialize( Database(), FE_Manager(), vset );

    // 4. forming unique root Region called "All Elements" as well as default computational domain called "Model"
    const bool withNeighborConnectivity(true);
    const bool valid_master_region = this->CreateNonUniqueMasterRegionFromRootNode( withNeighborConnectivity );
    assert( valid_master_region );
    const bool place_copy_in_unique_regions( (mesh_topology.ModelRegions()==0) );
    this->CopyRegion( "All Elements", "Model", place_copy_in_unique_regions );

    // 5. Testing with a flood-fill whether the model is contiguous
    //    if not Accumulate all will not have reached all the elements
    if ( Mesh().Elements() < vset.Elements() )
        throw csmp::Exception( FATAL_ERROR, "Model<dim>::Initialize",
                                            "Model appears to be fragmented. Are all regions connected?" );
    cout <<"\nModel<dim>::Initialize: ";
    cout <<"Mesh has been built successfully..." << endl;

    // 6. Associating supplied subregions with regions (model subdomains)
    this->FormRegionsFrom( mesh_topology );

    if ( mesh_topology.BoxShapedModel() ) {
         cout <<"\nModel<"<< dim <<">::Initialize: This model is box-shaped so that you can assign ";
         cout <<"boundary conditions in the standard way. "<< endl;
      }

    // 7. Forming Boundaries
    if ( create_boundaries )
      {
          bool box_shaped(this->BoxShaped());
          if ( box_shaped and box_shaped != !non_box_shaped_model )
            ErrorHandler::Instance().notice( WARNING, "Model<dim>::Initialize:",
                                            "while model contains all relevant box side boundaries it will be treated as non-box shaped." );
          
          if ( !non_box_shaped_model && box_shaped )
            {
                this->EstablishBoxBoundaries();
                // build boundary edges consisting of faces
                if ( dim == 3U ) {
                     const bool created_edges(this->EstablishEdgeBoundariesOfBoxShapedModel());
                     assert( created_edges );
                  }
            }
          else this->EstablishBoundaries();
      }
    else cout<<"\nModel<dim>::Initialize: CSMP boundaries disabled." << endl;

    // 8. Adding property storage to the Model
    InitializeLocalVariableStorage();
    UpdateSubdomainPropertyStorage();

    // 9. assigning properties to mesh
    InputVariablesFrom( vset );

    cout << "\n============================================================================";
    cout << "\nModel '"<< Name() <<"' has been established successfully!";
    cout << "\n============================================================================";
    cout << endl;
  
} // end Initialize (VSet / ModelTopology)




// TESTING
/*

set<BOX_BOUNDARY> node_flags;
const Region<dim>& model(this->Region("Model"));
for ( auto nit=model.NodesBegin(); nit!=model.NodesEnd(); nit++ )
      node_flags.insert( (*nit)->AtBoundary() );
*/
              

/**
    Initialises model from VSet. Very similar to Initialise(VSet,ModelTopology), but without
    the creation of regions other than 'Model'.
    
    @note this version of Initialise() is called when the model is build from ANSYS input data
*/
template<size_t dim>
void Model<dim>::Initialize( bool isoparametric_elements,
                             VSet<dim>& vset,
                             bool create_boundaries,
                             bool non_box_shaped_model )
{
    // 0. initializing the finite-element manager true=isoparametric
    fem_manager_.InitializeElements( dim,
                                     vset.OrderOfFiniteElementInterpolationFunctions(),
                                     isoparametric_elements );

    // 1. building the finite element mesh and property storage
    mesh_manager_.Initialize( Database(), FE_Manager(), vset );

    // 2. forming unique root Region called "All Elements"
    const bool withNeighborConnectivity(true);
    this->CreateNonUniqueMasterRegionFromRootNode( withNeighborConnectivity );
    const bool place_copy_in_unique_regions(true);
    this->CopyRegion( "All Elements", "Model", place_copy_in_unique_regions );

    // 3. Testing with a flood-fill whether the model is contiguous
    //    if not Accumulate all will not have reached all the elements
    if ( Mesh().Elements() < vset.Elements() )
        throw csmp::Exception( FATAL_ERROR, "Model<dim>::Initialize",
                                            "Model appears to be fragmented. Are all regions connected?" );

    cout <<"\nModel<dim>::Initialize(VSet): ";
    cout <<"Mesh has been built successfully..." << endl;

    // 4. Forming Boundaries
    if ( create_boundaries ) {
          if ( !non_box_shaped_model && this->BoxShaped() ) {
                this->EstablishBoxBoundaries();
                if ( dim == 3U ) this->EstablishEdgeBoundariesOfBoxShapedModel();
            }
          else this->EstablishBoundaries();
      }
    else cout<<"\nModel<dim>::Initialize: CSMP boundaries disabled." << endl;

    // 5. Adding potentially required property storage
    InitializeLocalVariableStorage();
    UpdateSubdomainPropertyStorage();

    // 6. assigning properties to mesh
    InputVariablesFrom( vset );

    cout << "\n================================================";
    cout << "\nModel has been established successfully!";
    cout << "\n================================================";
    cout << endl;
}




template<size_t dim>
const char* Model<dim>::Name() const
 {
    return model_name_.c_str();
 }


template<size_t dim>
void Model<dim>::Name( const char* new_name )
 {
    model_name_ = new_name;
 }




/**
     Instantiates FiniteVolumeStencilManager and initialises FiniteVolumeStencils 
     if correpoding pointers are NULL.
     
     delegates this step to MeshManager::InitializeFiniteVolumeStencils()
 
*/
template<size_t dim>
void Model<dim>::InstantiateFiniteVolumes()
{
  if(!fvStencilManager_)
    {
      fvStencilManager_ = new FiniteVolumeStencilManager<dim>();
      cout << "\nModel<dim>::InstantiateFiniteVolumeStencilManager: Created local FiniteVolumeStencilManager\n";
    }
  Mesh().InitializeFiniteVolumeStencils( Database(), fem_manager_, *fvStencilManager_ );
}


/// removes dynamically allocated finite volume stencil manager
template<size_t dim>
Model<dim>::~Model()
  {
    if(fvStencilManager_)
      delete fvStencilManager_;
  } // end



template<size_t dim>
void Model<dim>::Verbose( bool verbose )
{
    this->verbose_ = verbose;
}



template<size_t dim>
bool Model<dim>::Verbose()
{
    return this->verbose_;
}

template<size_t dim>
string Model<dim>::BinaryVsetFileName( const char* base_file_name )
{
  string fullFileName(base_file_name);
  fullFileName.append(".vset");
  return fullFileName;
}

template<size_t dim>
string Model<dim>::BinaryRegionsFileName( const char* base_file_name )
{
  string fullFileName(base_file_name);
  fullFileName.append("_regions.dat");
  return fullFileName;
}

template<size_t dim>
string Model<dim>::BinaryBoundariesFileName( const char* base_file_name )
{
  string fullFileName(base_file_name);
  fullFileName.append("_boundaries.dat");
  return fullFileName;
}

template<size_t dim>
string Model<dim>::BinarySplitBoundariesFileName( const char* base_file_name )
{
  string fullFileName(base_file_name);
  fullFileName.append("_splitboundaries.dat");
  return fullFileName;
}

template<size_t dim>
string Model<dim>::BinaryVariablesFileName( const char* base_file_name )
{
  string fullFileName(base_file_name);
  fullFileName.append("_variables.dat");
  return fullFileName;
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
template<size_t dim>
template<class Var>
void  Model<dim>::OutputVariableTo( const char* out_var, FEM_Data<Var>& data ) const
 {
    this->Region(this->MasterRegion().c_str()).OutputVariableTo( out_var, data );
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
template<size_t dim>
void Model<dim>::InputVariablesFrom( const VSet<dim>& vset )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if( !vset.DataEmpty() ) mesh_manager_.InputStoredVariablesFrom( Database(), vset );
     else ErrorHandler::Instance().notice( INFO, "Model<dim>::InputVariablesFrom:", "No properties found in VSet." );

     // properties and values stored on the model itself
     for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
       {
         // apart from the name string key in the map, PropertyData contains the most important variable specifications
         if ( (*pit).second.Placement() != MODEL ) continue;
         // some checks
         assert( Database().IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(Database().StorageKey( (*pit).first.c_str() ));
         assert( key.place == MODEL );
         assert( (*pit).second.Size() / key.dataDepth == 1U );
        
         switch( key.type )
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
               csmp_error.notice( ERROR, "Model<dim>::InputVariablesFrom:",
                                  (*pit).first, "type of Model variable not recognized.");
           }
      }
}




/**
  InputVariableFrom() lets you input variable data stored in a FEM_Data
  template class object to a Model variable. class T here is a place
  holder for the data type which may be a double64 or any basic CSP variable.

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
template<size_t dim>
template<class Var>
void  Model<dim>::InputVariableFrom( const char* input_prop, const FEM_Data<Var>& vdata )
 {
    const string master_region("All Elements");
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
    Checks whether finite-volume related properties are contained in the VSet
*/
template<size_t dim>
bool Model<dim>::VSetHasFiniteVolumeProperties( const VSet<dim>& vset ) const
 {
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit ) {
         if ( (*pit).second.Placement() == SECTOR_INTEGRATION_POINT ) return true;
         if ( (*pit).second.Placement() == FACET_INTEGRATION_POINT ) return true;
         if ( (*pit).second.Placement() == FACE_SECTOR_INTEGRATION_POINT ) return true;
         if ( (*pit).second.Placement() == FACE_FACET_INTEGRATION_POINT ) return true;
         if ( (*pit).second.Placement() == INTER_FACE_SECTOR_INTEGRATION_POINT ) return true;
         if ( (*pit).second.Placement() == INTER_FACE_FACET_INTEGRATION_POINT ) return true;
      }
   
    return false;
 }



// Renumbers nodes, elements, faces & interfaces. Nodes and Element/Face/Interface may have same values, but the latter may not.
template<size_t dim>
size_t Model<dim>::UpdateIndices() const
  {
    // elements and nodes
    this->Region( "Model" ).UpdateMemberIndexes();
    size_t runningIndex( this->Region( "Model" ).Elements() );

    // boundaries
    for( typename Model<dim>::splitBoundaryConstIterator bit( this->SplitBoundariesBegin() ); bit != this->SplitBoundariesEnd(); ++bit )
      {
        for( typename vector<InterFace<dim>*>::const_iterator face( bit->second.ElementsBegin() ); face !=  bit->second.ElementsEnd(); ++face )
          (*face)->Idx( runningIndex++ );
      }
    for( typename Model<dim>::boundaryConstIterator bit( this->BoundariesBegin() ); bit != this->BoundariesEnd(); ++bit )
      {
        for( typename vector<Face<dim>*>::const_iterator face( bit->second.ElementsBegin() ); face !=  bit->second.ElementsEnd(); ++face )
          (*face)->Idx( runningIndex++ );
      }
    return runningIndex;
  }



/**
  The copy-constructor permits to duplicate models which may come in handy
  if one wants to compare the results of slightly different computations
  at runtime. This has, for instance, the advantage that large datasets
  must not be written to file first before they can be compared.
*/
template<size_t dim>
Model<dim>& Model<dim>::operator=( const Model<dim>& model )
 {
    if ( &model != this )
      {
        database_                           = model.database_;
        fem_manager_                        = model.fem_manager_;
        mesh_manager_                       = model.mesh_manager_;
        this->uniqueGroupMap_               = model.uniqueGroupMap_;
        this->groupMap_                     = model.groupMap_;
        this->faceBoundaryMap_              = model.faceBoundaryMap_;
        this->interFaceSplitBoundaryMap_    = model.interFaceSplitBoundaryMap_;
        this->LVS( model.LVS() );
      }
    return *this;

 } // end assignment operator




 
/** Add a property at runtime
@warning Any runtime change in variables invalidates existing Index objects! They need to be refreshed.
@todo (1-F) Support for BOUNDARY... / SPLIT_BOUNDARY (create overload taking index)
*/
template<size_t dim>
csmp::Index  Model<dim>::CreateProperty( const char* new_prop,
                                         const char* unit,
                                         VARIABLE_TYPE vtype,
                                         PLACEMENT vplace,
                                         size_t vsize,
                                         double64 vmin,
                                         double64 vmax,
                                         string usage)
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( database_.IsDefined(new_prop) ) {
         csmp_error.notice( WARNING, "Model<dim>::CreateProperty:", new_prop, "property already exists." );
         return database_.StorageKey(new_prop);
      }

    size_t  prop_index = Database().VariableCount(vplace, vtype);
    csmp::Index  prop_key = database_.AddProperty( new_prop, unit, prop_index, vtype, vplace, vsize , vmin, vmax, usage );

    if ( vplace == NODE ) {
          csmp::Region<dim>&  gref(this->Region("Model"));
          for ( typename vector<csmp::Node<dim>*>::iterator
                eit=gref.NodesBegin(); eit!=gref.NodesEnd(); eit++ )
            (*eit)->AddProperty( prop_key );
      }
    else if ( vplace == ELEMENT || vplace == ELEMENT_INTEGRATION_POINT || vplace == SECTOR_INTEGRATION_POINT || vplace == FACET_INTEGRATION_POINT) {
          csmp::Region<dim>&  gref(this->Region("Model"));
          for ( typename vector<csmp::Element<dim>*>::iterator
                eit=gref.ElementsBegin(); eit!=gref.ElementsEnd(); eit++ )
            (*eit)->AddProperty( prop_key );
      }
    else if ( vplace == FACE || vplace == FACE_INTEGRATION_POINT || vplace == FACE_SECTOR_INTEGRATION_POINT || vplace == FACE_FACET_INTEGRATION_POINT ) {
         for ( typename map<std::string,csmp::Boundary<dim> >::iterator
               git = this->BoundariesBegin(); git != this->BoundariesEnd(); ++git )
           for ( typename vector<csmp::Face<dim>*>::iterator
                 eit=(*git).second.ElementsBegin(); eit!=(*git).second.ElementsEnd(); eit++ )
             (*eit)->AddProperty( prop_key );
      }
    else if ( vplace == INTER_FACE || vplace == INTER_FACE_INTEGRATION_POINT || vplace == INTER_FACE_SECTOR_INTEGRATION_POINT || vplace == INTER_FACE_FACET_INTEGRATION_POINT ) {
         for ( typename map<std::string,csmp::SplitBoundary<dim> >::iterator
               git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); ++git )
           for ( typename vector<csmp::InterFace<dim>*>::iterator
                 eit=(*git).second.ElementsBegin(); eit!=(*git).second.ElementsEnd(); eit++ )
             (*eit)->AddProperty( prop_key );
      }
    else if ( vplace == MODEL ) {
        LocalVariableStorage<dim, Model<dim> >::AddProperty( prop_key );
      }
    else if ( vplace == REGION ) {
         // each region has its unique property value even if they are overlapping!
         for ( typename map<string,csmp::Region<dim> >::iterator
               git=this->UniqueRegionsBegin(); git!=this->UniqueRegionsEnd(); git++ )
           (*git).second.AddProperty( prop_key );
         for ( typename map<string,csmp::Region<dim> >::iterator
               git=this->RegionsBegin(); git!=this->RegionsEnd(); git++ )
           (*git).second.AddProperty( prop_key );
      }
    else if ( vplace == BOUNDARY ) {
         for ( typename map<std::string,csmp::Boundary<dim> >::iterator
               git = this->BoundariesBegin(); git != this->BoundariesEnd(); git++ )
           (*git).second.AddProperty( prop_key );
         for ( typename map<std::string,csmp::SplitBoundary<dim> >::iterator
               git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); git++ )
           (*git).second.AddProperty( prop_key );
      }
    else
    throw csmp::Exception( ERROR, "Model<dim>::CreateProperty:",
                           new_prop, "placement not supported yet, nothing was done." );

    if( Verbose() )
        cout <<"\nModel<dim>::CreateProperty: successfully created the new property '"<< new_prop <<"'\n";
    return prop_key;

 } // end CreateProperty





/** Removes a property from the model. 
@warning All indices have to be updated. 
@todo (1-F) Support for new placements
@warning Existing Index objects may be invalidated after invoking this function
*/
template<size_t dim>
void  Model<dim>::DeleteProperty( const char* property )
 {

    if ( !database_.IsDefined(property) ) {
         throw csmp::Exception( INFO,
                                "Model<dim>::CreateProperty",
                                property,
                                "property does not exist" );
         return;
      }

    csmp::Index  prop_key = database_.StorageKey( property );
    csmp::Region<dim>&  gref(this->Region("Model"));

    if ( prop_key.place == NODE ) {
          for ( typename vector<csmp::Node<dim>*>::iterator
                eit=gref.NodesBegin(); eit!=gref.NodesEnd(); eit++ )
            (*eit)->DeleteProperty( prop_key );
          }
    else if ( prop_key.place == ELEMENT  || prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACET_INTEGRATION_POINT ) {
          for ( typename vector<csmp::Element<dim>*>::iterator
                eit=gref.ElementsBegin(); eit!=gref.ElementsEnd(); eit++ )
            (*eit)->DeleteProperty( prop_key );
      }
    else if ( prop_key.place == FACE || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT ) {
         for ( typename map<std::string,csmp::Boundary<dim> >::iterator
               git = this->BoundariesBegin(); git != this->BoundariesEnd(); ++git )
           for ( typename vector<csmp::Face<dim>*>::iterator
                 eit=(*git).second.ElementsBegin(); eit!=(*git).second.ElementsEnd(); eit++ )
             (*eit)->DeleteProperty( prop_key );
      }
    else if ( prop_key.place == INTER_FACE || prop_key.place == INTER_FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT || prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
         for ( typename map<std::string,csmp::SplitBoundary<dim> >::iterator
               git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); ++git )
           for ( typename vector<csmp::InterFace<dim>*>::iterator
                 eit=(*git).second.ElementsBegin(); eit!=(*git).second.ElementsEnd(); eit++ )
             (*eit)->DeleteProperty( prop_key );
      }
    else if ( prop_key.place == MODEL ) {
        LocalVariableStorage<dim, Model<dim> >::DeleteProperty( prop_key );
      }
    else if ( prop_key.place == REGION ) {
         // each region has its unique property value even if they are overlapping!
         for ( typename map<string,csmp::Region<dim> >::iterator
               git=this->UniqueRegionsBegin(); git!=this->UniqueRegionsEnd(); git++ )
           (*git).second.DeleteProperty( prop_key );
         for ( typename map<string,csmp::Region<dim> >::iterator
               git=this->RegionsBegin(); git!=this->RegionsEnd(); git++ )
           (*git).second.DeleteProperty( prop_key );
      }
    else if ( prop_key.place == BOUNDARY ) {
         for ( typename map<std::string,csmp::Boundary<dim> >::iterator
               git = this->BoundariesBegin(); git != this->BoundariesEnd(); git++ )
           (*git).second.DeleteProperty( prop_key );
         for ( typename map<std::string,csmp::SplitBoundary<dim> >::iterator
               git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); git++ )
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
template<size_t dim>
void  Model<dim>::ExtrapolateElementToNodeProperty( const char* eprop, const char* nprop, bool by_distance )
 {
    this->Region("Model").ExtrapolateElementToNodeProperty( eprop, nprop, by_distance );

 } // end ExtrapolateElementToNodeProperty


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
template<size_t dim>
void Model<dim>::ChangePropertyStatus( const char* input_prop, VARIABLE_FLAG new_status )
 {
    this->Region("Model").ChangePropertyStatus( input_prop, new_status, COMPLETE );

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
template<size_t dim>
void Model<dim>::ChangePropertyStatusWhere( const char* var,
                                            double64 vmin, double64 vmax,
                                            VARIABLE_FLAG new_status )
 {
    this->Region("Model").ChangePropertyStatusWhere( var, new_status, vmin, vmax );

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
template<size_t dim>
void Model<dim>::Accept( csmp::Visitor<dim>& v )
 {
    switch( v.ApplicationLevel() ) {
    // means direct access to elements and nodes without being able to fetch
    // data from regions or boundaries
    case MODEL: {
        csmp::Region<dim>&  mref(this->Region("Model"));
        v.Visit(this);
        switch( v.ApplicationTarget() ) {
        case ELEMENT: {
            for ( typename vector<csmp::Element<dim>*>::iterator
                  el_it=mref.ElementsBegin(); el_it!=mref.ElementsEnd(); el_it++ )
                (*el_it)->Accept( v );
        }
            return;
        case NODE: {
            for ( typename vector<csmp::Node<dim>*>::iterator
                  nd_it=mref.NodesBegin(); nd_it!=mref.NodesEnd(); nd_it++ )
                (*nd_it)->Accept( v );
        }
            return;
        case FACE: { // all faces contained in the model live inside of the boundaries
            for ( typename map<std::string,csmp::Boundary<dim> >::iterator
                  bit = this->BoundariesBegin(); bit != this->BoundariesEnd(); bit++ )
                for ( typename vector<Face<dim>*>::iterator
                      it=(*bit).second.ElementsBegin(); it!=(*bit).second.ElementsEnd(); it++ )
                    (*it)->Accept( v );
        }
            return;
        case INTER_FACE:
            for ( typename map<std::string,csmp::SplitBoundary<dim> >::iterator
                  bit = this->SplitBoundariesBegin(); bit != this->SplitBoundariesEnd(); bit++ )
                for ( typename vector<InterFace<dim>*>::iterator
                      it=(*bit).second.ElementsBegin(); it!=(*bit).second.ElementsEnd(); it++ )
                    (*it)->Accept( v );
            return;
            // access is granted to all boundaries of each type
        case SPLIT_BOUNDARY:
            for( typename Model<dim>::splitBoundaryIterator
                 it( this->SplitBoundariesBegin() ); it != this->SplitBoundariesEnd(); ++it )
                (*it).second.Accept(v);
            break;
        case BOUNDARY:
            for( typename Model<dim>::boundaryIterator
                 it( this->BoundariesBegin() ); it != this->BoundariesEnd(); ++it )
                (*it).second.Accept(v);
            break;
        default:
            throw csmp::Exception( WARNING, "Model<dim>::Accept:",
                                   "application target of visitor unresolved; nothing was done.");
        }
    }
        return;
    case REGION: { // processing all unique regions except for the "Model"
        typename std::map<std::string,csmp::Region<dim> >::iterator  rit = this->UniqueRegionsBegin();
        // if there is only the region "Model"
        if ( this->Regions() == 1U and (*rit).first == "Model" ) {
            (*rit).second.Accept( v );
            return;
        }
        else {
            // visitor is applied to all unique regions except for "Model"
            while ( rit!=this->UniqueRegionsEnd() ) {
                if ( (*rit).first != "Model" ) (*rit).second.Accept(v);
                rit++;
            }
        }
    }
        return;
    case SPLIT_BOUNDARY:
        for( typename Model<dim>::splitBoundaryIterator it( this->SplitBoundariesBegin() );
             it != this->SplitBoundariesEnd();
             ++it )
            (*it).second.Accept(v);
        break;
    case BOUNDARY:
        if ( v.ApplicationTarget() == FACE || v.ApplicationTarget() == BOUNDARY )
            for ( typename map<std::string,csmp::Boundary<dim> >::iterator
                  it = this->BoundariesBegin(); it != this->BoundariesEnd(); ++it )
                (*it).second.Accept( v );
        // TODO: check wether this is redundant
        else if ( v.ApplicationTarget() == INTER_FACE || v.ApplicationTarget() == SPLIT_BOUNDARY )
            for ( typename map<std::string,csmp::SplitBoundary<dim> >::iterator
                  it = this->SplitBoundariesBegin(); it != this->SplitBoundariesEnd(); ++it )
                (*it).second.Accept( v );
        return;
    default:
        throw csmp::Exception( WARNING, "Model<dim>::Accept",
                               "application level of visitor unresolved, nothing done...");
    }
} // end Accept




/**
  Allows the input of properties via the basic CSP variable types
  ScalarVariable, VectorVariable and TensorVariable. Any previous
  values are overwritten irrespective of the flags that the current
  values had. Thus, boundary conditions need to be reassigned after this
  method has been called.

  Note that the method checks the input variables against the ranges
  specified in the variables database '*.txt' file. If the values
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

  @todo (1-F) New placements not supported?
*/
template<size_t dim>
template<class T>
void Model<dim>::InputPropertyValue( const char* input_prop, const T& value )
  {
    csmp::Index prop_key = this->database_.StorageKey(input_prop);
    if ( prop_key.place == FACE || prop_key.place == BOUNDARY || prop_key.place == FACE_INTEGRATION_POINT ||
         prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT )
      {
      for (typename map<std::string,csmp::Boundary<dim> >::iterator
        it = this->BoundariesBegin(); it != this->BoundariesEnd(); ++it)
        (*it).second.InputPropertyValue( input_prop, value, COMPLETE );
      }
    else if ( prop_key.place == INTER_FACE || prop_key.place == SPLIT_BOUNDARY || prop_key.place == INTER_FACE_INTEGRATION_POINT ||
              prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
      {
        throw csmp::Exception( ERROR, "Model<dim>::InputPropertyValue", "SplitBoundary/InterFace properties not supported here yet" );
      }
    else if ( prop_key.place == MODEL )
      {
        this->Store( prop_key, value );
      }
    else if ( prop_key.place == REGION )
      {
      for (typename map<string,csmp::Region<dim> >::iterator
        it = this->UniqueRegionsBegin(); it != this->UniqueRegionsEnd(); ++it)
        (*it).second.InputPropertyValue( input_prop, value, COMPLETE );
      for (typename map<string,csmp::Region<dim> >::iterator
        it = this->RegionsBegin(); it != this->RegionsEnd(); ++it)
        (*it).second.InputPropertyValue( input_prop, value, COMPLETE );
      }
    else 
      this->Region("Model").InputPropertyValue( input_prop, value, COMPLETE );
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
template<size_t dim>
template<class T>
void Model<dim>::InputBoundaryValue( BOX_BOUNDARY boundary, const char* input_prop, const T& value )
 {
    csmp::Index  prop_key = database_.StorageKey(input_prop);

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryValue: input variable:", input_prop,
                                    "method has not been implemented for TensorVariables.");

    if ( prop_key.place != NODE )
      throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryValue: input variable:", input_prop,
                                    "method can only be applied to node variables.");

    csmp::Region<dim>&  super_group(this->Region("Model"));

    if ( isSide(boundary) ) {
         for ( typename vector<Node<dim>*>::iterator
               nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
           if ( belongsToSide( boundary, (*nit)->AtBoundary() ) ) (*nit)->Store( prop_key, value );
         return;
      }

    if ( isEdge(boundary) ) {
         for ( typename vector<Node<dim>*>::iterator
               nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
           if ( belongsToEdge( boundary, (*nit)->AtBoundary() ) ) (*nit)->Store( prop_key, value );
         return;
      }

    if ( isCorner(boundary) ) {
         for ( typename vector<Node<dim>*>::iterator
               nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
           if ( (*nit)->AtBoundary() == boundary ) (*nit)->Store( prop_key, value );
         return;
      }

    throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryValue: ",
                              "boundary flag could not be parsed.");

 } // end InputBoundaryValue


template void Model<1U>::InputBoundaryValue<ScalarVariable >( BOX_BOUNDARY, const char*, const ScalarVariable& );
template void Model<2U>::InputBoundaryValue<ScalarVariable >( BOX_BOUNDARY, const char*, const ScalarVariable& );
template void Model<3U>::InputBoundaryValue<ScalarVariable >( BOX_BOUNDARY, const char*, const ScalarVariable& );

template void Model<1U>::InputBoundaryValue<VectorVariable<1U> >( BOX_BOUNDARY, const char*, const VectorVariable<1U>& );
template void Model<2U>::InputBoundaryValue<VectorVariable<2U> >( BOX_BOUNDARY, const char*, const VectorVariable<2U>& );
template void Model<3U>::InputBoundaryValue<VectorVariable<3U> >( BOX_BOUNDARY, const char*, const VectorVariable<3U>& );







template<size_t dim>
void Model<dim>::InputBoundaryFlags( BOX_BOUNDARY boundary, const char* property, const std::vector<VARIABLE_FLAG>& flags )
 {
    csmp::Index  prop_key = database_.StorageKey(property);

    if ( flags.empty() )
      throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags: ",
                                "boundary flag vector is empty, nothing was done");

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags: ",
                                "method has not been implemented for TensorVariables");

    if ( prop_key.place != NODE )
      throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags: ",
                                "method can only be applied to node variables");

    csmp::Region<dim>&  super_group(this->Region("Model"));

    if ( prop_key.type == SCALAR ) {
          if ( isSide(boundary) ) {
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
                 if ( belongsToSide( boundary, (*nit)->AtBoundary() ) ) (*nit)->Status( prop_key, flags[0] );
               return;
            }

          if ( isEdge(boundary) ) {
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
                 if ( belongsToEdge( boundary, (*nit)->AtBoundary() ) ) (*nit)->Status( prop_key, flags[0] );
               return;
            }

          if ( isCorner(boundary) ) {
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
                 if ( (*nit)->AtBoundary() == boundary ) (*nit)->Status( prop_key, flags[0] );
               return;
            }
      }
    else { // vector variables
          if ( flags.size() != dim )
            throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags (vector variable): ",
                                      "boundary flag vector has the wrong number of entries");
          if ( isSide(boundary) ) {
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
                 if ( belongsToSide( boundary, (*nit)->AtBoundary() ) )
                   for ( size_t i=0U; i<dim; i++ ) (*nit)->Status( prop_key, i, flags[i] );
               return;
            }

          if ( isEdge(boundary) ) {
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
                 if ( belongsToEdge( boundary, (*nit)->AtBoundary() ) )
                   for ( size_t i=0U; i<dim; i++ ) (*nit)->Status( prop_key, i, flags[i] );
               return;
            }

          if ( isCorner(boundary) ) {
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ )
                 if ( (*nit)->AtBoundary() == boundary )
                   for ( size_t i=0U; i<dim; i++ ) (*nit)->Status( prop_key, i, flags[i] );
               return;
            }
      }
    throw csmp::Exception( ERROR, "Model<dim>::InputBoundaryFlags: ",
                              "boundary flag could not be handled.");
 } // end InputBoundaryFlags





/**
  Uses linear interpolation to interpolate boundary values across model boundary.

  Last argument is a vector because there may either be two or four boundary endpoints
  depending on whether a 2 or three-dimensional model is used.

  @attention This method is about to be deprecated.
*/
template<size_t dim>
void Model<dim>::InterpolateBoundaryValues( BOX_BOUNDARY side, const char* input_prop,
                                            const vector<ScalarVariable>& bvalues )
 {
    csmp::Index   prop_key = database_.StorageKey(input_prop);

    if ( prop_key.type != SCALAR )
      throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                                "this method can only be applied to scalar variables");

    if ( prop_key.place != NODE )
      throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                                "method can only be applied to node variables");

    // checking that the boundary flags are all consistent
    assert( !bvalues.empty() );
    ScalarVariable  res(bvalues[0]); // getting the variable flag
    for ( typename vector<ScalarVariable >::const_iterator
          vit=bvalues.begin(); vit!=bvalues.end(); vit++ )
      if ( (*vit).Flag() != res.Flag() )
        throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                                       "inconsistent flagging of scalar boundary values");

    csmp::Region<dim>&  super_group(this->Region("Model"));

    assert( bvalues.size() >= 2U );
    double64  v1 = bvalues[0].Value();
    double64  v2 = bvalues[1].Value();

    csmp::Point<dim>  xyz_min, xyz_max;
    MinMaxCoordinates( xyz_min, xyz_max );

    // 2D case
    if ( dim == 2U ) {
         assert( isSide(side) );
         boundaryMinMaxCoordinates( side, xyz_min, xyz_max );
         for ( typename vector<Node<dim>*>::iterator
               nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
               res = linearInterpolate( make_pair(xyz_min,v1), make_pair(xyz_max,v2), (*nit)->Coordinate() );
               (*nit)->Store( prop_key, res );
            }
         return;

      } // end 2D

    assert( !isCorner(side) );

    if ( isSide(side) ) {
        assert( bvalues.size() == 4U );
        double64  v3 = bvalues[2].Value();
        double64  v4 = bvalues[3].Value();

        switch( side )
          {
             case LEFT:   // YZ PLANE
               boundaryMinMaxCoordinates( LEFT, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case RIGHT:  // YZ PLANE
               boundaryMinMaxCoordinates( RIGHT, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res   = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case BACK:   // XY PLANE
               boundaryMinMaxCoordinates( BACK, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res   = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case FRONT:  // XY PLANE
               boundaryMinMaxCoordinates( FRONT, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res   = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case TOP:    // XZ PLANE
               boundaryMinMaxCoordinates( TOP, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res   = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case BOTTOM: // XZ PLANE
               boundaryMinMaxCoordinates( BOTTOM, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res   = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1, v2, v3, v4 );
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
         for ( typename vector<Node<dim>*>::iterator
               nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
              res = linearInterpolate( make_pair(xyz_min,v1), make_pair(xyz_max,v2), (*nit)->Coordinate() );
              (*nit)->Store( prop_key, res );
           }
      }

  } // end InterpolateBoundaryValues (scalar version)






///  vector variable version of previous method
template<size_t dim>
void Model<dim>::InterpolateBoundaryValues( BOX_BOUNDARY side, const char* input_prop,
                                            const vector<VectorVariable<dim> >& bvalues )
 {
    csmp::Index  prop_key = database_.StorageKey(input_prop);

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                                "method has not been implemented for tensor variables");

    if ( prop_key.place != NODE )
      throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                                "method can only be applied to node variables");

    // checking that the boundary flags are all consistent
    assert( !bvalues.empty() );
    VectorVariable<dim>  res(bvalues[0]); // getting the flags
    for ( typename vector<VectorVariable<dim> >::const_iterator
          vit=bvalues.begin(); vit!=bvalues.end(); vit++ )
      for ( size_t i=0U; i<dim; i++ )
        if ( (*vit).Flag(i) != res.Flag(i) )
          throw csmp::Exception( ERROR, "Model<dim>::InterpolateBoundaryValues: ",
                                   "inconsistent flagging of boundary vector variables");

    csmp::Region<dim>&  super_group(this->Region("Model"));

    assert( bvalues.size() >= 2U );
    VectorVariable<dim>  v1 = bvalues[0];
    VectorVariable<dim>  v2 = bvalues[1];

    csmp::Point<dim>  xyz_min, xyz_max;
    MinMaxCoordinates( xyz_min, xyz_max );

    // 2D case
    if ( dim == 2U ) {
         assert( isSide(side) );
         boundaryMinMaxCoordinates( side, xyz_min, xyz_max );
         for ( typename vector<Node<dim>*>::iterator
               nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
               res(0) = linearInterpolate( make_pair(xyz_min,v1[0]), make_pair(xyz_max,v2[0]), (*nit)->Coordinate() );
               res(1) = linearInterpolate( make_pair(xyz_min,v1[1]), make_pair(xyz_max,v2[1]), (*nit)->Coordinate() );
               (*nit)->Store( prop_key, res );
            }
         return;

      } // end 2D

    assert( !isCorner(side) );

    if ( isSide(side) ) {
        assert( bvalues.size() == 4U );
        VectorVariable<dim>  v3 = bvalues[2];
        VectorVariable<dim>  v4 = bvalues[3];

        switch( side )
          {
             case LEFT:   // YZ PLANE
               boundaryMinMaxCoordinates( LEFT, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res(0) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
                     res(1) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
                     res(2) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case RIGHT:  // YZ PLANE
               boundaryMinMaxCoordinates( RIGHT, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res(0) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
                     res(1) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
                     res(2) = bilinearInterpolate( 1, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case BACK:   // XY PLANE
               boundaryMinMaxCoordinates( BACK, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res(0) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
                     res(1) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
                     res(2) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case FRONT:  // XY PLANE
               boundaryMinMaxCoordinates( FRONT, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res(0) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
                     res(1) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
                     res(2) = bilinearInterpolate( 0, 1, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case TOP:    // XZ PLANE
               boundaryMinMaxCoordinates( TOP, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res(0) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
                     res(1) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
                     res(2) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
                     (*nit)->Store( prop_key, res );
                  }
               break;

             case BOTTOM: // XZ PLANE
               boundaryMinMaxCoordinates( BOTTOM, xyz_min, xyz_max );
               for ( typename vector<Node<dim>*>::iterator
                     nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
                     res(0) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[0], v2[0], v3[0], v4[0] );
                     res(1) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[1], v2[1], v3[1], v4[1] );
                     res(2) = bilinearInterpolate( 0, 2, xyz_min, xyz_max, (*nit)->Coordinate(), v1[2], v2[2], v3[2], v4[2] );
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
         for ( typename vector<Node<dim>*>::iterator
               nit=super_group.NodesBegin(); nit!=super_group.NodesEnd(); nit++ ) {
              res(0) = linearInterpolate( make_pair(xyz_min,v1[0]), make_pair(xyz_max,v2[0]), (*nit)->Coordinate() );
              res(1) = linearInterpolate( make_pair(xyz_min,v1[1]), make_pair(xyz_max,v2[1]), (*nit)->Coordinate() );
              res(2) = linearInterpolate( make_pair(xyz_min,v1[2]), make_pair(xyz_max,v2[2]), (*nit)->Coordinate() );
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
template<size_t dim>
void Model<dim>::CopyReplace( const char* from, const char* to )
 {
     this->Region("Model").CopyReplace( from, to );
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
template<size_t dim>
bool  Model<dim>::CopyGradientOfProperty_A_To_B( const char* prop_a, const char* prop_b )
 {
    return this->Region("Model").CopyGradientOfProperty_A_To_B( prop_a, prop_b );
 }  





/**
  InterpolateNodePropertyToElementProperty() interpolates node property to the
  'barycenter' of the triangle. The resulting
  value is different from the result obtained by applying the interrelation
  subclass NodeToElementProperty. The Calculate() method of the latter assigns
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
  InterpolateNodePropertyToElementProperty() will report an error and
  return without completing its task.
*/
template<size_t dim>
void  Model<dim>::InterpolateNodeToElementProperty( const char* nprop, const char* eprop , bool verbose)
 {
    this->Region("Model").InterpolateNodeToElementProperty( nprop, eprop );

    if (verbose){
        cout <<"\nModel<" << dim << ">::InterpolateNodeToElementProperty: ";
        cout <<"'" << nprop <<"' has been successfully interpolated to '"<< eprop <<"'." << endl;
    }

 } // end InterpolateNodeToElementProperty





template<size_t dim>
void  Model<dim>::InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* cpprop )
 {
    this->Region("Model").InterpolateNodeToIntegrationPointProperty( nprop, cpprop );

    cout <<"\nModel<" << dim << ">::InterpolateNodeToIntegrationPointProperty: ";
    cout <<"'" << nprop <<"' has been successfully interpolated to '"<< cpprop <<"'." << endl;

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
template<size_t dim>
void  Model<dim>::InterpolateIntegrationPointToElementProperty( const char* cprop, const char* eprop )
 {
    this->Region("Model").InterpolateIntegrationPointToElementProperty( cprop, eprop );

    cout <<"\nModel<" << dim << ">::InterpolateIntegrationPointToElementProperty: ";
    cout <<"'" << cprop <<"' has been successfully interpolated to '"<< eprop <<"'." << endl;

 } // end InterpolateIntegrationPointToElementProperty



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
template<size_t dim>
void  Model<dim>::ExtrapolateIntegrationPointToNodeProperty( const char* cprop, const char* nprop )
 {
    this->Region("Model").ExtrapolateIntegrationPointToNodeProperty( cprop, nprop );

    cout <<"\nModel<" << dim << ">::ExtrapolateIntegrationPointToNodeProperty: ";
    cout <<"'" << cprop <<"' has been successfully interpolated to '"<< nprop <<"'." << endl;

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
template<size_t dim>
void Model<dim>::MinMaxCoordinates( Point<dim>& xyz_min,
                                    Point<dim>& xyz_max ) const
 {
    this->Region("Model").MinMaxCoordinates( xyz_min, xyz_max );

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
template<size_t dim>
void Model<dim>::Out() const
 {
    cout <<"\n\n\n\nModel<"<< dim <<">::Out: ";
    database_.Out();
    fem_manager_.Out();
    mesh_manager_.Out();

    cout <<"\nunique Regions: ";
    for ( typename map<string,csmp::Region<dim> >::const_iterator
          gr_it=this->uniqueGroupMap_.begin(); gr_it!=this->uniqueGroupMap_.end(); gr_it++ )
      {
         cout <<"\n\t"<< (*gr_it).first;
         (*gr_it).second.Out();
      }

    if ( !this->groupMap_.empty() ) {
        cout <<"\n\n\nnon-unique Regions: ";
        for ( typename map<string,csmp::Region<dim> >::const_iterator
              gr_it=this->groupMap_.begin(); gr_it!=this->groupMap_.end(); gr_it++ )
          {
             cout <<"\n\t"<< (*gr_it).first;
             (*gr_it).second.Out();
          }
      }

     // boundaries
     if ( this->Boundaries() != 0U ) {
        cout <<"\n\n\nBoundaries: "<< endl;
        for (typename map<std::string,csmp::Boundary<dim> >::const_iterator
             it = this->BoundariesBegin(); it != this->BoundariesEnd(); it++ )
          {
             cout <<"\n\t"<< (*it).first;
             (*it).second.Out();
          }
      }

    // split boundaries
    if ( this->SplitBoundaries() != 0U ) {
        cout <<"\n\n\nSplitBoundaries: "<< endl;
        for ( typename map<std::string,csmp::SplitBoundary<dim> >::const_iterator
              it = this->SplitBoundariesBegin(); it != this->SplitBoundariesEnd(); it++ )
          {
             cout <<"\n\t"<< (*it).first;
             (*it).second.Out();
          }
      }

 } // end Out






/**
  AssignElementCharacteristicsTo() allows to assign a number of Element
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

  AssignElementCharacteristicsTo() can be used to:
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
template<size_t dim>
void Model<dim>::AssignElementCharacteristicsTo( const char* characteristic, const char* var )
 {
     this->Region("Model").AssignElementCharacteristicsTo( characteristic, var ) ;

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
template<size_t dim>
void Model<dim>::MoveNodeCoordinatesBy( const char* vector_variable )
 {
     this->Region("Model").MoveNodeCoordinatesBy( vector_variable );
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
template<size_t dim>
void Model<dim>::AssignNodeCoordinatesTo( const char* vector_variable )
 {
    this->Region("Model").AssignNodeCoordinatesTo( vector_variable );
 } 



template<size_t dim>
void Model<dim>::AssignNodeCoordinatesTo( const char* scalar_variable, char c )
 {
    this->Region("Model").AssignNodeCoordinatesTo( scalar_variable, c );
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
template<size_t dim>
void Model<dim>::Apply( PDE_Integrator<dim,csmp::Region>& problem, bool debug )
 {
    problem.IntegrateOver( this->Region("Model"), debug );

 } 


/// application to all boundaries
template<size_t dim>
void Model<dim>::Apply( PDE_Integrator<dim,csmp::Boundary>& problem, bool debug )
 {
    for ( typename map<std::string,csmp::Boundary<dim> >::iterator it = this->BoundariesBegin(); it != this->BoundariesEnd(); ++it )
      problem.IntegrateOver( (*it).second, debug );
 } 


/// application to all boundaries
template<size_t dim>
void Model<dim>::Apply( PDE_Integrator<dim,csmp::SplitBoundary>& problem, bool debug )
 {
    for ( typename map<std::string,csmp::SplitBoundary<dim> >::iterator it = this->SplitBoundariesBegin(); it != this->SplitBoundariesEnd(); ++it )
      problem.IntegrateOver( (*it).second, debug );

 } 


/**
     specific regions
*/
template<size_t dim>
void Model<dim>::Apply( PDE_Integrator<dim,csmp::Region>& problem, const char* region_name, bool debug )
 {
    problem.IntegrateOver( this->Region(region_name), debug );

 } 


/// PDE solution applied to a specific boundary
template<size_t dim>
void Model<dim>::Apply( PDE_Integrator<dim,csmp::Boundary>& problem, const std::string& boundary_name, bool debug )
 {
    problem.IntegrateOver( this->Boundary(boundary_name), debug );

 } // Apply (PDE_Integrator)



/// PDE solution applied to a specific split boundary
template<size_t dim>
void Model<dim>::Apply( PDE_Integrator<dim,csmp::SplitBoundary>& problem, const std::string& splitboundary_name, bool debug )
 {
    problem.IntegrateOver( this->SplitBoundary(splitboundary_name), debug );
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
template<size_t dim>
void Model<dim>::Apply( Interrelation<dim>& relation, const char* region )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( this->uniqueGroupMap_.empty() ) {
         csmp_error.notice( ERROR, "Model<dim>::Apply",
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
    if ( this->uniqueGroupMap_.find(region) != this->uniqueGroupMap_.end() ) {
         this->Region(region).Apply( relation );
         return;
      }

    for ( typename map<std::string,csmp::Region<dim> >::iterator
          it=this->uniqueGroupMap_.begin(); it!=this->uniqueGroupMap_.end(); it++ )
       (*it).second.Apply( relation );

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
template<size_t dim>
void Model<dim>::OutputVariableToScreen( const char* prop ) const
 {
     this->Region("Model").OutputVariableToScreen( prop );

 } // end OutputVariableToScreen



// this wants to be a lambda function in the next method
/// recovering and sorting to find minimum and maximum Eigen values
template<size_t dim>
void minMaxEigenValues( const TensorVariable<dim>& ts, double64& tmin, double64& tmax )
 {
    VectorVariable<dim>  evals;
    ts.EigenValues( evals );
    std::set<double64> min_max;
    for ( size_t i=0U; i<dim; i++ ) min_max.insert( evals[i] );
    tmin = (*min_max.begin());
    tmax = (*min_max.rbegin());
 }




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
template<size_t dim>
void Model<dim>::MinMaxOf( const char* prop, double64& vmin, double64& vmax ) const
 {
    ErrorHandler&  csmp_error(ErrorHandler::Instance());
 
    if ( !database_.IsDefined(prop) ) {
          csmp_error.notice( ERROR, "Model<dim>::MinMaxOf", prop, "is undefined; nothing could be done" );
           return;
      }
    csmp::Index  prop_key(database_.StorageKey(prop));
   
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
               ArrayVariable  a(prop_key.dataDepth);
               this->Read( prop_key, a );
               a.MinMax( vmin, vmax );
            }
          else if ( prop_key.type == FLAGGEDARRAY ) {
               FlaggedArrayVariable  fa(prop_key.dataDepth);
               this->Read( prop_key, fa );
               fa.MinMax( vmin, vmax );
            }
          return;
      }


    // standard properties present everywhere in the model
    if ( prop_key.place == NODE || prop_key.place == ELEMENT || prop_key.place == ELEMENT_INTEGRATION_POINT ||
         prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == SECTOR_INTEGRATION_POINT ) {
       this->Region("Model").MinMaxOf( prop_key, vmin, vmax );
       return;
    }

   
    // unique and non-unique regions 
    if ( prop_key.place == REGION ) {
         typename std::map<std::string,csmp::Region<dim> >::const_iterator  git(this->UniqueRegionsBegin());
         (*git).second.MinMaxOf( prop_key, vmin, vmax );
         double64  gmin(vmin), gmax(vmax);
         while ( git!=this->UniqueRegionsEnd() ) {
              (*git).second.MinMaxOf( prop_key, vmin, vmax );
              gmin = std::min(gmin,vmin);
              gmax = std::max(gmax,vmax);
              git++;
           }
         for ( typename std::map<std::string,csmp::Region<dim> >::const_iterator
               git=this->RegionsBegin(); git!=this->RegionsEnd(); git++ ) {
              (*git).second.MinMaxOf( prop_key, vmin, vmax );
              gmin = std::min(gmin,vmin);
              gmax = std::max(gmax,vmax);
           }
         vmin = gmin;
         vmax = gmax;
         return;
      }

   
    // boundaries
    if ( prop_key.place == BOUNDARY || prop_key.place == FACE || prop_key.place == FACE_INTEGRATION_POINT ||
         prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT ) {
         typename map<std::string,csmp::Boundary<dim> >::const_iterator  git( this->BoundariesBegin() );
         (*git).second.MinMaxOf( prop_key, vmin, vmax );
         double64  gmin(vmin), gmax(vmax);
         while ( git != this->BoundariesEnd() ) {
               (*git).second.MinMaxOf( prop_key, vmin, vmax );
               gmin = std::min(gmin,vmin);
               gmax = std::max(gmax,vmax);
               git++;
            }
         vmin = gmin;
         vmax = gmax;
         return;
      }


    // split boundaries
    if ( prop_key.place == SPLIT_BOUNDARY || prop_key.place == INTER_FACE || prop_key.place == INTER_FACE_INTEGRATION_POINT ||
         prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
        typename map<std::string,csmp::SplitBoundary<dim> >::const_iterator  git( this->SplitBoundariesBegin() );
        (*git).second.MinMaxOf( prop_key, vmin, vmax );
        double64  gmin(vmin), gmax(vmax);
        while ( git != this->SplitBoundariesEnd() ) {
            (*git).second.MinMaxOf( prop_key, vmin, vmax );
            gmin = std::min(gmin,vmin);
            gmax = std::max(gmax,vmax);
            git++;
          }
        vmin = gmin;
        vmax = gmax;
     }

 } // end MinMaxOf




/// Instantiates MODEL LocalVariables
template<size_t dim>
void Model<dim>::InitializeLocalVariableStorage()
  {
    this->ResizePropertyStorage( Database().LocalVariablesAt(MODEL) );
  }


/**
  Property storage for the geometric primitives that form the mesh is created
  before Regions, Boundaries, and / or SplitBoundaries are formed.

  This method fulfills the job at a later stage, updating only if such
  variables exist.
*/
template<size_t dim>
bool Model<dim>::UpdateSubdomainPropertyStorage()
 {
    bool  increased_storage(false);

    // 1. Regions
    if ( this->Database().VariableCount(REGION) > 0 ) {
        for ( typename map<string,csmp::Region<dim> >::iterator
              git=this->UniqueRegionsBegin(); git!=this->UniqueRegionsEnd(); git++ )
          (*git).second.ResizePropertyStorage( this->Database().LocalVariablesAt(REGION) );
        for ( typename map<string,csmp::Region<dim> >::iterator
              git=this->RegionsBegin(); git!=this->RegionsEnd(); git++ )
          (*git).second.ResizePropertyStorage( this->Database().LocalVariablesAt(REGION) );
        increased_storage = true;
      }

    // 2. Boundaries
    if ( this->Database().VariableCount(BOUNDARY) > 0 ) {
         for ( typename map<std::string,csmp::Boundary<dim> >::iterator
               git = this->BoundariesBegin(); git != this->BoundariesEnd(); git++ )
           (*git).second.ResizePropertyStorage( this->Database().LocalVariablesAt(BOUNDARY) );
      }

    // 3. SplitBoundaries
    if ( this->Database().VariableCount(SPLIT_BOUNDARY) > 0  ) {
      for ( typename map<std::string,csmp::SplitBoundary<dim> >::iterator
        git = this->SplitBoundariesBegin(); git != this->SplitBoundariesEnd(); git++ )
        (*git).second.ResizePropertyStorage( this->Database().LocalVariablesAt(SPLIT_BOUNDARY) );
      }

    return increased_storage;

 }  // end InitializeRegionProperties


 

// ================================================================================================================
//
//            NEW BINARY FILE INPUT / OUTPUT
//
// ================================================================================================================

  
/**
    writes entire model with associated properties / variables to CSMP native set of binary files.
*/
template<size_t dim>
void Model<dim>::OutputToBinaryFile( const char* file_string ) const
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // 0. IO & resetting of indices
    double64& model_time( ModelTime::Instance().modelTime );
    cout << "\nModel<" << dim <<">::OutputToBinaryFile: Saving model '"<< Name();
    cout <<"' at current time level, t = "<< model_time <<" secs." << endl;

    // 1. mesh output: creating a VSet including Face and InterFace objects
    VSet<dim> vset;
    const bool simplices_numbered_in_a_single_sequence(true);
    mesh_manager_.AssignUniqueNumbers( simplices_numbered_in_a_single_sequence );
    mesh_manager_.OutputMeshTo( vset );

    // 2. property output into VSet including Face and InterFace data
    mesh_manager_.OutputStoredVariablesTo( Database(), vset );

    // model properties
    map<string,Index>  properties;
    Database().ListProperties( MODEL, properties );

    // for all model properties
    for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
      {
         // setting the specifications for the property storage (no memory allocation yet)
         PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
         // for the given property type
         const size_t flag_capacity( (*pit).second.flagDepth );
         const size_t data_capacity( (*pit).second.dataDepth );
         // allocating memory to store the property flags and values
         data.Reserve( flag_capacity, data_capacity );

         switch( (*pit).second.type )
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
               csmp_error.notice( ERROR, "Model<dim>::OutputToBinaryFile:",
                                  (*pit).first, "type of Model variable not recognized.");
           }
         // storing the data in the VSet
         vset.AddData( (*pit).first.c_str(), data );
      }
   
    vset.OutputTo( BinaryVsetFileName(file_string).c_str(), model_time );

    // 3. regions: "All Elements" region first appending the unique and then the non-unique regions
    this->OutputAllRegionsToBinary( BinaryRegionsFileName(file_string).c_str() );
   
    // 4. boundaries "All Faces"
    this->OutputAllBoundariesToBinary( BinaryBoundariesFileName(file_string).c_str() );
   
    // 5. splitboundaries "AllInterFaces"
//    this->OutputAllSplitBoundariesToBinary( BinarySplitBoundariesFileName(file_string).c_str() );
   
    // 6. variable specifications through the database
    Database().BinaryOut( BinaryVariablesFileName(file_string).c_str() );
    cout << "\nModel<"<< dim <<">::OutputToBinaryFile: Output of model '"<< Name();
    cout <<"' to CSMP binaries completed successfully.\n\n";

 } // end OutputToBinaryFile

// TESTING
//for ( auto bit=this->BoundariesBegin(); bit!=this->BoundariesEnd(); ++bit )
//  (*bit).second.Out();





/** reads model written by OutputToBinaryFile() including all associated properties
*/
template<size_t dim>
void Model<dim>::InputFromBinaryFile( const char* model_name )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     // 1. The binary data are read into VSet
     double64& model_time( ModelTime::Instance().modelTime );
     VSet<dim>  vset;
     cout <<"\nModel<"<< dim;
     cout <<">::InputFromBinaryFile: Reading '"<< model_name;
     cout <<"' from VSet... " << endl;  
     vset.InputFrom( BinaryVsetFileName(model_name).c_str(), model_time ); 
    
     // 2. initializing the finite-element manager true=isoparametric
     fem_manager_.InitializeElements( dim, vset.OrderOfFiniteElementInterpolationFunctions(), true );
     cout << "\nModel<"<< dim <<">::InputFromBinaryFile: it is assumed that the model is based on 'isoparametric' finite elements.\n\n";

     // 3. rebuilds finite element mesh and associated property storage
     mesh_manager_.Reconstruct( database_, fem_manager_, vset );
   
     // 4. checking whether the FV stencils need to be initialised
     if ( vset.ContainsFiniteVolumeIntegrationPointData() )
       InstantiateFiniteVolumes();

     // 5. assigning properties to mesh (this reads in the properties output to file via Region::OutputTo(VSet) )
     if( !vset.DataEmpty() ) mesh_manager_.InputStoredVariablesFrom( Database(), vset );
     else ErrorHandler::Instance().notice( INFO, "Model<dim>::InputFromBinaryFile:", "No properties found in VSet." );

     // properties and values stored on the model itself
     for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
       {
         // apart from the name string key in the map, PropertyData contains the most important variable specifications
         if ( (*pit).second.Placement() != MODEL ) continue;
         // some checks
         assert( Database().IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(Database().StorageKey( (*pit).first.c_str() ));
         assert( key.place == MODEL );
         assert( (*pit).second.Size() / key.dataDepth == 1U );
        
         switch( key.type )
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
               csmp_error.notice( ERROR, "Model<dim>::InputFromBinaryFile:",
                                  (*pit).first, "type of Model variable not recognized.");
           }
      }

     // 6. reconstruction of the regions
     this->InputAllRegionsFromBinary( BinaryRegionsFileName(model_name).c_str() );
     // making sure that the computational region has been built
     if ( !this->ContainsRegion("Model") )
       throw csmp::Exception( ERROR, "Model<>::InputFromBinaryFile", "Root region 'Model' is not present." );

     // 7. reconstructing the boundaries (TODO: what if there are no boundaries?)
    this->InputAllBoundariesFromBinary( BinaryBoundariesFileName(model_name).c_str() );

     // 8. reconstructing the splitboundaries
//    this->InputSplitBoundariesFromBinary( BinarySplitBoundariesFileName(model_name).c_str() );

    cout << "\nModel<"<< dim <<">::InputFromBinaryFile: input from binaries (file set: "<< model_name <<") completed successfully.\n\n";

 } // end InputFromBinaryFile


// TESTING
//for ( auto bit=this->BoundariesBegin(); bit!=this->BoundariesEnd(); ++bit )
//  (*bit).second.Out();






// ================================================================================================================
//
//            GLOBAL FUNCTIONS INVOLVING MODEL
//
// ================================================================================================================




/**
    Calculates the x, y, z extent of the model and returns these lengths into its arguments.
    The return value is a string that contains the dimensions with explanations.
*/
std::string  boundingBox( const Model<3U>& sg, double64& dim_x, double64& dim_y, double64& dim_z )
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
    dimensions +=",  ";
    dimensions += "y-length: ";
    sprintf( num, "%lf", dim_y );
    dimensions += num;
    dimensions +=",  ";
    dimensions += "z-length: ";
    sprintf( num, "%lf", dim_z );
    dimensions += num;
    dimensions +=" (box-shaped model). ";
   
    return dimensions;

 } // end boundingBox




// explicit instantiations
template class RegionInterface<1U,Model>;
template class RegionInterface<2U,Model>;
template class RegionInterface<3U,Model>;

template class BoundaryInterface<1U,Model>;
template class BoundaryInterface<2U,Model>;
template class BoundaryInterface<3U,Model>;

template class SplitBoundaryInterface<1U,Model>;
template class SplitBoundaryInterface<2U,Model>;
template class SplitBoundaryInterface<3U,Model>;

template class Model<1U>;
template class Model<2U>;
template class Model<3U>;

} // end namespace csmp





