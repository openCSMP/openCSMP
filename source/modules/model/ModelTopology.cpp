#include "ModelTopology.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"
#include "Standard_IO_Handler.h"
#include "VSet.h"
#include "Box.h"
#include "ConsecutiveSequenceChecker.h"

using namespace std;

namespace csmp {

ModelTopology::ModelTopology( bool isoparametric_element_mesh )
 : model_name("not initialized"),  
   isoparametric_mesh(isoparametric_element_mesh)
 {
 }

ModelTopology::ModelTopology( const char* model_name,
                              bool isoparametric_element_mesh )
 : model_name(model_name),  
   isoparametric_mesh(isoparametric_element_mesh)
 {
 }

ModelTopology& ModelTopology::operator=( const ModelTopology& mt )
 {
    if ( &mt != this ) {
         model_name         = mt.model_name;
         model_regions      = mt.model_regions;
         isoparametric_mesh = mt.isoparametric_mesh;
      }
    return *this;
 }


ModelTopology::ModelTopology( const ModelTopology& mt )
 {
    *this = mt;
 }

ModelTopology::~ModelTopology()
 {
 }


/** Writes all the currently stored topological information to standard output, i.e. the screen.
*/
void  ModelTopology::Out() const
 {
    if ( model_regions.empty() ) {
         std::cout <<"\nModelTopology::Out: Topology of '"<< model_name;
         std::cout <<"' is not defined."<< std::endl;
         return;
      }

    std::cout <<"\nModelTopology::Out: Model: '"<< model_name <<"' with "<< Elements() <<" elements."<< std::endl;
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      {
         std::cout <<"\nRegion: '"<< (*it).first <<"' of ";
         for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
               sit!=(*it).second.first.end(); sit++ ) std::cout << (*sit) <<" ";
         if ( isoparametric_mesh )
           std::cout <<" isoparametric finite elements";
         else
           std::cout <<" finite elements";
         std::cout <<"\nNumber of elements in region: "<< (*it).second.second.size();
         std::cout <<"\nElement ID numbers: "<< std::endl;
         for ( std::vector<size_t>::const_iterator
               lit=(*it).second.second.begin(); lit!=(*it).second.second.end(); lit++ )
           std::cout << (*lit) <<" ";
         std::cout << std::endl;
      }

    std::cout << std::endl;

 } // end Out()


/**

Writes all the currently stored topological information to a user
specified text file. If no information is present when the method
is called, a message is printed to screen in stead of writing
an output file.

@param output_file The name of the textfile which will either be newly created or
overwritten.
*/
void  ModelTopology::Out( const char* output_file ) const
 {
    if ( model_regions.empty() ) {
         std::cout <<"\nModelTopology::Out: Topology of '"<< model_name;
         std::cout <<"' is not defined. No output is written to '"<< output_file <<"'"<< std::endl;
         return;
      }
    std::ofstream  ofs( output_file );

    ofs <<"\nModelTopology::Out: Model: '"<< model_name <<"'"<< std::endl;
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      {
         ofs <<"\nRegion: '"<< (*it).first <<"' of ";
         for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
               sit!=(*it).second.first.end(); sit++ ) ofs << (*sit) <<" ";
         ofs <<" finite elements";
         ofs <<"\nNumber of elements in region: "<< (*it).second.second.size();
         ofs <<"\nElement ID numbers: "<< std::endl;
         for ( std::vector<size_t>::const_iterator
               lit=(*it).second.second.begin(); lit!=(*it).second.second.end(); lit++ )
           ofs << (*lit) <<" ";
         ofs << std::endl;
      }

    ofs << std::endl;
 }





void   ModelTopology::ModelName( const char* name )
  {
     model_name = name;
  }
  
  
std::string ModelTopology::ModelName() const
 {
    return std::string(model_name);
 }

size_t  ModelTopology::ModelRegions() const
 {
    return model_regions.size();
 }



/**

Loops over the contained elements testing whether the model consists
entirely of line elements, i.e., linear, quadratic or cubic bars.  

@return If the model consists entirely of line elements, LineModel() returns
true, else it returns false.  

*/
bool ModelTopology::LineModel() const
 {
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
            sit!=(*it).second.first.end(); sit++ )
        if ( fem_specs::SurfaceElement( (*sit) ) || fem_specs::VolumeElement( (*sit) ) ) return false;
    return true;
 }


/**
 
If the model consists entirely of line and surface type elements it
is considered as a surface model. This is checked by the SurfaceModel().
Surface elements are triangles and quadrilaterals.  
 
@return bool reporting whether the model consists entirely of surfaces.

The method returns true if the model consists entirely of line and 
surface elements, else false is returned.  
*/
bool ModelTopology::SurfaceModel() const
 {
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
            sit!=(*it).second.first.end(); sit++ )
        if ( fem_specs::VolumeElement( (*sit) ) ) return false;
    return true;
 }




/**
 
The models is considered a solid if it contains volume elements, even 
if these do not enclose contiguous domains and are mixed with surfaces.
 

@return If the model contains volume elements, i.e., tetrahedra, hexahedra,
pyramids or prisms, the method returns true, else false.  
*/
bool ModelTopology::SolidModel() const
 {
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
            sit!=(*it).second.first.end(); sit++ )
        if ( fem_specs::VolumeElement( (*sit) ) ) return true;
      
    return false;
 } 



/**

The minimum spatial dimension of a model region is taken to be one if it
consists of line elements, and 2 or 3 if it is made up entirely of surface 
or volume elements, respectively. For regions that contain lines, surfaces
and or volumes, the elements of highest spatial dimensionality determine
the result. Still, this result is only a minimum estimate, because 
line and surface elements also exist in 3D space.  

@param region The method takes to model region name as input argument.

@return  The method returns the minimum spatial dimension of the target region.
For line elements this is 1, for surfaces 2, and for volumes 3.  

If the region cannot be found, the spatial dimension of the model will
be returned in stead.  

@section messages Messages

When the region cannot be found, the method will report an error.  
*/
size_t  ModelTopology::MinimumSpatialDimensionOfRegion( const char* region ) const
 {
    // 1. finding the target region in the current model
    std::string target_region(region);
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator  it(model_regions.find(target_region));
    
    if ( it == model_regions.end() )
      throw csmp::Exception( ERROR, "ModelTopology::MinimumSpatialDimensionOfRegion",
                            "Target region could not be found; returning dim of entire model." );
    
    // if the target region exists, its element types are used to establish
    // its minimum spatial dimension
    size_t min_dim(0U);
     
    for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
          sit!=(*it).second.first.end(); sit++ )
      min_dim = std::min( min_dim, fem_specs::MinimumSpatialDimension(*sit) );
      
    return min_dim;  

 } // end MinimumSpatialDimensionOfRegion




/**

A finite element mesh can either be treated as isoparametric (using
local interpolation functions and numeric integration) or as
standard. CSMP uses different element type names to differentiate
the two. Thus it is important to specify the desired approach.

@return This method returns whether the current model is treated as
isoparametric or not.
*/
bool ModelTopology::IsoparametricElements() const
 {
    return isoparametric_mesh;

 } // end IsoparametricElements




/**

Instructs the ModelTopology class to treat the contained finite elements
as isoparametric, i.e. as finite elements with local interpolation functions
and using numeric integration.
*/
void ModelTopology::TreatElementsAsIsoparametric()
 {
    isoparametric_mesh = true;
 }



/**
Returns the order of the finite-element interpolation functions of
the current finite-element model.

@return If all elements have the same interpolation order, it is returned
as an unsigned integer (1=linear, 2=quadratic, 3=cubic). If the order
varies 0 is returned.

@section implementation Implementation

Cubic interpolation is currently not tested for.

@section messages Messages

The method tests whether elements of different order appear in the mesh
(a no-no) and report this as an error.
 */
size_t ModelTopology::InterpolationOrder() const
 {
    bool linear_ipol(false), quadratic_ipol(false), cubic_ipol(false);

    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
            sit!=(*it).second.first.end(); sit++ ) {
           if ( fem_specs::LinearElement( fem_specs::CSMP_Type(*sit) ) )    linear_ipol    = true;
           if ( fem_specs::QuadraticElement( fem_specs::CSMP_Type(*sit) ) ) quadratic_ipol = true;
           if ( fem_specs::CubicElement( fem_specs::CSMP_Type(*sit) ) )     cubic_ipol     = true;
        }

    if ( linear_ipol && quadratic_ipol )
      std::cout <<"\nModelTopology::InterpolationOrder: Warning: \
                Model contains linear and quadratic elements at the same time !"<< std::endl;

    if ( quadratic_ipol && cubic_ipol )
      std::cout <<"\nModelTopology::InterpolationOrder: Warning: \
                Model contains cubic and quadratic elements at the same time !"<< std::endl;

    if ( linear_ipol && cubic_ipol )
      std::cout <<"\nModelTopology::InterpolationOrder: Warning: \
                Model contains cubic and linear elements at the same time !"<< std::endl;

    if ( linear_ipol )    return 1U;
    if ( quadratic_ipol ) return 2U;
    if ( cubic_ipol )     return 3U;

    std::cout <<"\nModelTopology::InterpolationOrder: ERROR: Interpolation order could not ";
    std::cout <<" be identified for one of the elements; returning 0."<< std::endl;

    return 0U;
 }

/**
 
The method checks the order of the finite element interpolation functions
of the current model.  

@return If all the element types in the current model use linear interpolation
the method returns true, else it returns false.  
 */
bool ModelTopology::LinearElementMesh() const
 {
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
            sit!=(*it).second.first.end(); sit++ )
        if ( !fem_specs::LinearElement( fem_specs::CSMP_Type(*sit) ) ) return false;
    return true;
 }

/**
 
The method checks the order of the finite element interpolation functions
of the current model.  

@return If all the element types in the current model use linear interpolation
the method returns true, else it returns false.  
*/
bool ModelTopology::QuadraticElementMesh() const
 {
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
            sit!=(*it).second.first.end(); sit++ )
        if ( !fem_specs::QuadraticElement( fem_specs::CSMP_Type(*sit) ) ) return false;
    return true;
 }


/**

Returns the names of ANSYS finite-element types used in the current
model to the supplied list of std::strings.

The list is not emptied before the element types are written to
it.

@return A set of integers or std::strings into which the current ANSYS element types will
be inserted.
*/
size_t  ModelTopology::FiniteElementTypes( std::set<std::string>& etypes ) const
 {
     for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
           it=model_regions.begin(); it!=model_regions.end(); it++ )
       for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
             sit!=(*it).second.first.end(); sit++ )
         etypes.insert( (*sit) );
     return etypes.size();
 }

size_t  ModelTopology::FiniteElementTypes( std::set<int32_t>& etypes ) const
 {
     for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
           it=model_regions.begin(); it!=model_regions.end(); it++ )
       for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
             sit!=(*it).second.first.end(); sit++ )
         etypes.insert( fem_specs::CSMP_Type(*sit) );
     return etypes.size();
 }



/** Cheanges the element type from CSMP_Type as specified in ANSYS_ElementSpecifications
    to user defined type.

@section arguments Input Arguments

Use, for example, to convert linear tetrahedra into quadratic tetrahedra. AP.  
*/
void ModelTopology::ChangeElementType( const std::string& old_element_type,
                                       std::string new_element_type )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  size_t         counter(0U);
  
  for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
  {
     std::set<std::string> new_set;
     for ( std::set<std::string>::iterator sit=(*it).second.first.begin();
            sit!=(*it).second.first.end(); sit++ )
       {
          if ( *sit == old_element_type ) {
               new_set.insert(new_element_type);
               counter++;
            }
          else
            new_set.insert(*sit);
       }       
     (*it).second.first = new_set;
  }
  
  if ( counter == 0U )
    csmp_error.notice( ERROR, "ModelTopology::ChangeElementType",
                      "topology did not contain requested element type (should be ANSYS type)" );
}






/**

Gets rid of all the model regions which only contain the target
element type. For poly-element-type model regions the method cannot
determine which of the elements are of the target type. Thus, they
are left untouched.

@param etype A character std::string that matches the element type that is to be removed.


@section messages Messages

The method reports back which regions are being eliminated.
*/
void  ModelTopology::EliminateElementType( const char* etype )
 {
    std::list<std::string>  to_cull;

    // erasing the element ids from 'elmt_ids' and flagging keys of 'model_regions'
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      if ( (*it).second.first.size() == 1U && (*(*it).second.first.begin()) == etype )
        {
            std::cout <<"\nModelTopology::EliminateElementType: ";
            std::cout <<" eliminating '"<< etype <<"' region '"<< (*it).first;
            std::cout <<"' from model topology '"<< model_name <<"'"<< std::endl;
            to_cull.push_back( (*it).first );
        }

    // erasing regions of etype from 'model regions'
    for ( std::list<std::string>::iterator
          dit=to_cull.begin(); dit!=to_cull.end(); dit++ ) model_regions.erase( (*dit) );

 } // end EliminateElementType


/**

Removes the enlisted element types from the current model topology,
irrespective whether this breaks the model up into spatially non-
contiguous domains.

EliminateElementTypes() will leave regions which consist of a range
of elements untouched, since it cannot identify which element IDs
correspond to the enlisted types.

@param etypes The method takes a list of std::strings as argument which define the ANSYS
element types.

@section application Application .

Use with utmost care !

RenumberElementsConsecutively() should be called after the element
elimination to re-establish contiguous element numbering in the
model. Corresponding operations must also be carried out on the
VSet that stores the element connectivity which has to be updated
when elements were removed.
*/
void  ModelTopology::EliminateElementTypes( const std::list<std::string>& etypes )
 {
    std::set<std::string>  existing_etypes;

    // checking what element types the model consists of
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
            sit!=(*it).second.first.end(); sit++ )
      existing_etypes.insert( (*sit) );

    // checking whether the element types which shall be eliminated do actually exist
    std::set<std::string>::const_iterator  sit;
    std::list<std::string>::const_iterator lit;

    for ( lit=etypes.begin(); lit!=etypes.end(); lit++ )
      if ( (sit=existing_etypes.find((*lit))) == existing_etypes.end() ) {
            std::cout <<"\n'"<< (*lit) <<"'"<< std::endl;
            throw csmp::Exception( INFO, "ModelTopology::EliminateElementTypes",
                                     "This element type is not contained in the model" );
        }

    // eliminating the instructed element types if they exist
    for ( lit=etypes.begin(); lit!=etypes.end(); lit++ )
      EliminateElementType( (*lit).c_str() );
 }

/** Removes all line elements from the current model.

@section application Application .

Use carefully ! - while the nodes of line elements are typically shared
with surface or volume elements, if single line wells are present in a
model, the geometry is disrupted.

RenumberElementsConsecutively() should be called after the element
elimination to re-establish contiguous element numbering in the
model. Corresponding operations must also be carried out on the
VSet that stores the element connectivity which has to be updated
when elements were removed.
*/
void  ModelTopology::EliminateLineElements()
 {
    std::list<std::string>  etypes;

    fem_specs::LineElements( etypes );

    for ( std::list<std::string>::const_iterator
          it=etypes.begin(); it!=etypes.end(); it++ )
      EliminateElementType( (*it).c_str() );
 }





/** Removes all surface elements from the current topology.

@section application Application .

It may be desirable to get rid of internal or external surfaces in a
model since their accumulation will add extra terms to the solution
of PDE integrals. If the surface element nodes are shared with volume
elements, the model topology will stay intact. Else the model will be
broken.

RenumberElementsConsecutively() should be called after the element
elimination to re-establish contiguous element numbering in the
model. Corresponding operations must also be carried out on the
VSet that stores the element connectivity which has to be updated
when elements were removed.
*/
void  ModelTopology::EliminateSurfaceElements()
 {
    std::list<std::string>  etypes;

    fem_specs::SurfaceElements( etypes );

    for ( std::list<std::string>::const_iterator
          it=etypes.begin(); it!=etypes.end(); it++ )
      EliminateElementType( (*it).c_str() );
 }




/** Method eliminates all volume elements from the current model topology.

@section application Application .

This operation typically greatly upsets a current model, but it gives
the option to do, for instance, a  fracture only model of a volumetric
mesh of a fractured rock. It is recommended, however, to control
this inside the meshing tool rather than by post-processing, because
a better and smaller mesh can typically be created if it is known
that only surfaces are required.

RenumberElementsConsecutively() should be called after the element
elimination to re-establish contiguous element numbering in the
model. Corresponding operations must also be carried out on the
VSet that stores the element connectivity which has to be updated
when elements were removed.
*/
void  ModelTopology::EliminateVolumeElements()
 {
    std::list<std::string>  etypes;

    fem_specs::VolumeElements( etypes );

    for (std::list<std::string>::const_iterator
           it=etypes.begin(); it!=etypes.end(); it++ )
      EliminateElementType( (*it).c_str() );
 }






// returns count of elements in unique regions
size_t  ModelTopology::Elements() const
 {
    //       region name          etypes-of-region       ids of elements in region
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
    it(model_regions.begin());

    size_t  nelements(0U);
    while ( it!=model_regions.end() ) {
         nelements += (*it).second.second.size();
         it++;
      }
    return nelements;
 }


void  ModelTopology::Elements( std::set<size_t>& eids ) const
 {
    //       region name          etypes-of-region       ids of elements in region
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
    it(model_regions.begin());

    while ( it!=model_regions.end() ) {
         for ( std::vector<size_t>::const_iterator
               eit=(*it).second.second.begin(); eit!=(*it).second.second.end(); eit++ )
           eids.insert( (*eit) );
         it++;
      }
 }

size_t  ModelTopology::ElementsOfRegion( const char* region ) const
 {
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator  it(model_regions.find(region));

    if ( it == model_regions.end() ) {
         throw csmp::Exception( INFO, "ModelTopology::ElementsOfRegion",
                                "Target region could not be found" );
         return 0U;
      }
    return (*it).second.second.size();
 }

// is the element with the specified ID (1..n) contained in the given
// region ?
bool  ModelTopology::IsWithinRegion( const char* region, size_t elmt_id ) const
 {
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator  it(model_regions.find(region));

    if ( model_regions.empty() ) {
         throw csmp::Exception( ERROR, "ModelTopology::IsWithinRegion",
                                    "The ModelTopology map of regions is empty" );
      }
    if ( it == model_regions.end() ) {
         throw csmp::Exception( ERROR, "ModelTopology::IsWithinRegion",
                                     region, " does not exist" );
         return false;
      }

    // assumes that the element ids of the target region are ordered consecutively
    if ( std::binary_search( (*it).second.second.begin(), (*it).second.second.end(), elmt_id ) ) return true;

    return false;
 }

/**

The method provides access to the element list of the target region
via a constant iterator.

@param region The name of the region the element IDs of which shall be recovered.

@return A constant iterator set to the beginning of the target list, if the
latter could be found. Otherwise an iterator to the end of the first
element list stored by the topology class is returned.

@section messages Messages

An ERROR is raised if the target region cannot be found.

*/
std::vector<size_t>::const_iterator  ModelTopology::ElementsOfRegionBegin( const char* region ) const
 {
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator  it(model_regions.find(region));

    if ( it == model_regions.end() ) {
         throw csmp::Exception( ERROR, "ModelTopology::ElementsOfRegionBegin",
                                    "Target region could not be found" );

         return (*model_regions.begin()).second.second.end();
      }
    return (*it).second.second.begin();
 }




/**

The method provides access to the element list of the target region
via a constant iterator.

@param region The name of the region the element IDs of which shall be recovered.

@return A constant iterator set to the end of the target list, if the
latter could be found. Otherwise an iterator to the end of the first
element list stored by the topology class is returned.

@section messages Messages

A ERROR is raised if the target region cannot be found.
*/
std::vector<size_t>::const_iterator  ModelTopology::ElementsOfRegionEnd( const char* region ) const
 {
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator  it(model_regions.find(region));

    if ( it == model_regions.end() ) {
         throw csmp::Exception( ERROR, "ModelTopology::ElementsOfRegionEnd",
                                      "Target region could not be found" );

         return (*model_regions.begin()).second.second.end();
      }
    return (*it).second.second.end();
 }

/**

Returns the names of the finite element types that constitute the region
of interest as std::strings.

@param region The name of the region that shall be investigated.

A set of unique finite-element type names of the elements that constitute
the regions of interest.

*/
void  ModelTopology::ElementTypesOfRegion( const char* region, std::set<std::string>& etypes ) const
 {
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator  it(model_regions.find(region));


    if ( it == model_regions.end() ) {
         throw csmp::Exception( ERROR, "ModelTopology::ElementTypeOfRegion",
                                    "Target region could not be found" );

         return;
      }

    if ( !etypes.empty() ) etypes.erase( etypes.begin(), etypes.end() );

    for ( std::set<std::string>::const_iterator sit=(*it).second.first.begin();
          sit!=(*it).second.first.end(); sit++ )
      etypes.insert( (*sit) );
 }



















/*
    Manipulation with Regions
 */



/**
 
ReduceToRegions() reads an ASCII (text) input file to determine the names
of the model regions that shall be retained by the current model 
toopology. The name of the input file is created from the user supplied
character std::string by appending '-regions.txt'.

The regions file consists of the following: 
 
1. headline with file name, model name etc.

2. properties whose specific values shall be listed after region name
(here the tab key separates the property names which can contain whitespace).
At least, the property 'permeability' must be specified.

If the file is used only to identify the regions of interest, the second
line and the property information in the later part of the file will
be ignored. Still the second line is treated as if it would contain property 
data. No region names can be specified in it. We recommend to write 
'no property data' into the second line to clarify if the file is only
used to select output regions. 

3. model regions followed by property values (after the last property value
in each line, arbitrary comments can be placed but will be ignored)
 

@param regions_file The name of the regions file for which an extension '-regions.txt' is
appended.  


@section implementation Implementation

This method does not recognize the standard comments yet.
Thus, if you put # in front of an item the method will still read the
line and try to interpret the following input.  

@section application Application .

In conjunction with input of mesh topologies such as ANSYS or FRED derived
meshes.  

@section messages Messages

The method reports if the input file cannot be opened or if one of the 
target regions is not contained in the model topology. In the latter case
a ERROR is raised. 
 
*/
bool  doesRegionsFileExist( const char* regions_file )
{
    std::string  file_name(regions_file);
    file_name += "-regions.txt";
    std::ifstream ifs(file_name.c_str());
    if( !ifs.good() )
        return false;
    return true;
}


void  readDesiredRegions( const char* regions_file,
                          std::set<std::string>& desired_regions )
{
    char         text_line[256];
    char*        token(0);
    const char* const delims =" ,\t,:,\n,\r";

    strcpy( text_line, regions_file );
    strcat( text_line, "-regions.txt" );
    std::string  file_name(text_line);
    std::ifstream ifs(file_name.c_str());

    if ( !ifs.is_open() )
      throw csmp::Exception( FATAL_ERROR,
                             "readDesiredRegions:", text_line,
                             "ASCII geometry input file could not be opened");

    // 1. reading and discarding file header
    ifs.getline( text_line, 500 );

    // 2. reading the names of the properties
    ifs.getline( text_line, 500 );

    // 3. reading the region descriptors from file
    while ( !ifs.eof() )
      {
         ifs.getline( text_line, 256 );
         if ( strlen(text_line) == 0 ) break;
         // only the first token is used which allows the user to add comments
         token = strtok( text_line, delims );
         if( token!= NULL )
            desired_regions.insert( token );
      }
    ifs.close();
}




/**

Copies the selected regions into the target model topology object. The
target object is emptied before the new regions are inserted.

The target model topology to which the selected regions shall be added is returned
into the second method argument.

@param regions The list of the regions that shall be copied.

@section messages Messages

If the target object is not empty, this will be reported. A ERROR
is raised if one of the selected regions does not exist in the
source object.
*/
void  ModelTopology::ExportSelectionTo( const std::list<std::string>& regions, ModelTopology& mt ) const
 {
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator  it;

    if ( !mt.model_regions.empty() ) {
         std::cout <<"\nModelTopology::ExportSelectionTo: ";
         std::cout <<"Erasing input topology object '"<< mt.ModelName() <<"'"<< std::endl;
         mt.Erase();
      }
    if ( &mt == this ) return;

    // giving the exported selection a name
    mt.model_name  = "subset of topology'";
    mt.model_name += model_name;

    for ( std::list<std::string>::const_iterator
          lit=regions.begin(); lit!=regions.end(); lit++ )
      if ( (it=model_regions.find((*lit))) == model_regions.end() )
         throw csmp::Exception( ERROR, "\nModelTopology::ExportSelectionTo",
                                    "Target region does not exist, e.g. ", (*lit).c_str() );
      else
      mt.AddRegion( (*it).first.c_str(), (*it).second.first, (*it).second.second );
 }

/** Out() returns writes all the current region names in the output list.

The method takes a list as an argument into which the region names
are inserted. If the list is not empty, its contents will be deleted
before the new names are inserted.

@section messages Messages

If the current model does not contain a topology a message is printed.
*/
void  ModelTopology::Out( std::list<std::string>& regions ) const
 {
    if ( model_regions.empty() ) {
         std::cout <<"\nModelTopology::Out: Topology of '"<< model_name;
         std::cout <<"' is not defined."<< std::endl;
         return;
      }
    if ( !regions.empty() ) regions.erase( regions.begin(), regions.end() );

    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      regions.push_back( (*it).first );
 }

/**

Contains() allows the user to test whether the current model contains
a certain region identified by the supplied name.

@param region Enter the name of the region that you are looking for.

@return The methods returns 'true' if the region can-, and 'false' if the region
cannot be found in the current model topology.

@section messages Messages

A CSMP error is printed if the region was found but is empty, i.e.
contains no elements.
*/
bool  ModelTopology::Contains( const char* region ) const
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( region == NULL ) {
         csmp_error.notice( WARNING, "ModelTopology::Contains",
                           "method was passed empty const char* std::string" );
         return false;
      }

    if ( model_regions.find(region) != model_regions.end() ) return true;

    return false;
 }




/**

Returns all the dynamically allocated memory to the operating system.
Consequently all the information about the current model is lost.

The name of the model topology is set to 'erased' to signify that
all data are gone.

@section application Application .

Use Erase() to free up memory once the construction of a model has been
completed, i.e. a Model object has been formed.
*/
void  ModelTopology::Erase()
 {
    model_regions.erase( model_regions.begin(), model_regions.end() );
    model_name = "erased";
 }


/** If the target region exists it is removed from the current model.

@param region The name of the region that shall be removed.

@section messages Messages

If the region could not be found in the current topology a warning
message is printed to stdout.
*/
void  ModelTopology::RemoveRegion( const char* region )
 {
    std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator  it(model_regions.find(region));

    if ( it == model_regions.end() ) {
         throw csmp::Exception( WARNING, "ModelTopology::RemoveRegion","Nothing was done." );
         return;
      }

    // removing the actual region
    model_regions.erase(region);
 }



void  ModelTopology::RemoveRegions( const std::set<std::string>& undesired_regions )
 {
    for ( std::set<std::string>::const_iterator
          uit=undesired_regions.begin(); uit!=undesired_regions.end(); uit++ )
      RemoveRegion( (*uit).c_str() );
 }



void  ModelTopology::ReduceToRegions( const char* regions_file )
 {
    std::set<std::string> desired_regions;
    readDesiredRegions( regions_file, desired_regions );
    ReduceToRegions( desired_regions );
 } // end ReduceToRegions




/**
    Eliminates all, but the desired regions from the model topology (element id- per region) container
*/
void  ModelTopology::ReduceToRegions( const std::set<std::string>& desired_regions )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    std::set<std::string>::const_iterator ritEnd = desired_regions.end();
    for( std::set<std::string>::const_iterator
         rit = desired_regions.begin(); rit != ritEnd; rit++ )
        if ( !Contains( (*rit).c_str() ) ) {
          std::string errorMessage("unrecognized region in region file");
          std::string token_message("token: ");
          token_message += (*rit);
          csmp_error.notice( ERROR, "ModelTopology::ReduceToRegions:",
                             token_message.c_str(), errorMessage.c_str() );
        }

    std::set<std::string> undesired_regions;
    // eliminating the undesired regions
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      if ( (desired_regions.find((*it).first)) == desired_regions.end() )
        undesired_regions.insert( (*it).first );
    RemoveRegions( undesired_regions );

    if( csmp_error.Verbose() ){
        std::cout <<"\nModelTopology::ReduceToRegions: Remaining regions:\n\n\t\t";
        for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
              it=model_regions.begin(); it!=model_regions.end(); it++ )
          std::cout << (*it).first <<"  ";
        std::cout << std::endl << std::endl;
    }
 }
 
 
 

/** Adds a new region of finite-elements to the current model topology.

@section arguments Input Arguments

The region must be defined by a unique name, a set of unique finite
element types that make up the region and the list of elements that
actually constitute the region in the greater finite element mesh.

@return the method reports 'true' if the region could be added to the
model topology, and 'false' if not.
 */
bool ModelTopology::AddRegion( const char* rname,
                               const std::set<std::string>& fem_types,
                               const std::vector<size_t>& elms )
 {
    std::pair<std::set<std::string>,std::vector<size_t> >  region(fem_types,elms);

    std::pair<std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator,bool>
    it = model_regions.insert( make_pair( rname, region ) );

    // if the region could not be inserted (probably because it exists already)
    if ( it.second == false ) return false;

    return true;
 }

void ModelTopology::AddRegionElementType( const char* rname,
                                          const std::string& fem_type )
 {
    const std::pair<std::set<std::string>,std::vector<size_t> >  empty_region;
    std::pair<std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator,bool>
    it = model_regions.insert( make_pair( rname, empty_region ) );
    (*it.first).second.first.insert( fem_type );
}

void ModelTopology::AddRegionElementTypes( const char* rname,
                                           const std::set<std::string>& fem_types )
 {
    const std::pair<std::set<std::string>,std::vector<size_t> >  empty_region;
    std::pair<std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator,bool>
    it = model_regions.insert( make_pair( rname, empty_region ) );
    (*it.first).second.first.insert( fem_types.begin(), fem_types.end() );
}

void ModelTopology::AddRegionElementId( const char* rname,
                                        size_t elm )
 {
    const std::pair<std::set<std::string>,std::vector<size_t> >  empty_region;
    std::pair<std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator,bool>
    it = model_regions.insert( make_pair( rname, empty_region ) );
    (*it.first).second.second.push_back( elm );
 }

void ModelTopology::AddRegionElementIds( const char* rname,
                                         const std::vector<size_t>& elms )
 {
    const std::pair<std::set<std::string>,std::vector<size_t> >  empty_region;
    std::pair<std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator,bool>
    it = model_regions.insert( make_pair( rname, empty_region ) );
    (*it.first).second.second.insert( (*it.first).second.second.end(), elms.begin(), elms.end() );
 }



bool ModelTopology::AddRegions( const std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >& unique_regions )
{
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    if( csmp_error.Verbose() &&  !unique_regions.empty() )
        std::cout <<"\nModelTopology::AddRegions: Exporting model subregions:\n\n\t";
    size_t  counter(0);
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          rit=unique_regions.begin(); rit!=unique_regions.end(); rit++ )
    {
        if( csmp_error.Verbose() )
        {
             std::cout << (*rit).first <<"->";
             for ( std::set<std::string>::const_iterator sit=(*rit).second.first.begin(); sit!=(*rit).second.first.end(); sit++ )
               std::cout << (*sit) <<"  ";
        }
        AddRegion( (*rit).first.c_str(), (*rit).second.first, (*rit).second.second );
        if( csmp_error.Verbose() )
        {
             if ( counter == 5 ) {
                   std::cout <<"\n\t";
                   counter=0;
               }
             counter++;
        }
    }
    if( csmp_error.Verbose() )
        std::cout << std::endl;
    return true;
}





/**

Collapses multiple regions with the same names into the first entry
so that these can be output into containers that require unique region
names.

This allows poly-element-type regions to be created that can then be
turned into CSMP groups.

The default is that elements of a lesser spatial dimension are not
incorporated into the unique region with the same type. Thus, volumetric
regions will not contain any surface elements unless the last method
argument is set to false. This largely reflects ANSYS's approach to
families and materials (used to assign properties to volumes). ANSYS
requires materials to belong to families with unique names.

@section arguments Input Arguments

The method requires a handle to the ANSYS finite element specifications
(1st argument), and the aforementioned boolean variable that determines
how duplicate region names shall be treated when they are coalesced so
that CSMP can form groups from them which must have unique names.

@return The method produces a ma map of unique regions, their element types,
and IDs of the elements that make up these regions.

@section application Application

The method is used to prepare the internally stored region information
for output to CSMP.

@section messages Messages

The method will report the names and element types of regions that were
excluded from the merged regions map.
*/

bool ModelTopology::AddRegions( const std::multimap<std::string,std::string>& object_specs,
                                const std::multimap<std::string,std::vector<size_t> >& object_elements )
 {
    if ( object_specs.empty() || object_elements.empty() )
      throw csmp::Exception( ERROR, "ModelTopology::AddRegions:", "Method did not receive any region data." );

    // 1. make a unique set of region names
    // ------------------------------------
    std::set<std::string>  unique_names;
    for ( std::multimap<std::string,std::string>::const_iterator
          it=object_specs.begin(); it!=object_specs.end(); it++ ) unique_names.insert( (*it).first );

    // 2. collapse region names and elements into new list
    // ---------------------------------------------------
    // merge all elements with the same region name together
    // --------------------------------------------------------------------
    // for each region that is unique by name
    const std::pair<std::set<std::string>,std::vector<size_t> >  empty_region;
    for ( std::set<std::string>::const_iterator
          sit=unique_names.begin(); sit!=unique_names.end(); sit++ ) {
          // get range iterators for the two multimaps with FE-names and element numbers
          std::multimap<std::string,std::string>::const_iterator eit1(object_specs.lower_bound( (*sit) ));
          std::multimap<std::string,std::string>::const_iterator eit2(object_specs.upper_bound( (*sit) ));
          std::multimap<std::string,std::vector<size_t> >::const_iterator lit1(object_elements.lower_bound( (*sit) ));
          std::multimap<std::string,std::vector<size_t> >::const_iterator lit2(object_elements.upper_bound( (*sit) ));

          // insert an empty unique new region in external 'model_regions' map
          // pair.second tests whether this region could indeed be inserted
          std::pair<std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator,bool>
          rit = model_regions.insert( std::make_pair( (*sit), empty_region ) );
          assert( rit.second == true );
          // check that object_specs and elmts for this region have the same ranges
          assert( distance( eit1, eit2 ) == distance( lit1, lit2 ) );
          // check that these ranges are not empty
          assert( distance( eit1, eit2 ) > 0 );

          // loop over the ranges, merging the element numbers and storing
          // the unique finite element types (if there is only one entry,
          // the while loop will just perform a single insertion)
          while ( eit1 != eit2 ) {
               // inserting finite element type information
               (*rit.first).second.first.insert( (*eit1).second ); // element types
               // inserting element numbers which are already unique
               (*rit.first).second.second.reserve( (*rit.first).second.second.size() + (*lit1).second.size() );
               copy( (*lit1).second.begin(), (*lit1).second.end(),
                     back_inserter((*rit.first).second.second) );
               eit1++;
               lit1++;
            }
        }

    return AddRegions( model_regions );

 } // end AddRegions





bool ModelTopology::AddRegionsWithEquidimensionalCheck( const multimap<string,string>& object_specs,
                                                        const multimap<string,vector<size_t> >& object_elements )
 {
    if ( object_specs.empty() || object_elements.empty() )
      throw csmp::Exception( WARNING, "ModelTopology::AddRegionsWithEquidimensionalCheck:", "Method did not receive any region data." );

    // 1. make a unique set of region names
    // ------------------------------------
    std::set<std::string>  unique_names;
    for ( std::multimap<std::string,std::string>::const_iterator
          it=object_specs.begin(); it!=object_specs.end(); it++ )
        unique_names.insert( (*it).first );

    // 2. collapse region names and elements into new list
    // ---------------------------------------------------
    // method enforces unique names for regions consisting of
    // volumes, surfaces, and line elements: If volume, surface and or line
    // elements belong to a family with the same name, the elements of the
    // lower spatial dimension will be ignored.
    // --------------------------------------------------------------------
    const std::pair<std::set<std::string>,std::vector<size_t> >  empty_region;
    for ( std::set<std::string>::const_iterator
          sit=unique_names.begin(); sit!=unique_names.end(); sit++ ) {
          // get range iterators for the two multimaps
          std::multimap<std::string,std::string>::const_iterator eit1(object_specs.lower_bound( (*sit) ));
          std::multimap<std::string,std::string>::const_iterator eit2(object_specs.upper_bound( (*sit) ));
          std::multimap<std::string,std::vector<size_t> >::const_iterator lit1(object_elements.lower_bound( (*sit) ));
          std::multimap<std::string,std::vector<size_t> >::const_iterator lit2(object_elements.upper_bound( (*sit) ));

          // form the unique new region
          std::pair<std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator,bool>
          rit = model_regions.insert( std::make_pair( (*sit), empty_region ) );
          assert( rit.second == true );
          // check that object_specs and elmts for this region have the same ranges
          assert( distance( eit1, eit2 ) == distance( lit1, lit2 ) );
          // check that these ranges are not empty
          assert( distance( eit1, eit2 ) > 0 );

          // getting the spatial dimension of the first region entry
          // in the region map, this dimension will be binding
          const size_t dim = fem_specs::MinimumSpatialDimension( (*eit1).second );

          // loop over the ranges, merging the element numbers and storing
          // the unique finite element types (if there is only one entry,
          // the while loop will just perform a single insertion)
          while ( eit1 != eit2 )  {
               // checking that the elements to be inserted have the
               // correct spatial dimension
               if ( fem_specs::MinimumSpatialDimension( (*eit1).second ) == dim )
                 {
                   // inserting finite element type information
                   (*rit.first).second.first.insert( (*eit1).second );
                   // inserting element numbers which are already unique
                   (*rit.first).second.second.reserve( (*rit.first).second.second.size() + (*lit1).second.size() );
                   copy( (*lit1).second.begin(), (*lit1).second.end(),
                         back_inserter((*rit.first).second.second) );
                 }
               else /*if( csmp_error.Verbose() )*/ { // if element has a different spatial dimension it is ignored;
                    std::cerr <<"\nModelTopology::AddRegionsWithEquidimensionalCheck: ";
                    std::cerr << (*lit1).second.size();
                    std::cerr <<" elements of type "<< (*eit1).second;
                    std::cerr <<" were excluded from region "<< (*sit) << std::endl;
                 }
               eit1++;
               lit1++;
            }
      }

    return AddRegions( model_regions );

 } // end AddRegionsWithEquidimensionalCheck



/**
    Method enforces unique names for regions consisting of
    volumes, surfaces, and line elements: If volume, surface and or line
    elements belong to element families with the same name, 
    the elements of lower spatial dimension will be ignored.
*/
template<size_t dim>
void ModelTopology::RemoveLowDimElementsFromRegions( csmp::VSet<dim>& vset )
 {
    size_t elmtdim;
    size_t elmtid;
    int32_t  elmttype;
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator
          rit=model_regions.begin(); rit!=model_regions.end(); rit++ )
    {
        /// find highest dimension for elements in current region
        elmtdim = 0U;
        std::set<std::string>::const_iterator etype_endit = (*rit).second.first.end();
        for( std::set<std::string>::const_iterator
             etype_it = (*rit).second.first.begin(); etype_it != etype_endit; ++etype_it )
            elmtdim = std::max( elmtdim, fem_specs::MinimumSpatialDimension( *etype_it ) );
        /// remove types of low dim elements
        std::vector<std::string> remove_types;
        for( std::set<std::string>::const_iterator
             etype_it = (*rit).second.first.begin(); etype_it != etype_endit; ++etype_it )
            if ( fem_specs::MinimumSpatialDimension( *etype_it ) != elmtdim )
                remove_types.push_back( *etype_it );
        /// remove ids of low dim elements
        std::vector<size_t> remove_eids;
        for ( size_t i = 0 ; i < (*rit).second.second.size(); ++i )
        {
            elmtid   = (*rit).second.second[ i ];
            elmttype = vset.ElementType( elmtid );
            // if element has a different spatial dimension it will be removed
            if ( fem_specs::MinimumSpatialDimension( static_cast<int8_t>(elmttype) ) != elmtdim )
                remove_eids.push_back( i );
        }
        if( !remove_eids.empty() || !remove_types.empty() )
        {
            //if( csmp_error.Verbose() ){
            std::cerr <<"\nModelTopology::RemoveLowDimElementsFromRegions: ";
            std::cerr << remove_eids.size();
            std::cerr <<" elements of types: ";
            for( size_t i = 0 ; i< remove_types.size(); ++i ){
                if( i != 0 ) std::cerr<< ", ";
                std::cerr << remove_types[ i ];
            }
            std::cerr <<" were excluded from region "<< (*rit).first << std::endl;
            //}
            for( size_t i = 0 ; i< remove_types.size(); ++i )
                (*rit).second.first.erase( remove_types[i] );
            for( size_t i = 0 ; i < remove_eids.size(); ++i )
                (*rit).second.second.erase( (*rit).second.second.begin() + ( remove_eids[i] - i ) );
        }
    }
 } // end RemoveLowDimElementsFromRegions

template void ModelTopology::RemoveLowDimElementsFromRegions( csmp::VSet<1U>& );
template void ModelTopology::RemoveLowDimElementsFromRegions( csmp::VSet<2U>& );
template void ModelTopology::RemoveLowDimElementsFromRegions( csmp::VSet<3U>& );


/**
	return names of regions in the model
*/
void ModelTopology::RegionNames( std::vector<std::string>& region_names ) const
{
	for ( map<std::string, std::pair<std::set<std::string>, std::vector<size_t> > >::const_iterator
		    it = model_regions.begin(); it != model_regions.end(); ++ it )
		region_names.push_back( (*it).first );
}

/*
    Region properties
*/

/**

Read material properties from '*-regions.txt' file and associate them with
the regions storage. The property values are going to be constant throughout
the regions to which they are assigned.

@section arguments Input Arguments

The method expects that the file name has the appendage and extension
'-regions.txt'.
*/
void  ModelTopology::PropertiesOfRegions( const char* regions_file,
                                          std::list<std::string>& properties,
                                          std::map<std::string,std::list<double> >& props ) const
 {
    char               text_line[500];
    char*              token(0);
    const char* const  delims ="\t,:,\n,\r";
    const char* const  delims2 =" ,\t,\n,\r";
    strcpy( text_line, regions_file );
    strcat( text_line, "-regions.txt" );
    std::ifstream        ifs(text_line);
    std::list<double>  prop_vals;
    double             val;
    std::string          region;

    if ( !ifs.is_open() )
      throw csmp::Exception( FATAL_ERROR, "ModelTopology::PropertiesOfRegions:", text_line,
                                      "ASCII geometry input file could not be opened");

    if ( !props.empty() )
      throw csmp::Exception( INFO, "ModelTopology::PropertiesOfRegions:", regions_file,
                            "Supplied non-empty regions map is erased");

    // 1. reading and discarding file header
    ifs.getline( text_line, 500 );

    // 2. reading the names of the properties
    ifs.getline( text_line, 500 );
    properties.push_back( strtok( text_line, delims ) );

    while ( (token=strtok( NULL, delims )) != NULL )
        properties.push_back( token );

    // 3. reading region descriptors and associated properties from file
    while ( !ifs.eof() )
      {
         ifs.getline( text_line, 256 );
         if ( strlen(text_line) == 0 ) break;
         // region name
         region = strtok( text_line, delims2 );
         // if the region is part of the current model its properties are read
         if ( Contains(region.c_str()) ) {
               for ( size_t i=0U; i<properties.size(); i++ ) {
                    val = atof( strtok( NULL, delims2 ) );
                    prop_vals.push_back( val );
                 }
               props[ region ] = prop_vals;
               prop_vals.erase( prop_vals.begin(), prop_vals.end() );
           }
         else
         throw csmp::Exception( ERROR, "\nModelTopology::PropertiesOfRegions",
                           "Region file contains unrecognized region, e.g.", region.c_str() );
      }
    ifs.close();

    // 4. screen output of read properties
    std::cout <<"\nModelTopology::PropertiesOfRegions: Input properties read from text file '";
    std::cout << regions_file <<"':\n\n\t";
    for ( std::list<std::string>::const_iterator
          pit=properties.begin(); pit!=properties.end(); pit++ ) std::cout << (*pit) <<" ";
    std::cout <<"\n\nValues of these properties for listed regions: "<< std::endl;
    for ( std::map<std::string,std::list<double> >::const_iterator
          it=props.begin(); it!=props.end(); it++ ) {
         std::cout <<"\t"<< (*it).first <<": ";
         for ( std::list<double>::const_iterator
               dit=(*it).second.begin(); dit!=(*it).second.end(); dit++ )
           std::cout << (*dit) <<"  ";
         std::cout << std::endl;
      }

 } // end PropertiesOfRegions




/**

AssignMaterialProperties() allows the user to interactively assign
properties to ANSYS models via the console. This is tedious and we
recommend to use a configuration file for this purpose.

The only properties that can be assigned are element properties, ie,
piecewise constant material properties.

The method adds the property data to the argument VSet.

@section arguments Input Arguments

The methods takes a const reference to the model topology object.

@section messages Messages

The user is prompted for region and property names.

*/
template<size_t dim>
void ModelTopology::AssignMaterialProperties( VSet<dim>& vset,
                                              const std::multimap<std::string,std::vector<size_t> >& object_elements )
 {
    std::string     prop_name;
    double        prop_val;
    std::set<std::string> box_boundaries;

    BoundariesOfBoxShapedModel( box_boundaries );

    std::cout <<"\nANSYS_Interface::AssignMaterialProperties: "<< std::endl;
    //  for the permeability (and the volume elements in the model)
    PropertyData  data( ELEMENT, SCALAR, dim );
    data.Resize( vset.Elements(), vset.Elements() );
   
    for ( std::multimap<std::string,std::vector<size_t> >::const_iterator
          it=object_elements.begin(); it!=object_elements.end(); it++ )
      {
          if ( MinimumSpatialDimensionOfRegion( ((*it).first.c_str()) ) ) {
               std::cout <<"\n\tEnter 'permeability' value for geometric object '"<< (*it).first <<"': ";
               std::cin  >> prop_val;
               // for each element
               for ( std::vector<size_t>::const_iterator
                     lit=(*it).second.begin(); lit!=(*it).second.end(); lit++ )
                 // data[ (*lit) ] = prop_val;
                 data.Value( (*lit), 0U ) = prop_val;
            }
         // putting the data into the vset
         vset.AddData( "permeability", data );
      }

    std::cout <<"\n\tEnter for how many (scalar) properties you would like to assign values to regions: ";
    int32_t assignments;
    std::cin >> assignments;
    if ( assignments == 0 ) return;

    // for each additional scalar property
    PropertyData  scdata( ELEMENT, SCALAR, dim );
    data.Resize( vset.Elements(), vset.Elements() );

    for ( int32_t i=0; i<assignments; i++ ) {
         std::cout <<"\n\tEnter property name: ";
         std::cin  >> prop_name;

         // for each object (except for boundaries and edges in box-shaped model)
         for ( std::multimap<std::string,std::vector<size_t> >::const_iterator
               it=object_elements.begin(); it!=object_elements.end(); it++ ) {
              // if it is not one of the standard model boundaries or edges
              if ( (box_boundaries.find((*it).first)) == box_boundaries.end() ) {
                   std::cout <<"\n\tEnter property value for geometric object '"<< (*it).first <<"': ";
                   std::cin  >> prop_val;
                   // for each element
                   for ( std::vector<size_t>::const_iterator
                         lit=(*it).second.begin(); lit!=(*it).second.end(); lit++ )
                     scdata.Value( (*lit), 0U ) = prop_val;
                }
           }
         // putting the data into the vset
         vset.AddData( prop_name.c_str(), scdata );
      }
 } // end AssignMaterialProperties


template void ModelTopology::AssignMaterialProperties( VSet<1U>&,const std::multimap<std::string,std::vector<size_t> >&);
template void ModelTopology::AssignMaterialProperties( VSet<2U>&,const std::multimap<std::string,std::vector<size_t> >&);
template void ModelTopology::AssignMaterialProperties( VSet<3U>&,const std::multimap<std::string,std::vector<size_t> >&);






/*
    Numbering
 */

/**
 
Tests whether the element numbers stored in ModelTopology are in a format
suitable for reorganising the VSet. The following tests are made:

- does the element numbering start with zero?
- is the largest element number equivalent to the number of elements-1?
- are the elements numbered consecutively?
- are the elemnt ids' within each Region consecutive.
- are there any duplicate elements?
 
The method assumes that that each element can only belong to a single
region.

@return Boolean. If any of the tests fails the method returns false.

@section messages Messages

The method complains if the element numbering does not start at 0
or the largest element number is not equal to nn-elements-1.  

*/
bool  ModelTopology::CheckElementNumbering() const
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    std::set<size_t> element_ids;

    // all numbers of elements from all the current regions are inserted into one set 
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          rit=model_regions.begin(); rit!=model_regions.end(); rit++ ) 
      {
         std::vector<size_t>::const_iterator litp1(++((*rit).second.second.begin()));
         const std::vector<size_t>::const_iterator litEnd( (*rit).second.second.end() );
         for ( std::vector<size_t>::const_iterator
               lit=(*rit).second.second.begin(); lit!=litEnd; lit++ ) 
           {
              // sequence is not consecutive if the next element has not got a number that is en+1
              if ( litp1 != (*rit).second.second.end() and (*lit+1U) != *litp1 )
              {
                  if( csmp_error.Verbose() )
                      csmp_error.notice( WARNING, "ModelTopology::CheckElementNumbering:", "Sequence of elemnt id's within the Region is not consecutive." );
                  return false;
              }
              // sequence is not consecutive if next element number cannot be inserted into it because it is non-unique
              std::pair<std::set<size_t>::iterator,bool> it=element_ids.insert(*lit);
              if ( it.second == false )
              {
                  if( csmp_error.Verbose() )
                      csmp_error.notice( WARNING, "ModelTopology::CheckElementNumbering:", "Sequence of elemnt id's within the Region is not consecutive." );
                  return false;
              }
              if( litp1 != litEnd ) ++litp1;
           }
      }

    // if the first element is not numbered zero
    if ( (*element_ids.begin()) != 0U ) {
          if( csmp_error.Verbose() ) {
                std::string  err_msg("first element number ");
                err_msg += to_string( (*max_element(element_ids.begin(),element_ids.end())) );
                err_msg +=" is not equal to zero.";
                csmp_error.notice( WARNING, "ModelTopology::CheckElementNumbering:", err_msg.c_str() );
            }
          return false;
      }
      
    // if the last element number is not equivalent to the total number of contained elements - 1
    if ( (*max_element(element_ids.begin(),element_ids.end())) != element_ids.size()-1U ) {
          if ( csmp_error.Verbose() ) {
                std::string  err_msg("largest element number ");
                err_msg += std::to_string( (*max_element(element_ids.begin(),element_ids.end())) );
                err_msg +="-1 is not equal to the total number of elements ";
                err_msg += std::to_string( element_ids.size() );
                csmp_error.notice( ERROR, "ModelTopology::CheckElementNumbering:", err_msg.c_str() );
            }
          return false;
      }
     
    return true;
 
 } // end CheckElementNumbering
 




/**
     Creates a new contiguous element numbering 0..n-1 and outputs an old-to-new mapping into its argument map.
     
     @param check_output - if true a consistency check is performed on the new mapping, using
     the ConsecutiveSequenceChecker.
*/
void  ModelTopology::CreateNewElementNumbers( std::map<size_t /* old-# */,size_t /* new-# */>& eid_mapping, bool check_output )
{
    // 1. assuming that the elements are already numbered correctly, this numbering only is output to the map
    if ( !eid_mapping.empty() ) eid_mapping.clear();
  
    size_t  eid(0U);
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          rit=model_regions.begin(); rit!=model_regions.end(); rit++ )
      // for all element IDs of each region
      for ( std::vector<size_t>::const_iterator
            lit=(*rit).second.second.begin(); lit!=(*rit).second.second.end(); ++lit )
        {
           // only if the element ID could be inserted the element counter is incremented            old  new
           std::pair<std::map<size_t,size_t>::iterator,bool> it=eid_mapping.insert( std::make_pair( *lit, eid ) );
           if ( it.second == true ) eid++;
        }
  
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    if ( csmp_error.Verbose() ) {
         std::cout <<"\nModelTopology::CreateNewElementNumbers: mapped "<< eid;
         std::cout <<" elements successfully to consecutive numbers."<< std::endl;
         std::cout.flush();
      }

    if ( check_output ) {
    // ConsecutiveSequenceChecker::Test_ConsecutiveSequenceChecker();
         const bool check_whether_max_value_is_size_minus1(true);
         if ( !ConsecutiveSequenceChecker::IsValueRangeUniqueAndBounded( eid_mapping, check_whether_max_value_is_size_minus1 ) )
           csmp_error.notice( WARNING, "ModelTopology::CreateNewElementNumbers:",
                            "the renumbered element range is not consecutive and unique; trying to fix this.");
      }

    // needs to be done only if something changed
    RenumberElements( eid_mapping );
  
} // end CreateNewElementNumbers







/**
 
The elements in each region are renumbered using the supplied number mapping. 
In debug mode, first a check is performed whether this numbering
is consecutive, starts with 0 and there are no duplicates.  

@param eid_mapping Map with ID correspondance between old and new IDs.

@section messages Messages

The method throws CSMP Exceptions.  
*/
void  ModelTopology::RenumberElements( const std::map<size_t,size_t>& eid_mapping )
 {
    //std::cerr <<"\nelement id mapping: size: "<< eid_mapping.size() <<":\n";
    //for ( auto it=eid_mapping.begin(); it!=eid_mapping.end(); ++it )
    //std::cerr << (*it).first <<"->"<< (*it).second <<" ";

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( eid_mapping.empty() ) {
         csmp_error.notice( WARNING, "ModelTopology::RenumberElements",
                           "The element number correspondance map is empty. Nothing was done.");
         return;
      }

#ifndef NDEBUG      
    std::set<size_t>  new_eids;
    
    // 0. checking whether the elements are consecutively numbered
    //    (these numbers are stored as values eid_map)
    //    making a set of the new element ID numbers and recording the old max element ID for 
    //    later comparison with new one
    size_t  new_min_id(NULL_IDX), new_max_id(0U);
    for ( std::map<size_t,size_t>::const_iterator
          it=eid_mapping.begin(); it!=eid_mapping.end(); it++ ) {
         new_eids.insert( (*it).second );   
         new_min_id = std::min( new_min_id, (*it).second );
         new_max_id = std::max( new_max_id, (*it).second );
      }

    if ( new_min_id != 0U )
      throw csmp::Exception( WARNING, "ModelTopology::RenumberElements",
                            "The new number range does not commence with zero. Nothing was done.");
                               
    if ( new_max_id != new_eids.size()-1U)
      throw csmp::Exception( WARNING, "ModelTopology::RenumberElements",
                            "The supplied maximum element number is not equal to N-elements-1. Nothing was done.");

    // is input range consecutive
    std::set<size_t>::const_iterator it2(++(new_eids.begin()));
    for ( std::set<size_t>::const_iterator it1=new_eids.begin(); it2!=new_eids.end(); it1++, it2++ )
      if ( (*it2) != ((*it1)+1U) ) 
        throw csmp::Exception( WARNING, "ModelTopology::RenumberElements",
                              "The supplied new element numbers are not consecutive. Nothing was done.");
#endif

    // 2. renumbering the elements
    std::vector<size_t>  new_region_eids;
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::iterator
          rit=model_regions.begin(); rit!=model_regions.end(); ++rit )
      {
         const size_t elements_to_renumber((*rit).second.second.size());
         // for each region
         new_region_eids.reserve(elements_to_renumber);
         // for all the element IDs of the region 
         for ( std::vector<size_t>::const_iterator
               lit=(*rit).second.second.begin(); lit!=(*rit).second.second.end(); ++lit ) {
              std::map<size_t,size_t>::const_iterator it=eid_mapping.find(*lit);
              // checking that the element ID was found
              assert( it != eid_mapping.end() );
              new_region_eids.push_back( (*it).second );
           }
         // writing the renumbered vector ModelTopology
         assert( new_region_eids.size() == elements_to_renumber );
         (*rit).second.second = new_region_eids;
         new_region_eids.clear();
      }
    if( csmp_error.Verbose() )
        std::cout <<"\nModelTopology::RenumberElements: Successfully re-established consecutive element number range."<< std::endl;
        
 } // end RenumberElements



/**
    Renumbers the element ids stored in ModelTopology and eliminates potentiallu unused
    elements from the argument VSet.
    
    @param check_range allows user to check the newly generated range again.
*/
template<size_t dim>
void  ModelTopology::RenumberElements(csmp::VSet<dim>& vset, bool check_range )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );

   std::map<size_t,size_t>  old_and_new_elmtids;
   CreateNewElementNumbers( old_and_new_elmtids, false );
   vset.ReduceTo( old_and_new_elmtids );
   old_and_new_elmtids.clear();

  if ( check_range ) {
        const bool check_whether_max_value_is_size_minus1(true);
        if ( !ConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive( old_and_new_elmtids, check_whether_max_value_is_size_minus1 ) )
          csmp_error.notice( ERROR, "ModelTopology::RenumberElements:", "failed to calculate consecutive new element idx range.");
    }
}

template void ModelTopology::RenumberElements( csmp::VSet<1U>&,bool );
template void ModelTopology::RenumberElements( csmp::VSet<2U>&,bool );
template void ModelTopology::RenumberElements( csmp::VSet<3U>&,bool );





/**
    Checks that the new numbers in the VSet form a consecutive range.
    
    test: assert( ConsecutiveSequenceChecker::Test_ConsecutiveSequenceChecker() );

*/
template<size_t dim>
bool ModelTopology::EstablishTopology( VSet<dim>& vset,
                                       const multimap<string,string>& object_specs,
                                       const multimap<string,vector<size_t> >& object_elements,
                                       bool require_unique_names_for_vol_surf_lines,
                                       bool interactive_property_assignment,
                                       bool correct_orientation_of_surface_elements,
                                       bool reassign_boundary_flags )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // 1. eliminating lower-dimensional regions that have the same name as higher-dimensional ones
    // -------------------------------------------------------------------------------------------
    // the chosen regions get added to the ModelTopology object
    if ( require_unique_names_for_vol_surf_lines )
        AddRegionsWithEquidimensionalCheck( object_specs, object_elements );
    else
        AddRegions( object_specs, object_elements );

    // 2. check element numbering
    // -------------------------------------------------------
    // checking that the new numbers form a consecutive range
    // test: assert( ConsecutiveSequenceChecker::Test_ConsecutiveSequenceChecker() );
    const bool check_whether_max_value_is_size_minus1(true);
    if ( !ConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive( object_elements, check_whether_max_value_is_size_minus1 ) )
      RenumberElements( vset, false );

    // 3. assign boundary flags to box-shaped model
    // --------------------------------------------
    if ( reassign_boundary_flags ) {
        if ( !AssignBoxShapedModelFlags(vset) )
          csmp_error.notice( WARNING, "ModelTopology::EstablishTopology:",
                            "although this claims to be a box-shaped model, a correct BOX_BOUNDARY flagging could not be established." );
     }
     
    // 4. correct surface mesh orientation
    // ---------------------------------------
    // SKM note: this functionality also deals with surface mesh in 3D meshes 
    if constexpr ( dim == 2U ) {
        if ( correct_orientation_of_surface_elements )
          if ( InterpolationOrder() <= 2 )
            // NB! - actually makes changes for most 2D ANSYS models
            vset.RenumberElementsCounterClockwise2D();
      }

    // 5. assigning properties to regions
    // ----------------------------------
    if ( interactive_property_assignment )
        AssignMaterialProperties( vset, object_elements );

    return true;
}

template bool ModelTopology::EstablishTopology( VSet<1U>&,const std::multimap<std::string,std::string>&,const std::multimap<std::string,std::vector<size_t> >&,bool,bool,bool,bool);
template bool ModelTopology::EstablishTopology( VSet<2U>&,const std::multimap<std::string,std::string>&,const std::multimap<std::string,std::vector<size_t> >&,bool,bool,bool,bool);
template bool ModelTopology::EstablishTopology( VSet<3U>&,const std::multimap<std::string,std::string>&,const std::multimap<std::string,std::vector<size_t> >&,bool,bool,bool,bool);





template<size_t dim>
bool ModelTopology::EstablishTopology( VSet<dim>& vset,
                                   bool require_unique_names_for_vol_surf_lines,
                                   bool correct_orientation_of_surface_elements,
                                   bool reassign_boundary_flags )
{
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    bool checks_passed(true);
  
    // 1. merge regions
    // -------------------------------------------------------
    if ( require_unique_names_for_vol_surf_lines )
        RemoveLowDimElementsFromRegions( vset );

    // 2. check numbering of keys and values in map
    // -------------------------------------------------------
    std::map<size_t,size_t>  old_and_new_elmtids;
    CreateNewElementNumbers( old_and_new_elmtids, false );
    bool check_whether_max_value_is_size_minus1( true );
    if ( !ConsecutiveSequenceChecker::IsKeyRangeOfUnsignedIntConsecutive( old_and_new_elmtids, check_whether_max_value_is_size_minus1 ) ) {
         csmp_error.notice( WARNING, "ModelTopology::EstablishTopology:", "input element number range is not consecutive.");
         // looking at the input  range
         std::cerr <<"\n\tfirst element-Idx stored in model topology: "<< (*old_and_new_elmtids.begin()).first;
         std::cerr <<"\n\tlast element-Idx stored in model topology: "<< (*old_and_new_elmtids.rbegin()).first <<"\n";
         // printing the new element range
         std::set<size_t>  range;
         for ( auto it=old_and_new_elmtids.begin(); it!=old_and_new_elmtids.end(); ++it ) range.insert( (*it).first );
         std::cerr <<"\n\telement Idx range (only printing the non-consective element numbers in the sequence):";
         size_t last_value(0U);
         for ( auto sit : range ) {
              if ( sit > 0U and sit != last_value+1U ) {
                   std::cerr <<"\n\t\t"<< last_value <<" sequence break "<< sit;
                }
              last_value = sit;
           }
         std::cerr <<"\n\n";
         checks_passed = false;
      }

    if ( !ConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive( old_and_new_elmtids, check_whether_max_value_is_size_minus1 ) ) {
         csmp_error.notice( WARNING, "ModelTopology::EstablishTopology:", "output (new) element number range is not consecutive.");
         RenumberElements( vset, check_whether_max_value_is_size_minus1=false ); // false=done already
         checks_passed = false;
      }
 
    // 3. assign boundary flags for box-shaped model
    // -----------------------------------------------
    if ( reassign_boundary_flags ) {
         if ( !AssignBoxShapedModelFlags(vset) )
           csmp_error.notice( ERROR, "ModelTopology::EstablishTopology:",
                             "although this claims to be a box-shaped model, a correct BOX_BOUNDARY flagging could not be established." );
      }
      
    // 4. correct surface mesh orientation
    // ---------------------------------------
    if constexpr ( dim == 2U )
      if ( correct_orientation_of_surface_elements )
        vset.RenumberElementsCounterClockwise2D();
      
    return checks_passed;
}

template bool ModelTopology::EstablishTopology( VSet<1U>&,bool,bool,bool);
template bool ModelTopology::EstablishTopology( VSet<2U>&,bool,bool,bool);
template bool ModelTopology::EstablishTopology( VSet<3U>&,bool,bool,bool);












/*
    Box Shaped Model related
*/


/**
A box-shaped model has walls that lie in the XY, XZ, and YZ planes of
the coordinate system. These wall must be identified by groups of
surface elements with the correct CSP names LEFT, RIGHT, BOTTOM,
TOP, FRONT and BACK. If these boundary identifying surfaces are present
the topology is interpreted as box shaped. A model can also be 
box shaped if the top boundary is not present but a surface exists
in its place that is called IRREGULAR.  

@return If the model contains the aforementioned boundary surfaces 'true' is
returned, else 'false'.  
*/
bool  ModelTopology::BoxShapedModel() const
 {
    std::map<std::string,bool>            box_boundaries;
    std::map<std::string,bool>::iterator  bit;
    
    // boundaries
    box_boundaries.insert( std::make_pair(std::string("BOTTOM"),false) );
    box_boundaries.insert( std::make_pair(std::string("LEFT"),false) );
    box_boundaries.insert( std::make_pair(std::string("RIGHT"),false) );
    box_boundaries.insert( std::make_pair(std::string("TOP"),false) );
    box_boundaries.insert( std::make_pair(std::string("FRONT"),false) );
    box_boundaries.insert( std::make_pair(std::string("BACK"),false) );

    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); ++it )
      if ( (bit=box_boundaries.find((*it).first)) != box_boundaries.end() )
        (*bit).second = true;
        
    // checking whether all boundary regions were found in the dataset
    for (  bit=box_boundaries.begin(); bit!=box_boundaries.end(); bit++ )
      if ( (*bit).second == false ) return false;
      
    return true;

 } // end BoxShapedModel


 bool  ModelTopology::RectangleShapedModel() const
 {
    // making a boundary map
    std::map<std::string,bool>  box_boundaries;
    box_boundaries.insert( std::make_pair(std::string("BOTTOM"),false) );
    box_boundaries.insert( std::make_pair(std::string("LEFT"),false) );
    box_boundaries.insert( std::make_pair(std::string("RIGHT"),false) );
    box_boundaries.insert( std::make_pair(std::string("TOP"),false) );
    
    std::map<std::string,bool>::iterator  bit;
    for ( std::map<std::string,std::pair<std::set<std::string>,std::vector<size_t> > >::const_iterator
          it=model_regions.begin(); it!=model_regions.end(); it++ )
      if ( (bit=box_boundaries.find((*it).first)) != box_boundaries.end() )
        (*bit).second = true;
        
    // checking whether all boundaries were found in the dataset
    for (  bit=box_boundaries.begin(); bit!=box_boundaries.end(); bit++ )
      if ( (*bit).second == false ) return false;
      
    return true;

 } // end RectangleShapedModel




/**
 
 Writes the standard boundary name std::string objects
 into the supplied set.

 An empty set of std::strings that will be initialized with the boundary
 names. The set is emptied before the data are stored if it already
 contains std::strings.
 */
 void BoundariesOfBoxShapedModel( std::set<std::string>& bs )
  {
     if ( !bs.empty() )
         bs.erase( bs.begin(), bs.end() );
     bs.insert("BOTTOM");
     bs.insert("LEFT");
     bs.insert("RIGHT");
     bs.insert("TOP");
     bs.insert("FRONT");
     bs.insert("BACK");
     bs.insert("BACK_BOTTOM");
     bs.insert("BACK_RIGHT");
     bs.insert("BACK_TOP");
     bs.insert("BACK_LEFT");
     bs.insert("BOTTOM_LEFT");
     bs.insert("BOTTOM_RIGHT");
     bs.insert("TOP_RIGHT");
     bs.insert("TOP_LEFT");
     bs.insert("FRONT_BOTTOM");
     bs.insert("FRONT_RIGHT");
     bs.insert("FRONT_TOP");
     bs.insert("FRONT_LEFT");
  }

 void BoundariesOfRectangleShapedModel( std::set<std::string>& bs )
  {
     if ( !bs.empty() ) bs.erase( bs.begin(), bs.end() );
     bs.insert("BOTTOM");
     bs.insert("LEFT");
     bs.insert("RIGHT");
     bs.insert("TOP");
  }



/**
    Overwrites all pre-existing boundary flags.
*/
 template<size_t dim>
 bool ModelTopology::AssignBoxShapedModelFlags( VSet<dim>& vset ) const
 {
     ErrorHandler& csmp_error ( ErrorHandler::Instance() );
     bool flagging_was_done_successfully(false);

     // flag boundaries if the model is box-shaped ( 3D )
     // ------------------------------------------------------
     if constexpr ( dim == 3U ) {
         if ( BoxShapedModel() )
           {
              if ( csmp_error.Verbose() )
                {
                   std::cout <<"\n\nModelTopology::AssignBoxShapedModelFlags: The model is 'box shaped'. ";
                   std::cout <<"Re-assigning boundary flags to box-shaped model..."<< std::endl;
                }
               // flagging the boundary nodes according to CSMP specs
               flagging_was_done_successfully = FlagNodesUsingBoundaryRegions( vset );
           }
       }
     // flag boundaries if the model is rectangle-shaped ( 2D )
     // -------------------------------------------------------
     else if constexpr ( dim == 2U ) {
           if ( RectangleShapedModel() )
             {
                 if ( csmp_error.Verbose() ) {
                      std::cout <<"\n\nModelTopology::AssignBoxShapedModelFlags: The model is 'rectangle shaped'. ";
                      std::cout <<"Building neighbor connectivity and assigning boundary flags in rectangle-shaped model..."<< std::endl;
                   }
                 // for all surface elements in the mesh which have bar element neighbors,
                 // assign appropriate boundary flags
                 flagging_was_done_successfully = FlagNodesUsingBoundaryRegions( vset );
             }
            return flagging_was_done_successfully;
         }
         
    return flagging_was_done_successfully;
    
 } // end AssignBoxShapedModelFlags

 template bool ModelTopology::AssignBoxShapedModelFlags( VSet<2U>& ) const;
 template bool ModelTopology::AssignBoxShapedModelFlags( VSet<3U>& ) const;





/**
       Deduces boundary flags from region names, including the basic BOX_BOUNDARY flags,
       All other boundary regions (regionts that contain the string "BOUNDARY" are flagged IRREGULAR.
       
       @author SKM
       @date 17/10/2021
*/
bool ModelTopology::FlagNodesUsingBoundaryRegions( VSet<2U>& vset ) const
  {
     set<string> boundary_regions;

     // 0. Extracting the names of valid BOUNDARY regions from file
     // -----------------------------------------------------------
     //  region name  etypes-of-region  elmt-ids for region
     //   map<string, pair<set<string>,vector<size_t> > >  model_regions;
     for ( auto& it : model_regions ) {
          const string region_name{ it.first };
          if ( region_name.find("BOUNDARY") != std::string::npos ||
               region_name.find("IRREGULAR") != std::string::npos ) {
               // if it is a lower dimensional region
               const auto etype = it.second.first.begin(); // .first.begin());
               if ( isLineElement( parseFiniteElementType(*etype) ) == true )
                 boundary_regions.insert( region_name );
            }
       }
     
     // 1. Making map of the target regions in which to search for the boundary elements
     deque<size_t>  nodes_left, nodes_right, nodes_bottom, nodes_top, nodes_irregular;
     std::vector<size_t>::const_iterator  bit;

     for ( bit=ElementsOfRegionBegin("BOTTOM");
           bit!=ElementsOfRegionEnd("BOTTOM"); bit++ )
       for ( vector<int64_t>::iterator nit=vset.PlistBegin(*bit); nit!=vset.PlistEnd(*bit); ++nit )
         nodes_bottom.push_back( (*nit) );
       
     for ( bit=ElementsOfRegionBegin("RIGHT");
           bit!=ElementsOfRegionEnd("RIGHT"); bit++ )
       for ( vector<int64_t>::iterator nit=vset.PlistBegin(*bit); nit!=vset.PlistEnd(*bit); ++nit )
         nodes_right.push_back( (*nit) );
       
     for ( bit=ElementsOfRegionBegin("TOP");
           bit!=ElementsOfRegionEnd("TOP"); bit++ )
       for ( vector<int64_t>::iterator nit=vset.PlistBegin(*bit); nit!=vset.PlistEnd(*bit); ++nit )
         nodes_top.push_back( (*nit) );
       
     for ( bit=ElementsOfRegionBegin("LEFT");
           bit!=ElementsOfRegionEnd("LEFT"); bit++ )
       for ( vector<int64_t>::iterator nit=vset.PlistBegin(*bit); nit!=vset.PlistEnd(*bit); ++nit )
         nodes_left.push_back( (*nit) );
 
      for ( auto& it : boundary_regions )
       for ( bit=ElementsOfRegionBegin(it.c_str());
             bit!=ElementsOfRegionEnd(it.c_str()); bit++ )
         for ( vector<int64_t>::iterator nit=vset.PlistBegin(*bit); nit!=vset.PlistEnd(*bit); ++nit )
           nodes_irregular.push_back( (*bit) );

     // making these containers unique
     sort( nodes_bottom.begin(), nodes_bottom.end() );
     nodes_bottom.erase( unique( nodes_bottom.begin(), nodes_bottom.end() ), nodes_bottom.end() );

     sort( nodes_right.begin(), nodes_right.end() );
     nodes_right.erase( unique( nodes_right.begin(), nodes_right.end() ), nodes_right.end() );

     sort( nodes_top.begin(), nodes_top.end() );
     nodes_top.erase( unique( nodes_top.begin(), nodes_top.end() ), nodes_top.end() );

     sort( nodes_left.begin(), nodes_left.end() );
     nodes_left.erase( unique( nodes_left.begin(), nodes_left.end() ), nodes_left.end() );

     sort( nodes_irregular.begin(), nodes_irregular.end() );
     nodes_irregular.erase( unique( nodes_irregular.begin(), nodes_irregular.end() ), nodes_irregular.end() );


     // 2. Going over the nodes, making BOX_BOUNDARY flag assignments
     // -------------------------------------------------------------
     // zapping all previous box boundary flags
     for ( auto bit=vset.BFlagsBegin(); bit!= vset.BFlagsEnd(); ++bit ) (*bit) = NOT;
     // (starting with the least specific regions so that their corners are overwritten by box boundary flags)
     for ( auto nit : nodes_irregular ) vset.AddBFlag( nit, IRREGULAR );
     for ( auto nit : nodes_bottom ) vset.AddBFlag( nit, BOTTOM );
     for ( auto nit : nodes_right ) vset.AddBFlag( nit, RIGHT );
     for ( auto nit : nodes_top ) vset.AddBFlag( nit, TOP );
     for ( auto nit : nodes_left ) vset.AddBFlag( nit, LEFT );


     // 3. finding the corners of a box-shaped model, if any
     // ----------------------------------------------------
     if ( !nodes_left.empty() && !nodes_bottom.empty() ) {
          vector<int64_t> cnr;
          set_intersection( nodes_left.begin(), nodes_left.end(),
                            nodes_bottom.begin(), nodes_bottom.end(),
                            back_inserter( cnr ) );
         
          assert( !cnr.empty() );
          vset.AddBFlag( cnr[0], CNR1 );
       }

     if ( !nodes_right.empty() && !nodes_bottom.empty() ) {
          vector<int64_t> cnr;
          set_intersection( nodes_right.begin(), nodes_right.end(),
                            nodes_bottom.begin(), nodes_bottom.end(),
                            back_inserter( cnr ) );
         
          assert( !cnr.empty() );
          vset.AddBFlag( cnr[0], CNR2 );
       }

     if ( !nodes_right.empty() && !nodes_top.empty() ) {
          vector<int64_t> cnr;
          set_intersection( nodes_right.begin(), nodes_right.end(),
                            nodes_top.begin(), nodes_top.end(),
                            back_inserter( cnr ) );
         
          assert( !cnr.empty() );
          vset.AddBFlag( cnr[0], CNR3 );
       }

     if ( !nodes_left.empty() && !nodes_top.empty() ) {
          vector<int64_t> cnr;
          set_intersection( nodes_left.begin(), nodes_left.end(),
                            nodes_top.begin(), nodes_top.end(),
                            back_inserter( cnr ) );
         
          assert( !cnr.empty() );
          vset.AddBFlag( cnr[0], CNR4 );
       }

     // 4. potential corners with an irregular boundary
     // -----------------------------------------------

     // if there is an irregular boundary on the model top instead of TOP
     if ( !nodes_irregular.empty() && nodes_top.empty() ) {
         if ( !nodes_left.empty() ) {
              vector<int64_t> cnr;
              set_intersection( nodes_irregular.begin(), nodes_irregular.end(),
                                nodes_left.begin(), nodes_left.end(),
                                back_inserter( cnr ) );
             
              vset.AddBFlag( cnr[0], CNR4 );
           }
         if ( !nodes_right.empty() ) {
              vector<int64_t> cnr;
              set_intersection( nodes_irregular.begin(), nodes_irregular.end(),
                                nodes_right.begin(), nodes_right.end(),
                                back_inserter( cnr ) );
             
              vset.AddBFlag( cnr[0], CNR3 );
           }
       }

     // if there is an irregular boundary on the model bottom instead of BOTTOM
     if ( !nodes_irregular.empty() && nodes_bottom.empty() ) {
         if ( !nodes_left.empty() ) {
              vector<int64_t> cnr;
              set_intersection( nodes_irregular.begin(), nodes_irregular.end(),
                                nodes_left.begin(), nodes_left.end(),
                                back_inserter( cnr ) );
             
              vset.AddBFlag( cnr[0], CNR1 );
           }
         if ( !nodes_right.empty() ) {
              vector<int64_t> cnr;
              set_intersection( nodes_irregular.begin(), nodes_irregular.end(),
                                nodes_right.begin(), nodes_right.end(),
                                back_inserter( cnr ) );
             
              vset.AddBFlag( cnr[0], CNR2 );
           }
       }

     return true;

  } // end FlagNodesUsingBoundaryRegions(2D)





 



 /**

 Infer_BOX_BOUNDARY_EdgesAndCornersFlagsFromSideFlags() 
 looks for families of surface
 elements in the supplied model topology which are named FRONT, BACK,
 LEFT, RIGHT, TOP and BOTTOM. If these can be found the model is accepted
 as box shaped and the missing edge and corner point boundary flags
 are assigned to the input VSet in addition of the ones for the nodes
 at the model surface.

 @section arguments Input Arguments

 The method needs a reference to the model topology and the
 corresponding VSet.

 @return Nodes in the argument VSet will be assigned the correct boundary flags
 if the model is genuinely box-shaped and the correct surface family
 names have been assigned.

 @section implementation Implementation

 Since the method intersects the bounding surfaces of the box-shaped model
 it gets confused if these have not been labeled according to CSMP
 conventions which require that x increases from LEFT to RIGHT, y
 increases from BOTTOM to TOP, and z increases from BACK to FRONT.
 */
 /// method only checks whether the 4 corners are present and report missing ones, returning false.
bool ModelTopology::Infer_BOX_BOUNDARY_EdgeAndCornerFlagsFromSideFlags( const VSet<2U>& vset ) const
 {
     ErrorHandler& csmp_error ( ErrorHandler::Instance() );
     csmp_error.notice( WARNING, "ModelTopology<2>::Infer_BOX_BOUNDARY_EdgeAndCornerFlagsFromSideFlags:",
                       "Rectangular model has no edges.");
     // checking for the presence of corners
     std::set<BOX_BOUNDARY> corners;
     for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); ++bit )
       if ( (*bit) == CNR1 || (*bit) == CNR2 || (*bit) == CNR3 || (*bit) == CNR4 )
         corners.insert( static_cast<BOX_BOUNDARY>(*bit) );
      
     // reporting potential errors 
     if ( corners.size() != 4U ) {
          if ( !corners.empty() ) {
               std::cerr <<"\n\tdetected corners: ";
               for ( auto cit=corners.begin(); cit!=corners.end(); ++cit )
                 std::cerr << parseBoundary( (*cit) ) <<" "; 
               std::cerr << std::endl;
            }
          csmp_error.notice( ERROR, "ModelTopology<2>::Infer_BOX_BOUNDARY_EdgeAndCornerFlagsFromSideFlags:",
                            "Some of the model corners are not flagged.");
          return false;
       }
     return true;
 }
 
 
 
 
 
 
bool ModelTopology::FlagNodesUsingBoundaryRegions( VSet<3U>& vset ) const
  {
     ErrorHandler& csmp_error ( ErrorHandler::Instance() );

     set<string> boundary_regions;

     // 0. Extracting the names of valid BOUNDARY regions from file
     // -----------------------------------------------------------
     //  region name  etypes-of-region  elmt-ids for region
     //   map<string, pair<set<string>,vector<size_t> > >  model_regions;
     for ( auto& it : model_regions ) {
          const string region_name{ it.first };
          if ( region_name.find("BOUNDARY") != std::string::npos ||
               region_name.find("IRREGULAR") != std::string::npos ) {
               // if it is a lower dimensional region
               const auto etype = it.second.first.begin(); // .first.begin());
               if ( isLineElement( parseFiniteElementType(*etype) ) == true )
                 boundary_regions.insert( region_name );
            }
       }

     // 1. Making node ID sets for each of the standard boundaries
     // box boundaries
     // expects "BOTTOM","LEFT","RIGHT","TOP","FRONT","BACK"
     deque<size_t>  bottom, right, left, top, front, back,
                    front_left, front_right, front_bottom, front_top,
                    bottom_left, bottom_right,
                    back_bottom, back_left, back_top, back_right,
                    top_left, top_right, irregular;

     for ( auto it=ElementsOfRegionBegin("BOTTOM");
           it!=ElementsOfRegionEnd("BOTTOM"); it++ ) {
           // accessing contiguous ranges of element ID's with it->size_t
           for ( auto vit=vset.PlistBegin(*it); vit!=vset.PlistEnd(*it); vit++ )
             bottom.push_back( (*vit) );
       }
     for ( auto it=ElementsOfRegionBegin("LEFT");
           it!=ElementsOfRegionEnd("LEFT"); it++ ) {
           for ( auto vit=vset.PlistBegin(*it); vit!=vset.PlistEnd(*it); vit++ )
             left.push_back( (*vit) );
       }
     for ( auto it=ElementsOfRegionBegin("RIGHT");
           it!=ElementsOfRegionEnd("RIGHT"); it++ ) {
           for ( auto vit=vset.PlistBegin(*it); vit!=vset.PlistEnd(*it); vit++ )
             right.push_back( (*vit) );
       }
     for ( auto it=ElementsOfRegionBegin("TOP");
           it!=ElementsOfRegionEnd("TOP"); it++ ) {
           for ( auto vit=vset.PlistBegin(*it); vit!=vset.PlistEnd(*it); vit++ )
             top.push_back( (*vit) );
       }
     for ( auto it=ElementsOfRegionBegin("FRONT");
           it!=ElementsOfRegionEnd("FRONT"); it++ ) {
           for ( auto vit=vset.PlistBegin(*it); vit!=vset.PlistEnd(*it); vit++ )
             front.push_back( (*vit) );
       }
     for ( auto it=ElementsOfRegionBegin("BACK");
           it!=ElementsOfRegionEnd("BACK"); it++ ) {
           for ( auto vit=vset.PlistBegin(*it); vit!=vset.PlistEnd(*it); vit++ )
             back.push_back( (*vit) );
       }
     for ( auto& it : boundary_regions )
       for ( auto bit=ElementsOfRegionBegin(it.c_str());
             bit!=ElementsOfRegionEnd(it.c_str()); bit++ )
         for ( vector<int64_t>::iterator nit=vset.PlistBegin(*bit); nit!=vset.PlistEnd(*bit); ++nit )
           irregular.push_back( (*bit) );

     // making these containers unique
     sort( bottom.begin(), bottom.end() );
     bottom.erase( unique( bottom.begin(), bottom.end() ), bottom.end() );

     sort( right.begin(), right.end() );
     right.erase( unique( right.begin(), right.end() ), right.end() );

     sort( top.begin(), top.end() );
     top.erase( unique( top.begin(), top.end() ), top.end() );

     sort( left.begin(), left.end() );
     left.erase( unique( left.begin(), left.end() ), left.end() );

     sort( front.begin(), front.end() );
     front.erase( unique( front.begin(), front.end() ), front.end() );

     sort( back.begin(), back.end() );
     back.erase( unique( back.begin(), back.end() ), back.end() );

     sort( irregular.begin(), irregular.end() );
     irregular.erase( unique( irregular.begin(), irregular.end() ), irregular.end() );


     // -------------------------------------------------------------------------------
     // intersecting the sides to identify the edges
     // -------------------------------------------------------------------------------
     // FRONT_LEFT
     insert_iterator<deque<size_t> >  fl_it(front_left,front_left.begin());
     set_intersection( front.begin(), front.end(), left.begin(), left.end(), fl_it );

     // FRONT_RIGHT
     insert_iterator<deque<size_t> >  fr_it(front_right,front_right.begin());
     set_intersection( front.begin(), front.end(), right.begin(), right.end(), fr_it );

     // FRONT_BOTTOM
     insert_iterator<deque<size_t> >  fb_it(front_bottom,front_bottom.begin());
     set_intersection( front.begin(), front.end(), bottom.begin(), bottom.end(), fb_it );

     // FRONT_TOP
     insert_iterator<deque<size_t> >  ft_it(front_top,front_top.begin());
     set_intersection( front.begin(), front.end(), top.begin(), top.end(), ft_it );

     // BOTTOM_LEFT
     insert_iterator<deque<size_t> >  bl_it(bottom_left,bottom_left.begin());
     set_intersection( bottom.begin(), bottom.end(), left.begin(), left.end(), bl_it );

     // BOTTOM_RIGHT
     insert_iterator<deque<size_t> >  br_it(bottom_right,bottom_right.begin());
     set_intersection( bottom.begin(), bottom.end(), right.begin(), right.end(), br_it );

     // BACK_BOTTOM
     insert_iterator<deque<size_t> >  bb_it(back_bottom,back_bottom.begin());
     set_intersection( back.begin(), back.end(), bottom.begin(), bottom.end(), bb_it );

     // BACK_LEFT
     insert_iterator<deque<size_t> >  bal_it(back_left,back_left.begin());
     set_intersection( back.begin(), back.end(), left.begin(), left.end(), bal_it );

     // TOP_LEFT
     insert_iterator<deque<size_t> >  tl_it(top_left,top_left.begin());
     set_intersection( top.begin(), top.end(), left.begin(), left.end(), tl_it );

     // BACK_TOP
     insert_iterator<deque<size_t> >  bt_it(back_top,back_top.begin());
     set_intersection( back.begin(), back.end(), top.begin(), top.end(), bt_it );

     // TOP_RIGHT
     insert_iterator<deque<size_t> >  tr_it(top_right,top_right.begin());
     set_intersection( top.begin(), top.end(), right.begin(), right.end(), tr_it );

     // BACK_RIGHT
     insert_iterator<deque<size_t> >  bar_it(back_right,back_right.begin());
     set_intersection( back.begin(), back.end(), right.begin(), right.end(), bar_it );

     // 2. Flagging the nodes on the sides according to the boundaries
     // zapping all previous box boundary flags
     for ( auto bit=vset.BFlagsBegin(); bit!= vset.BFlagsEnd(); ++bit ) (*bit) = NOT;
     // assigning "IRREGULAR","BOTTOM","LEFT","RIGHT","TOP","FRONT","BACK"
     for ( auto& it : irregular ) vset.AddBFlag( it, IRREGULAR_OUTSIDE );
     for ( auto& it : bottom ) vset.AddBFlag( it, BOTTOM_OUTSIDE );
     for ( auto& it : left ) vset.AddBFlag( it, LEFT_OUTSIDE );
     for ( auto& it : right ) vset.AddBFlag( it, RIGHT_OUTSIDE );
     for ( auto& it : top ) vset.AddBFlag( it, TOP_OUTSIDE );
     for ( auto& it : front ) vset.AddBFlag( it, FRONT_OUTSIDE );
     for ( auto& it : back ) vset.AddBFlag( it, BACK_OUTSIDE );
     // edges
     for ( auto& it : back_bottom ) vset.AddBFlag( it, BACK_BOTTOM );
     for ( auto& it : back_right ) vset.AddBFlag( it, BACK_RIGHT );
     for ( auto& it : back_top ) vset.AddBFlag( it, BACK_TOP );
     for ( auto& it : back_left ) vset.AddBFlag( it, BACK_LEFT );
     for ( auto& it : bottom_left ) vset.AddBFlag( it, BOTTOM_LEFT );
     for ( auto& it : bottom_right ) vset.AddBFlag( it, BOTTOM_RIGHT );
     for ( auto& it : top_right ) vset.AddBFlag( it, TOP_RIGHT );
     for ( auto& it : top_left ) vset.AddBFlag( it, TOP_LEFT );
     for ( auto& it : front_bottom ) vset.AddBFlag( it, FRONT_BOTTOM );
     for ( auto& it : front_right ) vset.AddBFlag( it, FRONT_RIGHT );
     for ( auto& it : front_top ) vset.AddBFlag( it, FRONT_TOP );
     for ( auto& it : front_left ) vset.AddBFlag( it, FRONT_LEFT );

     // Flagging the corner nodes

     // The -Z axis (backward) facing plane of the model
     // CNR1
     {
       deque<size_t> corner;
       insert_iterator<deque<size_t> >  cit(corner,corner.begin());
       set_intersection( back_left.begin(), back_left.end(),
                         back_bottom.begin(), back_bottom.end(), cit );

       if ( corner.empty() )
         throw csmp::Exception( ERROR, "ModelTopology<3>::FlagNodesUsingBoundaryRegions",
                                                             "CNR1 could not be identified");
       else vset.AddBFlag( (*corner.begin()), CNR_MIN );
     }
     
     // CNR2
     {
       deque<size_t> corner;
       insert_iterator<deque<size_t> >  cit(corner,corner.begin());
       set_intersection( back_bottom.begin(), back_bottom.end(),
                         back_right.begin(), back_right.end(), cit );
       if ( corner.empty() )
         throw csmp::Exception( ERROR, "ModelTopology<3>::FlagNodesUsingBoundaryRegions",
                                                             "CNR2 could not be identified");
       else vset.AddBFlag( (*corner.begin()), CNR_MIN_MAXX );
     }
     // CNR3
     {
       deque<size_t> corner;
       insert_iterator<deque<size_t> >  cit(corner,corner.begin());
       set_intersection( back_right.begin(), back_right.end(),
                         back_top.begin(), back_top.end(), cit );
       if ( corner.empty() )
         throw csmp::Exception( ERROR, "ModelTopology<3>::FlagNodesUsingBoundaryRegions",
                                                             "CNR3 could not be identified");
       else vset.AddBFlag( (*corner.begin()), CNR_MAX_MAXX );
     }
     // CNR4
     {
       deque<size_t> corner;
       insert_iterator<deque<size_t> >  cit(corner,corner.begin());
       set_intersection( back_left.begin(), back_left.end(),
                         back_top.begin(), back_top.end(), cit );
       if ( corner.empty() )
         throw csmp::Exception( ERROR, "ModelTopology<3>::FlagNodesUsingBoundaryRegions",
                                                             "CNR4 could not be identified");
       else vset.AddBFlag( (*corner.begin()), CNR_MAX_MINXZ );
     }
     // The Z axis (forward) facing plane of the model
     // CNR5
     {
       deque<size_t> corner;
       insert_iterator<deque<size_t> >  cit(corner,corner.begin());
       set_intersection( front_left.begin(), front_left.end(),
                         front_bottom.begin(), front_bottom.end(), cit );
       if ( corner.empty() )
         throw csmp::Exception( ERROR, "ModelTopology<3>::FlagNodesUsingBoundaryRegions",
                                                             "CNR5 could not be identified");
       else vset.AddBFlag( (*corner.begin()), CNR_MIN_MAXZ );
     }
     // CNR6
     {
       deque<size_t> corner;
       insert_iterator<deque<size_t> >  cit(corner,corner.begin());
       set_intersection( front_right.begin(), front_right.end(),
                         front_bottom.begin(), front_bottom.end(), cit );
       if ( corner.empty() )
         throw csmp::Exception( ERROR, "ModelTopology<3>::FlagNodesUsingBoundaryRegions",
                                                             "CNR6 could not be identified");
       else vset.AddBFlag( (*corner.begin()), CNR_MIN_MAXXZ );
     }
     // CNR7
     {
       deque<size_t> corner;
       insert_iterator<deque<size_t> >  cit(corner,corner.begin());
       set_intersection( front_right.begin(), front_right.end(),
                         front_top.begin(), front_top.end(), cit );
       if ( corner.empty() )
         throw csmp::Exception( ERROR, "ModelTopology<3>::FlagNodesUsingBoundaryRegions",
                                                             "CNR7 could not be identified");
       else vset.AddBFlag( (*corner.begin()), CNR_MAX );
     }
     // CNR8
     {
       deque<size_t> corner;
       insert_iterator<deque<size_t> >  cit(corner,corner.begin());
       set_intersection( front_left.begin(), front_left.end(),
                         front_top.begin(), front_top.end(), cit );
       if ( corner.empty() )
         throw csmp::Exception( ERROR, "ModelTopology::FlagNodesUsingBoundaryRegions",
                                                             "CNR8 could not be identified");
       else vset.AddBFlag( (*corner.begin()), CNR_MAX_MAXZ );
     }

     if ( csmp_error.Verbose() ) {
          cout <<"\nModelTopology<3>::FlagNodesUsingBoundaryRegions: ";
          cout <<"Assigned CSMP associated boundary flags to the nodes."<< std::endl;
       }

     return true;

  } // end FlagNodesUsingBoundaryRegions

// DEBUGGING
//cerr <<"\nModelTopology::FlagNodesUsingBoundaryRegions: boundary flags:\n";
//size_t counter{0};
//for ( auto it=vset.BFlagsBegin(); it!=vset.BFlagsEnd(); ++it, ++counter )
//  if ( static_cast<BOX_BOUNDARY>(*it) != NOT ) {
//       cerr <<" "<< counter <<":"<< parseBoundary( static_cast<BOX_BOUNDARY>(*it) );
//    }
//cerr << endl << endl;

 

 } // end namespace csmp

