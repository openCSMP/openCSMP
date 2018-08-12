#include "ANSYS_Model3D.h"
#include "ANSYS_Interface.h"
#include "Element.h"
#include "Region.h"
#include "Box.h"
#include "Exception.h"
#include "ModelTopology.h"
#include "ModelTime.h"

using namespace std;

namespace csmp {

/** Builds Model after it was constructed with the default constructor.
    
    The following steps are performed:
    
    1. initializes model topology and ansys interface classes from data contained in *.asc and *.dat files
    
    2. Construct the model from topology and vset containers initialised by ANSYS_Interface object

    The initialisation process involves the following steps:
    
    2.1 create Boundary objects from lower-dimensional regions whose name begins with the string BOUNDARY
    
    - if the model is box-shaped:
      - establish the BOX_BOUNDARY flagging
      - create corresponding Boundary objects with the box boundary names
      - create Boundary edge objects (in 3D) that represent the model edges
 
    @attention model time is initialised to zero

    @attention per default isoparametric is true and the connectivity information
    between the elements from ANSYS is not used, but this data is recreated
*/

	void ANSYS_Model3D::Initialize(bool isoparametric,
		const char* mesh_file_set,
		const char* regions_file_prefix,
		bool irregular_mesh,
		bool binary_input_file,
		bool use_regions_file,
		bool create_boundaries)
	{
		double64& model_time(ModelTime::Instance().modelTime);
		model_time = 0.;

		// -------------------------------------------------
		// initializing the empty Model from the ANSYS
		// data imported into a vset.
		// -------------------------------------------------
		try {
			VSet<3U>  vset;
			bool isoparametric_elements(isoparametric);

			ModelTopology   mesh_topology(isoparametric_elements);
			ANSYS_Interface mesh_interface(isoparametric_elements);

			// 0. reading the mesh from ANSYS-CSMP-input files
			mesh_interface.Read_ANSYS_Mesh(std::string(mesh_file_set), vset, mesh_topology, binary_input_file, irregular_mesh);

			// 1. preserving numbered node coordinates in a vector
			const size_t vertices(vset.Vertices());
			node_coords_.reserve(vertices);
			for (size_t i = 0U; i<vertices; ++i)
				node_coords_.emplace_back(Point<3U>(vset.Px(i), vset.Py(i), vset.Pz(i)));

			// 2. construct model based on obtained model topology and vset
			if (use_regions_file)
				Model<3U>::Initialize(regions_file_prefix,
					mesh_topology,
					vset,
					create_boundaries,
					irregular_mesh);
			else
				Model<3U>::Initialize(mesh_topology,
					vset,
					create_boundaries,
					irregular_mesh);

		}

		// -------------------------------------------------
		// catching all possible standard and csmp::Exceptions
		// -------------------------------------------------
		catch (bad_alloc& ba) {
			cout << "\nbad_alloc: Memory allocation error caused by: " << ba.what() << endl;
		}
		catch (bad_cast& ba) {
			cout << "\nbad_cast: Type casting error caused by: " << ba.what() << endl;
		}
		catch (bad_exception& ba) {
			cout << "\nbad_exception: Exception error caused by: " << ba.what() << endl;
		}
		catch (bad_typeid& ba) {
			cout << "\nbad_typeid: Type ID error caused by: " << ba.what() << endl;
		}
		catch (ios_base::failure& ba) {
			cout << "\nios_base::failure: Probable I/O error caused by: " << ba.what() << endl;
		}
		// standard logic errors
		catch (domain_error& ba) {
			cout << "\ndomain_error: Logic error caused by: " << ba.what() << endl;
		}
		catch (invalid_argument& ba) {
			cout << "\ninvalid_argument: Logic error caused by: " << ba.what() << endl;
		}
		catch (length_error& ba) {
			cout << "\nlength_error: Logic error caused by: " << ba.what() << endl;
		}
		catch (out_of_range& ba) {
			cout << "\nout_of_range: Logic error caused by: " << ba.what() << endl;
		}
		// runtime errors
		catch (overflow_error& ba) {
			cout << "\noverflow_error: Runtime error caused by: " << ba.what() << endl;
		}
		catch (range_error& ba) {
			cout << "\nrange_error: Runtime error caused by: " << ba.what() << endl;
		}
		catch (underflow_error& ba) {
			cout << "\nunderflow_error: Runtime error caused by: " << ba.what() << endl;
		}
		catch (Exception& ba) {
			cout << "\nException: Exception raised by: " << ba.What() << endl;
			cout << "\nDiagnostics:" << endl;
			ba.Out();
		}

	} // end Initialize

void ANSYS_Model3D::Initialize( const char* mesh_file_set,
                                const char* regions_file_prefix,
                                bool irregular_mesh,
                                bool binary_input_file,
                                bool use_regions_file,
                                bool create_boundaries )
{
  double64& model_time( ModelTime::Instance().modelTime );
  model_time = 0.;
  
  // -------------------------------------------------
  // initializing the empty Model from the ANSYS
  // data imported into a vset.
  // -------------------------------------------------
  try {
       VSet<3U>  vset;
       bool isoparametric_elements( true );

       ModelTopology   mesh_topology ( isoparametric_elements );
       ANSYS_Interface mesh_interface( isoparametric_elements );

       // 0. reading the mesh from ANSYS-CSMP-input files
       mesh_interface.Read_ANSYS_Mesh( std::string(mesh_file_set), vset, mesh_topology, binary_input_file, irregular_mesh );
       // ATTENTION (comment from SKM): Since ANSYS does not output the neighbour connectivity correctly,
       // the 'pfverts' neighbor container is zapped here so that VData does not think anymore that it has neighbor connectivity
       // later on this connectivity will be recreated inside of the Model where suitable machinery exists.
       vset.RemovePfverts();

       // 1. preserving numbered node coordinates in a vector
       const size_t vertices(vset.Vertices());
       node_coords_.reserve(vertices);
       for ( size_t i=0U; i<vertices; ++i )
         node_coords_.emplace_back( Point<3U>(vset.Px(i),vset.Py(i),vset.Pz(i)) );

       // 2. construct model based on obtained model topology and vset
       if ( use_regions_file )
         Model<3U>::Initialize( regions_file_prefix,
                                mesh_topology,
                                vset,
                                create_boundaries,
                                irregular_mesh );
       else
         Model<3U>::Initialize( mesh_topology,
                                vset,
                                create_boundaries,
                                irregular_mesh );

    }

  // -------------------------------------------------
  // catching all possible standard and csmp::Exceptions
  // -------------------------------------------------
  catch( bad_alloc& ba ) {
       cout <<"\nbad_alloc: Memory allocation error caused by: "<< ba.what() << endl;
    }
  catch( bad_cast& ba ) {
       cout <<"\nbad_cast: Type casting error caused by: "<< ba.what() << endl;
    }
  catch( bad_exception& ba ) {
       cout <<"\nbad_exception: Exception error caused by: "<< ba.what() << endl;
    }
  catch( bad_typeid& ba ) {
       cout <<"\nbad_typeid: Type ID error caused by: "<< ba.what() << endl;
    }
  catch( ios_base::failure& ba ) {
       cout <<"\nios_base::failure: Probable I/O error caused by: "<< ba.what() << endl;
    }
  // standard logic errors
  catch( domain_error& ba ) {
       cout <<"\ndomain_error: Logic error caused by: "<< ba.what() << endl;
    }
  catch( invalid_argument& ba ) {
       cout <<"\ninvalid_argument: Logic error caused by: "<< ba.what() << endl;
    }
  catch( length_error& ba ) {
       cout <<"\nlength_error: Logic error caused by: "<< ba.what() << endl;
     }
  catch( out_of_range& ba ) {
       cout <<"\nout_of_range: Logic error caused by: "<< ba.what() << endl;
    }
  // runtime errors
  catch( overflow_error& ba ) {
       cout <<"\noverflow_error: Runtime error caused by: "<< ba.what() << endl;
    }
  catch( range_error& ba ) {
       cout <<"\nrange_error: Runtime error caused by: "<< ba.what() << endl;
    }
  catch( underflow_error& ba ) {
       cout <<"\nunderflow_error: Runtime error caused by: "<< ba.what() << endl;
    }
  catch( Exception& ba ) {
       cout <<"\nException: Exception raised by: "<< ba.What() << endl;
       cout <<"\nDiagnostics:"<< endl;
       ba.Out();
    }

 } // end Initialize





// ==============================================================================
// TESTING of Initialise() method
// ==============================================================================
#ifdef DEBUG_ANSYS_MODEL3D_H
assert( old_and_new_elmtids.size() == vset.Elements() );
cout <<"\nprinting element mapping: old,new: ";
size_t  counter(0U);
for ( map<size_t,size_t>::const_iterator
      it=old_and_new_elmtids.begin(); it!=old_and_new_elmtids.end(); it++ )
  cout << (*it).first <<","<< (*it).second <<"  ";
cout << endl;

cout <<"\n\n\nPRINTING REMAINING REGIONS AND THEIR ELEMENT TYPES: "<< endl;
list<string>  regions;
mesh_topology.Out( regions );
for ( list<string>::const_iterator
      git=regions.begin(); git!=regions.end(); git++ ) {
     cout <<"\nregion: "<< *git <<" new e-number and (type): "<< endl;
      for ( vector<size_t>::const_iterator
            rit=mesh_topology.ElementsOfRegionBegin( (*git).c_str() );
            rit!=mesh_topology.ElementsOfRegionEnd( (*git).c_str() ); rit++ )
        {
//            map<size_t,size_t>::const_iterator eidx=old_and_new_elmtids.find()
            cout << *rit <<"("<< vset.ElementType( counter++ ) <<")  ";
        }
  }
cout << endl << endl;
cout <<"\ntesting node associations: "<< endl;
set<size_t>  nodes;
for ( size_t i=0U; i<vset.Elements(); i++ ) {
     cout << i <<"("<< vset.ElementType(i) <<"): nodes: ";
     for ( vector<size_t>::const_iterator
           pit=vset.PlistBegin(i); pit!=vset.PlistEnd(i); pit++ ) {
          cout << (*pit) <<" ";
          nodes.insert( (*pit) );
       }
     cout << endl;
  }
assert( nodes.size() == vset.Vertices() );
cout <<"\nnew node numbers: "<< endl;
for ( set<size_t>::const_iterator
      sit=nodes.begin(); sit!=nodes.end(); sit++ )
  cout << *sit <<" ";
cout << endl << endl;
#endif
// ==============================================================================
// END OF TESTING
// ==============================================================================


// inlines

std::vector<Point<3U> >::const_iterator ANSYS_Model3D::VerticesBegin() const { return node_coords_.begin(); }

std::vector<Point<3U> >::const_iterator ANSYS_Model3D::VerticesEnd() const { return node_coords_.end(); }







/**

  ANSYS_Model3D::ANSYS_Model3D(const char* icem_file_set,
                               const char* variable_file,
                               const char* regions_file_prefix,
                               bool irregular_mesh,
                               bool binary_file )
\par Description:
Default constructor of Model is called.
The model is built from the 'icem_file_set' '*.asc' and '*.dat',
variables file, and the configuration file prefix is used to read the regions file.
Uses the method Initialize.

Added by: Julian E. Mindel 16-03-2012

*/

ANSYS_Model3D::ANSYS_Model3D( const char* icem_file_set,
                              const char* regions_file_prefix,
                              const char* variable_file,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file,
                              bool create_boundaries)
 : Model<3U>( variable_file, false )
 {
    this->Name( icem_file_set );
    Initialize( icem_file_set,
                regions_file_prefix,
                irregular_mesh,
                binary_file,
                use_regions_file,
                create_boundaries);
 }

/**

ANSYS_Model3D::ANSYS_Model3D( const char* icem_file_set,
const char* variable_file )
\par Description:
Choose isoparametric or non-isoparametric!
Default constructor of Model is called. Then the model is build
from the 'icem_file_set' '*.asc' and '*.dat' files using the method
Initialize.
*/

ANSYS_Model3D::ANSYS_Model3D(bool isoparametric,
	const char* icem_file_set,
	const char* variable_file,
	bool irregular_mesh,
	bool binary_file,
	bool use_regions_file,
	bool create_boundaries)
	: Model<3U>(variable_file, false)
{
	this->Name(icem_file_set);
	Initialize(isoparametric,
		icem_file_set,
		icem_file_set,
		irregular_mesh,
		binary_file,
		use_regions_file,
		create_boundaries);
}
/**

  ANSYS_Model3D::ANSYS_Model3D( const char* icem_file_set,
                                const char* variable_file )
\par Description:

Default constructor of Model is called. Then the model is build
from the 'icem_file_set' '*.asc' and '*.dat' files using the method
Initialize.
*/

ANSYS_Model3D::ANSYS_Model3D( const char* icem_file_set,
                              const char* variable_file,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file,
                              bool create_boundaries)
 : Model<3U>( variable_file, false )
 {
    this->Name( icem_file_set );
    Initialize( icem_file_set,
                icem_file_set,
                irregular_mesh,
                binary_file,
                use_regions_file,
                create_boundaries);
 }




/**
@ Description:
Default constructor of Model is called.
The model is built from the 'icem_file_set' '*.asc' and '*.dat',
variables file, and the configuration file prefix is used to read the regions file.
Uses the method Initialize.

@author Julian E. Mindel 16-03-2012
*/
ANSYS_Model3D::ANSYS_Model3D( const char* icem_file_set,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file,
                              bool create_boundaries )
{
    this->Name(icem_file_set);
    Initialize( icem_file_set,
                icem_file_set,
                irregular_mesh,
                binary_file,
                use_regions_file,
                create_boundaries);
}



/// nothing needs to be done here
ANSYS_Model3D::~ANSYS_Model3D()
 {
 }
 
 
 
/** 
     Renumbers the nodes (0..n) as in the original ANSYS model.
     
     @return returns whether any changes in the numbering were made.
     
     @author SKM
     @date 6/12/2016
*/
bool ANSYS_Model3D::RestoreOriginalNodeNumbering( bool verbose )
 {
    // making a binary tree of the original node numbers, searchable for point coordinates
    map<Point<3U>,size_t>  original_node_numbers;
    for ( size_t i=0U; i<node_coords_.size(); ++i )
      original_node_numbers.insert( make_pair( node_coords_[i], i ) );
   
    // renumbering the nodes of the model consecutively
    bool first_call(true), made_changes(false);
    const auto nodesEnd(Mesh().NodesEnd());
    const auto onodesEnd(original_node_numbers.end());
    for ( auto nit=Mesh().NodesBegin(); nit!=nodesEnd; ++nit ) {
         auto onit( original_node_numbers.find( (*nit).Coordinate() ) );
         if ( onit != onodesEnd ) {
              if ( (*nit).Idx() != (*onit).second ) {
                   if ( first_call ) {
                        if ( verbose ) cout <<"\nANSYS_Model3D::RestoreOriginalNodeNumbering: changed indices of following nodes:";
                        first_call   = false;
                        made_changes = true;
                     }
                   if ( verbose ) cout <<"\n\t"<< (*nit).Idx() <<" -> "<< (*onit).second;
                }
              (*nit).Idx( (*onit).second );
           }
         else
           throw csmp::Exception( ERROR, "ANSYS_Model3D::RestoreOriginalNodeNumbering:",
                                 "node could not be identified; has it been newly created?" );
      }
   
   return made_changes;
   
 } // end RestoreOriginalNodeNumbering

 
 
 
 
 
 

} // end csmp
