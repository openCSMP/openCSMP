#include "EclipseModel.h"
#include "ModelTopology.h"
#include "Region.h"
#include "ModelTime.h"
#include "variableOperations.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

/// Model constructor with provided "variables_file.txt" file is used
EclipseModel::EclipseModel(EclipseModelSettings& settings,
	const string& model_name,
	const string& variables_file)
	: csmp::Model<3U>(variables_file.c_str()),
	eclipse_model_settings_(settings)
{
	this->Name(model_name.c_str());
	Initialize();
}


/// Default Model constructor with empty property data base is called
EclipseModel::EclipseModel(EclipseModelSettings& settings,
	const string& model_name)
	: csmp::Model<3U>(),
	eclipse_model_settings_(settings)
{
	this->Name(model_name.c_str());
	Initialize();
}


EclipseModel::~EclipseModel()
{
}


/**
Master function that BUILDS CSMP MODEL FROM ECLIPSE DATA
*/
void EclipseModel::Initialize()
{
	csmp::ErrorHandler& error_handler(csmp::ErrorHandler::Instance());

	double& model_time(csmp::ModelTime::Instance().modelTime);
	model_time = 0.;

	try {
		csmp::VSet<3U>  vset;
		bool isoparametric_elements(true);

		csmp::ModelTopology  mesh_topology(isoparametric_elements);
		
		// =====================================================================
		// 0. reads grid from ECLIPSE input files and converts into CSMP mesh
		// =====================================================================
		mesh_interface.SetProperties(eclipse_model_settings_.properties_);

		// KEY METHOD here
		mesh_interface.ReadFile( vset, mesh_topology,
                             eclipse_model_settings_.mesh_file_prefix_,
                             eclipse_model_settings_.exclude_inactive_cells_,
                             eclipse_model_settings_.tetra_mesh_);

		// =====================================================================
		// 1. selectively read properties of interest, adding them to VSET
		// =====================================================================
		// porosity, permeability, saturations
		set<string> vset_props;
		for ( auto epit = eclipse_model_settings_.properties_.begin(); epit != eclipse_model_settings_.properties_.end(); ++epit)
		{
			const csmp::Parameter& prop = (*epit).second;
			auto prop_it = vset_props.find(prop.name);
			if ( prop_it != vset_props.end())
			{
				if (!this->Database().IsDefined(prop.name.c_str()))
					// SKM FIX 
					this->Database().AddProperty( prop.name.c_str(), prop.notation.c_str(), prop.unit.c_str(),
						                            prop.key.type, prop.key.place,
                                        this->Database().VariableCount(prop.key.place, prop.key.type), prop.min, prop.max,
						                            prop.usage.c_str());
				vset_props.erase(prop_it);
			}
		}

		/// checking whether any undefined properties are left
		set<string> undefined_props;
		for ( auto prop_it = vset_props.begin(); prop_it != vset_props.end(); ++prop_it)
			if (!this->Database().IsDefined((*prop_it).c_str()))
				undefined_props.insert(*prop_it);
		if (!undefined_props.empty())
		{
			if (error_handler.Verbose())
			{
				string message = "VSet contains data for undefined properties: ";
				for ( auto prop_it = undefined_props.begin(); prop_it != undefined_props.end(); ++prop_it)
				{
					if (prop_it != undefined_props.begin())
						message += ", ";
					message += *prop_it;
				}
				message += " !!!";
				error_handler.Note(csmp::INFO, "EclipseModel<3U>::BuildModel", message.c_str());
			}
		}
    
		// =====================================================================
		// 2. construct CSMP model from obtained topology and mesh in vset
		// =====================================================================
		//    we won't use eclipse neighbor info since it includes neighbor information
		//    of elements of different dimensionality (i.e. e volumetric element has a surface element neighbors )
		//    all cells are lumped into the region "Eclipse Model" that is stored in the model topology
		csmp::Model<3U>::Initialize( mesh_topology, vset );

   }
	// ---------------------------------------------------
	// catching all possible standard and csmp::Exceptions
	// ---------------------------------------------------
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
	catch (csmp::Exception& ba) {
		cout << "\nException: Exception raised by: " << ba.What() << endl;
		cout << "\nDiagnostics:" << endl;
		ba.Out();
	}
} // end


/// existing special regions
template<class Container>
void EclipseModel::GetRegions(Container& data)
{
	mesh_interface.GetRegions(data);
}

template void EclipseModel::GetRegions(vector<string>&);
template void EclipseModel::GetRegions(list<string>&);
template void EclipseModel::GetRegions(set<string>&);


template<class Container>
void EclipseModel::GetFaults(Container& data)
{
	mesh_interface.GetFaults(data);
}

template void EclipseModel::GetFaults(vector<string>&);
template void EclipseModel::GetFaults(list<string>&);
template void EclipseModel::GetFaults(set<string>&);


template<class Container>
void EclipseModel::GetWells(Container& data)
{
	mesh_interface.GetWells(data);
}

template void EclipseModel::GetWells(vector<string>&);
template void EclipseModel::GetWells(list<string>&);
template void EclipseModel::GetWells(set<string>&);



void EclipseModel::AddWell(const string& well_name, const Point<3U>& well_start_point, const Point<3U>& well_end_point)
{
	mesh_interface.AddWell(well_name, well_start_point, well_end_point);
}



/// processing special regions
void EclipseModel::CreateBoundariesAroundFaults( bool keep_fault_regions )
{
	// create external boundaries from regions whose name contains 'BOUNDARY' or a box boundary identifier
  // ( the input regions are deleted)
  EstablishBoundariesFromRegions();

   // combine all regions imported under the category of faults into a single one call faults
	this->MergeRegions( faults_, "FAULTS");
	faults_.insert("FAULTS");
// 	this->CreateInternalBoundaryFrom( "FAULTS" );
}


void EclipseModel::CreateSplitBoundariesAroundFaults( bool delete_fault_regions )
{
	this->MergeRegions(faults_, "FAULTS");
	faults_.insert("FAULTS");
	this->CreateSplitBoundaryFrom("FAULTS");

	// create splitboundaries
	//for( set<string>::const_iterator
	//     rit = faults_.begin(); rit != faults_.end(); ++rit )
	//    this->InsertSplitBoundary( (*rit), delete_fault_regions );
	//faults_.insert("FAULTS");

}




/**
BOX flag nodes and elements of volumetric target region
This method assumes that the nodes of the mesh were previously flagged correctly.
*/
void EclipseModel::AssignBoxBoundaryFlagsWherePossible(const char* target_region)
{
	csmp::ErrorHandler& error_handler(csmp::ErrorHandler::Instance());

	auto& domain = this->Region(target_region);
	vector<uint32_t>        fnids;
	multimap<size_t, pair<BOX_BOUNDARY, Node<3U>*> >  boundary_nodes;
	vector<double>      nrml, nrml_right, nrml_left, nrml_top, nrml_bottom, nrml_front, nrml_back;
	Box                   box;
	double              minLength(0.71); // dot-product of 2 unit vectors at an angle >=45 degrees
	BOX_BOUNDARY          bflag(NOT);

	const uint32_t dim(3U);
	box.UnitNormalTo(BOTTOM, dim, nrml_bottom);
	box.UnitNormalTo(TOP, dim, nrml_top);
	box.UnitNormalTo(LEFT, dim, nrml_left);
	box.UnitNormalTo(RIGHT, dim, nrml_right);
	box.UnitNormalTo(FRONT, dim, nrml_front);
	box.UnitNormalTo(BACK, dim, nrml_back);

	cout << "\nAssignBoxBoundaryFlagsWherePossible: scanning hexahedral elements for boundary adffiliation...";
	for (auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it)
	{
		// ignore elements that are not hexahedra
		if ((*it)->FE_Type() != ISOPARAMETRIC_LINEAR_HEXAHEDRON) {
			cout << "\n\tignored: " << parseFiniteElementType((*it)->FE_Type());
			continue;
		}
		// idea: loop over the faces of the cell and where there is no neighbor
		// check in which direction the face normal is pointing, assign boundary flags accordingly
		// if the element has more than one face at the boundary, idenfify it as an edge or a corner
		for (auto i{0U}; i<(*it)->Faces(); ++i)
			if ((*it)->Neighbor(i) == nullptr) {
				// determining in which direction the face normal points
				(*it)->UnitNormalToFace(i, nrml);
				// projecting: perfect alignment would give dot-product equal 1, inclinations up to 37 degrees cos(37)~0.8 are tolerated
				if (dotProduct<3U>(nrml, nrml_bottom) >= minLength) bflag = BOTTOM;
				else if (dotProduct<3U>(nrml, nrml_top) >= minLength) bflag = TOP;
				else if (dotProduct<3U>(nrml, nrml_left) >= minLength) bflag = LEFT;
				else if (dotProduct<3U>(nrml, nrml_right) >= minLength) bflag = RIGHT;
				else if (dotProduct<3U>(nrml, nrml_front) >= minLength) bflag = FRONT;
				else if (dotProduct<3U>(nrml, nrml_back) >= minLength) bflag = BACK;
				// getting the nodes for flagging the faces
				for ( const auto& j : (*it)->FE()->NodesOfFace(i) ) {
					(*it)->N(j)->AtBoundary(bflag);
					// storing the nodes to determine which ones lie on EDGES (duplicates) or even corners (triplicates)
					//                   local node #            flag   local node #
					boundary_nodes.insert(make_pair(j, make_pair(bflag, (*it)->N(j))));
				}
			}
		// flagging elements with duplicate and triplicate boundary nodes accordingly
		for ( auto n{0U}; n<(*it)->Nodes(); ++n) {
			// dealing with any cases where there are multiple boundary flags
			// duplicates = edges
			if (boundary_nodes.count(n) == 2U) {
				// figuring out which edge we are on
				auto range = boundary_nodes.equal_range(n);
				bflag = whichEdge((*range.first).second.first, (*range.second).second.first);
				(*it)->N(n)->AtBoundary(bflag);
			}
			// triplicates = corners
			else if (boundary_nodes.count(n) == 3U) {
				// figuring out which corner we have found
				auto range = boundary_nodes.equal_range(n);
				auto it_2nd(range.first); it_2nd++;
				bflag = whichCorner((*range.first).second.first, (*it_2nd).second.first, (*range.second).second.first);
				(*it)->N(n)->AtBoundary(bflag);
			}
			else if (boundary_nodes.count(n) > 3U) { // potentially a hexahedron which sits at a model edge (7-boundary nodes)
				(*it)->Out();
				cerr << "\n\tdetected " << boundary_nodes.count(n) << " boundary flags for element " << (*it)->Idx();
				error_handler.Note(WARNING, "EclipseModel<3U>::AssignBoxBoundaryFlagsWherePossible:",
					"this may be a completely disconnected element.");
			}
		}

		// testing
		//cerr <<"\n("<< (*it)->Idx() <<"): ";
		//for ( auto n=0U; n<(*it)->Nodes(); ++n )
		//  cerr << parseBoundary( (*it)->N(n)->AtBoundary() ) <<" ";

		// resetting
		boundary_nodes.clear();
		bflag = NOT;
	}

} // end AssignBoxBoundaryFlagsWherePossible

} // end csmp
