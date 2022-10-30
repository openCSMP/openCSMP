#include "EclipseMeshInterface_Example.h"
#include "Region.h"

#include "CSMP_definitions.h"
#include "VTU_Interface.h"
#include "InputDataManager.h"

#include "EclipseModel.h"
#include "VTK_Interface.h"


using namespace std;

namespace csmp {

	void EclipseMeshInterface_Example::Specifications()
	{
		SetTitle("EclipseMeshInterface_Example");
		SetDifficulty(1);
		SetCategory("Software Interfaces");
		AddAuthor("Roman Manasipov");
		AddDescription("source in: EclipseMeshInterface_Example.cpp");
		AddDescription("Example show's how to use EclispeModel constructor.");
		AddRequirement("One has to have corresponding mesh files to be present.");
		AddRequirement("by default: NPD5.grdecl, Johansen Data Set ( http://www.sintef.no/projectweb/matmora/downloads/johansen/ )");
	}

	void EclipseMeshInterface_Example::Run()
	{
		// Model name
    string model_name;
    cout << "\nPlease enter the name of input model, or press ENTER to use the default model 'NPD5':" << endl;
    cin.ignore();
    getline(cin, model_name);
    if (model_name.length() == 0) model_name = "NPD5";

		// Load Model
		std::string variables_file("CSMP_Eclipse_example-variables.txt");

    // create a working directory with current example name, go into this directory, and copy input files into it.
    CopyInputFiles(model_name, variables_file);

		// Setup mesh
		std::string regions_file(model_name);
		bool exclude_inactive_cells(true);
		bool tetra_mesh(false);
		const bool create_boundaries(false);
		EclipseModelSettings settings(model_name);
		settings.MeshSetup(regions_file,
			exclude_inactive_cells,
			tetra_mesh,
			create_boundaries);
		// Setup properties
		// porosity
		std::string   porosity_name("porosity");
		VARIABLE_TYPE porosity_type(SCALAR);
		PLACEMENT     porosity_place(ELEMENT);
		settings.PoroPropertySetup(porosity_name,
			porosity_type,
			porosity_place);
		// permeability
		std::string   permeability_name("tensor permeability");
		VARIABLE_TYPE permeability_type(TENSOR);
		PLACEMENT     permeability_place(ELEMENT);
		std::string   permeability_unit("mD"); // other options m2,D
		settings.PermPropertySetup(permeability_name,
			permeability_type,
			permeability_place,
			permeability_unit);
		// rock type
		std::string   rocktype_name("rock type");
		VARIABLE_TYPE rocktype_type(SCALAR);
		PLACEMENT     rocktype_place(ELEMENT);
		settings.RockNumPropertySetup(rocktype_name,
			rocktype_type,
			rocktype_place);

		EclipseModel modelOut(settings, model_name, variables_file);

		// Get Fault Regions
		std::vector<std::string> faults;
		modelOut.GetFaults(faults);

		// Get Well Regions
		std::vector<std::string> wells;
		modelOut.GetWells(wells);

		// We are in the process of rewriting the Eclipse interface, and the
		// following part is not yet fully ported.  - AJB

		Region<3U>&  model_domain(modelOut.Region("Model"));
		model_domain.UpdateMemberIndexes();

		// 3. eliminating any potentially disfunctional elements / cells from the model
		// ----------------------------------------------------------------------------
		// elements that have a negative Jacobian determinant are assumed to be degenerate and flagged for deletion
		vector<uint32_t> degenerate_elements;
		int volume_e_removed(0U), surface_e_removed(0U), line_e_removed(0U);
		for (auto it = model_domain.CellsBegin(); it != model_domain.CellsEnd(); ++it) {
			// find broken elements
			// (an element is regarded as broken if the determinant of its Jacobian inverse is negative at least
			//  at one of the integration points
			bool broken_elmt(false);
			for (size_t ipoint = 0U; ipoint<(*it)->IntegrationPoints(); ++ipoint)
				if ((*it)->det_JINV_AtIntegrationPoint(ipoint) <= 0.) {
					broken_elmt = true;
					break;
				}
			if (broken_elmt) {
				if ((*it)->IsVolume()) volume_e_removed++;
				else if ((*it)->IsSurface()) surface_e_removed++;
				else if ((*it)->IsLine()) line_e_removed++;
				degenerate_elements.push_back((*it)->Idx());
				auto eclipseCoord = modelOut.EclipseCoordinates(*it);
				std::cerr << "Broken element at " << eclipseCoord.i << ' ' << eclipseCoord.j << ' ' << eclipseCoord.k << ' '
					<< parseFiniteElementType((*it)->FE()->ElementType())
					<< '\n';
#if 0
				for (auto nit = (*it)->NodesBegin(); nit != (*it)->NodesEnd(); ++nit) {
					std::cerr << (*nit)->Coordinate() << '\n';
				}
				std::cerr << "Done\n";
#endif
			}
		}
		if (volume_e_removed > 0 || surface_e_removed > 0 || line_e_removed > 0) {
			cout << "\nread_and_configure_ECLIPSE_model: removing degenerate elements:\n";
			cout << "\n\tvolume elements removed:  " << volume_e_removed;
			cout << "\n\tsurface elements removed: " << surface_e_removed;
			cout << "\n\tline elements removed:    " << line_e_removed;
		}

		VTK_Interface<3U>  vtk_output;
		vtk_output.OutputNodeDataToVTK(modelOut, "EclipseInterfaceExample", 0);

		cout << "\nEclipseMeshInterface_Example: That's it!\n";

    fs::current_path("../../example_inputs/");
	}


  void EclipseMeshInterface_Example::CopyInputFiles(std::string& model_name, std::string& variable_file) {
    //find the name of current example source file
    string example_name = GetExampleFileName(__FILE__);
    //create a working directory with the name of this example and go into it
    fs::create_directory("../example_outputs");
    fs::current_path("../example_outputs");
    if (fs::is_directory(example_name)) fs::remove_all(example_name); //if directory already exists, delete it
    fs::create_directory(example_name);
    fs::current_path(example_name);

    string input_directory = fs::current_path().parent_path().parent_path();
    input_directory += "/example_inputs/input_meshes/";

    //copy Eclipse model files into working directory
    string path = "../../example_inputs/input_meshes/";
    string name = model_name + ".grdecl";
    string file_name = path + name;
    if (fs::exists(file_name)) fs::copy(file_name, "./");
    else {
      string error_message = "\n\nError: file '";
      error_message += (name + "' does not exist in directory " + input_directory);
      error_message += (", example cannot run, please copy this file into this directory\n");
      throw std::runtime_error(error_message);
    }

    //copy variable file into working directory
    input_directory = fs::current_path().parent_path().parent_path();
    input_directory += "/example_inputs/variables_and_configuration_files/";

    path = "../../example_inputs/variables_and_configuration_files/";
    file_name = path + variable_file;
    if (fs::exists(file_name)) fs::copy(file_name, "./");
    else {
      string error_message = "\n\nError: file '";
      error_message += (variable_file + "' does not exist in directory " + input_directory);
      error_message += (", example cannot run, please copy this file into this directory\n");
      throw std::runtime_error(error_message);
    }
  }

} // csmp
