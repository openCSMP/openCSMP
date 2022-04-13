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
		SetCategory("Software Functionality");
		AddAuthor("Roman Manasipov");
		AddDescription("source in: EclipseMeshInterface_Example.cpp");
		AddDescription("Example show's how to use EclispeModel constructor.");
		AddRequirement("One has to have corresponding mesh files to be present.");
		AddRequirement("by default: NPD5.grdecl, Johansen Data Set ( http://www.sintef.no/projectweb/matmora/downloads/johansen/ )");
	}

	void EclipseMeshInterface_Example::Run()
	{
		// Model name
		std::string model_name;
		cout << "\nPlease input model (mesh file prefix) name: ";
		cin >> model_name;

		// Load Model
		std::string variables_file("CSMP_Eclipse_example-variables.txt");

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
				if ((*it)->IsVolumeElement()) volume_e_removed++;
				else if ((*it)->IsSurfaceElement()) surface_e_removed++;
				else if ((*it)->IsLineElement()) line_e_removed++;
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
	}

} // csmp
