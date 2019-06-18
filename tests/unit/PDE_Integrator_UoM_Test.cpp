#include "PDE_Integrator_UoM_Test.h"
#include "GaussJordan_Solver.h"



using namespace std;

namespace csmp {

  PDE_Integrator_UoM_Test::PDE_Integrator_UoM_Test(std::string modelName) {
    const bool       isoparametric(true);
    ANSYS_Interface  mesh_interface(isoparametric);  // true = isoparametric elements
    VSet<2U>         mesh_container;
    ModelTopology    mesh_topology(isoparametric);   // true = isoparametric elements

    string mesh_name(modelName);
    const bool binary_file(true);
    const bool irregular_mesh(false); // irregular_boundary! BOX SHAPE or IRREGULAR SHAPE
    mesh_interface.Read_ANSYS_Mesh(mesh_name.c_str(), mesh_container, mesh_topology, binary_file, irregular_mesh);
    cout << "Finished reading mesh..." << endl;
    cout << "Building Model..." << endl;
    model = new Model<2U>(mesh_topology, mesh_container, "PDE_Integrate_UoM_2phase-variables.txt");
   
    Reset();

    // Create pde-operators
    pressureLHS = new NumIntegral_dNT_op_dN_dV<2U, Element<2U> >(model->Database(),
      "permeability",
      "fluid pressure",
      "fluid pressure");

    sourceVolume = new NumIntegral_NT_op_N_dV<2U, Element<2U> >(model->Database(),
      "fluid volume source", "fluid pressure");

    sourcePoint = new PointSource_rhsop<2U, Element<2U>>(model->Database(),
      "nodal fluid point source", "fluid pressure");
    
    gravityTerm = new NumIntegral_dNT_op_dV<2U, Element<2U>>(model->Database(),
      "gravity vector",
      "fluid pressure");

    //fluid_velocity = new VelocityAndVolumeFlux<2U, Element<2U>>(*model, "total mobility", "porosity", "fluid pressure", true, "velocity");

    temperatureLHS = new NumIntegral_dNT_op_dN_dV<2U, Element<2U> >(model->Database(),
      "permeability",
      "temperature",
      "temperature");

    sourceHeat = new NumIntegral_NT_op_N_dV<2U, Element<2U>>(model->Database(),
      "thermal volume source", "temperature");

  }

  PDE_Integrator_UoM_Test::~PDE_Integrator_UoM_Test() {

    if (model != nullptr) {
      delete model;
    }
    if (sourceVolume != nullptr) {
      delete sourceVolume;
    }
    if (sourcePoint != nullptr) {
      delete sourcePoint;
    }
    if (gravityTerm != nullptr) {
      delete gravityTerm;
    }
    if (pressureLHS != nullptr) {
      delete pressureLHS;
    }
    if (temperatureLHS != nullptr) {
      delete temperatureLHS;
    }
    if (sourceHeat != nullptr) {
      delete sourceHeat;
    }
    if (fluid_velocity != nullptr) {
      delete fluid_velocity;
    }
    /*
    if (pde_validate != nullptr) {
      delete pde_validate;
    }
    if (pde_test != nullptr) {
      delete pde_test;
    }
    */
  }

  void PDE_Integrator_UoM_Test::Reset() {
    Region<2U> region = model->Region("Model");
    size_t elementNum = region.Elements();

    // Set values on nodes and elements
    Index pressureKey = model->Database().StorageKey("fluid pressure");
    Index pressureValKey = model->Database().StorageKey("fluid pressure previous");
    Index temperatureKey = model->Database().StorageKey("temperature");
    Index fluidVolumeKey = model->Database().StorageKey("fluid volume source");
    Index fluidNodeSourceKey = model->Database().StorageKey("nodal fluid point source");
    Index heatSourceKey = model->Database().StorageKey("thermal volume source");
    Index gravityVectorKey = model->Database().StorageKey("gravity vector");

    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      (*nIter)->Store(pressureKey, makeScalar(PLAIN, (*nIter)->Idx()));
      (*nIter)->Store(pressureValKey, makeScalar(PLAIN, (*nIter)->Idx()));
      (*nIter)->Store(temperatureKey, makeScalar(PLAIN, (*nIter)->Idx()));
      (*nIter)->Store(fluidNodeSourceKey, makeScalar(PLAIN, (*nIter)->Idx()));
    }

    for (auto eIter = region.ElementsBegin(); eIter != region.ElementsEnd(); ++eIter) {
      (*eIter)->Store(fluidVolumeKey, makeScalar(PLAIN, (*eIter)->Idx()));
      (*eIter)->Store(heatSourceKey, makeScalar(PLAIN, (*eIter)->Idx()) / elementNum);
      (*eIter)->Store(gravityVectorKey, makeVector(PLAIN, PLAIN, (*eIter)->BaryCenter()[0], (*eIter)->BaryCenter()[1]));
    }


    model->InputPropertyValue("permeability", makeScalar(PLAIN, 1.));
    model->InputPropertyValue("total mobility", makeScalar(PLAIN, 1.)); // conductivity
                                                                        // Set boundary conditions
    model->InputBoundaryValue(BOTTOM, "fluid pressure", makeScalar(DIRICH, 1.0));
    model->InputBoundaryValue(BOTTOM, "fluid pressure previous", makeScalar(DIRICH, 1.0));
    model->InputBoundaryValue(BOTTOM, "temperature", makeScalar(DIRICH, 1.0));

  
  }

  void PDE_Integrator_UoM_Test::run()
  {
    //=======================================
    // test single variable
    //=======================================

    // test scalar variable
    TestSingleVariable();
    TestOutputSingleVariable();
    // test vector variable

    // test array variable

    // test flagged array variable

    // test tensor variable




    //=======================================
    // test multiple single variable
    //=======================================

    // test 2 scalar variables
    TestTwoScalarVariables();

    // test scalar variable and vector variable

    // test vector variable and vector variable

  }

  void PDE_Integrator_UoM_Test::TestSingleVariable() {
    // 0. prepare
    Region<2U>& region = model->Region("Model");
    region.RenumberNodes();

    GaussJordan_Solver solver;

    pde_validate = new PDE_Integrator_UoM_Mock<2U, Region>(solver);
    pde_test = new PDE_Integrator_UoM_Mock<2U, Region>(solver);

    pde_validate->Add(pressureLHS);
    pde_validate->Add(sourceVolume);
    pde_validate->Add(sourcePoint);
    pde_validate->Add(gravityTerm);
    pde_test->Add(pressureLHS);
    pde_test->Add(sourceVolume);
    pde_test->Add(sourcePoint);
    pde_test->Add(gravityTerm);

    const std::vector<size_t>& DOF_indexes = pde_test->GetDOFIndex();

    // 1. test matrix establish and enumerate DOFs
    pde_validate->EstablishMatrixSetup(region);
    pde_test->EstablishMatrixSetupTest(region);

    size_t pressureDirchletDOFs(0);
    Index pressureKey = model->Database().StorageKey("fluid pressure");
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) == DIRICH) {
        pressureDirchletDOFs += 1;
      }
    }

    _test((pde_validate->GetRH()->size()) == (pde_test->GetRH()->size()));

    pde_test->EnumerateAndFixMatrixSize(region);
    _test(pde_validate->GetRH()->size() == DOF_indexes.size());
    _test(pde_validate->GetRH()->size() == pde_test->GetRH()->size() + pressureDirchletDOFs);
    _test(pde_validate->GetG()->Rows() == pde_test->GetG()->Rows() + pressureDirchletDOFs);
    _test(pde_validate->GetG()->Cols() == pde_test->GetG()->Cols() + pressureDirchletDOFs);

    // 2. test accumulate
    pde_validate->Accumulate(region);
    pde_test->AccumulateTest(region);

    // 3. test AssignInitialConditions & LateAccumulate
    pde_validate->TimeIncrement(0.1);
    pde_test->TimeIncrement(0.1);
    pde_validate->AssignInitialConditions(region);
    pde_test->AssignInitialConditionsTest(region);
	_test(pde_validate->Transient() == true);
	_test(pde_test->Transient() == true);
	pde_validate->LateAccumulate(region);
	pde_test->LateAccumulateTest(region);

    for (size_t i(0); i < DOF_indexes.size(); ++i) {
      if (DOF_indexes[i] != NULL_IDX) {
        // _test conductance matrix
        for (size_t j(0); j < DOF_indexes.size(); ++j) {
          if (DOF_indexes[j] != NULL_IDX) {
            _test(pde_test->GetG()->At(DOF_indexes[i], DOF_indexes[j]) == pde_validate->GetG()->At(i, j));
          }
        }
        // check load vector
        _test((*pde_test->GetRH())[DOF_indexes[i]] == (*pde_validate->GetRH())[i]);
      }
    }

    delete pde_validate;
    delete pde_test;
  }



void PDE_Integrator_UoM_Test::TestTwoScalarVariables() {
    Region<2U>& region = model->Region("Model");
    region.RenumberNodes();
    
    GaussJordan_Solver solver;

    pde_validate = new PDE_Integrator_UoM_Mock<2U, Region>(solver);
    pde_test = new PDE_Integrator_UoM_Mock<2U, Region>(solver);

    pde_validate->Add(pressureLHS);
    pde_validate->Add(sourceVolume);
    pde_validate->Add(sourcePoint);
    pde_validate->Add(gravityTerm);
    pde_validate->Add(temperatureLHS);
    pde_validate->Add(sourceHeat);

    pde_test->Add(pressureLHS);
    pde_test->Add(sourceVolume);
    pde_test->Add(sourcePoint);
    pde_test->Add(gravityTerm);
    pde_test->Add(temperatureLHS);
    pde_test->Add(sourceHeat);

    const std::vector<size_t>& DOF_indexes = pde_test->GetDOFIndex();
    // 1. test matrix establish and enumerate DOFs
    pde_validate->EstablishMatrixSetup(region);
    pde_test->EstablishMatrixSetupTest(region);
    size_t dirchletDOFs(0);
    Index pressureKey = model->Database().StorageKey("fluid pressure");
    Index temperatureKey = model->Database().StorageKey("temperature");
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) == DIRICH) {
        dirchletDOFs += 1;
      }
      if ((*nIter)->Status(temperatureKey) == DIRICH) {
        dirchletDOFs += 1;
      }
    }

    _test((pde_validate->GetRH()->size()) == (pde_test->GetRH()->size()));
    pde_test->EnumerateAndFixMatrixSize(region);
    _test(pde_validate->GetRH()->size() == DOF_indexes.size());
    _test(pde_validate->GetRH()->size() == pde_test->GetRH()->size() + dirchletDOFs);
    _test(pde_validate->GetG()->Rows() == pde_test->GetG()->Rows() + dirchletDOFs);
    _test(pde_validate->GetG()->Cols() == pde_test->GetG()->Cols() + dirchletDOFs);

    // 2. test accumulate
    pde_validate->Accumulate(region);
    pde_test->AccumulateTest(region);

	// 3. test AssignInitialConditions & LateAccumulate
	pde_validate->TimeIncrement(0.1);
	pde_test->TimeIncrement(0.1);
	pde_validate->AssignInitialConditions(region);
	pde_test->AssignInitialConditionsTest(region);
	_test(pde_validate->Transient() == true);
	_test(pde_test->Transient() == true);
	pde_validate->LateAccumulate(region);
	pde_test->LateAccumulateTest(region);
    for (size_t i(0); i < DOF_indexes.size(); ++i) {
      if (DOF_indexes[i] != NULL_IDX) {
        // _test conductance matrix
        for (size_t j(0); j < DOF_indexes.size(); ++j) {
          if (DOF_indexes[j] != NULL_IDX) {
            _test(pde_test->GetG()->At(DOF_indexes[i], DOF_indexes[j]) == pde_validate->GetG()->At(i, j));
          }
        }
        // check load vector
        _test((*pde_test->GetRH())[DOF_indexes[i]] == (*pde_validate->GetRH())[i]);
      }
    }

    delete pde_validate;
    delete pde_test;
  }


  void PDE_Integrator_UoM_Test::TestOutputSingleVariable() {
    Region<2U>& region = model->Region("Model");
    region.RenumberNodes();
    Reset();

    Index pressureKey = model->Database().StorageKey("fluid pressure");
    Index pressureValKey = model->Database().StorageKey("fluid pressure previous");
    
    GaussJordan_Solver solver;

    pde_test = new PDE_Integrator_UoM_Mock<2U, Region>(solver);
    pde_test->Add(pressureLHS);
    pde_test->Add(sourceVolume);    
    pde_test->EstablishMatrixSetupTest(region);
    pde_test->EnumerateAndFixMatrixSize(region);
    std::map<size_t, double64> result;
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) != DIRICH) {
        result[(*nIter)->Idx()] = (*nIter)->Read(pressureKey);
      } 
    }

    std::vector<double64>* valid_x = pde_test->GetX();
    size_t idx(0);
    for (auto& it : result) {
      (*valid_x)[idx] = it.second;
      ++idx;
    }

    pde_test->OutputResultsTest(region);
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      _test((*nIter)->Read(pressureKey) == (*nIter)->Read(pressureValKey));     
    }
    
    delete pde_test;
  }

} // end namespace csmp
