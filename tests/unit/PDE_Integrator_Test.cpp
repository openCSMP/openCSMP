#include "PDE_Integrator_Test.h"
#include "PDE_Integrator.h"

#include "VSet.h"
#include "ModelTopology.h"
#include "vsetMakers.h"

#include "GaussJordan_Solver.h"
#include "Model.h"
#include "Region.h"

// TODO: distinguish this UNIT TEST from an INTEGRATION TEST
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PointSource_rhsop.h"
#include "NumIntegral_dNT_op_dV.h"
#include "VelocityAndVolumeFlux.h"
    // thermal equation
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"

using namespace std;

namespace csmp {

PDE_Integrator_Test::PDE_Integrator_Test( Model<2U>& model )
    : model_{model}
  {
  }

  PDE_Integrator_Test::~PDE_Integrator_Test()
  {
  }




void PDE_Integrator_Test::TestAssembly()
 {
    // 0. creates test model with 4 elements and 6 Face objects for the box boundaries
    VSet<2U> vset;
    ModelTopology topo = test_CreateSimplestPolyElement2DModel( vset );
    const bool treat_regions_as_boundaries{false};
    Model<2U> model( topo, vset, "CSMP-1phase-variables.txt", treat_regions_as_boundaries );
    model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1.) );
    
    // 1. setting up PDE_Integrator for simple case with a single scalar variable
    // --------------------------------------------------------------------------
    GaussJordan_Solver        solver;
    PDE_Integrator<2U,Region> pde_integrator(solver);

    // Create pde-operators
    LHS_FixedValueMatrix<2U>  lhs( model.Database(), "permeability", "fluid pressure", "fluid pressure", 1. );
    RHS_FixedValueMatrix<2U>  rhs( model.Database(), "fluid volume source", "fluid pressure", 1. );

    // assign them to pde integrator
    pde_integrator.Add(&lhs);
    pde_integrator.Add(&rhs);
    
    const bool debug{true};
    model.Apply( pde_integrator, debug );
    
    // compare matrix with expected matrix
    // TODO: do testing here; not sure how to get to matrix and vector
    
    
    // 2. setting up PDE_Integrator for case with a single vector solution variable
    // ----------------------------------------------------------------------------
    
    
    // 3. setting up PDE_Integrator for coupled system of 2 scalars
    // ----------------------------------------------------------------------------
    

    // 4. setting up PDE_Integrator for coupled system of 1 scalar + 1 vector variable
    // -------------------------------------------------------------------------------
    
    cout <<"\nPDE_Integrator_Test::TestAssembly: finished test."<< endl;
    
 } // end TestAssembly


    // Create pde-operators
//    NumIntegral_dNT_op_dN_dV<2U> pressureLHS(model_.Database(), "permeability", "fluid pressure", "fluid pressure");
//    NumIntegral_NT_op_N_dV<2U> sourceVolume(model_.Database(), "fluid volume source", "fluid pressure");
//    pde_integrator.Add(&pressureLHS);
//    pde_integrator.Add(&sourceVolume);





void PDE_Integrator_Test::Reset()
  {
    Region<2U>& region = model_.Region("Model");
    size_t elementNum = region.Cells();

    // Set values on nodes and elements
    Index pressureKey = model_.Database().StorageKey("fluid pressure");
    Index pressureValKey = model_.Database().StorageKey("fluid pressure previous");
    Index temperatureKey = model_.Database().StorageKey("temperature");
    Index fluidVolumeKey = model_.Database().StorageKey("fluid volume source");
    Index fluidNodeSourceKey = model_.Database().StorageKey("nodal fluid point source");
    Index heatSourceKey = model_.Database().StorageKey("thermal volume source");
    Index gravityVectorKey = model_.Database().StorageKey("gravity vector");

    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      (*nIter)->Store(pressureKey, makeScalar(PLAIN, (*nIter)->Idx()));
      (*nIter)->Store(pressureValKey, makeScalar(PLAIN, (*nIter)->Idx()));
      (*nIter)->Store(temperatureKey, makeScalar(PLAIN, (*nIter)->Idx()));
      (*nIter)->Store(fluidNodeSourceKey, makeScalar(PLAIN, (*nIter)->Idx()));
    }

    for (auto eIter = region.CellsBegin(); eIter != region.CellsEnd(); ++eIter) {
      (*eIter)->Store(fluidVolumeKey, makeScalar(PLAIN, (*eIter)->Idx()));
      (*eIter)->Store(heatSourceKey, makeScalar(PLAIN, (*eIter)->Idx()) / elementNum);
      (*eIter)->Store(gravityVectorKey, makeVector(PLAIN, PLAIN, (*eIter)->BaryCenter()[0], (*eIter)->BaryCenter()[1]));
    }


    model_.InputPropertyValue("permeability", makeScalar(PLAIN, 1.));
    model_.InputPropertyValue("total mobility", makeScalar(PLAIN, 1.)); // conductivity
                                                                        // Set boundary conditions
    model_.InputBoundaryValue(BOTTOM, "fluid pressure", makeScalar(DIRICH, 1.0));
    model_.InputBoundaryValue(BOTTOM, "fluid pressure previous", makeScalar(DIRICH, 1.0));
    model_.InputBoundaryValue(BOTTOM, "temperature", makeScalar(DIRICH, 1.0));

  
  }

  void PDE_Integrator_Test::run()
  {
    //=======================================
    // test single variable
    //=======================================
    TestAssembly();

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

/*
  void PDE_Integrator_Test::TestSingleVariable() {
    // 0. prepare
    Region<2U>& region = model_.Region("Model");
    region.RenumberNodes();

    GaussJordan_Solver solver;

    // Create pde-operators
    NumIntegral_dNT_op_dN_dV<2U> pressureLHS(model_.Database(),
      "permeability",
      "fluid pressure",
      "fluid pressure");

    NumIntegral_NT_op_N_dV<2U> sourceVolume(model_.Database(),
      "fluid volume source", "fluid pressure");

    PointSource_rhsop<2U> sourcePoint(model_.Database(),
      "nodal fluid point source", "fluid pressure");
    
    NumIntegral_dNT_op_dV<2U> gravityTerm(model_.Database(),
      "gravity vector",
      "fluid pressure");

    //fluid_velocity = new VelocityAndVolumeFlux<2U, Element<2U>>(*model, "total mobility", "porosity", "fluid pressure", true, "velocity");

    NumIntegral_dNT_op_dN_dV<2U> temperatureLHS(model_.Database(),
      "permeability",
      "temperature",
      "temperature");

    NumIntegral_NT_op_N_dV<2U> sourceHeat(model_.Database(),
      "thermal volume source", "temperature");


    pde_reference_ = new PDE_Integrator<2U,Region>(solver);
    pde_test_      = new PDE_Integrator<2U,Region>(solver);

    pde_reference_->Add(&pressureLHS);
    pde_reference_->Add(&sourceVolume);
    pde_reference_->Add(&sourcePoint);
    pde_reference_->Add(&gravityTerm);
    pde_test_->Add(&pressureLHS);
    pde_test_->Add(&sourceVolume);
    pde_test_->Add(&sourcePoint);
    pde_test_->Add(&gravityTerm);

    const vector<size_t>& DOF_indexes = pde_test_->DOF_indexes_;

    // 1. test matrix establish and enumerate DOFs
    pde_reference_->EstablishMatrixSetup(region);
    pde_test_->EstablishMatrixSetup(region);

    size_t pressureDirchletDOFs(0);
    Index pressureKey = model_.Database().StorageKey("fluid pressure");
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) == DIRICH) {
        pressureDirchletDOFs += 1;
      }
    }

    _test((pde_reference_->GetRH()->size()) == (pde_test_->GetRH()->size()));

    pde_test_->EnumerateAndFixMatrixSize(region);
    _test(pde_reference_->GetRH()->size() == DOF_indexes.size());
    _test(pde_reference_->GetRH()->size() == pde_test_->GetRH()->size() + pressureDirchletDOFs);
    _test(pde_reference_->GetG()->Rows() == pde_test_->GetG()->Rows() + pressureDirchletDOFs);
    _test(pde_reference_->GetG()->Cols() == pde_test_->GetG()->Cols() + pressureDirchletDOFs);

    // 2. test accumulate
    pde_reference_->Accumulate(region);
    pde_test_->Accumulate(region);

    // 3. test AssignInitialConditions & LateAccumulate
    pde_reference_->TimeIncrement(0.1);
    pde_test_->TimeIncrement(0.1);
    pde_reference_->AssignInitialConditions(region);
    pde_test_->AssignInitialConditions(region);
	_test(pde_reference_->Transient() == true);
	_test(pde_test_->Transient() == true);
	pde_reference_->LateAccumulate(region);
	pde_test_->LateAccumulate(region);

    for (size_t i(0); i < DOF_indexes.size(); ++i) {
      if (DOF_indexes[i] != NULL_IDX) {
        // _test conductance matrix
        for (size_t j(0); j < DOF_indexes.size(); ++j) {
          if (DOF_indexes[j] != NULL_IDX) {
            _test(pde_test_->GetG()->At(DOF_indexes[i], DOF_indexes[j]) == pde_reference_->GetG()->At(i, j));
          }
        }
        // check load vector
        _test((*pde_test_->GetRH())[DOF_indexes[i]] == (*pde_reference_->GetRH())[i]);
      }
    }

    delete pde_reference_;
    delete pde_test;

  } // end method
*/



#if 0

void PDE_Integrator_Test::TestTwoScalarVariables() {
    Region<2U>& region = model_.Region("Model");
    region.RenumberNodes();
    
    GaussJordan_Solver solver;

    pde_reference_ = new PDE_Integrator_<2U, Region>(solver);
    pde_test = new PDE_Integrator<2U, Region>(solver);

    pde_reference_->Add(pressureLHS);
    pde_reference_->Add(sourceVolume);
    pde_reference_->Add(sourcePoint);
    pde_reference_->Add(gravityTerm);
    pde_reference_->Add(temperatureLHS);
    pde_reference_->Add(sourceHeat);

    pde_test_->Add(pressureLHS);
    pde_test_->Add(sourceVolume);
    pde_test_->Add(sourcePoint);
    pde_test_->Add(gravityTerm);
    pde_test_->Add(temperatureLHS);
    pde_test_->Add(sourceHeat);

    const std::vector<size_t>& DOF_indexes = pde_test_->GetDOFIndex();
    // 1. test matrix establish and enumerate DOFs
    pde_reference_->EstablishMatrixSetup(region);
    pde_test_->EstablishMatrixSetupTest(region);
    size_t dirchletDOFs(0);
    Index pressureKey = model_.Database().StorageKey("fluid pressure");
    Index temperatureKey = model_.Database().StorageKey("temperature");
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) == DIRICH) {
        dirchletDOFs += 1;
      }
      if ((*nIter)->Status(temperatureKey) == DIRICH) {
        dirchletDOFs += 1;
      }
    }

    _test((pde_reference_->GetRH()->size()) == (pde_test_->GetRH()->size()));
    pde_test_->EnumerateAndFixMatrixSize(region);
    _test(pde_reference_->GetRH()->size() == DOF_indexes.size());
    _test(pde_reference_->GetRH()->size() == pde_test_->GetRH()->size() + dirchletDOFs);
    _test(pde_reference_->GetG()->Rows() == pde_test_->GetG()->Rows() + dirchletDOFs);
    _test(pde_reference_->GetG()->Cols() == pde_test_->GetG()->Cols() + dirchletDOFs);

    // 2. test accumulate
    pde_reference_->Accumulate(region);
    pde_test_->AccumulateTest(region);

	// 3. test AssignInitialConditions & LateAccumulate
	pde_reference_->TimeIncrement(0.1);
	pde_test_->TimeIncrement(0.1);
	pde_reference_->AssignInitialConditions(region);
	pde_test_->AssignInitialConditionsTest(region);
	_test(pde_reference_->Transient() == true);
	_test(pde_test_->Transient() == true);
	pde_reference_->LateAccumulate(region);
	pde_test_->LateAccumulateTest(region);
    for (size_t i(0); i < DOF_indexes.size(); ++i) {
      if (DOF_indexes[i] != NULL_IDX) {
        // _test conductance matrix
        for (size_t j(0); j < DOF_indexes.size(); ++j) {
          if (DOF_indexes[j] != NULL_IDX) {
            _test(pde_test_->GetG()->At(DOF_indexes[i], DOF_indexes[j]) == pde_reference_->GetG()->At(i, j));
          }
        }
        // check load vector
        _test((*pde_test_->GetRH())[DOF_indexes[i]] == (*pde_reference_->GetRH())[i]);
      }
    }

    delete pde_reference_;
    delete pde_test;
  }


  void PDE_Integrator_Test::TestOutputSingleVariable() {
    Region<2U>& region = model_.Region("Model");
    region.RenumberNodes();
    Reset();

    Index pressureKey = model_.Database().StorageKey("fluid pressure");
    Index pressureValKey = model_.Database().StorageKey("fluid pressure previous");
    
    GaussJordan_Solver solver;

    pde_test = new PDE_Integrator_UoM_Mock<2U, Region>(solver);
    pde_test_->Add(pressureLHS);
    pde_test_->Add(sourceVolume);
    pde_test_->EstablishMatrixSetupTest(region);
    pde_test_->EnumerateAndFixMatrixSize(region);
    std::map<size_t, double> result;
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      if ((*nIter)->Status(pressureKey) != DIRICH) {
        result[(*nIter)->Idx()] = (*nIter)->Read(pressureKey);
      } 
    }

    std::vector<double>* valid_x = pde_test_->GetX();
    size_t idx(0);
    for (auto& it : result) {
      (*valid_x)[idx] = it.second;
      ++idx;
    }

    pde_test_->OutputResultsTest(region);
    for (auto nIter = region.NodesBegin(); nIter != region.NodesEnd(); ++nIter) {
      _test((*nIter)->Read(pressureKey) == (*nIter)->Read(pressureValKey));     
    }
    
    delete pde_test;
  }

#endif

} // end namespace csmp
