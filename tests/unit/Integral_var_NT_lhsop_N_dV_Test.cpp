#include "Integral_var_NT_lhsop_N_dV_Test.h"
#include "vsetMakers.h"
#include "VSetConverter.h"

using namespace std;

namespace csmp {

Integral_var_NT_lhsop_N_dV_Test::Integral_var_NT_lhsop_N_dV_Test( bool verbose ) : tol_(0.001), verbose_(verbose)
  {
    VSet<2U>    mesh_container;
    test_Create_TrianglePatch_VSet( mesh_container );

    // Building Region object from ANSYS data files
    string mesh_name("triangle_patch");
    if ( verbose_ ) cout <<"\nIntegral_var_NT_lhsop_N_dV_Test: Building Model..."<<endl;
    sg_= new Model<2U>( mesh_container, "CSMP-2phase-variables.txt" );

    // Set values on nodes
    sg_->InputPropertyValue("fluid pressure", makeScalar(PLAIN,1.));
    sg_->InputPropertyValue("diffusivity", makeScalar(PLAIN,1.));
    sg_->InputPropertyValue("permeability", makeScalar(PLAIN,1.));
    sg_->InputPropertyValue("nodal fluid volume source", makeScalar(PLAIN,1.));
  }



  Integral_var_NT_lhsop_N_dV_Test::~Integral_var_NT_lhsop_N_dV_Test()
  {
      if (sg_!=NULL)  delete sg_;
  }


  void Integral_var_NT_lhsop_N_dV_Test::run() {
    valueTest();
    compareConsistentTest();
    compareLumpedTest();
    lumpedTest();
    rowSumTest();
  }

  void Integral_var_NT_lhsop_N_dV_Test::valueTest() {
    // create object
    Integral_var_NT_lhsop_N_dV<2U,Element<2U> > integral(sg_->Database(),
                                                     "diffusivity",
                                                     "fluid pressure",
                                                     "fluid pressure",
                                                     "nodal fluid volume source");

    // Set values
    vector<double> mobility;
    mobility.push_back(1.0);
    mobility.push_back(2.0);
    mobility.push_back(3.0);
    mobility.push_back(4.0);
    mobility.push_back(5.0);
    setNodeVariable(mobility, "nodal fluid volume source");

    vector<double> conductivity;
    conductivity.push_back(1.0);
    conductivity.push_back(2.0);
    conductivity.push_back(3.0);
    conductivity.push_back(4.0);
    setElementVariable(conductivity, "diffusivity");

    SparseMatrix sm;
    calculateGlobalMatrix(sm, integral);

    /* Matlab output
    (output from assemblyNonlinear.m, nonlinMat.m)
    1.0e+005 *

    1.2333    0.1833         0    0.5625    0.8333
    0.1833    1.3333    0.5000         0    0.8167
         0    0.5000    1.5167    0.2375    0.8500
    0.5625         0    0.2375    1.8500    0.9750
    0.8333    0.8167    0.8500    0.9750    3.9833
   */

    _equal(sm.At(0,0), 1.2333*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(0,1), 0.1833*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(0,2), 0             , tol_*1.0e+04);
    _equal(sm.At(0,3), 0.5625*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(0,4), 0.8333*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(1,0), 0.1833*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(1,1), 1.3333*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(1,2), 0.5000*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(1,3), 0             , tol_*1.0e+04);
    _equal(sm.At(1,4), 0.8167*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(2,0), 0             , tol_*1.0e+04);
    _equal(sm.At(2,1), 0.5000*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(2,2), 1.5167*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(2,3), 0.2375*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(2,4), 0.8500*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(3,0), 0.5625*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(3,1), 0             , tol_*1.0e+04);
    _equal(sm.At(3,2), 0.2375*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(3,3), 1.8500*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(3,4), 0.9750*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(4,0), 0.8333*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(4,1), 0.8167*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(4,2), 0.8500*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(4,3), 0.9750*1.0e+05, tol_*1.0e+04);
    _equal(sm.At(4,4), 3.9833*1.0e+05, tol_*1.0e+04);

    //sm.Out();
  }

  void Integral_var_NT_lhsop_N_dV_Test::compareConsistentTest() {
    compareTest(false);
  }

  void Integral_var_NT_lhsop_N_dV_Test::compareLumpedTest() {
    compareTest(true);
  }

  void Integral_var_NT_lhsop_N_dV_Test::compareTest(bool lumped) {
    Integral_var_NT_lhsop_N_dV<2U,Element<2U> > integral(sg_->Database(),
                                            "diffusivity",
                                            "fluid pressure",
                                            "fluid pressure",
                                            "nodal fluid volume source");
    integral.LumpedFormulation(lumped);
    Integral_NT_lhsop_N_dV<2U,Element<2U> > simple(sg_->Database(),
                                      "diffusivity",
                                      "fluid pressure",
                                      "fluid pressure");
    simple.LumpedFormulation(lumped);

    // Set values
    vector<double> mobility;
    mobility.push_back(1.0);
    mobility.push_back(1.0);
    mobility.push_back(1.0);
    mobility.push_back(1.0);
    mobility.push_back(1.0);
    setNodeVariable(mobility, "total mobility");
    
    vector<double> conductivity;
    conductivity.push_back(1.0);
    conductivity.push_back(2.0);
    conductivity.push_back(3.0);
    conductivity.push_back(4.0);
    setElementVariable(conductivity, "diffusivity");
    
    SparseMatrix sm_int;
    SparseMatrix sm_simp;
    calculateGlobalMatrix(sm_int, integral);
    calculateGlobalMatrix(sm_simp, simple);
    
    for (auto i = 0; i < sm_int.Rows(); ++i) {
      for (size_t j = 0; j < sm_int.Cols(); ++j) {
        _equal(sm_int.At(i, j), sm_simp.At(i, j), tol_);
      }
    }
  }

  void Integral_var_NT_lhsop_N_dV_Test::lumpedTest() {
    Integral_var_NT_lhsop_N_dV<2U,Element<2U> > lumped(sg_->Database(),
                                                   "diffusivity",
                                                   "fluid pressure",
                                                   "fluid pressure",
                                                   "nodal fluid volume source");
    lumped.LumpedFormulation(true);

    // Set values
    vector<double> mobility;
    mobility.push_back(1.0);
    mobility.push_back(1.0);
    mobility.push_back(2.0);
    mobility.push_back(2.0);
    mobility.push_back(5.0);
    setNodeVariable(mobility, "nodal fluid volume source");
    
    vector<double> conductivity;
    conductivity.push_back(1.0);
    conductivity.push_back(2.0);
    conductivity.push_back(3.0);
    conductivity.push_back(4.0);
    setElementVariable(conductivity, "diffusivity");
    
    /* Matlab output (assemblyLumped.m)
  1.0e+005 *

    2.3542         0         0         0         0
         0    2.1667         0         0         0
         0         0    2.3542         0         0
         0         0         0    2.5625         0
         0         0         0         0    6.4792
   */

    SparseMatrix sm;
    calculateGlobalMatrix(sm, lumped);

    _equal(sm.At(0,0), 2.3542e+05, tol_*1.0e+04);
    _equal(sm.At(1,1), 2.1667e+05, tol_*1.0e+04);
    _equal(sm.At(2,2), 2.3542e+05, tol_*1.0e+04);
    _equal(sm.At(3,3), 2.5625e+05, tol_*1.0e+04);
    _equal(sm.At(4,4), 6.4792e+05, tol_*1.0e+04);

  }

  void Integral_var_NT_lhsop_N_dV_Test::rowSumTest() {
    // create object
    Integral_var_NT_lhsop_N_dV<2U,Element<2U> > lumped(sg_->Database(),
                                                   "diffusivity",
                                                   "fluid pressure",
                                                   "fluid pressure",
                                                   "nodal fluid volume source");
    lumped.LumpedFormulation(true);
    Integral_var_NT_lhsop_N_dV<2U,Element<2U> > consistent(sg_->Database(),
                                                       "diffusivity",
                                                       "fluid pressure",
                                                       "fluid pressure",
                                                       "nodal fluid volume source");
    consistent.LumpedFormulation(false);

    const size_t dof = sg_->Region("Model").Nodes();

    SparseMatrix sm_lumped;
    SparseMatrix sm_consistent;
    calculateGlobalMatrix(sm_lumped, lumped);
    calculateGlobalMatrix(sm_consistent, consistent);

    std::vector<double> vec(dof);
    std::fill(vec.begin(), vec.end(), 0.0);
    for (auto i = 0; i < dof; ++i) {
      for (size_t j = 0; j < dof; ++j) {
        vec[i] += sm_consistent.At(i,j);
      }
    }

    for (auto i = 0; i < dof; ++i) {
      _equal(sm_lumped.At(i,i), vec[i], tol_*1.0e+03);
      for (size_t j = 0; j < dof; ++j) {
        if (i != j) _equal(sm_lumped.At(i, j), 0.0, tol_);
      }
    }

    //cout << "Lumped formulation" << endl;
    //sm_lumped.Out();
    //cout << "Consistent formulation" << endl;
    //sm_consistent.Out();
    //cout << "Diagonal sum" << endl;
    //for (auto i = 0; i < dof; ++i) {
    //  cout << vec[i] << endl;
    //}
  }

  void Integral_var_NT_lhsop_N_dV_Test::showNodeVariable(const char* var_name) {
    const csmp::Index key(sg_->Database().StorageKey(var_name));
    unsigned int i = 0;
    Region<2>& model_domain( sg_->Region("Model") );
    for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it) {
        if ( verbose_ ) cout << "Node " << i << ": " << (*it)->Read(key) << endl;
        ++i;
      }
  }

  void Integral_var_NT_lhsop_N_dV_Test::setNodeVariable(vector<double>& var, const char* var_name) {
    const csmp::Index key(sg_->Database().StorageKey(var_name));
    unsigned int i = 0;
    Region<2>& model_domain( sg_->Region("Model") );
    // TODO: Test fails here because there are many more nodes in the model than supplied to setNodeVariable
    assert( var.size() == model_domain.Nodes() );
    for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it) {
        (*it)->Store( key, makeScalar(PLAIN, var[i++]) );
      }
  }

  void Integral_var_NT_lhsop_N_dV_Test::setElementVariable(vector<double>& var, const char* var_name)
   {
    const csmp::Index key(sg_->Database().StorageKey(var_name));
    unsigned int i = 0;
    Region<2>& model_domain( sg_->Region("Model") );
    assert( var.size() == model_domain.Elements() );
    for ( auto it = model_domain.ElementsBegin(); it != model_domain.ElementsEnd(); ++it) {
       (*it)->Store(key,ScalarVariable(PLAIN, var[i++]) );
      }
  }

  void Integral_var_NT_lhsop_N_dV_Test::calculateGlobalMatrix( SparseMatrix& sm, MathOperatorLHS<2U>& oper ) {
    Region<2>& model_domain( sg_->Region("Model") );
    sm.Resize( model_domain.Nodes() );
    sm.Zero();
    
    for ( auto it = model_domain.ElementsBegin(); it != model_domain.ElementsEnd(); ++it ) {
        oper.GetOperands( *(*it));
        oper.ComputeContribution(*(*it));
        oper.AssignToGlobal(*(*it), sm);
      }
  }

} // csmp
