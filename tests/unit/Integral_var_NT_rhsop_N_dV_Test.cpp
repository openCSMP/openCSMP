#include "Integral_var_NT_rhsop_N_dV_Test.h"

using namespace std;

namespace csmp {

Integral_var_NT_rhsop_N_dV_Test::Integral_var_NT_rhsop_N_dV_Test( bool verbose )
 : tol_(0.001), verbose_(verbose), sg_(nullptr)
{
  const bool       isoparametric(true);
  ANSYS_Interface  mesh_interface(isoparametric);  // true = isoparametric elements
  VSet<2U>         mesh_container;
  ModelTopology    mesh_topology(isoparametric);   // true = isoparametric elements

  // Building Region object from ANSYS data files
  if ( verbose_ ) cout <<"Reading mesh..."<<endl;
  string mesh_name("pde_integrator_test");
  const bool binary_file( true );
  mesh_interface.Read_ANSYS_Mesh( mesh_name.c_str(), mesh_container, mesh_topology, binary_file, true );
  if ( verbose_ ) {
      cout <<"Finished reading mesh..."<<endl;
      cout <<"Building Model..."<<endl;
    }
  sg_= new Model<2U> ( mesh_topology, mesh_container, "CSMP-2phase-variables.txt");

    // Set values on nodes
  sg_->InputPropertyValue("fluid pressure", makeScalar(PLAIN,1.));
  sg_->InputPropertyValue("diffusivity", makeScalar(PLAIN,1.));
  sg_->InputPropertyValue("permeability", makeScalar(PLAIN,1.));
  sg_->InputPropertyValue("total mobility", makeScalar(PLAIN,1.));
}


Integral_var_NT_rhsop_N_dV_Test::~Integral_var_NT_rhsop_N_dV_Test()
{
    delete sg_;
}

void Integral_var_NT_rhsop_N_dV_Test::run() {
  valueTest();
  compareConsistentTest();
  compareLumpedTest();
}

void Integral_var_NT_rhsop_N_dV_Test::valueTest() {
  // create object
  Integral_var_NT_rhsop_N_dV<2U,Element<2U> > integral(sg_->Database(),
                                                    "diffusivity",
                                                    "fluid pressure",
                                                    "fluid pressure",
                                                    "total mobility",
                                                    2.0);

  const size_t dof = sg_->Region("Model").Nodes();
                                                    
    // Set values
    vector<double> mobility;
    mobility.push_back(1.0);
    mobility.push_back(2.0);
    mobility.push_back(3.0);
    mobility.push_back(4.0);
    mobility.push_back(5.0);
    setNodeVariable(mobility, "mobility");
    
    vector<double> conductivity;
    conductivity.push_back(1.0);
    conductivity.push_back(2.0);
    conductivity.push_back(3.0);
    conductivity.push_back(4.0);
    setElementVariable(conductivity, "conductivity1");
    
        vector<double> pressure;
    pressure.push_back(1.0);
    pressure.push_back(1.0);
    pressure.push_back(1.0);
    pressure.push_back(1.0);
    pressure.push_back(1.0);
    setNodeVariable(pressure, "pressure");
    
    std::vector<double> rhs_vec(dof);
    calculateGlobalRHS(rhs_vec, integral);

    /* Matlab output
    (output from assemblyNonlinear.m, nonlinMat.m)
    1.0e+005 *

    1.2333    0.1833         0    0.5625    0.8333
    0.1833    1.3333    0.5000         0    0.8167
         0    0.5000    1.5167    0.2375    0.8500
    0.5625         0    0.2375    1.8500    0.9750
    0.8333    0.8167    0.8500    0.9750    3.9833
    
    A*p
    2.8124
    2.8333
    3.1042
    3.6250
    7.4583
    
   */

  _equal(rhs_vec[0], 2.0*2.8124e+05, tol_*1.0e+05);
  _equal(rhs_vec[1], 2.0*2.8333e+05, tol_*1.0e+05);
  _equal(rhs_vec[2], 2.0*3.1042e+05, tol_*1.0e+05);
  _equal(rhs_vec[3], 2.0*3.6250e+05, tol_*1.0e+05);
  _equal(rhs_vec[4], 2.0*7.4583e+05, tol_*1.0e+05);
}

void Integral_var_NT_rhsop_N_dV_Test::compareConsistentTest() {
  compareTest(false);
}

void Integral_var_NT_rhsop_N_dV_Test::compareLumpedTest() {
  compareTest(true);
}

void Integral_var_NT_rhsop_N_dV_Test::compareTest(bool lumped) {
  Integral_var_NT_lhsop_N_dV<2U,Element<2U> > lhs(sg_->Database(),
                                     "diffusivity",
                                     "fluid pressure",
                                     "fluid pressure",
                                     "total mobility");
  lhs.LumpedFormulation(lumped);
  Integral_var_NT_rhsop_N_dV<2U,Element<2U> > rhs(sg_->Database(),
                                     "diffusivity",
                                     "fluid pressure",
                                     "fluid pressure",
                                     "total mobility");
  rhs.LumpedFormulation(lumped);
  const size_t dof = sg_->Region("Model").Nodes();
                                              
// Set values
    vector<double> mobility;
    mobility.push_back(1.0);
    mobility.push_back(2.0);
    mobility.push_back(3.0);
    mobility.push_back(4.0);
    mobility.push_back(5.0);
    setNodeVariable(mobility, "total mobility");
    
    vector<double> conductivity;
    conductivity.push_back(1.0);
    conductivity.push_back(2.0);
    conductivity.push_back(3.0);
    conductivity.push_back(4.0);
    setElementVariable(conductivity, "diffusivity");
    
    vector<double> pressure;
    pressure.push_back(0.0);
    pressure.push_back(0.0);
    pressure.push_back(1.0);
    pressure.push_back(1.0);
    pressure.push_back(0.5);
    setNodeVariable(pressure, "fluid pressure");
    
    // Calculate lhs
    SparseMatrix sm;
    calculateGlobalMatrix(sm, lhs);
    std::vector<double> lhs_vec(dof);
    std::fill(lhs_vec.begin(), lhs_vec.end(), 0.0);
    for (size_t i = 0; i < dof; ++i) {
      for (size_t j = 0; j < dof; ++j) {
        lhs_vec[i] += sm(i, j) * pressure[j];
      }
    }
    
  // Calculate rhs
    std::vector<double> rhs_vec(dof);
    calculateGlobalRHS(rhs_vec, rhs);
    
    // Test for equality
    assert(lhs_vec.size() == rhs_vec.size());
    for (size_t i = 0; i < lhs_vec.size(); ++i) {
      _equal(lhs_vec[i], rhs_vec[i], tol_);
    } 
}


void Integral_var_NT_rhsop_N_dV_Test::showNodeVariable(const char* var_name) {
    //std::deque<Node<2U> >::iterator it;
    csmp::Index key(sg_->Database().StorageKey(var_name));
    ScalarVariable sc;
        
    unsigned int i = 0;
    for (vector<Node<2U>*>::const_iterator it = sg_->Region("Model").NodesBegin(); it != sg_->Region("Model").NodesEnd(); ++it) {
        (*it)->Read(key,sc);
        if ( verbose_ ) cout << "Node " << i << ": " << sc() << endl;
        ++i;
    }
}

void Integral_var_NT_rhsop_N_dV_Test::setNodeVariable(vector<double>& var, const char* var_name) {
    //std::deque<Node<2U> >::iterator it;
    csmp::Index key(sg_->Database().StorageKey(var_name));

    unsigned int i = 0;
    for (vector<Node<2U>*>::iterator it = sg_->Region("Model").NodesBegin(); it != sg_->Region("Model").NodesEnd(); ++it) {
        (*it)->Store(key,ScalarVariable(PLAIN, static_cast<double>(var[i])));
        ++i;
    }
}

void Integral_var_NT_rhsop_N_dV_Test::setElementVariable(vector<double>& var, const char* var_name) {
  //std::deque<Element<2U> >::iterator it;
  csmp::Index key(sg_->Database().StorageKey(var_name));
  
  unsigned int i = 0;
  for (vector<Element<2U>*>::const_iterator it = sg_->Region("Model").ElementsBegin(); it != sg_->Region("Model").ElementsEnd(); ++it) {
    (*it)->Store(key,ScalarVariable(PLAIN, static_cast<double>(var[i])));
              
    ++i;
  }
}

void Integral_var_NT_rhsop_N_dV_Test::calculateGlobalMatrix(SparseMatrix& sm, MathOperatorLHS<2U>& oper) {
    sm.Zero();
    sm.Resize(sg_->Region("Model").Nodes());
    
    //std::deque<Element<2U> >::iterator it;
    for (vector<Element<2U>*>::const_iterator it = sg_->Region("Model").ElementsBegin(); it != sg_->Region("Model").ElementsEnd(); ++it) {
        oper.GetOperands( *(*it));
        oper.ComputeContribution(*(*it));
        oper.AssignToGlobal(*(*it), sm);
    }
    
}

void Integral_var_NT_rhsop_N_dV_Test::calculateGlobalRHS(vector<double>& rhs, MathOperatorRHS<2U>& oper) {
    //std::deque<Element<2U> >::iterator it;
    for (vector<Element<2U>*>::const_iterator it = sg_->Region("Model").ElementsBegin(); it != sg_->Region("Model").ElementsEnd(); ++it) {
        oper.GetOperands(*(*it));
        oper.ComputeContribution(*(*it));
        oper.AssignToGlobal(*(*it), rhs);
    }
    
}

} // csmp
