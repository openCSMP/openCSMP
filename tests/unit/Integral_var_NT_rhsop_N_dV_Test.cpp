#include "Integral_var_NT_rhsop_N_dV_Test.h"
#include "TRIANGLE_Interface.h"

using namespace std;

namespace csmp {

/**
  Test by Adrian Burri (ETHZ, 2004), by comparison with Matlab model.
  Uses Triangle iso.1 mesh called 'triangle_patch" as input. Mesh is a square with a central node.
  Mesh consists out of 4 triangles.
 */
Integral_var_NT_rhsop_N_dV_Test::Integral_var_NT_rhsop_N_dV_Test( bool verbose )
 : tol_(0.001), verbose_(verbose), sg_(nullptr)
{
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

  // Building Region object from ANSYS data files
  const string mesh_name("triangle_patch");
  if ( verbose_ ) cout <<"\nIntegral_var_NT_rhsop_N_dV_Test: Building Model..."<<endl;
  sg_= new Model<2U>( mesh_container, "CSMP-2phase-variables.txt" );
  sg_->Name(mesh_name.c_str());

  sg_->InputPropertyValue("permeability", makeScalar(PLAIN,1.));

  sg_->CreateProperty("conductivity", "HK", "m/s", SCALAR, ELEMENT ); // scalar-element
  sg_->InputPropertyValue("conductivity", makeScalar(PLAIN,1.));

  sg_->InputPropertyValue("fluid pressure", makeScalar(PLAIN,1.));
  
  sg_->CreateProperty("mobility", "lambda","none"); // default scalar-node
  sg_->InputPropertyValue("mobility", makeScalar(PLAIN,1.));
  sg_->Out();
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
  Integral_var_NT_rhsop_N_dV<2U> integral(sg_->Database(),
                                                    "conductivity",
                                                    "fluid pressure",
                                                    "fluid pressure",
                                                    "mobility",
                                                    2.0);

    setNodeVariable(    { 1.0, 2.0, 3.0, 4.0, 5.0 }, "mobility"       );
    setElementVariable( { 1.0, 2.0, 3.0, 4.0 },       "conductivity"   );
    setNodeVariable(    { 1.0, 1.0, 1.0, 1.0, 1.0 }, "fluid pressure" );

    const size_t dof{ sg_->Region("Model").Nodes() };
    std::vector<double> rhs_vec( dof, 0.0 );
    calculateGlobalRHS( rhs_vec, integral );

    if ( verbose_ )
      {
        cout << "\nvalueTest: RHS vector:\n";
        for ( size_t i{0U}; i < dof; ++i )
          cout << "  RHS[" << i << "] = " << rhs_vec[i] << "\n";
      }

    // -----------------------------------------------------------------------
    // Analytical reference values
    // -----------------------------------------------------------------------
    // With pressure={1,1,1,1,1}, RHS[j] = prefactor * sum_k M_jk
    // where M_jk is the variable mass matrix with conductivity={1,2,3,4}
    // and mobility={1,2,3,4,5}.
    //
    // From Integral_var_NT_lhsop_N_dV_Test::valueTest we know the assembled
    // matrix for these exact inputs. The row sums are:
    //   row 0: sum of (0,0),(0,1),(0,3),(0,4) entries
    //   row 1: sum of (1,0),(1,1),(1,2),(1,4) entries
    //   etc.
    //
    // Using M_jk^(e) = E_e * (A/60) * (2*op_j + 2*op_k + op_l):
    // c = A/60 = 62500/60 = 1041.667
    //
    // Row sum for node j = sum over all connected k of M_jk
    // These were computed in Integral_var_NT_lhsop_N_dV_Test::valueTest:
    //   row 0: M_00 + M_01 + M_03 + M_04
    //   row 1: M_10 + M_11 + M_12 + M_14
    //   row 2: M_21 + M_22 + M_23 + M_24
    //   row 3: M_30 + M_32 + M_33 + M_34
    //   row 4: M_40 + M_41 + M_42 + M_43 + M_44
    // -----------------------------------------------------------------------

    constexpr double c{ 62500.0 / 60.0 };
    constexpr double prefactor{ 2.0 };
    const double t{ tol_ * c };

    // node 0: e0={0,1,4} E=1, e3={3,0,4} E=4
    // M_00 = c*(1*(6*1+2*2+2*5) + 4*(6*1+2*4+2*5)) = c*(1*20+4*24) = c*116
    // M_01 = c*(1*(2*1+2*2+1*5)) = c*11
    // M_03 = c*(4*(2*3+2*1+1*5)) = c*60  [j=3,k=0 in e3]
    // M_04 = c*(1*(2*1+2*5+1*2) + 4*(2*1+2*5+1*4)) = c*(14+64) = c*78
    const double row0 = c*( 116.0 + 11.0 + 60.0 + 78.0 );

    // node 1: e0={0,1,4} E=1, e1={1,2,4} E=2
    // M_11 = c*(1*(6*2+2*1+2*5) + 2*(6*2+2*3+2*5)) = c*(24+56) = c*80
    // M_10 = M_01 = c*11
    // M_12 = c*(2*(2*2+2*3+1*5)) = c*30
    // M_14 = c*(1*(2*2+2*5+1*1) + 2*(2*2+2*5+1*3)) = c*(15+34) = c*49
    const double row1 = c*( 80.0 + 11.0 + 30.0 + 49.0 );

    // node 2: e1={1,2,4} E=2, e2={2,3,4} E=3
    // M_22 = c*(2*(6*3+2*2+2*5) + 3*(6*3+2*4+2*5)) = c*(64+108) = c*172
    // M_21 = M_12 = c*30
    // M_23 = c*(3*(2*3+2*4+1*5)) = c*57
    // M_24 = c*(2*(2*3+2*5+1*2) + 3*(2*3+2*5+1*4)) = c*(36+60) = c*96
    const double row2 = c*( 172.0 + 30.0 + 57.0 + 96.0 );

    // node 3: e2={2,3,4} E=3, e3={3,0,4} E=4
    // M_33 = c*(3*(6*4+2*3+2*5) + 4*(6*4+2*1+2*5)) = c*(120+144) = c*264
    // M_30 = M_03 = c*60
    // M_32 = M_23 = c*57
    // M_34 = c*(3*(2*4+2*5+1*3) + 4*(2*4+2*5+1*1)) = c*(63+76) = c*139
    const double row3 = c*( 264.0 + 60.0 + 57.0 + 139.0 );

    // node 4: all 4 elements
    // M_44 = c*(1*(6*5+2*1+2*2)+2*(6*5+2*2+2*3)+3*(6*5+2*3+2*4)+4*(6*5+2*4+2*1))
    //      = c*(36+80+132+160) = c*408
    // M_40 = M_04 = c*78
    // M_41 = M_14 = c*49
    // M_42 = M_24 = c*96
    // M_43 = M_34 = c*139
    const double row4 = c*( 408.0 + 78.0 + 49.0 + 96.0 + 139.0 );

    _equal( rhs_vec[0], prefactor * row0, t );
    _equal( rhs_vec[1], prefactor * row1, t );
    _equal( rhs_vec[2], prefactor * row2, t );
    _equal( rhs_vec[3], prefactor * row3, t );
    _equal( rhs_vec[4], prefactor * row4, t );
}



void Integral_var_NT_rhsop_N_dV_Test::compareConsistentTest() {
  compareTest(false);
}



void Integral_var_NT_rhsop_N_dV_Test::compareLumpedTest() {
  compareTest(true);
}


void Integral_var_NT_rhsop_N_dV_Test::compareTest(bool lumped)
{
  Integral_var_NT_lhsop_N_dV<2U> lhs(sg_->Database(),
                                     "conductivity",
                                     "fluid pressure",
                                     "fluid pressure",
                                     "mobility");
  lhs.LumpedFormulation(lumped);
  Integral_var_NT_rhsop_N_dV<2U> rhs(sg_->Database(),
                                     "conductivity",
                                     "fluid pressure",
                                     "fluid pressure",
                                     "mobility");
  rhs.LumpedFormulation(lumped);
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
    setElementVariable(conductivity, "conductivity");
    
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
    const csmp::Index    key(sg_->Database().StorageKey(var_name));
    ScalarVariable sc;
    Region<2U>&    domain{ sg_->Region("Model") };
        
    unsigned int i = 0; // BARF!
    for ( auto it = domain.NodesBegin(); it != domain.NodesEnd(); ++it) {
        (*it)->Read(key,sc);
        if ( verbose_ ) cout << "Node " << i << ": " << sc() << endl;
        ++i;
    }
}

void Integral_var_NT_rhsop_N_dV_Test::setNodeVariable( const vector<double>& var, const char* var_name) {
    //std::deque<Node<2U> >::iterator it;
    const csmp::Index key(sg_->Database().StorageKey(var_name));
    Region<2U>& domain{ sg_->Region("Model") };

    unsigned int i = 0;
    
    for ( auto it = domain.NodesBegin(); it != domain.NodesEnd(); ++it ) {
        (*it)->Store(key,ScalarVariable(PLAIN, static_cast<double>(var[i])));
        ++i;
    }
}

void Integral_var_NT_rhsop_N_dV_Test::setElementVariable( const vector<double>& var, const char* var_name) {
  //std::deque<Element<2U> >::iterator it;
  const csmp::Index key(sg_->Database().StorageKey(var_name));
  Region<2U>& domain{ sg_->Region("Model") };
  
  unsigned int i = 0;
  for (auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it) {
    (*it)->Store(key,ScalarVariable(PLAIN, static_cast<double>(var[i])));
              
    ++i;
  }
}

void Integral_var_NT_rhsop_N_dV_Test::calculateGlobalMatrix(SparseMatrix& sm, MathOperatorLHS<2U>& oper) {
    sm.Zero();
    sm.Resize(sg_->Region("Model").Nodes());
    Region<2U>& domain{ sg_->Region("Model") };
    
    //std::deque<Element<2U> >::iterator it;
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it) {
        oper.GetOperands( *(*it));
        oper.ComputeContribution(*(*it));
        oper.AssignToGlobal(*(*it), sm);
    }
    
}

void Integral_var_NT_rhsop_N_dV_Test::calculateGlobalRHS(vector<double>& rhs, MathOperatorRHS<2U>& oper) {
    Region<2U>& domain{ sg_->Region("Model") };
    //std::deque<Element<2U> >::iterator it;
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it) {
        oper.GetOperands(*(*it));
        oper.ComputeContribution(*(*it));
        oper.AssignToGlobal(*(*it), rhs);
    }
    
}

} // csmp
