#ifndef SplitBoundaryTest_JK_H
#define SplitBoundaryTest_JK_H

#include "Test.h"
#include "SplitBoundary.h"
#include "Model.h"

using namespace std;
using namespace csmp;


namespace csmp
{

/** SplitBoundaryTest_JK : 2D Test Case
  =================================
  Mesh:       rectangle_2SBs
  Model:      2D
  Test:       Split boundaries
  BC:         TBC
  Criterion:
  =================================
  */

class SplitBoundaryTest_JK : public Test
{
public:
    SplitBoundaryTest_JK ();
    ~SplitBoundaryTest_JK();
    
    void run() override final;
    
private:

    void Test_CreateSplitBoundaryFrom();
    void Test_CreateSplitBoundaryFrom_Simplified();
    
    const char* geometry_name_;
    const char* config_file_name_;
    const char* vars_name_;
    const char* system_file_list_name_;
    const char* recipes_file_list_name_;
    list<string> output_props_;

    enum{dim=2};

    Model<dim>*  model_ = nullptr;

    // //! CVFEM instance
    // Coupled_reservoir_well_CVFEM_PHX_Scheme<dim>* CVFEM_PHX;

    // //! Permeability Visitor
    // PermeabilityVisitor<dim>* permeability_visitor;

    string
        output_prefix_,
        output_name_,
        restart_file_,
        model_ID_path_,
        output_path;
};

} // csmp

#endif // SplitBoundaryTest_JK
