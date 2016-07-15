#include "MathOperatorLHS_Test.h"
//#include "MathOperatorLHS.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_dNT_dN_dV.h"


namespace csmp {



MathOperatorLHS_Test::MathOperatorLHS_Test()
: model_ (new Model1D<1U>( "MathOperatorLHS_Test", "CSMP-2phase-variables.txt", 0.01, 100U )),
  database_ (model_->Database())

{
  fTolerance = 1.e-20;
}


MathOperatorLHS_Test::~MathOperatorLHS_Test()
{
    if(model_!=NULL)
        delete model_;
}


void MathOperatorLHS_Test::run()
{

    MathOperatorLHS_Ctor();
    MathOperatorLHS_CopyCtor();
    MathOperatorLHS_Equal();
    MathOperatorLHS_Name();
    MathOperatorLHS_OperandName();
    MathOperatorLHS_BasicOperandName();
    MathOperatorLHS_TestFunctionName();
    MathOperatorLHS_Add();
    MathOperatorLHS_Multiply();
    MathOperatorLHS_AddLater();
    MathOperatorLHS_MultiplyWithTimeIncrement();
    MathOperatorLHS_MultiplyBy();
    MathOperatorLHS_LumpedFormulation();
    MathOperatorLHS_BasicOffset();
    MathOperatorLHS_TestOffset();
    MathOperatorLHS_ApplicationCycles();
    MathOperatorLHS_ApplicationCycle();
    MathOperatorLHS_MaterialOperand();
    MathOperatorLHS_BasicOperand();
    MathOperatorLHS_TestFunctionOperand();


}


void MathOperatorLHS_Test::MathOperatorLHS_Ctor()
{
    NumIntegral_dNT_dN_dV<1U> MOLHS1( database_, "fluid pressure", "hydrostatic pressure" );
    NumIntegral_dNT_op_dN_dV<1U> MOLHS2( database_, "permeability", "fluid pressure", "hydrostatic pressure" );
    csmp::Index ph_key = database_.StorageKey("hydrostatic pressure");
    csmp::Index pf_key = database_.StorageKey("fluid pressure");
    csmp::Index perm_key = database_.StorageKey("permeability");

    _test( MOLHS1.Name() == "NumIntegral_dNT_dN_dV: Operand: 'fluid pressure', 'hydrostatic pressure" );
    _test( MOLHS2.Name() == "NumIntegral_dNT_op_dN_dV: Operand: 'permeability', 'fluid pressure', 'hydrostatic pressure" );

    _test( MOLHS2.MaterialOperandName() == "permeability" );

    _test( MOLHS1.BasicOperandName() == "fluid pressure" );
    _test( MOLHS2.BasicOperandName() == "fluid pressure" );

    _test( MOLHS1.TestOperandName() == "hydrostatic pressure" );
    _test( MOLHS2.TestOperandName() == "hydrostatic pressure" );
}


void MathOperatorLHS_Test::MathOperatorLHS_CopyCtor()
{
    
    string op1("permeability");
    string bop1("fluid pressure");
    string top1("hydrostatic pressure");
    std::string name1("MOLHS1");
    bool add_accumulate = true;
    bool multiply_accumulate = false;
    bool add_accumulate_later = false;
    bool lump_matrices = true;
    bool time_multiply = true;
    size_t application_cycles = 5U;
    size_t application_cycle = 3U;
    double64 factor = 1.75;
    size_t basic_offset = 2U;
    size_t test_offset = 4U;
    csmp::Index perm_key = database_.StorageKey(op1.c_str());
    csmp::Index pf_key = database_.StorageKey(bop1.c_str());
    csmp::Index ph_key = database_.StorageKey(top1.c_str());
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS1( database_, op1.c_str(), bop1.c_str(), top1.c_str() );
    
    MOLHS1.Name("MOLHS1", bop1.c_str(), top1.c_str());
    MOLHS1.AddAccumulate();
    MOLHS1.LumpedFormulation(lump_matrices);
    MOLHS1.MultiplyWithTimeIncrement(time_multiply);
    MOLHS1.ApplicationCycles(application_cycles);
    MOLHS1.ApplicationCycle(application_cycle);
    MOLHS1.ApplicationCycles(application_cycles);
    MOLHS1.MultiplyBy(factor);
    MOLHS1.BasicOperandOffset(basic_offset);
    MOLHS1.TestOperandOffset(test_offset);

    NumIntegral_dNT_op_dN_dV<1U> MOLHS2(MOLHS1);
    
    _test( MOLHS2.Name() == MOLHS1.Name() );
    _test( MOLHS2.MaterialOperandName() == op1 );
    _test( MOLHS2.BasicOperandName() == bop1 );
    _test( MOLHS2.TestOperandName() == top1 );
    _test( MOLHS2.Add() == add_accumulate );
    _test( MOLHS2.Multiply() == multiply_accumulate );
    _test( MOLHS2.AddLater() == add_accumulate_later );
    _test( MOLHS2.LumpedFormulation() == lump_matrices );
    _test( MOLHS2.MultiplyWithTimeIncrement() == time_multiply );
    _test( MOLHS2.ApplicationCycles() == application_cycles );
    _test( MOLHS2.ApplicationCycle() == application_cycle );
    _equal( MOLHS2.MultiplyBy(), factor, fTolerance );
    _test( MOLHS2.BasicOperandOffset() == basic_offset );
    _test( MOLHS2.TestOperandOffset() == test_offset );
    _test( MOLHS2.MaterialOperandKey() == perm_key );
    _test( MOLHS2.BasicOperandKey() == pf_key );
    _test( MOLHS2.TestOperandKey() == ph_key );

}


void MathOperatorLHS_Test::MathOperatorLHS_Equal()
{
    string op1("permeability");
    string bop1("fluid pressure");
    string top1("hydrostatic pressure");
    string op2("diffusivity");
    string bop2("hydrostatic pressure");
    string top2("fluid pressure");
    std::string name1("MOLHS1");
    bool add_accumulate = true;
    bool multiply_accumulate = false;
    bool add_accumulate_later = false;
    bool lump_matrices = true;
    bool time_multiply = true;
    size_t application_cycles = 5U;
    size_t application_cycle = 3U;
    double64 factor = 1.75;
    size_t basic_offset = 2U;
    size_t test_offset = 4U;
    csmp::Index perm_key = database_.StorageKey(op1.c_str());
    csmp::Index pf_key = database_.StorageKey(bop1.c_str());
    csmp::Index ph_key = database_.StorageKey(top1.c_str());
    
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS1( database_, op1.c_str(), bop1.c_str(), top1.c_str() );
    NumIntegral_dNT_op_dN_dV<1U> MOLHS2( database_, op2.c_str(), bop2.c_str(), top2.c_str() );
    
    MOLHS1.Name("MOLHS1", bop1.c_str(), top1.c_str());
    MOLHS1.AddAccumulate();
    MOLHS1.LumpedFormulation(lump_matrices);
    MOLHS1.MultiplyWithTimeIncrement(time_multiply);
    MOLHS1.ApplicationCycles(application_cycles);
    MOLHS1.ApplicationCycle(application_cycle);
    MOLHS1.ApplicationCycles(application_cycles);
    MOLHS1.MultiplyBy(factor);
    MOLHS1.BasicOperandOffset(basic_offset);
    MOLHS1.TestOperandOffset(test_offset);

    MOLHS2 = MOLHS1;
    
    _test( MOLHS2.Name() == MOLHS1.Name() );
    _test( MOLHS2.MaterialOperandName() == op1 );
    _test( MOLHS2.BasicOperandName() == bop1 );
    _test( MOLHS2.TestOperandName() == top1 );
    _test( MOLHS2.Add() == add_accumulate );
    _test( MOLHS2.Multiply() == multiply_accumulate );
    _test( MOLHS2.AddLater() == add_accumulate_later );
    _test( MOLHS2.LumpedFormulation() == lump_matrices );
    _test( MOLHS2.MultiplyWithTimeIncrement() == time_multiply );
    _test( MOLHS2.ApplicationCycles() == application_cycles );
    _test( MOLHS2.ApplicationCycle() == application_cycle );
    _equal( MOLHS2.MultiplyBy(), factor, fTolerance );
    _test( MOLHS2.BasicOperandOffset() == basic_offset );
    _test( MOLHS2.TestOperandOffset() == test_offset );
    _test( MOLHS2.MaterialOperandKey() == perm_key );
    _test( MOLHS2.BasicOperandKey() == pf_key );
    _test( MOLHS2.TestOperandKey() == ph_key );

}


void MathOperatorLHS_Test::MathOperatorLHS_Name()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");

    NumIntegral_dNT_op_dN_dV<1U> MOLHS1( database_, op.c_str(), bop.c_str(), top.c_str() );
    NumIntegral_dNT_op_dN_dV<1U> MOLHS2( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    _test( MOLHS1.Name() == "NumIntegral_dNT_op_dN_dV: Operand: 'permeability', 'fluid pressure', 'hydrostatic pressure" );
    _test( MOLHS2.Name() == "NumIntegral_dNT_op_dN_dV: Operand: 'permeability', 'fluid pressure', 'hydrostatic pressure" );
    
    MOLHS1.Name("MOLHS1", bop.c_str(), top.c_str());
    MOLHS2.Name("MOLHS2", op.c_str(), bop.c_str(), top.c_str());
    
    std::string name1;
    std::string name2;
    
    name1 = "MOLHS1";
    name1 += ": Operand: '";
    name1 += bop;
    name1 += "', '";
    name1 += top;
    name2 = "MOLHS2";
    name2 += ": Operand: '";
    name2 += op;
    name2 += "', '";
    name2 += bop;
    name2 += "', '";
    name2 += top;

    _test( MOLHS1.Name() == name1 );
    _test( MOLHS2.Name() == name2 );


}


void MathOperatorLHS_Test::MathOperatorLHS_OperandName()
{
    NumIntegral_dNT_op_dN_dV<1U> MOLHS2( database_, "permeability", "fluid pressure", "hydrostatic pressure" );

    _test( MOLHS2.MaterialOperandName() == "permeability" );

}


void MathOperatorLHS_Test::MathOperatorLHS_BasicOperandName()
{
    NumIntegral_dNT_dN_dV<1U> MOLHS1( database_, "fluid pressure", "hydrostatic pressure" );
    NumIntegral_dNT_op_dN_dV<1U> MOLHS2( database_, "permeability", "fluid pressure", "hydrostatic pressure" );

    _test( MOLHS1.BasicOperandName() == "fluid pressure" );
    _test( MOLHS2.BasicOperandName() == "fluid pressure" );

}


void MathOperatorLHS_Test::MathOperatorLHS_TestFunctionName()
{
    NumIntegral_dNT_dN_dV<1U> MOLHS1( database_, "fluid pressure", "hydrostatic pressure" );
    NumIntegral_dNT_op_dN_dV<1U> MOLHS2( database_, "permeability", "fluid pressure", "hydrostatic pressure" );

    _test( MOLHS1.TestOperandName() == "hydrostatic pressure" );
    _test( MOLHS2.TestOperandName() == "hydrostatic pressure" );

}


void MathOperatorLHS_Test::MathOperatorLHS_Add()
{
    NumIntegral_dNT_dN_dV<1U> MOLHS( database_, "fluid pressure", "hydrostatic pressure" );
    
    MOLHS.AddAccumulate();
    
    _test( MOLHS.Add() == true );
    _test( MOLHS.Multiply() == false );
    _test( MOLHS.AddLater() == false );

}


void MathOperatorLHS_Test::MathOperatorLHS_Multiply()
{
    NumIntegral_dNT_dN_dV<1U> MOLHS( database_, "fluid pressure", "hydrostatic pressure" );
    
    MOLHS.MultiplyAccumulate();
    
    _test( MOLHS.Add() == false );
    _test( MOLHS.Multiply() == true );
    _test( MOLHS.AddLater() == false );

}


void MathOperatorLHS_Test::MathOperatorLHS_AddLater()
{
    NumIntegral_dNT_dN_dV<1U> MOLHS( database_, "fluid pressure", "hydrostatic pressure" );
    
    MOLHS.AddAccumulateLater();
    
    _test( MOLHS.Add() == false );
    _test( MOLHS.Multiply() == false );
    _test( MOLHS.AddLater() == true );

}


void MathOperatorLHS_Test::MathOperatorLHS_MultiplyWithTimeIncrement()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    bool time_multiply = true;
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    _test( MOLHS.MultiplyWithTimeIncrement() == false );
    
    MOLHS.MultiplyWithTimeIncrement(time_multiply);
    
    _test( MOLHS.MultiplyWithTimeIncrement() == time_multiply );

}


void MathOperatorLHS_Test::MathOperatorLHS_MultiplyBy()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    double64 factor = 1.75;    
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    _equal( MOLHS.MultiplyBy(), 1., fTolerance );
    
    MOLHS.MultiplyBy(factor);
    
    _equal( MOLHS.MultiplyBy(), factor, fTolerance );

}


void MathOperatorLHS_Test::MathOperatorLHS_LumpedFormulation()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    bool lump_matrices = true;    
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    _test( MOLHS.LumpedFormulation() == false );
    
    MOLHS.LumpedFormulation(lump_matrices);
    
    _test( MOLHS.LumpedFormulation() == lump_matrices );

}



void MathOperatorLHS_Test::MathOperatorLHS_BasicOffset()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    size_t basic_offset = 2U;
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    _test( MOLHS.BasicOperandOffset() == 0U );
    
    MOLHS.BasicOperandOffset(basic_offset);

    _test( MOLHS.BasicOperandOffset() == basic_offset );

}


void MathOperatorLHS_Test::MathOperatorLHS_TestOffset()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    size_t test_offset = 3U;
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    _test( MOLHS.TestOperandOffset() == 0U );
    
    MOLHS.TestOperandOffset(test_offset);

    _test( MOLHS.TestOperandOffset() == test_offset );

}



void MathOperatorLHS_Test::MathOperatorLHS_ApplicationCycles()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    size_t application_cycles = 5U;
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );

    _test( MOLHS.ApplicationCycles() == 1U );
    
    MOLHS.ApplicationCycles(application_cycles);

    _test( MOLHS.ApplicationCycles() == application_cycles );

}


void MathOperatorLHS_Test::MathOperatorLHS_ApplicationCycle()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    size_t application_cycle = 3U;
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    _test( MOLHS.ApplicationCycle() == 0U );
    
    MOLHS.ApplicationCycle(application_cycle);
    
    _test( MOLHS.ApplicationCycle() == application_cycle );

}


void MathOperatorLHS_Test::MathOperatorLHS_MaterialOperand()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    csmp::Index perm_key = database_.StorageKey(op.c_str());
    
    _test( MOLHS.MaterialOperandKey() == perm_key );
    
}


void MathOperatorLHS_Test::MathOperatorLHS_BasicOperand()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    csmp::Index pf_key = database_.StorageKey(bop.c_str());
    
    _test( MOLHS.BasicOperandKey() == pf_key );

}


void MathOperatorLHS_Test::MathOperatorLHS_TestFunctionOperand()
{
    string op("permeability");
    string bop("fluid pressure");
    string top("hydrostatic pressure");
    
    NumIntegral_dNT_op_dN_dV<1U> MOLHS( database_, op.c_str(), bop.c_str(), top.c_str() );
    
    csmp::Index ph_key = database_.StorageKey(top.c_str());
    
    _test( MOLHS.TestOperandKey() == ph_key );
    
}


} // end namespace csmp
