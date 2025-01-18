#include "MathOperatorRHS_Test.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_SetRHS_to_One.h"

using namespace std;

namespace csmp {

MathOperatorRHS_Test::MathOperatorRHS_Test()
: model_ (new Model1D<1U>( "MathOperatorRHS_Test", "CSMP-2phase-variables.txt", 0.01, 100U )),
  database_ (model_->Database())

{
  fTolerance = 1.e-20;
}


MathOperatorRHS_Test::~MathOperatorRHS_Test()
{
   delete model_;
}


void MathOperatorRHS_Test::run()
{
    cout <<"\nMathOperatorRHS_Test::run: running test...\n";
    MathOperatorRHS_Ctor();
    MathOperatorRHS_CopyCtor();
    MathOperatorRHS_Equal();
    MathOperatorRHS_Name();
    MathOperatorRHS_OperandName();
    MathOperatorRHS_BasicOperandName();
    MathOperatorRHS_TestFunctionName();
    MathOperatorRHS_Subtract();
    MathOperatorRHS_Add();
    MathOperatorRHS_Multiply();
    MathOperatorRHS_AddLater();
    MathOperatorRHS_SubtractLater();
    MathOperatorRHS_MultiplyWithTimeIncrement();
    MathOperatorRHS_MultiplyBy();
    MathOperatorRHS_LumpedFormulation();
    MathOperatorRHS_ApplicationCycles();
    MathOperatorRHS_ApplicationCycle();
    MathOperatorRHS_Offset();
    MathOperatorRHS_MaterialOperand();
    MathOperatorRHS_BasicOperand();
    MathOperatorRHS_TestFunctionOperand();
}


void MathOperatorRHS_Test::MathOperatorRHS_Ctor()
{
    NumIntegral_SetRHS_to_One<1U> MORHS1( database_, "fluid pressure" );
    NumIntegral_NT_op_N_dV<1U> MORHS2( database_, "permeability", "fluid pressure" );
    csmp::Index pf_key = database_.StorageKey("fluid pressure");
    csmp::Index perm_key = database_.StorageKey("permeability");
    
    _test( MORHS1.Name() == "NumIntegral_SetRHS_to_One: Operand: 'fluid pressure" );
    _test( MORHS2.Name() == "NumIntegral_NT_op_N_dV: Operand: 'permeability', 'fluid pressure" );

    _test( MORHS2.MaterialOperandName() == "permeability" );

    _test( MORHS1.BasicOperandName() == "fluid pressure" );
    _test( MORHS2.BasicOperandName() == "fluid pressure" );

    _test( MORHS1.TestOperandName() == "fluid pressure" );
    _test( MORHS2.TestOperandName() == "fluid pressure" );

}



/** testing

        MTRL                 = mo.MTRL; // basic Operand storage
        RHS                  = mo.RHS;   // solution vector to be accumulated
        ID                   = mo.ID;     // node-ID & global constraint points vector
        DERIV                = mo.DERIV;
        IPOL                 = mo.IPOL;
        SC                   = mo.SC;
        VC                   = mo.VC;
        TS                   = mo.TS;
*/
void MathOperatorRHS_Test::MathOperatorRHS_CopyCtor()
{
   std::string op1("permeability");
    std::string top1("fluid pressure");
    std::string name1("MORHS1");
    bool add_accumulate = true;
    bool subtract_accumulate = false;
    bool multiply_accumulate = false;
    bool add_accumulate_later = false;
    bool subtract_accumulate_later = false;
    bool lump_matrices = true;
    bool time_multiply = true;
    uint32_t application_cycles = 5U;
    uint32_t application_cycle = 3U;
    double factor = 1.75;
    size_t offset = 2U;
    csmp::Index perm_key = database_.StorageKey(op1.c_str());
    csmp::Index pf_key = database_.StorageKey(top1.c_str());
    
    
    NumIntegral_NT_op_N_dV<1U> MORHS1( database_, op1.c_str(), top1.c_str() );
    
    MORHS1.Name("MORHS1", top1.c_str());
    MORHS1.AddAccumulate();
    MORHS1.LumpedFormulation(lump_matrices);
    MORHS1.MultiplyWithTimeIncrement(time_multiply);
    MORHS1.ApplicationCycles(application_cycles);
    MORHS1.ApplicationCycle(application_cycle);
    MORHS1.ApplicationCycles(application_cycles);
    MORHS1.MultiplyBy(factor);
    MORHS1.TestOperandOffset(offset);
    
    NumIntegral_NT_op_N_dV<1U> MORHS2(MORHS1);
    
    _test( MORHS2.Name() == MORHS1.Name() );
    _test( MORHS2.MaterialOperandName() == op1 );
    _test( MORHS2.TestOperandName() == top1 );
    _test( MORHS2.Add() == add_accumulate );
    _test( MORHS2.Multiply() == multiply_accumulate );
    _test( MORHS2.Subtract() == subtract_accumulate );
    _test( MORHS2.AddLater() == add_accumulate_later );
    _test( MORHS2.SubtractLater() == subtract_accumulate_later );
    _test( MORHS2.LumpedFormulation() == lump_matrices );
    _test( MORHS2.MultiplyWithTimeIncrement() == time_multiply );
    _test( MORHS2.ApplicationCycles() == application_cycles );
    _test( MORHS2.ApplicationCycle() == application_cycle );
    _equal( MORHS2.MultiplyBy(), factor, fTolerance );
    _test( MORHS2.TestOperandOffset() == offset );
    _test( MORHS2.MaterialOperandKey() == perm_key );
    _test( MORHS2.BasicOperandKey() == pf_key );
    _test( MORHS2.TestOperandKey() == pf_key );

}


void MathOperatorRHS_Test::MathOperatorRHS_Equal()
{
    string op1("permeability");
    string top1("fluid pressure");
    string op2("diffusivity");
    string top2("hydrostatic pressure");
    std::string name1("MORHS1");
    bool add_accumulate = true;
    bool subtract_accumulate = false;
    bool multiply_accumulate = false;
    bool add_accumulate_later = false;
    bool subtract_accumulate_later = false;
    bool lump_matrices = true;
    bool time_multiply = true;
    uint32_t application_cycles = 5U;
    uint32_t application_cycle = 3U;
    double factor = 1.75;
    size_t offset = 2U;
    csmp::Index perm_key = database_.StorageKey(op1.c_str());
    csmp::Index pf_key = database_.StorageKey(top1.c_str());
    
    
    NumIntegral_NT_op_N_dV<1U> MORHS1( database_, op1.c_str(), top1.c_str() );
    NumIntegral_NT_op_N_dV<1U> MORHS2( database_, op2.c_str(), top2.c_str() );
    
    MORHS1.Name("MORHS1", top1.c_str());
    MORHS1.AddAccumulate();
    MORHS1.LumpedFormulation(lump_matrices);
    MORHS1.MultiplyWithTimeIncrement(time_multiply);
    MORHS1.ApplicationCycles(application_cycles);
    MORHS1.ApplicationCycle(application_cycle);
    MORHS1.ApplicationCycles(application_cycles);
    MORHS1.MultiplyBy(factor);
    MORHS1.TestOperandOffset(offset);

    MORHS2 = MORHS1;
    
    _test( MORHS2.Name() == MORHS1.Name() );
    _test( MORHS2.MaterialOperandName() == op1 );
    _test( MORHS2.TestOperandName() == top1 );
    _test( MORHS2.Add() == add_accumulate );
    _test( MORHS2.Multiply() == multiply_accumulate );
    _test( MORHS2.Subtract() == subtract_accumulate );
    _test( MORHS2.AddLater() == add_accumulate_later );
    _test( MORHS2.SubtractLater() == subtract_accumulate_later );
    _test( MORHS2.LumpedFormulation() == lump_matrices );
    _test( MORHS2.MultiplyWithTimeIncrement() == time_multiply );
    _test( MORHS2.ApplicationCycles() == application_cycles );
    _test( MORHS2.ApplicationCycle() == application_cycle );
    _equal( MORHS2.MultiplyBy(), factor, fTolerance );
    _test( MORHS2.TestOperandOffset() == offset );
    _test( MORHS2.MaterialOperandKey() == perm_key );
    _test( MORHS2.BasicOperandKey() == pf_key );
    _test( MORHS2.TestOperandKey() == pf_key );

}




void MathOperatorRHS_Test::MathOperatorRHS_Name()
{
    string op("permeability");
    string top("fluid pressure");

    NumIntegral_NT_op_N_dV<1U> MORHS1( database_, op.c_str(), top.c_str() );
    NumIntegral_NT_op_N_dV<1U> MORHS2( database_, op.c_str(), top.c_str() );
    
    _test( MORHS1.Name() == "NumIntegral_NT_op_N_dV: Operand: 'permeability', 'fluid pressure" );
    _test( MORHS2.Name() == "NumIntegral_NT_op_N_dV: Operand: 'permeability', 'fluid pressure" );
    
    MORHS1.Name("MORHS1", top.c_str());
    MORHS2.Name("MORHS2", op.c_str(), top.c_str());
    
    std::string name1;
    std::string name2;
    
    name1 = "MORHS1";
    name1 += ": Operand: '";
    name1 += top;
    name2 = "MORHS2";
    name2 += ": Operand: '";
    name2 += op;
    name2 += "', '";
    name2 += top;

    _test( MORHS1.Name() == name1 );
    _test( MORHS2.Name() == name2 );


}



void MathOperatorRHS_Test::MathOperatorRHS_OperandName()
{
    NumIntegral_NT_op_N_dV<1U> MORHS2( database_, "permeability", "fluid pressure" );

    _test( MORHS2.MaterialOperandName() == "permeability" );

}


void MathOperatorRHS_Test::MathOperatorRHS_BasicOperandName()
{
    NumIntegral_SetRHS_to_One<1U> MORHS1( database_, "fluid pressure" );
    NumIntegral_NT_op_N_dV<1U> MORHS2( database_, "permeability", "fluid pressure" );

    _test( MORHS1.BasicOperandName() == "fluid pressure" );
    _test( MORHS2.BasicOperandName() == "fluid pressure" );

}


void MathOperatorRHS_Test::MathOperatorRHS_TestFunctionName()
{
    NumIntegral_SetRHS_to_One<1U> MORHS1( database_, "fluid pressure" );
    NumIntegral_NT_op_N_dV<1U> MORHS2( database_, "permeability", "fluid pressure" );

    _test( MORHS1.TestOperandName() == "fluid pressure" );
    _test( MORHS2.TestOperandName() == "fluid pressure" );

}


void MathOperatorRHS_Test::MathOperatorRHS_Subtract()
{
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, "permeability", "fluid pressure" );
    
    MORHS.SubtractAccumulate();
    
    _test( MORHS.Add() == false );
    _test( MORHS.Subtract() == true );
    _test( MORHS.Multiply() == false );
    _test( MORHS.AddLater() == false );
    _test( MORHS.SubtractLater() == false );
    
}


void MathOperatorRHS_Test::MathOperatorRHS_Add()
{
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, "permeability", "fluid pressure" );
    
    MORHS.AddAccumulate();
    
    _test( MORHS.Add() == true );
    _test( MORHS.Subtract() == false );
    _test( MORHS.Multiply() == false );
    _test( MORHS.AddLater() == false );
    _test( MORHS.SubtractLater() == false );

}


void MathOperatorRHS_Test::MathOperatorRHS_Multiply()
{
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, "permeability", "fluid pressure" );
    
    MORHS.MultiplyAccumulate();
    
    _test( MORHS.Add() == false );
    _test( MORHS.Subtract() == false );
    _test( MORHS.Multiply() == true );
    _test( MORHS.AddLater() == false );
    _test( MORHS.SubtractLater() == false );

}


void MathOperatorRHS_Test::MathOperatorRHS_AddLater()
{
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, "permeability", "fluid pressure" );
    
    MORHS.AddAccumulateLater();
    
    _test( MORHS.Add() == false );
    _test( MORHS.Subtract() == false );
    _test( MORHS.Multiply() == false );
    _test( MORHS.AddLater() == true );
    _test( MORHS.SubtractLater() == false );

}


void MathOperatorRHS_Test::MathOperatorRHS_SubtractLater()
{
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, "permeability", "fluid pressure" );
    
    MORHS.SubtractAccumulateLater();
    
    _test( MORHS.Add() == false );
    _test( MORHS.Subtract() == false );
    _test( MORHS.Multiply() == false );
    _test( MORHS.AddLater() == false );
    _test( MORHS.SubtractLater() == true );

}


void MathOperatorRHS_Test::MathOperatorRHS_MultiplyWithTimeIncrement()
{
    string op("permeability");
    string top("fluid pressure");
    bool time_multiply = true;
    
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, op.c_str(), top.c_str() );
    
    _test( MORHS.MultiplyWithTimeIncrement() == false );
    
    MORHS.MultiplyWithTimeIncrement(time_multiply);
    
    _test( MORHS.MultiplyWithTimeIncrement() == time_multiply );

}


void MathOperatorRHS_Test::MathOperatorRHS_MultiplyBy()
{
    string op("permeability");
    string top("fluid pressure");
    double factor = 1.75;    
    
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, op.c_str(), top.c_str() );
    
    _equal( MORHS.MultiplyBy(), 1., fTolerance );
    
    MORHS.MultiplyBy(factor);
    
    _equal( MORHS.MultiplyBy(), factor, fTolerance );

}


void MathOperatorRHS_Test::MathOperatorRHS_LumpedFormulation()
{
    string op("permeability");
    string top("fluid pressure");
    bool lump_matrices = true;    
    
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, op.c_str(), top.c_str() );
    
    _test( MORHS.LumpedFormulation() == false );
    
    MORHS.LumpedFormulation(lump_matrices);
    
    _test( MORHS.LumpedFormulation() == lump_matrices );

}


void MathOperatorRHS_Test::MathOperatorRHS_ApplicationCycles()
{
    string op("permeability");
    string top("fluid pressure");
    uint32_t application_cycles = 5U;
    
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, op.c_str(), top.c_str() );
    
    _test( MORHS.ApplicationCycles() == 1U );
    
    MORHS.ApplicationCycles(application_cycles);

    _test( MORHS.ApplicationCycles() == application_cycles );

}


void MathOperatorRHS_Test::MathOperatorRHS_ApplicationCycle()
{
    string op("permeability");
    string top("fluid pressure");
    uint32_t application_cycle = 3U;
    
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, op.c_str(), top.c_str() );
    
    _test( MORHS.ApplicationCycle() == 0U );
    
    MORHS.ApplicationCycle(application_cycle);
    
    _test( MORHS.ApplicationCycle() == application_cycle );

}


void MathOperatorRHS_Test::MathOperatorRHS_Offset()
{
    string op("permeability");
    string top("fluid pressure");
    size_t offset = 2U;
    
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, op.c_str(), top.c_str() );
    
    _test( MORHS.TestOperandOffset() == 0U );
    
    MORHS.TestOperandOffset(offset);

    _test( MORHS.TestOperandOffset() == offset );

}


void MathOperatorRHS_Test::MathOperatorRHS_MaterialOperand()
{
    string op("permeability");
    string top("fluid pressure");
    
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, op.c_str(), top.c_str() );
    
    csmp::Index perm_key = database_.StorageKey(op.c_str());
    
    _test( MORHS.MaterialOperandKey() == perm_key );
    
}


void MathOperatorRHS_Test::MathOperatorRHS_BasicOperand()
{
    string op("permeability");
    string top("fluid pressure");
    
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, op.c_str(), top.c_str() );
    
    csmp::Index pf_key = database_.StorageKey(top.c_str());
    
    _test( MORHS.BasicOperandKey() == pf_key );

}


void MathOperatorRHS_Test::MathOperatorRHS_TestFunctionOperand()
{
    string op("permeability");
    string top("fluid pressure");
    
    NumIntegral_NT_op_N_dV<1U> MORHS( database_, op.c_str(), top.c_str() );
    
    csmp::Index pf_key = database_.StorageKey(top.c_str());
    
    _test( MORHS.TestOperandKey() == pf_key );
    
}

} // end namespace csmp
