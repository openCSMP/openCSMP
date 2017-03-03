#include "MathOperatorRHS.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim>
MathOperatorRHS<dim>::MathOperatorRHS()
 :  name_("unspecified RHS-operator"),
    IDT(3),
    MTRL(13,DenseMatrix<DM_MIN>(2,2)),
    factor_(1.),
    add_accumulate_(true),
    subtract_accumulate_(false),
    add_accumulate_later_(false),
    subtract_accumulate_later_(false),
    multiply_accumulate_(false),
    lump_matrices_(false),
    application_cycles_(1),
    application_cycle_(0),
    time_multiply_(false),
    time_divide_(false)
 {
 }


template<size_t dim>
MathOperatorRHS<dim>::MathOperatorRHS( const PropertyDatabase<dim>& pref,
                                       const char* test )
  : name_("unspecified RHS-operator"),
    top(make_pair(pref.Parameter(test),0)),
    IDT(3),
    MTRL(13,DenseMatrix<DM_MIN>(2,2)),
    DERIV(2,3),
    IPOL(3),
    SC(3),
    VC(3),
    TS(3),
    AR(3),
    FR(3),
    factor_(1.),
    add_accumulate_(true),
    subtract_accumulate_(false),
    add_accumulate_later_(false),
    subtract_accumulate_later_(false),
    multiply_accumulate_(false),
    lump_matrices_(false),
    application_cycles_(1),
    application_cycle_(0),
    time_multiply_(false),
    time_divide_(false)
 {
 }


template<size_t dim>
MathOperatorRHS<dim>::MathOperatorRHS( const PropertyDatabase<dim>& pref, 
                                       const char* oper,
                                       const char* test )
  : name_("unspecified RHS-operator"),
    op(pref.Parameter(oper)),
    top(make_pair(pref.Parameter(test),0)),
    IDT(3),
    MTRL(13,DenseMatrix<DM_MIN>(2,2)),
    DERIV(2,3),
    IPOL(3),
    SC(3),
    VC(3),
    TS(3),
    AR(3),
    FR(3),
    factor_(1.),
    add_accumulate_(true),
    subtract_accumulate_(false),
    add_accumulate_later_(false),
    subtract_accumulate_later_(false),
    multiply_accumulate_(false),
    lump_matrices_(false),
    application_cycles_(1),
    application_cycle_(0),
    time_multiply_(false),
    time_divide_(false)

{
}



template<size_t dim>
MathOperatorRHS<dim>::MathOperatorRHS( const MathOperatorRHS<dim>& mo )
 :  name_                       (mo.name_),
    op                          (mo.op),    // basic operand
    top                         (mo.top),   // testfunction operand
    RHS                         (mo.RHS),   // solution vector to be accumulated
    IDT                         (mo.IDT),   // node-ID & global constraint points vector
    MTRL                        (mo.MTRL),  // basic Operand storage
    DERIV                       (mo.DERIV),
    IPOL                        (mo.IPOL),
    SC                          (mo.SC),
    VC                          (mo.VC),
    TS                          (mo.TS),
    AR                          (mo.AR),
    FR                          (mo.FR),
    factor_                     (mo.factor_),
    add_accumulate_             (mo.add_accumulate_),
    subtract_accumulate_        (mo.subtract_accumulate_),
    add_accumulate_later_       (mo.add_accumulate_later_),
    subtract_accumulate_later_  (mo.subtract_accumulate_later_),
    multiply_accumulate_        (mo.multiply_accumulate_),
    lump_matrices_              (mo.lump_matrices_),
    application_cycles_         (mo.application_cycles_),
    application_cycle_          (mo.application_cycle_),
    time_multiply_              (mo.time_multiply_),
    time_divide_                (mo.time_divide_)
{
}



template<size_t dim>
MathOperatorRHS<dim>& MathOperatorRHS<dim>::operator=( const MathOperatorRHS<dim>& mo )
 {
    if ( this != &mo ) {
        name_                       = mo.name_;
        op                          = mo.op;        // basic operand
        top                         = mo.top;       // testfunction operand
        RHS                         = mo.RHS;       // solution vector to be accumulated
        IDT                         = mo.IDT;       // node-ID & global constraint points vector
        MTRL                        = mo.MTRL;      // basic Operand storage
        DERIV                       = mo.DERIV;
        IPOL                        = mo.IPOL;
        SC                          = mo.SC;
        VC                          = mo.VC;
        TS                          = mo.TS;
        AR                          = mo.AR;
        FR                          = mo.FR;
        add_accumulate_             = mo.add_accumulate_;
        subtract_accumulate_        = mo.subtract_accumulate_;
        add_accumulate_later_       = mo.add_accumulate_later_;
        subtract_accumulate_later_  = mo.subtract_accumulate_later_;
        multiply_accumulate_        = mo.multiply_accumulate_;
        lump_matrices_              = mo.lump_matrices_;
        application_cycles_         = mo.application_cycles_;
        application_cycle_          = mo.application_cycle_;
        time_multiply_              = mo.time_multiply_;
        time_divide_                = mo.time_divide_;
        factor_                     = mo.factor_;
      }
    return *this;
 }

template<size_t dim>
MathOperatorRHS<dim>::~MathOperatorRHS()
{
}


/// ===================================  Operand Functions ========================================================================

template<size_t dim>
void MathOperatorRHS<dim>::Name( const char* s, const char* topname )
 {
    name_  = s;
    name_ += ": Operand: '";
    name_ += topname;
 }



template<size_t dim>
void MathOperatorRHS<dim>::Name( const char* s, const char* operand_name, const char* topname )
 {
    name_  = s;
    name_ += ": Operand: '";
    name_ += operand_name;
    name_ += "', '";
    name_ += topname;
 }

template<size_t dim>
std::string  MathOperatorRHS<dim>::Name() const
{
    return name_;
}


template<size_t dim>
void MathOperatorRHS<dim>::Out(std::ostream& os) const
{
    op.Out(os);
    top.first.Out(os);
}




// MATERIAL PROPERTY OPERAND

template<size_t dim>
const csmp::Parameter&   MathOperatorRHS<dim>::MaterialOperand() const
{
    return op;
}

template<size_t dim>
const csmp::Index&   MathOperatorRHS<dim>::MaterialOperandKey() const
{
    return op.key;
}

template<size_t dim>
std::string  MathOperatorRHS<dim>::MaterialOperandName() const
{
    return op.name;
}

template<size_t dim>
csmp::VARIABLE_TYPE   MathOperatorRHS<dim>::MaterialOperandType() const
{
    return op.key.type;
}

template<size_t dim>
csmp::PLACEMENT   MathOperatorRHS<dim>::MaterialOperandPlacement() const
{
    return op.key.place;
}

template<size_t dim>
size_t   MathOperatorRHS<dim>::MaterialOperandDataDepth() const
{
    return op.key.dataDepth;
}


// BASIC FUNCTION OPERAND

template<size_t dim>
const csmp::Parameter&  MathOperatorRHS<dim>::BasicOperand() const
{
    return top.first;
}

template<size_t dim>
const csmp::Index&   MathOperatorRHS<dim>::BasicOperandKey() const
{
    return top.first.key;
}

template<size_t dim>
std::string  MathOperatorRHS<dim>::BasicOperandName() const
{
    return top.first.name;
}

template<size_t dim>
csmp::VARIABLE_TYPE   MathOperatorRHS<dim>::BasicOperandType() const
{
    return top.first.key.type;
}

template<size_t dim>
csmp::PLACEMENT   MathOperatorRHS<dim>::BasicOperandPlacement() const
{
    return top.first.key.place;
}

template<size_t dim>
size_t   MathOperatorRHS<dim>::BasicOperandDataDepth() const
{
    return top.first.key.dataDepth;
}

template<size_t dim>
void MathOperatorRHS<dim>::BasicOperandOffset( size_t os )
{
    top.second = os;
}

template<size_t dim>
size_t MathOperatorRHS<dim>::BasicOperandOffset() const
{
    return top.second;
}

// TEST FUNCTION OPERAND

template<size_t dim>
const csmp::Parameter&   MathOperatorRHS<dim>::TestOperand() const
{
    return top.first;
}

template<size_t dim>
const csmp::Index&   MathOperatorRHS<dim>::TestOperandKey() const
{
    return top.first.key;
}

template<size_t dim>
std::string  MathOperatorRHS<dim>::TestOperandName() const
{
    return top.first.name;
}

template<size_t dim>
csmp::VARIABLE_TYPE   MathOperatorRHS<dim>::TestOperandType() const
{
    return top.first.key.type;
}

template<size_t dim>
csmp::PLACEMENT   MathOperatorRHS<dim>::TestOperandPlacement() const
{
    return top.first.key.place;
}

template<size_t dim>
size_t   MathOperatorRHS<dim>::TestOperandDataDepth() const
{
    return top.first.key.dataDepth;
}

template<size_t dim>
void MathOperatorRHS<dim>::TestOperandOffset( size_t os )
{
    top.second = os;
}

template<size_t dim>
size_t MathOperatorRHS<dim>::TestOperandOffset() const
{
    return top.second;
}





/// ===================================  Accumulation Process Settings ========================================================================



template<size_t dim>
void MathOperatorRHS<dim>::AddAccumulate()
{
   add_accumulate_            = true;
   subtract_accumulate_       = false;
   add_accumulate_later_      = false;
   subtract_accumulate_later_ = false;
   multiply_accumulate_       = false;

}

template<size_t dim>
void MathOperatorRHS<dim>::SubtractAccumulate()
{
   add_accumulate_            = false;
   subtract_accumulate_       = true;
   add_accumulate_later_      = false;
   subtract_accumulate_later_ = false;
   multiply_accumulate_       = false;

}

template<size_t dim>
void MathOperatorRHS<dim>::AddAccumulateLater()
{
   add_accumulate_            = false;
   subtract_accumulate_       = false;
   add_accumulate_later_      = true;
   subtract_accumulate_later_ = false;
   multiply_accumulate_       = false;

}

template<size_t dim>
void MathOperatorRHS<dim>::SubtractAccumulateLater()
{
   add_accumulate_            = false;
   subtract_accumulate_       = false;
   add_accumulate_later_      = false;
   subtract_accumulate_later_ = true;
   multiply_accumulate_       = false;

}

template<size_t dim>
void MathOperatorRHS<dim>::MultiplyAccumulate()
{
   add_accumulate_            = false;
   subtract_accumulate_       = false;
   add_accumulate_later_      = false;
   subtract_accumulate_later_ = false;
   multiply_accumulate_       = true;
}

template<size_t dim>
void MathOperatorRHS<dim>::LumpedFormulation( bool lumped )
{
    lump_matrices_=lumped;
}

template<size_t dim>
void  MathOperatorRHS<dim>::ApplicationCycles( size_t c )
{
    application_cycles_=c;
}

template<size_t dim>
void  MathOperatorRHS<dim>::ApplicationCycle( size_t c )
{
    application_cycle_=c;
}

template<size_t dim>
void MathOperatorRHS<dim>::MultiplyBy( double64 integral_mult_factor )
{
    factor_ = integral_mult_factor;
}

template<size_t dim>
void MathOperatorRHS<dim>::MultiplyWithTimeIncrement( bool multiply )
{
    time_multiply_  =   multiply;
}

template<size_t dim>
void MathOperatorRHS<dim>::DivideByTimeIncrement( bool divide )
{
    time_divide_  =   divide;
}

template<size_t dim>
 bool MathOperatorRHS<dim>::Add() const
{
    return add_accumulate_;
}

template<size_t dim>
 bool MathOperatorRHS<dim>::Subtract() const
{
    return subtract_accumulate_;
}

template<size_t dim>
 bool MathOperatorRHS<dim>::AddLater() const
{
    return add_accumulate_later_;
}

template<size_t dim>
 bool MathOperatorRHS<dim>::SubtractLater() const
{
    return subtract_accumulate_later_;
}

template<size_t dim>
bool MathOperatorRHS<dim>::Multiply() const
{
    return multiply_accumulate_;
}

template<size_t dim>
 size_t  MathOperatorRHS<dim>::ApplicationCycles() const
{
    return application_cycles_;
}

template<size_t dim>
 size_t  MathOperatorRHS<dim>::ApplicationCycle() const
{
    return application_cycle_;
}

template<size_t dim>
 bool  MathOperatorRHS<dim>::LumpedFormulation() const
{
    return lump_matrices_;
}

template<size_t dim>
double64   MathOperatorRHS<dim>::MultiplyBy() const
{
    return factor_;
}

template<size_t dim>
 bool MathOperatorRHS<dim>::MultiplyWithTimeIncrement() const
{
    return time_multiply_;
}

 template<size_t dim>
  bool MathOperatorRHS<dim>::DivideByTimeIncrement() const
 {
     return time_divide_;
 }


/// ===================================  Main Computational Procedure Functions ========================================================================


// interpolation of property if isoparametric elements are used
template<size_t dim>
void MathOperatorRHS<dim>::PropertyAtIntegrationPoint( const Element<dim>& e_ref,
                                                       const csmp::Index& idx,
                                                       size_t ip,
                                                       DenseMatrix<DM_MIN>& M )
 {
    if ( !e_ref.UsesLocalCoordinates() )
       throw csmp::Exception( CSMP_FATAL_ERROR,
                              "MathOperatorRHS<dim>::PropertyValueAtIntegrationPoint",
                              "Current element does not support numerical integration.");

     // 0. If the property is a constraint point variable a one-to-one mapping
     //    can be performed; no interpolation is needed
     if ( idx.place == ELEMENT_INTEGRATION_POINT ) {
          if ( idx.type == SCALAR ) {
               M.AssignToDiagonal( dim, e_ref.Read( ip, idx ) );
            }
          else if ( idx.type == VECTOR ) {
               e_ref.Read( ip, idx, VC[0] );
               M.AssignToDiagonal( VC[0] );
            }
          else if ( idx.type == TENSOR ) {
               e_ref.Read( ip, idx, TS[0] );
               M = TS[0];
           }
          else if ( idx.type == ARRAY )
            {
              AR[0].Resize( idx.dataDepth );
              e_ref.Read( ip, idx, AR[0] );
              M.AssignToDiagonal( AR[0] );
            }
          else if ( idx.type == FLAGGEDARRAY )
            {
              FR[0].Resize( idx.dataDepth );
              e_ref.Read( ip, idx, FR[0] );
              M.AssignToDiagonal( FR[0] );
            }
          return;
       }

     // 1. If the property is an element variable it is constant over the
     //    element, such that no interpolation is needed
     if ( idx.place == ELEMENT || idx.place == REGION )
       {
          if ( idx.type == SCALAR ) {
               M.AssignToDiagonal( dim, e_ref.Read(  idx ) );
            }
          else if ( idx.type == VECTOR ) {
               e_ref.Read(  idx, VC[0] );
               M.AssignToDiagonal( VC[0] );
            }
          else if ( idx.type == TENSOR ) {
               e_ref.Read(  idx, TS[0] );
               M = TS[0];
           }
          else if ( idx.type == ARRAY )
            {
               AR[0].Resize( idx.dataDepth );
               e_ref.Read(  idx, AR[0] );
               M.AssignToDiagonal( AR[0] );
            }
          else if ( idx.type == FLAGGEDARRAY )
            {
               FR[0].Resize( idx.dataDepth );
               e_ref.Read(  idx, FR[0] );
               M.AssignToDiagonal( FR[0] );
            }
          return;
       }

     // 2. If the property is a node property, it must be interpolated to the
     //    integration point
     if ( idx.place == NODE )
       {
          // interpolating properties
          if ( idx.type == SCALAR ) {
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, SC[0] );
               M.AssignToDiagonal( dim, SC[0] );
            }
          else if ( idx.type == VECTOR ) {
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, VC[0] );
               M.AssignToDiagonal( VC[0] );
            }
          else if ( idx.type == TENSOR ) {
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, TS[0] );
               M = TS[0];
           }
          else if ( idx.type == ARRAY )
            {
               AR[0].Resize( idx.dataDepth );
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, AR[0] );
               M.AssignToDiagonal( AR[0] );
            }
          else if ( idx.type == FLAGGEDARRAY )
            {
               FR[0].Resize( idx.dataDepth );
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, FR[0] );
               M.AssignToDiagonal( FR[0] );
            }
          return;
       }

     throw csmp::Exception( CSMP_FATAL_ERROR,
                            "MathOperatorRHS<dim>::PropertyValueAtIntegrationPoint",
                            "This method cannot interpolate IntegrationPoint variables.");

 } // end PropertyValueAtIntegrationPoint( Element )

template<size_t dim>
void MathOperatorRHS<dim>::PropertyAtIntegrationPoint( const Face<dim>& e_ref,
                                                       const csmp::Index& idx,
                                                       size_t ip,
                                                       DenseMatrix<DM_MIN>& M )
 {
     if ( !e_ref.UsesLocalCoordinates() )
       throw csmp::Exception( CSMP_FATAL_ERROR,
                              "MathOperatorRHS<dim>::PropertyValueAtIntegrationPoint",
                              "Current element does not support numerical integration.");

     // 0. If the property is a constraint point variable a one-to-one mapping
     //    can be performed; no interpolation is needed
     if ( idx.place == FACE_INTEGRATION_POINT ) {
          if ( idx.type == SCALAR ) {
               M.AssignToDiagonal( dim, e_ref.Read( ip, idx ) );
            }
          else if ( idx.type == VECTOR ) {
               e_ref.Read( ip, idx, VC[0] );
               M.AssignToDiagonal( VC[0] );
            }
          else if ( idx.type == TENSOR ) {
               e_ref.Read( ip, idx, TS[0] );
               M = TS[0];
           }
          else if ( idx.type == ARRAY )
            {
              AR[0].Resize( idx.dataDepth );
              e_ref.Read( ip, idx, AR[0] );
              M.AssignToDiagonal( AR[0] );
            }
          else if ( idx.type == FLAGGEDARRAY )
            {
              FR[0].Resize( idx.dataDepth );
              e_ref.Read( ip, idx, FR[0] );
              M.AssignToDiagonal( FR[0] );
            }
          return;
       }

     // 1. If the property is an element variable it is constant over the
     //    element, such that no interpolation is needed
     if ( idx.place == FACE || idx.place == BOUNDARY )
       {
          if ( idx.type == SCALAR ) {
               M.AssignToDiagonal( dim, e_ref.Read(  idx ) );
            }
          else if ( idx.type == VECTOR ) {
               e_ref.Read(  idx, VC[0] );
               M.AssignToDiagonal( VC[0] );
            }
          else if ( idx.type == TENSOR ) {
               e_ref.Read(  idx, TS[0] );
               M = TS[0];
           }
          else if ( idx.type == ARRAY )
            {
               AR[0].Resize( idx.dataDepth );
               e_ref.Read(  idx, AR[0] );
               M.AssignToDiagonal( AR[0] );
            }
          else if ( idx.type == FLAGGEDARRAY )
            {
               FR[0].Resize( idx.dataDepth );
               e_ref.Read(  idx, FR[0] );
               M.AssignToDiagonal( FR[0] );
            }
          return;
       }

     // 2. If the property is a node property, it must be interpolated to the
     //    integration point
     if ( idx.place == NODE )
       {
          // interpolating properties
          if ( idx.type == SCALAR ) {
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, SC[0] );
               M.AssignToDiagonal( dim, SC[0] );
            }
          else if ( idx.type == VECTOR ) {
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, VC[0] );
               M.AssignToDiagonal( VC[0] );
            }
          else if ( idx.type == TENSOR ) {
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, TS[0] );
               M = TS[0];
           }
          else if ( idx.type == ARRAY )
            {
               AR[0].Resize( idx.dataDepth );
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, AR[0] );
               M.AssignToDiagonal( AR[0] );
            }
          else if ( idx.type == FLAGGEDARRAY )
            {
               FR[0].Resize( idx.dataDepth );
               e_ref.PropertyValueAtIntegrationPoint(  idx, ip, FR[0] );
               M.AssignToDiagonal( FR[0] );
            }
          return;
       }

     throw csmp::Exception( CSMP_FATAL_ERROR,
                            "MathOperatorRHS<dim>::PropertyValueAtIntegrationPoint",
                           "This method cannot interpolate IntegrationPoint variables.");

 } // end PropertyValueAtIntegrationPoint( Face )






/**
GetOperands fills a vector of material property matrices with the 
required values for later computation. These matrices always have 
the dimensions spatial-dimension^2 and they will hold either scalar,
vector or tensor properties, depending on what kind of property the 
Operand is.  

When the property is an element property, it will be put into the
first vector entry MTRL[0].*/
template<size_t dim>
void MathOperatorRHS<dim>::GetOperands( Element<dim>& e_ref )
 {
    // if operand property is an element property
   if ( MaterialOperandPlacement() == ELEMENT || MaterialOperandPlacement() == REGION )
      {
         MTRL.resize(1U);
         if ( MaterialOperandType() == SCALAR )
           MTRL[0].AssignToDiagonal( dim, e_ref.Read(  MaterialOperandKey() ) );
         else if ( MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e_ref.Read(  MaterialOperandKey(), vc );
              MTRL[0].AssignToDiagonal( vc );
           }
         else if ( MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e_ref.Read(  MaterialOperandKey(), ts );
              MTRL[0] = ts;
           }
         else if ( MaterialOperandType() == ARRAY )
           {
              ArrayVariable ar( MaterialOperandDataDepth() );
              e_ref.Read(  MaterialOperandKey(), ar );
              MTRL[0].AssignToDiagonal( ar );
           }
         else if ( MaterialOperandType() == FLAGGEDARRAY )
           {
              FlaggedArrayVariable fr( MaterialOperandDataDepth() );
              e_ref.Read(  MaterialOperandKey(), fr );
              MTRL[0].AssignToDiagonal( fr );
           }
      }
    else // if the operand is placed on the constraint-points
      {
         if ( MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ) {
               MTRL.resize( e_ref.FE()->IntegrationPoints() );
               const size_t n_integration_points( e_ref.FE()->IntegrationPoints() );
               for ( size_t i=0U; i<n_integration_points; i++ )
                  {
                     if ( MaterialOperandType() == SCALAR )
                       MTRL[i].AssignToDiagonal( dim, e_ref.Read( i, MaterialOperandKey() ) );
                     else if ( MaterialOperandType() == VECTOR ) {
                          VectorVariable<dim>  vc;
                          e_ref.Read( i, MaterialOperandKey(), vc );
                          MTRL[i].AssignToDiagonal( vc );
                       }
                     else if ( MaterialOperandType() == TENSOR ) {
                          TensorVariable<dim>  ts;
                          e_ref.Read( i, MaterialOperandKey(), ts );
                          MTRL[i] = ts;
                       }
                     else if ( MaterialOperandType() == ARRAY )
                       {
                          ArrayVariable ar( MaterialOperandDataDepth() );
                          e_ref.Read( i, MaterialOperandKey(), ar );
                          MTRL[0].AssignToDiagonal( ar );
                       }
                     else if ( MaterialOperandType() == FLAGGEDARRAY )
                       {
                          FlaggedArrayVariable fr( MaterialOperandDataDepth() );
                          e_ref.Read( i, MaterialOperandKey(), fr );
                          MTRL[0].AssignToDiagonal( fr );
                       }
                 }
           }

         // if the operand is placed on the node
         else if ( MaterialOperandPlacement() == NODE ) {
             const size_t n_integration_points( e_ref.FE()->IntegrationPoints() );
             for ( size_t i=0U; i<n_integration_points; i++ )
                  PropertyAtIntegrationPoint(  e_ref, MaterialOperandKey(), i, MTRL[i] );
           }
         else
             throw csmp::Exception( CSMP_FATAL_ERROR,
                                    "MathOperatorRHS<dim>::GetOperands(Element):",
                                    "Face based operands cannot be accumulated with this method");
      }

 } // end GetOperands(Element)

template<size_t dim>
void MathOperatorRHS<dim>::GetOperands( Face<dim>&  e_ref )
 {
    // if operand property is an element property
   if ( MaterialOperandPlacement() == FACE || MaterialOperandPlacement() == BOUNDARY )
      {
         MTRL.resize(1U);
         if ( MaterialOperandType() == SCALAR )
           MTRL[0].AssignToDiagonal( dim, e_ref.Read(  MaterialOperandKey() ) );
         else if ( MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e_ref.Read(  MaterialOperandKey(), vc );
              MTRL[0].AssignToDiagonal( vc );
           }
         else if ( MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e_ref.Read(  MaterialOperandKey(), ts );
              MTRL[0] = ts;
           }
         else if ( MaterialOperandType() == ARRAY )
           {
              ArrayVariable ar( MaterialOperandDataDepth() );
              e_ref.Read(  MaterialOperandKey(), ar );
              MTRL[0].AssignToDiagonal( ar );
           }
         else if ( MaterialOperandType() == FLAGGEDARRAY )
           {
              FlaggedArrayVariable fr( MaterialOperandDataDepth() );
              e_ref.Read(  MaterialOperandKey(), fr );
              MTRL[0].AssignToDiagonal( fr );
           }
      }
    else // if the operand is placed on the constraint-points
      {
         if ( MaterialOperandPlacement() == FACE_INTEGRATION_POINT ) {
               MTRL.resize( e_ref.FE()->IntegrationPoints() );
               const size_t n_integration_points( e_ref.FE()->IntegrationPoints() );
               for ( size_t i=0U; i<n_integration_points; i++ )
                  {
                     if ( MaterialOperandType() == SCALAR )
                       MTRL[i].AssignToDiagonal( dim, e_ref.Read( i, MaterialOperandKey() ) );
                     else if ( MaterialOperandType() == VECTOR ) {
                          VectorVariable<dim>  vc;
                          e_ref.Read( i, MaterialOperandKey(), vc );
                          MTRL[i].AssignToDiagonal( vc );
                       }
                     else if ( MaterialOperandType() == TENSOR ) {
                          TensorVariable<dim>  ts;
                          e_ref.Read( i, MaterialOperandKey(), ts );
                          MTRL[i] = ts;
                       }
                     else if ( MaterialOperandType() == ARRAY )
                       {
                          ArrayVariable ar( MaterialOperandDataDepth() );
                          e_ref.Read( i, MaterialOperandKey(), ar );
                          MTRL[0].AssignToDiagonal( ar );
                       }
                     else if ( MaterialOperandType() == FLAGGEDARRAY )
                       {
                          FlaggedArrayVariable fr( MaterialOperandDataDepth() );
                          e_ref.Read( i, MaterialOperandKey(), fr );
                          MTRL[0].AssignToDiagonal( fr );
                       }
                 }
           }

         // if the operand is placed on the node
         else if ( MaterialOperandPlacement() == NODE ) {
             const size_t n_integration_points( e_ref.FE()->IntegrationPoints() );
             for ( size_t i=0U; i<n_integration_points; i++ )
                  PropertyAtIntegrationPoint(  e_ref, MaterialOperandKey(), i, MTRL[i] );
           }
         else
             throw csmp::Exception( CSMP_FATAL_ERROR,
                                    "MathOperatorRHS<dim>::GetOperands(Face):",
                                    "Element based operands cannot be accumulated with this method");
      }

 } // end GetOperands(Face)

template<size_t dim>
void MathOperatorRHS<dim>::GetOperands( InterFace<dim>& f )
 {
    std::cerr <<"\nMathOperatorRHS<dim>::GetOperands(InterFace): ";
    std::cerr <<" Overload to get LHS operands from interface: "<< f.Idx() << std::endl;
    throw invalid_argument("MathOperatorRHS<dim>::GetOperands(InterFace)");
 } // end GetOperands(InterFace)


template<size_t dim>
void MathOperatorRHS<dim>::WriteOperands( Element<dim>& e )
 {
    std::cerr <<"\nMathOperatorRHS<dim>::WriteOperands(Element): ";
    std::cerr <<" Overload to write LHS operands from element: "<< e.Idx() << std::endl;
    throw invalid_argument("MathOperatorRHS<dim>::WriteOperands(Element)");
 } // end WriteOperands(Element)


template<size_t dim>
void MathOperatorRHS<dim>::WriteOperands( Face<dim>& f )
 {
    std::cerr <<"\nMathOperatorRHS<dim>::WriteOperands(Face): ";
    std::cerr <<" Overload to write LHS operands from face: "<< f.Idx() << std::endl;
    throw invalid_argument("MathOperatorRHS<dim>::WriteOperands(Face)");
 } // end WriteOperands(Face)

template<size_t dim>
void MathOperatorRHS<dim>::WriteOperands( InterFace<dim>& f )
 {
    std::cerr <<"\nMathOperatorRHS<dim>::WriteOperands(InterFace): ";
    std::cerr <<" Overload to write LHS operands from interface: "<< f.Idx() << std::endl;
    throw invalid_argument("MathOperatorRHS<dim>::WriteOperands(InterFace)");
 } // end WriteOperands(InterFace)

template<size_t dim>
void MathOperatorRHS<dim>::ComputeContribution( Element<dim>& e )
 {
    std::cerr <<"\nMathOperatorRHS<dim>::ComputeContribution(Element): ";
    std::cerr <<" Overload to calculate LHS contribution from Element: "<< e.Idx() << std::endl;
    throw invalid_argument("MathOperatorRHS<dim>::ComputeContribution(Element)");
 }// end ComputeContribution(Element)

template<size_t dim>
void MathOperatorRHS<dim>::ComputeContribution( Face<dim>& f )
 {
    std::cerr <<"\nMathOperatorRHS<dim>::ComputeContribution(Face): ";
    std::cerr <<" Overload to calculate LHS constribution from face: "<< f.Idx() << std::endl;
    throw invalid_argument("MathOperatorRHS<dim>::ComputeContribution(Face)");
 }// end ComputeContribution(Face)

template<size_t dim>
void MathOperatorRHS<dim>::ComputeContribution( InterFace<dim>& f )
 {
    std::cerr <<"\nMathOperatorRHS<dim>::ComputeContribution(InterFace): ";
    std::cerr <<" Overload to calculate LHS constribution from interface: "<< f.Idx() << std::endl;
    throw invalid_argument("MathOperatorRHS<dim>::ComputeContribution(InterFace)");
 }// end ComputeContribution(InterFace)



template<size_t dim>
void MathOperatorRHS<dim>::MultiplyWithTimeFactor( double64 dt )
{
    for ( size_t i=0; i<RHS.size(); i++ )
        RHS[i] *= dt;
}


/// AssignToGlobal matrix functions

template<size_t dim>
void MathOperatorRHS<dim>::AssignToGlobal( const Element<dim>& e, vector<double64>& rhs )
{
    // map local to global indexes for test(basic) operands
    const size_t nodes(e.Nodes());
    IDT.resize( nodes );
    for ( size_t i=0U; i<nodes; i++ )
       IDT[i] = e.N(i)->Idx();

    if ( TestOperandType() != SCALAR )
      transformNodeIndexVector( dim, TestOperandKey(), IDT );

    const size_t dof(IDT.size());

    for ( size_t i=0U; i<dof; i++ )
       IDT[i] += this->TestOperandOffset();

   // perform assignment from local matrix to global matrix
   if ( multiply_accumulate_ )
     for ( size_t i=0U; i<dof; i++ )
       rhs[ IDT[i] ] *= RHS[i] * factor_;
   
   else if ( add_accumulate_ || add_accumulate_later_ )
     for ( size_t i=0U; i<dof; i++ )
       rhs[ IDT[i] ] += RHS[i] * factor_;

   else if ( subtract_accumulate_ || subtract_accumulate_later_ )
     for ( size_t i=0U; i<dof; i++ )
       rhs[ IDT[i] ] -= RHS[i] * factor_;
   else
       throw csmp::Exception( CSMP_ERROR,
                              "MathOperatorRHS<dim>::AssignToGlobal(Element)",
                              "accumulation instructions could not be parsed.");
}

template<size_t dim>
void MathOperatorRHS<dim>::AssignToGlobal( const Face<dim>& e, vector<double64>& rhs )
{
    // map local to global indexes for test(basic) operands
    const size_t nodes(e.Nodes());

    IDT.resize( nodes );
    for ( size_t i=0U; i<nodes; i++ )
       IDT[i] = e.N(i)->Idx();

    if ( TestOperandType() != SCALAR )
      transformNodeIndexVector( dim, TestOperandKey(), IDT );

    const size_t dof(IDT.size());

    for ( size_t i=0U; i<dof; i++ )
       IDT[i] += this->TestOperandOffset();

   // perform assignment from local matrix to global matrix

   if ( multiply_accumulate_ )
     for ( size_t i=0U; i<dof; i++ )
       rhs[ IDT[i] ] *= RHS[i] * factor_;

   else if ( add_accumulate_ || add_accumulate_later_ )
     for ( size_t i=0U; i<dof; i++ )
       rhs[ IDT[i] ] += RHS[i] * factor_;

   else if ( subtract_accumulate_ || subtract_accumulate_later_ )
     for ( size_t i=0U; i<dof; i++ )
       rhs[ IDT[i] ] -= RHS[i] * factor_;
   else
       throw csmp::Exception( CSMP_ERROR,
                              "MathOperatorRHS<dim>::AssignToGlobal(Face)",
                              "accumulation instructions could not be parsed.");
}

template<size_t dim>
void MathOperatorRHS<dim>::AssignToGlobal( const InterFace<dim>& f, vector<double64>& rhs )
{
    std::cerr <<"\nMathOperatorRHS<dim>::AssignToGlobal(InterFace): ";
    std::cerr <<" Overload to assign RHS  local entries to global matrix: "<< f.Idx() << std::endl;
    throw invalid_argument("MathOperatorRHS<dim>::AssignToGlobal(InterFace)");
}


/// Transform Node Index Vector

/// Roman, 2013: Added Array Variable and Falgged Array Variable Index Accessor

/**

If an output variable from the element is a VECTOR,TENSOR,ARRAY or FLAGGEDARRAY property,
a mapping is applied such that the number of entries that were read at the element are
expanded out into the N vector. This accomodates the additional degrees of freedom
in the global solution matrix.

@param N The result is returned into the input vector. The final size
of this vector<size_t> will be

nodes x dim for VECTOR,
nodes x dim x dim for TENSOR,
nodes x ArrayLength for ARRAY and FLAGGEDARRAY

@section implementation Implementation

If the output variable is a vector property a mapping is applied such that
the number of entries that were read are expanded out into
the N vector. This works independently of the placement of the variable
looping backwards such that no temporary storage for the variables is
needed.

@section application Application

Used in the assembly of global solution matrices where the dependent variable
is a vector or tensor property.
*/

void transformNodeIndexVector( size_t dim, const csmp::Index& idx, vector<size_t>& N )
 {
    assert( !N.empty() );
    assert( idx.type != SCALAR );

    // creating a running index for decrementation
    const int32  oldNsize_m1(static_cast<int32>(N.size()) - 1);

    // VECTOR variables
    if ( idx.type == VECTOR ) {
          // new elements are set to zero
          N.resize( N.size() * dim, 0U ); // new size of N vector
          for ( int32 i=oldNsize_m1, k=static_cast<int32>(N.size())-1; i>=0; i-- )
            for ( size_t j=0U; j<dim; j++ )
              N[ static_cast<size_t>(k--) ] = ((N[ static_cast<size_t>(i) ]+1) * dim - j) - 1;
          return;
      }

    // TENSOR variables
    if ( idx.type == TENSOR ) {
        const size_t  dim2(dim*dim);
        N.resize( N.size() * dim2, 0U );
        for ( int32 i=oldNsize_m1, k=static_cast<int32>(N.size())-1; i>=0; i-- )
          for ( size_t j=0; j<dim2; j++ )
            N[ static_cast<size_t>(k--) ] = ((N[ static_cast<size_t>(i) ]+1) * dim2 - j) - 1;
    }

    // ARRAY or FLAGGEDARRAY variables
    if ( idx.type == ARRAY || idx.type == FLAGGEDARRAY ) {
          const long64 length( idx.dataDepth );
          // new elements are set to zero
          N.resize( N.size() * length, 0U ); // new size of N vector ( length of array * number of nodes )
          for ( long64 i=oldNsize_m1, k=static_cast<long64>(N.size())-1; i>=0; i-- )
            // cycle for one node through all the array variable components
            for ( size_t j=0U; j<length; j++ )
              N[ static_cast<size_t>(k--) ] = ((N[ static_cast<size_t>(i) ]+1) * length - j) - 1;
          return;
      }

} // end transformNodeIndexVector

template class MathOperatorRHS<1U>;
template class MathOperatorRHS<2U>;
template class MathOperatorRHS<3U>;

} // end namespace csp

