#include "MathOperatorRHS.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Exception.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
  MathOperatorRHS<dim,CELL>::MathOperatorRHS()
    : name_("unspecified RHS-operator"),
      MTRL(13, DenseMatrix<DM_MIN>(2, 2)),
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




template<uint32_t dim, template<uint32_t> class CELL>
  MathOperatorRHS<dim,CELL>::MathOperatorRHS( const PropertyDatabase<dim>& pref,
                                              const char* test)
    : name_("unspecified RHS-operator"),
      top(make_pair(pref.Parameter(test), 0)),
      MTRL(13, DenseMatrix<DM_MIN>(2, 2)),
      DERIV(2, 3),
      IPOL(3),
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




template<uint32_t dim, template<uint32_t> class CELL>
  MathOperatorRHS<dim,CELL>::MathOperatorRHS( const PropertyDatabase<dim>& pref,
                                              const char* oper,
                                              const char* test )
    : name_("unspecified RHS-operator"),
      op(pref.Parameter(oper)),
      top(make_pair(pref.Parameter(test), 0)),
      MTRL(13, DenseMatrix<DM_MIN>(2, 2)),
      DERIV(2, 3),
      IPOL(3),
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



template<uint32_t dim, template<uint32_t> class CELL>
  MathOperatorRHS<dim,CELL>::MathOperatorRHS(const MathOperatorRHS<dim,CELL>& mo)
    : name_(mo.name_),
      op(mo.op),    // basic operand
      top(mo.top),   // testfunction operand
      RHS(mo.RHS),   // solution vector to be accumulated
      MTRL(mo.MTRL),  // basic Operand storage
      DERIV(mo.DERIV),
      IPOL(mo.IPOL),
      factor_(mo.factor_),
      add_accumulate_(mo.add_accumulate_),
      subtract_accumulate_(mo.subtract_accumulate_),
      add_accumulate_later_(mo.add_accumulate_later_),
      subtract_accumulate_later_(mo.subtract_accumulate_later_),
      multiply_accumulate_(mo.multiply_accumulate_),
      lump_matrices_(mo.lump_matrices_),
      application_cycles_(mo.application_cycles_),
      application_cycle_(mo.application_cycle_),
      time_multiply_(mo.time_multiply_),
      time_divide_(mo.time_divide_)
  {
  }



template<uint32_t dim, template<uint32_t> class CELL>
  MathOperatorRHS<dim,CELL>& MathOperatorRHS<dim,CELL>::operator=(const MathOperatorRHS<dim,CELL>& mo)
  {
    if (this != &mo) {
      name_ = mo.name_;
      op = mo.op;        // basic operand
      top = mo.top;       // testfunction operand
      RHS = mo.RHS;       // solution vector to be accumulated
      MTRL = mo.MTRL;      // basic Operand storage
      DERIV = mo.DERIV;
      IPOL = mo.IPOL;
      add_accumulate_ = mo.add_accumulate_;
      subtract_accumulate_ = mo.subtract_accumulate_;
      add_accumulate_later_ = mo.add_accumulate_later_;
      subtract_accumulate_later_ = mo.subtract_accumulate_later_;
      multiply_accumulate_ = mo.multiply_accumulate_;
      lump_matrices_ = mo.lump_matrices_;
      application_cycles_ = mo.application_cycles_;
      application_cycle_ = mo.application_cycle_;
      time_multiply_ = mo.time_multiply_;
      time_divide_ = mo.time_divide_;
      factor_ = mo.factor_;
    }
    return *this;
  }




  /// ===================================  Operand Functions ========================================================================

template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::Name(const char* s, const char* topname) noexcept
  {
    name_ = s;
    name_ += ": Operand: '";
    name_ += topname;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::Name(const char* s, const char* operand_name, const char* topname) noexcept
  {
    name_ = s;
    name_ += ": Operand: '";
    name_ += operand_name;
    name_ += "', '";
    name_ += topname;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  string  MathOperatorRHS<dim,CELL>::Name() const noexcept
  {
    return name_;
  }






template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::Out() const
  {
    op.Out();
    top.first.Out();
  }




  // MATERIAL PROPERTY OPERAND

template<uint32_t dim, template<uint32_t> class CELL>
  const csmp::Parameter&   MathOperatorRHS<dim,CELL>::MaterialOperand() const noexcept
  {
    return op;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  const csmp::Index&   MathOperatorRHS<dim,CELL>::MaterialOperandKey() const noexcept
  {
    return op.key;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  string  MathOperatorRHS<dim,CELL>::MaterialOperandName() const noexcept
  {
    return op.name;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  csmp::VARIABLE_TYPE   MathOperatorRHS<dim,CELL>::MaterialOperandType() const noexcept
  {
    return op.key.type;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  csmp::PLACEMENT   MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() const noexcept
  {
    return op.key.place;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  uint32_t   MathOperatorRHS<dim,CELL>::MaterialOperandDataDepth() const noexcept
  {
    return op.key.dataDepth;
  }


  // BASIC FUNCTION OPERAND




template<uint32_t dim, template<uint32_t> class CELL>
  const csmp::Parameter&  MathOperatorRHS<dim,CELL>::BasicOperand() const noexcept
  {
    return top.first;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  const csmp::Index&   MathOperatorRHS<dim,CELL>::BasicOperandKey() const noexcept
  {
    return top.first.key;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  string  MathOperatorRHS<dim,CELL>::BasicOperandName() const noexcept
  {
    return top.first.name;
  }





template<uint32_t dim, template<uint32_t> class CELL>
  csmp::VARIABLE_TYPE   MathOperatorRHS<dim,CELL>::BasicOperandType() const noexcept
  {
    return top.first.key.type;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  csmp::PLACEMENT   MathOperatorRHS<dim,CELL>::BasicOperandPlacement() const noexcept
  {
    return top.first.key.place;
  }





template<uint32_t dim, template<uint32_t> class CELL>
  uint32_t   MathOperatorRHS<dim,CELL>::BasicOperandDataDepth() const noexcept
  {
    return top.first.key.dataDepth;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::BasicOperandOffset(size_t os) noexcept
  {
    top.second = os;
  }




template<uint32_t dim, template<uint32_t> class CELL>
size_t MathOperatorRHS<dim,CELL>::BasicOperandOffset() const noexcept
  {
    return top.second;
  }

  // TEST FUNCTION OPERAND


template<uint32_t dim, template<uint32_t> class CELL>
  const csmp::Parameter&   MathOperatorRHS<dim,CELL>::TestOperand() const noexcept
  {
    return top.first;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  const csmp::Index&   MathOperatorRHS<dim,CELL>::TestOperandKey() const noexcept
  {
    return top.first.key;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  string  MathOperatorRHS<dim,CELL>::TestOperandName() const noexcept
  {
    return top.first.name;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  csmp::VARIABLE_TYPE   MathOperatorRHS<dim,CELL>::TestOperandType() const noexcept
  {
    return top.first.key.type;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  csmp::PLACEMENT   MathOperatorRHS<dim,CELL>::TestOperandPlacement() const noexcept
  {
    return top.first.key.place;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  uint32_t   MathOperatorRHS<dim,CELL>::TestOperandDataDepth() const noexcept
  {
    return top.first.key.dataDepth;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::TestOperandOffset(size_t os) noexcept
  {
    top.second = os;
  }



template<uint32_t dim, template<uint32_t> class CELL>
size_t MathOperatorRHS<dim,CELL>::TestOperandOffset() const noexcept
  {
    return top.second;
  }





  /// ===================================  Accumulation Process Settings ========================================================================



template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::AddAccumulate() noexcept
  {
    add_accumulate_ = true;
    subtract_accumulate_ = false;
    add_accumulate_later_ = false;
    subtract_accumulate_later_ = false;
    multiply_accumulate_ = false;

  }




template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::SubtractAccumulate() noexcept
  {
    add_accumulate_ = false;
    subtract_accumulate_ = true;
    add_accumulate_later_ = false;
    subtract_accumulate_later_ = false;
    multiply_accumulate_ = false;

  }




template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::AddAccumulateLater() noexcept
  {
    add_accumulate_ = false;
    subtract_accumulate_ = false;
    add_accumulate_later_ = true;
    subtract_accumulate_later_ = false;
    multiply_accumulate_ = false;

  }




template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::SubtractAccumulateLater() noexcept
  {
    add_accumulate_ = false;
    subtract_accumulate_ = false;
    add_accumulate_later_ = false;
    subtract_accumulate_later_ = true;
    multiply_accumulate_ = false;

  }




template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::MultiplyAccumulate() noexcept
  {
    add_accumulate_ = false;
    subtract_accumulate_ = false;
    add_accumulate_later_ = false;
    subtract_accumulate_later_ = false;
    multiply_accumulate_ = true;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::LumpedFormulation(bool lumped) noexcept
  {
    lump_matrices_ = lumped;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  void  MathOperatorRHS<dim,CELL>::ApplicationCycles(uint32_t c) noexcept
  {
    application_cycles_ = c;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  void  MathOperatorRHS<dim,CELL>::ApplicationCycle(uint32_t c) noexcept
  {
    application_cycle_ = c;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::MultiplyBy(double integral_mult_factor) noexcept
  {
    factor_ = integral_mult_factor;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::MultiplyWithTimeIncrement(bool multiply) noexcept
  {
    time_multiply_ = multiply;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::DivideByTimeIncrement(bool divide) noexcept
  {
    time_divide_ = divide;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  bool MathOperatorRHS<dim,CELL>::Add() const noexcept
  {
    return add_accumulate_;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  bool MathOperatorRHS<dim,CELL>::Subtract() const noexcept
  {
    return subtract_accumulate_;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  bool MathOperatorRHS<dim,CELL>::AddLater() const noexcept
  {
    return add_accumulate_later_;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  bool MathOperatorRHS<dim,CELL>::SubtractLater() const noexcept
  {
    return subtract_accumulate_later_;
  }




template<uint32_t dim, template<uint32_t> class CELL>
  bool MathOperatorRHS<dim,CELL>::Multiply() const noexcept
  {
    return multiply_accumulate_;
  }



template<uint32_t dim, template<uint32_t> class CELL>
uint32_t  MathOperatorRHS<dim,CELL>::ApplicationCycles() const noexcept
  {
    return application_cycles_;
  }



template<uint32_t dim, template<uint32_t> class CELL>
uint32_t  MathOperatorRHS<dim,CELL>::ApplicationCycle() const noexcept
  {
    return application_cycle_;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  bool  MathOperatorRHS<dim,CELL>::LumpedFormulation() const noexcept
  {
    return lump_matrices_;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  double   MathOperatorRHS<dim,CELL>::MultiplyBy() const noexcept
  {
    return factor_;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  bool MathOperatorRHS<dim,CELL>::MultiplyWithTimeIncrement() const noexcept
  {
    return time_multiply_;
  }



template<uint32_t dim, template<uint32_t> class CELL>
  bool MathOperatorRHS<dim,CELL>::DivideByTimeIncrement() const noexcept
  {
    return time_divide_;
  }


  /// ===================================  Main Computational Procedure Functions ========================================================================


  // interpolation of property if isoparametric elements are used
template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorRHS<dim,CELL>::PropertyAtIntegrationPoint( const CELL<dim>& e_ref,
                                                            const csmp::Index& idx,
                                                            uint32_t ip,
                                                            DenseMatrix<DM_MIN>& M )
  {
    if (!e_ref.UsesLocalCoordinates())
      throw csmp::Exception(FATAL_ERROR,
        "MathOperatorRHS<dim,CELL>::PropertyValueAtIntegrationPoint",
        "Current element does not support numerical integration.");

    // 0. If the property is a constraint point variable a one-to-one mapping
    //    can be performed; no interpolation is needed
    if (idx.place == ELEMENT_INTEGRATION_POINT) {
      if (idx.type == SCALAR) {
          M.AssignToDiagonalAndZeroOffDiagonal( dim, e_ref.Read(ip,idx) );
        }
      else if (idx.type == VECTOR) {
          // TODO: change LocalVariableStorage so that vector gets returned directly
          VectorVariable<dim> vc;
          e_ref.Read(ip, idx, vc );
          M.AssignToDiagonal( vc );
        }
      else if (idx.type == TENSOR) {
          TensorVariable<dim> ts;
          e_ref.Read(ip, idx, ts );
          M = ts;
        }
      else if (idx.type == ARRAY)
        {
           ArrayVariable ar( idx.dataDepth );
           e_ref.Read(ip, idx, ar );
           M.AssignToDiagonal( ar );
        }
      else if (idx.type == FLAGGEDARRAY)
        {
          FlaggedArrayVariable fr( idx.dataDepth );
          e_ref.Read(ip, idx, fr );
          M.AssignToDiagonal( fr );
        }
      return;
    }

    // 1. If the property is an element variable it is constant over the
    //    element, such that no interpolation is needed
    if (idx.place == ELEMENT || idx.place == REGION)
      {
        if (idx.type == SCALAR) {
            M.AssignToDiagonalAndZeroOffDiagonal(dim, e_ref.Read(idx));
          }
        else if (idx.type == VECTOR) {
            VectorVariable<dim> vc;
            e_ref.Read(idx, vc );
            M.AssignToDiagonal( vc );
          }
        else if (idx.type == TENSOR) {
            TensorVariable<dim> ts;
            e_ref.Read(idx, ts );
            M = ts;
          }
        else if (idx.type == ARRAY)
          {
            ArrayVariable ar( idx.dataDepth );
            e_ref.Read(idx, ar );
            M.AssignToDiagonal( ar );
          }
        else if (idx.type == FLAGGEDARRAY)
          {
            FlaggedArrayVariable fr( idx.dataDepth );
            e_ref.Read(idx, fr );
            M.AssignToDiagonal( fr );
          }
        return;
      }

    // 2. If the property is a node property, it must be interpolated to the
    //    integration point
    if (idx.place == NODE)
      {
        // interpolating properties
        if (idx.type == SCALAR) {
            ScalarVariable sc;
            e_ref.PropertyValueAtIntegrationPoint(idx, ip, sc );
            M.AssignToDiagonal(dim, sc );
          }
        else if (idx.type == VECTOR) {
            VectorVariable<dim> vc;
            e_ref.PropertyValueAtIntegrationPoint(idx, ip, vc );
            M.AssignToDiagonal( vc );
          }
        else if (idx.type == TENSOR) {
            TensorVariable<dim> ts;
            e_ref.PropertyValueAtIntegrationPoint(idx, ip, ts );
            M = ts;
          }
        else if (idx.type == ARRAY)
          {
            ArrayVariable ar( idx.dataDepth );
            e_ref.PropertyValueAtIntegrationPoint(idx, ip, ar );
            M.AssignToDiagonal( ar );
          }
        else if (idx.type == FLAGGEDARRAY)
          {
            FlaggedArrayVariable fr( idx.dataDepth );
            e_ref.PropertyValueAtIntegrationPoint(idx, ip, fr );
            M.AssignToDiagonal( fr );
          }
        return;
      }

    throw csmp::Exception(FATAL_ERROR,
      "MathOperatorRHS<dim,CELL>::PropertyValueAtIntegrationPoint",
      "This method cannot interpolate IntegrationPoint variables.");

  } // end PropertyValueAtIntegrationPoint( CELL )





 /**
 GetOperands fills a vector of material property matrices with the
 required values for later computation. These matrices always have
 the dimensions spatial-dimension^2 and they will hold either scalar,
 vector or tensor properties, depending on what kind of property the
 Operand is.

 When the property is an element property, it will be put into the
 first vector entry MTRL[0].*/
template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::GetOperands( const CELL<dim>& e_ref )
  {
    // if operand property is an element property
    if ( MaterialOperandPlacement() == ELEMENT || MaterialOperandPlacement() == FACE || MaterialOperandPlacement() == REGION )
      {
        MTRL.resize(1U);
        if (MaterialOperandType() == SCALAR)
            MTRL[0].AssignToDiagonalAndZeroOffDiagonal(dim, e_ref.Read(MaterialOperandKey()));
          else if (MaterialOperandType() == VECTOR) {
            VectorVariable<dim>  vc;
            e_ref.Read(MaterialOperandKey(), vc);
            MTRL[0].AssignToDiagonal(vc);
          }
        else if (MaterialOperandType() == TENSOR) {
            TensorVariable<dim>  ts;
            e_ref.Read(MaterialOperandKey(), ts);
            MTRL[0] = ts;
          }
        else if (MaterialOperandType() == ARRAY)
          {
            ArrayVariable ar(MaterialOperandDataDepth());
            e_ref.Read(MaterialOperandKey(), ar);
            MTRL[0].AssignToDiagonal(ar);
          }
        else if (MaterialOperandType() == FLAGGEDARRAY)
          {
            FlaggedArrayVariable fr(MaterialOperandDataDepth());
            e_ref.Read(MaterialOperandKey(), fr);
            MTRL[0].AssignToDiagonal(fr);
          }
      }
    else // if the operand is placed on the constraint-points
      {
        const auto n_integration_points{ e_ref.IntegrationPoints() };

        if (MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT || MaterialOperandPlacement() == FACE_INTEGRATION_POINT ) {
          MTRL.resize( n_integration_points );
          for (auto i{0U}; i < n_integration_points; i++)
          {
            if (MaterialOperandType() == SCALAR)
              MTRL[i].AssignToDiagonalAndZeroOffDiagonal(dim, e_ref.Read(i, MaterialOperandKey()));
            else if (MaterialOperandType() == VECTOR) {
              VectorVariable<dim>  vc;
              e_ref.Read(i, MaterialOperandKey(), vc);
              MTRL[i].AssignToDiagonal(vc);
            }
            else if (MaterialOperandType() == TENSOR) {
              TensorVariable<dim>  ts;
              e_ref.Read(i, MaterialOperandKey(), ts);
              MTRL[i] = ts;
            }
            else if (MaterialOperandType() == ARRAY)
            {
              ArrayVariable ar(MaterialOperandDataDepth());
              e_ref.Read(i, MaterialOperandKey(), ar);
              MTRL[0].AssignToDiagonal(ar);
            }
            else if (MaterialOperandType() == FLAGGEDARRAY)
            {
              FlaggedArrayVariable fr(MaterialOperandDataDepth());
              e_ref.Read(i, MaterialOperandKey(), fr);
              MTRL[0].AssignToDiagonal(fr);
            }
          }
        }

      // if the operand is placed on the node
      else if (MaterialOperandPlacement() == NODE) {
          MTRL.resize(n_integration_points);
          for (auto i{0U}; i < n_integration_points; i++)
            PropertyAtIntegrationPoint(e_ref, MaterialOperandKey(), i, MTRL[i]);
        }
      else
        throw csmp::Exception(FATAL_ERROR,
          "MathOperatorRHS<dim,CELL>::GetOperands(Element):",
          "InterFace based operands cannot be accumulated with this method");
    }

  } // end GetOperands(Element)


template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::WriteOperands(CELL<dim>& e)
  {
    cerr << "\nMathOperatorRHS<dim,CELL>::WriteOperands(CELL): ";
    cerr << " Overload to write LHS operands from element: " << e.Idx() << endl;
    throw invalid_argument("MathOperatorRHS<dim,CELL>::WriteOperands(Element)");
  } // end WriteOperands(Element)





template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorRHS<dim,CELL>::ComputeContribution( const CELL<dim>& e )
  {
    cerr << "\nMathOperatorRHS<dim,CELL>::ComputeContribution(CELL): ";
    cerr << " Overload to calculate LHS contribution from CELL: " << e.Idx() << endl;
    throw invalid_argument("MathOperatorRHS<dim,CELL>::ComputeContribution(CELL)");
  }// end ComputeContribution(Element)




template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorRHS<dim,CELL>::MultiplyWithTimeFactor( double dt ) noexcept
  {
    for (uint32_t i = 0; i < RHS.size(); i++)
      RHS[i] *= dt;
  }


  /// AssignToGlobal matrix functions

template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorRHS<dim,CELL>::AssignToGlobal(const CELL<dim>& e, vector<double>& rhs )
  {
		// vectors for mapping local to global node indices for test operand
		vector<size_t> IDT( e.Nodes() ); // ii
    for ( uint32_t n{0U}; n<e.Nodes(); ++n )
      IDT[n] = e.N(n)->Idx();

    if (TestOperandType() != SCALAR)
      transformNodeIndexVector( TestOperandKey(), IDT );

    const size_t dof(IDT.size());

    for (uint32_t i{0U}; i < dof; i++)
      IDT[i] += this->TestOperandOffset();

    // perform assignment from local matrix to global matrix
    if (multiply_accumulate_)
      for (uint32_t i{0U}; i < dof; i++)
        rhs[IDT[i]] *= RHS[i] * factor_;

    else if (add_accumulate_ || add_accumulate_later_)
      for (uint32_t i{0U}; i < dof; i++)
        rhs[IDT[i]] += RHS[i] * factor_;

    else if (subtract_accumulate_ || subtract_accumulate_later_)
      for (uint32_t i{0U}; i < dof; i++)
        rhs[IDT[i]] -= RHS[i] * factor_;
    else
      throw csmp::Exception(ERROR,
        "MathOperatorRHS<dim,CELL>::AssignToGlobal(Element)",
        "accumulation instructions could not be parsed.");
  }




  /// Transform Node Index Vector

  /// Roman, 2013: Added Array Variable and Falgged Array Variable Index Accessor


template<uint32_t dim, template<uint32_t> class CELL>
void  MathOperatorRHS<dim,CELL>::AssignToGlobal(const CELL<dim>& e, vector<double>& rhs, const vector<size_t>& DOF_indexes )
 {
// TODO: get static assert to pass by moving functionality for SplitBoundaries to subclasses
//    static_assert( is_same<CELL<dim>,InterFace<dim>>::value, "Override this method for InterFace is subclass");
    // map local to global indexes for test and basic operands
		// vectors for mapping local to global node indices for test and basic operands
		// vectors for mapping local to global node indices for test operand
		vector<size_t> IDT( e.Nodes() ); // ii
    for ( uint32_t n{0U}; n<e.Nodes(); ++n )
      IDT[n] = e.N(n)->Idx();

    if (TestOperandType() != SCALAR)
      transformNodeIndexVector( TestOperandKey(), IDT );

    const size_t dof(IDT.size());

    for (uint32_t i{0U}; i < dof; i++) {
      IDT[i] += this->TestOperandOffset();
       IDT[i] = DOF_indexes[IDT[i]]; //-> to the global index 
    }

    // perform assignment from local matrix to global matrix
    if (multiply_accumulate_) {
      for (uint32_t i{0U}; i < dof; i++)
        if (IDT[i] != NULL_IDX) {
          rhs[IDT[i]] *= RHS[i] * factor_;
        }
    }
    else if (add_accumulate_ || add_accumulate_later_){
      for (uint32_t i{0U}; i < dof; i++)
        if (IDT[i] != NULL_IDX) {
          rhs[IDT[i]] += RHS[i] * factor_;
        }
    }
    else if (subtract_accumulate_ || subtract_accumulate_later_) {
      for (uint32_t i{0U}; i < dof; i++)
        if (IDT[i] != NULL_IDX) {
          rhs[IDT[i]] -= RHS[i] * factor_;
        }
    }
    else {
      throw csmp::Exception(ERROR,
        "MathOperatorRHS<dim,CELL>::AssignToGlobal(CELL)",
        "accumulation instructions could not be parsed.");
    }
  }
  



  /**

  If an output variable from the element is a VECTOR, TENSOR, ARRAY or FLAGGEDARRAY property,
  a mapping is applied such that the number of entries that were read at the element are
  expanded out into the N vector. This accomodates the additional degrees of freedom
  in the global solution matrix.

  @param idx Test function operand
  @param N The result is returned into the input vector. The final size
  of this vector<uint32_t> will be

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

void transformNodeIndexVector(const csmp::Index& idx, std::vector<size_t>& N) noexcept {
    if (N.empty() || idx.type == SCALAR) return;

    const size_t num_nodes = N.size();
    const size_t length = static_cast<size_t>(idx.dataDepth);

    // 1. Resize upfront
    N.resize(num_nodes * length);

    // 2. Work backwards to avoid overwriting data we haven't read yet
    // Use ptrdiff_t to safely handle the decrement down to 0
    for (auto i = static_cast<std::ptrdiff_t>(num_nodes) - 1; i >= 0; --i) {
        size_t original_val = N[static_cast<size_t>(i)];
        
        for (size_t j = 0; j < length; ++j) {
            // This formula maps Node Index -> Component Index
            // Example: Node 5, Length 3 -> Indices 15, 16, 17
            N[static_cast<size_t>(i) * length + j] = (original_val * length) + j;
        }
    }
}


template class MathOperatorRHS<1U>;
template class MathOperatorRHS<2U>;
template class MathOperatorRHS<3U>;

template class MathOperatorRHS<1U,Face>;
template class MathOperatorRHS<2U,Face>;
template class MathOperatorRHS<3U,Face>;

// TODO: this instantiation should not be necessary
template class MathOperatorRHS<1U,InterFace>;
template class MathOperatorRHS<2U,InterFace>;
template class MathOperatorRHS<3U,InterFace>;


} // end namespace csmp

