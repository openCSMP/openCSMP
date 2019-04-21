#include "MathOperatorLHS.h"
#include "MathOperatorRHS.h"
#include "Element.h"
#include "InterFace.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

  template<size_t dim>
  MathOperatorLHS<dim>::MathOperatorLHS()
    : name_("unspecified LHS-operator"),
    IDT(3),
    IDB(3),
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



  template<size_t dim>
  MathOperatorLHS<dim>::MathOperatorLHS(const PropertyDatabase<dim>& pref,
    const char* basic,
    const char* test)
    : name_("unspecified LHS-operator"),
    bop(make_pair(pref.Parameter(basic), 0)),
    top(make_pair(pref.Parameter(test), 0)),
    IDT(3),
    IDB(3),
    MTRL(13, DenseMatrix<DM_MIN>(2, 2)),
    DERIV(2, 3),
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
  MathOperatorLHS<dim>::MathOperatorLHS(const PropertyDatabase<dim>& pref,
    const char* oper,
    const char* basic,
    const char* test)
    : name_("unspecified LHS-operator"),
    op(pref.Parameter(oper)),
    bop(make_pair(pref.Parameter(basic), 0)),
    top(make_pair(pref.Parameter(test), 0)),
    IDT(3),
    IDB(3),
    MTRL(13, DenseMatrix<DM_MIN>(2, 2)),
    DERIV(2, 3),
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
  MathOperatorLHS<dim>::MathOperatorLHS(const MathOperatorLHS<dim>& mo)
    : name_(mo.name_),
    op(mo.op),        // operand
    bop(mo.bop),       // basic operand
    top(mo.top),       // testfunction operand
    IDT(mo.IDT),       // node-ID & global constraint points vector
    IDB(mo.IDB),       // node-ID & global constraint points vector
    LHS(mo.LHS),
    MTRL(mo.MTRL),      // Operand storage
    DERIV(mo.DERIV),
    IPOL(mo.IPOL),
    SC(mo.SC),
    VC(mo.VC),
    TS(mo.TS),
    AR(mo.AR),
    FR(mo.FR),
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




  template<size_t dim>
  MathOperatorLHS<dim>& MathOperatorLHS<dim>::operator=(const MathOperatorLHS<dim>& mo)
  {
    if (this != &mo)
    {
      name_ = mo.name_;
      op = mo.op;    // operand
      bop = mo.bop;   // basic operand
      top = mo.top;   // testfunction operand
      IDT = mo.IDT;   // node-ID & global constraint points vector
      IDB = mo.IDB;   // node-ID & global constraint points vector
      LHS = mo.LHS;
      MTRL = mo.MTRL;  // Operand storage
      DERIV = mo.DERIV;
      IPOL = mo.IPOL;
      SC = mo.SC;
      VC = mo.VC;
      TS = mo.TS;
      AR = mo.AR;
      FR = mo.FR;
      factor_ = mo.factor_;
      add_accumulate_ = mo.add_accumulate_;
      add_accumulate_later_ = mo.add_accumulate_later_;
      multiply_accumulate_ = mo.multiply_accumulate_;
      lump_matrices_ = mo.lump_matrices_;
      application_cycles_ = mo.application_cycles_;
      application_cycle_ = mo.application_cycle_;
      time_multiply_ = mo.time_multiply_;
      time_divide_ = mo.time_divide_;
    }
    return *this;
  }


  template<size_t dim>
  MathOperatorLHS<dim>::~MathOperatorLHS()
  {
  }


  /// ===================================  Operand Functions ========================================================================

  template<size_t dim>
  void MathOperatorLHS<dim>::Name(const char* s,
    const char* bopname,
    const char* topname)
  {
    name_ = s;
    name_ += ": Operand: '";
    name_ += bopname;
    name_ += "', '";
    name_ += topname;
  }


  template<size_t dim>
  void MathOperatorLHS<dim>::Name(const char* s,
    const char* opname,
    const char* bopname,
    const char* topname)
  {
    name_ = s;
    name_ += ": Operand: '";
    name_ += opname;
    name_ += "', '";
    name_ += bopname;
    name_ += "', '";
    name_ += topname;
  }

  template<size_t dim>
  std::string  MathOperatorLHS<dim>::Name() const
  {
    return name_;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::Out()  const
  {
    op.Out();
    bop.first.Out();
    top.first.Out();
  }

  // MATERIAL PROPERTY OPERAND

  template<size_t dim>
  const csmp::Parameter&   MathOperatorLHS<dim>::MaterialOperand() const
  {
    return op;
  }

  template<size_t dim>
  const csmp::Index&   MathOperatorLHS<dim>::MaterialOperandKey() const
  {
    return op.key;
  }

  template<size_t dim>
  std::string  MathOperatorLHS<dim>::MaterialOperandName() const
  {
    return op.name;
  }

  template<size_t dim>
  csmp::VARIABLE_TYPE   MathOperatorLHS<dim>::MaterialOperandType() const
  {
    return op.key.type;
  }

  template<size_t dim>
  csmp::PLACEMENT   MathOperatorLHS<dim>::MaterialOperandPlacement() const
  {
    return op.key.place;
  }

  template<size_t dim>
  size_t   MathOperatorLHS<dim>::MaterialOperandDataDepth() const
  {
    return op.key.dataDepth;
  }


  // BASIC FUNCTION OPERAND

  template<size_t dim>
  const csmp::Parameter&   MathOperatorLHS<dim>::BasicOperand() const
  {
    return bop.first;
  }

  template<size_t dim>
  const csmp::Index&   MathOperatorLHS<dim>::BasicOperandKey() const
  {
    return bop.first.key;
  }

  template<size_t dim>
  std::string  MathOperatorLHS<dim>::BasicOperandName() const
  {
    return bop.first.name;
  }

  template<size_t dim>
  csmp::VARIABLE_TYPE   MathOperatorLHS<dim>::BasicOperandType() const
  {
    return bop.first.key.type;
  }

  template<size_t dim>
  csmp::PLACEMENT   MathOperatorLHS<dim>::BasicOperandPlacement() const
  {
    return bop.first.key.place;
  }

  template<size_t dim>
  size_t   MathOperatorLHS<dim>::BasicOperandDataDepth() const
  {
    return bop.first.key.dataDepth;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::BasicOperandOffset(size_t os)
  {
    bop.second = os;
  }

  template<size_t dim>
  size_t MathOperatorLHS<dim>::BasicOperandOffset()
    const
  {
    return bop.second;
  }


  // TEST FUNCTION OPERAND

  template<size_t dim>
  const csmp::Parameter&  MathOperatorLHS<dim>::TestOperand() const
  {
    return top.first;
  }

  template<size_t dim>
  const csmp::Index&  MathOperatorLHS<dim>::TestOperandKey() const
  {
    return top.first.key;
  }

  template<size_t dim>
  std::string  MathOperatorLHS<dim>::TestOperandName() const
  {
    return top.first.name;
  }

  template<size_t dim>
  csmp::VARIABLE_TYPE   MathOperatorLHS<dim>::TestOperandType() const
  {
    return top.first.key.type;
  }

  template<size_t dim>
  csmp::PLACEMENT   MathOperatorLHS<dim>::TestOperandPlacement() const
  {
    return top.first.key.place;
  }

  template<size_t dim>
  size_t   MathOperatorLHS<dim>::TestOperandDataDepth() const
  {
    return top.first.key.dataDepth;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::TestOperandOffset(size_t os)
  {
    top.second = os;
  }

  template<size_t dim>
  size_t MathOperatorLHS<dim>::TestOperandOffset() const
  {
    return top.second;
  }





  /// ===================================  Accumulation Process Settings ========================================================================


  template<size_t dim>
  void MathOperatorLHS<dim>::AddAccumulate()
  {
    add_accumulate_ = true;
    subtract_accumulate_ = false;
    add_accumulate_later_ = false;
    subtract_accumulate_later_ = false;
    multiply_accumulate_ = false;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::SubtractAccumulate()
  {
    add_accumulate_ = false;
    subtract_accumulate_ = true;
    add_accumulate_later_ = false;
    subtract_accumulate_later_ = false;
    multiply_accumulate_ = false;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::AddAccumulateLater()
  {
    add_accumulate_ = false;
    subtract_accumulate_ = false;
    add_accumulate_later_ = true;
    subtract_accumulate_later_ = false;
    multiply_accumulate_ = false;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::SubtractAccumulateLater()
  {
    add_accumulate_ = false;
    subtract_accumulate_ = false;
    add_accumulate_later_ = false;
    subtract_accumulate_later_ = true;
    multiply_accumulate_ = false;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::MultiplyAccumulate()
  {
    add_accumulate_ = false;
    subtract_accumulate_ = false;
    add_accumulate_later_ = false;
    subtract_accumulate_later_ = false;
    multiply_accumulate_ = true;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::LumpedFormulation(bool lumped)
  {
    lump_matrices_ = lumped;
  }

  template<size_t dim>
  void  MathOperatorLHS<dim>::ApplicationCycles(size_t c)
  {
    application_cycles_ = c;
  }

  template<size_t dim>
  void  MathOperatorLHS<dim>::ApplicationCycle(size_t c)
  {
    application_cycle_ = c;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::MultiplyBy(double64 multiplication_factor)
  {
    factor_ = multiplication_factor;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::MultiplyWithTimeIncrement(bool multiply)
  {
    time_multiply_ = multiply;
  }

  template<size_t dim>
  void MathOperatorLHS<dim>::DivideByTimeIncrement(bool divide)
  {
    time_divide_ = divide;
  }



  template<size_t dim>
  bool MathOperatorLHS<dim>::Add() const
  {
    return add_accumulate_;
  }

  template<size_t dim>
  bool MathOperatorLHS<dim>::Subtract() const
  {
    return subtract_accumulate_;
  }

  template<size_t dim>
  bool MathOperatorLHS<dim>::AddLater() const
  {
    return add_accumulate_later_;
  }

  template<size_t dim>
  bool MathOperatorLHS<dim>::SubtractLater() const
  {
    return subtract_accumulate_later_;
  }

  template<size_t dim>
  bool MathOperatorLHS<dim>::Multiply() const
  {
    return multiply_accumulate_;
  }

  template<size_t dim>
  bool  MathOperatorLHS<dim>::LumpedFormulation() const
  {
    return lump_matrices_;
  }

  template<size_t dim>
  size_t MathOperatorLHS<dim>::ApplicationCycles() const
  {
    return application_cycles_;
  }

  template<size_t dim>
  size_t  MathOperatorLHS<dim>::ApplicationCycle() const
  {
    return application_cycle_;
  }

  template<size_t dim>
  double64   MathOperatorLHS<dim>::MultiplyBy() const
  {
    return factor_;
  }

  template<size_t dim>
  bool MathOperatorLHS<dim>::MultiplyWithTimeIncrement() const
  {
    return time_multiply_;
  }

  template<size_t dim>
  bool MathOperatorLHS<dim>::DivideByTimeIncrement() const
  {
    return time_divide_;
  }


  /// ===================================  Accumulation Process Settings ========================================================================

  /// interpolation of property if isoparametric elements are used
  template<size_t dim>
  void MathOperatorLHS<dim>::PropertyAtIntegrationPoint(const Element<dim>& e_ref,
    const csmp::Index&  idx,
    size_t ip, DenseMatrix<DM_MIN>& M)
  {
    if (!e_ref.UsesLocalCoordinates())
      throw csmp::Exception(FATAL_ERROR,
        "MathOperatorLHS<dim>::PropertyValueAtIntegrationPoint",
        "Current element does not support numerical integration.");

    // 0. If the property is a constraint point variable a one-to-one mapping
    //    can be performed; no interpolation is needed
    if (idx.place == ELEMENT_INTEGRATION_POINT) {
      if (idx.type == SCALAR) {
        M.AssignToDiagonal(dim, e_ref.Read(ip, idx));
      }
      else if (idx.type == VECTOR) {
        e_ref.Read(ip, idx, VC[0]);
        M.AssignToDiagonal(VC[0]);
      }
      else if (idx.type == TENSOR) {
        e_ref.Read(ip, idx, TS[0]);
        M = TS[0];
      }
      else if (idx.type == ARRAY)
      {
        AR[0].Resize(idx.dataDepth);
        e_ref.Read(ip, idx, AR[0]);
        M.AssignToDiagonal(AR[0]);
      }
      else if (idx.type == FLAGGEDARRAY)
      {
        FR[0].Resize(idx.dataDepth);
        e_ref.Read(ip, idx, FR[0]);
        M.AssignToDiagonal(FR[0]);
      }
      return;
    }

    // 1. If the property is an element variable it is constant over the
    //    element, such that no interpolation is needed
    if (idx.place == ELEMENT || idx.place == REGION)
    {
      if (idx.type == SCALAR) {
        M.AssignToDiagonal(dim, e_ref.Read(idx));
      }
      else if (idx.type == VECTOR) {
        e_ref.Read(idx, VC[0]);
        M.AssignToDiagonal(VC[0]);
      }
      else if (idx.type == TENSOR) {
        e_ref.Read(idx, TS[0]);
        M = TS[0];
      }
      else if (idx.type == ARRAY)
      {
        AR[0].Resize(idx.dataDepth);
        e_ref.Read(idx, AR[0]);
        M.AssignToDiagonal(AR[0]);
      }
      else if (idx.type == FLAGGEDARRAY)
      {
        FR[0].Resize(idx.dataDepth);
        e_ref.Read(idx, FR[0]);
        M.AssignToDiagonal(FR[0]);
      }
      return;
    }

    // 2. If the property is a node property, it must be interpolated to the
    //    integration point
    if (idx.place == NODE)
    {
      // interpolating properties
      if (idx.type == SCALAR) {
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, SC[0]);
        M.AssignToDiagonal(dim, SC[0]);
      }
      else if (idx.type == VECTOR) {
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, VC[0]);
        M.AssignToDiagonal(VC[0]);
      }
      else if (idx.type == TENSOR) {
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, TS[0]);
        M = TS[0];
      }
      else if (idx.type == ARRAY)
      {
        AR[0].Resize(idx.dataDepth);
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, AR[0]);
        M.AssignToDiagonal(AR[0]);
      }
      else if (idx.type == FLAGGEDARRAY)
      {
        FR[0].Resize(idx.dataDepth);
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, FR[0]);
        M.AssignToDiagonal(FR[0]);
      }
      return;
    }

    throw csmp::Exception(FATAL_ERROR,
      "MathOperatorLHS<dim>::PropertyValueAtIntegrationPoint",
      "This method cannot interpolate IntegrationPoint variables.");

  } // end PropertyValueAtIntegrationPoint( Element )






  template<size_t dim>
  void MathOperatorLHS<dim>::PropertyAtIntegrationPoint(const Face<dim>& e_ref,
    const csmp::Index&  idx,
    size_t ip, DenseMatrix<DM_MIN>& M)
  {
    if (!e_ref.UsesLocalCoordinates())
      throw csmp::Exception(FATAL_ERROR,
        "MathOperatorLHS<dim>::PropertyValueAtIntegrationPoint",
        "Current element does not support numerical integration.");

    // 0. If the property is a constraint point variable a one-to-one mapping
    //    can be performed; no interpolation is needed
    if (idx.place == FACE_INTEGRATION_POINT) {
      if (idx.type == SCALAR) {
        M.AssignToDiagonal(dim, e_ref.Read(ip, idx));
      }
      else if (idx.type == VECTOR) {
        e_ref.Read(ip, idx, VC[0]);
        M.AssignToDiagonal(VC[0]);
      }
      else if (idx.type == TENSOR) {
        e_ref.Read(ip, idx, TS[0]);
        M = TS[0];
      }
      else if (idx.type == ARRAY)
      {
        AR[0].Resize(idx.dataDepth);
        e_ref.Read(ip, idx, AR[0]);
        M.AssignToDiagonal(AR[0]);
      }
      else if (idx.type == FLAGGEDARRAY)
      {
        FR[0].Resize(idx.dataDepth);
        e_ref.Read(ip, idx, FR[0]);
        M.AssignToDiagonal(FR[0]);
      }
      return;
    }

    // 1. If the property is an element variable it is constant over the
    //    element, such that no interpolation is needed
    if (idx.place == FACE || idx.place == BOUNDARY)
    {
      if (idx.type == SCALAR) {
        M.AssignToDiagonal(dim, e_ref.Read(idx));
      }
      else if (idx.type == VECTOR) {
        e_ref.Read(idx, VC[0]);
        M.AssignToDiagonal(VC[0]);
      }
      else if (idx.type == TENSOR) {
        e_ref.Read(idx, TS[0]);
        M = TS[0];
      }
      else if (idx.type == ARRAY)
      {
        AR[0].Resize(idx.dataDepth);
        e_ref.Read(idx, AR[0]);
        M.AssignToDiagonal(AR[0]);
      }
      else if (idx.type == FLAGGEDARRAY)
      {
        FR[0].Resize(idx.dataDepth);
        e_ref.Read(idx, FR[0]);
        M.AssignToDiagonal(FR[0]);
      }
      return;
    }

    // 2. If the property is a node property, it must be interpolated to the
    //    integration point
    if (idx.place == NODE)
    {
      // interpolating properties
      if (idx.type == SCALAR) {
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, SC[0]);
        M.AssignToDiagonal(dim, SC[0]);
      }
      else if (idx.type == VECTOR) {
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, VC[0]);
        M.AssignToDiagonal(VC[0]);
      }
      else if (idx.type == TENSOR) {
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, TS[0]);
        M = TS[0];
      }
      else if (idx.type == ARRAY)
      {
        AR[0].Resize(idx.dataDepth);
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, AR[0]);
        M.AssignToDiagonal(AR[0]);
      }
      else if (idx.type == FLAGGEDARRAY)
      {
        FR[0].Resize(idx.dataDepth);
        e_ref.PropertyValueAtIntegrationPoint(idx, ip, FR[0]);
        M.AssignToDiagonal(FR[0]);
      }
      return;
    }

    throw csmp::Exception(FATAL_ERROR,
      "MathOperatorLHS<dim>::PropertyValueAtIntegrationPoint",
      "This method cannot interpolate IntegrationPoint variables.");

  } // end PropertyValueAtIntegrationPoint( Face )


    /**

    GetOperands fills a vector of material property matrices with the
    required values for later computation. These matrices always have
    the dimensions spatial-dimension^2 and they will hold either scalar,
    vector or tensor properties, depending on what kind of property the
    Operand is.

    When the property is an element property, it will be put into the
    first vector entry MTRL[0].
    */
  template<size_t dim>
  void MathOperatorLHS<dim>::GetOperands(Element<dim>& e_ref)
  {
    // if operand property is an element property
    if (MaterialOperandPlacement() == ELEMENT || MaterialOperandPlacement() == REGION)
    {
      MTRL.resize(1U);
      if (MaterialOperandType() == SCALAR)
        MTRL[0].AssignToDiagonal(dim, e_ref.Read(MaterialOperandKey()));
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
      if (MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT) {
        MTRL.resize(e_ref.FE()->IntegrationPoints());
        const size_t n_integration_points(e_ref.FE()->IntegrationPoints());
        for (size_t i = 0U; i<n_integration_points; i++)
        {
          if (MaterialOperandType() == SCALAR)
            MTRL[i].AssignToDiagonal(dim, e_ref.Read(i, MaterialOperandKey()));
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
        const size_t n_integration_points(e_ref.FE()->IntegrationPoints());
        for (size_t i = 0U; i<n_integration_points; i++)
          PropertyAtIntegrationPoint(e_ref, MaterialOperandKey(), i, MTRL[i]);
      }
      else
        throw csmp::Exception(FATAL_ERROR,
          "MathOperatorLHS<dim>::GetOperands",
          "Face based operands cannot be accumulated with this method");
    }

  } // end GetOperands(Element)

  template<size_t dim>
  void MathOperatorLHS<dim>::GetOperands(Face<dim>&  e_ref)
  {
    // if operand property is an element property
    if (MaterialOperandPlacement() == FACE || MaterialOperandPlacement() == BOUNDARY)
    {
      MTRL.resize(1U);
      if (MaterialOperandType() == SCALAR)
        MTRL[0].AssignToDiagonal(dim, e_ref.Read(MaterialOperandKey()));
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
      if (MaterialOperandPlacement() == FACE_INTEGRATION_POINT) {
        MTRL.resize(e_ref.FE()->IntegrationPoints());
        const size_t n_integration_points(e_ref.FE()->IntegrationPoints());
        for (size_t i = 0U; i<n_integration_points; i++)
        {
          if (MaterialOperandType() == SCALAR)
            MTRL[i].AssignToDiagonal(dim, e_ref.Read(i, MaterialOperandKey()));
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
        const size_t n_integration_points(e_ref.FE()->IntegrationPoints());
        for (size_t i = 0U; i<n_integration_points; i++)
          PropertyAtIntegrationPoint(e_ref, MaterialOperandKey(), i, MTRL[i]);
      }
      else
        throw csmp::Exception(FATAL_ERROR,
          "MathOperatorLHS<dim>::GetOperands",
          "Element based operands cannot be accumulated with this method");
    }

  } // end GetOperands(Face)

  template<size_t dim>
  void MathOperatorLHS<dim>::GetOperands(InterFace<dim>& f)
  {
    std::cerr << "\nMathOperatorLHS<dim>::GetOperands: ";
    std::cerr << " Overload to get LHS operands from interface: " << f.Idx() << std::endl;
    throw invalid_argument("MathOperatorLHS<dim>::GetOperands(InterFace)");
  } // end GetOperands(InterFace)


  template<size_t dim>
  void MathOperatorLHS<dim>::WriteOperands(Element<dim>& e)
  {
    std::cerr << "\nMathOperatorLHS<dim>::WriteOperands: ";
    std::cerr << " Overload to write LHS operands from element: " << e.Idx() << std::endl;
    throw invalid_argument("MathOperatorLHS<dim>::WriteOperands(Element)");
  } // end WriteOperands(Element)


  template<size_t dim>
  void MathOperatorLHS<dim>::WriteOperands(Face<dim>& f)
  {
    std::cerr << "\nMathOperatorLHS<dim>::WriteOperands: ";
    std::cerr << " Overload to write LHS operands from face: " << f.Idx() << std::endl;
    throw invalid_argument("MathOperatorLHS<dim>::WriteOperands(Face)");
  } // end WriteOperands(Face)

  template<size_t dim>
  void MathOperatorLHS<dim>::WriteOperands(InterFace<dim>& f)
  {
    std::cerr << "\nMathOperatorLHS<dim>::WriteOperands: ";
    std::cerr << " Overload to write LHS operands from interface: " << f.Idx() << std::endl;
    throw invalid_argument("MathOperatorLHS<dim>::WriteOperands(InterFace)");
  } // end WriteOperands(InterFace)

  template<size_t dim>
  void MathOperatorLHS<dim>::ComputeContribution(Element<dim>& e)
  {
    std::cerr << "\nMathOperatorLHS<dim>::ComputeContribution: ";
    std::cerr << " Overload to calculate LHS contribution from Element: " << e.Idx() << std::endl;
    throw invalid_argument("MathOperatorLHS<dim>::ComputeContribution(Element)");
  }// end ComputeContribution(Element)

  template<size_t dim>
  void MathOperatorLHS<dim>::ComputeContribution(Face<dim>& f)
  {
    std::cerr << "\nMathOperatorLHS<dim>::ComputeContribution: ";
    std::cerr << " Overload to calculate LHS constribution from face: " << f.Idx() << std::endl;
    throw invalid_argument("MathOperatorLHS<dim>::ComputeContribution(Face)");
  }// end ComputeContribution(Face)

  template<size_t dim>
  void MathOperatorLHS<dim>::ComputeContribution(InterFace<dim>& f)
  {
    std::cerr << "\nMathOperatorLHS<dim>::ComputeContribution: ";
    std::cerr << " Overload to calculate LHS constribution from interface: " << f.Idx() << std::endl;
    throw invalid_argument("MathOperatorLHS<dim>::ComputeContribution(InterFace)");
  }// end ComputeContribution(InterFace)


  template<size_t dim>
  inline void MathOperatorLHS<dim>::MultiplyWithTimeFactor(double64 dt)
  {
    LHS *= dt;
  }

  /// AssignToGlobal matrix function
  template<size_t dim>
  void MathOperatorLHS<dim>::AssignToGlobal(const Element<dim>& e, SparseMatrix& G, std::vector<double64>& pivotVector, const  std::vector<size_t>& DOF_indexes)
  {

    // map local to global indexes for test and basic operands
    
    IDT.resize(e.Nodes()); // ii
    IDB.resize(e.Nodes()); // jj
    for (size_t i = 0U; i<e.Nodes(); i++) {
      IDT[i] = e.N(i)->Idx();       
      IDB[i] = IDT[i];
    }

    if (TestOperandType() != SCALAR)
      transformNodeIndexVector(dim, TestOperandKey(), IDT);

    for (size_t i = 0U; i < IDT.size(); i++) {
      IDT[i] += this->TestOperandOffset();
      IDT[i] = DOF_indexes[IDT[i]];
    }


    if (BasicOperandType() != SCALAR)
      transformNodeIndexVector(dim, BasicOperandKey(), IDB);

    for (size_t i = 0U; i < IDB.size(); i++) {
      IDB[i] += this->BasicOperandOffset();
      IDB[i] = DOF_indexes[IDB[i]];
    }
    // perform assignment from local matrix to global matrix

    if (multiply_accumulate_)
    {
		throw csmp::Exception(ERROR,
			"MathOperatorLHS<dim>::AssignToGlobal(Element):",
			"Multiply Accumulate has not been supported yet");

        //if (IDT[i] != NULL_IDX) {
        //  for (size_t j = 0U; j<LHS.Cols(); j++)
        //    if (IDB[j] != NULL_IDX) {
        //      G.MultiplyEntryWith(IDT[i],
        //        IDB[j],
        //        LHS(i, j) * factor_);
        //    }
        //}
    }
    else if (add_accumulate_ || add_accumulate_later_)
    {
		for (size_t i = 0U; i < LHS.Rows(); i++) {
			if (IDT[i] != NULL_IDX) {			
				for (size_t j = 0U; j < LHS.Cols(); j++) {
					if (IDB[j] == NULL_IDX) {
						pivotVector[IDT[i]] -= LHS(i, j) * e.N(j)->Read(TestOperandKey());  // notice the sign
					}
					else {
						G.Add(IDT[i],
							IDB[j],
							LHS(i, j) * factor_);
					}					
				}			
			}
		}

      /*for (size_t i = 0U; i<LHS.Rows(); i++)
        if (IDT[i] != NULL_IDX) {
          for (size_t j = 0U; j < LHS.Cols(); j++)
            if (IDB[j] != NULL_IDX) {
              G.Add(IDT[i],
                IDB[j],
                LHS(i, j) * factor_);
            }
        }*/
    }
    else if (subtract_accumulate_ || subtract_accumulate_later_)
    {
		for (size_t i = 0U; i < LHS.Rows(); i++) {
			if (IDT[i] != NULL_IDX) {
				for (size_t j = 0U; j < LHS.Cols(); j++) {
					if (IDB[j] == NULL_IDX) {
						pivotVector[IDT[i]] += LHS(i, j) * e.N(j)->Read(TestOperandKey());  // notice the sign
					}
					else {
						G.Add(IDT[i],
							IDB[j],
							-LHS(i, j) * factor_);
					}
				}
			}
		}

      /*for (size_t i = 0U; i<LHS.Rows(); i++)
        if (IDT[i] != NULL_IDX) {
          for (size_t j = 0U; j < LHS.Cols(); j++)
            if (IDB[j] != NULL_IDX) {
              G.Add(IDT[i],
                IDB[j],
                -LHS(i, j) * factor_);
            }
        }*/
    }
    else
      throw csmp::Exception(ERROR,
        "MathOperatorLHS<dim>::AssignToGlobal(Element):",
        "accumulation instructions could not be parsed.");

  } // end AssignToGlobal (Element)

  template<size_t dim>
  void MathOperatorLHS<dim>::AssignToGlobal(const Face<dim>& e, SparseMatrix& G, std::vector<double64>& pivotVector, const  std::vector<size_t>& DOF_indexes)
  {
    // map local to global indexes for test and basic operands
    IDT.resize(e.Nodes());
    IDB.resize(e.Nodes());
    for (size_t i = 0U; i<e.Nodes(); i++) {
      IDT[i] = e.N(i)->Idx();
      IDB[i] = IDT[i];
    }

    if (TestOperandType() != SCALAR)
      transformNodeIndexVector(dim, TestOperandKey(), IDT);

    for (size_t i = 0U; i < IDT.size(); i++) {
      IDT[i] += this->TestOperandOffset();
      IDT[i] = DOF_indexes[IDT[i]];
    }


    if (BasicOperandType() != SCALAR)
      transformNodeIndexVector(dim, BasicOperandKey(), IDB);

    for (size_t i = 0U; i < IDB.size(); i++) {
      IDB[i] += this->BasicOperandOffset();
      IDB[i] = DOF_indexes[IDB[i]];
    }


    // perform assignment from local matrix to global matrix

    if (multiply_accumulate_)
    {
		throw csmp::Exception(ERROR,
			"MathOperatorLHS<dim>::AssignToGlobal(face):",
			"Multiply Accumulate has not been supported yet");
      /*for (size_t i = 0U; i<LHS.Rows(); i++)
        if (IDT[i] != NULL_IDX) {
          for (size_t j = 0U; j<LHS.Cols(); j++)
            if (IDB[j] != NULL_IDX) {
              G.MultiplyEntryWith(IDT[i],
                IDB[j],
                LHS(i, j) * factor_);
            }
        }*/
    }
    else if (add_accumulate_ || add_accumulate_later_)
    {
		for (size_t i = 0U; i < LHS.Rows(); i++) {
			if (IDT[i] != NULL_IDX) {
				for (size_t j = 0U; j < LHS.Cols(); j++) {
					if (IDB[j] == NULL_IDX) {
						pivotVector[IDT[i]] -= LHS(i, j) * e.N(j)->Read(TestOperandKey());  // notice the sign
					}
					else {
						G.Add(IDT[i],
							IDB[j],
							LHS(i, j) * factor_);
					}
				}
			}
		}

      /*for (size_t i = 0U; i<LHS.Rows(); i++)
        if (IDT[i] != NULL_IDX) {
          for (size_t j = 0U; j < LHS.Cols(); j++)
            if (IDB[j] != NULL_IDX) {
              G.Add(IDT[i],
                IDB[j],
                LHS(i, j) * factor_);
            }
        }*/
    }
    else if (subtract_accumulate_ || subtract_accumulate_later_)
    {
		for (size_t i = 0U; i < LHS.Rows(); i++) {
			if (IDT[i] != NULL_IDX) {
				for (size_t j = 0U; j < LHS.Cols(); j++) {
					if (IDB[j] == NULL_IDX) {
						pivotVector[IDT[i]] += LHS(i, j) * e.N(j)->Read(TestOperandKey());  // notice the sign
					}
					else {
						G.Add(IDT[i],
							IDB[j],
							-LHS(i, j) * factor_);
					}
				}
			}
		}
      /*for (size_t i = 0U; i<LHS.Rows(); i++)
        if (IDT[i] != NULL_IDX) {
          for (size_t j = 0U; j < LHS.Cols(); j++)
            if (IDB[j] != NULL_IDX) {
              G.Add(IDT[i],
                IDB[j],
                -LHS(i, j) * factor_);
            }
        }*/
    }
    else
      throw csmp::Exception(ERROR,
        "MathOperatorLHS<dim>::AssignToGlobal(Element):",
        "accumulation instructions could not be parsed.");
  
  } // end AssignToGlobal (Face)

  template<size_t dim>
  void MathOperatorLHS<dim>::AssignToGlobal(const InterFace<dim>& f, SparseMatrix& G, std::vector<double64>& pivotVector, const  std::vector<size_t>& DOF_indexes)
  {
    std::cerr << "\nMathOperatorLHS<dim>::AssignToGlobal(InterFace): ";
    std::cerr << " Overload to assign LHS local entries to global matrix: " << f.Idx() << std::endl;
    throw invalid_argument("MathOperatorLHS<dim>::AssignToGlobal(InterFace)");
  } // end AssignToGlobal (InterFace)


  template<size_t dim>
  void MathOperatorLHS<dim>::AssignToGlobal(const Element<dim>& e, SparseMatrix& G)
  {
    // map local to global indexes for test and basic operands

    IDT.resize(e.Nodes());
    IDB.resize(e.Nodes());
    for (size_t i = 0U; i<e.Nodes(); i++) {
      IDT[i] = e.N(i)->Idx();
      IDB[i] = IDT[i];
    }

    if (TestOperandType() != SCALAR)
      transformNodeIndexVector(dim, TestOperandKey(), IDT);

    for (size_t i = 0U; i<IDT.size(); i++)
      IDT[i] += this->TestOperandOffset();

    if (BasicOperandType() != SCALAR)
      transformNodeIndexVector(dim, BasicOperandKey(), IDB);

    for (size_t i = 0U; i<IDB.size(); i++)
      IDB[i] += this->BasicOperandOffset();

    // perform assignment from local matrix to global matrix

    if (multiply_accumulate_)
    {
      for (size_t i = 0U; i<LHS.Rows(); i++)
        for (size_t j = 0U; j<LHS.Cols(); j++)
          G.MultiplyEntryWith(IDT[i],
            IDB[j],
            LHS(i, j) * factor_);
    }
    else if (add_accumulate_ || add_accumulate_later_)
    {
      for (size_t i = 0U; i<LHS.Rows(); i++)
        for (size_t j = 0U; j<LHS.Cols(); j++)
          G.Add(IDT[i],
            IDB[j],
            LHS(i, j) * factor_);
    }
    else if (subtract_accumulate_ || subtract_accumulate_later_)
    {
      for (size_t i = 0U; i<LHS.Rows(); i++)
        for (size_t j = 0U; j<LHS.Cols(); j++)
          G.Add(IDT[i],
            IDB[j],
            -LHS(i, j) * factor_);
    }
    else
      throw csmp::Exception(ERROR,
        "MathOperatorLHS<dim>::AssignToGlobal(Element):",
        "accumulation instructions could not be parsed.");

  } // end AssignToGlobal (Element)

  template<size_t dim>
  void MathOperatorLHS<dim>::AssignToGlobal(const Face<dim>& e, SparseMatrix& G)
  {
    // map local to global indexes for test and basic operands

    IDT.resize(e.Nodes());
    IDB.resize(e.Nodes());
    for (size_t i = 0U; i<e.Nodes(); i++) {
      IDT[i] = e.N(i)->Idx();
      IDB[i] = IDT[i];
    }

    if (TestOperandType() != SCALAR)
      transformNodeIndexVector(dim, TestOperandKey(), IDT);

    for (size_t i = 0U; i<IDT.size(); i++)
      IDT[i] += this->TestOperandOffset();

    if (BasicOperandType() != SCALAR)
      transformNodeIndexVector(dim, BasicOperandKey(), IDB);

    for (size_t i = 0U; i<IDB.size(); i++)
      IDB[i] += this->BasicOperandOffset();

    // perform assignment from local matrix to global matrix

    if (multiply_accumulate_)
    {
      for (size_t i = 0U; i<LHS.Rows(); i++)
        for (size_t j = 0U; j<LHS.Cols(); j++)
          G.MultiplyEntryWith(IDT[i],
            IDB[j],
            LHS(i, j) * factor_);
    }
    else if (add_accumulate_ || add_accumulate_later_)
    {
      for (size_t i = 0U; i<LHS.Rows(); i++)
        for (size_t j = 0U; j<LHS.Cols(); j++)
          G.Add(IDT[i],
            IDB[j],
            LHS(i, j) * factor_);
    }
    else if (subtract_accumulate_ || subtract_accumulate_later_)
    {
      for (size_t i = 0U; i<LHS.Rows(); i++)
        for (size_t j = 0U; j<LHS.Cols(); j++)
          G.Add(IDT[i],
            IDB[j],
            -LHS(i, j) * factor_);
    }
    else
      throw csmp::Exception(ERROR,
        "MathOperatorLHS<dim>::AssignToGlobal(Face):",
        "accumulation instructions could not be parsed.");

  } // end AssignToGlobal (Face)

  template<size_t dim>
  void MathOperatorLHS<dim>::AssignToGlobal(const InterFace<dim>& f, SparseMatrix& G)
  {
    std::cerr << "\nMathOperatorLHS<dim>::AssignToGlobal(InterFace): ";
    std::cerr << " Overload to assign LHS local entries to global matrix: " << f.Idx() << std::endl;
    throw invalid_argument("MathOperatorLHS<dim>::AssignToGlobal(InterFace)");
  } // end AssignToGlobal (InterFace)


  template class MathOperatorLHS<1U>;
  template class MathOperatorLHS<2U>;
  template class MathOperatorLHS<3U>;

} // end namespace csmp
