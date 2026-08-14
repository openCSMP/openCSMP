#include "MathOperatorLHS.h"
#include "MathOperatorRHS.h"
#include "Element.h"
#include "InterFace.h"
#include "Face.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
MathOperatorLHS<dim,CELL>::MathOperatorLHS()
  : name_("unspecified LHS-operator"),
    MTRL(13, DenseMatrix<dim>(dim,dim)),
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
MathOperatorLHS<dim,CELL>::MathOperatorLHS( const PropertyDatabase<dim>& pref,
                                            const char* basic,
                                            const char* test)
		: name_("unspecified LHS-operator"),
      bop(make_pair(pref.Parameter(basic), 0)),
      top(make_pair(pref.Parameter(test), 0)),
      MTRL(13, DenseMatrix<dim>(dim,dim)),
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
MathOperatorLHS<dim,CELL>::MathOperatorLHS( const PropertyDatabase<dim>& pref,
                                            const char* oper,
                                            const char* basic,
                                            const char* test )
		: name_("unspecified LHS-operator"),
      op(pref.Parameter(oper)),
      bop(make_pair(pref.Parameter(basic), 0)),
      top(make_pair(pref.Parameter(test), 0)),
      MTRL(13, DenseMatrix<dim>(dim,dim)),
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
MathOperatorLHS<dim,CELL>::MathOperatorLHS( const MathOperatorLHS<dim,CELL>& mo )
		: name_(mo.name_),
      op(mo.op),        // operand
      bop(mo.bop),       // basic operand
      top(mo.top),       // testfunction operand
      LHS(mo.LHS),
      MTRL(mo.MTRL),      // Operand storage
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
MathOperatorLHS<dim,CELL>& MathOperatorLHS<dim,CELL>::operator=( const MathOperatorLHS<dim,CELL>& mo )
	{
		if (this != &mo)
      {
        name_ = mo.name_;
        op = mo.op;    // operand
        bop = mo.bop;   // basic operand
        top = mo.top;   // testfunction operand
        LHS = mo.LHS;
        MTRL = mo.MTRL;  // Operand storage
        DERIV = mo.DERIV;
        IPOL = mo.IPOL;
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




	/// ===================================  Operand Functions ========================================================================

template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::Name( const char* s,
                                      const char* bopname,
                                      const char* topname) noexcept
	{
		name_ = s;
		name_ += ": Operand: '";
		name_ += bopname;
		name_ += "', '";
		name_ += topname;
	}


template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::Name( const char* s,
                                      const char* opname,
                                      const char* bopname,
                                      const char* topname ) noexcept
	{
		name_ = s;
		name_ += ": Operand: '";
		name_ += opname;
		name_ += "', '";
		name_ += bopname;
		name_ += "', '";
		name_ += topname;
	}



template<uint32_t dim, template<uint32_t> class CELL>
string MathOperatorLHS<dim,CELL>::Name() const noexcept
	{
		return name_;
	}



template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::Out()  const
	{
		op.Out();
		bop.first.Out();
		top.first.Out();
	}


	// MATERIAL PROPERTY OPERAND

template<uint32_t dim, template<uint32_t> class CELL>
const csmp::Parameter&   MathOperatorLHS<dim,CELL>::MaterialOperand() const noexcept
	{
		return op;
	}



template<uint32_t dim, template<uint32_t> class CELL>
const csmp::Index&   MathOperatorLHS<dim,CELL>::MaterialOperandKey() const noexcept
	{
		return op.key;
	}



template<uint32_t dim, template<uint32_t> class CELL>
string  MathOperatorLHS<dim,CELL>::MaterialOperandName() const noexcept
	{
		return op.name;
	}



template<uint32_t dim, template<uint32_t> class CELL>
csmp::VARIABLE_TYPE   MathOperatorLHS<dim,CELL>::MaterialOperandType() const noexcept
	{
		return op.key.type;
	}



template<uint32_t dim, template<uint32_t> class CELL>
csmp::PLACEMENT   MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() const noexcept
	{
		return op.key.place;
	}



template<uint32_t dim, template<uint32_t> class CELL>
uint32_t   MathOperatorLHS<dim,CELL>::MaterialOperandDataDepth() const noexcept
	{
		return op.key.dataDepth;
	}


	// BASIC FUNCTION OPERAND

template<uint32_t dim, template<uint32_t> class CELL>
	const csmp::Parameter&   MathOperatorLHS<dim,CELL>::BasicOperand() const noexcept
	{
		return bop.first;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	const csmp::Index&   MathOperatorLHS<dim,CELL>::BasicOperandKey() const noexcept
	{
		return bop.first.key;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	string  MathOperatorLHS<dim,CELL>::BasicOperandName() const noexcept
	{
		return bop.first.name;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	csmp::VARIABLE_TYPE   MathOperatorLHS<dim,CELL>::BasicOperandType() const noexcept
	{
		return bop.first.key.type;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	csmp::PLACEMENT   MathOperatorLHS<dim,CELL>::BasicOperandPlacement() const noexcept
	{
		return bop.first.key.place;
	}

template<uint32_t dim, template<uint32_t> class CELL>
uint32_t   MathOperatorLHS<dim,CELL>::BasicOperandDataDepth() const noexcept
	{
		return bop.first.key.dataDepth;
	}

template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::BasicOperandOffset(size_t os) noexcept
	{
		bop.second = os;
	}

template<uint32_t dim, template<uint32_t> class CELL>
size_t MathOperatorLHS<dim,CELL>::BasicOperandOffset() const noexcept
	{
		return bop.second;
	}


	// TEST FUNCTION OPERAND

template<uint32_t dim, template<uint32_t> class CELL>
	const csmp::Parameter&  MathOperatorLHS<dim,CELL>::TestOperand() const noexcept
	{
		return top.first;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	const csmp::Index&  MathOperatorLHS<dim,CELL>::TestOperandKey() const noexcept
	{
		return top.first.key;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	string  MathOperatorLHS<dim,CELL>::TestOperandName() const noexcept
	{
		return top.first.name;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	csmp::VARIABLE_TYPE   MathOperatorLHS<dim,CELL>::TestOperandType() const noexcept
	{
		return top.first.key.type;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	csmp::PLACEMENT   MathOperatorLHS<dim,CELL>::TestOperandPlacement() const noexcept
	{
		return top.first.key.place;
	}

template<uint32_t dim, template<uint32_t> class CELL>
uint32_t  MathOperatorLHS<dim,CELL>::TestOperandDataDepth() const noexcept
	{
		return top.first.key.dataDepth;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::TestOperandOffset(size_t os) noexcept
	{
		top.second = os;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	size_t MathOperatorLHS<dim,CELL>::TestOperandOffset() const noexcept
	{
		return top.second;
	}





	/// ===================================  Accumulation Process Settings ========================================================================


template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::AddAccumulate() noexcept
	{
		add_accumulate_ = true;
		subtract_accumulate_ = false;
		add_accumulate_later_ = false;
		subtract_accumulate_later_ = false;
		multiply_accumulate_ = false;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::SubtractAccumulate() noexcept
	{
		add_accumulate_ = false;
		subtract_accumulate_ = true;
		add_accumulate_later_ = false;
		subtract_accumulate_later_ = false;
		multiply_accumulate_ = false;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::AddAccumulateLater() noexcept
	{
		add_accumulate_ = false;
		subtract_accumulate_ = false;
		add_accumulate_later_ = true;
		subtract_accumulate_later_ = false;
		multiply_accumulate_ = false;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::SubtractAccumulateLater() noexcept
	{
		add_accumulate_ = false;
		subtract_accumulate_ = false;
		add_accumulate_later_ = false;
		subtract_accumulate_later_ = true;
		multiply_accumulate_ = false;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::MultiplyAccumulate() noexcept
	{
		add_accumulate_ = false;
		subtract_accumulate_ = false;
		add_accumulate_later_ = false;
		subtract_accumulate_later_ = false;
		multiply_accumulate_ = true;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::LumpedFormulation(bool lumped) noexcept
	{
		lump_matrices_ = lumped;
	}


template<uint32_t dim, template<uint32_t> class CELL>
	void  MathOperatorLHS<dim,CELL>::ApplicationCycles( uint32_t c) noexcept
	{
		application_cycles_ = c;
	}


template<uint32_t dim, template<uint32_t> class CELL>
	void  MathOperatorLHS<dim,CELL>::ApplicationCycle( uint32_t c) noexcept
	{
		application_cycle_ = c;
	}


template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::MultiplyBy(double multiplication_factor) noexcept
	{
		factor_ = multiplication_factor;
	}


template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::MultiplyWithTimeIncrement(bool multiply) noexcept
	{
		time_multiply_ = multiply;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::DivideByTimeIncrement(bool divide) noexcept
	{
		time_divide_ = divide;
	}



template<uint32_t dim, template<uint32_t> class CELL>
bool MathOperatorLHS<dim,CELL>::Add() const noexcept
	{
		return add_accumulate_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
bool MathOperatorLHS<dim,CELL>::Subtract() const noexcept
	{
		return subtract_accumulate_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::AddLater() const noexcept
	{
		return add_accumulate_later_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::SubtractLater() const noexcept
	{
		return subtract_accumulate_later_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::Multiply() const noexcept
	{
		return multiply_accumulate_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool  MathOperatorLHS<dim,CELL>::LumpedFormulation() const noexcept
	{
		return lump_matrices_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
uint32_t MathOperatorLHS<dim,CELL>::ApplicationCycles() const noexcept
	{
		return application_cycles_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
uint32_t  MathOperatorLHS<dim,CELL>::ApplicationCycle() const noexcept
	{
		return application_cycle_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	double   MathOperatorLHS<dim,CELL>::MultiplyBy() const noexcept
	{
		return factor_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::MultiplyWithTimeIncrement() const noexcept
	{
		return time_multiply_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::DivideByTimeIncrement() const noexcept
	{
		return time_divide_;
	}


	/// ===================================  Accumulation Process Settings ========================================================================

	/// interpolation of property if isoparametric elements are used
template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::PropertyAtIntegrationPoint( const CELL<dim>& e_ref,
                                                            const csmp::Index&  idx,
                                                            uint32_t ip, DenseMatrix<dim>& M )
	{
		if (!e_ref.UsesLocalCoordinates())
			throw csmp::Exception(FATAL_ERROR,
				"MathOperatorLHS<dim,CELL>::PropertyValueAtIntegrationPoint",
				"Current element does not support numerical integration.");

		// 0. If the property is a constraint point variable a one-to-one mapping
		//    can be performed; no interpolation is needed
		if (idx.place == ELEMENT_INTEGRATION_POINT) {
        if (idx.type == SCALAR) {
          M.AssignToDiagonalAndZeroOffDiagonal(dim, e_ref.Read(ip, idx));
        }
			else if (idx.type == VECTOR) {
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
          ArrayVariable ar(idx.dataDepth);
          e_ref.Read(ip, idx, ar );
          M.AssignToDiagonal( ar );
        }
			else if (idx.type == FLAGGEDARRAY)
        {
          FlaggedArrayVariable fr(idx.dataDepth);
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
          ArrayVariable ar(idx.dataDepth);
          e_ref.Read(idx, ar );
          M.AssignToDiagonal( ar );
        }
			else if (idx.type == FLAGGEDARRAY)
        {
          FlaggedArrayVariable fr(idx.dataDepth);
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
          ArrayVariable ar(idx.dataDepth);
          e_ref.PropertyValueAtIntegrationPoint(idx, ip, ar );
          M.AssignToDiagonal( ar );
        }
			else if (idx.type == FLAGGEDARRAY)
        {
          FlaggedArrayVariable fr(idx.dataDepth);
          e_ref.PropertyValueAtIntegrationPoint(idx, ip, fr );
          M.AssignToDiagonal( fr );
        }
			return;
		}

		throw csmp::Exception(FATAL_ERROR,
			"MathOperatorLHS<dim,CELL>::PropertyValueAtIntegrationPoint",
			"This method cannot interpolate IntegrationPoint variables.");

	} // end PropertyValueAtIntegrationPoint( Element )








	  /**

	  GetOperands fills a vector of material property matrices with the
	  required values for later computation. These matrices always have
	  the dimensions spatial-dimension^2 and they will hold either scalar,
	  vector or tensor properties, depending on what kind of property the
	  Operand is.

	  When the property is an element property, it will be put into the
	  first vector entry MTRL[0].
*/
template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::GetOperands( const CELL<dim>& e_ref )
 {
    // if operand property is an element property
    if ( MaterialOperandPlacement() == ELEMENT || MaterialOperandPlacement() == FACE || MaterialOperandPlacement() == REGION )
      {
        MTRL.resize(1U);
        switch ( MaterialOperandType() )
          {
            case SCALAR:
              MTRL[0].AssignToDiagonalAndZeroOffDiagonal(dim, e_ref.Read(MaterialOperandKey()));
              break;
            case VECTOR:
              {
                VectorVariable<dim>  vc;
                e_ref.Read(MaterialOperandKey(), vc);
                MTRL[0].AssignToDiagonal(vc);
              }
              break;
            case TENSOR:
              {
                TensorVariable<dim>  ts;
                e_ref.Read(MaterialOperandKey(), ts);
                MTRL[0] = ts;
              }
              break;
            default:
              throw csmp::Exception(ERROR,
                "MathOperatorLHS<dim>::GetOperands",
                "Unsupported MaterialOperandType for ELEMENT/FACE/REGION placement");
          }
      }
    else // if the operand is placed on the constraint-points
      {
        const auto n_integration_points{ e_ref.IntegrationPoints() };

        if ( MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT || MaterialOperandPlacement() == FACE_INTEGRATION_POINT )
          {
            MTRL.resize( n_integration_points );
            for ( uint32_t i{0U}; i < n_integration_points; i++ )
              {
                switch ( MaterialOperandType() )
                  {
                    case SCALAR:
                      MTRL[i].AssignToDiagonalAndZeroOffDiagonal(dim, e_ref.Read(i, MaterialOperandKey()));
                      break;
                    case VECTOR:
                      {
                        VectorVariable<dim>  vc;
                        e_ref.Read(i, MaterialOperandKey(), vc);
                        MTRL[i].AssignToDiagonal(vc);
                      }
                      break;
                    case TENSOR:
                      {
                        TensorVariable<dim>  ts;
                        e_ref.Read(i, MaterialOperandKey(), ts);
                        MTRL[i] = ts;
                      }
                      break;
                    default:
                      throw csmp::Exception(ERROR,
                        "MathOperatorLHS<dim>::GetOperands",
                        "Unsupported MaterialOperandType for ELEMENT_INTEGRATION_POINT/FACE_INTEGRATION_POINT placement");
                  }
              }
          }
        // if the operand is placed on the node
        else if ( MaterialOperandPlacement() == NODE )
          {
            for ( uint32_t i{0U}; i < n_integration_points; i++ )
              PropertyAtIntegrationPoint(e_ref, MaterialOperandKey(), i, MTRL[i]);
          }
        else
          throw csmp::Exception(FATAL_ERROR,
            "MathOperatorLHS<dim>::GetOperands",
            "InterFace based operands cannot be accumulated with this method");
      }

 } // end GetOperands(CELL)







template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::WriteOperands( CELL<dim>& e )
	{
		cerr << "\nMathOperatorLHS<dim,CELL>::WriteOperands: ";
		cerr << " Overload to write LHS operands from element: " << e.Idx() << endl;
		throw invalid_argument("MathOperatorLHS<dim,CELL>::WriteOperands(Element)");
	} // end WriteOperands(Element)



template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::ComputeContribution( const CELL<dim>& e )
	{
		cerr << "\nMathOperatorLHS<dim>::ComputeContribution: ";
		cerr << " Overload to calculate LHS contribution from Element: " << e.Idx() << endl;
		throw invalid_argument("MathOperatorLHS<dim>::ComputeContribution(Element)");
	}// end ComputeContribution(Element)




/**
 * @brief Accumulates element LHS matrix into global SparseMatrix
 *        WITHOUT Dirichlet elimination.
 *
 * No pivot vector is needed — all DOFs including constrained ones
 * are assembled directly into G.
 *
 * Rows    -> test  operand DOFs (IDT)
 * Columns -> basic operand DOFs (IDB)
 *
 * @attention Do not use this legacy method unless you do not want full elimination of Dirichlet boundary conditions
 */
template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::AssignToGlobal( const CELL<dim>& e,
                                                SparseMatrix& G )
{
    const auto n_nodes{ e.Nodes() };

    // --- Build IDB (basic/column DOFs) and IDT (test/row DOFs) ---
    vector<size_t> IDB( n_nodes );
    vector<size_t> IDT( n_nodes );
    for ( uint32_t n{0U}; n < n_nodes; ++n )
        IDB[n] = IDT[n] = e.N(n)->Idx();

    if ( BasicOperandType() != SCALAR )
        transformNodeIndexVector( BasicOperandKey(), IDB );
    for ( uint32_t i{0U}; i < IDB.size(); ++i )
        IDB[i] += this->BasicOperandOffset();

    if ( TestOperandType() != SCALAR )
        transformNodeIndexVector( TestOperandKey(), IDT );
    for ( uint32_t i{0U}; i < IDT.size(); ++i )
        IDT[i] += this->TestOperandOffset();

#ifdef DEBUG
    csmp::ErrorHandler& csmp_err( csmp::ErrorHandler::Instance() );
    if ( LHS.Cols() != IDB.size() ) {
        cout << "\nbasic operand '" << BasicOperandName()
             << "' (" << parseType(BasicOperandType()) << ")"
             << "; element nodes " << n_nodes
             << "; LHS dimensions: " << LHS.Rows() << " x " << LHS.Cols() << "\n";
        csmp_err.Note( ERROR, "MathOperatorLHS<dim,CELL>::AssignToGlobal:",
                       "number of matrix COLUMNS does not match nodes x nodal DOFs." );
        return;
    }
    if ( LHS.Rows() != IDT.size() ) {
        cout << "\ntest operand '" << TestOperandName()
             << "' (" << parseType(TestOperandType()) << ")"
             << "; element nodes " << n_nodes
             << "; LHS dimensions: " << LHS.Rows() << " x " << LHS.Cols() << "\n";
        csmp_err.Note( ERROR, "MathOperatorLHS<dim,CELL>::AssignToGlobal:",
                       "number of matrix ROWS does not match nodes x nodal DOFs." );
        return;
    }
#endif

    // --- Scatter: rows from IDT, columns from IDB ---
    if ( multiply_accumulate_ )
    {
        for ( uint32_t i{0U}; i < LHS.Rows(); ++i )
            for ( uint32_t j{0U}; j < LHS.Cols(); ++j )
                G.MultiplyEntryWith( IDT[i], IDB[j], LHS(i,j) * factor_ );
    }
    else if ( add_accumulate_ || add_accumulate_later_ )
    {
        for ( uint32_t i{0U}; i < LHS.Rows(); ++i )
            for ( uint32_t j{0U}; j < LHS.Cols(); ++j )
                G.Add( IDT[i], IDB[j], LHS(i,j) * factor_ );
    }
    else if ( subtract_accumulate_ || subtract_accumulate_later_ )
    {
        for ( uint32_t i{0U}; i < LHS.Rows(); ++i )
            for ( uint32_t j{0U}; j < LHS.Cols(); ++j )
                G.Add( IDT[i], IDB[j], -LHS(i,j) * factor_ );
    }
    else
        throw csmp::Exception( ERROR,
                               "MathOperatorLHS<dim>::AssignToGlobal(SparseMatrix, no elimination):",
                               "Accumulation mode not recognised." );

} // end AssignToGlobal (SparseMatrix, no Dirichlet elimination)





// SKM fixed for coupling operators 1/7/2024
template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::AssignToGlobal( const CELL<dim>& e,
                                                 SparseMatrix& G,
                                                 vector<double>& pivotVector,
                                                 const vector<size_t>& DOF_indexes )
{
    const auto n_nodes{ e.Nodes() };

    // --- Build IDB and IDT ---
    vector<size_t> IDB( n_nodes );
    vector<size_t> IDT( n_nodes );
    for ( uint32_t n{0U}; n < n_nodes; ++n )
        IDB[n] = IDT[n] = e.N(n)->Idx();

    if ( BasicOperandType() != SCALAR )
        transformNodeIndexVector( BasicOperandKey(), IDB );
    for ( uint32_t i{0U}; i < IDB.size(); ++i ) {
        IDB[i] += this->BasicOperandOffset();
        IDB[i]  = DOF_indexes[ IDB[i] ];
    }

    if ( TestOperandType() != SCALAR )
        transformNodeIndexVector( TestOperandKey(), IDT );
    for ( uint32_t i{0U}; i < IDT.size(); ++i ) {
        IDT[i] += this->TestOperandOffset();
        IDT[i]  = DOF_indexes[ IDT[i] ];
    }

#ifdef DEBUG
    csmp::ErrorHandler& csmp_err( csmp::ErrorHandler::Instance() );
    if ( LHS.Cols() != IDB.size() ) {
        cout << "\nbasic operand '" << BasicOperandName()
             << "' (" << parseType(BasicOperandType()) << ")"
             << "; element nodes " << n_nodes
             << "; LHS dimensions: " << LHS.Rows() << " x " << LHS.Cols() << "\n";
        csmp_err.Note( ERROR, "MathOperatorLHS<dim,CELL>::AssignToGlobal:",
                       "number of matrix COLUMNS does not match nodes x nodal DOFs." );
        return;
    }
    if ( LHS.Rows() != IDT.size() ) {
        cout << "\ntest operand '" << TestOperandName()
             << "' (" << parseType(TestOperandType()) << ")"
             << "; element nodes " << n_nodes
             << "; LHS dimensions: " << LHS.Rows() << " x " << LHS.Cols() << "\n";
        csmp_err.Note( ERROR, "MathOperatorLHS<dim,CELL>::AssignToGlobal:",
                       "number of matrix ROWS does not match nodes x nodal DOFs." );
        return;
    }
#endif

    // --- Prescribed values of the BASIC operand ---
    vector<double> top_variable_node_values( IDB.size() );
    ReadBasicOperandValues( e, n_nodes, top_variable_node_values );

    // --- Scatter ---
    if ( multiply_accumulate_ )
    {
        throw csmp::Exception( ERROR,
                               "MathOperatorLHS<dim>::AssignToGlobal(SparseMatrix):",
                               "Multiply Accumulate not supported yet" );
    }
    else if ( add_accumulate_ || add_accumulate_later_ )
    {
        assert( IDB.size() == LHS.Cols() );
        assert( IDT.size() == LHS.Rows() );

        for ( uint32_t i{0U}; i < LHS.Rows(); ++i ) {
            if ( IDT[i] == NULL_IDX ) continue;
            for ( uint32_t j{0U}; j < LHS.Cols(); ++j ) {
                if ( IDB[j] == NULL_IDX )
                    pivotVector[ IDT[i] ] -= LHS(i,j) * top_variable_node_values[j] * factor_;
                else
                    G.Add( IDT[i], IDB[j], LHS(i,j) * factor_ );
            }
        }
    }
    else if ( subtract_accumulate_ || subtract_accumulate_later_ )
    {
        assert( IDB.size() == LHS.Cols() );
        assert( IDT.size() == LHS.Rows() );

        for ( uint32_t i{0U}; i < LHS.Rows(); ++i ) {
            if ( IDT[i] == NULL_IDX ) continue;
            for ( uint32_t j{0U}; j < LHS.Cols(); ++j ) {
                if ( IDB[j] == NULL_IDX )
                    pivotVector[ IDT[i] ] += LHS(i,j) * top_variable_node_values[j] * factor_;
                else
                    G.Add( IDT[i], IDB[j], LHS(i,j) * factor_ );
            }
        }
    }
    else
        throw csmp::Exception( ERROR,
                               "MathOperatorLHS<dim>::AssignToGlobal(SparseMatrix):",
                               "Accumulation mode not recognised." );

} // end AssignToGlobal (SparseMatrix)





/**
    Accumulation without elimination of Dirichlet degrees of freedom
*/
template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::AssignToGlobal( const CELL<dim>& e,
                                                CompressedRowMatrix& G,
                                                vector<double>& pivotVector,
                                                const vector<size_t>& DOF_indexes )
{
    const auto n_nodes{ e.Nodes() };

    // --- Build IDB (basic/column DOF indices) and IDT (test/row DOF indices) ---
    vector<size_t> IDB( n_nodes );
    vector<size_t> IDT( n_nodes );
    for ( uint32_t n{0U}; n < n_nodes; ++n )
        IDB[n] = IDT[n] = e.N(n)->Idx();

    if ( BasicOperandType() != SCALAR )
        transformNodeIndexVector( BasicOperandKey(), IDB );
    for ( uint32_t i{0U}; i < IDB.size(); ++i ) {
        IDB[i] += this->BasicOperandOffset();
        IDB[i]  = DOF_indexes[ IDB[i] ];
    }

    if ( TestOperandType() != SCALAR )
        transformNodeIndexVector( TestOperandKey(), IDT );
    for ( uint32_t i{0U}; i < IDT.size(); ++i ) {
        IDT[i] += this->TestOperandOffset();
        IDT[i]  = DOF_indexes[ IDT[i] ];
    }

#ifdef DEBUG
    csmp::ErrorHandler& csmp_err( csmp::ErrorHandler::Instance() );
    if ( LHS.Cols() != IDB.size() ) {
        cout << "\nbasic operand '" << BasicOperandName()
             << "' (" << parseType(BasicOperandType()) << ")"
             << "; element nodes " << n_nodes
             << "; LHS dimensions: " << LHS.Rows() << " x " << LHS.Cols() << "\n";
        csmp_err.Note( ERROR, "MathOperatorLHS<dim,CELL>::AssignToGlobal:",
                       "number of matrix COLUMNS does not match nodes x nodal DOFs." );
        return;
    }
    if ( LHS.Rows() != IDT.size() ) {
        cout << "\ntest operand '" << TestOperandName()
             << "' (" << parseType(TestOperandType()) << ")"
             << "; element nodes " << n_nodes
             << "; LHS dimensions: " << LHS.Rows() << " x " << LHS.Cols() << "\n";
        csmp_err.Note( ERROR, "MathOperatorLHS<dim,CELL>::AssignToGlobal:",
                       "number of matrix ROWS does not match nodes x nodal DOFs." );
        return;
    }
#endif

    // --- Prescribed values of the BASIC operand (column variable) ---
    // Sized to IDB.size(): one entry per column DOF.
    // Used to modify the pivot vector when a basic DOF is Dirichlet-constrained.
    vector<double> top_variable_node_values( IDB.size() );
    ReadBasicOperandValues( e, n_nodes, top_variable_node_values );

    // --- Scatter LHS into global matrix and pivot vector ---
    if ( multiply_accumulate_ )
    {
        throw csmp::Exception( ERROR,
                               "MathOperatorLHS<dim>::AssignToGlobal(CompressedRowMatrix):",
                               "Multiply Accumulate not supported yet" );
    }
    else if ( add_accumulate_ || add_accumulate_later_ )
    {
        assert( IDB.size() == LHS.Cols() );
        assert( IDT.size() == LHS.Rows() );

        for ( uint32_t i{0U}; i < LHS.Rows(); ++i ) {
            if ( IDT[i] == NULL_IDX ) continue;
            for ( uint32_t j{0U}; j < LHS.Cols(); ++j ) {
                if ( IDB[j] == NULL_IDX )
                    pivotVector[ IDT[i] ] -= LHS(i,j) * top_variable_node_values[j] * factor_;
                else // if there is a pre-allocated space in the CRM
                    G.AddIf( IDT[i], IDB[j], LHS(i,j) * factor_ );
            }
        }
    }
    else if ( subtract_accumulate_ || subtract_accumulate_later_ )
    {
        assert( IDB.size() == LHS.Cols() );
        assert( IDT.size() == LHS.Rows() );

        for ( uint32_t i{0U}; i < LHS.Rows(); ++i ) {
            if ( IDT[i] == NULL_IDX ) continue;
            for ( uint32_t j{0U}; j < LHS.Cols(); ++j ) {
                if ( IDB[j] == NULL_IDX )
                    pivotVector[ IDT[i] ] += LHS(i,j) * top_variable_node_values[j] * factor_;
                else
                    G.AddIf( IDT[i], IDB[j], LHS(i,j) * factor_ );
            }
        }
    }
    else
        throw csmp::Exception( ERROR,
                               "MathOperatorLHS<dim>::AssignToGlobal(CompressedRowMatrix):",
                               "Accumulation mode not recognised." );

} // end AssignToGlobal (CompressedRowMatrix)



template class MathOperatorLHS<1U>;
template class MathOperatorLHS<2U>;
template class MathOperatorLHS<3U>;

template class MathOperatorLHS<1U,Face>;
template class MathOperatorLHS<2U,Face>;
template class MathOperatorLHS<3U,Face>;

// TODO: this instantiation is only necessary if FE computations are done on InterFace objects
template class MathOperatorLHS<1U,InterFace>;
template class MathOperatorLHS<2U,InterFace>;
template class MathOperatorLHS<3U,InterFace>;

} // end namespace csmp
