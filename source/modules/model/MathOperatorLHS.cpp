#include "MathOperatorLHS.h"
#include "MathOperatorRHS.h"
#include "Element.h"
#include "InterFace.h"
#include "Face.h"
#include "PropertyDatabase.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
MathOperatorLHS<dim,CELL>::MathOperatorLHS()
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



template<uint32_t dim, template<uint32_t> class CELL>
MathOperatorLHS<dim,CELL>::MathOperatorLHS( const PropertyDatabase<dim>& pref,
                                            const char* basic,
                                            const char* test)
		: name_("unspecified LHS-operator"),
		bop(make_pair(pref.Parameter(basic), 0)),
		top(make_pair(pref.Parameter(test), 0)),
		IDT(3),
		IDB(3),
		MTRL(13, DenseMatrix<DM_MIN>(2, 2) ),
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
		IDT(3),
		IDB(3),
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
MathOperatorLHS<dim,CELL>::MathOperatorLHS( const MathOperatorLHS<dim,CELL>& mo )
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
			IDT = mo.IDT;   // node-ID & global constraint points vector
			IDB = mo.IDB;   // node-ID & global constraint points vector
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


template<uint32_t dim, template<uint32_t> class CELL>
MathOperatorLHS<dim,CELL>::~MathOperatorLHS()
	{
	}


	/// ===================================  Operand Functions ========================================================================

template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::Name( const char* s,
                                      const char* bopname,
                                      const char* topname)
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
                                      const char* topname )
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
string MathOperatorLHS<dim,CELL>::Name() const
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
const csmp::Parameter&   MathOperatorLHS<dim,CELL>::MaterialOperand() const
	{
		return op;
	}



template<uint32_t dim, template<uint32_t> class CELL>
const csmp::Index&   MathOperatorLHS<dim,CELL>::MaterialOperandKey() const
	{
		return op.key;
	}



template<uint32_t dim, template<uint32_t> class CELL>
string  MathOperatorLHS<dim,CELL>::MaterialOperandName() const
	{
		return op.name;
	}



template<uint32_t dim, template<uint32_t> class CELL>
csmp::VARIABLE_TYPE   MathOperatorLHS<dim,CELL>::MaterialOperandType() const
	{
		return op.key.type;
	}



template<uint32_t dim, template<uint32_t> class CELL>
csmp::PLACEMENT   MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() const
	{
		return op.key.place;
	}



template<uint32_t dim, template<uint32_t> class CELL>
uint32_t   MathOperatorLHS<dim,CELL>::MaterialOperandDataDepth() const
	{
		return op.key.dataDepth;
	}


	// BASIC FUNCTION OPERAND

template<uint32_t dim, template<uint32_t> class CELL>
	const csmp::Parameter&   MathOperatorLHS<dim,CELL>::BasicOperand() const
	{
		return bop.first;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	const csmp::Index&   MathOperatorLHS<dim,CELL>::BasicOperandKey() const
	{
		return bop.first.key;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	string  MathOperatorLHS<dim,CELL>::BasicOperandName() const
	{
		return bop.first.name;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	csmp::VARIABLE_TYPE   MathOperatorLHS<dim,CELL>::BasicOperandType() const
	{
		return bop.first.key.type;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	csmp::PLACEMENT   MathOperatorLHS<dim,CELL>::BasicOperandPlacement() const
	{
		return bop.first.key.place;
	}

template<uint32_t dim, template<uint32_t> class CELL>
uint32_t   MathOperatorLHS<dim,CELL>::BasicOperandDataDepth() const
	{
		return bop.first.key.dataDepth;
	}

template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::BasicOperandOffset(size_t os)
	{
		bop.second = os;
	}

template<uint32_t dim, template<uint32_t> class CELL>
size_t MathOperatorLHS<dim,CELL>::BasicOperandOffset()
		const
	{
		return bop.second;
	}


	// TEST FUNCTION OPERAND

template<uint32_t dim, template<uint32_t> class CELL>
	const csmp::Parameter&  MathOperatorLHS<dim,CELL>::TestOperand() const
	{
		return top.first;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	const csmp::Index&  MathOperatorLHS<dim,CELL>::TestOperandKey() const
	{
		return top.first.key;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	string  MathOperatorLHS<dim,CELL>::TestOperandName() const
	{
		return top.first.name;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	csmp::VARIABLE_TYPE   MathOperatorLHS<dim,CELL>::TestOperandType() const
	{
		return top.first.key.type;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	csmp::PLACEMENT   MathOperatorLHS<dim,CELL>::TestOperandPlacement() const
	{
		return top.first.key.place;
	}

template<uint32_t dim, template<uint32_t> class CELL>
uint32_t  MathOperatorLHS<dim,CELL>::TestOperandDataDepth() const
	{
		return top.first.key.dataDepth;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::TestOperandOffset(size_t os)
	{
		top.second = os;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	size_t MathOperatorLHS<dim,CELL>::TestOperandOffset() const
	{
		return top.second;
	}





	/// ===================================  Accumulation Process Settings ========================================================================


template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::AddAccumulate()
	{
		add_accumulate_ = true;
		subtract_accumulate_ = false;
		add_accumulate_later_ = false;
		subtract_accumulate_later_ = false;
		multiply_accumulate_ = false;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::SubtractAccumulate()
	{
		add_accumulate_ = false;
		subtract_accumulate_ = true;
		add_accumulate_later_ = false;
		subtract_accumulate_later_ = false;
		multiply_accumulate_ = false;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::AddAccumulateLater()
	{
		add_accumulate_ = false;
		subtract_accumulate_ = false;
		add_accumulate_later_ = true;
		subtract_accumulate_later_ = false;
		multiply_accumulate_ = false;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::SubtractAccumulateLater()
	{
		add_accumulate_ = false;
		subtract_accumulate_ = false;
		add_accumulate_later_ = false;
		subtract_accumulate_later_ = true;
		multiply_accumulate_ = false;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::MultiplyAccumulate()
	{
		add_accumulate_ = false;
		subtract_accumulate_ = false;
		add_accumulate_later_ = false;
		subtract_accumulate_later_ = false;
		multiply_accumulate_ = true;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::LumpedFormulation(bool lumped)
	{
		lump_matrices_ = lumped;
	}


template<uint32_t dim, template<uint32_t> class CELL>
	void  MathOperatorLHS<dim,CELL>::ApplicationCycles( uint32_t c)
	{
		application_cycles_ = c;
	}


template<uint32_t dim, template<uint32_t> class CELL>
	void  MathOperatorLHS<dim,CELL>::ApplicationCycle( uint32_t c)
	{
		application_cycle_ = c;
	}


template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::MultiplyBy(double multiplication_factor)
	{
		factor_ = multiplication_factor;
	}


template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::MultiplyWithTimeIncrement(bool multiply)
	{
		time_multiply_ = multiply;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	void MathOperatorLHS<dim,CELL>::DivideByTimeIncrement(bool divide)
	{
		time_divide_ = divide;
	}



template<uint32_t dim, template<uint32_t> class CELL>
bool MathOperatorLHS<dim,CELL>::Add() const
	{
		return add_accumulate_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
bool MathOperatorLHS<dim,CELL>::Subtract() const
	{
		return subtract_accumulate_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::AddLater() const
	{
		return add_accumulate_later_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::SubtractLater() const
	{
		return subtract_accumulate_later_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::Multiply() const
	{
		return multiply_accumulate_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool  MathOperatorLHS<dim,CELL>::LumpedFormulation() const
	{
		return lump_matrices_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
uint32_t MathOperatorLHS<dim,CELL>::ApplicationCycles() const
	{
		return application_cycles_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
uint32_t  MathOperatorLHS<dim,CELL>::ApplicationCycle() const
	{
		return application_cycle_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	double   MathOperatorLHS<dim,CELL>::MultiplyBy() const
	{
		return factor_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::MultiplyWithTimeIncrement() const
	{
		return time_multiply_;
	}

template<uint32_t dim, template<uint32_t> class CELL>
	bool MathOperatorLHS<dim,CELL>::DivideByTimeIncrement() const
	{
		return time_divide_;
	}


	/// ===================================  Accumulation Process Settings ========================================================================

	/// interpolation of property if isoparametric elements are used
template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::PropertyAtIntegrationPoint( const CELL<dim>& e_ref,
                                                            const csmp::Index&  idx,
                                                            uint32_t ip, DenseMatrix<DM_MIN>& M )
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
 // TODO: put static assert in here as needed
		// if operand property is an element property
		if (MaterialOperandPlacement() == ELEMENT || MaterialOperandPlacement() == REGION)
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
        
        if (MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT) {
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
				for (auto i{0U}; i < n_integration_points; i++)
					PropertyAtIntegrationPoint(e_ref, MaterialOperandKey(), i, MTRL[i]);
			}
			else
				throw csmp::Exception(FATAL_ERROR,
					"MathOperatorLHS<dim>::GetOperands",
					"Face based operands cannot be accumulated with this method");
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






template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::MultiplyWithTimeFactor( double dt )
	{
		LHS *= dt;
	}




	/// AssignToGlobal matrix function
template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::AssignToGlobal( const CELL<dim>& e, SparseMatrix& G,
                                                vector<double>& pivotVector, const vector<size_t>& DOF_indexes )
	{
		// map local to global indexes for test and basic operands
		IDT.resize(e.Nodes()); // ii
		IDB.resize(e.Nodes()); // jj
		for (auto i{0U}; i < e.Nodes(); i++) {
			IDT[i] = e.N(i)->Idx();
			IDB[i] = IDT[i];
		}

		if (TestOperandType() != SCALAR)
			transformNodeIndexVector(dim, TestOperandKey(), IDT);

		for (auto i{0U}; i < IDT.size(); i++) {
			IDT[i] += this->TestOperandOffset();
			IDT[i] = DOF_indexes[IDT[i]];
		}

		if (BasicOperandType() != SCALAR)
			transformNodeIndexVector(dim, BasicOperandKey(), IDB);

		for (auto i{0U}; i < IDB.size(); i++) {
			IDB[i] += this->BasicOperandOffset();
			IDB[i] = DOF_indexes[IDB[i]];
		}

		// get (component of) local value of a local node 
		vector<double> nodal_values(IDB.size());
		{
			if (this->TestOperandType() == SCALAR) {
				for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
					nodal_values[nIdx] = e.N(nIdx)->Read(this->TestOperandKey());
				}
			}
			else if (this->TestOperandType() == VECTOR) {
				VectorVariable<dim> var;
				for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
					for (auto i = 0; i < dim; ++i) {
						e.N(nIdx)->Read(this->TestOperandKey(), var);
						nodal_values[nIdx*dim + i] = var.Component(i);
					}
				}
			}
			else if (this->TestOperandType() == TENSOR) {
				uint32_t dim2 = dim*dim;
				TensorVariable<dim> var;
				for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
					e.N(nIdx)->Read(this->TestOperandKey(), var);
					for (auto i = 0; i < dim; ++i) {
						for (auto j = 0; j < dim; ++j) {
							nodal_values[nIdx*dim2 + i*dim + j] = var.Component(i*dim + j);
						}
					}
				}
			}
			else if (this->TestOperandType() == ARRAY) {
				const auto datadepth(TestOperandKey().dataDepth);
				ArrayVariable var;
				for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
					e.N(nIdx)->Read(this->TestOperandKey(), var);
					for (auto i = 0; i < datadepth; ++i) {
						nodal_values[nIdx*datadepth + i] = var.Component(i);
					}
				}
			}
			else if (this->TestOperandType() == FLAGGEDARRAY) {
				const auto datadepth(TestOperandKey().dataDepth);
				FlaggedArrayVariable var;
				for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
					e.N(nIdx)->Read(this->TestOperandKey(), var);
					for (auto i = 0; i < datadepth; ++i) {
						nodal_values[nIdx*datadepth + i] = var.Component(i);
					}
				}
			}
			else {
				throw csmp::Exception(ERROR,
					"MathOperatorLHS<dim>::AssignToGlobal(Element):",
					"Undefined variable type");
			}
		};

		// perform assignment from local matrix to global matrix

		if (multiply_accumulate_)
		{
			throw csmp::Exception(ERROR,
				"MathOperatorLHS<dim>::AssignToGlobal(Element):",
				"Multiply Accumulate has not been supported yet");

			//if (IDT[i] != NULL_IDX) {
			//  for (size_t j{0U}; j<LHS.Cols(); j++)
			//    if (IDB[j] != NULL_IDX) {
			//      G.MultiplyEntryWith(IDT[i],
			//        IDB[j],
			//        LHS(i, j) * factor_);
			//    }
			//}
		}
		else if ( add_accumulate_ || add_accumulate_later_)
		{
			for (auto i{0U}; i < LHS.Rows(); i++) {
				if (IDT[i] != NULL_IDX) {
					for (auto j{0U}; j < LHS.Cols(); j++) {
						if ( IDB[j] == NULL_IDX ) {
							pivotVector[IDT[i]] -= LHS(i, j) * nodal_values[j] * factor_;  // notice the sign e.N(j / this->TestOperandOffset())->Read(this->TestOperand())
						}
						else {
							G.Add(IDT[i],
								IDB[j],
								LHS(i, j) * factor_);
						}
					}
				}
			}

			/*for (auto i{0U}; i<LHS.Rows(); i++)
			  if (IDT[i] != NULL_IDX) {
				for (size_t j{0U}; j < LHS.Cols(); j++)
				  if (IDB[j] != NULL_IDX) {
					G.Add(IDT[i],
					  IDB[j],
					  LHS(i, j) * factor_);
				  }
			  }*/
		}
		else if (subtract_accumulate_ || subtract_accumulate_later_)
		{
			for (auto i{0U}; i < LHS.Rows(); i++) {
				if (IDT[i] != NULL_IDX) {
					for (auto j{0U}; j < LHS.Cols(); j++) {
						if (IDB[j] == NULL_IDX) {
							pivotVector[IDT[i]] += LHS(i, j) * nodal_values[j] * factor_;   // notice the sign LHS(i, j) * e.N(j)->Read(TestOperandKey());
						}
						else {
							G.Add(IDT[i],
								IDB[j],
								-LHS(i, j) * factor_);
						}
					}
				}
			}

			/*for (auto i{0U}; i<LHS.Rows(); i++)
			  if (IDT[i] != NULL_IDX) {
				for (size_t j{0U}; j < LHS.Cols(); j++)
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




  /// AssignToGlobal matrix function
  template<uint32_t dim, template<uint32_t> class CELL>
  void MathOperatorLHS<dim,CELL>::AssignToGlobal( const CELL<dim>& e, CompressedRowMatrix& G,
                                                  vector<double>& pivotVector, const vector<size_t>& DOF_indexes )
  {
    // map local to global indexes for test and basic operands
    IDT.resize(e.Nodes()); // ii
    IDB.resize(e.Nodes()); // jj
    for (auto i{0U}; i < e.Nodes(); i++) {
      IDT[i] = e.N(i)->Idx();
      IDB[i] = IDT[i];
    }

    if (TestOperandType() != SCALAR)
      transformNodeIndexVector(dim, TestOperandKey(), IDT);

    for (auto i{0U}; i < IDT.size(); i++) {
      IDT[i] += this->TestOperandOffset();
      IDT[i] = DOF_indexes[IDT[i]];
    }

    if (BasicOperandType() != SCALAR)
      transformNodeIndexVector(dim, BasicOperandKey(), IDB);

    for (auto i{0U}; i < IDB.size(); i++) {
      IDB[i] += this->BasicOperandOffset();
      IDB[i] = DOF_indexes[IDB[i]];
    }

    // get (component of) local value of a local node
    vector<double> nodal_values(IDB.size());
    {
      if (this->TestOperandType() == SCALAR) {
        for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
          nodal_values[nIdx] = e.N(nIdx)->Read(this->TestOperandKey());
        }
      }
      else if (this->TestOperandType() == VECTOR) {
        VectorVariable<dim> var;
        for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
          for (auto i = 0; i < dim; ++i) {
            e.N(nIdx)->Read(this->TestOperandKey(), var);
            nodal_values[nIdx*dim + i] = var.Component(i);
          }
        }
      }
      else if (this->TestOperandType() == TENSOR) {
        uint32_t dim2 = dim*dim;
        TensorVariable<dim> var;
        for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
          e.N(nIdx)->Read(this->TestOperandKey(), var);
          for (auto i = 0; i < dim; ++i) {
            for (auto j = 0; j < dim; ++j) {
              nodal_values[nIdx*dim2 + i*dim + j] = var.Component(i*dim + j);
            }
          }
        }
      }
      else if (this->TestOperandType() == ARRAY) {
        const auto datadepth(TestOperandKey().dataDepth);
        ArrayVariable var;
        for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
          e.N(nIdx)->Read(this->TestOperandKey(), var);
          for (auto i = 0; i < datadepth; ++i) {
            nodal_values[nIdx*datadepth + i] = var.Component(i);
          }
        }
      }
      else if (this->TestOperandType() == FLAGGEDARRAY) {
        const auto datadepth(TestOperandKey().dataDepth);
        FlaggedArrayVariable var;
        for (auto nIdx = 0; nIdx < e.Nodes(); ++nIdx) {
          e.N(nIdx)->Read(this->TestOperandKey(), var);
          for (auto i = 0; i < datadepth; ++i) {
            nodal_values[nIdx*datadepth + i] = var.Component(i);
          }
        }
      }
      else {
        throw csmp::Exception(ERROR,
                              "MathOperatorLHS<dim>::AssignToGlobal(Element):",
                              "Undefined variable type");
      }
    };

    // perform assignment from local matrix to global matrix

    if (multiply_accumulate_)
    {
      throw csmp::Exception(ERROR,
                            "MathOperatorLHS<dim>::AssignToGlobal(Element):",
                            "Multiply Accumulate has not been supported yet");

      //if (IDT[i] != NULL_IDX) {
      //  for (size_t j{0U}; j<LHS.Cols(); j++)
      //    if (IDB[j] != NULL_IDX) {
      //      G.MultiplyEntryWith(IDT[i],
      //        IDB[j],
      //        LHS(i, j) * factor_);
      //    }
      //}
    }
    else if ( add_accumulate_ || add_accumulate_later_)
    {
      for (auto i{0U}; i < LHS.Rows(); i++) {
        if (IDT[i] != NULL_IDX) {
          for (auto j{0U}; j < LHS.Cols(); j++) {
            if ( IDB[j] == NULL_IDX ) {
              pivotVector[IDT[i]] -= LHS(i, j) * nodal_values[j] * factor_;  // notice the sign e.N(j / this->TestOperandOffset())->Read(this->TestOperand())
            }
            else {
              G.Add(IDT[i],
                    IDB[j],
                    LHS(i, j) * factor_);
            }
          }
        }
      }

      /*for (auto i{0U}; i<LHS.Rows(); i++)
        if (IDT[i] != NULL_IDX) {
        for (size_t j{0U}; j < LHS.Cols(); j++)
          if (IDB[j] != NULL_IDX) {
          G.Add(IDT[i],
            IDB[j],
            LHS(i, j) * factor_);
          }
        }*/
    }
    else if (subtract_accumulate_ || subtract_accumulate_later_)
    {
      for (auto i{0U}; i < LHS.Rows(); i++) {
        if (IDT[i] != NULL_IDX) {
          for (auto j{0U}; j < LHS.Cols(); j++) {
            if (IDB[j] == NULL_IDX) {
              pivotVector[IDT[i]] += LHS(i, j) * nodal_values[j] * factor_;   // notice the sign LHS(i, j) * e.N(j)->Read(TestOperandKey());
            }
            else {
              G.Add(IDT[i],
                    IDB[j],
                    -LHS(i, j) * factor_);
            }
          }
        }
      }

      /*for (auto i{0U}; i<LHS.Rows(); i++)
        if (IDT[i] != NULL_IDX) {
        for (size_t j{0U}; j < LHS.Cols(); j++)
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

	

template<uint32_t dim, template<uint32_t> class CELL>
void MathOperatorLHS<dim,CELL>::AssignToGlobal( const CELL<dim>& e, SparseMatrix& G )
	{
		// map local to global indexes for test and basic operands

		IDT.resize(e.Nodes());
		IDB.resize(e.Nodes());
		for (auto i{0U}; i < e.Nodes(); i++) {
			IDT[i] = e.N(i)->Idx();
			IDB[i] = IDT[i];
		}

		if (TestOperandType() != SCALAR)
			transformNodeIndexVector(dim, TestOperandKey(), IDT);

		for (auto i{0U}; i < IDT.size(); i++)
			IDT[i] += this->TestOperandOffset();

		if (BasicOperandType() != SCALAR)
			transformNodeIndexVector(dim, BasicOperandKey(), IDB);

		for (auto i{0U}; i < IDB.size(); i++)
			IDB[i] += this->BasicOperandOffset();

		// perform assignment from local matrix to global matrix

		if (multiply_accumulate_)
		{
			for (auto i{0U}; i < LHS.Rows(); i++)
				for (auto j{0U}; j < LHS.Cols(); j++)
					G.MultiplyEntryWith(IDT[i],
						IDB[j],
						LHS(i, j) * factor_);
		}
		else if (add_accumulate_ || add_accumulate_later_)
		{
			for (auto i{0U}; i < LHS.Rows(); i++)
				for (auto j{0U}; j < LHS.Cols(); j++)
					G.Add(IDT[i],
						IDB[j],
						LHS(i, j) * factor_);
		}
		else if (subtract_accumulate_ || subtract_accumulate_later_)
		{
			for (auto i{0U}; i < LHS.Rows(); i++)
				for (auto j{0U}; j < LHS.Cols(); j++)
					G.Add(IDT[i],
						IDB[j],
						-LHS(i, j) * factor_);
		}
		else
			throw csmp::Exception(ERROR,
				"MathOperatorLHS<dim>::AssignToGlobal(Element):",
				"accumulation instructions could not be parsed.");

	} // end AssignToGlobal (Element)



template class MathOperatorLHS<1U>;
template class MathOperatorLHS<2U>;
template class MathOperatorLHS<3U>;

template class MathOperatorLHS<1U,Face>;
template class MathOperatorLHS<2U,Face>;
template class MathOperatorLHS<3U,Face>;

// TODO: this instantiation should not be necessary
template class MathOperatorLHS<1U,InterFace>;
template class MathOperatorLHS<2U,InterFace>;
template class MathOperatorLHS<3U,InterFace>;

} // end namespace csmp
