#include "PDE_Integrator_CRM.h"

using namespace std;

namespace csmp {

	/*
	******************************************************************************
	PAERT 1: The following functions need modification
	******************************************************************************
	*/
	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::AssignEssentialConditions(const SIMPLICIAL_COMPLEX<dim>& gref)
	{
		if (!setup_established_)
			throw csmp::Exception(ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignEssentialConditions",
				"please call EstablishMatrixSetup() prior to this method.");

		if (basic_operands_.empty())
			throw csmp::Exception(ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignEssentialConditions",
				"No (basic) operands have been specified...");

		for (operandsConstIterator it = test_operands_.begin(); it != test_operands_.end(); it++)
		{
			typename vector<csmp::Node<dim>*>::const_iterator  niter(gref.NodesBegin());
			csmp::Index prop_key = (*it).first.key;
			size_t      offset = (*it).second;

			if (prop_key.place != NODE)
				throw csmp::Exception(ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignEssentialConditions: CRM:",
					"So far no conditions are assigned to elements, faces, segments");

			switch (prop_key.type)
			{
			case SCALAR:
				while (niter != gref.NodesEnd()) {
					if ((*niter)->Status(prop_key) == DIRICH)
						dirich_.push_back(Entry((*niter)->Idx() + offset , scale_factor_ * (*niter)->Read(prop_key)));
					niter++;
				}
				break;
			case VECTOR: {
				VectorVariable<dim>  vc;
				while (niter != gref.NodesEnd()) {
					(*niter)->Read(prop_key, vc);
					for (size_t i = 0U; i<dim; i++)
						if (vc.Flag(i) == DIRICH) {
							size_t  position = (*niter)->Idx() * dim + i + offset;
							cout << "FIX it --> PDE_Integrator_CRM 1\n"; getchar();
							//G_.ZeroRow(position);
							//G_.Add(position, position, scale_factor_);
							rh_[position] = scale_factor_ * vc(i);
						}
					niter++;
				}
			}
						 break;
			case TENSOR: {
				TensorVariable<dim>  ts;
				while (niter != gref.NodesEnd()) {
					(*niter)->Read(prop_key, ts);
					for (size_t i = 0U; i<dim; i++)
						if (ts.Flag(i) == DIRICH)
							for (size_t j = 0U; j<dim; j++)
							{
								size_t  position = (*niter)->Idx() * dim2_ + i * dim + j + offset;
								cout << "FIX it --> PDE_Integrator_CRM 2\n"; getchar();
								//G_.ZeroRow(position);
								//G_.Add(position, position, scale_factor_);
								rh_[position] = scale_factor_ * ts(i, j);
							}
					++niter;
				}
			}
						 break;
			case ARRAY: {
				ArrayVariable  ar(prop_key.dataDepth);
				while (niter != gref.NodesEnd()) {
					(*niter)->Read(prop_key, ar);
					if (ar.Flag() == DIRICH)
					{
						for (size_t i = 0U; i<prop_key.dataDepth; i++)
						{
							size_t  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
							cout << "FIX it --> PDE_Integrator_CRM 3\n"; getchar();
							//G_.ZeroRow(position);
							//G_.Add(position, position, scale_factor_);
							rh_[position] = scale_factor_ * ar(i);
						}
					}
					niter++;
				}
			}
						break;
			case FLAGGEDARRAY: {
				FlaggedArrayVariable  ar(prop_key.dataDepth);
				while (niter != gref.NodesEnd()) {
					(*niter)->Read(prop_key, ar);
					for (size_t i = 0U; i<prop_key.dataDepth; i++)
						if (ar.Flag(i) == DIRICH) {
							size_t  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
							cout << "FIX it --> PDE_Integrator_CRM 4\n"; getchar();
							//G_.ZeroRow(position);
							//G_.Add(position, position, scale_factor_);
							rh_[position] = scale_factor_ * ar(i);
						}
					niter++;
				}
			}
							   break;
			default:
				throw csmp::Exception(FATAL_ERROR,
					"PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignEssentialConditions",
					"Variable type not recognised by this method");
			}
		} // end for

	} // end AssignEssentialConditions

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::IntegrateOver(SIMPLICIAL_COMPLEX<dim>& domain, bool debug)
	{
		// 1. configure algorithm
		clock_t loc_t = clock();
		EstablishMatrixSetup(domain);
		
		
		// 2. Accumulation: Note that the conditions that pertain to the group must be input !
		Accumulate(domain);

		// 3. If the computation is transient initial conditions must be input into the righthand vector
		if (Transient() == true) { cout << "Transient in PDE_Integrator_CRM should be fixed\n"; getchar(); AssignInitialConditions(domain); }

		// 4. If the computation is transient initial conditions must be input into the righthand vector
		if (Transient() == true) { cout << "LateAccumulate in PDE_Integrator_CRM should be fixed\n"; getchar();  LateAccumulate(domain); }

		// 5. assign conditions like Dirichlet or Neumann boundary conditions etc.
		AssignEssentialConditions(domain);

		// 5.1 Remove Dirichelet Condition from Matrix, Modify RHS and Convert data structure JV (in CRM G_) to ia, ja, a
		G_.Set_Dirichelet_RHS_CRM(rh_ , dirich_);

		loc_t = clock() - loc_t;
		cout << "\n\n CLOCK Accumulate crm: \t" << loc_t << "\n\n";

		// 6. diagnostics
		if (debug) {
			OutputGlobals();
			// OutputInput();
		}

		loc_t = clock();
		// 7. invert global matrix
		x_.resize(rh_.size());
		Solve();
		
		//7.5 insert Dirichelet value in solution (map To Global)
		G_.mapToGlobal(x_, dirich_);

		loc_t = clock() - loc_t;
		cout << "\n\n CLOCK (Solve and) get Global Solution crm: \t" << loc_t << "\n\n";

		// 8. write results back into Model
		OutputResults(domain);

		// 9. Calculation of result-dependent properties                                 
		PostProcess(domain);

	} // end IntegrateOver

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::EstablishMatrixSetup(const SIMPLICIAL_COMPLEX<dim>& gref)
	{
		// -------------------------------------------------------------------
		// 0. If the algorithm is just re-used, (and has not been reset by
		//    the user, the righthand vector and
		//    the solution vector are zeroed and nothing else is done.
		// -------------------------------------------------------------------
		// * Change: x is not reset anymore
		if (setup_established_ && (!basic_operands_.empty() && !test_operands_.empty()))
		{
			if (rh_.size() > 0) fill(rh_.begin(), rh_.end(), 0.);
			if (G_.RowsJV()  > 0 && !retain_matrix_) G_.Resize(rh_.size());
		}
		// -------------------------------------------
		// 1. determine basic sizes for G, x, rh
		// -------------------------------------------
		dof_per_node_ = 0U;
		target_.nodes = gref.Nodes();
		target_.elements = gref.Elements();
		// ------------------------------------------------------
		// 2. if A.-setup for first time or if rebuild is necessary:
		//    checking basic operands and comparing MathOperators
		// ------------------------------------------------------
		const Index  unspecified;
		// lefthand MathOperators first
		// ----------------------------
		for (typename map<string, MathOperatorLHS<dim>*>::iterator
			lhs_it = lhs_operators_.begin(); lhs_it != lhs_operators_.end(); lhs_it++)
		{
			// making list of unique basic operands
			Index pkey = (*lhs_it).second->BasicOperandKey();
			if (verbose_)
				cout << "\nFor: '" << (*lhs_it).first << "' PDE operator (CRM) is added to lefthand term list." << endl;

			// checking whether the intended variables exist in the database
			if (unspecified != pkey)
				basic_operands_[(*lhs_it).second->BasicOperand()] = 0U;
			else
				throw csmp::Exception(WARNING,
					"PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
					"lefthand basic operand not found");
			// test function operands are picked up when the righthandside is accumulated
			// since they must also be present in there
		}
		// righthand MathOperators
		// -----------------------
		for (typename map<string, MathOperatorRHS<dim>*>::iterator
			rhs_it = rhs_operators_.begin(); rhs_it != rhs_operators_.end(); rhs_it++)
		{
			// making a list of unique test operands
			Index pkey = (*rhs_it).second->TestOperandKey();
			if (verbose_)
				cout << "\nFor: '" << (*rhs_it).first << "' PDE operator (CRM) is added to righthand term list." << endl;

			if (unspecified != pkey)
				test_operands_[(*rhs_it).second->TestOperand()] = 0U;
			else
				throw csmp::Exception(ERROR, "PDE_Integrator_CRM<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
					"righthand test operand not found");
		}
		// ---------------------------------------------------------------------------------------
		// 3. Testing: for each righthand operand there must be a basic or test operand on the LHS
		// ---------------------------------------------------------------------------------------
		for (typename map<string, MathOperatorLHS<dim>*>::iterator
			lhs_it = lhs_operators_.begin(); lhs_it != lhs_operators_.end(); lhs_it++)
			if (basic_operands_.find(((*lhs_it).second->BasicOperand())) == basic_operands_.end() &&
				test_operands_.find(((*lhs_it).second->TestOperand())) == test_operands_.end())
			{
				cout << "\nPDE_Integrator_CRM<" << dim << ">::EstablishMatrixSetup: ";
				cout << "There is no lefthand operand corresponding to righthand operand. ";
				cout << "\nThe system of equations is undefined. ";
				cout << "\nCreate corresponding LHS basic or test Operand for: ";
				cout << (*lhs_it).first << endl;
				throw csmp::Exception(ERROR, "PDE_Integrator_CRM<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
					"lefthand basic or test operand missing");
			}
		// ----------------------------------------------------------------
		//   4. Offsets are assigned to test Operands
		//      indicating positions, i,j in solution matrix.
		// ----------------------------------------------------------------
		size_t  offset(0U);
		for (operandsIterator iter = test_operands_.begin(); iter != test_operands_.end(); iter++)
		{
			// calculating the degrees of freedom per node
			if ((*iter).first.key.type == SCALAR)
				dof_per_node_ += 1U;
			else if ((*iter).first.key.type == VECTOR)
				dof_per_node_ += dim;
			else if ((*iter).first.key.type == TENSOR)
				dof_per_node_ += dim2_;
			else if ((*iter).first.key.type == ARRAY)
				dof_per_node_ += (*iter).first.key.dataDepth;
			else if ((*iter).first.key.type == FLAGGEDARRAY)
				dof_per_node_ += (*iter).first.key.dataDepth;

			// the matrix/righthand offsets are stored with the operands
			// offset starts out as zero.
			switch ((*iter).first.key.place)
			{
			case NODE:
				(*iter).second = offset;
				switch ((*iter).first.key.type) {
				case SCALAR: offset += target_.nodes;
					break;
				case VECTOR: offset += target_.nodes * dim;
					break;
				case TENSOR: offset += target_.nodes * dim2_;
					break;
				case ARRAY:  offset += target_.nodes * (*iter).first.key.dataDepth;
					break;
				case FLAGGEDARRAY:  offset += target_.nodes * (*iter).first.key.dataDepth;
				}
				break;
			case ELEMENT:
				cout << "\nPDE_Integrator_CRM<" << dim << ">::EstablishMatrixSetup: ";
				cout << " Test operand is an Element variable. Is this intended ?" << endl;
				(*iter).second = offset;
				switch ((*iter).first.key.type) {
				case SCALAR: offset += target_.elements;
					break;
				case VECTOR: offset += target_.elements * dim;
					break;
				case TENSOR: offset += target_.elements * dim2_;
					break;
				case ARRAY:  offset += target_.elements * (*iter).first.key.dataDepth;
					break;
				case FLAGGEDARRAY:  offset += target_.elements * (*iter).first.key.dataDepth;
				}
				break;
			default:
				throw csmp::Exception(FATAL_ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
					"Test operand placement unresolved");
			}
		}

		// --------------------------------------------------------------------------
		// 5. communicating offsets to MathOperators
		// --------------------------------------------------------------------------
		operandsIterator  iter;

		for (typename map<string, MathOperatorLHS<dim>*>::iterator
			lhs_it = lhs_operators_.begin(); lhs_it != lhs_operators_.end(); lhs_it++)
		{
			// since basic and test operand offsets must be the same, but only the test operands
			// have been assigned an offset basic operand offsets are derived from test operand offsets
			if ((iter = test_operands_.find((*lhs_it).second->BasicOperand())) != test_operands_.end()) {
				(*lhs_it).second->BasicOperandOffset((*iter).second);
				if ((iter = basic_operands_.find((*lhs_it).second->BasicOperand())) != basic_operands_.end())
					(*iter).second = (*lhs_it).second->BasicOperandOffset();
			}
			else
				throw csmp::Exception(FATAL_ERROR, "PDE_Integrator_CRM<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
					"LHS basic operand matrix placement i unresolved");

			if ((iter = test_operands_.find((*lhs_it).second->TestOperand())) != test_operands_.end())
				(*lhs_it).second->TestOperandOffset((*iter).second);
			else
				throw csmp::Exception(FATAL_ERROR, "PDE_Integrator_CRM<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
					"LHS test function operand matrix placement j unresolved");
		}
		for (typename map<string, MathOperatorRHS<dim>*>::iterator
			rhs_it = rhs_operators_.begin(); rhs_it != rhs_operators_.end(); rhs_it++)
		{
			if ((iter = test_operands_.find((*rhs_it).second->TestOperand())) != test_operands_.end())
				(*rhs_it).second->TestOperandOffset((*iter).second);
			else
				throw csmp::Exception(ERROR, "PDE_Integrator_CRM<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
				(*rhs_it).first.c_str(), "RHS operand vector^T placement i unresolved...");
		}

		// ----------------------------------------------------------------------
		// 6. Resizing 'G' and righthand vector 'rh' which is initialised to zero
		// ----------------------------------------------------------------------
		// * Change: resize x as well
		G_.Resize(offset);
		dirich_.reserve(gref.PerimeterNodes());
		rh_.resize(offset);
		vector<double64>(rh_).swap(rh_);
		fill(rh_.begin(), rh_.end(), 0.);
		x_.reserve(gref.InteriorNodes());
		vector<double64>(x_).swap(x_);
		setup_established_ = true;

	} // end EstablishMatrixSetup()

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void  PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::OutputGlobals(int32 precision)
	{}

	/*
	******************************************************************************
	PAERT 2: The following functions are copied from coressponding PDE_Integrator
	******************************************************************************
	*/

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void  PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::OutputResults(SIMPLICIAL_COMPLEX<dim>& gref)
	{
		Index   prop_key;
		size_t  offset;

		for (operandsIterator
			it = basic_operands_.begin(); it != basic_operands_.end(); it++)
		{
			typename vector<Node<dim>*>::iterator  gfirst(gref.NodesBegin());
			prop_key = (*it).first.key;
			offset = (*it).second;

			if (prop_key.place != NODE)
				throw csmp::Exception(ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::OutputResults(Model)",
					"only nodal properties can be output by this method.");
			switch (prop_key.type)
			{
			case SCALAR:
				while (gfirst != gref.NodesEnd()) {
					const double64 sc = x_[(*gfirst)->Idx() + offset];
					(*gfirst)->Store(prop_key, makeScalar((*gfirst)->Status(prop_key), sc));
					gfirst++;
				}
				break;
			case VECTOR: {
				VectorVariable<dim>  vc;
				while (gfirst != gref.NodesEnd()) {
					(*gfirst)->Read(prop_key, vc);
					for (size_t i = 0U; i<dim; i++) vc(i) = x_[(*gfirst)->Idx() * dim + i + offset];
					(*gfirst)->Store(prop_key, vc);
					gfirst++;
				}
			}
						 break;
			case TENSOR: {
				TensorVariable<dim>  ts;
				while (gfirst != gref.NodesEnd()) {
					(*gfirst)->Read(prop_key, ts);
					for (size_t i = 0U; i<dim; i++)
						for (size_t k = 0U; k<dim; k++)
							ts(i, k) = x_[(*gfirst)->Idx() * dim2_ + i * dim + k + offset];
					(*gfirst)->Store(prop_key, ts);
					gfirst++;
				}
			}
						 break;
			case ARRAY: {
				ArrayVariable  ar(prop_key.dataDepth);
				while (gfirst != gref.NodesEnd()) {
					(*gfirst)->Read(prop_key, ar);
					for (size_t i = 0U; i<prop_key.dataDepth; i++) ar(i) = x_[(*gfirst)->Idx() * prop_key.dataDepth + i + offset];
					(*gfirst)->Store(prop_key, ar);
					gfirst++;
				}
			}
						break;
			case FLAGGEDARRAY: {
				FlaggedArrayVariable  ar(prop_key.dataDepth);
				while (gfirst != gref.NodesEnd()) {
					(*gfirst)->Read(prop_key, ar);
					for (size_t i = 0U; i<prop_key.dataDepth; i++) ar(i) = x_[(*gfirst)->Idx() * prop_key.dataDepth + i + offset];
					(*gfirst)->Store(prop_key, ar);
					gfirst++;
				}
			}
							   break;
			default:
				throw csmp::Exception(ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::OutputResults(Model)",
					"Output to ARRAY type variables is not supported by this method yet.");

			} // end switch(type)

		} // end for

	} // end OutputResults

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void  PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::LateAccumulate(const SIMPLICIAL_COMPLEX<dim>& gref)
	{
		// accumulating as late addition into the righhand vector 'rhs'
		// ------------------------------------------------------------
		for (typename map<string, MathOperatorRHS<dim>*>::const_iterator
			it_rhs = rhs_operators_.begin(); it_rhs != rhs_operators_.end(); it_rhs++)
			if ((*it_rhs).second->AddLater() || (*it_rhs).second->SubtractLater())
				for (typename vector<typename SIMPLICIAL_COMPLEX<dim>::Simplex*>::const_iterator
					git = gref.ElementsBegin(); git != gref.ElementsEnd(); git++)
				{
					(*it_rhs).second->GetOperands(*(*git));
					(*it_rhs).second->ComputeContribution(*(*git));
					if ((*it_rhs).second->MultiplyWithTimeIncrement())
						(*it_rhs).second->MultiplyWithTimeFactor(time_increment_);
					if ((*it_rhs).second->DivideByTimeIncrement())
						(*it_rhs).second->MultiplyWithTimeFactor(1. / time_increment_);
					(*it_rhs).second->AssignToGlobal(*(*git), rh_);
				}
	} // end Late Accumulate

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void  PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::PostProcess(const SIMPLICIAL_COMPLEX<dim>& gref)

	{
		if (postpro_operators_.empty()) return;

		for (typename map<string, MathOperatorLHS<dim>*>::iterator
			it = postpro_operators_.begin(); it != postpro_operators_.end(); it++)
		{
			for (size_t i = 1; i <= (*it).second->ApplicationCycles(); i++)
			{
				// setting application cylce such that it can be used by PDE operator
				(*it).second->ApplicationCycle(i);
				if (verbose_) cout << "\nPDE_Integrator<" << dim;
				//              if (verbose_) cout <<">::PostProcess: Computing: "<< (*it).first <<" in region'"<< gref.Name() <<"'\n";
				if (verbose_) cout << ">::PostProcess: Computing: " << (*it).first << "\n";
				for (typename vector<typename SIMPLICIAL_COMPLEX<dim>::Simplex*>::const_iterator
					git = gref.ElementsBegin(); git != gref.ElementsEnd(); git++)
				{
					(*it).second->GetOperands(*(*git));
					(*it).second->ComputeContribution(*(*git));
					(*it).second->WriteOperands(*(*git));
				}
			}
		}
	} // end PostProcess

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::AssignInitialConditions(const SIMPLICIAL_COMPLEX<dim>& gref)

	{
		size_t                  i, j;
		size_t                  position, offset;
		ScalarVariable          sc;
		VectorVariable<dim>     vc;
		TensorVariable<dim>     ts;
		Index                   prop_key;

		if (!setup_established_)
			throw csmp::Exception(ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignInitialConditions",
				"please call EstablishMatrixSetup() prior to this method.");

		if (basic_operands_.empty()) {
			throw csmp::Exception(ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignInitialConditions",
				"No basic operands have been specified...");
			return;
		}

		for (operandsConstIterator
			it = test_operands_.begin(); it != test_operands_.end(); it++)
		{
			prop_key = (*it).first.key;
			offset = (*it).second;
			typename vector<csmp::Node<dim>*>::const_iterator  niter(gref.NodesBegin());

			if (prop_key.place != NODE) {
				throw csmp::Exception(WARNING, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignInitialConditions",
					"So far no conditions are assigned to elements, faces, segments");
				return;
			}

			switch (prop_key.type)
			{
			case SCALAR:
				while (niter != gref.NodesEnd()) {
					// here the number in the new list is required to facilitate input
					// into the size-restricted computation matrix
					(*niter)->Read(prop_key, sc);
					position = (*niter)->Idx() + offset;
					rh_[position] *= sc();
					niter++;
				}
				break;
			case VECTOR:
				while (niter != gref.NodesEnd()) {
					(*niter)->Read(prop_key, vc);
					for (i = 0; i<dim; i++) {
						position = (*niter)->Idx() * dim + i + offset;
						rh_[position] *= vc(i);
					}
					niter++;
				}
				break;
			case TENSOR:
				while (niter != gref.NodesEnd()) {
					(*niter)->Read(prop_key, ts);
					for (i = 0; i<dim; i++)
						for (j = 0; j<dim; j++) {
							position = (*niter)->Idx() * dim2_ + i * dim + j + offset;
							rh_[position] *= ts(i, j);
						}
					niter++;
				}
				break;
			case ARRAY: {
				ArrayVariable  ar(prop_key.dataDepth);
				while (niter != gref.NodesEnd()) {
					(*niter)->Read(prop_key, ar);
					for (i = 0; i<prop_key.dataDepth; i++) {
						position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
						rh_[position] *= ar(i);
					}
					niter++;
				}
			}
						break;
			case FLAGGEDARRAY: {
				FlaggedArrayVariable far(prop_key.dataDepth);
				while (niter != gref.NodesEnd()) {
					(*niter)->Read(prop_key, far);
					for (i = 0; i<prop_key.dataDepth; i++) {
						position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
						rh_[position] *= far(i);
					}
					niter++;
				}
			}
			}
		}

	} // end AssignInitialConditions

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::Accumulate(const SIMPLICIAL_COMPLEX<dim>& gref)
	{
		// setting up the index mapping from global to local node ID numbers
		gref.RenumberNodes();
		// accumulating into the sparse matrix 'G'
		// ---------------------------------------
		for (typename map<string, MathOperatorLHS<dim>*>::iterator
			it_lhs = lhs_operators_.begin(); it_lhs != lhs_operators_.end(); it_lhs++)
			if (!(*it_lhs).second->AddLater() && !(*it_lhs).second->SubtractLater())
				for (typename vector<typename SIMPLICIAL_COMPLEX<dim>::Simplex*>::const_iterator
					git = gref.ElementsBegin(); git != gref.ElementsEnd(); git++)
				{
					(*it_lhs).second->GetOperands(*(*git));
					(*it_lhs).second->ComputeContribution(*(*git));
					if ((*it_lhs).second->MultiplyWithTimeIncrement())
						(*it_lhs).second->MultiplyWithTimeFactor(time_increment_);
					if ((*it_lhs).second->DivideByTimeIncrement())
						(*it_lhs).second->MultiplyWithTimeFactor(1. / time_increment_);
					(*it_lhs).second->AssignToGlobal(*(*git), G_);
				}

		// accumulating via multiplication into the righthand vector 'rhs'
		// --------------------------------------------------------------
		for (typename map<string, MathOperatorRHS<dim>*>::iterator
			it_rhs = rhs_operators_.begin(); it_rhs != rhs_operators_.end(); it_rhs++)
			if (!(*it_rhs).second->AddLater() && !(*it_rhs).second->SubtractLater())
				for (typename vector<typename SIMPLICIAL_COMPLEX<dim>::Simplex*>::const_iterator
					git = gref.ElementsBegin(); git != gref.ElementsEnd(); git++)
				{
					(*it_rhs).second->GetOperands(*(*git));
					(*it_rhs).second->ComputeContribution(*(*git));
					if ((*it_rhs).second->MultiplyWithTimeIncrement())
						(*it_rhs).second->MultiplyWithTimeFactor(time_increment_);
					if ((*it_rhs).second->DivideByTimeIncrement())
						(*it_rhs).second->MultiplyWithTimeFactor(1. / time_increment_);
					(*it_rhs).second->AssignToGlobal(*(*git), rh_);
				}
	} // end Accumulate

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::ScaleEssentialConditions(double64 scale_factor)
	{
		scale_factor_ = scale_factor;

	} // end ScaleEssentialConditions

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void  PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::FirstGuess(const vector<double64>& guess) {
		assert(guess.size() == x_.size());
		copy(guess.begin(), guess.end(), x_.begin());
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void  PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::Reset(bool delete_math_operators)
	{
		if (delete_math_operators) {
			lhs_operators_.erase(lhs_operators_.begin(), lhs_operators_.end());
			rhs_operators_.erase(rhs_operators_.begin(), rhs_operators_.end());
			basic_operands_.erase(basic_operands_.begin(), basic_operands_.end());
			test_operands_.erase(test_operands_.begin(), test_operands_.end());
			postpro_operators_.erase(postpro_operators_.begin(), postpro_operators_.end());
			setup_established_ = false;
		}

		// restoring defaults
		time_increment_ = 0.;
		retain_matrix_ = false;
		target_.nodes = target_.elements = 0U;
	} // end Reset

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::Solve()
	{
		solver_->Solve(G_, rh_, x_, dof_per_node_);
	} // end Solve

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void  PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::SolutionVector(vector<double64>& sol) const {
		sol.resize(x_.size());
		vector<double64>(sol).swap(sol);
		copy(x_.begin(), x_.end(), sol.begin());
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::ListMathOperatorsLHS() const
	{
		typename map<string, MathOperatorLHS<dim>*>::const_iterator it;
		for (it = lhs_operators_.begin(); it != lhs_operators_.end(); it++)
		{ cout << (*it).first << ":  "; (*it).second->Out(); }
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::ListMathOperatorsRHS() const
	{
		typename map<string, MathOperatorRHS<dim>*>::const_iterator it;
		for (it = rhs_operators_.begin(); it != rhs_operators_.end(); it++)
		{ cout << (*it).first << ":  ";	(*it).second->Out(); }
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::Add(MathOperatorLHS<dim>* op)
	{
		// The name for the algorithm is combined out of its operands
		lhs_operators_[op->Name()] = op;
		// force update during next application
		setup_established_ = false;
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::Add(MathOperatorRHS<dim>* op)
	{
		rhs_operators_[op->Name()] = op;
		// force update during next application
		setup_established_ = false;
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::AddPostProcess(MathOperatorLHS<dim>* op)
	{
		postpro_operators_[op->Name()] = op;
		// force update during next application
		setup_established_ = false;
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::~PDE_Integrator_CRM() { if (newed_Solver_object) delete solver_; }

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::RetainGlobalSolutionMatrix(bool retain) { retain_matrix_ = retain; }

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::SetSolver(Solver* new_solver) {
		assert(new_solver != NULL);
		if (solver_ != new_solver) {
			delete solver_;
			solver_ = new_solver;
		}
		newed_Solver_object = false;
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	Solver* PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::GetSolver() const { return solver_; }

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::AdjustSolverSettings()
	{
		ErrorHandler&  csmp_error(ErrorHandler::Instance());
		csmp_error.notice(WARNING, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AdjustSolverSettings",
			"this call to the method did nothing, the method is defined only in the subclasses.");
	} // end AdjustSolverSettings

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	bool  PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::Transient() const { return !(time_increment_ < numeric_limits<double64>::epsilon()); }

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	void   PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::TimeIncrement(double64 dt) { time_increment_ = dt; }

	// Constructors
	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::PDE_Integrator_CRM()
		:
#ifdef CSMP_WITH_SAMG_SOLVER
		solver_(new SAMG_Solver()),
#else
		/// add extra functionality for alternative solver if needed
		solver_(new CSMP_DEFAULT_LINEAR_SOLVER()),
#endif
		dim2_(dim*dim),
		dof_per_node_(0),
		setup_established_(false),
		retain_matrix_(false),
		newed_Solver_object(true),
		time_increment_(0.),
		scale_factor_(1.),
		verbose_(true)
	{
		target_.nodes = target_.elements = 0U;
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::PDE_Integrator_CRM(Solver& solver)
		: solver_(&solver),
		dim2_(dim*dim),
		dof_per_node_(0),
		setup_established_(false),
		retain_matrix_(false),
		newed_Solver_object(false),
		time_increment_(0.),
		scale_factor_(1.),
		verbose_(true)
	{
		target_.nodes = target_.elements = 0U;
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::PDE_Integrator_CRM(Solver* solver)
		: solver_(solver),
		dim2_(dim*dim),
		dof_per_node_(0),
		setup_established_(false),
		retain_matrix_(false),
		newed_Solver_object(false),
		time_increment_(0.),
		scale_factor_(1.),
		verbose_(true)
	{
		assert(solver != NULL); target_.nodes = target_.elements = 0U;
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::PDE_Integrator_CRM(const PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>& a)
		:
#ifdef CSMP_WITH_SAMG_SOLVER
		solver_(new SAMG_Solver()),
#else
		/// add extra functionality for alternative solver if needed
		solver_(new CSMP_DEFAULT_LINEAR_SOLVER()),
#endif
		rh_(a.rh_),
		x_(a.x_),
		G_(a.G_),
		lhs_operators_(a.lhs_operators_),
		rhs_operators_(a.rhs_operators_),
		basic_operands_(a.basic_operands_),
		test_operands_(a.test_operands_),
		postpro_operators_(a.postpro_operators_),
		dim2_(dim*dim),
		dof_per_node_(a.dof_per_node_),
		setup_established_(a.setup_established_),
		retain_matrix_(a.retain_matrix_),
		newed_Solver_object(true),
		scale_factor_(1.),
		time_increment_(a.time_increment_),
		target_(a.target_),
		verbose_(a.verbose_)
	{
		assert(solver_ != NULL);
		cout << "\nPDE_Integrator: copy constructor: ";
		cout << "New algorithm uses new instance of solver." << endl;
	}

	template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
	PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>& PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>::operator=(const PDE_Integrator_CRM<dim, SIMPLICIAL_COMPLEX>& a)
	{
		if (&a != this) {
			G_ = a.G_;
			rh_ = a.rh_;
			x_ = a.x_;

			if (a.solver_ != NULL) {
				assert(solver_ != NULL or solver_ == 0);
				if (newed_Solver_object) {
					delete solver_;
					/// SKM fix:  this is not a clean solution. One should bring over solver from other integrator instance
#ifdef CSMP_WITH_SAMG_SOLVER
					solver_ = new SAMG_Solver();
#else
					/// add extra functionality for alternative solver if needed
					solver_ = new CSMP_DEFAULT_LINEAR_SOLVER();
#endif
				}
				else solver_ = a.solver_;
			}

			lhs_operators_ = a.lhs_operators_;
			rhs_operators_ = a.rhs_operators_;
			postpro_operators_ = a.postpro_operators_;
			basic_operands_ = a.basic_operands_;
			test_operands_ = a.test_operands_;
			time_increment_ = a.time_increment_;
			retain_matrix_ = a.retain_matrix_;
			setup_established_ = a.setup_established_;
			target_ = a.target_;
			dof_per_node_ = a.dof_per_node_;
			scale_factor_ = a.scale_factor_;
			verbose_ = a.verbose_;
		}
		return *this;
	}
	// end of Constructors

	template class PDE_Integrator_CRM<1U, Region>;
	template class PDE_Integrator_CRM<2U, Region>;
	template class PDE_Integrator_CRM<3U, Region>;
	
	template class PDE_Integrator_CRM<1U, Boundary>;
	template class PDE_Integrator_CRM<2U, Boundary>;
	template class PDE_Integrator_CRM<3U, Boundary>;

	template class PDE_Integrator_CRM<1U, SplitBoundary>;
	template class PDE_Integrator_CRM<2U, SplitBoundary>;
	template class PDE_Integrator_CRM<3U, SplitBoundary>;
	
} // end namespace csmp
