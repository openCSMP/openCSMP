#include "PDE_Integrator_UoM.h"
#include "NimbleRegion.h"

#if defined(_OPENMP )
#include "omp.h"
#include "FiniteElementManager.h"
#endif

using namespace std;

namespace csmp {

template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::PDE_Integrator_UoM( Solver& solver )
 : PDE_Integrator<dim,COMPUTATION_DOMAIN>(&solver)
{
}



template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::~PDE_Integrator_UoM()
 {
    // are we sure that any newed solver object in the base class is deleted ?
 }


template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
void PDE_Integrator_UoM<dim, COMPUTATION_DOMAIN>::EstablishMatrixSetup(const COMPUTATION_DOMAIN<dim>& gref)
{
    // -------------------------------------------------------------------
    // 0. If the algorithm is just re-used, (and has not been reset by
    //    the user, the righthand vector and
    //    the solution vector are zeroed and nothing else is done.
    // -------------------------------------------------------------------
    // * Change: x is not reset anymore
    if (this->setup_established_ &&
      (!this->basic_operands_.empty() && !this->test_operands_.empty()))
      {
         if (this->rh_.size() > 0) fill( this->rh_.begin(), this->rh_.end(), 0.);
#if defined(_OPENMP )
         for (size_t tid = 0; tid < omp_get_max_threads(); tid++) {
           if (thread_rh_[tid].size() > 0) fill(thread_rh_[tid].begin(), thread_rh_[tid].end(), 0.);
      }
#endif
      if (this->G_.Rows() > 0 && !this->retain_matrix_) {
         this->G_.Erase();
         this->G_.Resize(this->rh_.size());
#if defined(_OPENMP )
         for (size_t tid = 0; tid < omp_get_max_threads(); tid++) {
              this->thread_G_[tid].Erase();
              this->thread_G_[tid].Resize(rh_.size());
           }
#endif
        }
    }

    // -------------------------------------------
    // 1. determine basic sizes for G, x, rh
    // -------------------------------------------
    this->dof_per_node_ = 0U;
    this->target_.nodes = gref.Nodes();
    this->target_.elements = gref.Elements();


    // ------------------------------------------------------
    // 2. if A.-setup for first time or if rebuild is necessary:
    //    checking basic operands and comparing MathOperators
    // ------------------------------------------------------
    const Index  unspecified;

    // lefthand MathOperators first
    // ----------------------------
    for ( typename map<string, MathOperatorLHS<dim>*>::iterator
          lhs_it = this->lhs_operators_.begin(); lhs_it != this->lhs_operators_.end(); lhs_it++ )
      {
      // making list of unique basic operands
      Index pkey = (*lhs_it).second->BasicOperandKey();
      if (verbose_)
        cout << "\nFor: '" << (*lhs_it).first << "' PDE operator is added to lefthand term list." << endl;

      // checking whether the intended variables exist in the database
      if (unspecified != pkey)
        this->basic_operands_[(*lhs_it).second->BasicOperand()] = 0U;
      else
        throw csmp::Exception(WARNING,
          "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
          "lefthand basic operand not found");
      // test function operands are picked up when the righthandside is accumulated
      // since they must also be present in there
    }

    // righthand MathOperators
    // -----------------------
    for ( typename map<string, MathOperatorRHS<dim>*>::iterator
          rhs_it = this->rhs_operators_.begin(); rhs_it != this->rhs_operators_.end(); rhs_it++ )
      {
        // making a list of unique test operands
        Index pkey = (*rhs_it).second->TestOperandKey();
        if (verbose_)
          cout << "\nFor: '" << (*rhs_it).first << "' PDE operator is added to righthand term list." << endl;

        if (unspecified != pkey)
          this->test_operands_[(*rhs_it).second->TestOperand()] = 0U;
        else
          throw csmp::Exception(ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
            "righthand test operand not found");
      }

    // ---------------------------------------------------------------------------------------
    // 3. Testing: for each righthand operand there must be a basic or test operand on the LHS
    // ---------------------------------------------------------------------------------------
    for ( typename map<string, MathOperatorLHS<dim>*>::iterator
          lhs_it = this->lhs_operators_.begin(); lhs_it != this->lhs_operators_.end(); lhs_it++ )
      if ( this->basic_operands_.find(((*lhs_it).second->BasicOperand())) == this->basic_operands_.end() &&
           this->test_operands_.find(((*lhs_it).second->TestOperand())) == this->test_operands_.end())
        {
          cout << "\nPDE_Integrator<" << dim << ">::EstablishMatrixSetup: ";
          cout << "There is no lefthand operand corresponding to righthand operand. ";
          cout << "\nThe system of equations is undefined. ";
          cout << "\nCreate corresponding LHS basic or test Operand for: ";
          cout << (*lhs_it).first << endl;
          throw csmp::Exception(ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
            "lefthand basic or test operand missing");
        }

    // ----------------------------------------------------------------
    //   4. Offsets are assigned to test Operands
    //      indicating positions, i,j in solution matrix.
    // ----------------------------------------------------------------
    size_t  offset(0U);
    for ( typename PDE_Integrator<dim,COMPUTATION_DOMAIN>::operandsIterator
          iter = this->test_operands_.begin(); iter != this->test_operands_.end(); iter++ )
      {
        // calculating the degrees of freedom per node
        if ((*iter).first.key.type == SCALAR)
          this->dof_per_node_ += 1U;
        else if ((*iter).first.key.type == VECTOR)
          this->dof_per_node_ += dim;
        else if ((*iter).first.key.type == TENSOR)
          this->dof_per_node_ += this->dim2_;
        else if ((*iter).first.key.type == ARRAY)
          this->dof_per_node_ += (*iter).first.key.dataDepth;
        else if ((*iter).first.key.type == FLAGGEDARRAY)
          this->dof_per_node_ += (*iter).first.key.dataDepth;

        // the matrix/righthand offsets are stored with the operands
        // offset starts out as zero.
        switch ((*iter).first.key.place)
          {
            case NODE:
              (*iter).second = offset;
              switch ((*iter).first.key.type) {
                  case SCALAR: offset += this->target_.nodes;
                    break;
                  case VECTOR: offset += this->target_.nodes * dim;
                    break;
                  case TENSOR: offset += this->target_.nodes * this->dim2_;
                    break;
                  case ARRAY:  offset += this->target_.nodes * (*iter).first.key.dataDepth;
                    break;
                  case FLAGGEDARRAY:  offset += this->target_.nodes * (*iter).first.key.dataDepth;
                }
              break;
            case ELEMENT:
              cout << "\nPDE_Integrator<" << dim << ">::EstablishMatrixSetup: ";
              cout << " Test operand is an Element variable. Is this intended ?" << endl;
              (*iter).second = offset;
              switch ((*iter).first.key.type) {
                case SCALAR: offset += this->target_.elements;
                  break;
                case VECTOR: offset += this->target_.elements * dim;
                  break;
                case TENSOR: offset += this->target_.elements * this->dim2_;
                  break;
                case ARRAY:  offset += this->target_.elements * (*iter).first.key.dataDepth;
                  break;
                case FLAGGEDARRAY:  offset += this->target_.elements * (*iter).first.key.dataDepth;
                }
              break;
            default:
              throw csmp::Exception(FATAL_ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                "Test operand placement unresolved");
          }
     }

   // --------------------------------------------------------------------------
   // 5. communicating offsets to MathOperators
   // --------------------------------------------------------------------------
   typename PDE_Integrator<dim,COMPUTATION_DOMAIN>::operandsIterator  iter;

   for ( typename map<string,MathOperatorLHS<dim>*>::iterator
         lhs_it=this->lhs_operators_.begin(); lhs_it!=this->lhs_operators_.end(); lhs_it++ )
     {
       // since basic and test operand offsets must be the same, but only the test operands
       // have been assigned an offset basic operand offsets are derived from test operand offsets
       if ( (iter=this->test_operands_.find((*lhs_it).second->BasicOperand())) != this->test_operands_.end() ) {
           (*lhs_it).second->BasicOperandOffset( (*iter).second );
           if ( (iter=this->basic_operands_.find((*lhs_it).second->BasicOperand())) != this->basic_operands_.end() )
               (*iter).second = (*lhs_it).second->BasicOperandOffset();
       }
       else
           throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                                  "LHS basic operand matrix placement i unresolved");

        if ( (iter=this->test_operands_.find((*lhs_it).second->TestOperand())) != this->test_operands_.end() )
          (*lhs_it).second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                                       "LHS test function operand matrix placement j unresolved");
     }
   for ( typename map<string,MathOperatorRHS<dim>*>::iterator
         rhs_it=this->rhs_operators_.begin(); rhs_it!=this->rhs_operators_.end(); rhs_it++ )
     {
        if ( (iter=this->test_operands_.find((*rhs_it).second->TestOperand())) != this->test_operands_.end() )
          (*rhs_it).second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                          (*rhs_it).first.c_str(), "RHS operand vector^T placement i unresolved...");
     }

   // ----------------------------------------------------------------------
   // 6. Resizing 'G' and righthand vector 'rh' which is initialised to zero
   // ----------------------------------------------------------------------
   this->G_.Resize( offset );
   this->rh_.resize( offset );
   vector<double64>( this->rh_ ).swap( this->rh_ );
   fill( this->rh_.begin(), this->rh_.end(), 0. );
   this->x_.resize( offset );
   vector<double64>( this->x_ ).swap( this->x_ );
#if defined(_OPENMP )
   for (size_t tid = 0 ; tid < omp_get_max_threads() ; tid++){
       this->thread_G_[tid].Resize(offset);
       this->thread_rh_[tid].resize(offset);
       vector<double64>( thread_rh_[tid] ).swap( thread_rh_[tid] );
       fill( thread_rh_[tid].begin(), thread_rh_[tid].end(), 0. );
   }
#endif
   this->setup_established_ = true;

  } // end EstablishMatrixSetup()




/**
    Needed to overwrite base clase method because DOF_indexes_ yield different positions in matrices and right-hand vector.
*/
template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
void PDE_Integrator_UoM<dim, COMPUTATION_DOMAIN>::AssignInitialConditions(const COMPUTATION_DOMAIN<dim>& gref)
  {
    size_t                  i, j;
    size_t                  position, offset;
    ScalarVariable          sc;
    VectorVariable<dim>     vc;
    TensorVariable<dim>     ts;
    Index                   prop_key;

    if (!this->setup_established_)
      throw csmp::Exception(ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::AssignInitialConditions",
                                   "please call EstablishMatrixSetup() prior to this method.");

    if (this->basic_operands_.empty() )
      throw csmp::Exception(ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::AssignInitialConditions",
                                   "No basic operands have been specified...");

    for ( typename PDE_Integrator<dim,COMPUTATION_DOMAIN>::operandsConstIterator
          it = this->test_operands_.begin(); it != this->test_operands_.end(); it++ )
      {
        prop_key = (*it).first.key;
        offset = (*it).second;
        typename vector<csmp::Node<dim>*>::const_iterator  niter(gref.NodesBegin());

        if (prop_key.place != NODE )
          throw csmp::Exception( WARNING, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::AssignInitialConditions",
                                "So far no conditions are assigned to elements, faces, segments");

        switch (prop_key.type)
          {
          case SCALAR:
            while (niter != gref.NodesEnd()) {
              // here the number in the new list is required to facilitate input
              // into the size-restricted computation matrix
              (*niter)->Read(prop_key, sc);
              position = (*niter)->Idx() + offset;
              position = DOF_indexes_[position];
              if (position != NULL_IDX) {
                this->rh_[position] *= sc();
              }
              
              niter++;
            }
            break;
          case VECTOR:
            while (niter != gref.NodesEnd()) {
              (*niter)->Read(prop_key, vc);
              for (i = 0; i < dim; i++) {
                position = (*niter)->Idx() * dim + i + offset;
                position = DOF_indexes_[position];
                if (position != NULL_IDX) {
                  this->rh_[position] *= vc(i);
                }
              }
              niter++;
            }
            break;
          case TENSOR:
            while (niter != gref.NodesEnd()) {
              (*niter)->Read(prop_key, ts);
              for (i = 0; i < dim; i++)
                for (j = 0; j < dim; j++) {
                  position = (*niter)->Idx() * this->dim2_ + i * dim + j + offset;
                  position = DOF_indexes_[position];
                  if (position != NULL_IDX) {
                    this->rh_[position] *= ts(i, j);
                  }

                }
              niter++;
            }
            break;
          case ARRAY: {
            ArrayVariable  ar(prop_key.dataDepth);
            while (niter != gref.NodesEnd()) {
              (*niter)->Read(prop_key, ar);
              for (i = 0; i < prop_key.dataDepth; i++) {
                position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                position = DOF_indexes_[position];
                if (position != NULL_IDX) {
                  this->rh_[position] *= ar(i);
                }
              }
              niter++;
            }
          }
            break;
          case FLAGGEDARRAY: {
            FlaggedArrayVariable far(prop_key.dataDepth);
            while (niter != gref.NodesEnd()) {
              (*niter)->Read(prop_key, far);
              for (i = 0; i < prop_key.dataDepth; i++) {
                position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                position = DOF_indexes_[position];
                if (position != NULL_IDX) {
                  this->rh_[position] *= far(i);
                }
              }
              niter++;
            }
          }
      }
    }

  } // end AssignInitialConditions
  
  
  
  
  
template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
void PDE_Integrator_UoM<dim, COMPUTATION_DOMAIN>::Accumulate(const COMPUTATION_DOMAIN<dim>& gref)
  {
    // setting up the index mapping from global to local node ID numbers
#if !defined(_OPENMP)
    // accumulating into the compressed row matrix 'G'
    // ---------------------------------------
    for (typename map<string, MathOperatorLHS<dim>*>::iterator
      it_lhs = this->lhs_operators_.begin(); it_lhs != this->lhs_operators_.end(); it_lhs++)
      if (!(*it_lhs).second->AddLater() && !(*it_lhs).second->SubtractLater())
        for (typename vector<typename COMPUTATION_DOMAIN<dim>::CellType*>::const_iterator
          git = gref.ElementsBegin(); git != gref.ElementsEnd(); git++)
        {
          (*it_lhs).second->GetOperands(*(*git));
          (*it_lhs).second->ComputeContribution(*(*git));
          if ((*it_lhs).second->MultiplyWithTimeIncrement())
            (*it_lhs).second->MultiplyWithTimeFactor( this->time_increment_);
          if ((*it_lhs).second->DivideByTimeIncrement())
            (*it_lhs).second->MultiplyWithTimeFactor( 1. / this->time_increment_);
          //  (*it_lhs).second->AssignToGlobal(*(*git), G_);
          (*it_lhs).second->AssignToGlobal(*(*git), this->G_, pivotVector_, DOF_indexes_); // luat changed here

        }

    // accumulating via multiplication into the righthand vector 'rhs'
    // --------------------------------------------------------------
    for (typename map<string, MathOperatorRHS<dim>*>::iterator
      it_rhs = this->rhs_operators_.begin(); it_rhs != this->rhs_operators_.end(); it_rhs++)
      if (!(*it_rhs).second->AddLater() && !(*it_rhs).second->SubtractLater())
        for (typename vector<typename COMPUTATION_DOMAIN<dim>::CellType*>::const_iterator
          git = gref.ElementsBegin(); git != gref.ElementsEnd(); git++)
        {
          (*it_rhs).second->GetOperands(*(*git));
          (*it_rhs).second->ComputeContribution(*(*git));
          if ((*it_rhs).second->MultiplyWithTimeIncrement())
            (*it_rhs).second->MultiplyWithTimeFactor(this->time_increment_);
          if ((*it_rhs).second->DivideByTimeIncrement())
            (*it_rhs).second->MultiplyWithTimeFactor(1. / this->time_increment_);
          //(*it_rhs).second->AssignToGlobal(*(*git), rh_);
          (*it_rhs).second->AssignToGlobal(*(*git), this->rh_, DOF_indexes_); // luat changed here
        }
#else
    // accumulating into the compressed row matrix 'G' for each thread
    // ---------------------------------------
#pragma omp parallel
    {
      FiniteElement* fe_tmp;
      size_t tid = omp_get_thread_num();
      for (typename map<string, MathOperatorLHS<dim>*>::iterator
        it_lhs = this->thread_lhs_operators_[tid].begin(); it_lhs != this->thread_lhs_operators_[tid].end(); it_lhs++) {
        if (!(*it_lhs).second->AddLater() && !(*it_lhs).second->SubtractLater()) {
#pragma omp for
          for (int32 e = 0; e < gref.Elements(); e++)
          {
            //                    cout<<"element: "<<e<<endl;
            typename COMPUTATION_DOMAIN<dim>::CellType* eit = gref.E(e);
            fe_tmp = eit->FE(); //save old pointer.
                                // change pointer here
            eit->Assign(femgrs_[tid].E(eit->FE_Type()));

            // do needed work.
            (*it_lhs).second->GetOperands((*eit));
            (*it_lhs).second->ComputeContribution((*eit));
            if ((*it_lhs).second->MultiplyWithTimeIncrement())
              (*it_lhs).second->MultiplyWithTimeFactor(time_increment_);
            if ((*it_lhs).second->DivideByTimeIncrement())
              (*it_lhs).second->MultiplyWithTimeFactor(1. / time_increment_);
            (*it_lhs).second->AssignToGlobal((*eit), thread_G_[tid]);

            // put it back here
            eit->Assign(fe_tmp);
          }
        }
      }
    }


#pragma omp parallel
    {
      FiniteElement* fe_tmp;
      size_t tid = omp_get_thread_num();
      // accumulating via multiplication into the righhand vector 'rhs' for each thread
      // --------------------------------------------------------------
      for (typename map<string, MathOperatorRHS<dim>*>::iterator
        it_rhs = this->thread_rhs_operators_[tid].begin(); it_rhs != this->thread_rhs_operators_[tid].end(); it_rhs++) {
        if (!(*it_rhs).second->AddLater() && !(*it_rhs).second->SubtractLater()) {
#pragma omp for
          for (int32 e = 0; e < gref.Elements(); e++)
          {
            typename COMPUTATION_DOMAIN<dim>::CellType* eit = gref.E(e);
            fe_tmp = eit->FE(); //save old pointer.
                                // change pointer here
            eit->Assign(femgrs_[tid].E(eit->FE_Type()));

            // do needed work.
            (*it_rhs).second->GetOperands((*eit));
            (*it_rhs).second->ComputeContribution((*eit));
            if ((*it_rhs).second->MultiplyWithTimeIncrement())
              (*it_rhs).second->MultiplyWithTimeFactor(time_increment_);
            if ((*it_rhs).second->DivideByTimeIncrement())
              (*it_rhs).second->MultiplyWithTimeFactor(1. / time_increment_);

            (*it_rhs).second->AssignToGlobal((*eit), thread_rh_[tid]);

            // put it back here
            eit->Assign(fe_tmp);
          }
        }
      }
    }

    for (size_t tid = 0; tid < omp_get_max_threads(); tid++)
      this->G_ += this->thread_G_[tid];

    for (size_t tid = 0; tid < omp_get_max_threads(); tid++)
      for (size_t i = 0; i < this->rh_.size(); i++)
        this->rh_[i] = this->rh_[i] + this->thread_rh_[tid][i];

    //    G_.Out();
    //    for (size_t i = 0; i<rh_.size();i++)
    //        cout<<"rh:["<<i<<"]: "<<rh_[i]<<endl;
    //    exit(1);

#endif
  } // end Accumulate

  template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
  void PDE_Integrator_UoM<dim, COMPUTATION_DOMAIN>::EnumerateAndFixMatrixSize(const COMPUTATION_DOMAIN<dim>& gref) {
    DOF_indexes_.resize(this->rh_.size());
    fill(DOF_indexes_.begin(), DOF_indexes_.end(), 0);

    if (!this->setup_established_)
      throw csmp::Exception(ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::AssignEssentialConditions",
        "please call EstablishMatrixSetup() prior to this method.");

    if (this->basic_operands_.empty())
      throw csmp::Exception(ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::AssignEssentialConditions",
        "No (basic) operands have been specified...");

    //  Enumerate must be applied AFTER established process.
    //  and BEFORE the accumulated proecess.
    //  When you assemble vector or tensor variables, you also have the option
    //  of only assembling one of their components. To do this just set the
    //  components that you do not want to assemble to DBL_MAX.
    // -------------------------------------------------------------
    size_t DOF(0);

    for ( typename PDE_Integrator<dim,COMPUTATION_DOMAIN>::operandsConstIterator
          it = this->test_operands_.begin(); it != this->test_operands_.end(); it++)
        {
          typename vector<csmp::Node<dim>*>::const_iterator  niter(gref.NodesBegin());
          csmp::Index prop_key = (*it).first.key;
          size_t      offset = (*it).second;

          if (prop_key.place != NODE)
            throw csmp::Exception(ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::AssignEssentialConditions",
              "So far no conditions are assigned to elements, faces, segments");

          size_t  position(0);
          switch (prop_key.type)
          {
          case SCALAR:
            while (niter != gref.NodesEnd()) {
              position = (*niter)->Idx() + offset;
              if ((*niter)->Status(prop_key) == DIRICH)
              {
                DOF_indexes_[position] = NULL_IDX;
              }
              else {
                DOF_indexes_[position] = DOF;
                DOF = DOF + 1;
              }
              niter++;
            }
            break;
          case VECTOR: {
            VectorVariable<dim>  vc;
            while (niter != gref.NodesEnd()) {
              (*niter)->Read(prop_key, vc);
              for (size_t i = 0U; i < dim; i++) {
                position = (*niter)->Idx() * dim + i + offset;
                if (vc.Flag(i) == DIRICH) {
                  DOF_indexes_[position] = NULL_IDX;
                }
                else {
                  DOF_indexes_[position] = DOF;
                  DOF = DOF + 1;
                }
              }
              niter++;
            }
          }
                       break;
          case TENSOR: {
            TensorVariable<dim>  ts;
            while (niter != gref.NodesEnd()) {
              (*niter)->Read(prop_key, ts);
              for (size_t i = 0U; i < dim; i++) {
                if (ts.Flag(i) == DIRICH) {
                  for (size_t j = 0U; j < dim; j++) {
                    position = (*niter)->Idx() * this->dim2_ + i * dim + j + offset;
                    DOF_indexes_[position] = NULL_IDX;
                  }
                }
                else {
                  for (size_t j = 0U; j < dim; j++) {
                    position = (*niter)->Idx() * this->dim2_ + i * dim + j + offset;
                    DOF_indexes_[position] = DOF;
                    DOF = DOF + 1;
                  }
                }
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
                for (size_t i = 0U; i < prop_key.dataDepth; i++) {
                  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                  DOF_indexes_[position] = NULL_IDX;
                }
              }
              else {
                for (size_t i = 0U; i < prop_key.dataDepth; i++) {
                  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                  DOF_indexes_[position] = DOF;
                  DOF = DOF + 1;
                }
              }// end if
              niter++;
            } // end while
          }
                      break;
          case FLAGGEDARRAY: {
            FlaggedArrayVariable  ar(prop_key.dataDepth);
            while (niter != gref.NodesEnd()) {
              (*niter)->Read(prop_key, ar);
              for (size_t i = 0U; i < prop_key.dataDepth; i++)
                if (ar.Flag(i) == DIRICH) {
                  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                  DOF_indexes_[position] = NULL_IDX;
                }
                else {
                  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                  DOF_indexes_[position] = DOF;
                  DOF = DOF + 1;
                }
                niter++;
            }
          }
         break;
          default:
            throw csmp::Exception(FATAL_ERROR,
              "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::AssignEssentialConditions",
              "Variable type not recognised by this method");
          }
        } // end for

    this->G_.Resize(DOF);
    this->rh_.resize(DOF);
    vector<double64>(this->rh_).swap(this->rh_);
    fill(this->rh_.begin(), this->rh_.end(), 0.);
    this->x_.resize(DOF);
    
    // LUAT's new code
    pivotVector_.resize(DOF);
	  fill(pivotVector_.begin(), pivotVector_.end(), 0.);
    vector<double64>(this->x_).swap(this->x_);

 } // end
 


template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator_UoM<dim, COMPUTATION_DOMAIN>::LateAccumulate(const COMPUTATION_DOMAIN<dim>& gref)
  {
    // accumulating as late addition into the righhand vector 'rhs'
    // ------------------------------------------------------------
    for (typename map<string, MathOperatorRHS<dim>*>::const_iterator
      it_rhs = this->rhs_operators_.begin(); it_rhs != this->rhs_operators_.end(); it_rhs++ )
      if ((*it_rhs).second->AddLater() || (*it_rhs).second->SubtractLater())
        for (typename vector<typename COMPUTATION_DOMAIN<dim>::CellType*>::const_iterator
             git = gref.ElementsBegin(); git != gref.ElementsEnd(); git++)
          {
            (*it_rhs).second->GetOperands(*(*git));
            (*it_rhs).second->ComputeContribution(*(*git));
            if ((*it_rhs).second->MultiplyWithTimeIncrement())
              (*it_rhs).second->MultiplyWithTimeFactor( this->time_increment_ );
            if ((*it_rhs).second->DivideByTimeIncrement())
              (*it_rhs).second->MultiplyWithTimeFactor(1. / this->time_increment_ );
            // LUAT's new code
            (*it_rhs).second->AssignToGlobal(*(*git), this->rh_, DOF_indexes_ ); // Luat edited here
          }

  } // end Late Accumulate
  



/**
    Overriding method of the base class because a mapping is required from the smaller solution vector
    back onto the computational domain.
*/
  template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
  void PDE_Integrator_UoM<dim, COMPUTATION_DOMAIN>::OutputResults(COMPUTATION_DOMAIN<dim>& gref)
  {
    Index   prop_key;
    size_t  offset;
    size_t  position;
    for ( typename PDE_Integrator<dim,COMPUTATION_DOMAIN>::operandsIterator
          it = this->basic_operands_.begin(); it != this->basic_operands_.end(); it++ )
      {
        typename vector<Node<dim>*>::iterator  gfirst(gref.NodesBegin());
        prop_key = (*it).first.key;
        offset = (*it).second;

        if (prop_key.place != NODE)
          throw csmp::Exception(ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::OutputResults(Model)",
            "only nodal properties can be output by this method.");
        switch (prop_key.type)
        {
        case SCALAR:
          while (gfirst != gref.NodesEnd()) {
            position = (*gfirst)->Idx() + offset;
            position = DOF_indexes_[position];
            if (position != NULL_IDX) {
              const double64 sc = this->x_[position];
              (*gfirst)->Store(prop_key, makeScalar((*gfirst)->Status(prop_key), sc));
            }
            gfirst++;
          }
          break;
        case VECTOR: {
          VectorVariable<dim>  vc;
          while (gfirst != gref.NodesEnd()) {
            (*gfirst)->Read(prop_key, vc);
            for (size_t i = 0U; i < dim; i++) {
              position = (*gfirst)->Idx() * dim + i + offset;
              position = DOF_indexes_[position];
              if (position != NULL_IDX) {
                vc(i) = this->x_[position];
              }
            }
            (*gfirst)->Store(prop_key, vc);
            gfirst++;
          }
        }
      break;
      
      case TENSOR: {
        TensorVariable<dim>  ts;
        while (gfirst != gref.NodesEnd()) {
          (*gfirst)->Read(prop_key, ts);
          for (size_t i = 0U; i < dim; i++)
            for (size_t k = 0U; k < dim; k++) {
              position = (*gfirst)->Idx() * this->dim2_ + i * dim + k + offset;
              position = DOF_indexes_[position];
              if (position != NULL_IDX) {
                ts(i, k) = this->x_[position];
              }
            }
              
          (*gfirst)->Store(prop_key, ts);
          gfirst++;
        }
      }
                   break;
      case ARRAY: {
        ArrayVariable  ar(prop_key.dataDepth);
        while (gfirst != gref.NodesEnd()) {
          (*gfirst)->Read(prop_key, ar);
          for (size_t i = 0U; i < prop_key.dataDepth; i++) {
            position = (*gfirst)->Idx() * prop_key.dataDepth + i + offset;
            position = DOF_indexes_[position];
            if (position != NULL_IDX) {
              ar(i) = this->x_[position];
            }
          }
          (*gfirst)->Store(prop_key, ar);
          gfirst++;
        }
      }
      break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable  ar(prop_key.dataDepth);
        while (gfirst != gref.NodesEnd()) {
          (*gfirst)->Read(prop_key, ar);
          for (size_t i = 0U; i < prop_key.dataDepth; i++) {
            position = (*gfirst)->Idx() * prop_key.dataDepth + i + offset;
            position = DOF_indexes_[position];
            if (position != NULL_IDX) {
              ar(i) = this->x_[position];
            }
          }
          (*gfirst)->Store(prop_key, ar);
          gfirst++;
        }
      }
      break;
      default:
        throw csmp::Exception(ERROR, "PDE_Integrator_UoM<dim,COMPUTATION_DOMAIN>::OutputResults(Model)",
          "Output to ARRAY type variables is not supported by this method yet.");

      } // end switch(type)

    } // end for

  } // end OutputResults




/**
    New aspect as compared with the base class: reduction of matrix and rhs size.
*/
template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
void PDE_Integrator_UoM<dim, COMPUTATION_DOMAIN>::IntegrateOver( COMPUTATION_DOMAIN<dim>& domain )
  {
    //0. renumber nodes
    domain.RenumberNodes();
    
    // 1. configure algorithm
    EstablishMatrixSetup(domain);
    // => this is importance, since in dynamic changing of DIRICHLET BCs i.e: coupling and decoupling process
	  EnumerateAndFixMatrixSize(domain);
    
    // 2. Accumulation: Note that the conditions that pertain to the group must be input !                                 
    Accumulate(domain);

    // 3. If the computation is transient initial conditions must be input into the righthand vector
    if ( this->Transient() == true ) AssignInitialConditions(domain);

    // 4. If the computation is transient initial conditions must be input into the righthand vector
    if ( this->Transient() == true ) LateAccumulate(domain);

	  // 5. assign conditions like Dirichlet conditions etc.
    // LUAT's extra method
	  AssignEssentialConditions(domain);

    // 6. invert global matrix
    this->Solve();

    // 7. write results back into Model
    OutputResults(domain);

    // 8. Calculation of result-dependent properties                                 
    this->PostProcess(domain);

} // end IntegrateOver



 
/**
    ??? - document
*/
template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
void PDE_Integrator_UoM<dim, COMPUTATION_DOMAIN>::AssignEssentialConditions(const COMPUTATION_DOMAIN<dim>& domain) {
	  for (size_t i(0); i < this->rh_.size(); ++i) {
		     this->rh_[i] += pivotVector_[i];
	    }
 }


template class PDE_Integrator_UoM<1U, Region>;
template class PDE_Integrator_UoM<2U, Region>;
template class PDE_Integrator_UoM<3U, Region>;

template class PDE_Integrator_UoM<1U, Boundary>;
template class PDE_Integrator_UoM<2U, Boundary>;
template class PDE_Integrator_UoM<3U, Boundary>;

template class PDE_Integrator_UoM<1U, SplitBoundary>;
template class PDE_Integrator_UoM<2U, SplitBoundary>;
template class PDE_Integrator_UoM<3U, SplitBoundary>;

template class PDE_Integrator_UoM<1U,NimbleRegion>;
template class PDE_Integrator_UoM<2U,NimbleRegion>;
template class PDE_Integrator_UoM<3U,NimbleRegion>;

} // end namespace csmp
