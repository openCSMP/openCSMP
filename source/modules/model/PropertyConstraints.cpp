#include "PropertyConstraints.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Node.h"
#include "Exception.h"
#include "PropertyDatabase.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

PropertyConstraints::PropertyConstraints()
 {
 }
 

PropertyConstraints::PropertyConstraints( const PropertyConstraints& cr )
 :  criteria_(cr.criteria_),
    check_list_(cr.check_list_),
    vector_length_check_(cr.vector_length_check_),
    one_node_only_(cr.one_node_only_),
    nodal_average_(cr.nodal_average_)
 {
 }
 
template<uint32_t dim>
PropertyConstraints::PropertyConstraints( const PropertyDatabase<dim>& database, const char* prop_name, double pmin, double pmax )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
  
    if ( database.IsDefined(prop_name) == false ) {
         csmp_error.Note( ERROR, "PropertyConstraints (ctor)", prop_name, "is undefined; constraint was ignored" );
         return;
      }
 
    // map<string,pair<double,double> >
    criteria_[ prop_name ] = pair<double,double>(pmin,pmax);
    
    // map<Index,pair<double,double> >
    check_list_.insert( make_pair( database.StorageKey(prop_name), make_pair(pmin,pmax) ) );
 }

template PropertyConstraints::PropertyConstraints( const PropertyDatabase<3>&, const char*, double, double );
template PropertyConstraints::PropertyConstraints( const PropertyDatabase<2>&, const char*, double, double );
template PropertyConstraints::PropertyConstraints( const PropertyDatabase<1>&, const char*, double, double );




PropertyConstraints& PropertyConstraints::operator=( const PropertyConstraints& cr )
 {
    if ( &cr != this ) {
         criteria_                = cr.criteria_;
         check_list_              = cr.check_list_;
         vector_length_check_     = cr.vector_length_check_;
         one_node_only_           = cr.one_node_only_;
         nodal_average_           = cr.nodal_average_;
      }
    return *this;
 }
 
 

bool PropertyConstraints::WithIndexes() const
 {
    return !check_list_.empty();
 }
  
    
void PropertyConstraints::SatisfyConstraintsForAtLeastOneNode( bool satisfy )
 {
    one_node_only_ = satisfy;
 }
 
 
void PropertyConstraints::SatisfyConstraintsForNodalAverage( bool satisfy )
 {
    nodal_average_ = satisfy;
 }


  

bool PropertyConstraints::AddConstraint( const char* prop_name, double pmin, double pmax )
 {
    string  property(prop_name);
    
    pair<map<string,pair<double,double>,less<string> >::iterator,bool>
      it = criteria_.insert( make_pair(property,make_pair(pmin,pmax)) );
    
    return it.second;
 }


void PropertyConstraints::ChangeConstraint( const char* prop_name, double pmin, double pmax )
 {
    map<string,pair<double,double> >::iterator it;
    map<Index,pair<double,double> >::iterator  cit;

    if ( (it=criteria_.find(prop_name)) != criteria_.end() ) {
         for ( it=criteria_.begin(), cit=check_list_.begin(); it!=criteria_.end(); it++, cit++ )
           if ( (*it).first == prop_name ) {
                (*it).second.first  = (*cit).second.first   = pmin;
                (*it).second.second = (*cit).second.second  = pmax;
             }
      }
    else
    throw csmp::Exception( ERROR, "PropertyConstraints::ChangeConstraint", "Constraint does not exist");
 }
 


template<uint32_t dim>
bool PropertyConstraints::InitializePropertyIndices( const PropertyDatabase<dim>& pref )
  {
    if ( criteria_.empty() ) {
        throw csmp::Exception( ERROR, "PropertyConstraints::InitializePropertyIndices",
          "No criteria have been defined so far");
        return false;
      }

    if ( !check_list_.empty() && check_list_.size() == criteria_.size() ) return true; 

    for ( const auto& it : criteria_ )
      check_list_[ pref.StorageKey( it.first.c_str() ) ] = make_pair( it.second.first, it.second.second );

    return true;
  }

template bool PropertyConstraints::InitializePropertyIndices( const PropertyDatabase<1U>& );
template bool PropertyConstraints::InitializePropertyIndices( const PropertyDatabase<2U>& );
template bool PropertyConstraints::InitializePropertyIndices( const PropertyDatabase<3U>& );




template<uint32_t dim>
void PropertyConstraints::DeleteConstraint( const PropertyDatabase<dim>& pref, const char* prop_name )
  {
     criteria_.erase( prop_name );
     check_list_.erase( pref.StorageKey(prop_name) );
  }

template void PropertyConstraints::DeleteConstraint( const PropertyDatabase<1U>&, const char* );
template void PropertyConstraints::DeleteConstraint( const PropertyDatabase<2U>&, const char* );
template void PropertyConstraints::DeleteConstraint( const PropertyDatabase<3U>&, const char* );
  

/**
 
For all constraints the checks are performed until a check fails.
Since the constraints scalar ranges, vector and tensor variables are checked 
for their individual components unless specific check instructions were
specified. If so, for instance, vectors can be checked for  
their length and tensors for their Eigenvalues.  

@return whether the element satisfies the user-supplied constraints.
 */
template<uint32_t dim, template<uint32_t> class CELL>
bool PropertyConstraints::CheckConstraints( const CELL<dim>* e, Index& failed_upon ) const
 {  
    if ( vector_length_check_ ) {
         throw csmp::Exception( FATAL_ERROR, "PropertyConstraints::CheckConstraints",
                                      "Tensor Eigenvalue check not implemented yet");
         return false;
      }      

    
    // in this process, the actual variable ranges could be collected into the criteria map
    for ( typename map<csmp::Index,pair<double,double> >::const_iterator 
          it=check_list_.begin(); it!=check_list_.end(); it++ )
      {
         if ( vector_length_check_ && (*it).first.place )
           return VectorLengthCheck( e, (*it).first, (*it).second.first, (*it).second.second );

         switch( (*it).first.place ) {
              case NODE:
                   if ( (*it).first.type == SCALAR ) {
                        ScalarVariable sc;
                        for ( uint32_t i{0U}; i<e->Nodes(); i++ ) {
                            e->N(i)->Read( (*it).first, sc );
                            if ( (*it).second.first  > sc() ||
                                 (*it).second.second < sc() ) {
                                 failed_upon = (*it).first;
                                 return false;
                              }
                          }
                     }
                   else if ( (*it).first.type == VECTOR ) {
                        VectorVariable<dim> vc;
                        for ( uint32_t i{0U}; i<e->Nodes(); i++ ) {
                            e->N(i)->Read( (*it).first, vc );
                            if ( vector_length_check_ ) {
                              if ( (*it).second.first  > vc.Length() ||
                                   (*it).second.second < vc.Length() ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                              }
                            else
                            for ( uint32_t j{0U}; j<dim; j++ )
                              if ( (*it).second.first  > vc(j) ||
                                   (*it).second.second < vc(j) ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                          }
                     }
                   else if ( (*it).first.type == TENSOR ) {
                        TensorVariable<dim> ts;
                        for ( auto i{0U}; i<e->Nodes(); i++ ) {
                            e->N(i)->Read( (*it).first, ts );
                            for ( uint32_t j{0U}; j<dim; j++ )
                              for ( uint32_t k=0; k<dim; k++ )
                                if ( (*it).second.first  > ts(j,k) ||
                                     (*it).second.second < ts(j,k) ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                          }
                     }
                break;
              case ELEMENT_INTEGRATION_POINT:
                   if ( (*it).first.type == SCALAR ) {
                        ScalarVariable sc;
                        for ( uint32_t i{0U}; i<e->IntegrationPoints(); i++ ) {
                            e->Read( i, (*it).first, sc );
                            if ( (*it).second.first  > sc() ||
                                 (*it).second.second < sc() ) {
                                 failed_upon = (*it).first;
                                 return false;
                              }
                          }
                     }
                   else if ( (*it).first.type == VECTOR ) {
                        VectorVariable<dim> vc;
                        for ( uint32_t i{0U}; i<e->IntegrationPoints(); i++ ) {
                            e->Read( i, (*it).first, vc );
                            if ( vector_length_check_ ) {
                              if ( (*it).second.first  > vc.Length() ||
                                   (*it).second.second < vc.Length() ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                              }
                            else
                            for ( uint32_t j{0U}; j<dim; j++ )
                              if ( (*it).second.first  > vc(j) ||
                                   (*it).second.second < vc(j) ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                          }
                     }
                   else if ( (*it).first.type == TENSOR ) {
                        TensorVariable<dim> ts;
                        for ( auto i{0U}; i<e->IntegrationPoints(); i++ ) {
                            e->Read( i, (*it).first, ts );
                            for ( uint32_t j{0U}; j<dim; j++ )
                              for ( uint32_t k=0; k<dim; k++ )
                                if ( (*it).second.first  > ts(j,k) ||
                                     (*it).second.second < ts(j,k) ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                          }
                     }
                break;
              case ELEMENT:
                   if ( (*it).first.type == SCALAR ) {
                        ScalarVariable  sc;
                        e->Read( (*it).first, sc );
                        if ( (*it).second.first  > sc() ||
                             (*it).second.second < sc() ) {
                                 failed_upon = (*it).first;
                                 return false;
                          }
                     }
                   else if ( (*it).first.type == VECTOR ) {
                        VectorVariable<dim>  vc;
                        e->Read( (*it).first, vc );
                        if ( vector_length_check_ ) {
                           if ( (*it).second.first  > vc.Length() ||
                               (*it).second.second < vc.Length() ) {
                                failed_upon = (*it).first;
                                return false;
                             }
                          }
                        else
                        for ( uint32_t j{0U}; j<dim; j++ )
                          if ( (*it).second.first  > vc(j) ||
                               (*it).second.second < vc(j) ) {
                                failed_upon = (*it).first;
                                return false;
                            }
                     }
                   else if ( (*it).first.type == TENSOR ) {
                        TensorVariable<dim> ts;
                        e->Read( (*it).first, ts );
                        for ( uint32_t j{0U}; j<dim; j++ )
                          for ( uint32_t k=0; k<dim; k++ )
                            if ( (*it).second.first  > ts(j,k) ||
                                 (*it).second.second < ts(j,k) ) {
                                 failed_upon = (*it).first;
                                 return false;
                              }
                     }
                break;
              default:
                   throw csmp::Exception( ERROR, "PropertyConstraints::CheckConstraints",
                                              "Placement of constraint variable could not be identified");
           }
      }
    return true;
    
 } // end CheckConstraints
  



/**
*/
template<uint32_t dim, template<uint32_t> class CELL>
bool PropertyConstraints::CheckConstraints( const CELL<dim>* e ) const
 {  
    if ( one_node_only_ )
      return CheckSingleNodeConstraints( e );

    if ( nodal_average_ )
      return CheckNodeAverageConstraints( e );

    // in this process, the actual variable ranges could be collected into the criteria map
    for ( const auto& it : check_list_ )
      {
         if ( vector_length_check_ && it.first.type == VECTOR )
           return VectorLengthCheck( e, it.first, it.second.first, it.second.second );

         switch( it.first.place ) {
              case NODE:
                   for ( uint32_t i{0U}; i<e->Nodes(); i++ )
                     if ( !e->N(i)->IsWithinRange( it.first, it.second.first, it.second.second ) )
                       return false;
                break;
              case ELEMENT_INTEGRATION_POINT:
                   for ( uint32_t i{0U}; i<e->IntegrationPoints(); i++ )
                     if ( !e->IsWithinRange( i, it.first, it.second.first, it.second.second ) )
                       return false;
                break;
              case ELEMENT:
                   if ( !e->IsWithinRange( it.first, it.second.first, it.second.second ) )
                     return false;
                break;
              default:
                   throw csmp::Exception( ERROR, "PropertyConstraints::CheckConstraints",
                                              "Placement of constraint variable could not be identified");
           }
      }
      
    return true;
    
 } // end CheckConstraints







uint32_t  PropertyConstraints::Constraints() const { return static_cast<uint32_t>(criteria_.size()); }
 
 
void PropertyConstraints::CheckLengthOfVectorVariables( bool check ) { vector_length_check_=check; }
 
 
 
    
void PropertyConstraints::Erase()
 {
    criteria_.clear();
    check_list_.clear();
 }
 




template<uint32_t dim, template<uint32_t> class CELL>
bool PropertyConstraints::CheckSingleNodeConstraints( const CELL<dim>* e ) const
{
    for ( const auto& it : check_list_ )
      {
         // vector length check is a special case — handled separately
         if ( vector_length_check_ && it.first.type == VECTOR )
           return VectorLengthCheck( e, it.first, it.second.first, it.second.second );

         switch ( it.first.place ) {
              case NODE:
                 {
                   // count nodes that fall outside the constraint range.
                   // IsWithinRange handles scalar, vector, and tensor types
                   // internally — PropertyConstraints does not need to
                   // distinguish between them here.
                   uint32_t out_of_range{ 0U };
                   for ( uint32_t i{ 0U }; i < e->Nodes(); ++i )
                     if ( !e->N(i)->IsWithinRange( it.first,
                                                   it.second.first,
                                                   it.second.second ) )
                       ++out_of_range;

                   // fail only if NO node satisfies the constraint
                   if ( out_of_range >= e->Nodes() ) return false;
                 }
                break;

              case ELEMENT_INTEGRATION_POINT:
                   // single-node mode does not apply to integration point
                   // variables — use the strict all-points check
                   for ( uint32_t i{ 0U }; i < e->IntegrationPoints(); ++i )
                     if ( !e->IsWithinRange( i, it.first,
                                             it.second.first,
                                             it.second.second ) )
                       return false;
                break;

              case ELEMENT:
                   if ( !e->IsWithinRange( it.first,
                                           it.second.first,
                                           it.second.second ) )
                     return false;
                break;

              default:
                   throw csmp::Exception( ERROR,
                       "PropertyConstraints::CheckSingleNodeConstraints",
                       "Placement of constraint variable could not be identified" );
           }
      }
    return true;
}
 
 



template<uint32_t dim, template<uint32_t> class CELL>
bool PropertyConstraints::CheckNodeAverageConstraints( const CELL<dim>* e ) const
{
    for ( const auto& it : check_list_ )
      {
         // vector length check is a special case — handled separately
         if ( vector_length_check_ && it.first.type == VECTOR )
           return VectorLengthCheck( e, it.first, it.second.first, it.second.second );

         switch ( it.first.place ) {
              case NODE:
                 {
                   // PropertyValueAtBaryCenter interpolates nodal values to
                   // the barycentre using shape functions. IsWithinRange on
                   // the resulting variable handles scalar, vector, and tensor
                   // types internally — no type dispatch needed here.
                   if ( it.first.type == SCALAR ) {
                        ScalarVariable sc;
                        e->PropertyValueAtBaryCenter( it.first, sc );
                        if ( !sc.IsWithinRange( it.second.first,
                                                it.second.second ) ) return false;
                     }
                   else if ( it.first.type == VECTOR ) {
                        VectorVariable<dim> vc;
                        e->PropertyValueAtBaryCenter( it.first, vc );
                        if ( !vc.IsWithinRange( it.second.first,
                                                it.second.second ) ) return false;
                     }
                   else if ( it.first.type == TENSOR ) {
                        TensorVariable<dim> ts;
                        e->PropertyValueAtBaryCenter( it.first, ts );
                        if ( !ts.IsWithinRange( it.second.first,
                                                it.second.second ) ) return false;
                     }
                 }
                break;

              case ELEMENT_INTEGRATION_POINT:
                   for ( uint32_t i{ 0U }; i < e->IntegrationPoints(); ++i )
                     if ( !e->IsWithinRange( i, it.first,
                                             it.second.first,
                                             it.second.second ) )
                       return false;
                break;

              case ELEMENT:
                   if ( !e->IsWithinRange( it.first,
                                           it.second.first,
                                           it.second.second ) )
                     return false;
                break;

              default:
                   throw csmp::Exception( ERROR,
                       "PropertyConstraints::CheckNodeAverageConstraints",
                       "Placement of constraint variable could not be identified" );
           }
      }
    return true;
}




template<uint32_t dim, template<uint32_t> class CELL>
bool PropertyConstraints::VectorLengthCheck( const CELL<dim>* e,
                                             const csmp::Index& idx,
                                             double vmin, double vmax ) const
{
    VectorVariable<dim> vc;

    if ( idx.type != VECTOR )
        throw csmp::Exception( ERROR, "PropertyConstraints::VectorLengthCheck",
                               "the variable index must define a VectorVariable" );

    if ( nodal_average_ )
        throw csmp::Exception( ERROR, "PropertyConstraints::VectorLengthCheck",
                               "'nodal average' and 'vector_length_check' are mutually exclusive switches" );

    if ( idx.place == NODE ) {
        if ( one_node_only_ ) {
            for ( uint32_t i{ 0U }; i < e->Nodes(); ++i ) {
                e->N(i)->Read( idx, vc );
                if ( vc.Length() >= vmin && vc.Length() <= vmax ) return true;
            }
            return false;
        }
        else {
            for ( uint32_t i{ 0U }; i < e->Nodes(); ++i ) {
                e->N(i)->Read( idx, vc );
                if ( vc.Length() < vmin || vc.Length() > vmax ) return false;
            }
        }
    }
    else if ( idx.place == ELEMENT_INTEGRATION_POINT ) {
        for ( uint32_t i{ 0U }; i < e->IntegrationPoints(); ++i ) {
            e->Read( i, idx, vc );
            if ( vc.Length() < vmin || vc.Length() > vmax ) return false;
        }
    }
    else if ( idx.place == ELEMENT ) {
        e->Read( idx, vc );
        if ( vc.Length() < vmin || vc.Length() > vmax ) return false;
    }

    return true;
}



void PropertyConstraints::Out() const
 {
    cout <<"\nPropertyConstraints::Out: "<< endl;
    if ( vector_length_check_ )      cout <<"\tSet to check the length of vector variables"<< endl;
    if ( one_node_only_ )            cout <<"\tIf at least one node matches criteria, constraints are satisfied"<< endl;
    if ( nodal_average_ )            cout <<"\tAll nodal variables are averaged"<< endl;
    if ( check_list_.empty() )       cout <<"\tIndex data are not established yet."<< endl;
    else                            cout <<"\tIndex data have been established."<< endl;

    cout <<"\tAssigned property constraints, and their ranges:";
    for ( auto crit=criteria_.begin(); crit!=criteria_.end(); crit++ )
      cout <<"\n\t\t'"<< (*crit).first <<"' range: "<< (*crit).second.first <<" to "<< (*crit).second.second;
    
    cout << endl; 
 }

 // 1d
template bool PropertyConstraints::CheckSingleNodeConstraints<1U>( const Element<1>* ) const;
template bool PropertyConstraints::CheckNodeAverageConstraints<1U>( const Element<1>* ) const;
template bool PropertyConstraints::VectorLengthCheck<1U>( const Element<1>*, const csmp::Index&, double, double ) const;
template bool PropertyConstraints::CheckConstraints<1U>( const Element<1>* ) const;
template bool PropertyConstraints::CheckConstraints<1U>( const Element<1>*, Index& ) const;

// 2d
template bool PropertyConstraints::CheckSingleNodeConstraints<2U>( const Element<2>* ) const;
template bool PropertyConstraints::CheckNodeAverageConstraints<2U>( const Element<2>* ) const;
template bool PropertyConstraints::VectorLengthCheck<2U>( const Element<2>*, const csmp::Index&, double, double ) const;
template bool PropertyConstraints::CheckConstraints<2U>( const Element<2>* ) const;
template bool PropertyConstraints::CheckConstraints<2U>( const Element<2>*, Index& ) const;

// 3d
template bool PropertyConstraints::CheckSingleNodeConstraints<3U>( const Element<3>* ) const;
template bool PropertyConstraints::CheckNodeAverageConstraints<3U>( const Element<3>* ) const;
template bool PropertyConstraints::VectorLengthCheck<3U>( const Element<3>*, const csmp::Index&, double, double ) const;
template bool PropertyConstraints::CheckConstraints<3U>( const Element<3>* ) const;
template bool PropertyConstraints::CheckConstraints<3U>( const Element<3>*, Index& ) const;


 // 1d
template bool PropertyConstraints::CheckSingleNodeConstraints<1U>( const Face<1>* ) const;
template bool PropertyConstraints::CheckNodeAverageConstraints<1U>( const Face<1>* ) const;
template bool PropertyConstraints::VectorLengthCheck<1U>( const Face<1>*, const csmp::Index&, double, double ) const;
template bool PropertyConstraints::CheckConstraints<1U>( const Face<1>* ) const;
template bool PropertyConstraints::CheckConstraints<1U>( const Face<1>*, Index& ) const;

// 2d
template bool PropertyConstraints::CheckSingleNodeConstraints<2U>( const Face<2>* ) const;
template bool PropertyConstraints::CheckNodeAverageConstraints<2U>( const Face<2>* ) const;
template bool PropertyConstraints::VectorLengthCheck<2U>( const Face<2>*, const csmp::Index&, double, double ) const;
template bool PropertyConstraints::CheckConstraints<2U>( const Face<2>* ) const;
template bool PropertyConstraints::CheckConstraints<2U>( const Face<2>*, Index& ) const;

// 3d
template bool PropertyConstraints::CheckSingleNodeConstraints<3U>( const Face<3>* ) const;
template bool PropertyConstraints::CheckNodeAverageConstraints<3U>( const Face<3>* ) const;
template bool PropertyConstraints::VectorLengthCheck<3U>( const Face<3>*, const csmp::Index&, double, double ) const;
template bool PropertyConstraints::CheckConstraints<3U>( const Face<3>* ) const;
template bool PropertyConstraints::CheckConstraints<3U>( const Face<3>*, Index& ) const;


 // 1d
template bool PropertyConstraints::CheckSingleNodeConstraints<1U>( const InterFace<1>* ) const;
template bool PropertyConstraints::CheckNodeAverageConstraints<1U>( const InterFace<1>* ) const;
template bool PropertyConstraints::VectorLengthCheck<1U>( const InterFace<1>*, const csmp::Index&, double, double ) const;
template bool PropertyConstraints::CheckConstraints<1U>( const InterFace<1>* ) const;
template bool PropertyConstraints::CheckConstraints<1U>( const InterFace<1>*, Index& ) const;

// 2d
template bool PropertyConstraints::CheckSingleNodeConstraints<2U>( const InterFace<2>* ) const;
template bool PropertyConstraints::CheckNodeAverageConstraints<2U>( const InterFace<2>* ) const;
template bool PropertyConstraints::VectorLengthCheck<2U>( const InterFace<2>*, const csmp::Index&, double, double ) const;
template bool PropertyConstraints::CheckConstraints<2U>( const InterFace<2>* ) const;
template bool PropertyConstraints::CheckConstraints<2U>( const InterFace<2>*, Index& ) const;

// 3d
template bool PropertyConstraints::CheckSingleNodeConstraints<3U>( const InterFace<3>* ) const;
template bool PropertyConstraints::CheckNodeAverageConstraints<3U>( const InterFace<3>* ) const;
template bool PropertyConstraints::VectorLengthCheck<3U>( const InterFace<3>*, const csmp::Index&, double, double ) const;
template bool PropertyConstraints::CheckConstraints<3U>( const InterFace<3>* ) const;
template bool PropertyConstraints::CheckConstraints<3U>( const InterFace<3>*, Index& ) const;


} // end namespace csp











