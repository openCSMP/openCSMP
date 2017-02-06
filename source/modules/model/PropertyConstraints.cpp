#include "PropertyConstraints.h"
#include "Node.h"
#include "Exception.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {

PropertyConstraints::PropertyConstraints()
 : vector_length_check(false),
   one_node_only(false),
   nodal_average(false)
 {
 }
 

PropertyConstraints::PropertyConstraints( const PropertyConstraints& cr )
 {
   *this = cr;
 }
 

PropertyConstraints::PropertyConstraints( const char* prop_name, double64 pmin, double64 pmax )
 : vector_length_check(false),
   one_node_only(false),
   nodal_average(false)
 {
    criteria[ prop_name ] = pair<double64,double64>(pmin,pmax);
 }
 

PropertyConstraints& PropertyConstraints::operator=( const PropertyConstraints& cr )
 {
    if ( &cr != this ) {
         criteria                = cr.criteria;
         check_list              = cr.check_list;
         vector_length_check     = cr.vector_length_check;
         one_node_only           = cr.one_node_only;
         nodal_average           = cr.nodal_average;
      }
    return *this;
 }
 

PropertyConstraints::~PropertyConstraints()
 {
 }
 

bool PropertyConstraints::WithIndexes() const
 {
    return !check_list.empty();
 }
  
    
void PropertyConstraints::SatisfyConstraintsForAtLeastOneNode( bool satisfy )
 {
    one_node_only = satisfy;
 }
 
 
void PropertyConstraints::SatisfyConstraintsForNodalAverage( bool satisfy )
 {
    nodal_average = satisfy;
 }


  

bool PropertyConstraints::AddConstraint( const char* prop_name, double64 pmin, double64 pmax )
 {
    pair<map<string,pair<double64,double64>,less<string> >::iterator,bool>  it;
    string  property(prop_name);
    
    it=criteria.insert( make_pair(property,make_pair(pmin,pmax)) );
    
    return it.second;
 }


void PropertyConstraints::ChangeConstraint( const char* prop_name, double64 pmin, double64 pmax )
 {
    map<string,pair<double64,double64> >::iterator     it;
    map<Index,pair<double64,double64> >::iterator  cit;

    if ( (it=criteria.find(prop_name)) != criteria.end() ) {
         for ( it=criteria.begin(), cit=check_list.begin(); it!=criteria.end(); it++, cit++ )
           if ( (*it).first == prop_name ) {
                (*it).second.first  = (*cit).second.first   = pmin;
                (*it).second.second = (*cit).second.second  = pmax;
             }
      }
    else
    throw csmp::Exception( INFO, "PropertyConstraints::ChangeConstraint", "Constraint did not exist, but was added");
    criteria[ prop_name ] = make_pair(pmin,pmax);
 }
 
    

/**
 
For all constraints the checks are performed until a check fails.
Since the constraints scalar ranges, vector and tensor variables are checked 
for their individual components unless specific check instructions were
specified. If so, for instance, vectors can be checked for  
their length and tensors for their Eigenvalues.  

@return whether the element satisfies the user-supplied constraints.
 */
template<size_t dim>
bool PropertyConstraints::CheckConstraints( const Element<dim>& e, 
                                            Index& failed_upon ) const
 {  
    if ( vector_length_check ) {
         throw csmp::Exception( FATAL_ERROR, "PropertyConstraints::CheckConstraints",
                                      "Tensor Eigenvalue check not implemented yet");
         return false;
      }      

    
    // in this process, the actual variable ranges could be collected into the criteria map
    for ( typename map<csmp::Index,pair<double64,double64> >::const_iterator 
          it=check_list.begin(); it!=check_list.end(); it++ )
      {
         if ( vector_length_check && (*it).first.place )
           return VectorLengthCheck( e, (*it).first, (*it).second.first, (*it).second.second );

         switch( (*it).first.place ) {
              case NODE:
                   if ( (*it).first.type == SCALAR ) {
                        ScalarVariable sc;
                        for ( size_t i=0; i<e.Nodes(); i++ ) {
                            e.N(i)->Read( (*it).first, sc );
                            if ( (*it).second.first  > sc() ||
                                 (*it).second.second < sc() ) {
                                 failed_upon = (*it).first;
                                 return false;
                              }
                          }
                     }
                   else if ( (*it).first.type == VECTOR ) {
                        VectorVariable<dim> vc;
                        for ( size_t i=0; i<e.Nodes(); i++ ) {
                            e.N(i)->Read( (*it).first, vc );
                            if ( vector_length_check ) {
                              if ( (*it).second.first  > vc.Length() ||
                                   (*it).second.second < vc.Length() ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                              }
                            else
                            for ( size_t j=0; j<dim; j++ )
                              if ( (*it).second.first  > vc(j) ||
                                   (*it).second.second < vc(j) ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                          }
                     }
                   else if ( (*it).first.type == TENSOR ) {
                        TensorVariable<dim> ts;
                        for ( size_t i=0; i<e.Nodes(); i++ ) {
                            e.N(i)->Read( (*it).first, ts );
                            for ( size_t j=0; j<dim; j++ )
                              for ( size_t k=0; k<dim; k++ )
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
                        for ( size_t i=0; i<e.IntegrationPoints(); i++ ) {
                            e.Read( i, (*it).first, sc );
                            if ( (*it).second.first  > sc() ||
                                 (*it).second.second < sc() ) {
                                 failed_upon = (*it).first;
                                 return false;
                              }
                          }
                     }
                   else if ( (*it).first.type == VECTOR ) {
                        VectorVariable<dim> vc;
                        for ( size_t i=0; i<e.IntegrationPoints(); i++ ) {
                            e.Read( i, (*it).first, vc );
                            if ( vector_length_check ) {
                              if ( (*it).second.first  > vc.Length() ||
                                   (*it).second.second < vc.Length() ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                              }
                            else
                            for ( size_t j=0; j<dim; j++ )
                              if ( (*it).second.first  > vc(j) ||
                                   (*it).second.second < vc(j) ) {
                                    failed_upon = (*it).first;
                                    return false;
                                 }
                          }
                     }
                   else if ( (*it).first.type == TENSOR ) {
                        TensorVariable<dim> ts;
                        for ( size_t i=0; i<e.IntegrationPoints(); i++ ) {
                            e.Read( i, (*it).first, ts );
                            for ( size_t j=0; j<dim; j++ )
                              for ( size_t k=0; k<dim; k++ )
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
                        e.Read( (*it).first, sc );
                        if ( (*it).second.first  > sc() ||
                             (*it).second.second < sc() ) {
                                 failed_upon = (*it).first;
                                 return false;
                          }
                     }
                   else if ( (*it).first.type == VECTOR ) {
                        VectorVariable<dim>  vc;
                        e.Read( (*it).first, vc );
                        if ( vector_length_check ) {
                           if ( (*it).second.first  > vc.Length() ||
                               (*it).second.second < vc.Length() ) {
                                failed_upon = (*it).first;
                                return false;
                             }
                          }
                        else
                        for ( size_t j=0; j<dim; j++ )
                          if ( (*it).second.first  > vc(j) ||
                               (*it).second.second < vc(j) ) {
                                failed_upon = (*it).first;
                                return false;
                            }
                     }
                   else if ( (*it).first.type == TENSOR ) {
                        TensorVariable<dim> ts;
                        e.Read( (*it).first, ts );
                        for ( size_t j=0; j<dim; j++ )
                          for ( size_t k=0; k<dim; k++ )
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
template<size_t dim>
bool PropertyConstraints::CheckConstraints( const Element<dim>& e ) const
 {  
    if ( one_node_only )
      return CheckSingleNodeConstraints( e );

    if ( nodal_average )
      return CheckNodeAverageConstraints( e );

    // in this process, the actual variable ranges could be collected into the criteria map
    for ( typename map<Index,pair<double64,double64> >::const_iterator
          it=check_list.begin(); it!=check_list.end(); it++ )
      {
         if ( vector_length_check && (*it).first.place )
           return VectorLengthCheck( e, (*it).first, (*it).second.first, (*it).second.second );

         switch( (*it).first.place ) {
              case NODE:
                   for ( size_t i=0; i<e.Nodes(); i++ )
                     if ( !e.N(i)->IsWithinRange( (*it).first, (*it).second.first, (*it).second.second ) )
                       return false;
                break;
              case ELEMENT_INTEGRATION_POINT:
                   for ( size_t i=0; i<e.IntegrationPoints(); i++ )
                     if ( !e.IsWithinRange( i, (*it).first, (*it).second.first, (*it).second.second ) )
                       return false;
                break;
              case ELEMENT:
                   if ( !e.IsWithinRange( (*it).first, (*it).second.first, (*it).second.second ) )
                     return false;
                break;
              default:
                   throw csmp::Exception( ERROR, "PropertyConstraints::CheckConstraints",
                                              "Placement of constraint variable could not be identified");
           }
      }
      
    return true;
    
 } // end CheckConstraints







size_t  PropertyConstraints::Constraints() const { return criteria.size(); }
 
void PropertyConstraints::CheckLengthOfVectorVariables( bool check ) { vector_length_check=check; }
 
    
void PropertyConstraints::Erase()
 {
    criteria.erase( criteria.begin(), criteria.end() );
    check_list.erase( check_list.begin(), check_list.end() );
 }
 


template<size_t dim>
bool PropertyConstraints::CheckSingleNodeConstraints( const Element<dim>& e ) const
 {
    size_t i, counter;
    
    for ( typename map<Index,pair<double64,double64> >::const_iterator
          it=check_list.begin(); it!=check_list.end(); it++ )
      {
         if ( vector_length_check && (*it).first.place )
           return VectorLengthCheck( e, (*it).first, (*it).second.first, (*it).second.second );

         switch( (*it).first.place ) {
              case NODE:
                   for ( counter=i=0U; i<e.Nodes(); i++ ) {
                     if ( !e.N(i)->IsWithinRange( (*it).first, (*it).second.first, (*it).second.second ) ) 
                          counter++;
                     }   
                   if ( counter >= e.Nodes()-1U ) return false; 
                break;
              case ELEMENT_INTEGRATION_POINT:
                   for ( size_t i=0; i<e.IntegrationPoints(); i++ )
                     if ( !e.IsWithinRange( i, (*it).first, (*it).second.first, (*it).second.second ) )
                       return false;
                break;
              case ELEMENT:
                   if ( !e.IsWithinRange( (*it).first, (*it).second.first, (*it).second.second ) )
                     return false;
                break;
              default:
                   throw csmp::Exception( ERROR, "PropertyConstraints::CheckSingleNodeConstraints",
                                              "Placement of constraint variable could not be identified");
           }
      }
    return true;
 }
 
 
template<size_t dim>
bool PropertyConstraints::CheckNodeAverageConstraints( const Element<dim>& e ) const
 {
    for ( typename map<Index,pair<double64,double64> >::const_iterator
          it=check_list.begin(); it!=check_list.end(); it++ )
      {
         if ( vector_length_check && (*it).first.place )
           return VectorLengthCheck( e, (*it).first, (*it).second.first, (*it).second.second );

         csmp::Index index = (*it).first;
         
         switch( (*it).first.place ) {
              case NODE:
                   if ( (*it).first.type == SCALAR ) {
                        ScalarVariable sc;
                        e.PropertyValueAtBaryCenter( index, sc );
                        if ( !sc.IsWithinRange( (*it).second.first, (*it).second.second ) ) return false; 
                     }
                   else if ( (*it).first.type == VECTOR ) {
                        VectorVariable<dim> vc;
                        e.PropertyValueAtBaryCenter( index, vc );
                        if ( !vc.IsWithinRange( (*it).second.first, (*it).second.second ) ) return false; 
                     }
                   else if ( (*it).first.type == TENSOR ) {
                        TensorVariable<dim> ts;
                        e.PropertyValueAtBaryCenter( index, ts );
                        if ( !ts.IsWithinRange( (*it).second.first, (*it).second.second ) ) return false; 
                     }
                break;
              case ELEMENT_INTEGRATION_POINT:
                   for ( size_t i=0; i<e.IntegrationPoints(); i++ )
                     if ( !e.IsWithinRange( i, (*it).first, (*it).second.first, (*it).second.second ) )
                       return false;
                break;
              case ELEMENT:
                   if ( !e.IsWithinRange( (*it).first, (*it).second.first, (*it).second.second ) )
                     return false;
                break;
              default:
                   throw csmp::Exception( ERROR, "PropertyConstraints::CheckNodeAverageConstraints",
                                              "Placement of constraint variable could not be identified");
           }
      }
    return true;
 }




template<size_t dim>
bool PropertyConstraints::VectorLengthCheck( const Element<dim>& e, 
                                             const csmp::Index& idx,
                                             double64 vmin, double64 vmax ) const
 {
    VectorVariable<dim>  vc;
 
    if ( nodal_average ) {
         throw csmp::Exception( FATAL_ERROR, "PropertyConstraints::VectorLengthCheck",
                                      "'nodal average' and 'vector_length_check' are mutually exclusive switches");
         return false;
      }
 
    if ( idx.place == NODE ) {
         if ( one_node_only ) {
              for ( size_t i=0; i<e.Nodes(); i++ ) {
                   e.N(i)->Read( idx, vc );
                   if ( vc.IsWithinRange( vmin, vmax ) ) return true; 
                }
              return false;
           }
         else
         for ( size_t i=0; i<e.Nodes(); i++ ) {
              e.N(i)->Read( idx, vc );
              if ( !vc.IsWithinRange( vmin, vmax ) ) return false; 
           }
      }
    else if ( idx.place == ELEMENT_INTEGRATION_POINT ) {
         for ( size_t i=0; i<e.IntegrationPoints(); i++ ) {
              e.Read( i, idx, vc );
              if ( !vc.IsWithinRange( vmin, vmax ) ) return false; 
           }
      }
    else if ( idx.place == ELEMENT ) {
              e.Read( idx, vc );
              if ( !vc.IsWithinRange( vmin, vmax ) ) return false; 
      }
      
    return true;
 }



void PropertyConstraints::Out(std::ostream& os) const
 {
    os <<"\nPropertyConstraints::Out: "<< endl;
    if ( vector_length_check )      os <<"\tSet to check the length of vector variables"<< endl;
    if ( one_node_only )            os <<"\tIf at least one node matches criteria, constraints are satisfied"<< endl;
    if ( nodal_average )            os <<"\tAll nodal variables are averaged"<< endl;
    if ( check_list.empty() )       os <<"\tIndex data are not established yet."<< endl;
    else                            os <<"\tIndex data have been established."<< endl;

    map<string,pair<double64,double64> >::const_iterator  crit;       
    
    os <<"\tAssigned property constraints, and their ranges:";
    for ( crit=criteria.begin(); crit!=criteria.end(); crit++ )
      os <<"\n\t\t'"<< (*crit).first <<"' range: "<< (*crit).second.first <<" to "<< (*crit).second.second;
    
    os << endl; 
 }

 // 1d
template bool PropertyConstraints::CheckSingleNodeConstraints<1U>(
                                              const Element<1U>& e ) const;

template bool PropertyConstraints::CheckNodeAverageConstraints<1U>(
                                              const Element<1U>& e ) const;

template bool PropertyConstraints::VectorLengthCheck<1U>( 
                                              const Element<1U>& e,
                                              const csmp::Index& idx, double64 vmin, double64 vmax ) const;

template bool PropertyConstraints::CheckConstraints<1U>( 
                                             const Element<1U>& e ) const;

template bool PropertyConstraints::CheckConstraints<1U>( 
                                             const Element<1U>& e, Index& idx ) const;

// 2d
template bool PropertyConstraints::CheckSingleNodeConstraints<2U>(
                                              const Element<2U>& e ) const;

template bool PropertyConstraints::CheckNodeAverageConstraints<2U>(
                                              const Element<2U>& e ) const;

template bool PropertyConstraints::VectorLengthCheck<2U>( 
                                              const Element<2U>& e,
                                              const csmp::Index& idx, double64 vmin, double64 vmax ) const;

template bool PropertyConstraints::CheckConstraints<2U>( 
                                             const Element<2U>& e ) const;

template bool PropertyConstraints::CheckConstraints<2U>( 
                                             const Element<2U>& e, Index& idx ) const;

// 3d
template bool PropertyConstraints::CheckSingleNodeConstraints<3U>(
                                              const Element<3U>& e ) const;
    
template bool PropertyConstraints::CheckNodeAverageConstraints<3U>( 
                                              const Element<3U>& e ) const;
    
template bool PropertyConstraints::VectorLengthCheck<3U>( 
                                              const Element<3U>& e, 
                                              const csmp::Index& idx, double64 vmin, double64 vmax ) const;

template bool PropertyConstraints::CheckConstraints<3U>( 
                                             const Element<3U>& e ) const;

template bool PropertyConstraints::CheckConstraints<3U>( 
                                             const Element<3U>& e, Index& idx ) const;


} // end namespace csp











