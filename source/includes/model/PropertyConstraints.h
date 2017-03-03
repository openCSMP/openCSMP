#ifndef PROPERTY_CONSTRAINTS_H
#define PROPERTY_CONSTRAINTS_H

#include "CSMP_definitions.h"
#include "Element.h"

namespace csmp {

/**
    @brief PropertyConstraints permits to combine multiple criteria on the 
    basis of which non-unique regions can be formed.
    
    In the terminology of this class, a criterion is called a constraint.
    This class stores user-defined constraints that may include, for instance, 
    the magnitude of a vector.
    
    @note look in the region interface to see how to build regions
    using PropertyConstraints.
    
    @author SKM
    @date 15/6/1999
*/
class PropertyConstraints {
  public:
    PropertyConstraints();
    PropertyConstraints( const PropertyConstraints& cr );
    PropertyConstraints( const char* prop_name, double64 pmin, double64 pmax );
    PropertyConstraints& operator=( const PropertyConstraints& cr );
    ~PropertyConstraints();
    bool    WithIndexes() const;
    size_t  Constraints() const;
    void    CheckLengthOfVectorVariables( bool check );
    
    // default: void SatisfyConstraintsForAllNodes( bool satisfy );
    void SatisfyConstraintsForAtLeastOneNode( bool satisfy );
    void SatisfyConstraintsForNodalAverage( bool satisfy );

    bool AddConstraint( const char* prop_name, double64 pmin, double64 pmax );
    void ChangeConstraint( const char* prop_name, double64 pmin, double64 pmax );

    template<size_t dim>
    bool InitializePropertyIndices( const PropertyDatabase<dim>& pref ); 

    template<size_t dim>
    void DeleteConstraint( const PropertyDatabase<dim>& pref, const char* prop_name );
    
    template<size_t dim>
    bool CheckConstraints( const Element<dim>& e ) const;
    
    template<size_t dim>
    bool CheckConstraints( const Element<dim>& e, Index& idx ) const;
    
    void Erase();
    void Out() const { Out(std::cout); }
    void Out(std::ostream& os) const;
  
  private:
    std::map<std::string,std::pair<double64,double64> >  criteria;
    std::map<Index,std::pair<double64,double64> >        check_list;
    bool                                                   vector_length_check;
    bool                                                   one_node_only;
    bool                                                   nodal_average;

    template<size_t dim>
    bool CheckSingleNodeConstraints( const Element<dim>& e ) const;
    
    template<size_t dim>
    bool CheckNodeAverageConstraints( const Element<dim>& e ) const;
    
    template<size_t dim>
    bool VectorLengthCheck( const Element<dim>& e, 
                            const csmp::Index& idx, double64 vmin, double64 vmax ) const;
};

<<<<<<< HEAD
=======
// templatized functions
template<size_t dim>
bool PropertyConstraints::InitializePropertyIndices( const PropertyDatabase<dim>& pref )
  {
    if ( criteria.empty() ) {
      throw csmp::Exception( CSMP_WARNING, "PropertyConstraints::InitializePropertyIndices",
        "No criteria have been defined so far");
      return false;
      }

    if ( !check_list.empty() && check_list.size() == criteria.size() ) return true; 

    std::map<std::string,std::pair<double64,double64> >::const_iterator  it;

    for ( it=criteria.begin(); it!=criteria.end(); it++ )
      check_list[ pref.StorageKey( (*it).first.c_str() ) ] =
      std::make_pair((*it).second.first,(*it).second.second);

    return true;
  }

template<size_t dim>
void PropertyConstraints::DeleteConstraint( const PropertyDatabase<dim>& pref, const char* prop_name )
  {
  criteria.erase( prop_name );
  check_list.erase( pref.StorageKey(prop_name) );
  }


>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3
} // csmp

#endif

