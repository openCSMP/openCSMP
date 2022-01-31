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
    PropertyConstraints( const char* prop_name, double pmin, double pmax );
    PropertyConstraints& operator=( const PropertyConstraints& cr );
    ~PropertyConstraints();
    bool    WithIndexes() const;
    size_t  Constraints() const;
    void    CheckLengthOfVectorVariables( bool check );
    
    // default: void SatisfyConstraintsForAllNodes( bool satisfy );
    void SatisfyConstraintsForAtLeastOneNode( bool satisfy );
    void SatisfyConstraintsForNodalAverage( bool satisfy );

    bool AddConstraint( const char* prop_name, double pmin, double pmax );
    void ChangeConstraint( const char* prop_name, double pmin, double pmax );

    template<size_t dim>
    bool InitializePropertyIndices( const PropertyDatabase<dim>& pref ); 

    template<size_t dim>
    void DeleteConstraint( const PropertyDatabase<dim>& pref, const char* prop_name );
    
    template<size_t dim>
    bool CheckConstraints( const Element<dim>& e ) const;
    
    template<size_t dim>
    bool CheckConstraints( const Element<dim>& e, Index& idx ) const;
    
    void Erase();
    void Out() const;
  
  private:
    std::map<std::string,std::pair<double,double> >  criteria;
    std::map<Index,std::pair<double,double> >        check_list;
    bool                                                   vector_length_check;
    bool                                                   one_node_only;
    bool                                                   nodal_average;

    template<size_t dim>
    bool CheckSingleNodeConstraints( const Element<dim>& e ) const;
    
    template<size_t dim>
    bool CheckNodeAverageConstraints( const Element<dim>& e ) const;
    
    template<size_t dim>
    bool VectorLengthCheck( const Element<dim>& e, 
                            const csmp::Index& idx, double vmin, double vmax ) const;
};

} // csmp

#endif

