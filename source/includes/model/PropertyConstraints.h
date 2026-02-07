#ifndef PROPERTY_CONSTRAINTS_H
#define PROPERTY_CONSTRAINTS_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class Node;
template<uint32_t> class PropertyDatabase;

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

    template<uint32_t dim>
    PropertyConstraints( const PropertyDatabase<dim>&, const char* prop_name, double pmin, double pmax );

    PropertyConstraints( const PropertyConstraints& cr );
    ~PropertyConstraints() = default;
    
    PropertyConstraints& operator=( const PropertyConstraints& cr );
    
    bool     WithIndexes() const;
    uint32_t Constraints() const;
    void     CheckLengthOfVectorVariables( bool check );
    
    bool AddConstraint( const char* prop_name, double pmin, double pmax );
    void ChangeConstraint( const char* prop_name, double pmin, double pmax );

    template<uint32_t dim>
    bool InitializePropertyIndices( const PropertyDatabase<dim>& );

    template<uint32_t dim>
    void DeleteConstraint( const PropertyDatabase<dim>&, const char* prop_name );
    
    template<uint32_t dim, template<uint32_t> class CELL>
    bool CheckConstraints( const CELL<dim>* ) const;
    
    template<uint32_t dim, template<uint32_t> class CELL>
    bool CheckConstraints( const CELL<dim>* , Index& idx ) const;
    
    // default: void SatisfyConstraintsForAllNodes( bool satisfy );
    void SatisfyConstraintsForAtLeastOneNode( bool satisfy );
    void SatisfyConstraintsForNodalAverage( bool satisfy );

    void Erase();
    
    void Out() const;
  
  private:
    template<uint32_t dim, template<uint32_t> class CELL>
    bool CheckSingleNodeConstraints( const CELL<dim>* ) const;
    
    template<uint32_t dim, template<uint32_t> class CELL>
    bool CheckNodeAverageConstraints( const CELL<dim>* ) const;
    
    template<uint32_t dim, template<uint32_t> class CELL>
    bool VectorLengthCheck( const CELL<dim>*,
                            const csmp::Index&, double vmin, double vmax ) const;

    std::map<std::string,std::pair<double,double> >  criteria_;
    std::map<Index,std::pair<double,double> >        check_list_;
    bool                                             vector_length_check_ = true;  ///< for VectorVariable uses length as constraint
    bool                                             one_node_only_       = false; ///< will include elements even if only a single node falls into range
    bool                                             nodal_average_       = false; ///< will use the average of node properties to apply constraints
};

} // csmp

#endif

