#ifndef Partitioner__
#define Partitioner__

#include "CSP_definitions.h"
#include "CSP_VSet.h"
 

namespace csp {

template<typename fT, stl_index dim>
class Partitioner
  {
    public:
    
    // Constructor:
    Partitioner ( const VSet<fT,dim>& vs, int n );
	
    // Partition functions
    void MetisPartition( std::vector<std::vector<stl_index> >& element_vector, bool element_partition=true );
    void GeometricPartition( std::vector<std::vector<stl_index> > element_vector ) const;
    void FastPartition( std::vector<std::vector<stl_index> >& element_vector, bool three_D = false );

    // Access functions
    void GetSharedNodes( unsigned int proc, std::vector<std::pair<int,int> >& shared_n );
    void GetHaloElements( unsigned int proc, std::vector<int>& halo_el_eids );
    void GetHaloNodes(std::vector<std::set<std::pair<stl_index,int> > >&  inhalo, std::vector<std::set<std::pair<stl_index,int> > >&  outhalo );
    void GetOuterhaloNodes(std::vector<std::set<std::pair<stl_index,int> > >&  outhalo );

    private:
    void SharedNodes();
    void AddHaloElements(std::vector<std::vector<stl_index> >& element_vector);
    void ErroneousElements(std::vector<std::vector<stl_index> >& element_vector);
    void RemoveHaloElements(std::vector<std::vector<stl_index> >& element_vector);
    void CreateHaloElements(const std::vector<std::vector<stl_index> >& element_vector);
    void CreateHaloNodes();
    void CreateNodeVector(std::vector<std::vector<stl_index> >& element_vector);
    void RemoveDoubleEntries(std::vector<std::vector<stl_index> >& element_vector);
    void RemoveDoubleOuterhaloEntries();
    void PartitionInfoToScreen(std::vector<std::vector<stl_index> >& element_vector);

    VSet<fT,dim>	vset;
    int n_parts, nn, ne;
    std::vector<int> epart, npart;
    std::vector<std::set<std::pair<stl_index,int> > >  halo_elements;
    std::vector<std::set<stl_index> >                  node_vector;
    std::vector<std::set<std::pair<stl_index,int> > >  shared_nodes;
    std::vector<std::set<std::pair<stl_index,int> > >  innerhalo;
    std::vector<std::set<std::pair<stl_index,int> > >  outerhalo;
   }; // end class

} // csp

#endif
