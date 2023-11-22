#ifndef CSP_ParallelSuperGroup_H
#define CSP_ParallelSuperGroup_H

#include "SuperGroup2.h"
#include "CSP_VSet.h"

namespace csp {

template<typename fT,stl_index dim>
class ParallelSuperGroup : public SuperGroup<fT,dim> {

  private:
   // MPI variables
   int n_procs, my_rank;

   // Added member Functions
   void FlagCornerElements();
   bool FlagIrregularTopBoundary( fT right_left_tolerance );

    std::map<std::string,Group<fT,dim> >  group_list;
    PropertyDatabase                      phys_vars;
    MemoryManager<fT,dim>                 property_collection;
    MeshManager<fT,dim>                   local_mesh;
    FiniteElementManager                  fem_manager;

  public:

    // Constructors
    ~ParallelSuperGroup();
    ParallelSuperGroup( const ParallelSuperGroup<fT,dim>& sg );
    explicit ParallelSuperGroup( VSet<fT,dim>& vset, const char* var_file, bool isoparametric=false ); 
    ParallelSuperGroup& operator=( const ParallelSuperGroup& sg );

    bool EnlistBoundaryElements( SG_BOUNDARY side, std::list<stl_index>& belmts   ) const;
    void AddCornersTo( SG_BOUNDARY side, PLACEMENT place, std::list<stl_index>& bdata ) const;
    void AssignBoundaryValuesParallel( SG_BOUNDARY b, const char* prop, DATA_STYLE bcond, fT bmin, fT bmax  );

};

} // end namespace csp


#endif 



