#ifndef VSetConnectivity_h__
#define VSetConnectivity_h__

#include "CSP_definitions.h"
#include "SuperGroup.h"
namespace csp {

template<typename fT, stl_index dim>
 class VSetConnectivity {
  public:
      VSetConnectivity();
      VSetConnectivity(std::vector<VSet<fT,dim> >& , std::vector<stl_index>& );
      ~VSetConnectivity();
      VSetConnectivity( const VSetConnectivity& );

      // function to remove halo-halo connectivity
      void RemoveHaloToHaloConnections();

      void Initialize( std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >&,
                       std::vector<stl_index>&,
                       std::vector<std::set<std::pair<stl_index,int> > >& );

      void SwitchOuterhaloDirichlet();
      void CheckOuterhaloValues(SuperGroup<csp_float,dim>&, const char*);

      // Acces functions
      void GetNodeConnectivity(std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >& nc)
       {
         nc = node_corresp_;
         return;
       }
      void GetOuterhalo(std::vector<std::set<std::pair<stl_index,int> > >& oh)
       {
         oh = outerhalo_;
         return;
       }
      void GetFirstOuterhalo(std::vector<stl_index>& foh)
       {
         foh = first_outerhalo_;
         return;
       }
      stl_index GetFirstOuterhalo(int rank)
       {
         return first_outerhalo_[rank];
       }
      stl_index CommunicationVolume(int rank, int rank2)
       {
         std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >::
         const_iterator fit=node_corresp_.find(std::make_pair(rank,rank2));
         if(fit!=node_corresp_.end())
           return (*fit).second.size();
         else
           return 0;
       }

  private:
      std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > > node_corresp_;
      std::vector<stl_index>                            first_outerhalo_;
      std::vector<std::set<std::pair<stl_index,int> > > outerhalo_;

      void FirstOuterhalo(std::vector<VSet<fT,dim> >& vsets);
      void CreateCorrespondanceMap(std::vector<VSet<fT,dim> >& vsets);
      void CreateOuterhalo(std::vector<VSet<fT,dim> >& vsets);

 };

} // END CSP NAMESPACE 
#endif
