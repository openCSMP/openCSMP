#ifndef HAVE_MPI_CXX
#define HAVE_MPI_CXX
#endif

#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

#include "mpi.h"
#include "VSetConnectivity.h"

using namespace std;

namespace csp {

template<typename fT, stl_index dim>
VSetConnectivity<fT,dim>::VSetConnectivity()
 {
 }

template<typename fT, stl_index dim>
VSetConnectivity<fT,dim>::VSetConnectivity(std::vector<VSet<fT,dim> >& vsets, std::vector<stl_index>& frstOH )
: first_outerhalo_(frstOH)
 {
  CreateCorrespondanceMap(vsets);
  RemoveHaloToHaloConnections();
//  FirstOuterhalo(vsets);
  CreateOuterhalo(vsets);

  cout<<"\nVSetConnectivity<fT,dim>::VSetConnectivity: Connectivities between VSets determined succesfully."<< endl;
 }

template<typename fT, stl_index dim>
VSetConnectivity<fT,dim>::~VSetConnectivity()
 {
 }

template<typename fT, stl_index dim>
VSetConnectivity<fT,dim>::VSetConnectivity( const VSetConnectivity& sp )
 {
    *this = sp;
 }

template<typename fT, stl_index dim>
void VSetConnectivity<fT,dim>::CreateOuterhalo(std::vector<VSet<fT,dim> >& vsets)
 {
   outerhalo_.resize(vsets.size());
   for(stl_index proc=0;proc!=vsets.size();proc++)
    {
     outerhalo_[proc].clear();
     for(stl_index n_id=first_outerhalo_[proc];n_id!=(vsets[proc].Vertices()+1U);n_id++)
      for(stl_index proc2=0;proc2!=vsets.size();proc2++)
       if(proc2!=proc)
       {
         std::pair<stl_index,stl_index> pr;
         pr.first=proc;
         pr.second=proc2;
         std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >
           ::const_iterator mit = node_corresp_.find(pr);
         std::vector<stl_index> temp;
         temp.resize((*mit).second.size());
         for(unsigned int i=0;i!=temp.size();i++)
          temp[i]=(*mit).second[i].first;
         std::vector<stl_index>::const_iterator fit=find(temp.begin(),temp.end(),n_id);
         if(fit!=temp.end())
          {
            std::pair<stl_index,stl_index> temp2;
            temp2.first = n_id;
            temp2.second= proc2;
            outerhalo_[proc].insert(temp2);
            break;
          }
       }
    }
 }

template<typename fT, stl_index dim>
void VSetConnectivity<fT,dim>::CreateCorrespondanceMap(std::vector<VSet<fT,dim> >& vsets)
 {

int my_rank = MPI::COMM_WORLD.Get_rank();

  std::vector<std::pair<stl_index,stl_index> > corresp_nodes;
  // Loop over all vsets
  for(stl_index proc=0;proc!=vsets.size();proc++)
   for(stl_index proc2=0;proc2!=vsets.size();proc2++)
    if(proc2!=proc)
     {
      corresp_nodes.clear();
      for(stl_index n=0;n!=vsets[proc].Vertices();n++)
       for(stl_index n2=0;n2!=vsets[proc2].Vertices();n2++)
        if(vsets[proc].Px(n) == vsets[proc2].Px(n2)) // if same x
         {
          if(dim == 1U) corresp_nodes.push_back(std::make_pair((n+1U),(n2+1U)));
          else if (vsets[proc].Py(n)== vsets[proc2].Py(n2))
           {
             if(dim == 2U) corresp_nodes.push_back(std::make_pair((n+1U),(n2+1U)));
             else if (vsets[proc].Pz(n)== vsets[proc2].Pz(n2))
                corresp_nodes.push_back(std::make_pair((n+1U),(n2+1U)));
           }
         }
      if (!corresp_nodes.empty()) node_corresp_.insert(std::make_pair(std::make_pair(proc,proc2),corresp_nodes));
     }
 }


template<typename fT, stl_index dim>
void VSetConnectivity<fT,dim>::FirstOuterhalo(std::vector<VSet<fT,dim> >& vsets)
 {
   first_outerhalo_.resize(vsets.size());
   // Loop over VSets
   for(unsigned int proc=0;proc!=vsets.size();proc++)
    {
      first_outerhalo_[proc]=vsets[proc].Vertices()+1;
      for (std::map<stl_index,int32>::const_iterator bf=vsets[proc].BFlagsBegin(); bf!=vsets[proc].BFlagsEnd(); bf++)
       if((*bf).second==-1 && (*bf).first < first_outerhalo_[proc])// if IRREGULAR_OUTSIDE
        {
          // find in corresp_map
          for (std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >::iterator
            it=node_corresp_.begin();it!=node_corresp_.end();it++)
            if((*it).first.first==proc)
             {
               // write 1st entry of vector of pairs to temporary vector
               std::vector<stl_index> temp;
               temp.resize((*it).second.size());
               for(unsigned int j=0;j!=(*it).second.size();j++)
                 temp[j]=(*it).second[j].first;
               // find if ID is in temp
               std::vector<stl_index>::const_iterator fit=find(temp.begin(),temp.end(),(*bf).first);
               if(fit!=temp.end())
                 first_outerhalo_[proc]=(*bf).first;
               // this should work as innerhalo can never be flagged IRREGULAR OUTSIDE
             }
        }
    }
 }

template<typename fT, stl_index dim>
void VSetConnectivity<fT,dim>::Initialize ( std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >& c_map,
                                    std::vector<stl_index>& first_r_halo,
                                    std::vector<std::set<std::pair<stl_index,int> > >& r_halo)
 {
    if(!node_corresp_.empty())    node_corresp_.clear();
    if(!first_outerhalo_.empty()) first_outerhalo_.clear();
    if(!outerhalo_.empty())       outerhalo_.clear();

    node_corresp_ = c_map;
    first_outerhalo_ = first_r_halo;
    outerhalo_ = r_halo;

    cout<<"\nVSetConnectivity::Initialize: Succesfull initialization"<< endl;
 }

template<typename fT, stl_index dim>
void VSetConnectivity<fT,dim>::RemoveHaloToHaloConnections()
 {
   // Adapt correpondance map such that outerhalo-outerhalo entries are removed
   // OR: entries with the local node being a outerhalo node and the
   // neighbornode also being a outerhalo node need to be removed
   for (std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >::iterator
        it=node_corresp_.begin();it!=node_corresp_.end();it++)
    for (std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >::iterator
         it2=node_corresp_.begin();it2!=node_corresp_.end();it2++)
      if((*it2).first.first==(*it).first.second && (*it2).first.second==(*it).first.first) // if neighbours
        //loop over shared node vector
   	for (std::vector<std::pair<stl_index,stl_index> >::iterator i=(*it).second.begin();i!=(*it).second.end();i++) 
  	  if((*i).first >= first_outerhalo_[(*it).first.first])  
          // (*i) is outerhalo node of partition (*it).first.first
	   {
	     bool wrong_entry(true);
	     for (std::vector<std::pair<stl_index,stl_index> >::iterator j=(*it2).second.begin();j!=(*it2).second.end();j++)
	       if((*j).first <  first_outerhalo_[(*it2).first.first] && (*i).first == (*j).second)  
               // (*j).first is innerhalo node of partition (*it2).first.first
               // so this connectivity should be kept!
		  {
		    wrong_entry = false;
		    break;
		  }
	     if(wrong_entry)
	      {
cout<<"\nVSetConnectivity::RemoveHaloToHaloConnections: Node "<<(*i).first<<" is outerhalo of partition "<<(*it).first.first<<\
" and partition "<<(*it).first.second<<"..Removing connectivity."<< endl;
	       (*it).second.erase(i);
	       i--;
	      }
	   }

 }

template<typename fT, stl_index dim>
void VSetConnectivity<fT,dim>::CheckOuterhaloValues(SuperGroup<csp_float,dim>& sg, const char* prop)
 {
   // this checks if the value of "prop" on the outerhalo of 
   // supergroup "sg" correspond with values on the innerhalo
   // of its neighbors

   // This is a master/slave setup: needs to be called in PARALLEL
   // also for each node a MPI communication step is performed
   // This is thus CPU INTENSIVE!
   // Program aborts when difference in value is found!

   const PropertyDatabase&             pref  = sg.ReferencePropertyDatabase();
   const MemoryManager<csp_float,dim>& pmem  = sg.ReferencePropertyStorage();
   const MeshManager<csp_float,dim>&   pmesh = sg.ReferenceMesh();
   Index prop_key                            = pref.StorageKey(prop);

   // MPI variables
   int n_procs = MPI::COMM_WORLD.Get_size();
   int my_rank = MPI::COMM_WORLD.Get_rank();
   int tag(0);
   int master, slave;
   std::vector<csp_float> local_vector, remote_vector;

   for (std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >::const_iterator
         it=node_corresp_.begin();it!=node_corresp_.end();it++)
    {
      master = (*it).first.first;
      slave  = (*it).first.second;
      if(my_rank==master)
       {
         // resize vectors
         local_vector.resize((*it).second.size());
         remote_vector.resize((*it).second.size());

         // receive data from slave
         MPI::COMM_WORLD.Recv(&remote_vector[0],remote_vector.size(),MPI::DOUBLE,(*it).first.second,tag);

         // create local data vector
         for (unsigned int j=0;j!=(*it).second.size();j++)
           local_vector[j] = pmem.Read((*it).second[j].first,NODE,prop_key.index);

         assert(local_vector==remote_vector);
         for (unsigned int j=0;j!=local_vector.size();j++)
          if(local_vector[j]!=remote_vector[j])
           {
             cout<<"\n\nRank: "<<my_rank<<": VSetConnectivity::CheckOuterhaloValues "<< endl;
             cout<<"Node "<<(*it).second[j].first<<" has incorrect "<<prop<<" defined on boundary:"<< endl;
             cout<<"Partition: "<<(*it).first.first<<": value = "<<local_vector[j]<< endl;
             pmesh.N((*it).second[j].first-1).Out();
             cout<<"Partition: "<<(*it).first.second<<": value = "<<remote_vector[j]<< endl;
           }
       }
      else if (my_rank==slave)
       {
         // resize vector
         local_vector.resize((*it).second.size());

         // create local data vector
         for (unsigned int j=0;j!=(*it).second.size();j++)
           local_vector[j] = pmem.Read((*it).second[j].second,NODE,prop_key.index);

         // send data to master
         MPI::COMM_WORLD.Send(&local_vector[0],local_vector.size(),MPI::DOUBLE,master,tag);
       }
    }
 }


template class VSetConnectivity<csp_float,1U>;
template class VSetConnectivity<csp_float,2U>;
template class VSetConnectivity<csp_float,3U>;

} // END CSP NAMESPACE
