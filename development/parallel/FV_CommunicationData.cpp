#ifndef HAVE_MPI_CXX
#define HAVE_MPI_CXX
#endif

#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

#include "mpi.h"
#include "FV_CommunicationData.h"

using namespace std;

namespace csp {

extern csp::ErrorHandler           skm_err;

template<typename fT, stl_index dim>
FV_CommunicationData<fT,dim>::~FV_CommunicationData()
 {
 }
 
template<typename fT, stl_index dim>
FV_CommunicationData<fT,dim>::FV_CommunicationData( const FV_CommunicationData& fv_comm )
: pmem(fv_comm.pmem),
  pref(fv_comm.pref)
 {
    *this = fv_comm;
 }

template<typename fT, stl_index dim>
FV_CommunicationData<fT,dim>::FV_CommunicationData( SuperGroup<fT,dim>& sg,
                                                    SAMGp_CommunicationData<fT,dim>& samgp_comm,
                                                    int& my_id )
 : pmem(sg.ReferencePropertyStorage()),
   pref(sg.ReferencePropertyDatabase()),
   tag(0), property(PLAIN,0), my_rank(my_id),
   iranksnd(samgp_comm.iranksnd),
   ipts(samgp_comm.ipts),
   isndlist(samgp_comm.isndlist),
   irankrec(samgp_comm.irankrec),
   iptr(samgp_comm.iptr),
   ireclist(samgp_comm.ireclist)
{
   cout <<"\nFV_CommunicationData<fT,dim>::FV_CommunicationData: Constructing from SAMGp_CommunicationData object." << endl;

   // reserve size in data vectors
   isnddata.resize(isndlist.size());
   irecdata.resize(ireclist.size());
   mpi_req.resize (irankrec.size());

   cout <<"\nFV_CommunicationData<fT,dim>::FV_CommunicationData: Successful construction.(rank = "<< my_rank <<")"<< endl;

}


// Constructor with input correpondance map (output of partitionVSet) and rank of processor
template<typename fT, stl_index dim>
FV_CommunicationData<fT,dim>::FV_CommunicationData( SuperGroup<fT,dim>& sg, VSetConnectivity<fT,dim>& VSetConn, int& my_id )
 : pmem(sg.ReferencePropertyStorage()),
   pref(sg.ReferencePropertyDatabase()),
   tag(0), property(PLAIN,0), my_rank(my_id)
 {
   cout <<"\nFV_CommunicationData<fT,dim>::FV_CommunicationData: Constructing from VSetConnectivity object." << endl;

   // Resize all vectors to zero
   iranksnd.resize(0);
   ipts.resize(0);
   isndlist.resize(0);
   irankrec.resize(0);
   iptr.resize(0);
   ireclist.resize(0);

   // Push_back first entry of ipts and iptr
   ipts.push_back(1);
   iptr.push_back(1);

   std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > > correspondance_map;
   std::vector<stl_index>                            first_outerhalo;
   VSetConn.GetNodeConnectivity(correspondance_map);
   VSetConn.GetFirstOuterhalo(first_outerhalo);

   stl_index partition(0);
   std::set<std::pair<stl_index,stl_index> > irec_temp;

   std::map<std::pair<stl_index,stl_index>,std::vector<std::pair<stl_index,stl_index> > >::const_iterator it;
   for (it=correspondance_map.begin();it!=correspondance_map.end();it++, partition++)
    {
      // Determining send variables: isndlist, iranksnd, ipts
      if((*it).first.first==my_rank && !(*it).second.empty()) // entry contains send-variables
          {
	    if((*it).second[0].first < first_outerhalo[my_rank])
	      iranksnd.push_back((*it).first.second);
	    if((*it).second[(*it).second.size()-1U].first >= first_outerhalo[my_rank])
	      irankrec.push_back((*it).first.second);
	    // Loop over vector in map
	    irec_temp.clear();
   	    for (unsigned int i=0;i!=(*it).second.size();i++)
	     {
   	       if ((*it).second[i].first < first_outerhalo[my_rank])
	         isndlist.push_back((*it).second[i].first);
	       else
	         irec_temp.insert(std::make_pair((*it).second[i].second,(*it).second[i].first));
	     }
  	    // Write irec_temp (2nd entry) into ireclist
	    for (std::set<std::pair<stl_index,stl_index> >::const_iterator tit=irec_temp.begin();tit!=irec_temp.end();tit++)
	      ireclist.push_back((*tit).second);

	    ipts.push_back(isndlist.size()+1);
	    iptr.push_back(ireclist.size()+1);
          }
     }

   // reserve size in data vectors
   isnddata.resize(isndlist.size());
   irecdata.resize(ireclist.size());
   mpi_req.resize (irankrec.size());

   cout<<"\nFV_CommunicationData::FV_CommunicationData: Succesfull construction (rank = "<< my_rank <<")"<< endl;

 }

template<typename fT, stl_index dim>
FV_CommunicationData<fT,dim>&  FV_CommunicationData<fT,dim>::operator=( const FV_CommunicationData& sp )
 {
    if ( &sp != this ) {
	     iranksnd  = sp.iranksnd;
	     ipts      = sp.ipts;
	     isndlist  = sp.isndlist;
	     irankrec  = sp.irankrec;
	     iptr      = sp.iptr;
	     ireclist  = sp.ireclist;
      }
    return *this;
 }

template<typename fT, stl_index dim>
void FV_CommunicationData<fT,dim>::UpdateHaloNodes(const char* prop, bool blocking_send_recv, bool verbose )
 {
    cout<<"\nFV_CommunicationData<fT,dim>::UpdateHaloNodes: Rank "<<MPI::COMM_WORLD.Get_rank();
    if(blocking_send_recv) cout << " BLOCKING ";
    else                  cout<< " NON-BLOCKING ";
    cout<<" Send/Recv outerhalo "<<prop<<"-data."<< endl;


    Index prop_key = pref.StorageKey(prop);

    // Unique tag for communication step
    tag++;
    if(tag==1000)tag=0;

    // A Non-blocking MPI send-receive call
    // This method should be faster
    if(!blocking_send_recv)
     {
       // All process post non-blocking receive to all neighbors
       for (unsigned int i=0;i!=irankrec.size();i++)
        {
         mpi_req[i] = MPI::COMM_WORLD.Irecv(&irecdata[iptr[i]-1],(iptr[i+1]-iptr[i]),MPI::DOUBLE,irankrec[i],tag);
        }

       // Now update data vector with data to be send:
       for (unsigned int i=0;i!=isndlist.size();i++)
         isnddata[i]=pmem.Read(isndlist[i],NODE,prop_key.index);

       // All process post non-blocking send to all neighbors
       for (unsigned int i=0;i!=iranksnd.size();i++)
        {
         MPI::COMM_WORLD.Isend(&isnddata[ipts[i]-1],(ipts[i+1]-ipts[i]),MPI::DOUBLE,iranksnd[i],tag);
         if(verbose) cout<<"\nFV_CommunicationData<fT,dim>::UpdateHaloNodes: Rank "<<MPI::COMM_WORLD.Get_rank()<<" send "<<(ipts[i+1]-ipts[i])<<" halo '"<<prop<<"' values to "<<iranksnd[i]<< endl;
        }

       // Check if irecdata is received
       for (unsigned int i=0;i!=mpi_req.size();i++)
        {
          mpi_req[i].Wait();
          if(verbose) cout<<"\nFV_CommunicationData<fT,dim>::UpdateHaloNodes: Rank "<<MPI::COMM_WORLD.Get_rank()<<" received "<<(iptr[i+1]-iptr[i])<<" halo '"<<prop<<"' values from "<<irankrec[i]<< endl;
        }
     }
    // A Blocking MPI send-receive call
    // This method will be slower
    else
     {
       // First update data vector with data to be send:
       for (unsigned int i=0;i!=isndlist.size();i++)
         isnddata[i]=pmem.Read(isndlist[i],NODE,prop_key.index);

       // Loop over number of procs
       for (unsigned int send_proc=0;send_proc!=MPI::COMM_WORLD.Get_size();send_proc++)
        {
          // Determine receive proc:
          if (send_proc!=MPI::COMM_WORLD.Get_rank())
           {
            for (unsigned int i=0;i!=irankrec.size();i++)
             if(irankrec[i]==send_proc)
              {
              MPI::COMM_WORLD.Recv(&irecdata[iptr[i]-1],(iptr[i+1]-iptr[i]),MPI::DOUBLE,send_proc,tag);
              if(verbose)cout<<"\nFV_CommunicationData<fT,dim>::UpdateHaloNodes: Rank "<<MPI::COMM_WORLD.Get_rank()<<" received "<<(iptr[i+1]-iptr[i])<<" halo '"<<prop<<"' values from "<<send_proc<< endl;
              }
           }
          else
           {
            for (unsigned int i=0;i!=iranksnd.size();i++)
             {
              MPI::COMM_WORLD.Send(&isnddata[ipts[i]-1],(ipts[i+1]-ipts[i]),MPI::DOUBLE,iranksnd[i],tag);
              if(verbose)cout<<"\nFV_CommunicationData<fT,dim>::UpdateHaloNodes: Rank "<<MPI::COMM_WORLD.Get_rank()<<" send "<<(ipts[i+1]-ipts[i])<<" halo '"<<prop<<"' values to "<<iranksnd[i]<< endl;
             }
           }
        }
     }

    // Now update variables with received data
    for (unsigned int i=0;i!=irecdata.size();i++)
     {
       property = irecdata[i];
       pmem.Store(ireclist[i],NODE,prop_key.index,property);
     }
 }


template<typename fT, stl_index dim>
void FV_CommunicationData<fT,dim>::ConstantValueOnOuterhalo(const char* prop, csp_float val )
 {
    cout<<"\nFV_CommunicationData<fT,dim>::ConstantValueOnOuterhalo: Rank "<<MPI::COMM_WORLD.Get_rank()<<" Storing constant value of '"<<prop<<"' on outerhalo."<< endl;
    Index prop_key = pref.StorageKey(prop);

    // Now update variables with received data
    for (unsigned int i=0;i!=ireclist.size();i++)
     {
       property = val;
       pmem.Store(ireclist[i],NODE,prop_key.index,property);
     }

 }

template<typename fT, stl_index dim>
void FV_CommunicationData<fT,dim>::ConstantValueOnInnerhalo(const char* prop, csp_float val )
 {
    cout<<"\nFV_CommunicationData<fT,dim>::ConstantValueOnOuterhalo: Rank "<<MPI::COMM_WORLD.Get_rank()<<" Storing constant value of '"<<prop<<"' on outerhalo."<< endl;
    Index prop_key = pref.StorageKey(prop);

    // Now update variables with received data
    for (unsigned int i=0;i!=isndlist.size();i++)
     {
       property = val;
       pmem.Store(isndlist[i],NODE,prop_key.index,property);
     }

 }


template<typename fT, stl_index dim>
void FV_CommunicationData<fT,dim>::SetOuterhaloDirich(const char* prop)
 {
   cout<<"\nFV_CommunicationData<fT,dim>::SetOuterhaloDirich: Rank "<<MPI::COMM_WORLD.Get_rank()<<" flagged "<<ireclist.size()<<" "<<prop<<" variables"<< endl;

   Index prop_key = pref.StorageKey(prop);
   // loop over outerhalo node id's
   for(std::vector<int>::const_iterator it=ireclist.begin();it!=ireclist.end();it++)
      pmem.Status((*it),prop_key,DIRICH);
 }

template<typename fT, stl_index dim>
void FV_CommunicationData<fT,dim>::SetOuterhaloPlain(const char* prop)
 {
   cout<<"\nFV_CommunicationData<fT,dim>::SetOuterhaloPlain: Rank "<<MPI::COMM_WORLD.Get_rank()<<" flagged "<<ireclist.size()<<" "<<prop<<" variables"<< endl;
   Index prop_key = pref.StorageKey(prop);
   // loop over outerhalo node id's
   for(std::vector<int>::const_iterator it=ireclist.begin();it!=ireclist.end();it++)
      pmem.Status((*it),prop_key,PLAIN);
 }

template<typename fT, stl_index dim>
void FV_CommunicationData<fT,dim>::Out()
 {
    // Output to screen
    cout <<"\nSAMGp_CommunicationData::Out(): Rank "<<my_rank<<": Dumping data to screen:"<< endl;

    cout <<"\nRank "<<my_rank<<": SEND VARIABLES:"<< endl;
    cout <<"\nFV_CommunicationData::nshalo = "<<nshalo()<< endl;
    cout <<"\nFV_CommunicationData::isndlist: "<< endl;
    for(unsigned int i=0;i!=isndlist.size();i++) cout <<"\t"<< isndlist[i];
    cout <<"\nFV_CommunicationData::npsnd = "<<npsnd()<< endl;
    cout <<"\nFV_CommunicationData::iranksnd: "<< endl;
    for(unsigned int i=0;i!=iranksnd.size();i++) cout <<"\t"<< iranksnd[i];
    cout <<"\nFV_CommunicationData::ipts: "<< endl;
    for(unsigned int i=0;i!=ipts.size();i++) cout <<"\t"<< ipts[i];

    cout <<"\nRank "<<my_rank<<": RECEIVE VARIABLES:"<< endl;
    cout <<"\nFV_CommunicationData::nrhalo = "<<nrhalo()<< endl;
    cout <<"\nFV_CommunicationData::ireclist: "<< endl;
    for(unsigned int i=0;i!=ireclist.size();i++) cout <<"\t"<< ireclist[i];
    cout <<"\nFV_CommunicationData::nprec = "<<nprec()<< endl;
    cout <<"\nFV_CommunicationData::irankrec: "<< endl;
    for(unsigned int i=0;i!=irankrec.size();i++) cout <<"\t"<< irankrec[i];
    cout <<"\nFV_CommunicationData::iptr: "<< endl;
    for(unsigned int i=0;i!=iptr.size();i++) cout <<"\t"<< iptr[i];

 }


template class FV_CommunicationData<csp_float,1U>;
template class FV_CommunicationData<csp_float,2U>;
template class FV_CommunicationData<csp_float,3U>;


} // END CSP NAMESPACE
/*
COMMENT

iranksnd - vector of integers
         - contains rankids of processors to which data should be send
         - iranksnd.size()==npsnd

ipts - vector of integers
     - needed to determine which entries of isndlist have to be send to which processor
     - ipts.size() == iranksnd.size()+1
     - first entry is always one, next entry is i+1, with isndlist[i] the first element which
       need to be send to the next negihboring processor, etc..

isndlist - vector of integers
         - contains column numbers of entries outside the partition, but in the same column range,
           which contain non-zeros. For each neighboring processor only one column number needs to 
           be stored. As the same column could be send to several neighboring processors, this list
           can have several entries with the same value. The local column numbers are stored, so it
           should only contain numbers between the 1 and nnu_local

irankrec - vector of integers
         - contains rankids of processors from which data should be received
         - iranksnd.size()==nprec 

iptr - vector of integers
     - needed to determine which entries of ireclist have to be send to which processor
     - iptr.size() == irankrec.size()+1
     - first entry is always one, next entry is i+1, with ireclist[i] the first element which
       need to be received from the next negihboring processor, etc..

ireclist - vector of integers
         - contains column numbers of entries outside the partition, but in the same row range,
           which contain non-zeros. For each neighboring processor only one column number needs to 
           be stored. Since a column number can only relate to one neighboring processor, this list 
           can NOT have several entries with the same value. Local column numbers are stored. These
           are counted, starting from nnu_local+1 contiguitively. So this list will always have 
           opeenvolgende nummers.
*/
