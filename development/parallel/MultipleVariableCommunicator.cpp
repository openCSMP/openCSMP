#ifndef HAVE_MPI_CXX
#define HAVE_MPI_CXX
#endif

#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif


#include "mpi.h"
#include "MultipleVariableCommunicator.h"

using namespace std;

namespace csp {

extern csp::ErrorHandler           skm_err;

template<typename fT, stl_index dim>
MultipleVariableCommunicator<fT,dim>::~MultipleVariableCommunicator()
 {
 }

template<typename fT, stl_index dim>
MultipleVariableCommunicator<fT,dim>::MultipleVariableCommunicator( SuperGroup<fT,dim>& sg,
                                                            FV_CommunicationData<fT,dim>& FVC,
                                                            std::vector<const char*>& props )
 : pmem(sg.ReferencePropertyStorage()),
   pref(sg.ReferencePropertyDatabase()),
   FV_comm (FVC),
   my_rank(MPI::COMM_WORLD.Get_rank()),
   properties(props),
   tag(0)
{
   // reserve sizes in data vectors
   isnddata.resize( FV_comm.iranksnd.size() );
   irecdata.resize( FV_comm.irankrec.size() );
   for (unsigned int i=0;i!=isnddata.size();i++)
     isnddata[i].resize( (FV_comm.ipts[i+1]-FV_comm.ipts[i]) * props.size() );
   for (unsigned int i=0;i!=irecdata.size();i++)
     irecdata[i].resize( (FV_comm.iptr[i+1]-FV_comm.iptr[i]) * props.size() );

   cout <<"\nMultipleVariableCommunicator::MultipleVariableCommunicator: Rank = "<< my_rank <<": Created communicator for properies:\n"<< endl;
   for(std::vector<const char*>::const_iterator it=properties.begin();it!=properties.end();it++)
     cout <<"\t'"<<(*it)<<"'"<< endl;
}

template<typename fT, stl_index dim>
void MultipleVariableCommunicator<fT,dim>::UpdateHaloNodes()
 {

    csp_float start = MPI::Wtime();

    // Unique tag for communication step
    tag++;
    if(tag==1000)tag=0;

    // First update data vector with data to be send:
    unsigned int i(0);
    unsigned int k(0);
    for (std::vector<std::vector<csp_float> >::iterator it=isnddata.begin();it!=isnddata.end();it++,i++)
     {
       k = 0;
       for (unsigned int j=0;j!=properties.size();j++)
        {
          prop_key = pref.StorageKey(properties[j]);
          for (unsigned int l = (FV_comm.ipts[i]-1); l!=(FV_comm.ipts[i+1]-1); l++, k++)
            (*it)[k] = pmem.Read(FV_comm.isndlist[l],NODE,prop_key.index);
        }
     }
/*
    i=0;
    for (std::vector<std::vector<csp_float> >::iterator it=isnddata.begin();it!=isnddata.end();it++,i++)
     {
       cout<<"\nisnddata ["<<i<<"]: for processor :"<<FV_comm.iranksnd[i]<< endl;
       for(unsigned int j=0;j!=(*it).size();j++)
         cout<<j<<":"<< (*it)[j]<< endl;
     }
*/
    // Now send-recv data vectors
    for (unsigned int send_proc=0;send_proc!=MPI::COMM_WORLD.Get_size();send_proc++)
     {
       // Determine receive proc:
       if (send_proc!=MPI::COMM_WORLD.Get_rank())
        {
         for (unsigned int i=0;i!=FV_comm.irankrec.size();i++)
          if(FV_comm.irankrec[i]==send_proc)
           {
            MPI::COMM_WORLD.Recv(&irecdata[i][0],irecdata[i].size(),MPI::DOUBLE,send_proc,tag);
           }
        }
       else
        {
         for (unsigned int i=0;i!=FV_comm.iranksnd.size();i++)
          {
           MPI::COMM_WORLD.Send(&isnddata[i][0],isnddata[i].size(),MPI::DOUBLE,FV_comm.iranksnd[i],tag);
          }
        }
     }
/*
    i=0;
    for (std::vector<std::vector<csp_float> >::iterator it=irecdata.begin();it!=irecdata.end();it++,i++)
     {
       cout<<"\nirecdata ["<<i<<"]: from processor :"<<FV_comm.irankrec[i]<< endl;
       for(unsigned int j=0;j!=(*it).size();j++)
         cout<<j<<":"<< (*it)[j]<< endl;
     }
*/

    // Now update variables with received data
    i=0;
    for (std::vector<std::vector<csp_float> >::iterator it=irecdata.begin();it!=irecdata.end();it++,i++)
     {
       k=0;
       for (unsigned int j=0;j!=properties.size();j++)
        {
         prop_key = pref.StorageKey(properties[j]);
          for (unsigned int l = (FV_comm.iptr[i]-1); l!=(FV_comm.iptr[i+1]-1); l++, k++)
           {
            property = (*it)[k];
            pmem.Store(FV_comm.ireclist[l],NODE,prop_key.index,property);
           }
        }
     }

    cout<<"\nMultipleVariableCommunicator::UpdateHaloNodes: Communicated all variables in "<<MPI::Wtime()-start<<" s."<< endl;

 }

template class MultipleVariableCommunicator<csp_float,1U>;
template class MultipleVariableCommunicator<csp_float,2U>;
template class MultipleVariableCommunicator<csp_float,3U>;


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
