#include "SAMGp_CommunicationData.h"

using namespace std;

namespace csp {

extern csp::ErrorHandler           skm_err;

template<typename fT, stl_index dim>
SAMGp_CommunicationData<fT,dim>::SAMGp_CommunicationData()
 {
 }
 
template<typename fT, stl_index dim>
SAMGp_CommunicationData<fT,dim>::~SAMGp_CommunicationData()
 {
 }
 
template<typename fT, stl_index dim>
SAMGp_CommunicationData<fT,dim>::SAMGp_CommunicationData( const SAMGp_CommunicationData& sp )
 {
    *this = sp;
 }

// Constructor with input correpondance map (output of partitionVSet) and rank of processor
template<typename fT, stl_index dim>
SAMGp_CommunicationData<fT,dim>::SAMGp_CommunicationData( VSetConnectivity<fT,dim>& VSetConn, int& my_rank )
: my_rank_(my_rank)
 {
   cout<<"\nSAMGp_CommunicationData::SAMGp_CommunicationData: Constructing from VSetConnectivity object (rank = "<< my_rank <<")"<< endl;

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

   cout<<"\nSAMGp_CommunicationData::SAMGp_CommunicationData: Succesfull construction (rank = "<< my_rank <<")"<< endl;

 }


template<typename fT, stl_index dim>
SAMGp_CommunicationData<fT,dim>&  SAMGp_CommunicationData<fT,dim>::operator=( const SAMGp_CommunicationData& sp )
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
void SAMGp_CommunicationData<fT,dim>::Out()
 {
    // Output to screen
    cout <<"\nSAMGp_CommunicationData::Out(): Rank "<<my_rank_<<": Dumping data to screen:"<< endl;

    cout <<"\nRank "<<my_rank_<<": SEND VARIABLES:"<< endl;
    cout <<"\nSAMGp_CommunicationData::nshalo = "<<nshalo()<< endl;
    cout <<"\nSAMGp_CommunicationData::isndlist: "<< endl;
    for(unsigned int i=0;i!=isndlist.size();i++) cout <<"\t"<< isndlist[i];
    cout <<"\nSAMGp_CommunicationData::npsnd = "<<npsnd()<< endl;
    cout <<"\nSAMGp_CommunicationData::iranksnd: "<< endl;
    for(unsigned int i=0;i!=iranksnd.size();i++) cout <<"\t"<< iranksnd[i];
    cout <<"\nSAMGp_CommunicationData::ipts: "<< endl;
    for(unsigned int i=0;i!=ipts.size();i++) cout <<"\t"<< ipts[i];

    cout <<"\nRank "<<my_rank_<<": RECEIVE VARIABLES:"<< endl;
    cout <<"\nSAMGp_CommunicationData::nrhalo = "<<nrhalo()<< endl;
    cout <<"\nSAMGp_CommunicationData::ireclist: "<< endl;
    for(unsigned int i=0;i!=ireclist.size();i++) cout <<"\t"<< ireclist[i];
    cout <<"\nSAMGp_CommunicationData::nprec = "<<nprec()<< endl;
    cout <<"\nSAMGp_CommunicationData::irankrec: "<< endl;
    for(unsigned int i=0;i!=irankrec.size();i++) cout <<"\t"<< irankrec[i];
    cout <<"\nSAMGp_CommunicationData::iptr: "<< endl;
    for(unsigned int i=0;i!=iptr.size();i++) cout <<"\t"<< iptr[i];

 }


template<typename fT, stl_index dim>
void SAMGp_CommunicationData<fT,dim>::DumpToTextFile(const char* fname)
 {
   cout<<"\nSAMGp_CommunicationData::DumpToTextFile: Dumping communication data to text file "<<fname<< endl;

   std::ofstream ofs(fname);
   ofs <<"# This file contains data required for communication for the SAMGp solver" << endl;
   ofs <<"# Do NOT alter this file.";

   ofs <<"\n# iranksnd:";
   ofs <<"\n"<<iranksnd.size();
   for (unsigned int i=0;i!=iranksnd.size();i++)
     ofs <<"\n"<<iranksnd[i];

   ofs <<"\n# ipts:";
   ofs <<"\n"<<ipts.size();
   for (unsigned int i=0;i!=ipts.size();i++)
     ofs <<"\n"<<ipts[i];

   ofs <<"\n# isndlist:";
   ofs <<"\n"<<isndlist.size();
   for (unsigned int i=0;i!=isndlist.size();i++)
     ofs <<"\n"<<isndlist[i];

   ofs <<"\n# irankrec:";
   ofs <<"\n"<<irankrec.size();
   for (unsigned int i=0;i!=irankrec.size();i++)
     ofs <<"\n"<<irankrec[i];

   ofs <<"\n# iptr:";
   ofs <<"\n"<<iptr.size();
   for (unsigned int i=0;i!=iptr.size();i++)
     ofs <<"\n"<<iptr[i];

   ofs <<"\n# ireclist:";
   ofs <<"\n"<<ireclist.size();
   for (unsigned int i=0;i!=ireclist.size();i++)
     ofs <<"\n"<<ireclist[i];

   ofs <<"\n# End of file.";
 }

template<typename fT, stl_index dim>
bool SAMGp_CommunicationData<fT,dim>::InitializeFromTextFile(const char* fname, int my_rank )
 {
   cout<<"\nSAMGp_CommunicationData::InitializeFromTextFile: Initialise from file "<<fname<<" (rank = "<< my_rank <<")"<< endl;

   my_rank_ = my_rank;

   // clear all vectors first
   iranksnd.clear();
   ipts.clear();
   isndlist.clear();
   irankrec.clear();
   iptr.clear();
   ireclist.clear();

   // open file
   assert (fname != NULL);
   char filename[200];
   strcpy(filename, fname);
   std::ifstream ifs;
   ifs.open(filename);
   if(!ifs.is_open())
    {
      skm_err.notice(CSP_ERROR, "SAMGp_CommunicationData<fT,dim>::InitializeFromTextFile", 
                                "Text file seems to be missing", fname);
      return false;
    }

   char* pch;
   char text_line[256];
   int size;
   char* delim ="\n";
   // read in header
   ifs.getline(text_line, 256);
   ifs.getline(text_line, 256);

   // read in iranksnd
   ifs.getline(text_line, 256);
   ifs.getline(text_line, 256);
   pch = strtok(text_line, delim);
   size = atoi(pch);
   iranksnd.resize(size);
   for (unsigned int i=0;i!=iranksnd.size();i++)
    {
      ifs.getline(text_line, 256);
      pch = strtok(text_line, delim);
      iranksnd[i]=atoi(pch);
    }

   // read in ipts
   ifs.getline(text_line, 256);
   ifs.getline(text_line, 256);
   pch = strtok(text_line, delim);
   size = atoi(pch);
   ipts.resize(size);
   for (unsigned int i=0;i!=ipts.size();i++)
    {
      ifs.getline(text_line, 256);
      pch = strtok(text_line, delim);
      ipts[i]=atoi(pch);
    }

   // read in isndlist
   ifs.getline(text_line, 256);
   ifs.getline(text_line, 256);
   pch = strtok(text_line, delim);
   size = atoi(pch);
   isndlist.resize(size);
   for (unsigned int i=0;i!=isndlist.size();i++)
    {
      ifs.getline(text_line, 256);
      pch = strtok(text_line, delim);
      isndlist[i]=atoi(pch);
    }

   // read in irankrec
   ifs.getline(text_line, 256);
   ifs.getline(text_line, 256);
   pch = strtok(text_line, delim);
   size = atoi(pch);
   irankrec.resize(size);
   for (unsigned int i=0;i!=irankrec.size();i++)
    {
      ifs.getline(text_line, 256);
      pch = strtok(text_line, delim);
      irankrec[i]=atoi(pch);
    }

   // read in iptr
   ifs.getline(text_line, 256);
   ifs.getline(text_line, 256);
   pch = strtok(text_line, delim);
   size = atoi(pch);
   iptr.resize(size);
   for (unsigned int i=0;i!=iptr.size();i++)
    {
      ifs.getline(text_line, 256);
      pch = strtok(text_line, delim);
      iptr[i]=atoi(pch);
    }

   // read in ireclist
   ifs.getline(text_line, 256);
   ifs.getline(text_line, 256);
   pch = strtok(text_line, delim);
   size = atoi(pch);
   ireclist.resize(size);
   for (unsigned int i=0;i!=ireclist.size();i++)
    {
      ifs.getline(text_line, 256);
      pch = strtok(text_line, delim);
      ireclist[i]=atoi(pch);
    }

   ifs.getline(text_line, 256);

   return true;

 }


template class SAMGp_CommunicationData<csp_float,1U>;
template class SAMGp_CommunicationData<csp_float,2U>;
template class SAMGp_CommunicationData<csp_float,3U>;

} // END CSP NAMESPACE
/*
COMMENT

iranksnd	- vector of integers
			- contains rankids of processors to which data should be send
			- iranksnd.size()==npsnd

ipts		- vector of integers
			- needed to determine which entries of isndlist have to be send to which processor
			- ipts.size() == iranksnd.size()+1
			- first entry is always one, next entry is i+1, with isndlist[i] the first element which
			  need to be send to the next negihboring processor, etc..

isndlist	- vector of integers
			- contains column numbers of entries outside the partition, but in the same column range,
              which contain non-zeros. For each neighboring processor only one column number needs to 
              be stored. As the same column could be send to several neighboring processors, this list
              can have several entries with the same value. The local column numbers are stored, so it
              should only contain numbers between the 1 and nnu_local

irankrec	- vector of integers
			- contains rankids of processors from which data should be received
			- iranksnd.size()==nprec 

iptr		- vector of integers
			- needed to determine which entries of ireclist have to be send to which processor
			- iptr.size() == irankrec.size()+1
			- first entry is always one, next entry is i+1, with ireclist[i] the first element which
			  need to be received from the next negihboring processor, etc..

ireclist	- vector of integers
			- contains column numbers of entries outside the partition, but in the same row range,
              which contain non-zeros. For each neighboring processor only one column number needs to 
              be stored. Since a column number can only relate to one neighboring processor, this list 
              can NOT have several entries with the same value. Local column numbers are stored. These
              are counted, starting from nnu_local+1 contiguitively. So this list will always have 
              opeenvolgende nummers.
*/
