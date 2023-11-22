#include "ParallelIO.h"

using namespace std;

namespace csp {

extern csp::ErrorHandler           skm_err;

template<typename fT, stl_index dim>
ParallelIO<fT,dim>::ParallelIO(int n_procs, int my_rank)
: n_procs_(n_procs), my_rank_(my_rank)
 {
   if(my_rank_ >= n_procs_)
     skm_err.notice( CSP_ERROR, "ParallelIO::ParallelIO",
                                "\nProcessor rank is larger than total number of processors used in simulation!!.");

 }
 
template<typename fT, stl_index dim>
ParallelIO<fT,dim>::~ParallelIO()
 {
 }

template<typename fT, stl_index dim>
bool ParallelIO<fT,dim>::ReadinVSets(const char* fname, std::vector<VSet<fT,dim> >& vsets, std::vector<stl_index>& first_outerhalo, csp_float& g_time)
 {
   bool read_vset(false);
   vsets.resize(n_procs_);
   first_outerhalo.resize(n_procs_);

   for(int i=0;i!=n_procs_;i++)
    {
      // create name local input file
      std::string  filename(fname);
      char rank[200];
      sprintf(rank, "%ld", static_cast<long>(i));
      filename += "_";
      filename += rank;

      // Read in the VSet file of rank i
      csp_float year(31536000.0);
      stl_index frstohalo;
      cout <<"\nParallelIO::ReadinVSets: Rank "<<my_rank_<<": Restarting from VSet file, reading file "<<filename<<".vset"<< endl;
      read_vset = false;
      read_vset = vsets[i].ParallelInputFrom(filename.c_str(), g_time, frstohalo);
      if(!read_vset)
       {
         // Error message to screen
         cout <<"\nParallelIO::ReadinVSets: Rank "<<my_rank_<<": Could not read VSet from file "<<filename<<".vset, terminating..."<< endl;
         skm_err.notice( CSP_ERROR, "ParallelIO::ReadinVSets",
                                    "\nCould not read VSet file.");

         break;
       }
      else
       {
         first_outerhalo[i] = frstohalo;
         cout <<"\nParallelIO::ReadinVSets: Rank "<<my_rank_<<": Succesfully read file "<<filename<<".vset."<< endl;
         cout <<"\nParallelIO::ReadinVSets: Rank "<<my_rank_<<": Restarting simulation at "<<g_time/year<<" years "<< endl;
         cout <<"\nParallelIO::ReadinVSets: Rank "<<my_rank_<<": First outerhalo "<<frstohalo<< endl;
       }
    }
   return read_vset;
 }


template<typename fT, stl_index dim>
bool ParallelIO<fT,dim>::ReadinLocalVSet(const char* fname, VSet<fT,dim>& loc_vset, stl_index& first_outerhalo, csp_float& g_time)
 {
   bool read_vset(false);

   // create name local input file
   std::string  filename(fname);
   char rank[200];
   sprintf(rank, "%ld", static_cast<long>(my_rank_));
   filename += "_";
   filename += rank;

   // Read in the VSet file of rank i
   csp_float year(31536000.0);
   stl_index frstohalo;
   cout <<"\nParallelIO::ReadinVSets: Rank "<<my_rank_<<": Restarting from VSet file, reading file "<<filename<<".vset"<< endl;
   read_vset = false;
   read_vset = loc_vset.ParallelInputFrom(filename.c_str(), g_time, first_outerhalo);
   if(!read_vset)
    {
      // Error message to screen
      cout <<"\nParallelIO::ReadinLocalVSet: Rank "<<my_rank_<<": Could not read VSet from file "<<filename<<".vset, terminating..."<< endl;
      skm_err.notice( CSP_ERROR, "ParallelIO::ReadinLocalVSet",
                                 "\nCould not read VSet file.");
    }
   else
    {
      cout <<"\nParallelIO::ReadinVSets: Rank "<<my_rank_<<": Succesfully read file "<<filename<<".vset."<< endl;
      cout <<"\nParallelIO::ReadinVSets: Rank "<<my_rank_<<": Restarting simulation at "<<g_time/year<<" years "<< endl;
     }

   return read_vset;
 }

template<typename fT, stl_index dim>
void ParallelIO<fT,dim>::PipeOutputToFile(const char* fname)
 {

   // create name local input file
   std::string  filename(fname);
   char rank[200];
   sprintf(rank, "%ld", static_cast<long>(my_rank_));
   filename += "_";
   filename += rank;

   stdout_filename = filename;

   cout<<"-------------------------------------------------"<< endl;
   cout <<"Rank "<<my_rank_<<" of "<<n_procs_<<": std output piped to file "<<filename.c_str()<< endl;
   cout<<"-------------------------------------------------"<< endl;
   cout <<"\n"<< endl;

   cout.flush();
   // Reopen standard stream, and associate with file
   freopen(filename.c_str(),"w",stdout);

   // Store header
   cout<<"-------------------------------------------------"<< endl;
   cout <<"This file stores all screen output of rank "<<my_rank_<<" of "<<n_procs_<< endl;
   cout<<"-------------------------------------------------"<< endl;
   cout <<"\n"<< endl;

 }



template<typename fT, stl_index dim>
void ParallelIO<fT,dim>::ClearOutputFile()
 {

   // This function clears the std output file and starts a new one with the same name
   // It should only be applied when PipeOutputToFile(const char*) has been called!

   // clear file
   if(remove(stdout_filename.c_str())!=0)
    {
      cout <<"ParallelIO<fT,dim>::ClearOutputFile(): Rank "<<my_rank_<<" of "<<n_procs_<<": file "<<stdout_filename.c_str()<<" could not be deleted, does it exist? \nNothing is done..."<< endl;
      skm_err.notice( WARNING, "ParallelIO::ClearOutputFile",
                               "\nCould not delete output file.");
    }
   else
    {
      cout <<"ParallelIO<fT,dim>::ClearOutputFile(): Rank "<<my_rank_<<" of "<<n_procs_<<": succesfully deleted  "<<stdout_filename.c_str()<<"...creating new empty file."<< endl;
      freopen(stdout_filename.c_str(),"w",stdout);
    }
     

 }



template<typename fT, stl_index dim>
void ParallelIO<fT,dim>::RemoveScreenOutput()
 {
   cout<<"-------------------------------------------------"<< endl;
   cout <<"Stopping all std output for rank "<<my_rank_<<" of "<<n_procs_<<"\n"<< endl;
   cout<<"-------------------------------------------------"<< endl;

   // Reopen standard stream, and associate /dev/null
   freopen("/dev/null","w",stdout);

 }

template<typename fT, stl_index dim>
void ParallelIO<fT,dim>::SaveVSetFiles(SuperGroup<fT,dim>& sg, fT global_time, const char* fname, stl_index first_outerhalo)
 {
   // create name local input file
   std::string  filename(fname);
   char rank[200];
   sprintf(rank, "%ld", static_cast<long>(my_rank_));
   filename += "_";
   filename += rank;

   // screen output
   cout<<"\nParallelIO::SaveVSetFiles: Rank "<<my_rank_<<": Writing data to VSet file "<<filename<<".vset."<< endl;

   // output to vset file
   sg.OutputTo(output_vset);
   output_vset.ParallelOutputTo(filename.c_str(),global_time, first_outerhalo);
 }

template<typename fT, stl_index dim>
void ParallelIO<fT,dim>::WriteToTextFile(SAMGp_CommunicationData<fT,dim>& SAMGp_Comm, const char* fname)
 {
   // create name local input file
   std::string  filename(fname);
   char rank[200];
   sprintf(rank, "%ld", static_cast<long>(my_rank_));
   filename += "_";
   filename += rank;
   filename += "_";
   filename += "SAMGp_Comm.txt";

   // screen output
   cout<<"\nParallelIO::WriteToTextFile: Rank "<<my_rank_<<": Writing SAMGp_CommunicationData to text file "<<filename<< endl;

   // output to txt file
   SAMGp_Comm.DumpToTextFile(filename.c_str());
 }

template<typename fT, stl_index dim>
bool ParallelIO<fT,dim>::ReadFromTextFile(SAMGp_CommunicationData<fT,dim>& SAMGp_Comm, const char* fname)
 {
   // create name local input file
   std::string  filename(fname);
   char rank[200];
   sprintf(rank, "%ld", static_cast<long>(my_rank_));
   filename += "_";
   filename += rank;
   filename += "_";
   filename += "SAMGp_Comm.txt";

   // screen output
   cout<<"\nParallelIO::ReadFromTextFile: Rank "<<my_rank_<<": Initializing SAMGp_CommunicationData from text file "<<filename<< endl;

   // input from txt file
   bool reading_file = SAMGp_Comm.InitializeFromTextFile(filename.c_str(), my_rank_);
   if(!reading_file) 
     skm_err.notice( CSP_ERROR, "ParallelIO::ReadFromTextFile",
                                "\nProblem reading SAMGp_Communication file.");

   return reading_file;
 }

template<typename fT, stl_index dim>
void ParallelIO<fT,dim>::OutputVTKFiles(SuperGroup<fT,dim>& superg, const char* prop, long timestep)
 {
  csp::String variable(prop);
  variable.ReplaceWhiteSpaceBy('_');
  std::string  filename(variable);
  char rank[200];
  sprintf(rank, "%ld", static_cast<long>(my_rank_));
  filename += "_";
  filename += rank;
  filename += "_";

  VTK_Interface<fT,dim> vtk;

  cout<<"\nParallelIO::OutputVTKFiles: Rank "<<my_rank_<<": Outputting "<<prop<<" data to "<<filename<<".vtk."<< endl;

  vtk.OutputDataToVTK( superg, filename.c_str(), prop, timestep );

 }



//template class ParallelIO<csp_float,1U>;
template class ParallelIO<csp_float,2U>;
template class ParallelIO<csp_float,3U>;


} // END CSP NAMESPACE
