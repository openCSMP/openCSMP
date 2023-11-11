#ifndef ParallelIO_h__
#define ParallelIO_h__

#include "CSP_definitions.h"
#include "SuperGroup.h"
#include "CSP_VSet.h"
#include "SAMGp_CommunicationData.h"
#include "CSP_VTK_Interface.h"
namespace csp {

template<typename fT, stl_index dim>
 class ParallelIO
  {
  public:
    ParallelIO(int n_procs, int my_rank);
    ~ParallelIO();

    bool ReadinVSets(const char* , std::vector<VSet<fT,dim> >& , std::vector<stl_index>&, csp_float& );
    bool ReadinLocalVSet(const char* , VSet<fT,dim>& , stl_index&, csp_float& );
    void SaveVSetFiles(SuperGroup<fT,dim>& , fT , const char* , stl_index );
    void WriteToTextFile(SAMGp_CommunicationData<fT,dim>& SAMGp_Comm, const char* fname);
    bool ReadFromTextFile(SAMGp_CommunicationData<fT,dim>& SAMGp_Comm, const char* fname);
    void PipeOutputToFile(const char* fname);
    void RemoveScreenOutput();
    void ClearOutputFile();
    void PipeOutputToScreen();


    void OutputVTKFiles(SuperGroup<fT,dim>&, const char*, long);

  private:
    VSet<fT,dim> output_vset;
    int n_procs_;
    int my_rank_;
    std::string stdout_filename;

  };

} // END CSP NAMESPACE 
#endif
