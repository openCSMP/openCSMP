#ifndef ECLIPSE_MODEL_H
#define ECLIPSE_MODEL_H

#include "EclipseInterface.h"
#include "Model.h"

namespace csmp {


/**

@brief Input interface which converts a corner-point (hexahedral cell) grid from 
Schlumberger's reservoir simulator Eclipse into a CSMP model object in which
the cells are represented by hexahedra. Degenerate cells are converted to other
element types like prisms and tetrahedra.

@author R. Manasipov
@author Stephan Matthai (refactoring in progress)
@date 2014, 2016

*/
template<size_t dim>
class EclipseModel : public csmp::Model<dim>
{
  public:

      /// CSMP model construction from Eclipse files; property database is created from "variables_file.txt"
      EclipseModel( const std::string& model_name,
                    const std::string& variables_file );

      /// CSMP model construction from mesh files; empty property database
      EclipseModel( const std::string& model_name );

      virtual ~EclipseModel();

      /// MOST IMPORTANT METHOD!- but should be private (TODO: move EclipseModelSettings out of class)
      void BuildModel();

      /// name of the model that will be used when creating output files
      void        Name( const char* );
      const char* Name() const;
 
      /// access to the settings of the Eclipse interface
      EclipseModelSettings&  EclipseModelSetup();

      /// existing special regions
      template<class Container>  void GetRegions( Container& data );
  
      template<class Container>  void GetFaults( Container& data );
  
      template<class Container>  void GetWells ( Container& data );

      /// processing special regions
      void CreateBoundariesAroundFaults( bool keep_fault_regions = false );
      void CreateSplitBoundariesAroundFaults( bool delete_fault_regions = false );

  private:

      // PROPS and other specs from RUNSPECS file
      EclipseModelSettings      eclipse_model_settings_;
      std::string               model_name_;
      std::set<std::string>     regions_;
      std::set<std::string>     faults_;
      std::set<std::string>     wells_;
};

} // csmp

#endif
