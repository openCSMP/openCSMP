#ifndef SPLIT_BOUNDARY_INTERFACE_TEST_H
#define SPLIT_BOUNDARY_INTERFACE_TEST_H

#include "Test.h"
#include "SplitBoundaryInterface.h"

namespace csmp {

template<size_t> class Region;
template<size_t> class Face;
template<size_t> class InterFace;
template<size_t> class Model;


/**
    SplitBoundaryInterface_Test unit test tests all functionality related to the creation and destruction of SplitBoundary objects.
    Refer to SplitBoundaryInterface_Test for tests of the functionality of the SplitBoundary itself.

     @author Revised testst originally conceived by Roman Manasipov in 2013,
     @author refactored by Junchul Kim (2019).
     @date ported to 2016 version by SKM.
 */
  class SplitBoundaryInterface_Test : public Test
    {
    public:
      explicit SplitBoundaryInterface_Test( bool verbose=false ) : verbose_(verbose) {}
    
      virtual void run();

      // testing split boundary creation methods 
      
      template<size_t dim>
      void Test_splitboundary_between_regions( const std::string& model_name );

      template<size_t dim>
      void Test_splitboundary_around_regions( const std::string& model_name );

      template<size_t dim>
      void Detect_and_create_splitboundaries( const std::string& model_name );

      template<size_t dim>
      void Detect_and_create_splitboundaries_from_constructor( const std::string& model_name );

      template<size_t dim>
      void NodeParents( const csmp::Region<dim>& region, size_t minParentCount = 2 );

      template<size_t dim>
      void ElementNodes( const csmp::Region<dim>& region );

      template<size_t dim>
      void CheckRemovedLowDimParents( csmp::Boundary<dim>& boundary );

      /// tests whether all InterFace elements of the SplitBoundary have both higher-dimensional neighbors
      template<size_t dim>
      bool NoNeighborNull( const csmp::InterFace<dim>& interFace );

      /// tests whether the splitted nodes share parent elements, located on different sides of interfaces
      template<size_t dim>
      void TestSplitNodeAssignment( const csmp::Model<dim>& );

      template<size_t dim>
      void TestUnitNormals( csmp::Model<dim>&, const std::string& test_name  );

      template<size_t dim>
      void VisualiseSplitBoundaries( csmp::Model<dim>&, const std::string&test_name  );

      template<size_t dim>
      void PullApartSplitboundaries( Model<dim>& model, std::vector<std::string>& fractures, double64 displacement );

      template<size_t dim>
      void PullApartSplitboundaries( Model<dim>& model, double64 displacement );

      template<size_t dim>
      void LoadModel( const std::string& model_name );

      template<size_t dim>
      void LoadModel( const std::string& model_name, std::vector<std::string>& regions );

      template<size_t dim>
      void LoadContiguousModel( const std::string& model_name, std::vector<std::string>& fractures );

      template<size_t dim>
      void EstablishContiguousRegionsList( Model<dim>& model,
                                           const std::set<std::string>& fractures_basic_set,
                                           std::set<std::string>& fractures );

      /// reading fault modeling input data
      void InputFromFile( const char* file_name,
                          std::set<std::string>& fractures_basic_set );

      void InputFromFile( const char* file_name,
                          std::vector<std::string>& fractures );

      /// writing fault modeling input data into text file
      void OutputToFile( const char* file_name,
                         const std::vector<std::string>& fractures );
                         
      void OutputToFile( const char* file_name,
                         const std::set<std::string>& fractures );
      private:
        const bool verbose_;
  };

  } // csmp

#endif /* SPLIT_BOUNDARY_INTERFACE_TEST_H */
