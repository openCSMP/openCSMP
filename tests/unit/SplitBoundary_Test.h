#ifndef SPLIT_BOUNDARY_TEST_H
#define SPLIT_BOUNDARY_TEST_H

#include "Test.h"
#include "SplitBoundary.h"

namespace csmp {

/**
    SplitBoundary_Test of the functionality of individual SpitBoundary objects
    as opposed to their creation (see SplitBoundaryInterface_Test) or that of the InterFace
    objects that the SplitBoundary objects consist of (see InterFace_Test).
    
    @attention this also performs test on SplitBoundary objects created in ANSYS
    
    @section Tests Performed Here (2D and 3D)
    
    - do all interfaces have higher dimensional neighbors
    - do correspoding nodes coincide
    - is there only a single node along the perimeter of a split boundary
    - are all of the above operations still valid for a SplitBoundary read from binary file?


     @author Revised testst originally conceived by Roman Manasipov in 2013,
     @author refactored by Junchul Kim (2019).
     @date ported to 2016 version by SKM.
 */
  class SplitBoundary_Test : public Test
    {
    public:
      explicit SplitBoundary_Test( bool verbose=false ) : verbose_(verbose) {}
    
      virtual void run();

      template<size_t dim>
      void test_splitboundary_between_regions( const std::string& model_name );

      template<size_t dim>
      void test_splitboundary_around_regions( const std::string& model_name );

      template<size_t dim>
      void detect_and_create_splitboundaries( const std::string& model_name );

      template<size_t dim>
      void detect_and_create_splitboundaries_from_constructor( const std::string& model_name );

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
      void VisualiseSplitBoundaries( csmp::Model<dim>&, const std::string&test_name  );

      template<size_t dim>
      void TestUnitNormals( csmp::Model<dim>&, const std::string& test_name  );

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
      void etablishContiguosRegionsList( Model<dim>& model,
                                         const std::set<std::string>& fractures_basic_set,
                                         std::set<std::string>& fractures );

      /// reading fault modeling input data
      void inputFromFile( const char* file_name,
                          std::set<std::string>& fractures_basic_set );

      void inputFromFile( const char* file_name,
                          std::vector<std::string>& fractures );

      /// writing fault modeling input data into text file
      void outputToFile( const char* file_name,
                         const std::vector<std::string>& fractures );
      void outputToFile( const char* file_name,
                         const std::set<std::string>& fractures );
      private:
        const bool verbose_;
  };

  } // csmp

#endif
